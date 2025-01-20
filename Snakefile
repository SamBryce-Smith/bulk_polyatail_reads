import pandas as pd
import os
from typing import List

sample_tbl = pd.read_csv(config["sample_table"])

assert all(col in sample_tbl.columns for col in ["sample_name", "bam"])

if config["strandedness_from_sample_table"]:
    assert "strandedness" in sample_tbl.columns, f"'strandedness' column must be present in sample table"

# Check that each input sample has a BAM index (.bai)
for sample_name, bam in zip(sample_tbl["sample_name"], sample_tbl["bam"]):
    assert os.path.exists(os.path.join(bam + ".bai")), f".bai index file does not exist at same location as input BAM file for sample - {sample_name} - inferred path: {bam + '.bai'}"

SAMPLES = sample_tbl["sample_name"].tolist()

# get dict mapping pooling strategy ('pool_id') to corresponding sample names
assert isinstance(config["pooling_methods"], list)
pool2sample = {}

for pooling_method in config["pooling_methods"]:
    if pooling_method == "all":
        pool2sample[pooling_method] = sample_tbl["sample_name"].tolist()

    else:
        # standardise object type (facilitates id construction + column checking)
        pooling_method = [pooling_method] if isinstance(pooling_method, str) else pooling_method
        
        for col in pooling_method:
            assert col in sample_tbl.columns

        grpd = sample_tbl.groupby(pooling_method)
        for key, df in grpd:
            
            # single column groups = group key returned as string. Convert to iterable type so zipping doesn't unpack string
            key = [key] if isinstance(key,str) else key

            # column1__value___column2__value___n
            pool_id = "___".join(["__".join(col2val) for col2val in zip(pooling_method, key)])
            samples = df["sample_name"].tolist()
            pool2sample[pool_id] = samples


assert config["pas_filter_method"] in "vlasenok,two_class_simple,vlasenok_two_class".split(",")


out_dir = config["output_dir"]
output_prefix = config["output_prefix"]

if not os.path.exists(out_dir):
    os.system(f"mkdir -p {out_dir}")

if len(SAMPLES) <= 25:
    sys.stderr.write(f"Basenames for input BAM files - {', '.join(SAMPLES)}\n")

sys.stderr.write(f"Inferred groups for pooling samples: {', ' .join(pool2sample.keys())}\n")


wildcard_constraints:
    sample = "|".join(SAMPLES),
    pool_id = "|".join(pool2sample.keys())

#
localrules: make_pas_stats_table

def per_sample_targets(sample_names: List[str]) -> List[str]:
    beds = expand(os.path.join(out_dir, "pas_clusters", "per_sample", config["pas_filter_method"], ".".join(["{sample}", config["pas_clusters_filename_prefix"], "bed"])), sample=SAMPLES)
    # add stats file to targets and return
    return beds + [os.path.join(out_dir, "pas_clusters", "per_sample", "all_samples.patrs.valid_stats.tsv")]

rule all:
    input: 
        expand(os.path.join(out_dir, "{sample}" + output_prefix + ".global_softclip_tail_counts.tsv"), sample=SAMPLES),
        os.path.join(out_dir, "all_samples.pas.count_matrix.tsv"),
        expand(os.path.join(out_dir, "pas_clusters", "{pool_id}", config["pas_filter_method"], config["pas_clusters_filename_prefix"] + ".bed"), pool_id=pool2sample.keys()),
        per_sample_targets(SAMPLES) if config["cluster_per_sample"] else []


rule extract_polya_reads:
    input:
        bam=lambda wildcards: sample_tbl.loc[sample_tbl["sample_name"] == wildcards.sample, "bam"]

    output: 
        parquet=directory(os.path.join(out_dir, "{sample}" + output_prefix + ".parquet")),
        aln_stats=os.path.join(out_dir, "{sample}" + output_prefix + ".alignment_stats.tsv"),
        sclip_stats=os.path.join(out_dir, "{sample}" + output_prefix + ".global_softclip_tail_counts.tsv")

    params:
        out_prefix=os.path.join(out_dir, "{sample}" + output_prefix),
        strandedness=lambda wildcards: sample_tbl.loc[sample_tbl["sample_name"] == wildcards.sample, "strandedness"].iloc[0] if config["strandedness_from_sample_table"] else config["strandedness"],
        no_output_bam="--no-output-bam" if config["no_output_bam"] else "",
        min_clip_len=config["min_softclip_length"],
        min_fracA=config["min_fracA"]

    conda:
        "env_bulk_polya_reads.yaml"

    log:
        stdout = os.path.join(out_dir, "logs", "{sample}.extract_polya_reads.log.stdout.txt"),
        stderr = os.path.join(out_dir, "logs", "{sample}.extract_polya_reads.log.stderr.txt")
    benchmark:
        os.path.join(out_dir, "benchmark", "{sample}.extract_polya_reads.benchmark.txt")

    shell:
        """
        python scripts/extract_polya_reads.py \
        -b {input.bam} \
        -s {params.strandedness} \
        -m {params.min_clip_len} \
        -f {params.min_fracA} \
        {params.no_output_bam} \
        -o {params.out_prefix} \
        1> {log.stdout} \
        2> {log.stderr}
        """


if config["count_per_sample"]:
    
    rule count_pas_per_sample:
        input:
            parquet = os.path.join(out_dir, "{sample}" + output_prefix + ".parquet"),
            bed = config["pas_bed"]

        output:
            count_mtx = temp(os.path.join(out_dir, "pas_counts_per_sample", "{sample}.pas.count_matrix.tsv"))

        params:
            pas_window = config["pas_window_size"],
            parquet_suffix = config["output_prefix"] + ".parquet", # not modifiable by scripts/extract_polya_reads.py
            output_prefix = os.path.join(out_dir, "pas_counts_per_sample", "{sample}.pas")

        conda:
            "env_bulk_polya_reads.yaml"

        log:
            stdout = os.path.join(out_dir, "logs", "{sample}.count_pas.log.stdout.txt"),
            stderr = os.path.join(out_dir, "logs", "{sample}.count_pas.log.stderr.txt")

        benchmark:
            os.path.join(out_dir, "benchmark", "{sample}.count_pas.benchmark.txt")

        shell:
            """
            python scripts/count_pas.py \
            -i {input.parquet} \
            -b {input.bed} \
            -w {params.pas_window} \
            -s {params.parquet_suffix} \
            -o {params.output_prefix} \
            1> {log.stdout} \
            2> {log.stderr}
            """ 


    rule merge_pas_counts:
        input:
            mtxs = lambda wildcards: expand(os.path.join(out_dir, "pas_counts_per_sample", "{sample}.pas.count_matrix.tsv"), sample=SAMPLES),

        output:
            count_mtx = os.path.join(out_dir, "all_samples.pas.count_matrix.tsv")

        conda:
            "env_bulk_polya_reads.yaml"

        log:
            stdout = os.path.join(out_dir, "logs", "merge_pas_counts.log.stdout.txt"),
            stderr = os.path.join(out_dir, "logs", "merge_pas_counts.log.stderr.txt")

        benchmark:
            os.path.join(out_dir, "benchmark", "merge_pas_counts.benchmark.txt")

        shell:
            """
            python scripts/merge_pas_counts.py \
            --o {output} \
            {input.mtxs}
            """

else:
    # run counting script passing all samples at once
    rule count_pas:
        input:
            parquets = lambda wildcards: expand(os.path.join(out_dir, "{sample}" + output_prefix + ".parquet"), sample=SAMPLES),
            bed = config["pas_bed"]

        output:
            count_mtx = os.path.join(out_dir, "all_samples.pas.count_matrix.tsv")

        params:
            pas_window = config["pas_window_size"],
            parquet_suffix = config["output_prefix"] + ".parquet", # not modifiable by scripts/extract_polya_reads.py
            output_prefix = os.path.join(out_dir, "all_samples.pas")

        conda:
            "env_bulk_polya_reads.yaml"

        log:
            stdout = os.path.join(out_dir, "logs", "count_pas.log.stdout.txt"),
            stderr = os.path.join(out_dir, "logs", "count_pas.log.stderr.txt")

        benchmark:
            os.path.join(out_dir, "benchmark", "count_pas.benchmark.txt")

        shell:
            """
            python scripts/count_pas.py \
            -i {input.parquets} \
            -b {input.bed} \
            -w {params.pas_window} \
            -s {params.parquet_suffix} \
            -o {params.output_prefix} \
            1> {log.stdout} \
            2> {log.stderr}
            """ 


rule cluster_pas:
    input:
        parquets = lambda wildcards: expand(os.path.join(out_dir, "{sample}" + output_prefix + ".parquet"), sample=pool2sample[wildcards.pool_id])

    output:
        os.path.join(out_dir, "pas_clusters", "{pool_id}", config["pas_filter_method"], config["pas_clusters_filename_prefix"] + ".bed")

    params:
        filter_method = config["pas_filter_method"],
        cluster_window = config["pas_cluster_window"],
        parquet_suffix = config["output_prefix"] + ".parquet", # not modifiable by scripts/extract_polya_reads.py
        output_prefix = os.path.join(out_dir, "pas_clusters", "{pool_id}", config["pas_filter_method"], config["pas_clusters_filename_prefix"]),
        long_overhang_min_length = config["long_overhang_min_length"],
        long_overhang_fraction_a = config["long_overhang_fraction_a"],
        short_overhang_length_range = config["short_overhang_length_range"],
        short_overhang_fraction_a = config["short_overhang_fraction_a"],
        entropy_cutoff = config["entropy_cutoff"]

    conda:
        "env_bulk_polya_reads.yaml"

    log:
        stdout = os.path.join(out_dir, "logs", "cluster_pas." + "{pool_id}" + ".log.stdout.txt"),
        stderr = os.path.join(out_dir, "logs", "cluster_pas." + "{pool_id}" + ".log.stderr.txt")

    benchmark:
        os.path.join(out_dir, "benchmark", "cluster_pas." + "{pool_id}" + ".benchmark.txt")

    shell:
        """
        python scripts/cluster_pas.py \
        -i {input.parquets} \
        -f {params.filter_method} \
        -w {params.cluster_window} \
        -s {params.parquet_suffix} \
        --long-overhang-min-length {params.long_overhang_min_length} \
        --short-overhang-length-range {params.short_overhang_length_range} \
        --long-overhang-fraction-a {params.long_overhang_fraction_a} \
        --short-overhang-fraction-a {params.short_overhang_fraction_a} \
        --entropy-cutoff {params.entropy_cutoff} \
        -o {params.output_prefix} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule cluster_pas_per_sample:
    input:
        parquet=rules.extract_polya_reads.output.parquet

    output:
        os.path.join(out_dir, "pas_clusters", "per_sample", config["pas_filter_method"], ".".join(["{sample}", config["pas_clusters_filename_prefix"], "bed"]))

    params:
        filter_method = config["pas_filter_method"],
        cluster_window = config["pas_cluster_window"],
        parquet_suffix = config["output_prefix"] + ".parquet", # not modifiable by scripts/extract_polya_reads.py
        output_prefix = os.path.join(out_dir, "pas_clusters", "per_sample", config["pas_filter_method"], ".".join(["{sample}", config["pas_clusters_filename_prefix"]])),
        long_overhang_min_length = config["long_overhang_min_length"],
        long_overhang_fraction_a = config["long_overhang_fraction_a"],
        short_overhang_length_range = config["short_overhang_length_range"],
        short_overhang_fraction_a = config["short_overhang_fraction_a"],
        entropy_cutoff = config["entropy_cutoff"]

    conda:
        "env_bulk_polya_reads.yaml"

    log:
        stdout = os.path.join(out_dir, "logs", "cluster_pas.per_sample." + "{sample}" + ".log.stdout.txt"),
        stderr = os.path.join(out_dir, "logs", "cluster_pas.per_sample." + "{sample}" + ".log.stderr.txt")

    benchmark:
        os.path.join(out_dir, "benchmark", "cluster_pas.per_sample." + "{sample}" + ".benchmark.txt")

    shell:
        """
        python scripts/cluster_pas.py \
        -i {input.parquet} \
        -f {params.filter_method} \
        -w {params.cluster_window} \
        -s {params.parquet_suffix} \
        --long-overhang-min-length {params.long_overhang_min_length} \
        --short-overhang-length-range {params.short_overhang_length_range} \
        --long-overhang-fraction-a {params.long_overhang_fraction_a} \
        --short-overhang-fraction-a {params.short_overhang_fraction_a} \
        --entropy-cutoff {params.entropy_cutoff} \
        -o {params.output_prefix} \
        1> {log.stdout} \
        2> {log.stderr}
        """

rule make_pas_stats_table:
    '''Create TSV mapping samples to stats bed (total alignments) and clusters bed (valid PATRs)'''
    input:
        # Expand wildcards to get all sample files
        stats=expand(os.path.join(out_dir, "{sample}" + output_prefix + ".alignment_stats.tsv"),
                    sample=SAMPLES),
        bed=expand(os.path.join(out_dir, "pas_clusters", "per_sample", config["pas_filter_method"],
                               ".".join(["{sample}", config["pas_clusters_filename_prefix"], "bed"])),
                  sample=SAMPLES)
    output:
        tsv=os.path.join(out_dir, "pas_clusters", "per_sample", "patr_bed2stats.tsv")
    
    run:
        import csv
        import os
        with open(output.tsv, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['sample_name', 'stats_path', 'bed_path'])
            
            for sample, stats_path, bed_path in zip(SAMPLES, input.stats, input.bed):
                writer.writerow([sample, stats_path, bed_path])


rule per_sample_pas_stats:
    input:
        tsv=rules.make_pas_stats_table.output.tsv
    output:
        os.path.join(out_dir, "pas_clusters", "per_sample", "all_samples.patrs.valid_stats.tsv")
    log:
        stdout = os.path.join(out_dir, "logs", "per_sample_pas_stats.log.stdout.txt"),
        stderr = os.path.join(out_dir, "logs", "per_sample_pas_stats.log.stderr.txt")
    benchmark:
        os.path.join(out_dir, "benchmark", "per_sample_pas_stats.benchmark.txt")
    shell:
        """
        python scripts/per_sample_pas_stats.py \
        {input.tsv} {output} \
        1> {log.stdout} \
        2> {log.stderr}
        """



