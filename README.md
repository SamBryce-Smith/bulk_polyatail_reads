# Extract reads with soft-clipped polyA tails from bulk RNA-seq

Snakemake pipeline to extract polyA-tail containiing reads (PATRs) from aligned bulk RNA-seq datasets and identify poly(A) site clusters. Broadly mimics workflow reported by [Vlasenok et al. 2023](https://doi.org/10.1093/nargab/lqad051) with more customisation options (e.g. rescuing/using shorter overhang lengths)

## Prerequesites

Input files:

- RNA-seq BAM files aligned with soft-clipping enabled. Indexes must be present in the same location (suffixed with .bai)
- BED file containing PAS to count (required for plumbing purposes, can use provided file `data/query_pas.bed` if not interested in specific set of PAS)

Software:

- Snakemake (tested with 7.32.4)
- Conda/mamba

Currently, you will need to manage initial dependencies yourself (this may change in the future). Alternatively, you can use the development environment `env_dev_bulk_polya_reads.yaml`:

```bash
# conda / mamba depending on preference
conda env create -f env_dev_bulk_polya_reads.yaml
```

Pipeline dependencies are managed using a conda environment.

## Configuration

- Generate a config YAML file using `config.minimal.yaml` or `config.full.yaml` as a template. All parameters are described using comments in the config file
- Generate a CSV sample table to map BAM files to sample names and groups. The minimal columns are `sample_name` (unique identifier for sample) and `bam` (path to BAM file). A complete example can be found at [example_sample_table.csv](example_sample_table.csv).

## Usage

Dry run using example data:

```bash
snakemake -n -p --configfile <config_file.yaml>
```

Local run:

```bash
snakemake -n -p --configfile <config_file.yaml> --use-conda --cores <n>
```

replacing `<config_file.yaml>` with the path to the config file (`config.minimal.yaml` or `config.full.yaml` for complete functionality) and `<n>` with the number of cores for parallel processing of samples.

UCL CS cluster users can use the `submit.sh` bash script to submit the pipeline to for remote processing on the cluster. See `bash submit.sh -h` for details on how to use.

## Output

```bash
# output of config.full.yaml (manually cleaned up .parquet, benchmark and log directories)
$ tree -L 4 example_output_full/
example_output_full/
├── NT_19074709_S20.polya_tail_reads.alignment_stats.tsv
├── NT_19074709_S20.polya_tail_reads.global_softclip_tail_counts.tsv
├── NT_19074709_S20.polya_tail_reads.parquet
├── NT_19074717_S21.polya_tail_reads.alignment_stats.tsv
├── NT_19074717_S21.polya_tail_reads.global_softclip_tail_counts.tsv
├── NT_19074717_S21.polya_tail_reads.parquet
├── TDP43_19065403_S23.polya_tail_reads.alignment_stats.tsv
├── TDP43_19065403_S23.polya_tail_reads.global_softclip_tail_counts.tsv
├── TDP43_19065403_S23.polya_tail_reads.parquet
├── TDP43_19065407_S25.polya_tail_reads.alignment_stats.tsv
├── TDP43_19065407_S25.polya_tail_reads.global_softclip_tail_counts.tsv
├── TDP43_19065407_S25.polya_tail_reads.parquet
├── all_samples.pas.count_matrix.tsv
├── benchmark
├── duplicate_TDP43_19065407_S25.polya_tail_reads.alignment_stats.tsv
├── duplicate_TDP43_19065407_S25.polya_tail_reads.global_softclip_tail_counts.tsv
├── duplicate_TDP43_19065407_S25.polya_tail_reads.parquet
├── logs
└── pas_clusters
    ├── all
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── batch__batch1
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── batch__batch1___condition__NT
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── batch__batch1___condition__TDP43KD
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── batch__batch2
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── batch__batch2___condition__NT
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── batch__batch2___condition__TDP43KD
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── condition__NT
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    ├── condition__TDP43KD
    │   └── two_class_simple
    │       └── all_samples.polya_clusters.bed
    └── per_sample
        ├── all_samples.patrs.valid_stats.tsv
        ├── patr_bed2stats.tsv
        └── two_class_simple
            ├── NT_19074709_S20.all_samples.polya_clusters.bed
            ├── NT_19074717_S21.all_samples.polya_clusters.bed
            ├── TDP43_19065403_S23.all_samples.polya_clusters.bed
            ├── TDP43_19065407_S25.all_samples.polya_clusters.bed
            └── duplicate_TDP43_19065407_S25.all_samples.polya_clusters.bed

36 directories, 108 files
```

- `<sample_name>.polya_tail_reads.parquet` - per-sample parquet files containing extracted soft-clipped read statistics (length, sequence, A content, genomic coordinates etc.)
- `all_samples.pas.count_matrix.tsv` - PAS x sample count matrix for input BED file of target PAS (e.g. `data/query_pas.bed`)
- `benchmarks`, `logs` - directories storing per-rule, per invocation runtime statistics (benchmark) and STDOUT/STDERR logs (logs) from Snakemake's benchmark and log directives
- `pas_clusters` - BED6 files of inferred PAS clusters after filtering for valid poly(A)-tail containing reads. Score field contains raw counts
  - Subdirectories are for different pooling methods specified in config file, named as `key__value` and separated by `___` for multiple groupings
  - `per_sample` contains per-sample PAS clusters (if `cluster_per_sample: True`)
    - `all_samples.patrs.valid_stats.tsv` - summary statistics per-sample valid poly(A)-tail containing reads after filtering criteria
    - `patr_bed2stats.tsv` - intermediate table containing paths to alignment statistics and PATR cluster BEDs per-sample

## Assorted Notes

### Limitations

- All extraction code runs single threaded and is not particularly well optimised. The computation times are therefore dependent on the read depth of the samples (and generally have a lot of room for improvement).
- Several papers have reported alignment errors confounding PATR extraction (Vlasenok et al. 2023, Moon et al. 2023). I haven't implemented anything to account for this. This is essential if one wanted to use this pipeline as basis for global PAS definition (as opposed to my use case for validation), particularly for use of shorter length overhangs.

### Expected disk usage for full datasets

Reporting total disk space usage individually for all parquet files in a directory (.parquets are the hungriest outputs in terms of disk space)

```bash
du -h --summarize *.parquet
```

Using the Seddighi i3Neurons, found output size ranged from 250 MB to 441 MB (default chromosome partitioning, outputting all soft-clipped reads length >=1 and any A content). These are pretty high depth RNA-seq samples, so take as a rough guide of expected disk usage per-sample. More conservative thresholds for overhang length or tail content can reduce the size of these files quite substantially.

Cluster BED files with the standard filtering criteria tend to be in the low MB range (e.g. BED with 62,376 clusters was 3.3 MB)
