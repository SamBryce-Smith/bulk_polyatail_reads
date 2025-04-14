# Extract reads with soft-clipped polyA tails from bulk RNA-seq

[![DOI](https://zenodo.org/badge/837292050.svg)](https://doi.org/10.5281/zenodo.15210305)

Snakemake pipeline to extract polyA-tail containiing reads (PATRs) from aligned bulk RNA-seq datasets and identify poly(A) site clusters. Broadly mimics workflow reported by [Vlasenok et al. 2023](https://doi.org/10.1093/nargab/lqad051) with more customisation options (e.g. rescuing/using shorter overhang lengths)

## Prerequesites

Input files:

- RNA-seq BAM files aligned with soft-clipping enabled. Indexes must be present in the same location (suffixed with .bai)
- BED file containing PAS to count (required for plumbing purposes, can use provided file if not interested in specific set of PAS)

Software:

- Snakemake
- Conda/mamba

Currently, you will need to manage initial dependencies yourself (this may change in the future). Pipeline dependencies are managed using a conda environment.

## Configuration

- Generate a config YAML file using `config.yaml` as a template. All parameters are described using comments in the config file
- Generate a CSV sample table to map BAM files to sample names and groups. A complete example can be found at [example_sample_table.csv](example_sample_table.csv). The minimal columns are `sample_name` (unique identifier for sample) and `bam` (path to BAM file)

## Usage

### TODO: provide example data for testing/initial run

Dry run:

```bash
snakemake -n -p --configfile <config_file.yaml> 
```

Local run:

```bash
snakemake -n -p --configfile <config_file.yaml> --use-conda --cores <n>
```

replacing `<config_file.yaml>` with the path to the config file and `<n>` with the number of cores for parallel processing of samples.

UCL CS cluster users can use the `submit.sh` bash script to submit the pipeline to for remote processing on the cluster. See `bash submit.sh -h` for details on how to use.

## Output

```bash
bulk_polya_reads
├── NT_19074709_S20.polya_tail_reads.alignment_stats.tsv
├── NT_19074709_S20.polya_tail_reads.global_softclip_tail_counts.tsv
├── NT_19074709_S20.polya_tail_reads.parquet
├── TDP43_19065407_S25.polya_tail_reads.alignment_stats.tsv
├── TDP43_19065407_S25.polya_tail_reads.global_softclip_tail_counts.tsv
├── TDP43_19065407_S25.polya_tail_reads.parquet
├── all_samples.pas.count_matrix.tsv
├── benchmark
├── logs
├── pas_clusters
│   ├── all
│   │   └── two_class_simple
│   │       └── all_samples.polya_clusters.bed
│   ├── batch__seddighi_i3_cortical_batch1___condition__NT
│   │   └── two_class_simple
│   │       └── all_samples.polya_clusters.bed
│   ├── batch__seddighi_i3_cortical_batch1___condition__TDP43KD
│   │   └── two_class_simple
│   │       └── all_samples.polya_clusters.bed
│   ├── batch__seddighi_i3_cortical_batch2___condition__NT
│   │   └── two_class_simple
│   │       └── all_samples.polya_clusters.bed
│   ├── batch__seddighi_i3_cortical_batch2___condition__TDP43KD
│   │   └── two_class_simple
│   │       └── all_samples.polya_clusters.bed
│   ├── condition__NT
│   │   └── two_class_simple
│   │       └── all_samples.polya_clusters.bed
│   ├── condition__TDP43KD
│   │   └── two_class_simple
│   │       └── all_samples.polya_clusters.bed
```


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
