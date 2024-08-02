#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


if [[ ( $@ == "--help") ||  $@ == "-h" ]]; then
    echo "Usage: bash submit.sh CONFIG_FILE RUN_NAME"
    echo "CONFIG_FILE - Path to config file for run "
    echo "RUN_NAME - Optional argument to name run. Config file for run will be copied to folder containing cluster log files (submissions/<date><time>/) with run name prefixed"
    echo "-h/--help - print this help message and exit"
    exit 0
fi

if [ "$2" != "" ]; then
    RUN_NAME="$2"
else
    RUN_NAME="bulk_polya_reads"
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}
cp $1 ${FOLDER}/${RUN_NAME}_${1}

snakemake \
--configfile $1 \
--use-conda \
--jobscript cluster_qsub.sh \
--cluster-config cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 50 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
--keep-going
