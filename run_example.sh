#!/bin/bash

[ -z $SLURM_JOBID ] && is_on_slurm=false || is_on_slurm=true

$is_on_slurm || script=$0
$is_on_slurm && script=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | cut -d" " -f1)
marple_pgt_dir=$(cd $(dirname "$script"); pwd)

"$marple_pgt_dir"/reset_example.sh
cd "$marple_pgt_dir"/example
../src/reads_to_cds_concat.sh --trim no */*.fastq
../src/cds_concat_to_tree_imgs.sh \
    --start ../data/12_samples_276_genes_cds.fa.gz \
    --out_dir tree \
    --name new_samples \
    */*_cds_concat.fasta
