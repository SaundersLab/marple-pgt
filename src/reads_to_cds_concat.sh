#!/bin/bash

[ -z $SLURM_JOBID ] && is_on_slurm=false || is_on_slurm=true

$is_on_slurm || script=$0
$is_on_slurm && script=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | cut -d" " -f1)
src_dir=$(cd $(dirname "$script"); pwd)

# Bash strict mode
set -euo pipefail

marple_pgt_dir=$( dirname "$src_dir" )
source "$marple_pgt_dir"/marple_pgt_miniconda/bin/activate marple-pgt
python3 "$src_dir"/reads_to_cds_concat.py $@
