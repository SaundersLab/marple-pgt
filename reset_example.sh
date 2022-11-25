#!/bin/bash

[ -z $SLURM_JOBID ] && is_on_slurm=false || is_on_slurm=true

$is_on_slurm || script=$0
$is_on_slurm && script=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | cut -d" " -f1)
marple_pgt_dir=$(cd $(dirname "$script"); pwd)

# Bash strict mode
set -euo pipefail

[ -d "$marple_pgt_dir"/example ] && rm -r "$marple_pgt_dir"/example
tar -xzf "$marple_pgt_dir"/.example_backup.tar.gz 
