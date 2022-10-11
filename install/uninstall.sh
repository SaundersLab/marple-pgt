#!/bin/bash

# Bash strict mode
set -euo pipefail

install_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
marple_pgt_dir=$(dirname $install_dir)

if [ -f "$marple_pgt_dir"/marple_pgt_miniconda/bin/activate ] ; then
    source "$marple_pgt_dir"/marple_pgt_miniconda/bin/activate
    for i in $(seq ${CONDA_SHLVL}); do conda deactivate ; done
    conda env remove -n marple-pgt
fi

[ -f "$install_dir"/marple_pgt_miniconda.sh ] && rm "$install_dir"/marple_pgt_miniconda.sh
[ -d "$marple_pgt_dir"/marple_pgt_miniconda ] && rm -rf "$marple_pgt_dir"/marple_pgt_miniconda
