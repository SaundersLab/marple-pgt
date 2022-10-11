#!/bin/bash

# Bash strict mode
set -euo pipefail

src_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
marple_pgt_dir=$( dirname "$src_dir" )
source "$marple_pgt_dir"/marple_pgt_miniconda/bin/activate marple-pgt
python3 "$src_dir"/reads_to_consensus.py $@
