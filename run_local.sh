#!/bin/bash
# Run the whole pipeline on one machine, one library after another (no cluster needed).
# Activate the junctionseq conda environment first.
#
#   bash run_local.sh <project config.json> <output dir>
set -eo pipefail
if [ $# -lt 2 ]; then sed -n '2,6p' "$0"; exit 1; fi
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python "$HERE/jseq.py" check "$1"
for lib in $(python "$HERE/jseq.py" list-samples "$1"); do
  python "$HERE/jseq.py" run-one "$1" "$2" --library "$lib"
done
python "$HERE/jseq.py" summarize "$1" "$2"
