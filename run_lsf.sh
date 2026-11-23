#!/bin/bash
# Run the whole pipeline on an LSF cluster (written for the Janelia cluster).
#
#   bash run_lsf.sh <project config.json> <output dir> [LSF project code]
#
# Submits one array job that counts every library in the sample sheet (merging read pairs with
# fastp first when r2 is given), and one job that waits for it and writes the tables, the tile
# QC report and the heatmap to <output dir>/results/. Safe to run from a login node: it only
# submits jobs.
#
# Settings, through environment variables:
#   JSEQ_ENV       conda environment name or path           (default: junctionseq)
#   JSEQ_CONDA_SH  path to conda's etc/profile.d/conda.sh   (default: found from $CONDA_EXE)
#   JSEQ_SLOTS     slots per counting job                   (default: 4; peak memory is < 1 GB)
#   JSEQ_WALL      wall time per counting job               (default: 2:00)
#   JSEQ_FONT_DIR  folder with Arial .ttf files for the heatmap (optional)
set -eo pipefail
if [ $# -lt 2 ]; then sed -n '2,20p' "$0"; exit 1; fi
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CFG=$(readlink -f "$1")
OUT=$(readlink -m "$2")
PROJ=${3:-${JSEQ_LSF_PROJECT:-}}
ENV=${JSEQ_ENV:-junctionseq}
SLOTS=${JSEQ_SLOTS:-4}
WALL=${JSEQ_WALL:-2:00}
if [ -z "$JSEQ_CONDA_SH" ]; then
  if [ -n "$CONDA_EXE" ]; then JSEQ_CONDA_SH=$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh; fi
fi
[ -f "$JSEQ_CONDA_SH" ] || { echo "cannot find conda.sh; set JSEQ_CONDA_SH"; exit 1; }

python3 -I "$HERE/jseq.py" check "$CFG"          # stops here if a FASTQ is missing
N=$(python3 -I "$HERE/jseq.py" list-samples "$CFG" | wc -l)
mkdir -p "$OUT/logs"
TAG=jseq_$(basename "$CFG" .json)_$$
P=(); [ -n "$PROJ" ] && P=(-P "$PROJ")
SETUP="source '$JSEQ_CONDA_SH' && conda activate '$ENV' && export OMP_NUM_THREADS=\$LSB_DJOB_NUMPROC MKL_NUM_THREADS=\$LSB_DJOB_NUMPROC OPENBLAS_NUM_THREADS=\$LSB_DJOB_NUMPROC JSEQ_FONT_DIR='${JSEQ_FONT_DIR:-}'"

SUB=$(bsub -J "${TAG}[1-$N]" -n "$SLOTS" -W "$WALL" "${P[@]}" -o "$OUT/logs/count_%J_%I.log" \
  "$SETUP && python '$HERE/jseq.py' run-one '$CFG' '$OUT' --index \$LSB_JOBINDEX")
echo "$SUB"
ARRAY=$(sed -n 's/^Job <\([0-9]*\)>.*/\1/p' <<<"$SUB")
[ -n "$ARRAY" ] || { echo "could not read the array job id from: $SUB"; exit 1; }

# numdone(,*) requires EVERY element to have finished successfully. With ended() instead, the
# summary would run after a failed element and could report whatever counts were already on disk
# from an earlier run. If the summary job stays PEND after the array finishes, an element failed:
# look in $OUT/logs/count_*.log, fix it, and rerun this script.
bsub -J "${TAG}_summary" -w "numdone($ARRAY,*)" -n 2 -W 0:30 "${P[@]}" \
  -o "$OUT/logs/summary_%J.log" \
  "$SETUP && python '$HERE/jseq.py' summarize '$CFG' '$OUT'"
echo "submitted $N counting jobs (array $ARRAY) and one summary job (name ${TAG})."
echo "results will be in $OUT/results/ once every counting job has succeeded."
