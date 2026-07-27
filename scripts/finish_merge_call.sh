#!/usr/bin/env bash
#
# Re-do the AutoGVP merge_call step that OOM-killed.
#
# merge_call concatenates the 100 per-chunk tables and sorts them by genomic
# position. The pipeline used `csvtk sort`, which holds the whole table in
# memory and ran out at full-cohort scale. This does the same thing with a
# disk-backed GNU sort: bounded memory, runs on a normal node.
#
# All 100 chunk files already exist, so this is all that is left to finish the run.
# No conda env or modules needed: this is awk + GNU sort only.
#
set -euo pipefail

WORK=/DCEG/CGF/Bioinformatics/Production/Dongjing/PLP_prediction_workflow
CHUNKDIR=$WORK/output_SR0467-047/autogvp
OUTDIR=$WORK/output_SR0467-047/merge_call
THREADS=8
SORT_MEM=8G

# sort spills here; needs roughly the size of the biggest output free (point at
# node-local scratch if you have it, e.g. /lscratch/$SLURM_JOB_ID, for speed)
export TMPDIR=$WORK/TMP
mkdir -p "$OUTDIR" "$TMPDIR"

merge_one() {
  local suffix=$1 out=$2
  local first="$CHUNKDIR/00000-autogvp-annotated-$suffix.tsv"

  # column positions of chr/start/ref/alt (header is identical across chunks)
  local h chrc startc refc altc
  h=$(head -1 "$first")
  chrc=$( printf '%s\n' "$h" | tr '\t' '\n' | grep -nxF 'chr'   | cut -d: -f1)
  startc=$(printf '%s\n' "$h" | tr '\t' '\n' | grep -nxF 'start' | cut -d: -f1)
  refc=$( printf '%s\n' "$h" | tr '\t' '\n' | grep -nxF 'ref'   | cut -d: -f1)
  altc=$( printf '%s\n' "$h" | tr '\t' '\n' | grep -nxF 'alt'   | cut -d: -f1)
  [ -n "$chrc" ] && [ -n "$startc" ] && [ -n "$refc" ] && [ -n "$altc" ] \
    || { echo "ERROR: chr/start/ref/alt not all found in $first" >&2; return 1; }

  # concat (one header) | tag header so it floats to the top | external sort by
  # chr (version) then start (numeric) | drop the tag | append vid = chr:start:ref:alt
  awk -F'\t' 'FNR==1 && NR>1 {next} {print}' "$CHUNKDIR"/000{00..99}-autogvp-annotated-$suffix.tsv \
    | awk -F'\t' 'NR==1{print "0\t"$0; next}{print "1\t"$0}' \
    | LC_ALL=C sort -t$'\t' -k1,1n -k$((chrc+1)),$((chrc+1))V -k$((startc+1)),$((startc+1))n \
        -S "$SORT_MEM" -T "$TMPDIR" --parallel="$THREADS" \
    | cut -f2- \
    | awk -F'\t' -v OFS='\t' -v c="$chrc" -v s="$startc" -v r="$refc" -v a="$altc" \
        'NR==1{print $0,"vid"; next}{print $0,$c":"$s":"$r":"$a}' \
    > "$out"

  echo "$suffix -> $out  ($(($(wc -l < "$out") - 1)) rows)"
}

# discover whatever chunk variants exist (full, abridged, ...) and merge each
shopt -s nullglob
for f in "$CHUNKDIR"/00000-autogvp-annotated-*.tsv; do
  b=$(basename "$f"); sfx=${b#00000-autogvp-annotated-}; sfx=${sfx%.tsv}
  merge_one "$sfx" "$OUTDIR/SR0467-047.autogvp_${sfx}.tsv"
done

echo "Done. Final tables in $OUTDIR/"
