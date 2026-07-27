#!/usr/bin/env bash
#
# Extract Pathogenic / Likely pathogenic rows from one family's 2caller file.
# One family per Slurm array task. Pure awk, no modules needed.
#
set -euo pipefail
M=/DCEG/Projects/WGS/Build/SR0467-047/ensemble/AutoGVP/by_family

mapfile -t FILES < <(ls "$M"/*.autogvp_abridged.2caller.tsv | sort)
i=${SLURM_ARRAY_TASK_ID:-0}
[ "$i" -lt "${#FILES[@]}" ] || { echo "index $i >= ${#FILES[@]} files" >&2; exit 1; }

f=${FILES[$i]}
fam=$(basename "$f" .autogvp_abridged.2caller.tsv)
out="$M/$fam.PLP.2caller.tsv"

col=$(head -1 "$f" | tr '\t' '\n' | grep -nxF autogvp_call | cut -d: -f1)
[ -n "$col" ] || { echo "autogvp_call not found in $f" >&2; exit 1; }

awk -F'\t' -v c="$col" '
  NR==1 { print; next }
  $c=="Pathogenic" || $c=="Likely_pathogenic" || $c=="Pathogenic/Likely_pathogenic"
' "$f" > "$out"

echo "$fam: $(($(wc -l < "$out")-1)) P/LP variants -> $out"
