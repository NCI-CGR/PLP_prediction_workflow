#!/usr/bin/env bash
#
# Split the AutoGVP result table into per-family files. For each family it writes
# two TSVs, both with each member's full genotype in the ensemble VCF FORMAT layout
#   ...annotation... | ensemble_FILTER | FORMAT | <member1> | <member2> | ...
#
#   <family>.<tag>.all.tsv       all variants the family carries (any caller count)
#   <family>.<tag>.2caller.tsv   only variants called by >=2 callers (twoCallers/threeCallers)
#
# One family per Slurm array task. Reuses the concatenated genotype VCF if present.
# bcftools comes from the module (no conda active, so it loads cleanly).
#
set -euo pipefail
module load bcftools

FAM_DIR=/DCEG/Projects/WGS/Build/SR0467-047/GenCompass_Processing/family_annotations/samples_in_fam
ENS_DIR=/DCEG/Projects/WGS/Build/SR0467-047/ensemble
REF=/DCEG/CGF/Bioinformatics/Production/data/hg38/Homo_sapiens_assembly38.fasta
AUTOGVP=/DCEG/CGF/Bioinformatics/Production/Dongjing/PLP_prediction_workflow/output_SR0467-047/merge_call/SR0467-047.autogvp_abridged.tsv
OUTDIR=/DCEG/Projects/WGS/Build/SR0467-047/ensemble/AutoGVP/by_family
THREADS=8
GENO_ALL=$OUTDIR/SR0467-047.allchr.GT.vcf.gz

mkdir -p "$OUTDIR"
reader() { case "$AUTOGVP" in *.gz) zcat "$AUTOGVP";; *) cat "$AUTOGVP";; esac; }

if [ ! -f "$GENO_ALL" ]; then
  ls "$ENS_DIR"/SR0467-047.chr*.all_callers_merged_genotypes.vcf.gz | sort -V > "$OUTDIR/chr_vcfs.txt"
  bcftools concat -f "$OUTDIR/chr_vcfs.txt" -Oz --threads "$THREADS" -o "$GENO_ALL"
  tabix -p vcf "$GENO_ALL"
fi

IDCOL=$(reader | head -1 | tr '\t' '\n' | grep -nxF 'variant_ids' | cut -d: -f1 || true)
[ -n "${IDCOL:-}" ] || { echo "ERROR: variant_ids not found in $AUTOGVP" >&2; exit 1; }
TAG=$(basename "$AUTOGVP" .tsv); TAG=${TAG#SR0467-047.}

process_family() {
  local FAM_LST=$1
  local FAM; FAM=$(basename "$FAM_LST" .lst); FAM=${FAM%_members}
  local allf="$OUTDIR/$FAM.$TAG.all.tsv"
  local twof="$OUTDIR/$FAM.$TAG.2caller.tsv"
  { [ -s "$allf" ] && [ -s "$twof" ]; } && { echo "[$(date +%T)] $FAM done, skip"; return 0; }
  echo "[$(date +%T)] $FAM"

  # one scan: all variants the family carries (only star alleles removed),
  # capturing FILTER (col 7), FORMAT (col 9) and each member's full value
  bcftools view -S "$FAM_LST" --force-samples -e 'ALT="*"' -Ou "$GENO_ALL" \
    | bcftools norm -m-both -Ou \
    | bcftools norm -f "$REF" -Ou \
    | bcftools view -i 'GT[*]="alt"' -Ov \
    | awk -F'\t' -v OFS='\t' '
        /^##/ { next }
        /^#CHROM/ { s="variant_ids" OFS "ensemble_FILTER" OFS "FORMAT"; for(i=10;i<=NF;i++) s=s OFS $i; print s; next }
        { v=$1":"$2":"$4":"$5; s=v OFS $7 OFS $9; for(i=10;i<=NF;i++) s=s OFS $i; print s }
      ' > "$OUTDIR/$FAM.gt.tsv"

  # File 1: all carried variants, annotation + genotypes
  awk -F'\t' -v OFS='\t' -v c="$IDCOL" '
    NR==FNR { if(FNR==1){h=$0;sub(/^[^\t]*\t/,"",h);gthdr=h;next} key=$1;rest=$0;sub(/^[^\t]*\t/,"",rest);gt[key]=rest;next }
    FNR==1 { print $0, gthdr; next }
    ($c in gt) { print $0, gt[$c] }
  ' "$OUTDIR/$FAM.gt.tsv" <(reader) > "$allf"

  # File 2: keep only rows called by two or three callers (derived from File 1, no rescan)
  local fcol; fcol=$(head -1 "$allf" | tr '\t' '\n' | grep -nxF ensemble_FILTER | cut -d: -f1)
  awk -F'\t' -v c="$fcol" 'NR==1 || $c ~ /twoCallers|threeCallers/' "$allf" > "$twof"

  rm -f "$OUTDIR/$FAM.gt.tsv"
  echo "   all: $(($(wc -l < "$allf")-1))   >=2 callers: $(($(wc -l < "$twof")-1))"
}

mapfile -t FAMS < <(ls "$FAM_DIR"/*.lst | sort)
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  i=$SLURM_ARRAY_TASK_ID
  [ "$i" -lt "${#FAMS[@]}" ] || { echo "index $i >= ${#FAMS[@]}" >&2; exit 1; }
  process_family "${FAMS[$i]}"
else
  for f in "${FAMS[@]}"; do process_family "$f"; done
fi
