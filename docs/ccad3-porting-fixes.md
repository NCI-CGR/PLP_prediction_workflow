# Getting the AutoGVP P/LP workflow running on CCAD3 (Slurm)

Branch: `fix/ccad3-slurm-portability`

This branch fixes the problems I ran into while running the pipeline on CCAD3
for the full SR0467-047 cohort (246 IBMS Families WGS samples, hg38, about 86.7
million variant sites). Each section below is one problem: what went wrong, why,
and the exact change that fixes it. 

---

## 1. The pipeline only looked at one sample instead of all 246

**What I saw.** The run produced about 5.27 million variants instead of the
roughly 86.7 million I expected. It turned out only the first sample had been
classified; the other 245 were missing.

**Why.** The prep step was supposed to throw away the per-sample genotypes and
keep just the list of variant positions. It used `bcftools annotate -x FORMAT`
to do that, but that command does not actually remove the genotypes. Because the
genotypes were still there, the next tool (InterVar's converter) quietly kept
only the variants where the very first sample had a change, and everyone else
was dropped.

**The fix (`workflow/Snakefile`, `prep` rule).** Use `bcftools view -G`, which
really does drop the genotypes, so the file becomes a clean positions-only list
that covers every sample:

```
# before
... | bcftools annotate -Oz -x FORMAT,^INFO/AC,^INFO/AF,^INFO/AN -o {output.vcf}

# after
... | bcftools annotate -Ou -x ^INFO/AC,^INFO/AF,^INFO/AN \
    | bcftools view -G -Oz --threads {threads} -o {output.vcf}
```

To check it worked, a prep output file should have 8 columns:
`zcat <prep>.vcf.gz | grep -m1 -v '^#' | awk '{print NF}'` should print `8`.

**Commit:** `fix(prep): drop genotypes with view -G so all samples are classified`

---

## 2. bcftools would not load when a conda environment was active

**What I saw.** The compute jobs failed with `bcftools: command not found`. The
log said it was trying to unload `anaconda` and gave up:
`Run conda deactivate before unloading`.

**Why.** On CCAD3, loading the `bcftools` module automatically tries to unload
`anaconda`, and it refuses to do that while a conda environment is switched on.
Our pipeline runs inside the `AutoGVP` conda environment, and Slurm passes that
environment down to every job, so the module could never load. I also cannot
just unload anaconda ourselves, because that is where `conda` itself comes from,
so removing it breaks `conda activate`.

**The fix (`environment.yaml`).** Stop relying on the module and get bcftools
straight from the conda environment instead (plus htslib, which provides `bgzip`
and `tabix`). Once bcftools is in the environment, it does not matter that the
module fails to load:

```
dependencies:
  - bioconda::bcftools>=1.19
  - bioconda::htslib
```

(I did this live with `conda install -n AutoGVP -c bioconda bcftools htslib`.)

Small scripts that need bcftools but not the pipeline can still use
`module load bcftools`, as long as no conda environment is turned on in that
shell.

**Commit:** `fix(env): provide bcftools/htslib via conda to avoid module/conda conflict on CCAD3`

---

## 3. The final merge step ran out of memory

**What I saw.** The `merge_call` step was killed for using too much memory (more
than 180 GB) on the full cohort. It had worked fine before, but only because
back then it was merging one sample's worth of data.

**Why.** It used `csvtk sort`, which loads the whole table into memory to sort
it. That is fine for a few million rows but not for ~86.7 million.

**The fix (`workflow/Snakefile`, `merge_call` rule).** Swap `csvtk sort` for the
regular Unix `sort`, which sorts using disk instead of holding everything in
memory, so it stays within a small memory budget. The rest of the step stays the
same. Chromosome is column 1 and position is column 2 in these files, and the
little `awk` trick keeps the header row on top while sorting:

```
# before
| csvtk sort -t -k chr:N -k start:n -E \

# after
| awk -F'\t' 'NR==1{print "0\t"$0; next}{print "1\t"$0}' \
| LC_ALL=C sort -t$'\t' -k1,1n -k2,2V -k3,3n -S 4G -T "$TMPDIR" --parallel={threads} \
| cut -f2- \
```

After this change, the merge no longer needs a big-memory node. There is also a
standalone version of this in `scripts/finish_merge_call.sh` if you ever need to
redo just the merge without rerunning the whole pipeline.

**Commit:** `fix(merge_call): use disk-backed external sort to avoid OOM at cohort scale`

---

## 4. Some jobs needed more time to finish

**What I saw.** Now that the pipeline processes the whole cohort (the broken run
only ever reached these steps with one sample's data), several steps were at
risk of hitting their time limit and being killed.

**The fix (`config/cluster_config.yaml`).** Give the slower steps more wall-clock
time:

| step          | before      | after       |
|---------------|-------------|-------------|
| vep_ann       | 5 days      | 7 days      |
| intervar_ann  | 1 day       | 3 days      |
| autopvs1      | 12 hours    | 2 days      |
| annovar_ann   | 1 day       | 2 days      |
| autogvp       | 5 days      | 7 days      |
| merge_call    | 2 days      | 4 days      |

Check these against what each partition actually allows first:
`sinfo -o "%P %l"`. Thanks to fix #3, `merge_call` can now run on the normal
`defq` partition instead of the big-memory one.

**Commit:** `chore(cluster): raise wall-time limits for the full-cohort run`

---

## 5. Helper scripts

Small standalone scripts added under `scripts/`. They run as Slurm array jobs and
get bcftools from the conda environment:

- `finish_merge_call.sh` - rebuild the merged tables from the 100 chunk files
  using the memory-safe sort from fix #3, without rerunning the pipeline.
- `split_autogvp_by_family.sh` - split the results into one file per family,
  using the genotype VCF to figure out which variants each family carries. Each
  row includes every family member's full genotype, and it writes an
  all-variants file plus a version limited to variants called by at least two
  callers.
- `filter_plp_by_family.sh` - pull just the Pathogenic / Likely pathogenic rows
  for each family, giving small files that are easy to review.

**Commit:** `feat(scripts): add per-family split and P/LP extraction utilities`

