#########################################
# expect snakemake rule
# input: 
#    intervar=
#    annovar=
#    snpeff=
# output:
#    txt=
#    slim=


# Create a rscript to call P/LP from the annotation files
# rm(list=ls(all=T))
# sessionInfo()

### By default
# snpeff and annovar keep the original chr ID
# but intervar removes chr prefix
# therefore, the merge process is fine for Genomel and gnomad 
# but failed for ukb data
# I need revise the R code to drop chr prefix if exists.

# library(pacman)
library(tidyverse) #1.3.1
library(bigreadr)  #0.2.4
library(rlang) #1.0.5

intervar_fn <-  snakemake@input[["intervar"]]
annovar_fn <- snakemake@input[["annovar"]]
snpeff_fn <- snakemake@input[["snpeff"]]

# intervar_fn <- "/Volumes/data/ukbb/dnanexus/ukb_batch2/intervar/00000.hg38_multianno.txt.intervar"
# annovar_fn <- "/Volumes/data/ukbb/dnanexus/ukb_batch2/annovar/00000.hg38_multianno.txt"
# snpeff_fn <- "/Volumes/data/ukbb/dnanexus/ukb_batch2/snpeff/00000.vcf"



add_vid <- function(df, new_vid="variant_id", chr_pos_ref_alt_columns=c("Chr", "Start", "Ref", "Alt") ){
  var_enq <- rlang::sym(new_vid)
  rv <- df %>% 
    unite(!!var_enq, all_of(chr_pos_ref_alt_columns), sep=':', remove=F) %>% 
    mutate(!!var_enq := sub("^chr", "", !!var_enq))
  
  return(rv)
}

# tidyverse bigreadr

anno_impact <- bigreadr::fread2(sprintf("grep -v '^##' %s", snpeff_fn), sep="\t", fill=T, check.names=T, select=c("ID", "INFO")) %>% tidyr::extract(INFO, into=c("ANN"), regex=";ANN=(.+)$") %>% 
  separate(ANN, into=c(NA, NA, "impact"), sep='\\|', remove=F, extra="drop") %>%
  mutate(ID=sub("^chr", "", ID)) %>%
  dplyr::select(vid=ID, impact)

anno_intervar <- bigreadr::fread2(intervar_fn) %>%
  dplyr::rename(Chr='#Chr') %>%
  add_vid(new_vid ="vid") %>% 
  mutate(intervar=str_match(`InterVar: InterVar and Evidence`, "InterVar: (.+) PVS1") %>% .[,2])

anno_annovar <- bigreadr::fread2(annovar_fn, sep="\t", fill=T, check.names=T, na.strings='.') %>%  
  rename_at(vars(ends_with(".1")), ~sub("\\.1$", ".genome", .)) %>% 
  add_vid(new_vid ="vid")

recode_freq <- function(v){
  v[v=="."] <- "-1"
  return(as.numeric(v))
}

psum <- function(..., na.rm=FALSE) {
  x <- list(...)
  rowSums(matrix(unlist(x), ncol=length(x)), na.rm=na.rm)
}


# source("~/Rstat/my.util.r")
# comm <- compare(anno_impact$vid, anno_intervar$vid)
# comm <- compare(anno_impact$vid, anno_annovar$vid)

# dim(anno_annovar)
# dim(anno_intervar) #duplicate 
# dim(anno_impact) # 25656

# 
# anno_intervar_cnt <- anno_intervar%>% count(vid, intervar)
# anno_intervar_cnt %>% head
# anno_intervar %>% filter(vid == "1:29396:C:T" )
# vid Chr Start   End Ref Alt   Ref.Gene Func.refGene ExonicFunc.refGene Gene.ensGene avsnp147 AAChange.ensGene
# 1 1:29396:C:T   1 29396 29396   C   T MIR1302-10     upstream                  .  MIR1302-2HG        .                .
# 2 1:29396:C:T   1 29396 29396   C   T MIR1302-11     upstream                  .  MIR1302-2HG        .                .
# 3 1:29396:C:T   1 29396 29396   C   T  MIR1302-2     upstream                  .  MIR1302-2HG        .                .
# 4 1:29396:C:T   1 29396 29396   C   T  MIR1302-9     upstream                  .  MIR1302-2HG        .                .
# 5 1:29396:C:T   1 29396 29396   C   T     WASH7P     upstream                  .  MIR1302-2HG        .                .
# anno_annovar %>% filter(vid == "1:29396:C:T" )
# vid Chr Start   End Ref Alt Func.refGene                                     Gene.refGene GeneDetail.refGene
# 1 1:29396:C:T   1 29396 29396   C   T     upstream MIR1302-10;MIR1302-11;MIR1302-2;MIR1302-9;WASH7P           dist=970
# ExonicFunc.refGene AAChange.refGene Func.knownGene Gene.knownGene GeneDetail.knownGene ExonicFunc.knownGene AAChange.knownGene
# 1               <NA>             <NA> ncRNA_intronic         WASH7P                 <NA>                 <NA>               <NA>
#   Func.ensGene Gene.ensGene GeneDetail.ensGene ExonicFunc.ensGene AAChange.ensGene ExAC_nontcga_ALL ExAC_nontcga_AFR ExAC_nontcga_AMR
# 1     upstream  MIR1302-2HG  


# add impact from snpEff and InterVar from anno_intervar
# freq_cols <- c("Freq_gnomAD_genome_ALL",	"Freq_esp6500siv2_all",	"Freq_1000g2015aug_all")
freq_cols_annovar <- c("AF_popmax", "non_cancer_AF_popmax")

# anno_merged =  anno_annovar %>% 
#                left_join(anno_impact, by="vid") %>%
#                left_join(anno_intervar %>% filter(!duplicated(vid)) %>% select(vid, intervar), by="vid") %>%
#                # define popmax_freq
#                mutate(across(all_of(freq_cols_annovar), ~recode_freq(.x))) %>% 
#                rowwise() %>% 
#                mutate(popmax_freq = max(c_across(all_of(freq_cols_annovar)))) %>%
#                ungroup

# anno_merged =  anno_annovar %>% 
#   left_join(anno_impact, by="vid") %>%
#   left_join(anno_intervar %>% filter(!duplicated(vid)) %>% select(vid, intervar), by="vid") %>%
#   # define popmax_freq
#   mutate(across(all_of(freq_cols_annovar), ~recode_freq(.x))) %>% 
#   mutate(popmax_freq = pmax(AF_popmax, non_cancer_AF_popmax)  ) 


anno_merged =  anno_annovar %>% 
  left_join(anno_impact, by="vid") %>%
  left_join(anno_intervar %>% filter(!duplicated(vid)) %>% select(vid, intervar), by="vid") %>%
  # define popmax_freq
  mutate(across(all_of(freq_cols_annovar), ~recode_freq(.x))) %>% 
  mutate(popmax_freq = pmax(!!! syms(freq_cols_annovar )) %>% coalesce(0))



spliceai_cols <- c("DS_AG", "DS_AL", "DS_DG", "DS_DL")
ACCEPTED_CLNREVSTAT <- c("criteria_provided,_multiple_submitters,_no_conflicts", "practice_guideline", "reviewed_by_expert_panel")

plp <- anno_merged %>% 
  tidyr::extract(spliceai_filtered, into=spliceai_cols, regex="DS_AG=([^;]+);DS_AL=([^;]+);DS_DG=([^;]+);DS_DL=([^;]+);", remove=F, convert=T) %>%
  mutate(missense= psum( MetaSVM_pred=='D', as.numeric(REVEL_score) > 0.75, CADD_phred > 30, BayesDel_addAF_pred == 'D', na.rm=T) %>% coalesce(0) ) %>%
  mutate(splice= psum(dbscSNV_RF_SCORE > 0.6, abs(dpsi_zscore) > 2, pmax(!!! syms(spliceai_cols), na.rm=T) > 0.5, na.rm=T ) %>% coalesce(0) ) %>%
  mutate(
    PLP.clinvar = grepl("pathogenic", CLNSIG, ignore.case = T) & CLNREVSTAT %in%  ACCEPTED_CLNREVSTAT,
    PLP.intervar = grepl("pathogenic", intervar, ignore.case = T),
    PLP.impactHigh = (impact == "HIGH") %>% coalesce(F), 
    PLP.genomel = ((impact == 'HIGH') | (impact == "MODERATE" & Polyphen2_HDIV_pred == 'D')) %>% coalesce(F),
    PLP.m_sp = missense >= 3 | splice>=2 | missense + splice>=3,
    PLP.jung = PLP.clinvar | PLP.intervar | PLP.impactHigh | (intervar=="Uncertain significance" & PLP.m_sp)
  )

# QC
# head(plp %>% filter(PLP.m_sp) %>% select(MetaSVM_pred, REVEL_score, CADD_phred, BayesDel_addAF_pred, dbscSNV_RF_SCORE, dpsi_zscore, all_of(spliceai_cols) , missense, splice, PLP.m_sp ))

write_delim(plp, file=snakemake@output[["txt"]], quote="none", delim="\t")

plp %>% select(vid, Gene.refGene, popmax_freq, starts_with("PLP")) %>%
  write_delim( file=snakemake@output[["slim"]], quote="none", delim="\t")
