# Snakemake workflow to predict P/LP (pathogenic/likely pathogenic variants)

## Introduction 

Assuming the input data are a list of VCF files (for examples, multiple-sample VCF files split by chromosomes or genes), this workflow runs annovar, intervar and snpEff to predict P/LP variants in two ways: 1. Clinvar+Intervar (CI); 2. Jung's criteria as descried in [this paper](https://doi.org/10.1016/j.cancergen.2020.10.002) (as I learned from Jung Kim  ). The first method is straightforward, that is, the variant is predicted by as P/LP by either InterVar or ClinVar (with the CLINSTAT tag as "criteria_provided,_multiple_submitters,_no_conflicts",  "practice_guideline",  or "reviewed_by_expert_panel").  The second method is described in the [paper](https://doi.org/10.1016/j.cancergen.2020.10.002).  I named it "Jung's criteria" as I learned it from Dr. Jung Kim (and I have no idea how to call it ).  

Assuming the input data are a list of VCF files (for examples, multiple-sample VCF files split by chromosomes or genes), this workflow runs annovar, intervar and snpEff to predict P/LP variants in two ways: 1. Clinvar+Intervar (CI); 2. Jung's criteria as descried in [this paper](https://doi.org/10.1016/j.cancergen.2020.10.002). The first method is straightforward, that is, the variant is predicted by as P/LP by either InterVar or ClinVar (with the CLINSTAT tag as "criteria_provided,_multiple_submitters,_no_conflicts", "practice_guideline", or "reviewed_by_expert_panel").  The second method is a little complicated, as described in the [paper](https://doi.org/10.1016/j.cancergen.2020.10.002).  I named it "Jung's criteria" as I learned it from Dr. Jung Kim (and I have no idea how to call it :) ).  Briefly, it is a union set of Clinvar+InterVar and two other variant sets based on the variant annotations from InterVar, AnnoVar and snpEff:
1. variants with High impact (snpEff);
2. variants meeting all the following criteria;
   + Variants are predicted as VUS (variants of uncertain significance) by InterVar;
   + Variants meets "3/4", "2/3", or "3/7"
     + Predicted as ”damaging” missense mutations by 3 out of the 4 predictions;
     + Predicted as "damaging" splice site mutations by 2 out of 3 predictions;
     + Predicted as "damaging" missense or splice site mutations by 3 out of 7 predictions;
  
The overall process is as described in the diagram below: 
![](images/pathogenic_identificatioin.png)

Users may have a look of the R script [Call_patho.R]( workflow/scripts/Call_patho.R) to understand how P/LP variants are actually predicted.


There are several steps taking in this Snakemake workflow as outlined in [this diagram](./workflows/plp_dag_expand.pdf):
+ Merge input vcf files into one file; 
  + :bookmark: Only the first 8 columns of the VCF files are used so as to reduce the file size and facilitate the subsequent data processing.
+ Split vcf files into multiple parts to speed up the annotation process.  
+ Annotate each split part by *InterVar*, *Annovar* and *snpEff*.
+ Call P/LP variants using the R script [Call_patho.R]( workflow/scripts/Call_patho.R).
  + There are two output files generated in this step: *{chunk}.plp.txt* and *{chunk}.plp_slim.txt*
    +  All the essential  output from the variant annotators are kept in the *plp.txt* file.
    +  The plp_slim.txt file has only 8 columns:
    

| vid              | Gene.refGene   | popmax_freq | PLP.clinvar | PLP.intervar | PLP.impactHigh | PLP.genomel | PLP.m_sp | PLP.jung |
| ---------------- | -------------- | ----------- | ----------- | ------------ | -------------- | ----------- | -------- | -------- |
| 4:186317998:C:A  | ANKRD37;LRP2BP | 0           | FALSE       | FALSE        | FALSE          | FALSE       | FALSE    | FALSE    |
| 4:186318030:C:T  | ANKRD37;LRP2BP | 0           | FALSE       | FALSE        | FALSE          | FALSE       | FALSE    | FALSE    |
| 4:186318071:G:A  | ANKRD37;LRP2BP | 1.00E-04    | FALSE       | FALSE        | FALSE          | FALSE       | FALSE    | FALSE    |
| 4:186318186:G:A  | LRP2BP         | 0           | FALSE       | FALSE        | FALSE          | FALSE       | FALSE    | FALSE    |
| 4:186318261:G:A  | LRP2BP         | 8.00E-04    | FALSE       | FALSE        | FALSE          | FALSE       | FALSE    | FALSE    |
| 4:186318265:A:G  | LRP2BP         | 6.84E-05    | FALSE       | FALSE        | FALSE          | FALSE       | FALSE    | FALSE    |
| 4:186318282:T:TG | LRP2BP         | 0           | FALSE       | FALSE        | FALSE          | FALSE       | FALSE    | FALSE    |


+ Below is the description of the columns in the slim file

| Column ID      | Description                                                                                                      |
| -------------- | ---------------------------------------------------------------------------------------------------------------- |
| vid            | Variant ID in the format of CHROM:POS:REF:ALT                                                                    |
| Gene.refGene   | Gene ID                                                                                                          |
| popmax_freq    | Maximum AF among the subpopulations (derived from "AF_popmax" and "non_cancer_AF_popmax")                        |
| PLP.clinvar    | P/LP defined by ClinVar                                                                                          |
| PLP.intervar   | P/LP defined by InterVar                                                                                         |
| PLP.impactHigh | P/LP defined by high impact                                                                                      |
| PLP.genomel    | P/LP defined ***(impact == 'HIGH') \| (impact == "MODERATE" & Polyphen2_HDIV_pred == 'D')***                           |
| PLP.m_sp       | P/LP defined as Damaging missense and/or splite site mutations                                                   |
| PLP.jung       | P/LP defined by ***PLP.clinvar \| PLP.intervar \| PLP.impactHigh \| (intervar=="Uncertain significance" & PLP.m_sp)*** |