### I cannot use module load snakemake as /gpfs/gsfs12/users/zhuw10/ukbb/dnanexus/var_anno_snakemake/workflow/profiles/biowulf/slurm-submit.py cannot locate snakemake moduel properly
# slurm-submit.py uses /usr/bin/env pathon, which is different from the one used with snakemake module.  It might be the reason to have the error.
conda activate snakemake
snakemake --profile profiles/biowulf --jobs 100 --use-envmodules -R prep

snakemake --profile profiles/biowulf --jobs 100 --use-envmodules --latency-wait 60

squeue -u zhuw10

snakemake  --jobs 4 --use-envmodules --latency-wait 60

zcat /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2.vcf.gz | grep -v "^#"  > foo
mkdir /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2_splitted
split --numeric-suffixes=0 -n l/10 --suffix-length=5  --additional-suffix=".txt" foo  "$(dirname /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2_splitted/00000.txt)/"

ls -al /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2_splitted

tmpfile=$(mktemp /tmp/abc-script.XXXXXX); zcat /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2.vcf.gz | grep -v "^#"  > $tmpfile; split --numeric-suffixes=0 -n l/10 --suffix-length=5  --additional-suffix=".txt" $tmpfile  "$(dirname /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2_splitted/00000.txt)/"


tmpfile=$(mktemp /tmp/abc-script.XXXXXX); /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2.vcf.gz | grep -v "^#"  > $tmpfile; split --numeric-suffixes=0 -n l/10 --suffix-length=5  --additional-suffix=".txt" $tmpfile   "$(dirname /data/zhuw10/ukbb/dnanexus/batch2/ukb_batch2_splitted/00000.txt)/"

snakemake --list-code-changes

### rerun prep
snakemake -c4


module load bcftools
 bcftools norm -m-both -Ou --threads 4 /data/zhuw10/ukbb/dnanexus/pvcf2/UKBB_500K.BAP1.0.vcf.gz | bcftools norm -f /data/zhuw10/ukbb/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa | bcftools annotate -Oz -x ID  -I +"%CHROM:%POS:%REF:%ALT" --threads 4 | zcat | cut -f 1-8 > /data/zhuw10/ukbb/dnanexus/batch2/prep/UKBB_500K.BAP1.0.vcf

mkdir anno

# -nastring must be '.' when '-vcfinput' is specified
# output is: anno/genomel_28genes.hg38_multianno.vcf
# genomel_28genes.avinput
# genomel_28genes.hg38_multianno.txt

perl $ANNOVAR_HOME/table_annovar.pl --vcfinput ../../batch2/ukb_batch2_splitted/00009.vcf $ANNOVAR_DATA/hg38 --buildver hg38 --out anno/genomel_28genes --remove --protocol refGene,knownGene,ensGene,exac03nontcga,gnomad211_exome,gnomad211_genome,esp6500siv2_all,1000g2015aug_all,clinvar_20210123,dbnsfp41a,dbscsnv11,spliceai_filtered,spidex,cosmic92_coding --operation g,g,g,f,f,f,f,f,f,f,f,f,f,f --nastring '.'

grep -c -v "^#" ../../batch2/ukb_batch2_splitted/00009.vcf
1190

wc -l anno/genomel_28genes.avinput
1190 anno/genomel_28genes.avinput

wc -l anno/genomel_28genes.avinput

cut -f 11 anno/genomel_28genes.avinput | sort > anno/f1
cut -f 3 ../../batch2/ukb_batch2_splitted/00009.vcf | sort -u > anno/f2

chr16   67660036        chr16:67660036:C:T      C       T       48      .       AF=0.000675;AQ=48;ANNOVAR_DATE=2020-06-08;Func.
refGene=exonic;Gene.refGene=ACD;GeneDetail.refGene=.;ExonicFunc.refGene=nonsynonymous_SNV;AAChange.refGene=ACD:NM_001082486:exo
n2:c.G109A:p.D37N,ACD:NM_022914:exon2:c.G100A:p.D34N;Func.knownGene=exonic;Gene.knownGene=ACD;GeneDetail.knownGene=.;ExonicFunc
.knownGene=nonsynonymous_SNV;AAChange.knownGene=ACD:ENST00000219251.13:exon2:c.G100A:p.D34N,ACD:ENST00000602320.1:exon2:c.G100A
:p.D34N,ACD:ENST00000620338.4:exon2:c.G367A:p.D123N,ACD:ENST00000620761.6:exon2:c.G109A:p.D37N;Func.ensGene=exonic;Gene.ensGene
=ACD;GeneDetail.ensGene=.;ExonicFunc.ensGene=nonsynonymous_SNV;AAChange.ensGene=ACD:ENST00000219251.13:exon2:c.G100A:p.D34N,ACD
:ENST00000602320.1:exon2:c.G100A:p.D34N,ACD:ENST00000620338.4:exon2:c.G367A:p.D123N,ACD:ENST00000620761.6:exon2:c.G109A:p.D37N;
ExAC_nontcga_ALL=0.0004;ExAC_nontcga_AFR=0;ExAC_nontcga_AMR=0.0006;ExAC_nontcga_EAS=0;ExAC_nontcga_FIN=0;ExAC_nontcga_NFE=0.000
6;ExAC_nontcga_OTH=0;ExAC_nontcga_SAS=0;AF=0.0003;AF_popmax=0.0005;AF_male=0.0004;AF_female=0.0003;AF_raw=0.0003;AF_afr=0;AF_sa
s=0;AF_amr=0.0005;AF_eas=0;AF_nfe=0.0005;AF_fin=5.053e-05;AF_asj=0;AF_oth=0.0010;non_topmed_AF_popmax=0.0005;non_neuro_AF_popma
x=0.0005;non_cancer_AF_popmax=0.0006;controls_AF_popmax=0.0006;AF=0.0005;AF_popmax=0.0008;AF_male=0.0005;AF_female=0.0005;AF_ra
w=0.0005;AF_afr=0.0001;AF_sas=.;AF_amr=0;AF_eas=0;AF_nfe=0.0008;AF_fin=0;AF_asj=0;AF_oth=0.0009;non_topmed_AF_popmax=0.0009;non
_neuro_AF_popmax=0.0007;non_cancer_AF_popmax=.;controls_AF_popmax=0.0004;esp6500siv2_all=0.0005;1000g2015aug_all=.;CLNALLELEID=466553;CLNDN=Dyskeratosis_congenita,_autosomal_dominant_6;CLNDISDB=MONDO:MONDO:0014690,MedGen:C4225284,OMIM:616553;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;DamagePredCount=4.17;SIFT_pred=.;SIFT4G_pred=D;Polyphen2_HDIV_pred=P;Polyphen2_HVAR_pred=P;LRT_pred=N;MutationTaster_pred=N;MutationAssessor_pred=N;FATHMM_pred=.;PROVEAN_pred=.;VEST4_score=0.356;MetaSVM_pred=T;MetaLR_pred=T;M-CAP_pred=D;REVEL_score=0.081;MutPred_score=.;MVP_score=0.536;MPC_score=0.312;PrimateAI_pred=T;DEOGEN2_pred=T;BayesDel_addAF_pred=T;BayesDel_noAF_pred=T;ClinPred_pred=T;LIST-S2_pred=T;CADD_raw=3.477;CADD_phred=24.600;DANN_score=0.991;fathmm-MKL_coding_pred=D;fathmm-XF_coding_pred=D;Eigen-raw_coding=-0.087;Eigen-phred_coding=2.217;Eigen-PC-raw_

### Test intervar comamnd
module load intervar/2.1.3
# module load annovar/2020-06-08 # to update $ANNOVAR_DATA/hg38
# alternative soltuion is to use $ANNOVAR_DATA_CURRENT
LC_CTYPE=C InterVar --input_type=VCF -i /data/zhuw10/ukbb/dnanexus/batch2/splitted/00006.vcf -b hg38 -d $ANNOVAR_DATA_CURRENT/hg38 -o /data/zhuw10/ukbb/dnanexus/batch2/intervar/00006 
# Error: can't read the LOF genes file intervardb/PVS1.LOF.genes.hg38
Traceback (most recent call last):
  File "/usr/local/apps/intervar/2.1.3/bin/Intervar.py", line 2114, in <module>
    main()
  File "/usr/local/apps/intervar/2.1.3/bin/Intervar.py", line 2076, in main
    read_datasets()
  File "/usr/local/apps/intervar/2.1.3/bin/Intervar.py", line 440, in read_datasets
    strs = fh.read()
  File "/usr/local/Anaconda/envs/py3.8/lib/python3.8/codecs.py", line 322, in decode
    (result, consumed) = self._buffer_decode(data, self.errors, final)
UnicodeDecodeError: 'utf-8' codec can't decode byte 0xfc in position 6538: invalid start byte

### 
https://stackoverflow.com/questions/2276200/changing-default-encoding-of-python

    439         fh = open(paras['orpha'], "r")
    440         strs = fh.read()

cd ~/git/InterVar-2.2.1/

cat <<EOF >/data/zhuw10/ukbb/dnanexus/intervar.cmd 
./Intervar.py -i /data/zhuw10/ukbb/dnanexus/anno/genomel_28genes.avinput -b hg38 -d $ANNOVAR_DATA/hg38 -o /data/zhuw10/ukbb/dnanexus/anno/genomel_28genes.intervar.txt --table_annovar=$ANNOVAR_HOME/table_annovar.pl --convert2annovar=$ANNOVAR_HOME/convert2annovar.pl --annotate_variation=$ANNOVAR_HOME/annotate_variation.pl
EOF

./Intervar.py --input_type=VCF -i /data/zhuw10/ukbb/dnanexus/batch2/splitted/00006.vcf -b hg38 -d $ANNOVAR_DATA_CURRENT/hg38 -o /data/zhuw10/ukbb/dnanexus/batch2/intervar/00006 --table_annovar=$ANNOVAR_HOME/table_annovar.pl --convert2annovar=$ANNOVAR_HOME/convert2annovar.pl --annotate_variation=$ANNOVAR_HOME/annotate_variation.pl


### failed command from biowulf module
### InterVar
#!/bin/bash
$INTERVAR_BIN/Intervar.py -t $INTERVAR_DATA/intervardb --annotate_variation=$ANNOVAR_HOME/annotate_variation.pl --table_annovar=$ANNOVAR_HOME/table_annovar.pl --convert2annovar=$ANNOVAR_HOME/convert2annovar.pl $@

Notice: Your command of InterVar is ['/usr/local/apps/intervar/2.1.3/bin/Intervar.py', '-t', '/usr/local/apps/intervar/DOWNLOADS/InterVar/intervardb', '--annotate_variation=/usr/local/apps/annovar/2018-04-16/annotate_variation.pl', '--table_annovar=/usr/local/apps/annovar/2018-04-16/table_annovar.pl', '--convert2annovar=/usr/local/apps/annovar/2018-04-16/convert2annovar.pl', '--input_type=VCF', '-i', '/data/zhuw10/ukbb/dnanexus/batch2/splitted/00006.vcf', '-b', 'hg38', '-d', '/fdb/annovar/current/hg38', '-o', '/data/zhuw10/ukbb/dnanexus/batch2/intervar/00006']

### successful comamnd 
Notice: Your command of InterVar is ['./Intervar.py', '--input_type=VCF', '-i', '/data/zhuw10/ukbb/dnanexus/batch2/splitted/00006.vcf', '-b', 'hg38', '-d', '/fdb/annovar/current/hg38', '-o', '/data/zhuw10/ukbb/dnanexus/batch2/intervar/00006', '--table_annovar=/usr/local/apps/annovar/2018-04-16/table_annovar.pl', '--convert2annovar=/usr/local/apps/annovar/2018-04-16/convert2annovar.pl', '--annotate_variation=/usr/local/apps/annovar/2018-04-16/annotate_variation.pl']
INFO: The options are {'buildver': 'hg38', 'inputfile': '/data/zhuw10/ukbb/dnanexus/batch2/splitted/00006.vcf', 'inputfile_type': 'VCF', 'outfile': '/data/zhuw10/ukbb/dnanexus/batch2/intervar/00006', 'database_intervar': 'intervardb', 'lof_genes': 'intervardb/PVS1.LOF.genes.hg38', 'pm1_domain': 'intervardb/PM1_domains_with_benigns.hg38', 'mim2gene': 'intervardb/mim2gene.txt', 'mim_recessive': 'intervardb/mim_recessive.txt', 'mim_domin': 'intervardb/mim_domin.txt', 'mim_adultonset': 'intervardb/mim_adultonset.txt', 'mim_pheno': 'intervardb/mim_pheno.txt', 'mim_orpha': 'intervardb/mim_orpha.txt', 'orpha': 'intervardb/orpha.txt.utf8', 'knowngenecanonical': 'intervardb/knownGeneCanonical.txt.hg38', 'pp2_genes': 'intervardb/PP2.genes.hg38', 'bp1_genes': 'intervardb/BP1.genes.hg38', 'ps1_aa': 'intervardb/PS1.AA.change.patho.hg38', 'ps4_snps': 'intervardb/PS4.variants.hg38', 'bs2_snps': 'intervardb/BS2_hom_het.hg38', 'exclude_snps': 'intervardb/ext.variants.hg38', 'evidence_file': 'None', 'disorder_cutoff': '0.01', 'onetranscript': 'FALSE', 'otherinfo': 'TRUE', 'convert2annovar': '/usr/local/apps/annovar/2018-04-16/convert2annovar.pl', 'table_annovar': '/usr/local/apps/annovar/2018-04-16/table_annovar.pl', 'annotate_variation': '/usr/local/apps/annovar/2018-04-16/annotate_variation.pl', 'database_locat': '/fdb/annovar/current/hg38', 'database_names': 'refGene esp6500siv2_all 1000g2015aug avsnp147 dbnsfp42a clinvar_20210501 gnomad_genome dbscsnv11 rmsk ensGene knownGene', 'current_version': 'Intervar_20210727', 'public_dev': 'https://github.com/WGLab/InterVar/releases', 'skip_annovar': False} 

module load intervar/2.2.1
Intervar.py --input_type=VCF -i /data/zhuw10/ukbb/dnanexus/batch2/splitted/00005.vcf -b hg38 -d $ANNOVAR_DATA_CURRENT/hg38 -o /data/zhuw10/ukbb/dnanexus/batch2/intervar/00005 -t /home/zhuw10/git/InterVar-2.2.1/intervardb --table_annovar=$ANNOVAR_HOME/table_annovar.pl --convert2annovar=$ANNOVAR_HOME/convert2annovar.pl --annotate_variation=$ANNOVAR_HOME/annotate_variation.pl


### run for genomel
# make a new config.yaml
conda activate snakemake
snakemake --profile profiles/biowulf --jobs 400 --use-envmodules

module load annovar/2020-06-08
/home/zhuw10/git/InterVar-2.2.1/Intervar.py --input_type=VCF -i /data/zhuw10/ukbb/dnanexus/genomel_out/splitted/00003.vcf -b hg19 -d $ANNOVAR_DATA/hg19 -o /data/zhuw10/ukbb/dnanexus/genomel_out/intervar/00003 -t /home/zhuw10/git/InterVar-2.2.1/intervardb --table_annovar=$ANNOVAR_HOME/table_annovar.pl --convert2annovar=$ANNOVAR_HOME/convert2annovar.pl --annotate_variation=$ANNOVAR_HOME/annotate_variation.pl

slurm snakemake sinfo: error: If munged is up, restart with --num-threads=10