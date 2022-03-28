### Usage: Run_gatk_preprocessing.sh <sorted bam> <reference>
###: Author Madikay Senghore
###: Synposis: given a bam file. and a reference file, run the gatk snp calling pipeline
###: Date Wed 23 Mar 2022

### Load modules
module load bcftools/1.5-fasrc02
module load samtools/1.5-fasrc02
module load Java
module load gatk

i=$1
IFS='.' read -r -a f <<< "$i"

### ### Drop reads under 30 bases
#samtools view -h $i | awk 'length($10) > 30 || $1 ~ /^@/'  | samtools view -bS > $f.sf.bam

## classical mpileup and filter for know variant sites
bcftools view -M2 -v snps -i '%QUAL>=20' $f.pile.vcf > $f.confirmed_snps.vcf
gatk IndexFeatureFile -F $f.confirmed_snps.vcf

## Use gatk to add read groups to bam files
gatk AddOrReplaceReadGroups --INPUT $f.sf.bam --OUTPUT $f.s.f.rg.bam --RGLB library1 --RGPL illumina --RGPU unit1 --RGSM $f

### Add one liner to remove unmapped reads and reads less that 30 quality
samtools view -F 0x04 -q 30 -b $f.s.f.rg.bam > $f.s.f.rg.q30.bam

### Mark duplicates in the sorted bam files with read groups assigned
gatk --java-options "-Xmx16g" MarkDuplicates -I $f.s.f.rg.q30.bam -O $f.marked.bam -M $f.metrics

### Perform Base recalibration
gatk --java-options "-Xmx16g" BaseRecalibrator -I $f.marked.bam -R $2 --known-sites $f.confirmed_snps.vcf -O $f.recal.table

### Apply Base recalibration
gatk --java-options "-Xmx16g" ApplyBQSR -R $2 -I $f.marked.bam --bqsr-recal-file $f.recal.table -O $f.analysis.ready.bam

#In parallel run Haplotype caller
gatk --java-options "-Xmx16g" HaplotypeCaller -R $2 -I $f.analysis.ready.bam -O $f.g.vcf -ERC GVCF

# Genotype the vcf from haplottype caller
gatk --java-options "-Xmx4g" GenotypeGVCFs -R $2 -V $f.g.vcf -O $f.genotyped.vcf


### Use bcftools to filter biallelic SNP sites  where ALT allele is majority  (Threshold 80%) 
bcftools filter -i 'AD[0]/DP < 0.2  && QUAL>100 && DP>20' $f.genotyped.vcf | bcftools view -m2 -M2 -v snps -Ob -o $f.final.snps.bcf.gz
bcftools index $f.final.snps.bcf.gz
cat $2 | bcftools consensus $f.final.snps.bcf.gz > ../consensus/$f.consensus.fasta

#### Rename fasta file (Reference ID will vary)
sed -i "s/>CP000730\.1/>$f/g" ../consensus/$f.consensus.fasta; done
