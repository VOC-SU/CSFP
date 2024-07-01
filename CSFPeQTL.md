### 1 Identification of Single Nucleotide Polymorphisms (SNPs)

The genome reference is Prunus_persica_NCBIv2 and can be obtained from the source (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000346465.2)

##### 1.1 The whole genome DNA sequences (WGS) were quality trimmed using fastp

```shell
fastp -i WGS_read_1.gz -I WGS_read_2.gz -o WGS_read_qc_1.gz -O WGS_read_qc_2.gz
```

##### 1.2 The quality summary of reads was reported using FastqCount_v0.5

```shell
FastqCount_v0.5 WGS_read_1.gz > WGS_read_report.txt
```

##### 1.3 The genome index was built using BWA for DNA alignment

```shell
bwa index ref.fa ref.fa -a bwtsw
```

##### 1.4 DNA Reads were aligned to the reference using BWA,  samtools and picard

```shell
destdir=./bam
for fname in *_1.fastq.gz
do
base=${fname%_1.fastq.gz}
sm=${fname%%_*}
bwa mem -t 18 -M -R "@RG\tID:${sm}\tSM:${sm}\tPL:illumina\tLB:${sm}" ref.fa "$fname" "${base}_2.fastq.gz" | gzip -3 > "$destdir/${base}.sam.gz"
samtools fixmate -O bam $destdir/${base}.sam.gz $destdir/${base}_fixmate.bam
rm $destdir/${base}.sam.gz
samtools sort -@ 8 -O bam -o $destdir/${sm}_sorted.bam -T $destdir/${sm}_temp  $destdir/${base}_fixmate.bam
rm $destdir/${base}_fixmate.bam
java -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT="$destdir/${sm}_sorted.bam" OUTPUT="$destdir/${sm}_dedup.bam" METRICS_FILE="$destdir/${sm}_metrics.txt"
samtools index $destdir/${sm}_dedup.bam
rm $destdir/${sm}_sorted.bam
done
```

##### 1.5 Call SNP using samtools and bcftools

The name of 184 accessions were sorted in filename.txt 

The BAM files generated in step 1.5 were in the same directory.

```shell
samtools mpileup -q 20 -Q 15 -ugf ref.fa -b filename.txt|bcftools call -vmO z -o raw.vcf.gz
```

#### 1.6 SNP filitering using beagle and vcftools

```shell
vcftools --gzvcf raw.vcf.gz --remove-indels --recode --recode-INFO-all --out raw_01.vcf
awk '$5 !~ /([[:alpha:]])+,[[:alpha:]]/{print}' raw_01.vcf > raw_02.vcf
vcftools --vcf raw_02.vcf --max-missing 0.8 --minDP 3 --maf 0.05 --recode --recode-INFO-all --out raw_03.vcf
java -jar ~/bin/beagle.r1399.jar gt=raw_03.vcf out=qc.vcf nthreads=5
```

### 2 Genome-wide association studies using emmax

https://genome.sph.umich.edu/wiki/EMMAX

##### 2.1 Association analysis were conducted using EMMAX

```shell
plink --vcf Prunus_persica.raw.filter.snp.vcf --recode12 --output-missing-genotype 0 --transpose --out 184accessions_RC --allow-extra-chr &
emmax-kin-intel64 184accessions_RC -v -d 10 -o 184accessions_RC_kinship.BN.kinf & 
plink --vcf Prunus_persica.raw.filter.snp.vcf --make-bed --out Prunus_persica_cutscaffold --allow-extra-chr &
gcta64 --bfile Prunus_persica_cutscaffold --make-grm --autosome --out 184accessions_PCA &
gcta64 --grm 184accessions_PCA --pca --out 184accessions_PCA_2 &
cut -d " " -f 2,3,4,5 184accessions_PCA_2.eigenvec >184accessions_PCA_3.txt
awk '{print $1"\t"$1"\t1\t"$2"\t"$3"\t"$4}' 184accessions_PCA_3.txt > 1 && mv 1 PCA_final.txt
ls phenotype*.txt|while read id;do emmax-intel64 -v -d 10 -t Prunus_persica_cutscaffold -p $id -k 184accessions_RC_kinship.BN.kinf -c PCA_final.txt -o ${id}_GWAS;done 
# this step will generated GWAS result for each phenotype
```

##### 2.2 Calculating the best P value

```shell
vcftools --vcf Prunus_persica.raw.filter.snp.vcf --plink --out output
plink --noweb --file output --dog --transpose --recode12 --out output
plink --noweb --file output --dog --make-bed --out output
java -Xmx4g -jar /home/caoxm/soft/gec/gec.jar --no-web --effect-number --plink-binary output --genome --out P_Value 
# The best pvalue of this study is 2.73E-6
```

##### 2.3 Manhandun Plot for GWAS

The result of emmax were used for Manhandun Plot using R package CMplot

```shell
library(CMplot)
data <-read.table('phenotype_GWAS.ps')
head(data)
CMplot(data, plot.type='m', LOG10=TRUE, ylab = "-log10(p-value)", col=c("#266135", "#6a4b98"), threshold =c(2.73E-6), threshold.col=c("red"), threshold.lwd=c(1,1), threshold.lty=c(2,1), band = 1, cex = 0.3, signal.cex = 0.3, file.output=TRUE, file='jpg', dpi=600, height=NULL, width=NULL,  verbose=TRUE)
```

### 3 eQTL analysis

The gene expression values were formatted in a matrix named gene_expression.txt.

#### 3.1 Identification of genes expressed in more than 5% accession

```shell
perl ./gene_expression_filter.pl gene_expression.txt > QC_gene_expression.txt 
```

##### supplementary script 1: gene_expression_filter.pl

```shell
open (IN,"$ARGV[0]"); # read fpkm.txt
$head= <IN>; 
print $head;
$sample_num=0;
$sample_num_expressed=0;
while (<IN>) {
chomp;
@a=split("\t",$_);
$id=shift@a;
foreach $element(@a){
    if ($element > 1) { $sample_num_expressed++;}
        $sample_num ++;
}
    $raio= $sample_num_expressed/$sample_num;
    if ($raio > 0.05) { print $_."\n";
    }

$sample_num=0;
$sample_expressed=0;
```

#### 3.2 Normalization of gene expression values

```
Rscrpit ./qq_normal.r gene_expression.txt > QC_gene_expression.txt 
```

##### supplementary script 2: qq_normal.r

```r
setwd("./")
df=read.table("QC_gene_expression.txt ",header = T,sep = "\t",check.names = F,row.names = 1)
df_2=df
for (i in 1:nrow(df)) {
  df_2[i,] = qqnorm(df[i,])$x
}
write.table(df_2,"QC_gene_expression_normal.txt",quote = F,sep = "\t")
```

#### 3.3 eQTL using emmax

The gene expression of each gene was normalized, and the import format was changed to EMMAX.

```perl
perl ./emmax_import.pl QC_gene_expression_normal.txt
```

```shell
open (IN,"$ARGV[0]");
$head=<IN>;
chomp $head;
@head=split("\t",$head);
while ($line=<IN>) {
chomp $line;
@b=split("\t",$line);
$id=shift@b;
$i=0;
open (OUT,">${id}_time.txt");
foreach $element (@b) {
print OUT $head[$i]."\t";
print OUT $head[$i]."\t";
print OUT $element."\n";
$i++
}
}
```

##### Association analysis were similar to GWAS

#### 3.4 Calculating the best Pvalue

Calculating the best Pvalue was the same as in section 2.2. The best pvalue of this study is 2.71E-6

#### 3.5 Manhandun Plot for eQTL

Manhandun Plot for eQTL was the same as in section 2.3

#### 3.6 Identification of feature SNP

Identification of feature eSNP 

```shell
ls Prupe*ps|while read id; do awk '$4<2.71E-6{print $0}' $id > /data/eQTL/${id}_qc.txt;done &
ls *qc.txt|while read id; do cut -f 1 $id >> id.txt;done
cat id.txt|sort|uniq > 1 && mv 1 eQTL_Total_SNPs.txt
ls Prupe*.txt|while read id ;do perl qc.txt_to_bed.pl $id;done
find . -name "*" -type f -size 0c | xargs -n 1 rm -f
ls Prupe*bed|while read id; do /home/bedtools2/bin/mergeBed -i $id -d 20000 -c 1,4 -o count,collapse -delim ";" > ${id}_merge.txt;done
ls *merge.txt|while read id; do perl -ne 'print if /;.*?;.*?/' $id > 1 && mv 1 $id;done
find . -name "*" -type f -size 0c | xargs -n 1 rm -f
ls *merge.txt|while read id; do cut -f 5 $id > 1 && mv 1 $id;done
perl -pi -e 's/;/\n/g' *merge.txt
perl -pi -e 's/,/\t/g' *merge.txt
```

##### supplementary script 3: qc.txt_to_bed.pl

```perl
open (IN,"$ARGV[0]");
open (OUT,">$ARGV[0].bed");
while (<IN>) {
chomp;
@a=split("\t",$_);
$id=$a[0];
@loc_info=split(":",$id);
$chr=$loc_info[0];
$loc=$loc_info[1];
$pvalue=$a[3];
print OUT $chr."\t".$loc."\t".$loc."\t".$id.",".$pvalue."\n";
}
```

##### LD clump

```shell
ls *merge.txt|while read id; do cut -f 1 $id;done |sort|uniq > Uniq_SNPs.txt
vcftools --vcf Prunus_persica.raw.filter.snp.vcf --plink --out file
plink --noweb --file file --make-bed --out file 
plink --bfile file --extract Uniq_SNPs.txt --make-bed --out eQTL_snp
plink --bfile eQTL_snp --ld-window-kb 1000 --ld-window-r2 0.3 --ld-window 99999 --r2
awk 'BEGIN{OFS="\t"}{$1=$1;print}' plink.ld > 1 && mv 1 plink.ld
cut -f 3,6,7 plink.ld > snp-snp-r2.txt
perl eQTL_LD.pl 
awk '{print FILENAME"\t"$0}' *rm_dup.txt|perl -pi -e 's/.txt.ps_qc.txt.bed_merge.txt_rm_dup.txt//' > eQTL.txt
```

##### supplementary script 4: eQTL_LD.pl

```shell
open (IN,"snp-snp-r2.txt");
while (<IN>) {
chomp;
@a=split("\t",$_);
@arrary_temp= sort($a[0],$a[1]);
$pair= join("-",@arrary_temp);
#print $pair."\n";
$hash_pair{$pair}=1; 
}

my $DIR_PATH="./";
opendir DIR, ${DIR_PATH} or die "Can not open";
@filelist = readdir DIR;
foreach $file (@filelist) {
if ($file !~ /Prupe/) {next;}  
open (IN2,"$file"); 
open (OUT,">${file}_rm_dup.txt");
$first= <IN2>;
chomp $first;
@b=split("\t",$first);
$snp_temp= $b[0];
$snp_temp_pvalue= $b[1];
while ($line= <IN2>) {
chomp $line;
@c=split("\t",$line);
$snp= $c[0];
$snp_pvalue= $c[1];
@arrary_temp_2= sort($snp,$snp_temp);
$str= join("-",@arrary_temp_2);
if (exists $hash_pair{$str}  ) { 
if ($snp_pvalue  < $snp_temp_pvalue){
$snp_temp = $snp;
$snp_temp_pvalue= $snp_pvalue;}
}
else { 
print OUT $snp_temp."\t".$snp_temp_pvalue."\n";
$snp_temp = $snp;
$snp_temp_pvalue= $snp_pvalue;
}
}
print OUT $snp_temp."\t".$snp_temp_pvalue."\n";
}#
```

##### 3.7 LD block show of eQTL

```shell
LDBlockShow -InVCF ./ Prunus_persica.raw.filter.snp.vcf -OutPut Red_ld -Region 5:-668000-714000 -InGWAS ./ps_best/Red_GWAS.ps -InGFF ./Ppersica_298_v2.1.gene.gff3
```

#### 4 GRN construction

LD clump was the same as in section 3.6.

```shell
cut -f 1,2 snp-snp-r2.txt > snp-snp.txt
perl GRN.pl 
```

##### supplementary script 5: GRN.pl

```shell
open (IN2,"Gene_SNP_Pvalue.txt"); 
# Prupe.1G255500 Pp01:26557608    3.983784266e-07

while (<IN2>){
chomp;
@eQTL=split("\t",$_);
$snp= $eQTL[1];
$gene=$eQTL[0];
if (exists $hash{$snp}) {$hash{$snp}= $hash{$snp}.",".$gene;
}
else { $hash{$snp}= $gene;   
}
}
$temp= "SNP_A SNP_B";
chomp $temp;
@temp_arrary=split("\t",$temp);
open (IN,"snp-snp.txt"); 
while ($read=<IN>) { #
chomp $read;
@read_arrary=split("\t",$read);
foreach $case(@temp_arrary){
    if($case =~ m/$read_arrary[0]/)
    {  
$flag_1=1;} }
foreach $case(@temp_arrary){
    if($case =~ m/$read_arrary[1]/)
    {  
$flag_2=1;}  }
if ($flag_1 > 0 && $flag_2 ==0 ) {
push (@temp_arrary,$read_arrary[1]);
}
elsif ($flag_2 >0  && $flag_1 ==0){
push (@temp_arrary,$read_arrary[0]);
}
elsif ($flag_1 >0 && $flag_2 >0 ) {}
elsif ($flag_1 == 0 && $flag_2 == 0) {
#$str= join ( "\t", @temp_arrary);
#print $str."\n"; 
foreach $iterm (@temp_arrary) {
print $iterm.":".$hash{$iterm}.",";
}
print "\n";
@temp_arrary= @read_arrary;
$str="";
 }
$flag_1=0;
$flag_2=0;
} 
```

#### 5 Heritability estimation

Three data set of SNP were generated. The number of SNPs were consists.

```shell
for i in $(seq 1 100)  
do
gcta64 --bfile eQTL --extract snp_dataset1.txt --make-grm --out 01grmResult
gcta64 --bfile eQTL --extract snp_dataset2.txt --make-grm --out 02grmResult
gcta64 --bfile eQTL --extract snp_dataset3.txt --make-grm --out 03grmResult
gcta64 --reml --grm ../01grmResult --pheno ./phenotype/SI07an.txt --out ${i}_SI07an --reml-maxit 
10000 --reml-alg 2
done
```
