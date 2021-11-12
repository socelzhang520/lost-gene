## 1棉花多倍体丢失的基因 

####1.1 模拟重测序数据比对祖先种

#####1.1.1 生成模拟数据

A亚组基因组大小GhAt 1.4G , 那么1.4 *1000 *1000000 bp的序列， 如果需要覆盖50 X，需要150bp reads

需要模拟产生resds数目 450,000,000

使用bedtools模拟150 bp的读长 

```shell
$seqkit fx2tab -l GhAt.fa |cut -f 1,4 > Atgenome.txt
$bedtools random -l 150 -n 450000000 -g Atgenome.txt > Gh_At.bed
$fastaFromBed -bed ./Gh_At.bed -name -fi GhAt.fa -fo GhAt_simulation.fa
$seqkit fx2tab -l GbAt.fa |cut -f 1,4 > Atgenome.txt
$bedtools random -l 150 -n 400000000 -g Atgenome.txt > Gb_At.bed
$fastaFromBed -bed ./Gb_At.bed -name -fi GbAt.fa -fo GbAt_simulation.fa
```

##### 1.2分别把陆地棉和海岛棉比对亚洲棉基因组

```shell
$ bwa index Garboreum_Shixiya1_WHUv3.0rc.genome.standard.fa &
$ nohup bwa mem -t 18 -M -R "@RG\tID:GhAt\tSM:GhAt\tPL:illumina\tLB:GhAt" ~/genome/updata_cotton/Garboreum_Shixiya1_WHUv3.0rc.genome.standard.fa ./GhAt_simulation.fa > GhAt.sam &
$ nohup bwa mem -t 18 -M -R "@RG\tID:GbAt\tSM:GbAt\tPL:illumina\tLB:GbAt" ~/genome/updata_cotton/Garboreum_Shixiya1_WHUv3.0rc.genome.standard.fa ./GbAt_simulation.fa > GbAt.sam &
 
$for i in {1..9}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ~/genome/updata_cotton/Garboreum_Shixiya1_WHUv3.0rc.genome.standard.fa -b ./afilename.txt -r Chr0${i} |bcftools call -vmO z -o GaChr0${i}.vcf.gz";done > Call_Asub_snp.sh
$for i in {10..13}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ~/genome/updata_cotton/Garboreum_Shixiya1_WHUv3.0rc.genome.standard.fa -b ./afilename -r Chr${i} |bcftools call -vmO z -o GaChr${i}.vcf.gz";done >> Call_Asub_snp.sh
```

#####1.3 snp的注释

```shell
gffread -T -o Ga.gtf ./Garboreum_Shixiya1_WHUv3.0rc.gene.standard.gff3
cut -f 9 Ga.gtf |perl -pi -e 's/transcript_id/gene_id/' > table1.txt
cut -f 9 Ga.gtf |perl -pi -e 's/transcript_id/gene_name/' > table2.txt
paste -d " " Ga.gtf table1.txt table2.txt > 1 && mv 1 Ga.gtf 
gtfToGenePred -genePredExt Ga.gtf Ga_V2.1_refGene.txt

~/biosoftware/annovar/retrieve_seq_from_fasta.pl --format refGene -seqfile ./Garboreum_Shixiya1_WHUv3.0rc.genome.standard.fa Ga_V2.1_refGene.txt --outfile Ga_refGeneMrna.fa

perl ~/biosoftware/annovar/table_annovar.pl ./merge.vcf cottondb/ --vcfinput -out final -buildver Ga_V2.1 --protocol refGene --operation g
```

##### 1.4 找寻致命基因变异并统计

```shell
head final.Ga_V2.1_multianno.txt
grep -v "synonymous" final.Ga_V2.1_multianno.txt|grep -v "intergenic" |grep -v "intronic" > Ga_result.txt
```



##### 1.5 计算外显子变异的基因和类型

##### Gh

```shell
## for Gh
grep "exonic" Ga_result.txt |cut -f 1,2,3,6,7,9,23|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' > GhA_exon.var.txt
grep -v "\./\." GhA_exon.var.txt|grep -v "0/0"  |wc -l
13737
grep -v "\./\." GhA_exon.var.txt|grep -v "0/0"  |cut -f 6|sort|uniq -c

2919 frameshift deletion
2172 frameshift insertion
1724 nonframeshift deletion
1382 nonframeshift insertion
4609 stopgain
 925 stoploss
 
## For Gb
grep "exonic" Ga_result.txt |cut -f 1,2,3,6,7,9,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' > GbA_exon.var.txt
grep -v "\./\." GbA_exon.var.txt|grep -v "0/0"|wc -l
 13873
 
grep -v "\./\." GbA_exon.var.txt|grep -v "0/0"|cut -f 6|sort|uniq -c

3081 frameshift deletion
2163 frameshift insertion
1744 nonframeshift deletion
1405 nonframeshift insertion
4519 stopgain
 955 stoploss
```



##### 基因结构变异的基因

```shell
grep "exonic" Ga_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.'

### A亚组多倍体化丢失基因
grep "exonic" Ga_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.' |perl -ne 'print if /1\/1\t1\/1/' > Gb_Gh_both_var.txt
### Gh亚组多倍体化丢失基因
grep "exonic" Ga_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.' |perl -ne 'print if /1\/1\t0\/0/' > Gh_specific_var.txt

### Gb亚组多倍体丢失基因

grep "exonic" Ga_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.' |perl -ne 'print if /0\/0\t1\/1/' > Gb_specific_var.txt


```



统计的数目

```shell
cut -f 5 Gb_specific_var.txt|sort|uniq > Gb_specific.list
cut -f 5 Gh_specific_var.txt|sort|uniq > Gh_specific.list
wc -l *specific.list
```















陆地棉D基因组

GhDt 794M. GbDt 789 M

0.794 *1000 *1000000 bp序列，

如果需要覆盖50X ，需要150 bp的read 

270,000,000

#####分别比对雷蒙德氏棉

```shell
$seqkit fx2tab -l GhDt.fa |cut -f 1,4 > Dtgenome.txt
$bedtools random -l 150 -n 270000000 -g Dtgenome.txt > Gh_Dt.bed
$fastaFromBed -bed ./Gh_Dt.bed -name -fi GhDt.fa -fo GhDt_simulation.fa
$seqkit fx2tab -l GbDt.fa |cut -f 1,4 > Dtgenome.txt
$bedtools random -l 150 -n 270000000 -g Dtgenome.txt > Gb_Dt.bed
$fastaFromBed -bed ./Gb_Dt.bed -name -fi GbDt.fa -fo GbDt_simulation.fa
$nohup bwa mem -t 18 -M -R "@RG\tID:GhDt\tSM:GhDt\tPL:illumina\tLB:GhDt" ~/genome/updata_cotton/Graimondii_221.fa ./GhDt_simulation.fa > GhDt.sam &
$nohup bwa mem -t 18 -M -R "@RG\tID:GbDt\tSM:GbDt\tPL:illumina\tLB:GbDt" ~/genome/updata_cotton/Graimondii_221.fa ./GbDt_simulation.fa > GbDt.sam &
 $for i in {1..9}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ~/genome/updata_cotton/Graimondii_221.fa -b ./dfilename.txt -r Chr0${i} |bcftools call -vmO z -o GrChr0${i}.vcf.gz";done > Call_Dsub_snp.sh
 $for i in {10..13}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ~/genome/updata_cotton/Graimondii_221.fa -b ./dfilename -r Chr${i} |bcftools call -vmO z -o GrChr${i}.vcf.gz";done >> Call_Dsub_snp.sh
```



```shell
gffread -T -o Gr.gtf ./Graimondii_221_gene.gff3
cut -f 9 Gr.gtf |perl -pi -e 's/transcript_id/gene_id/' > table1.txt
cut -f 9 Gr.gtf |perl -pi -e 's/transcript_id/gene_name/' > table2.txt
paste -d " " Gr.gtf table1.txt table2.txt > 1 && mv 1 Gr.gtf 
rm table1.txt table2.txt
gtfToGenePred -genePredExt Gr.gtf Gr_refGene.txt

~/biosoftware/annovar/retrieve_seq_from_fasta.pl --format refGene -seqfile ./Graimondii_221.fa Gr_refGene.txt --outfile Gr_refGeneMrna.fa

perl ~/biosoftware/annovar/table_annovar.pl ./merge.vcf cottondb/ --vcfinput -out final -buildver Gr --protocol refGene --operation g
```



##### 筛差异基因

```shell
grep -v "synonymous" final.Gr_multianno.txt|grep -v "intergenic" |grep -v "intronic" > Gr_result.txt

grep "exonic" Gr_result.txt |cut -f 1,2,3,6,7,9,23|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' > GhD_exon.var.txt
grep -v "\./\." GhD_exon.var.txt|grep -v "0/0"|wc -l
 15116
grep -v "\./\." GhD_exon.var.txt|grep -v "0/0"|cut -f 6|sort|uniq -c

3052 frameshift deletion
2237 frameshift insertion
2528 nonframeshift deletion
2387 nonframeshift insertion
4171 stopgain
 741 stoploss



grep "exonic" Gr_result.txt |cut -f 1,2,3,6,7,9,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' > GbD_exon.var.txt
grep -v "\./\." GbD_exon.var.txt|grep -v "0/0"|wc -l
15110

grep -v "\./\." GbD_exon.var.txt|grep -v "0/0"|cut -f 6|sort|uniq -c
3145 frameshift deletion
2233 frameshift insertion
2494 nonframeshift deletion
2376 nonframeshift insertion
4124 stopgain
 738 stoploss
 
```



基因结构变异

```shell
grep "exonic" Gr_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.'

### A亚组多倍体化丢失基因
grep "exonic" Gr_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.' |perl -ne 'print if /1\/1\t1\/1/' > Gb_Gh_both_var.txt
### Gh亚组多倍体化丢失基因
grep "exonic" Gr_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.' |perl -ne 'print if /1\/1\t0\/0/' > Gh_specific_var.txt

### Gb亚组多倍体丢失基因

grep "exonic" Gr_result.txt |cut -f 1,2,3,6,7,9,23,24|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' |grep -v '\./\.' |perl -ne 'print if /0\/0\t1\/1/' > Gb_specific_var.txt



```









### 油菜多倍体丢失基因

AA和CC

#####分别比对AA基因组

232M 0.232 *1000  * 1000000 * 50/150= 77,000,000

```shell
seqkit fx2tab -l AACC_A.fa |cut -f 1,4 > Atgenome.txt
bedtools random -l 150 -n 77000000 -g Atgenome.txt > AACC_At.bed
fastaFromBed -bed ./AACC_At.bed -name -fi AACC_A.fa -fo AACC_A_simulation.fa
bwa index Brassica_rapa.Brapa_1.0.dna.toplevel.fa

nohup bwa mem -t 18 -M -R "@RG\tID:AA\tSM:AA\tPL:illumina\tLB:AA" ./Brassica_rapa.Brapa_1.0.dna.toplevel.fa ./AACC_A_simulation.fa > AACC_A.sam &

~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./Brassica_rapa.Brapa_1.0.dna.toplevel.fa -b filename.txt|bcftools call -vmO z -o AA.vcf.gz 

```

snp注释

```shell
gffread -T -o AA.gtf ./Brassica_rapa.Brapa_1.0.48.gff3
cut -f 9 AA.gtf |perl -pi -e 's/transcript_id/gene_id/' > table1.txt
cut -f 9 AA.gtf |perl -pi -e 's/transcript_id/gene_name/' > table2.txt
paste -d " " AA.gtf table1.txt table2.txt > 1 && mv 1 AA.gtf 
gtfToGenePred -genePredExt AA.gtf AA_refGene.txt

~/biosoftware/annovar/retrieve_seq_from_fasta.pl --format refGene -seqfile ./Brassica_rapa.Brapa_1.0.dna.toplevel.fa AA_refGene.txt --outfile AA_refGeneMrna.fa


perl ~/biosoftware/annovar/table_annovar.pl ./AA.vcf AAdb/ --vcfinput -out final -buildver AA --protocol refGene --operation g


```







#####分别比对CC基因组

395M

0. 395*1000 *1000000 * 50/150= 132,000,000

 ```shell
seqkit fx2tab -l AACC_C.fa |cut -f 1,4 > Ctgenome.txt
bedtools random -l 150 -n 132000000 -g Ctgenome.txt > AACC_Ct.bed
fastaFromBed -bed ./AACC_Ct.bed -name -fi AACC_C.fa -fo AACC_C_simulation.fa
bwa index Brassica_oleracea.BOL.dna.toplevel.fa
nohup bwa mem -t 18 -M -R "@RG\tID:CC\tSM:CC\tPL:illumina\tLB:CC" ./Brassica_oleracea.BOL.dna.toplevel.fa ./AACC_C_simulation.fa > AACC_C.sam &


nohup ~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./Brassica_oleracea.BOL.dna.toplevel.fa -b cfilename.txt|bcftools call -vmO z -o CC.vcf.gz &


 ```



```shell
gffread -T -o CC.gtf ./Brassica_oleracea.BOL.47.gff3
cut -f 9 CC.gtf |perl -pi -e 's/transcript_id/gene_id/' > table1.txt
cut -f 9 CC.gtf |perl -pi -e 's/transcript_id/gene_name/' > table2.txt
paste -d " " CC.gtf table1.txt table2.txt > 1 && mv 1 CC.gtf 
gtfToGenePred -genePredExt CC.gtf CC_refGene.txt

~/biosoftware/annovar/retrieve_seq_from_fasta.pl --format refGene -seqfile ./Brassica_oleracea.BOL.dna.toplevel.fa CC_refGene.txt --outfile CC_refGeneMrna.fa
perl ~/biosoftware/annovar/table_annovar.pl ./CC.vcf CCdb/ --vcfinput -out final -buildver CC --protocol refGene --operation g


grep -v "synonymous" final.CC_multianno.txt|grep -v "intergenic" |grep -v "intronic" > CC_result.txt

grep "exonic" CC_result.txt |cut -f 1,2,3,6,7,9,23|perl -pi -e 's/:.*?\t/\t/'|perl -pi -e 's/:\d.*?\n/\n/' > CC_exon.var.txt
```

####验证数据

```shell
seqkit fx2tab -l AA.fa |cut -f 1,4 > Agenome.txt
head -n 10 Agenome.txt  > 1 && mv 1 Agenome.txt
bedtools random -l 150 -n 132000000 -g Agenome.txt > A.bed
nohup fastaFromBed -bed ./A.bed -name -fi AA.fa -fo AA_simulation.fa &
nohup bwa index AACC_A.fa &
nohup bwa index AACC_C.fa &

bwa mem -t 18 -M -R "@RG\tID:CC\tSM:CC\tPL:illumina\tLB:CC" ./AACC_A.fa ./AA_simulation.fa > AACC_AA.sam 


### C
seqkit fx2tab -l CC.fa |cut -f 1,4 |head -n 9 > Cgenome.txt
perl -pi -e 's/ dna.*?\t/\t/' Cgenome.txt
bedtools random -l 150 -n 132000000 -g Cgenome.txt > C.bed
nohup fastaFromBed -bed ./C.bed -name -fi CC.fa -fo CC_simulation.fa &



bwa mem -t 18 -M -R "@RG\tID:AA\tSM:AA\tPL:illumina\tLB:AA" ./AACC_A.fa ./AA_simulation.fa > AACC_AA.sam& 
bwa mem -t 18 -M -R "@RG\tID:CC\tSM:CC\tPL:illumina\tLB:CC" ./AACC_C.fa ./CC_simulation.fa > AACC_CC.sam& 


```







## 验证，使用亚洲棉和雷蒙德



```shell
nohup bwa index -a bwtsw GhAt.fa &
seqkit fx2tab -l Ga.fa |cut -f 1,4 > Gagenome.txt
bedtools random -l 150 -n 450000000 -g Gagenome.txt > Ga.bed
nohup fastaFromBed -bed ./Ga.bed -name -fi Ga.fa -fo Ga_simulation.fa &

seqkit fx2tab -l Gr.fa |cut -f 1,4 > Grgenome.txt
nohup bedtools random -l 150 -n 450000000 -g Grgenome.txt > Gr.bed &
nohup fastaFromBed -bed ./Gr.bed -name -fi Gr.fa -fo Gr_simulation.fa &



nohup bwa index -a bwtsw GhDt.fa &
nohup bwa index -a bwtsw GbAt.fa &
nohup bwa index -a bwtsw GbDt.fa &


bwa mem -t 18 -M -R "@RG\tID:GhDt\tSM:Ga\tPL:illumina\tLB:GaGhAt" ./GhAt.fa ./Ga_simulation.fa > GavsGhAt.sam 

bwa mem -t 18 -M -R "@RG\tID:GhDt\tSM:Ga\tPL:illumina\tLB:GrGhDt" ./GhDt.fa ./Gr_simulation.fa > GrvsGhDt.sam 

ls GavsGhAt.sam.sorted.bam > afilename.txt

for i in {1..9}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./GhAt.fa -b ./afilename.txt -r A0${i} |bcftools call -vmO z -o A0${i}.vcf.gz";done > CallA_snp.sh
for i in {10..13}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./GhAt.fa -b ./afilename.txt -r A${i} |bcftools call -vmO z -o A0${i}.vcf.gz";done >> CallA_snp.sh



for i in {1..9}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./GhDt.fa -b ./dfilename.txt -r D0${i} |bcftools call -vmO z -o D0${i}.vcf.gz";done > CallD_snp.sh
for i in {10..13}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./GhDt.fa -b ./dfilename.txt -r D${i} |bcftools call -vmO z -o D${i}.vcf.gz";done >> CallD_snp.sh









bwa mem -t 18 -M -R "@RG\tID:GhDt\tSM:Ga\tPL:illumina\tLB:GaGbAt" ./GbAt.fa ../GhAD/Ga_simulation.fa > GavsGbAt.sam 

bwa mem -t 18 -M -R "@RG\tID:GhDt\tSM:Ga\tPL:illumina\tLB:GrGbDt" ./GbDt.fa ../GhAD/Gr_simulation.fa > GrvsGbDt.sam 



for i in {1..9}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./GbDt.fa -b ./dfilename.txt -r D0${i} |bcftools call -vmO z -o D0${i}.vcf.gz";done > CallD_snp.sh
for i in {10..13}; do echo -e "~/biosoftware/samtools-1.6/samtools mpileup -q 20 -Q 15 -ugf ./GbDt.fa -b ./dfilename.txt -r D${i} |bcftools call -vmO z -o D${i}.vcf.gz";done >> CallD_snp.sh


```

