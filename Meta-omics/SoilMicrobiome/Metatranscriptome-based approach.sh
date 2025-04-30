#### 1_Trimmomatic.sh
mkdir 2_trimmed
mkdir 2_unpaired

## Copy an adapter sequence information to the present working directory
cp /home/miniconda3/pkgs/trimmomatic-*/share/trimmomatic-*/adapters/TruSeq3-PE-2.fa .
mamba activate /run/media/root/data/Anaconda3_envs/trimmomatic
###
## TmM
for f in $(cat SampleList) # for each sample
do
    n=${f} # strip part of file name
    trimmomatic PE -threads 30 ./0_rawdata/${n}_1.fastq.gz ./0_rawdata/${n}_2.fastq.gz \
    ./1_adapter/PE2/${n}_1_trimmed.fastq.gz ./1_adapter/PE2/${n}_1_unpaired.fastq.gz ./1_adapter/PE2/${n}_2_trimmed.fastq.gz \
    ./1_adapter/PE2/${n}_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

##TmD
for f in $(cat SampleList2) # for each sample
do
    n=${f} # strip part of file name
    trimmomatic PE -threads 30 ./0_rawdata/${n}_1.fastq.gz ./0_rawdata/${n}_2.fastq.gz \
    ./1_adapter/PE2/${n}_1_trimmed.fastq.gz ./1_adapter/PE2/${n}_1_unpaired.fastq.gz ./1_adapter/PE2/${n}_2_trimmed.fastq.gz \
    ./1_adapter/PE2/${n}_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

trimmomatic PE -threads 10 ./0_rawdata/23R1PB_1.fastq.gz ./0_rawdata/23R1PB_2.fastq.gz \
    ./1_adapter/PE2/23R1PB_1_trimmed.fastq.gz ./1_adapter/PE2/23R1PB_1_unpaired.fastq.gz ./1_adapter/PE2/23R1PB_2_trimmed.fastq.gz \
    ./1_adapter/PE2/23R1PB_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#### Quality check using FastQC
mamba activate /run/media/root/data1/Anaconda3_envs/fastqc
for f in $(cat SampleList)  
do
    n=${f} 
    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/${n}_1_trimmed.fastq.gz
    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/${n}_2_trimmed.fastq.gz
done

for f in $(cat SampleList2)  
do
    n=${f} 
    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/${n}_1_trimmed.fastq.gz
    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/${n}_2_trimmed.fastq.gz
done

##### Filtering overrepresented sequences and G polymers using the bbduk.sh
for f in $(cat SampleList)  
do
    n=${f} 
    /run/media/root/data1/bbmap/bbduk.sh in1=./1_adapter/PE2/${n}_1_trimmed.fastq.gz in2=./1_adapter/PE2/${n}_2_trimmed.fastq.gz out1=./1_adapter/PE2/${n}_1_trimmed_filtered.fastq.gz out2=./1_adapter/PE2/${n}_2_trimmed_filtered.fastq.gz entropy=0.5 entropywindow=50 entropyk=5
done


##### Quality check after the filtering (G polymers are removed after the bbduk.sh step)
for f in $(cat SampleList)  
do
    n=${f} 
    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/${n}_1_trimmed_filtered.fastq.gz
    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/${n}_2_trimmed_filtered.fastq.gz
done

    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/23R1PB_1_trimmed_filtered.fastq.gz
    fastqc -t 30 -o 1_1_qualitycheck/PE2/ 1_adapter/PE2/23R1PB_2_trimmed_filtered.fastq.gz

##### SortMeRNA
mamba activate sortmerna
for f in $(cat SampleFullList)  
do
    n=${f}
    sortmerna --ref /run/media/root/data/DBforSortMeRNA/smr_v4.3_default_db.fasta \
              --reads 1_adapter/PE2/${n}_1_trimmed_filtered.fastq.gz --reads 1_adapter/PE2/${n}_2_trimmed_filtered.fastq.gz \
              --aligned 2_sortmerna/PairedOut/${n}_aligned_PE --other 2_sortmerna/PairedOut/${n}_unaligned_PE \
              --fastx T --out2 --paired_out
    rm -rf /root/sortmerna/run/kvdb/
done


#### 6-6-3-2. RNA-based: Expression of representative genes using salmon
#####  Bacteria
for f in $(cat /data/SongyiMetatranscriptome/SampleList2023) 
do
    n=${f}  
    mkdir -p /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}
done

for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)  
do
    n=${f}
    bwa mem -t 30 /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/gene/repre.bacteria.gene.fna \
                  /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_1_trimmed.fastq.gz /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_2_trimmed.fastq.gz > \
                  /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.sam
    samtools view -bS /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.sam > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.bam
    samtools sort /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.bam -o /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.sorted.bam
    samtools index /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.sorted.bam
    samtools flagstat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.sorted.bam > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.stat.txt
done

#### Using salmon
for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)
do
    n=${f}
    salmon quant -t /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/gene/repre.bacteria.gene.fna -l ISR \
                 -a /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/bwa.sorted.bam \
                 -o /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metatranscriptome/${n}/salmon \
                 -p 30 --meta
done

##### Fungi
for f in $(cat /data/SongyiMetatranscriptome/SampleList2023) 
do
    n=${f}  
    mkdir -p /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}
done

for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)  
do
    n=${f}
    bwa mem -t 30 /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/gene/repre.fungi.gene.fna \
                  /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_1_trimmed.fastq.gz /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_2_trimmed.fastq.gz > \
                  /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.sam
    samtools view -bS /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.sam > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.bam
    samtools sort /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.bam -o /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.sorted.bam
    samtools index /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.sorted.bam
    samtools flagstat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.sorted.bam > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.stat.txt
done

#### Using salmon
for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)
do
    n=${f}
    salmon quant -t /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/gene/repre.fungi.gene.fna -l ISR \
                 -a /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/bwa.sorted.bam \
                 -o /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metatranscriptome/${n}/salmon \
                 -p 30 --meta
done

