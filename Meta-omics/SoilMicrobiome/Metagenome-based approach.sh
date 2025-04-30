###### 1. Trimming raw reads using trimmomatic
conda activate
for f in $(cat SampleFullList)  # for each sample
do
    n=${f} 
    trimmomatic PE -threads 30 1_fastq/${n}_1.fastq.gz 1_fastq/${n}_2.fastq.gz \
    2_trimmed/PE2/${n}_1_trimmed.fastq.gz 2_trimmed/PE2/${n}_1_unpaired.fastq.gz 2_trimmed/PE2/${n}_2_trimmed.fastq.gz \
    2_trimmed/PE2/${n}_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

### 2. Quality check using FastQC
### 2-1. Before trimming
for f in $(cat SampleFullList)  
do
    n=${f} 
    fastqc -t 30 -o 1_fastq/QualityCheck 1_fastq/${n}_1.fastq.gz
    fastqc -t 30 -o 1_fastq/QualityCheck 1_fastq/${n}_2.fastq.gz
done

### 2-2. After trimming (none of overrepresented sequences were identified in the metagenomic data)
for f in $(cat SampleFullList)  
do
    n=${f} 
    fastqc -t 30 -o 2_trimmed/QualityCheck 2_trimmed/PE2/${n}_1_trimmed.fastq.gz
    fastqc -t 30 -o 2_trimmed/QualityCheck 2_trimmed/PE2/${n}_2_trimmed.fastq.gz
done

### 3. Assembly 
mkdir -p 4_assembly/megahit
### 3-1. Coassembly using megahit (default option with meta-large)
### Since an output directory is generated, do not need to make the output directory in advance.
### 3-1-1. TmD
megahit -1 ./2_trimmed/PE2/9D1PG_1_trimmed.fastq.gz,./2_trimmed/PE2/9D2PG_1_trimmed.fastq.gz,./2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz,./2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz,./2_trimmed/PE2/9D5PG_1_trimmed.fastq.gz,./2_trimmed/PE2/23D1PA_1_trimmed.fastq.gz,./2_trimmed/PE2/23D2PB_1_trimmed.fastq.gz,./2_trimmed/PE2/23D3PA_1_trimmed.fastq.gz,./2_trimmed/PE2/23D4PA_1_trimmed.fastq.gz,./2_trimmed/PE2/23D5PB_1_trimmed.fastq.gz \
        -2 ./2_trimmed/PE2/9D1PG_2_trimmed.fastq.gz,./2_trimmed/PE2/9D2PG_2_trimmed.fastq.gz,./2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz,./2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz,./2_trimmed/PE2/9D5PG_2_trimmed.fastq.gz,./2_trimmed/PE2/23D1PA_2_trimmed.fastq.gz,./2_trimmed/PE2/23D2PB_2_trimmed.fastq.gz,./2_trimmed/PE2/23D3PA_2_trimmed.fastq.gz,./2_trimmed/PE2/23D4PA_2_trimmed.fastq.gz,./2_trimmed/PE2/23D5PB_2_trimmed.fastq.gz \
        --presets meta-large --min-contig-len 300 -o 4_assembly/Coassembly/Dominant -t 20

### 3-1-2. TmM
megahit -1 ./2_trimmed/PE2/9D1NH_1_trimmed.fastq.gz,./2_trimmed/PE2/9D2NG_1_trimmed.fastq.gz,./2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz,./2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz,./2_trimmed/PE2/9D5NG_1_trimmed.fastq.gz,./2_trimmed/PE2/23D1NA_1_trimmed.fastq.gz,./2_trimmed/PE2/23D2NA_1_trimmed.fastq.gz,./2_trimmed/PE2/23D3NA_1_trimmed.fastq.gz,./2_trimmed/PE2/23D4NA_1_trimmed.fastq.gz,./2_trimmed/PE2/23D5NA_1_trimmed.fastq.gz \
        -2 ./2_trimmed/PE2/9D1NH_2_trimmed.fastq.gz,./2_trimmed/PE2/9D2NG_2_trimmed.fastq.gz,./2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz,./2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz,./2_trimmed/PE2/9D5NG_2_trimmed.fastq.gz,./2_trimmed/PE2/23D1NA_2_trimmed.fastq.gz,./2_trimmed/PE2/23D2NA_2_trimmed.fastq.gz,./2_trimmed/PE2/23D3NA_2_trimmed.fastq.gz,./2_trimmed/PE2/23D4NA_2_trimmed.fastq.gz,./2_trimmed/PE2/23D5NA_2_trimmed.fastq.gz \
        --presets meta-large --min-contig-len 300 -o 4_assembly/Coassembly/Minor -t 20

### 3-2. Individual assembly using megahit (default option with meta-large)
for f in $(cat SampleFullList) 
do
    n=${f}
    megahit -1 2_trimmed/PE2/${n}_1_trimmed.fastq.gz -2 2_trimmed/PE2/${n}_2_trimmed.fastq.gz --presets meta-large --min-contig-len 300 --continue -o 4_assembly/megahit_PE/MetaLarge/${n} -t 40
done

### 3-2. Change contig names (because of space in the seq names)
### 3-2-1. Replace spaces with underbars
for f in $(cat SampleFullList) 
do
    n=${f}
    sh RemoveSpaceInSeqNames.sh ./4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.fa > ./4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.removeSpace.fa
done

### 3-2-2. Remove redundant contigs using dedupe.sh
for f in $(cat SampleFullList) 
do
    n=${f}
    /run/media/root/data/bbmap/dedupe.sh in=./4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.removeSpace.fa out=./4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.deduplicated.fa
    rm -rf ./4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.removeSpace.fa
done

#### 3-2-3. Rename contigs
### 3-2-3-1. Check contig counts for each fasta file
for f in $(cat SampleFullList) 
do 
    grep -c ">" ./4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.deduplicated.fa
done

### 3-2-3-2. Rename contig headers for each fasta file
### 2022 samples
### TmD1
awk 'BEGIN{N=1;S=">TmD1IndAssemContig"}{if(/^>/){printf("%s%07d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D1PG/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D1PG/final.contigs.rename.fa
### TmD2
awk 'BEGIN{N=1;S=">TmD2IndAssemContig"}{if(/^>/){printf("%s%07d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D2PG/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D2PG/final.contigs.rename.fa
### TmD3
awk 'BEGIN{N=1;S=">TmD3IndAssemContig"}{if(/^>/){printf("%s%07d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D3PH/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D3PH/final.contigs.rename.fa
### TmD4
awk 'BEGIN{N=1;S=">TmD4IndAssemContig"}{if(/^>/){printf("%s%07d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D4PH/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D4PH/final.contigs.rename.fa
### TmD5
awk 'BEGIN{N=1;S=">TmD5IndAssemContig"}{if(/^>/){printf("%s%07d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D5PG/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D5PG/final.contigs.rename.fa
### TmM1
awk 'BEGIN{N=1;S=">TmM1IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D1NH/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D1NH/final.contigs.rename.fa
### TmM2
awk 'BEGIN{N=1;S=">TmM2IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D2NG/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D2NG/final.contigs.rename.fa
### TmM3
awk 'BEGIN{N=1;S=">TmM3IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D3NG/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D3NG/final.contigs.rename.fa
### TmM4
awk 'BEGIN{N=1;S=">TmM4IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D4NG/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D4NG/final.contigs.rename.fa
### TmM5
awk 'BEGIN{N=1;S=">TmM5IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/9D5NG/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/9D5NG/final.contigs.rename.fa

#### 2023 samples
grep -c ">" ./4_assembly/megahit_PE/MetaLarge/23DNCA/final.contigs.dedeplicated.fa
### Negative control
awk 'BEGIN{N=1;S=">Y23NegIndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23DNCA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23DNCA/final.contigs.rename.fa
### TmD1
awk 'BEGIN{N=1;S=">Y23TmD1IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D1PA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D1PA/final.contigs.rename.fa
### TmD2
awk 'BEGIN{N=1;S=">Y23TmD2IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D2PB/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D2PB/final.contigs.rename.fa
### TmD3
awk 'BEGIN{N=1;S=">Y23TmD3IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D3PA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D3PA/final.contigs.rename.fa
### TmD4
awk 'BEGIN{N=1;S=">Y23TmD4IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D4PA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D4PA/final.contigs.rename.fa
### TmD5
awk 'BEGIN{N=1;S=">Y23TmD5IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D5PB/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D5PB/final.contigs.rename.fa
### TmM1
awk 'BEGIN{N=1;S=">Y23TmM1IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D1NA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D1NA/final.contigs.rename.fa
### TmM2
awk 'BEGIN{N=1;S=">Y23TmM2IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D2NA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D2NA/final.contigs.rename.fa
### TmM3
awk 'BEGIN{N=1;S=">Y23TmM3IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D3NA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D3NA/final.contigs.rename.fa
### TmM4
awk 'BEGIN{N=1;S=">Y23TmM4IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D4NA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D4NA/final.contigs.rename.fa
### TmM5
awk 'BEGIN{N=1;S=">Y23TmM5IndAssemContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/megahit_PE/MetaLarge/23D5NA/final.contigs.deduplicated.fa > ./4_assembly/megahit_PE/MetaLarge/23D5NA/final.contigs.rename.fa

### Coassembled samples
### Dominant
awk 'BEGIN{N=1;S=">TmDContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/Coassembly/coassembled_contig_Dominant_deduplicated.fa > ./4_assembly/Coassembly/coassembled_contig_Dominant_rename.fa
### Minor
awk 'BEGIN{N=1;S=">TmMContig"}{if(/^>/){printf("%s%08d\n",S,N);N++}else{print $0}}' ./4_assembly/Coassembly/coassembled_contig_Minor_deduplicated.fa > ./4_assembly/Coassembly/coassembled_contig_Minor_rename.fa

### 4. Check quality of the assembled contigs using quast (without prepared reference genomes)
conda activate quastEnv
for f in $(cat SampleList)
do
    n=${f}
    python /home/miniconda3/envs/quastEnv/bin/metaquast.py --threads 5 \
    -o 4_assembly/QualityCheck/${n} \
    -l MegaHit,SPAdes 4_assembly/megahit/${n}/final.contigs.rename.fa 4_assembly/metaSpades/${n}/contigs.fasta
done
conda deactivate

### 5. Taxonomic assignment 
### Diamond using nrDB
### Individual assembly (in my work station)
#wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
#wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
#tar -xvzf taxdump.tar.gz
#wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

diamond makedb -p 20 --in nr.gz --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp -d nr231025.dmnd

#### 5-1. Individually assembled samples
for f in $(cat SampleFullList)  
do
    n=${f} 
    mkdir -p ./6_tax/Diamond/${n}
done

for f in $(cat SampleFullList) 
do
    n=${f}
    diamond blastx -d /data/DBforDiamond/nr231025.dmnd -q ./4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.rename.fa -o 6_tax/Diamond/${n}/contig_nr_diamond_tax.tsv -f 102
done

#### 5-2. Coassembled samples
###### 5-2-1. Split a fasta file into multiple files using seqkit
conda activate seqkit

seqkit split -s 707483 -O ./4_assembly/Coassembly/split/Dominant ./4_assembly/Coassembly/coassembled_contig_Dominant_rename.fa
seqkit split -s 818377 -O ./4_assembly/Coassembly/split/Minor ./4_assembly/Coassembly/coassembled_contig_Minor_rename.fa


###### 5-2-2. Run Diamond
for f in $(cat SampleCoassemblyList)  
do
    n=${f} 
    mkdir ./6_tax/Diamond/${n}
done

diamond blastx -d /data/DBforDiamond/nr231025.dmnd -q ./4_assembly/Coassembly/split/Dominant/coassembled_contig_Dominant_rename.part_001.fa -o ./6_tax/Diamond/Dominant/contig_part1_nr_diamond_tax.tsv -f 102
diamond blastx -d /data/DBforDiamond/nr231025.dmnd -q ./4_assembly/Coassembly/split/Dominant/coassembled_contig_Dominant_rename.part_002.fa -o ./6_tax/Diamond/Dominant/contig_part2_nr_diamond_tax.tsv -f 102
diamond blastx -d /data/DBforDiamond/nr231025.dmnd -q ./4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_001.fa -o ./6_tax/Diamond/Minor/contig_part1_nr_diamond_tax.tsv -f 102
diamond blastx -d /data/DBforDiamond/nr231025.dmnd -q ./4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_002.fa -o ./6_tax/Diamond/Minor/contig_part2_nr_diamond_tax.tsv -f 102
diamond blastx -d /data/DBforDiamond/nr231025.dmnd -q ./4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_003.fa -o ./6_tax/Diamond/Minor/contig_part3_nr_diamond_tax.tsv -f 102
diamond blastx -d /data/DBforDiamond/nr231025.dmnd -q ./4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_004.fa -o ./6_tax/Diamond/Minor/contig_part4_nr_diamond_tax.tsv -f 102

#### 5-2-2. Convert ncbi taxid to known lineages
### Accessed on 2024-07-15
ncbitax2lin --nodes-file /data/DBforMMSEQS2/tmp/14539138676643293848/taxonomy/nodes.dmp --names-file /data/DBforMMSEQS2/tmp/14539138676643293848/taxonomy/names.dmp

###### 6. Functional annotation at the community level ########
#### 6-1. Get bacterial and fungal contigs from each contig fasta file
### Bacterial and fungal contigs are divided based on the Diamond results using the R script 1_parse and get diamond tax result_updated.R.

#### 6-2. Structural annotation
##### 6-2-1. Individually assembled samples
#####  6-2-1-1. Prokka for bacterial contigs
mkdir -p ./7_annotation/ProteinClustering/{bacteria,fungi}/2_annotation
conda activate metaprokka
for f in $(cat SampleFullList) 
do
    n=${f}
    metaprokka --outdir ./7_annotation/ProteinClustering/bacteria/2_annotation/${n} --prefix metaprokka_res --locustag gene_${n} ./7_annotation/ProteinClustering/bacteria/1_contig/${n}_bacteria.contig.fa --cpus 20 --norrna
done

#####  6-2-1-2. Metaeuk for fungal contigs
conda activate metaeuk
## DB preparation
mmseqs databases UniRef90 /data/DBforMetaEuk/uniref90 /run/media/root/data1/DBforMetaEuk/tmp
# mmseqs databases NR /data/DBforMetaEuk/NR250110/nr /data/DBforMetaEuk/NR250110/tmp
# mmseqs databases UniProtKB /data/DBforMetaEuk/uniprot/uniprotKB /data/DBforMetaEuk/uniprot/tmp

## running
for f in $(cat SampleFullList) 
do
    n=${f}
    mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/${n}/tmp
done 

for f in $(cat SampleFullList) 
do
    n=${f}
    mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/NR/${n}/tmp
done 

for f in $(cat SampleFullList) 
do
    n=${f}
    mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/uniprot/${n}/tmp
done 

for f in $(cat SampleFullList) 
do
    n=${f}
    metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/${n}_fungi.contig.fa ./DBforMetaEuk/uniref90 ./7_annotation/ProteinClustering/fungi/2_annotation/${n}/output ./7_annotation/ProteinClustering/fungi/2_annotation/${n}/tmp
done

###### Annotation with nrDB
for f in $(cat SampleFullList) 
do
    n=${f}
    metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/${n}_fungi.contig.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/${n}/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/${n}/tmp
done

###### Annotation with uniprotKB
for f in $(cat SampleFullList) 
do
    n=${f}
    metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/${n}_fungi.contig.fa /data/DBforMetaEuk/uniprot/uniprotKB ./7_annotation/ProteinClustering/fungi/2_annotation/uniprot/${n}/output ./7_annotation/ProteinClustering/fungi/2_annotation/uniprot/${n}/tmp
done


conda activate seqkit
for f in $(cat SampleFullList) 
do
    n=${f}
    grep -c ">" ./7_annotation/ProteinClustering/fungi/1_contig/${n}_fungi.contig.fa
done

seqkit split -s 30000 -O ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG ./7_annotation/ProteinClustering/fungi/1_contig/9D1PG_fungi.contig.fa
seqkit split -s 30000 -O ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG ./7_annotation/ProteinClustering/fungi/1_contig/9D2PG_fungi.contig.fa
seqkit split -s 30000 -O ./7_annotation/ProteinClustering/fungi/1_contig/split/9D3PH ./7_annotation/ProteinClustering/fungi/1_contig/9D3PH_fungi.contig.fa
seqkit split -s 30000 -O ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1PA ./7_annotation/ProteinClustering/fungi/1_contig/23D1PA_fungi.contig.fa
seqkit split -s 30000 -O ./7_annotation/ProteinClustering/fungi/1_contig/split/23D3PA ./7_annotation/ProteinClustering/fungi/1_contig/23D3PA_fungi.contig.fa
seqkit split -s 30000 -O ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1NA ./7_annotation/ProteinClustering/fungi/1_contig/23D1NA_fungi.contig.fa


mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/{part1,part2,part3,part4,part5,part6,part7,part8,part9}/tmp
mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/{part1,part2,part3,part4,part5,part6,part7}/tmp
mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D3PH/{part1,part2,part3}/tmp
mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/{part1,part2,part3,part4,part5}/tmp
mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D3PA/{part1,part2,part3}/tmp
mkdir -p ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1NA/{part1,part2,part3}/tmp

metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_001.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part1/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part1/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_002.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part2/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part2/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_003.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part3/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part3/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_004.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part4/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part4/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_005.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part5/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part5/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_006.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part6/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part6/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_007.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part7/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part7/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_008.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part8/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part8/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D1PG/9D1PG_fungi.contig.part_009.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part9/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D1PG/part9/tmp

metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG/9D2PG_fungi.contig.part_001.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part1/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part1/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG/9D2PG_fungi.contig.part_002.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part2/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part2/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG/9D2PG_fungi.contig.part_003.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part3/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part3/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG/9D2PG_fungi.contig.part_004.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part4/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part4/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG/9D2PG_fungi.contig.part_005.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part5/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part5/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG/9D2PG_fungi.contig.part_006.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part6/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part6/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D2PG/9D2PG_fungi.contig.part_007.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part7/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D2PG/part7/tmp

metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D3PH/9D3PH_fungi.contig.part_001.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D3PH/part1/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D3PH/part1/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D3PH/9D3PH_fungi.contig.part_002.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D3PH/part2/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D3PH/part2/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/9D3PH/9D3PH_fungi.contig.part_003.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D3PH/part3/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/9D3PH/part3/tmp

metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1PA/23D1PA_fungi.contig.part_001.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part1/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part1/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1PA/23D1PA_fungi.contig.part_002.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part2/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part2/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1PA/23D1PA_fungi.contig.part_003.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part3/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part3/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1PA/23D1PA_fungi.contig.part_004.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part4/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part4/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1PA/23D1PA_fungi.contig.part_005.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part5/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1PA/part5/tmp

metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D3PA/23D3PA_fungi.contig.part_001.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D3PA/part1/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D3PA/part1/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D3PA/23D3PA_fungi.contig.part_002.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D3PA/part2/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D3PA/part2/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D3PA/23D3PA_fungi.contig.part_003.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D3PA/part3/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D3PA/part3/tmp

metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1NA/23D1NA_fungi.contig.part_001.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1NA/part1/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1NA/part1/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1NA/23D1NA_fungi.contig.part_002.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1NA/part2/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1NA/part2/tmp
metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/split/23D1NA/23D1NA_fungi.contig.part_003.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1NA/part3/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/23D1NA/part3/tmp


for f in $(cat SampleListLeftover) 
do
    n=${f}
    metaeuk easy-predict ./7_annotation/ProteinClustering/fungi/1_contig/${n}_fungi.contig.fa /data/DBforMMSEQS2/nrDB/nr ./7_annotation/ProteinClustering/fungi/2_annotation/NR/${n}/output ./7_annotation/ProteinClustering/fungi/2_annotation/NR/${n}/tmp
done


##### 6-2-2. Coassembled samples
conda activate metaprokka
metaprokka --outdir ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant --prefix metaprokka_res --locustag gene_TmDco ./4_assembly/Coassembly/Dominant/final.contigs.rename.bacteria.fa --cpus 40 --norrna

seqkit split -s 581925 -O ./4_assembly/Coassembly/Minor ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.fa
cat ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.part_004.fa ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.part_005.fa > ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.part_004_2.fa
metaprokka --outdir ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1 --prefix metaprokka_res --locustag gene_TmMcoP1 ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.part_001.fa --cpus 40 --norrna
metaprokka --outdir ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2 --prefix metaprokka_res --locustag gene_TmMcoP2 ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.part_002.fa --cpus 40 --norrna
metaprokka --outdir ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3 --prefix metaprokka_res --locustag gene_TmMcoP3 ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.part_003.fa --cpus 40 --norrna
metaprokka --outdir ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4 --prefix metaprokka_res --locustag gene_TmMcoP4 ./4_assembly/Coassembly/Minor/final.contigs.rename.bacteria.part_004_2.fa --cpus 40 --norrna

conda activate metaeuk
metaeuk easy-predict ./4_assembly/Coassembly/Dominant/final.contigs.rename.fungi.fa /data/DBforMetaEuk/uniref90 ./7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output ./7_annotation/ProteinClustering/fungi/2_annotation/Dominant/tmp
metaeuk easy-predict ./4_assembly/Coassembly/Minor/final.contigs.rename.fungi.fa /data/DBforMetaEuk/uniref90 ./7_annotation/ProteinClustering/fungi/2_annotation/Minor/output ./7_annotation/ProteinClustering/fungi/2_annotation/Minor/tmp


metaeuk easy-predict ./4_assembly/Coassembly/Dominant/final.contigs.rename.fungi.fa /data/DBforMetaEuk/uniprot/uniprotKB ./7_annotation/ProteinClustering/fungi/2_annotation/uniprot/Dominant/output ./7_annotation/ProteinClustering/fungi/2_annotation/uniprot/Dominant/tmp
metaeuk easy-predict ./4_assembly/Coassembly/Minor/final.contigs.rename.fungi.fa /data/DBforMetaEuk/uniprot/uniprotKB ./7_annotation/ProteinClustering/fungi/2_annotation/uniprot/Minor/output ./7_annotation/ProteinClustering/fungi/2_annotation/uniprot/Minor/tmp



#### 6-3. Rename fasta header to make consistency with id or TCS id in the gff files
#### 6-3-1. For bacteria
#### 6-3-1-1. Remove asterisks (non-amino acid bases) in the sequences
###### Individually assembled samples
for f in $(cat SampleFullList)
do
    n=${f}
    cat ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res.faa | perl -pe 's/\*//g' > ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res_noStar.faa
done
###### Coassembled samples
cat ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res.faa | perl -pe 's/\*//g' > ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res_noStar.faa

cat ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res.faa | perl -pe 's/\*//g' > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res_noStar.faa
cat ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res.faa | perl -pe 's/\*//g' > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res_noStar.faa
cat ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res.faa | perl -pe 's/\*//g' > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res_noStar.faa
cat ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res.faa | perl -pe 's/\*//g' > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res_noStar.faa

#### 6-3-1-2. Remove additional info in the fasta headers
###### Individually assembled samples
for f in $(cat SampleFullList) 
do
    n=${f}
    cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res_noStar.faa > ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res_rename.faa
    cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res.ffn > ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res_rename.ffn
done
###### Coassembled samples
cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res_noStar.faa > ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res_rename.faa
cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res.ffn > ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res_rename.ffn

cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res_noStar.faa > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res_rename.faa
cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res.ffn > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res_rename.ffn

cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res_noStar.faa > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res_rename.faa
cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res.ffn > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res_rename.ffn

cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res_noStar.faa > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res_rename.faa
cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res.ffn > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res_rename.ffn

cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res_noStar.faa > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res_rename.faa
cut -d" " -f1 ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res.ffn > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res_rename.ffn


cat ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res_rename.faa \
    ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res_rename.faa \
    ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res_rename.faa \
    ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res_rename.faa > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/combined/metaprokka_res_rename.faa

cat ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res_rename.ffn \
    ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res_rename.ffn \
    ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res_rename.ffn \
    ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res_rename.ffn > ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/combined/metaprokka_res_rename.ffn

#### 6-3-2. For fungi (more complicated than the bacterial case)
#### 6-3-2-1. Make reference name table using R (Get_name_table_for_changing_fasta_headers_for_fungi.R)

#### 6-3-2-2. Change header names using seqkit
###### Individually assembled samples
conda activate seqkit
for f in $(cat SampleFullList) 
do
    n=${f}
    seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output_rename.fas
    seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.codon.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.codon_rename.fas
done
###### Coassembled samples
seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output_rename.fas
seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.codon.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.codon_rename.fas

seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output_rename.fas
seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.codon.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.codon_rename.fas

##### 6-4. Get representative proteins using mmeseqs2 (protein clustering)
mkdir -p ./7_annotation/ProteinClustering/{bacteria,fungi}/{3_merge,4_cluster}/{gene,protein}
mkdir -p ./7_annotation/ProteinClustering/{bacteria,fungi}/4_cluster/{gene,protein}/tmp

######## 6-4-1. Concatenate individual and coassembly
### copy fasta files and change extension name: Individually assembled samples
for f in $(cat SampleFullList) 
do
    n=${f}
    cp ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res_rename.faa ./7_annotation/ProteinClustering/bacteria/3_merge/protein/metaprokka_res_${n}.faa
    cp ./7_annotation/ProteinClustering/fungi/2_annotation/${n}/output_rename.fas ./7_annotation/ProteinClustering/fungi/3_merge/protein/metaeuk_res_${n}.fas
done

### copy fasta files and change extension name: Coassembled samples
### bacteria
cp ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res_rename.faa ./7_annotation/ProteinClustering/bacteria/3_merge/protein/metaprokka_res_Dominant.faa
cp ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/combined/metaprokka_res_rename.faa ./7_annotation/ProteinClustering/bacteria/3_merge/protein/metaprokka_res_Minor.faa
### fungi
cp ./7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output_rename.fas ./7_annotation/ProteinClustering/fungi/3_merge/protein/metaeuk_res_Dominant.fas
cp ./7_annotation/ProteinClustering/fungi/2_annotation/Minor/output_rename.fas ./7_annotation/ProteinClustering/fungi/3_merge/protein/metaeuk_res_Minor.fas


### Merge all fasta files into a single one
### bacteria
cat ./7_annotation/ProteinClustering/bacteria/3_merge/protein/*.faa > ./7_annotation/ProteinClustering/bacteria/3_merge/protein/combined.bacteria.protein.faa
### fungi
cat ./7_annotation/ProteinClustering/fungi/3_merge/protein/*.fas > ./7_annotation/ProteinClustering/fungi/3_merge/protein/combined.fungi.protein.faa


######## 6-4-2. Cluster protein sequences using mmeseqs2
conda activate mmseqs2
#### bacteria
mmseqs easy-linclust ./7_annotation/ProteinClustering/bacteria/3_merge/protein/combined.bacteria.protein.faa ./7_annotation/ProteinClustering/bacteria/4_cluster/protein/clusterRes ./7_annotation/ProteinClustering/bacteria/4_cluster/protein/tmp \
                      --min-seq-id 0.95 -c 0.9 --threads 20 --cluster-mode 2 --cov-mode 1
#### fungi
mmseqs easy-linclust ./7_annotation/ProteinClustering/fungi/3_merge/protein/combined.fungi.protein.faa ./7_annotation/ProteinClustering/fungi/4_cluster/protein/clusterRes ./7_annotation/ProteinClustering/fungi/4_cluster/protein/tmp \
                     --min-seq-id 0.95 -c 0.9 --threads 20 --cluster-mode 2 --cov-mode 1

###### Individually assembled samples
conda activate seqkit
for f in $(cat SampleFullList) 
do
    n=${f}
    seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output_rename.fas
    seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.codon.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.codon_rename.fas
done
###### Coassembled samples
seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output_rename.fas
seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.codon.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.codon_rename.fas

seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output_rename.fas
seqkit replace -p "(.+)" -r '{kv}' -k /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.codon.fas > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.codon_rename.fas


######## 6-4-3. Extract genes corresponding to representative proteins from the merged gene sequence fasta file using seqkit
####### 6-4-3-1. concatenate identified gene fasta files
for i in $(cat SampleFullList)
do
    n=${f}
    cp ./7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res_rename.ffn ./7_annotation/ProteinClustering/bacteria/3_merge/gene/metaprokka_res_${n}.fna
    cp ./7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.codon_rename.fas ./7_annotation/ProteinClustering/fungi/3_merge/gene/metaeuk_res_${n}.fna
done

cp ./7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res_rename.ffn ./7_annotation/ProteinClustering/bacteria/3_merge/gene/metaprokka_res_Dominant.fna
cp ./7_annotation/ProteinClustering/bacteria/2_annotation/Minor/combined/metaprokka_res_rename.ffn ./7_annotation/ProteinClustering/bacteria/3_merge/gene/metaprokka_res_Minor.fna
cat ./7_annotation/ProteinClustering/bacteria/3_merge/gene/*.fna > ./7_annotation/ProteinClustering/bacteria/3_merge/gene/combined.bacteria.gene.fna

cp ./7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.codon_rename.fas ./7_annotation/ProteinClustering/fungi/3_merge/gene/metaeuk_res_Dominant.fas
cp ./7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.codon_rename.fas ./7_annotation/ProteinClustering/fungi/3_merge/gene/metaeuk_res_Minor.fas
cat ./7_annotation/ProteinClustering/fungi/3_merge/gene/*.fas > ./7_annotation/ProteinClustering/fungi/3_merge/gene/combined.fungi.gene.fna


####### 6-4-3-2. Extract genes corresponding to representative protein fasta files
conda activate seqkit
grep '^>' /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/clusterRes_rep_seq.fasta | awk '{print $1}' | sed 's/^>//' > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names
seqkit grep -n -f /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/3_merge/gene/combined.bacteria.gene.fna > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/gene/repre.bacteria.gene.fna

grep '^>' /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/clusterRes_rep_seq.fasta | awk '{print $1}' | sed 's/^>//' > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.names
seqkit grep -n -f /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.names /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/3_merge/gene/combined.fungi.gene.fna > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/gene/repre.fungi.gene.fna


##### 6-5. Find known orthologs using the EggNog-mapper (EggNog DB v.5.0)
mkdir -p ./7_annotation/ProteinClustering/{bacteria,fungi}/5_ortholog
conda activate eggnog
export EGGNOG_DATA_DIR=/data/DBforEggNog/
#### bacteria
emapper.py -m diamond --cpu 32 --itype proteins --allow_overlaps none --report_orthologs -i ./7_annotation/ProteinClustering/bacteria/4_cluster/protein/clusterRes_rep_seq.fasta --output_dir ./7_annotation/ProteinClustering/bacteria/5_ortholog --output eggnogRes
#### fungi
emapper.py -m diamond --cpu 32 --itype proteins --allow_overlaps none --report_orthologs -i ./7_annotation/ProteinClustering/fungi/4_cluster/protein/clusterRes_rep_seq.fasta --output_dir ./7_annotation/ProteinClustering/fungi/5_ortholog --output eggnogRes

####### 6-5-2. Parsing EggNog results (codes from Bin et al., 2023. Nat Commun)
#### Prepare GO Term data
python /data/ShellCodes/parse_go_obofile.py -i go-basic.obo -o go.tb

#### Parsing eggnog results (KEGG and GO terms)
python3 /data/ShellCodes/parse_eggNOG.py -i /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/5_ortholog/eggnogRes.emapper.annotations -g /data/SongyiMetagenome/7_annotation/Y2022/go.tb -o /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/5_ortholog
python3 /data/ShellCodes/parse_eggNOG.py -i /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/5_ortholog/eggnogRes.emapper.annotations -g /data/SongyiMetagenome/7_annotation/Y2022/go.tb -o /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/5_ortholog

##### 6-6. Estimation of gene abundance
####### 6-6-1. Get a single GTF file containing representative genes
######### 6-6-1-1. Export gene info from the originally annotated gff files (convert GFF to GTF)
######### bacteria
######### Individually assembled samples
mkdir -p /data/SongyiMetagenome/7_annotation/ProteinClustering/{bacteria,fungi}/6_abundance/{GTF,contig}
for f in $(cat SampleFullList)
do
    n=${f}
    sh /data/ShellCodes/prokkagff2gtf.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/2_annotation/${n}/metaprokka_res.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_${n}.gtf
done
######### Coassembled samples
sh /data/ShellCodes/prokkagff2gtf.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/2_annotation/Dominant/metaprokka_res.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Dominant.gtf
sh /data/ShellCodes/prokkagff2gtf.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part1/metaprokka_res.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part1.gtf
sh /data/ShellCodes/prokkagff2gtf.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part2/metaprokka_res.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part2.gtf
sh /data/ShellCodes/prokkagff2gtf.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part3/metaprokka_res.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part3.gtf
sh /data/ShellCodes/prokkagff2gtf.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/2_annotation/Minor/Part4/metaprokka_res.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part4.gtf

######### fungi
######### Individually assembled samples
for f in $(cat SampleFullList)
do
    n=${f}
    sh /data/ShellCodes/MetaEukgff2gtf_gene.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/${n}/output.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/metaeuk_res_${n}.gtf
done
######### Coassembled samples
sh /data/ShellCodes/MetaEukgff2gtf_gene.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Dominant/output.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/metaeuk_res_Dominant.gtf
sh /data/ShellCodes/MetaEukgff2gtf_gene.sh /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/2_annotation/Minor/output.gff > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/metaeuk_res_Minor.gtf

######### 6-6-1-2. Get GTF files containing only representative gene info 
mkdir -p /data/SongyiMetagenome/7_annotation/ProteinClustering/{bacteria,fungi}/6_abundance/GTF/filtered
######### bacteria
######### 1) Split protein list into each sample
for f in $(cat SampleFullList)
do
    n=${f}
    grep ${n} /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.${n}.names 
done
grep TmDco /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Dominant.names 

grep TmMcoP1 /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P1.names 
grep TmMcoP2 /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P2.names 
grep TmMcoP3 /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P3.names 
grep TmMcoP4 /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P4.names 

######### 2) Filter each GTF file
conda activate parallel
for f in $(cat SampleFullList)
do
    n=${f}
    cat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.${n}.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_${n}.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/metaprokka_res_${n}_filtered.gtf
done
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Dominant.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Dominant.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/metaprokka_res_Dominant_filtered.gtf

cat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P1.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part1.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/metaprokka_res_Minor_part1_filtered.gtf
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P2.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part2.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/metaprokka_res_Minor_part2_filtered.gtf
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P3.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part3.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/metaprokka_res_Minor_part3_filtered.gtf
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/protein/repre.prot.Minor.P4.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/metaprokka_res_Minor_part4.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/metaprokka_res_Minor_part4_filtered.gtf

######### 3) Merge filtered GTF files (since gene ids and annotation info were unique in each sample, use cat to merge GTF files)
cat ./7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/*.gtf > ./7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/rep.gene.bacteria.annot.gtf

######### fungi
conda activate parallel
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/*.gtf > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/combined.all.gene.fungi.gtf
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/combined.all.gene.fungi.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/rep.gene.fungi.annot.gtf

######### 1) Split protein list into each sample
for f in $(cat SampleFullList2)
do
    n=${f}
    grep ${n} /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.${n}.names 
done
grep TmDCon /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.Dominant.names 
grep TmMCon /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.names > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.Minor.names 

######### 2) Filter each GTF file
for f in $(cat SampleFullList2)
do
    n=${f}
    cat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.${n}.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/metaeuk_res_${n}.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/metaeuk_res_${n}_filtered.gtf
done
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.Dominant.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/metaeuk_res_Dominant.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/metaeuk_res_Dominant_filtered.gtf
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/protein/repre.prot.Minor.names | parallel -j 15 grep -w {} /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/metaeuk_res_Minor.gtf >> /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/metaeuk_res_Minor_filtered.gtf


######### 3) Merge filtered GTF files (since gene ids and annotation info were unique in each sample, use cat to merge GTF files)
cd ./7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/
cat *_filtered.gtf > rep.gene.fungi.annot.gtf
cd /data/SongyiMetagenome

####### 6-6-1-2. Convert GTF to SAF (for featureCounts)
awk '$3 == "CDS" {print $10 "\t" $1 "\t" $4 "\t" $5 "\t" $7}' /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/rep.gene.bacteria.annot.gtf | sed 's/gene_id "//;s/";//' > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/rep.gene.bacteria.annot.saf
echo -e "GeneID\tChr\tStart\tEnd\tStrand" | cat - /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/rep.gene.bacteria.annot.saf > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/temp && mv /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/temp /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/rep.gene.bacteria.annot.saf

awk '$3 == "CDS" {print $10 "\t" $1 "\t" $4 "\t" $5 "\t" $7}' /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/rep.gene.fungi.annot.gtf | sed 's/gene_id "//;s/";//' > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/rep.gene.fungi.annot.saf
echo -e "GeneID\tChr\tStart\tEnd\tStrand" | cat - /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/rep.gene.fungi.annot.saf > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/temp && mv /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/temp /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/rep.gene.fungi.annot.saf

####### 6-6-2. Extract contigs possessing representative genes
######## 6-6-2-1. Get a contig list
awk '{print $1}' /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/rep.gene.bacteria.annot.gtf |sort|uniq > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/bac.uniq.contig.with.rep.gene

awk '{print $1}' /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/rep.gene.fungi.annot.gtf |sort|uniq > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/fun.uniq.contig.with.rep.gene

######## 6-6-2-2. Get contig sequences from the originally assembled contig sequences
conda activate seqkit
######## bacteria
######## Individually assembled samples
for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f}
    seqkit grep --pattern-file /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/bac.uniq.contig.with.rep.gene /data/SongyiMetagenome/4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.rename.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/contig/contig.with.rep.gene_${n}.fa
done
######## Coassembled samples
seqkit grep --pattern-file /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/bac.uniq.contig.with.rep.gene /data/SongyiMetagenome/4_assembly/Coassembly/Dominant/final.contigs.rename.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/contig/contig.with.rep.gene_Dominant.fa
seqkit grep --pattern-file /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/GTF/filtered/bac.uniq.contig.with.rep.gene /data/SongyiMetagenome/4_assembly/Coassembly/Minor/final.contigs.rename.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/contig/contig.with.rep.gene_Minor.fa

####### Merge
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/contig/*.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/contig/combined.contigs.w.rep.genes.fa

######## fungi
######## Individually assembled samples
for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f}
    seqkit grep --pattern-file /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/fun.uniq.contig.with.rep.gene /data/SongyiMetagenome/4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.rename.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/contig/contig.with.rep.gene_${n}.fa
done
######## Coassembled samples
seqkit grep --pattern-file /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/fun.uniq.contig.with.rep.gene /data/SongyiMetagenome/4_assembly/Coassembly/Dominant/final.contigs.rename.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/contig/contig.with.rep.gene_Dominant.fa
seqkit grep --pattern-file /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/GTF/filtered/fun.uniq.contig.with.rep.gene /data/SongyiMetagenome/4_assembly/Coassembly/Minor/final.contigs.rename.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/contig/contig.with.rep.gene_Minor.fa

####### Merge
cat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/contig/*.fa > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/contig/combined.contigs.w.rep.genes.fa

####### 6-6-3. Mapping reads and estimating gene abundance 
for f in $(cat SampleFullList)  
do
    n=${f}
    mkdir -p ./7_annotation/ProteinClustering/{bacteria,fungi}/6_abundance/DNA/{1_mapping,2_count}/${n}
done

######## 6-6-3-1. Quantify gene abundance using salmon
###### Bacteria 
for f in $(cat /data/SongyiMetagenome/SampleFullList) 
do
    n=${f}  
    mkdir -p /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}
done


bwa index /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/gene/repre.bacteria.gene.fna


for f in $(cat /data/SongyiMetagenome/SampleFullList)  
do
    n=${f}
    bwa mem -t 30 /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/gene/repre.bacteria.gene.fna \
                  /data/SongyiMetagenome/2_trimmed/PE2/${n}_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/${n}_2_trimmed.fastq.gz > \
                  /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.sam
    samtools view -bS /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.sam > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.bam
    samtools sort /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.bam -o /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.sorted.bam
    samtools index /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.sorted.bam
    samtools flagstat /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.sorted.bam > /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.stat.txt
done

for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f}
    salmon quant -t /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/4_cluster/gene/repre.bacteria.gene.fna -l IU \
                 -a /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/bwa.sorted.bam \
                 -o /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/metagenome/${n}/salmon \
                 -p 30 --meta
done


###### Fungi
for f in $(cat /data/SongyiMetagenome/SampleFullList) 
do
    n=${f}  
    mkdir -p /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}
done


bwa index /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/gene/repre.fungi.gene.fna


for f in $(cat /data/SongyiMetagenome/SampleFullList)  
do
    n=${f}
    bwa mem -t 30 /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/4_cluster/gene/repre.fungi.gene.fna \
                  /data/SongyiMetagenome/2_trimmed/PE2/${n}_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/${n}_2_trimmed.fastq.gz > \
                  /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.sam
    samtools view -bS /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.sam > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.bam
    samtools sort /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.bam -o /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.sorted.bam
    samtools index /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.sorted.bam
    samtools flagstat /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.sorted.bam > /data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/6_abundance/metagenome/${n}/bwa.stat.txt
done
