##### Virome analysis from individually assembled and coassembled contigs
mamba create -n virsorter2 -c conda-forge -c bioconda virsorter=2
conda activate virsorter2
virsorter setup -d /var2/Hyun/230601_Songi_Metagenome/DBforVirsorter2/db2 -j 10
virsorter config --init-source --db-dir=/var2/Hyun/230601_Songi_Metagenome/DBforVirsorter2/db2
#### 1. retrieve DNA viruses from metagenome data (obtained from the individual assembly)
###### 1-1. VirSorter2
for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f} 
    mkdir -p ./DNAvirus/1_sorting/FromIndAssemble/1_1_virsorter/pass/${n}/tmp
done

for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f} 
    virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.rename.fa -w ./DNAvirus/1_sorting/FromIndAssemble/1_1_virsorter/pass/${n} --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromIndAssemble/1_1_virsorter/pass/${n}/tmp --min-length 1000 --min-score 0.5 -j 32 all
done

for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f} 
    virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/megahit_PE/MetaLarge/${n}/final.contigs.rename.fa -w ./DNAvirus/1_sorting/FromIndAssemble/1_1_virsorter/pass/${n} --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromIndAssemble/1_1_virsorter/pass/${n}/tmp --min-length 1000 --min-score 0.5 -j 16 all
done

###### 1-2. checkV
for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f} 
    mkdir ./DNAvirus/1_sorting/FromIndAssemble/1_2_checkv/${n}
done

for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f} 
    checkv end_to_end ./DNAvirus/1_sorting/FromIndAssemble/1_1_virsorter/pass/${n}/final-viral-combined.fa ./DNAvirus/1_sorting/FromIndAssemble/1_2_checkv/${n} -t 20 -d /data/ForestVirome/DBforCheckv/checkv-db-v1.5
done

###### 2. virome from co-assembled contigs
virsorter setup -d ./DBforVirsorter2/db -j 20

###### 2-1. VirSorter2
conda activate virsorter
mkdir -p ./DNAvirus/1_sorting/FromCoAssemble/Re/2_1_virsorter/pass

for f in $(cat /data/SongyiMetagenome/SampleCoassemblyList)
do
    n=${f} 
    mkdir ./DNAvirus/1_sorting/FromCoAssemble/Re/2_1_virsorter/pass/${n}
done

for f in $(cat /data/SongyiMetagenome/SampleCoassemblyList)
do
    n=${f} 
    virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/Coassembly/coassembled_contig_${n}_rename.fa -w ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/${n} --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/${n}/tmp --min-length 1000 --min-score 0.5 -j 20 all
done

virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/Coassembly/coassembled_contig_Dominant_rename.fa -w ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Dominant --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Dominant/tmp --min-length 1000 --min-score 0.5 -j 20 all
virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/Coassembly/coassembled_contig_Minor_rename.fa -w ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/tmp --min-length 1000 --min-score 0.5 -j 30 all

#### Splitted data for the Minor sample
virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_001.fa -w ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part1 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part1/tmp --min-length 1000 --min-score 0.5 -j 20 all
virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_002.fa -w ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part2 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part2/tmp --min-length 1000 --min-score 0.5 -j 20 all
virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_003.fa -w ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part3 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part3/tmp --min-length 1000 --min-score 0.5 -j 20 all
virsorter run --keep-original-seq -i /data/SongyiMetagenome/4_assembly/Coassembly/split/Minor/coassembled_contig_Minor_rename.part_004.fa -w ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part4 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae  --tmpdir ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part4/tmp --min-length 1000 --min-score 0.5 -j 20 all

#### Merge results
cat ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part1/final-viral-combined.fa \
    ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part2/final-viral-combined.fa \
    ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part3/final-viral-combined.fa \
    ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/Part4/final-viral-combined.fa > ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/combined/final-viral-combined.fa

###### 2-2. checkV
mkdir -p ./DNAvirus/1_sorting/FromCoAssemble/2_2_checkv/{Dominant,Minor}

checkv end_to_end ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Dominant/final-viral-combined.fa ./DNAvirus/1_sorting/FromCoAssemble/2_2_checkv/Dominant -t 20 -d ./DBforCheckv/checkv-db-v1.5
checkv end_to_end ./DNAvirus/1_sorting/FromCoAssemble/2_1_virsorter/pass/Minor/combined/final-viral-combined.fa ./DNAvirus/1_sorting/FromCoAssemble/2_2_checkv/Minor -t 40 -d ./DBforCheckv/checkv-db-v1.5

####### 3. filter qualified viral contigs based on checkv results
mkdir  /data/ForestVirome/DNAvirus/2_qc 
### R script "Filter viral contigs in R"
### qualified seqs will be saved in the directory /data/ForestVirome/DNAvirus/2_qc 

####### 4. vOTU generation
mkdir -p /data/ForestVirome/DNAvirus/3_cluster/{0_input,1_megablast,2_ANI,3_clustering}
cat /data/ForestVirome/DNAvirus/2_qc/*.fa > /data/ForestVirome/DNAvirus/3_cluster/0_input/combined.viral.genomes.fna

#mamba create -n blast bioconda::blast
conda activate blast
### 4-1. Perform all-vs-all BLAST using megablast utility, BLAST (v.2.15.0+) (From https://github.com/snayfach/MGV, accessed on Feb 23, 2024)
makeblastdb -in /data/ForestVirome/DNAvirus/3_cluster/0_input/combined.viral.genomes.fna -out /data/ForestVirome/DNAvirus/3_cluster/1_megablast/blastdb -dbtype nucl
blastn -query /data/ForestVirome/DNAvirus/3_cluster/0_input/combined.viral.genomes.fna -db /data/ForestVirome/DNAvirus/3_cluster/1_megablast/blastdb \
       -out /data/ForestVirome/DNAvirus/3_cluster/1_megablast/blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 19818 -perc_identity 90
### 4-2. Compute ANI from BLAST results
python3 /data/ForestVirome/MGV/ani_cluster/blastani.py -i /data/ForestVirome/DNAvirus/3_cluster/1_megablast/blast.tsv -o /data/ForestVirome/DNAvirus/3_cluster/2_ANI/ani.tsv
### 4-3. Perform centroid-based clustering
python3 /data/ForestVirome/MGV/ani_cluster/cluster.py --fna /data/ForestVirome/DNAvirus/3_cluster/0_input/combined.viral.genomes.fna --ani /data/ForestVirome/DNAvirus/3_cluster/2_ANI/ani.tsv --out /data/ForestVirome/DNAvirus/3_cluster/3_clustering/clusters.tsv --min_ani 95 --min_qcov 0 --min_tcov 85
#### 4-4. Get a fasta file consisting of representative vOTU sequences
### The first column of the clusters.tsv contains the representatives.
cut -f1 < /data/ForestVirome/DNAvirus/3_cluster/3_clustering/clusters.tsv > /data/ForestVirome/DNAvirus/3_cluster/3_clustering/repre_seqs.lst
conda activate seqtk
seqtk subseq /data/ForestVirome/DNAvirus/3_cluster/0_input/combined.viral.genomes.fna /data/ForestVirome/DNAvirus/3_cluster/3_clustering/repre_seqs.lst > /data/ForestVirome/DNAvirus/3_cluster/3_clustering/repre.votu.fna


### 5. Taxonomic and gene annotation using genomad (accessed Feb, 23, 2024)
#mamba create -n genomad -c conda-forge -c bioconda genomad
#conda activate genomad

#cd /data/ForestVirome/DBforgeNomad
#genomad download-database . ##db v.17

mkdir -p /data/ForestVirome/DNAvirus/4_annotation/{1_taxonomy,2_function}
conda activate genomad
genomad end-to-end --cleanup --splits 8 /data/ForestVirome/DNAvirus/3_cluster/3_clustering/repre.votu.fna /data/ForestVirome/DNAvirus/4_annotation/1_taxonomy/genomad /data/ForestVirome/DBforgeNomad/genomad_db

### Check the quality of the representative viral genomes again
checkv end_to_end /data/ForestVirome/DNAvirus/4_annotation/1_taxonomy/genomad/repre.votu_summary/repre.votu_virus.fna \
                  /data/ForestVirome/DNAvirus/4_annotation/1_taxonomy/genomad/repre.votu_summary/checkv/qc -t 20 -d /data/ForestVirome/DBforCheckv/checkv-db-v1.5


### 6. Host prediction using IPhoP (v1.3.3; accessed on Feb 28th, 2024)
### installation
mamba create -n iphop python=3.8 bioconda::iphop
### download DB
iphop download --db_dir /data/ForestVirome/DBforIPhoP/ -dbv iPHoP_db_Aug23_rw
### test running
wget https://bitbucket.org/srouxjgi/iphop/raw/d27b6bbdcd39a6a1cb8407c44ccbcc800d2b4f78/test/test_input_phages.fna
mkdir iphop_test_results
iphop predict --fa_file test_input_phages.fna --db_dir /data/ForestVirome/DBforIPhoP/iPHoP_db_Aug23_rw/ --out_dir iphop_test_results/test_input_phages_iphop

### 6-1. Running
mkdir -p ./DNAvirus/5_host_predic/output
iphop predict --fa_file /data/ForestVirome/DNAvirus/4_annotation/1_taxonomy/genomad/repre.votu_summary/repre.votu_virus.fna  --db_dir /data/ForestVirome/DBforIPhoP/Aug_2023_pub_rw --out_dir /data/ForestVirome/DNAvirus/5_host_predic/output


#### 6-2. host prediction using PHP
##### 6-2-1. decompress a host k-mer file (provided by the author)
tar -zxvf /data/ForestVirome/DNAvirus/5_host_predic/php/PHP/hostKmer_60105_kmer4.tar.gz

##### 6-2-2. split viral genome fasta to each files
conda activate seqkit
seqkit split /data/ForestVirome/DNAvirus/4_annotation/1_taxonomy/genomad/repre.votu_summary/repre.votu_virus.fna \
       -O /data/ForestVirome/DNAvirus/5_host_predic/php/input/fasta --by-id

##### 6-2-3. predict phage host
conda activate phagehostpredictor
python3 /data/ForestVirome/DNAvirus/5_host_predic/php/PHP/PHP.py \
        -v /data/ForestVirome/DNAvirus/5_host_predic/php/input/fasta  \
        -o /data/ForestVirome/DNAvirus/5_host_predic/php/output \
        -d /data/ForestVirome/DNAvirus/5_host_predic/php/PHP/host_kmer -n hostKmer_60105_kmer4


### 7. Viral abundance
### 7-1. DNA-based
### Mapped read count
coverm genome -d /data/ForestVirome/DNAvirus/6_abundance/1_input \
              -c /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23DNCA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23DNCA_2_trimmed.fastq.gz \
                 --min-covered-fraction 0 \
              -o /data/ForestVirome/DNAvirus/6_abundance/2_estimate/DNA/vmag.mapped.reads.tsv -p bwa-mem2 -m count

### RPKM
coverm genome -d /data/ForestVirome/DNAvirus/6_abundance/1_input \
              -c /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23DNCA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23DNCA_2_trimmed.fastq.gz \
                 --min-covered-fraction 0 \
              -o /data/ForestVirome/DNAvirus/6_abundance/2_estimate/DNA/vmag.mapped.reads.tsv -p bwa-mem2 -m rpkm



######## Gene abundance and expression
#### Gene abundance
for f in $(cat /data/SongyiMetagenome/SampleFullList) 
do
    n=${f}  
    mkdir -p /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}
done


bwa index /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.protein.meta.trueVirus.fna


for f in $(cat /data/SongyiMetagenome/SampleFullList)  
do
    n=${f}
    bwa mem -t 30 /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.protein.meta.trueVirus.fna \
                  /data/SongyiMetagenome/2_trimmed/PE2/${n}_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/${n}_2_trimmed.fastq.gz > \
                  /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.sam
    samtools view -bS /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.sam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.bam
    samtools sort /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.bam -o /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.sorted.bam
    samtools index /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.sorted.bam
    samtools flagstat /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.sorted.bam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.stat.txt
done

for f in $(cat /data/SongyiMetagenome/SampleFullList)
do
    n=${f}
    salmon quant -t /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.protein.meta.trueVirus.fna -l IU \
                 -a /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/bwa.sorted.bam \
                 -o /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/DNA/${n}/salmon \
                 -p 30 --meta
done

#### Gene expression
for f in $(cat /data/SongyiMetatranscriptome/SampleList2023) 
do
    n=${f}  
    mkdir -p /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}
done

for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)  
do
    n=${f}
    bwa mem -t 30 /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.protein.meta.trueVirus.fna \
                  /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_1_trimmed.fastq.gz /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_2_trimmed.fastq.gz > \
                  /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.sam
    samtools view -bS /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.sam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.bam
    samtools sort /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.bam -o /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.sorted.bam
    samtools index /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.sorted.bam
    samtools flagstat /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.sorted.bam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.stat.txt
done

#### Using salmon
for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)
do
    n=${f}
    salmon quant -t /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.protein.meta.trueVirus.fna -l ISR \
                 -a /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.sorted.bam \
                 -o /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/salmon \
                 -p 30 --meta
done


############ Find connection between phageome and bacrterial community
##### For the bacterial community, bacterial species abundances were calculated at the species level (merge contig abundances into species abunadnces)
###### get taxonomic information of representative prokaryotic contigs from the previously obtained diamond results
conda activate seqkit
seqkit seq -n /data/SongyiMetagenome_S2/1_contigs/2_clustered/clustered.id95_rep_seq.prok.fasta > /data/SongyiMetagenome_S2/1_contigs/2_clustered/repre.prokaryotic.contig.headers.txt

seqkit grep -f /data/SongyiMetagenome_S2/1_contigs/2_clustered/final.prokaryote.host.contig.list.txt /data/SongyiMetagenome_S2/1_contigs/2_clustered/clustered.id95_rep_seq.prok.fasta \
            -o /data/SongyiMetagenome_S2/1_contigs/2_clustered/prokaryote.matched.fasta


seqkit seq -n /data/SongyiMetagenome_S2/1_contigs/2_clustered/clustered.id95_rep_seq.euk.fasta > /data/SongyiMetagenome_S2/1_contigs/2_clustered/repre.eukaryotic.contig.headers.txt

seqkit grep -f /data/SongyiMetagenome_S2/1_contigs/2_clustered/final.eukaryote.host.contig.list.txt /data/SongyiMetagenome_S2/1_contigs/2_clustered/clustered.id95_rep_seq.fasta \
            -o /data/SongyiMetagenome_S2/1_contigs/2_clustered/eukaryote.matched.fasta

##### Abundance of prokaryotic contigs
#### Metagenome abundance
coverm contig -r /data/SongyiMetagenome_S2/1_contigs/2_clustered/prokaryote.matched.fasta \
              -c /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_2_trimmed.fastq.gz \
              -o /data/ForestVirome/DNAvirus/6_abundance/2_estimate/metag.rpkm.prokaryote.abundance_bwa.tsv -p bwa-mem -m rpkm -t 24


coverm contig -r /data/SongyiMetagenome_S2/1_contigs/2_clustered/prokaryote.matched.fasta \
              -c /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_2_trimmed.fastq.gz \
              -o /data/ForestVirome/DNAvirus/6_abundance/2_estimate/metag.read.prokaryote.abundance_bwa.tsv -p bwa-mem -m rpkm -t 24

### 8-1-1. Mapping short reads to representative viral contigs
for f in $(cat /data/SongyiMetagenome/SampleFullList)  
do
    n=${f}
    mkdir -p /data/ForestVirome/DNAvirus/5_abundance/{1_mapping,2_conut}/${n}
done

conda activate /data/Anaconda3_envs/samtoolsMamba

bwa index /data/ForestVirome/DNAvirus/3_cluster/3_clustering/repre.votu.fna
for f in $(cat /data/SongyiMetagenome/SampleFullList)  
do
    n=${f}
    bwa mem -t 10 ./DNAvirus/3_cluster/3_clustering/repre.votu.fna /data/SongyiMetagenome/2_trimmed/PE2/${n}_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/${n}_2_trimmed.fastq.gz > ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sam
    samtools view -bS ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sam > ./DNAvirus/5_abundance/1_mapping/${n}/bwa.bam
    samtools sort ./DNAvirus/5_abundance/1_mapping/${n}/bwa.bam -o ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sorted.bam
    samtools index ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sorted.bam
    samtools flagstat ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sorted.bam > ./DNAvirus/5_abundance/1_mapping/${n}/stat.txt
done

### 8-1-2. Get mapped contig length
conda activate bedtools
bedtools bamtobed -i ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sorted.bam > ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sorted.bed
awk '{print $1, $3-$2}' ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sorted.bed > ./DNAvirus/5_abundance/1_mapping/${n}/mapped_lengths.txt
awk '{contig_lengths[$1]+=$2} END {for (contig in contig_lengths) print contig, contig_lengths[contig]}' ./DNAvirus/5_abundance/1_mapping/${n}/mapped_lengths.txt > ./DNAvirus/5_abundance/1_mapping/${n}/total_mapped_lengths_per_contig.txt

#### 8-1-3. Get read counts mapped to each contig
for f in $(cat /data/SongyiMetagenome/SampleFullList)
do 
    n=${f}
    samtools view ./DNAvirus/5_abundance/1_mapping/${n}/bwa.sorted.bam | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c > ./DNAvirus/5_abundance/1_mapping/${n}/ReadCountMappedtoContigs.txt
done

samtools view bamfile.bam | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c > ./DNAvirus/5_abundance/1_mapping/DNA/${n}/ReadCountMappedtoContigs.txt
samtools view -F 0x100 -F 0x800 bamfile.bam | cut -f3 | sort | uniq -c > ./DNAvirus/5_abundance/1_mapping/DNA/${n}/ReadCountMappedtoContigs.txt

#### Considering strandedness
### Forward-stranded
samtools view -f 0x2 -F 0x100 -F 0x800 bamfile.bam | awk '{if ($2 & 0x40) { # First in pairif ($2 & 0x10) { print $3 }} else if ($2 & 0x80) { if (!($2 & 0x10)) { print $3 } }}' | sort | uniq -c > ./DNAvirus/5_abundance/1_mapping/RNA/${n}/ReadCountMappedtoContigs.txt
### reversed-stranded (we will use this.)
samtools view -f 0x2 -F 0x100 -F 0x800 bamfile.bam | awk '{if ($2 & 0x40) {if (!($2 & 0x10)) { print $3 }} else if ($2 & 0x80) {if ($2 & 0x10) { print $3 }}}' | sort | uniq -c > ./DNAvirus/5_abundance/1_mapping/RNA/${n}/ReadCountMappedtoContigs.txt

#### 8-1-4. Get the length of each viral contig
seqkit seq --name --length contigs.fasta

#### 8-1-5. Normalize read counts to TPM
for f in $(cat SampleFullList)
do
    n=${f}
    htseq-count -r pos -s no -t CDS -f bam ./7_annotation/ProteinClustering/bacteria/6_abundance/DNA/1_mapping/${n}/bwa.sorted.bam /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/0_gtf/rep.gene.bacteria.annot.gtf > ./7_annotation/ProteinClustering/bacteria/6_abundance/DNA/2_count/${n}/allcontigs.gene.count
    cut -f4,5,9 /data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/6_abundance/0_gtf/rep.gene.bacteria.annot.gtf | sed 's/gene_id //g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' >  ./7_annotation/ProteinClustering/bacteria/6_abundance/DNA/2_count/${n}/allcontigs.gene.genelengths
    /data/ShellCodes/tpm_table.py -n allcontigs.gene -c ./7_annotation/ProteinClustering/bacteria/6_abundance/DNA/2_count/${n}/allcontigs.gene.count -i <(echo -e "allcontigs.gene\t100") -l ./7_annotation/ProteinClustering/bacteria/6_abundance/DNA/2_count/${n}/allcontigs.gene.genelengths > ./7_annotation/ProteinClustering/bacteria/6_abundance/DNA/2_count/${n}/allcontigs.gene.tpm
done

htseq-count -r pos -s reverse -t CDS -f bam bamfile.bam map.gtf 

###### 8-2. Calculate genome coverages using CoverM (v.0.7.0; accessed May 9, 2024)
### Split genomes firit
mkdir -p /data/ForestVirome/DNAvirus/6_abundance/{1_input,2_estimate}
conda activate seqkit
seqkit split -i -O /data/ForestVirome/DNAvirus/6_abundance/1_input /data/ForestVirome/DNAvirus/3_cluster/3_clustering/repre.votu.trueVirus.fna 

### Estimate abundance using CoverM
#### relative abundance
coverm genome -d /data/ForestVirome/DNAvirus/6_abundance/1_input \
              -c /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_2_trimmed.fastq.gz \
              -o /data/ForestVirome/DNAvirus/6_abundance/2_estimate/metag.abundance.tsv

##### TPM (mapper: bwa)
coverm genome -d /data/ForestVirome/DNAvirus/6_abundance/1_input \
              -c /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4PH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5PG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D1NH_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D2NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D3NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D4NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/9D5NG_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4PA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5PB_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D1NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D2NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D3NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D4NA_2_trimmed.fastq.gz \
                 /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_1_trimmed.fastq.gz /data/SongyiMetagenome/2_trimmed/PE2/23D5NA_2_trimmed.fastq.gz \
              -o /data/ForestVirome/DNAvirus/6_abundance/2_estimate/metag.tpm.abundance_bwa.tsv -p bwa-mem2 -m tpm




################ AMG curation ##########
###### Target ORF
/data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.faa

#### DB preparation
##  VOG Feb 5 2025 version 228
cd /data/ForestVirome/VOG/v228
wget https://fileshare.lisc.univie.ac.at/vog/vog228/vog.hmm.tar.gz
tar -xvzf vog.hmm.tar.gz

### PHROGs
cd /data/ForestVirome/PHROGs
wget https://phrogs.lmge.uca.fr/downloads_from_website/HMM_phrog.tar.gz
tar -xvzf HMM_phrog.tar.gz

tar -xvzf MSA_phrogs.tar.gz

### pVOGs
cd /data/ForestVirome/pVOGs
wget https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
tar -xvzf AllvogHMMprofiles.tar.gz

###### Find viral hallmark genes from VOG, CheckV, and PHROGs
#### Clean headers
seqkit replace -p " #.*" -r "" /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.faa \
               -o /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa

##### VOG
find /data/ForestVirome/VOG/v228/hmm/ -name "*.hmm" | xargs cat > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/1_vog/combined_hmm.hmm
hmmpress /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/1_vog/combined_hmm.hmm
hmmsearch --tblout /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/1_vog/hmm_search_res.txt \
          --cpu 20 --noali /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/1_vog/combined_hmm.hmm /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa

##### CheckV
cat /data/ForestVirome/DBforCheckv/checkv-db-v1.5/hmm_db/checkv_hmms/*.hmm > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/2_checkv/combined_hmm.hmm
hmmpress /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/2_checkv/combined_hmm.hmm
hmmsearch --tblout /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/2_checkv/hmm_search_res.txt --cpu 20  --noali \
                   /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/2_checkv/combined_hmm.hmm /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa

##### pVOGs
mkdir /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs
cat /data/ForestVirome/pVOGs/AllvogHMMprofiles/*.hmm > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs/combined_hmm.hmm
hmmpress /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs/combined_hmm.hmm
hmmsearch --tblout /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs/hmm_search_res.txt --cpu 10  --noali \
                   /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs/combined_hmm.hmm /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa

###### Filtering based on the parameters e < 1e-5 & Bit Score > 50
awk '$5 <= 1e-5 && $6 >= 50' /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/1_vog/hmm_search_res.txt > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/1_vog/hmm_search_res_filtered.txt
awk '$5 <= 1e-5 && $6 >= 50' /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/2_checkv/hmm_search_res.txt > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/2_checkv/hmm_search_res_filtered.txt
awk '$5 <= 1e-5 && $6 >= 50' /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs/hmm_search_res.txt > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs/hmm_search_res_filtered.txt

cat /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/3_pVOGs/hmm_search_res_filtered.txt \
    /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/2_checkv/hmm_search_res_filtered.txt | sort -k1,1 -k6,6nr > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/4_merge/viral_hallmark_genes.txt
cut -f1 /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/4_merge/viral_hallmark_genes.txt | sort | uniq -c | sort -nr > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/4_merge/viral_gene_counts.txt

##### Define viral and non-viral ORFs
grep ">" /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa | sed 's/>//' > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/all_orfs.txt
awk '{print $1}' /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/4_merge/viral_hallmark_genes.txt | sort | uniq > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/viral_orfs.txt
seqkit grep -v -f /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/viral_orfs.txt /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa -o /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.faa

##### non-viral ORFs in viral contigs
grep ">" /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.faa | sed 's/>//' > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.txt
awk -F'_' '{OFS="_"; NF--; print}' /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/viral_orfs.txt  | sort | uniq > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/viral_contigs.txt

### obtain all ORFs in viral contigs
grep -F -f /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/viral_contigs.txt /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/all_orfs.txt > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/orfs_in_viral_contigs.txt
seqkit grep -f /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/orfs_in_viral_contigs.txt /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa \
            -o /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/orfs_in_viral_contigs.faa

grep -vF -f /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/viral_orfs.txt /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/orfs_in_viral_contigs.txt > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/amg_candidate_orfs.txt

seqkit grep -f  /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/amg_candidate_orfs.txt /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa -o  /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/AMG_candidate_orfs.faa


##### Kofamscan to search KEGG orthologs for each putative AMG
cd /data/DBforKOfamscan/Feb2025
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
tar -xvzf profiles.tar.gz
gunzip ko_list.gz

exec_annotation --cpu 32 -o /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/kofamscan_results.txt \
                --tmp-dir /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/tmp \
                -p /data/DBforKOfamscan/Feb2025/profiles -k /data/DBforKOfamscan/Feb2025/ko_list /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/AMG_candidate_orfs.faa

exec_annotation --cpu 32 -o /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/kofamscan_results_allNonviralORFs.txt \
                --tmp-dir /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/tmp \
                -p /data/DBforKOfamscan/Feb2025/profiles -k /data/DBforKOfamscan/Feb2025/ko_list  /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.faa

###### significant hits
awk '$5 <= 1e-5 && $4 >= 50' /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/kofamscan_results.txt > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/filtered_kofamscan_hits.txt
awk '$5 <= 1e-5 && $4 >= 50' /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/kofamscan_results_allNonviralORFs.txt > /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/6_kofam/filtered_kofamscan_results_allNonviralORFs.txt

###### PHROGs annotation
##### Using mmseqs2
#(1) Create a database with your fasta file :
mmseqs createdb /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/repre.votu.protein.cleanHeader.query

#(2) Compute the search and convert the results into a tab separated file :
mmseqs search /data/ForestVirome/DBforViralmmseqs/PHROGs/phrogs_mmseqs_db/phrogs_profile_db /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/repre.votu.protein.cleanHeader.query /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/phrogs.res /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/tmp -s 7
mmseqs createtsv /data/ForestVirome/DBforViralmmseqs/PHROGs/phrogs_mmseqs_db/phrogs_profile_db /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/repre.votu.protein.cleanHeader.query /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/phrogs.res /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/phrogs.results.tsv
mmseqs createtsv /data/ForestVirome/DBforViralmmseqs/PHROGs/phrogs_mmseqs_db/phrogs_profile_db /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/repre.votu.protein.cleanHeader.query /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/phrogs.res /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/7_phrogs/phrogs.result.fullheader.tsv --full-header




mmseqs databases UniRef90 /var2/Hyun/DBforMMseqs2/UniRef90/uniref90 /var2/Hyun/DBforMMseqs2/UniRef90/tmp 
mmseqs easy-search /var2/Hyun/ForestVirome/3_annotation/0_genes/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa /var2/Hyun/DBforMMseqs2/UniRef90/uniref90 \
                   /var2/Hyun/ForestVirome/3_annotation/1_uniref90/uniref90.res /var2/Hyun/ForestVirome/3_annotation/1_uniref90/tmp -s 7

mmseqs easy-search /var2/Hyun/ForestVirome/3_annotation/0_genes/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa /var2/Hyun/DBforMMseqs2/UniRef100/uniref100 \
                   /var2/Hyun/ForestVirome/3_annotation/1_uniref90/uniref100.coverage50.res /var2/Hyun/ForestVirome/3_annotation/1_uniref90/tmp -s 7 -c 0.5 --cov-mode 0

mmseqs easy-search /var2/Hyun/ForestVirome/3_annotation/0_genes/repre.votu.protein.meta.trueVirus.noStar.cleanHeader.faa /var2/Hyun/DBforMMseqs2/dbCAN/dbCAN \
                   /var2/Hyun/ForestVirome/3_annotation/1_uniref90/dbCAN.coverage50.res /var2/Hyun/ForestVirome/3_annotation/1_uniref90/tmp -s 7 -c 0.5 --cov-mode 0

sort annotated.ids.txt -o annotated.ids.sorted.txt
sort -k1,1 extracted_headers.txt -o extracted_headers.sorted.txt
join -t $'\t' annotated.ids.sorted.txt extracted_headers.sorted.txt > matched.txt


awk 'NR==FNR { ids[$1]; next } $1 in ids { print $1 "\t" $2 }' annotated_ids.txt extracted_headers.txt > matched_headers_descriptions.txt


##### dbcan (CAZyme)
mamba create -n dbcan python=3.8
conda activate dbcan
mamba install dbcan -c conda-forge -c bioconda


cd /data/DBforDBcan/db
wget http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08012023.tsv
mv fam-substrate-mapping-08012023.tsv fam-substrate-mapping.tsv

wget http://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa
makeblastdb -in PUL.faa -dbtype prot

wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.xlsx
mv dbCAN-PUL_12-12-2023.xlsx dbCAN-PUL.xlsx

wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz
tar xvf dbCAN-PUL.tar.gz
rm dbCAN-PUL.tar.gz

wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm
hmmpress dbCAN_sub.hmm

wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa
mv CAZyDB.07262023.fa CAZyDB.fa
diamond makedb --in CAZyDB.fa -d CAZy

wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt
mv dbCAN-HMMdb-V12.txt dbCAN.txt
hmmpress dbCAN.txt

wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa
diamond makedb --in tcdb.fa -d tcdb

wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm
hmmpress tf-1.hmm

wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm
hmmpress tf-2.hmm

wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm
hmmpress stp.hmm


dbcan_build --cpus 32 --db-dir /data/DBforDBcan/db --clean

cd /data/DBforDBcan/V13
wget https://bcb.unl.edu/dbCAN2/download/Databases/V13/CAZyDB.07142024.fa

wget https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt

wget https://bcb.unl.edu/dbCAN2/download/Databases/V13/hmmscan-parser.sh


** if you want to run dbCAN CAZyme annotation on your local linux computer, do the following:
** 1. download dbCAN-fam-HMMs.txt, hmmscan-parser.sh 
** 2. download HMMER 3.0 package [hmmer.org] and install it properly
** 3. format HMM db: hmmpress dbCAN-fam-HMMs.txt
** 4. run: hmmscan --domtblout yourfile.out.dm dbCAN-fam-HMMs.txt yourfile > yourfile.out
** 5. run: sh hmmscan-parser.sh yourfile.out.dm > yourfile.out.dm.ps (if alignment > 80aa, use E-value < 1e-5, otherwise use E-value < 1e-3; covered fraction of HMM > 0.3)
** 6. run: cat yourfile.out.dm.ps | awk '$5<1e-15&&$10>0.35' > yourfile.out.dm.ps.stringent (this allows you to get the same result as what is produced in our dbCAN2 webpage)
Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage
** About what E-value and Coverage cutoff thresholds you should use (in order to further parse yourfile.out.dm.ps file), we have done some evaluation analyses using arabidopsis, rice, Aspergillus nidulans FGSC A4, Saccharomyces cerevisiae S288c and Escherichia coli K-12 MG1655, Clostridium thermocellum ATCC 27405 and Anaerocellum thermophilum DSM 6725. Our suggestion is that for plants, use E-value < 1e-23 and coverage > 0.2; for bacteria, use E-value < 1e-18 and coverage > 0.35; and for fungi, use E-value < 1e-17 and coverage > 0.45.
** We have also performed evaluation for the five CAZyme classes separately, which suggests that the best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)
** On our dbCAN2 website, we use E-value < 1e-15 and coverage > 0.35, which is more stringent than the default ones in hmmscan-parser.sh


##### 
conda activate dbcan
run_dbcan --db_dir /data/DBforDBcan/db --tools all --out_dir /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/11_dbcan/output2 \
          --dia_eval 1e-5 --hmm_eval 1e-5 --tf_eval 1e-5 --stp_eval 1e-5 -evalue 1e-5 -hmmevalue 1e-5  \
          --out_pre nonviral_dbcan --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32 --stp_cpu 32 \
          /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.faa protein

###header for diamond
	1.	Query ID (qseqid): The ID of the query protein.
	2.	Subject ID (sseqid): The ID of the best matching CAZyme protein in the CAZy database.
	3.	Percentage of Identical Matches (pident): The percentage of identical matches in the alignment.
	4.	Alignment Length (length): The length of the alignment between the query and subject sequences.
	5.	Number of Mismatches (mismatch): The number of mismatches in the alignment.
	6.	Number of Gap Openings (gapopen): The number of gaps introduced in the alignment.
	7.	Start Position in Query (qstart): The start position of the alignment in the query sequence.
	8.	End Position in Query (qend): The end position of the alignment in the query sequence.
	9.	Start Position in Subject (sstart): The start position of the alignment in the subject sequence.
	10.	End Position in Subject (send): The end position of the alignment in the subject sequence.
	11.	E-value (evalue): The E-value of the alignment, indicating its statistical significance.
	12.	Bit Score (bitscore): A measure of the alignment’s quality, with higher scores indicating better matches.

##### Eggnog
conda activate eggnog
export EGGNOG_DATA_DIR=/data/DBforEggNog/
emapper.py -m diamond --cpu 48 --itype proteins --allow_overlaps none --report_orthologs \
               -i /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.faa \
               --output_dir /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/12_eggnog \
               --output EggNogRes --evalue 1e-05

##### COG
anvi-script-reformat-fasta /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.faa -o /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/14_COG/nonviral_orfs.simple.faa --simplify-names
anvi-gen-contigs-database -f /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/14_COG/nonviral_orfs.simple.faa -o  /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/14_COG/nonviral_orf.db -n phageAMG

conda activate mmseqs2
mmseqs createdb COGorg24.faa COG_mmseqs_db
mmseqs createindex COG_mmseqs_db tmp
mmseqs easy-search /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/5_sort/nonviral_orfs.faa /data/DBforMetaEuk/uniprot/uniprotKB /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/8_uniprotKB/uniprotKB.res /data/ForestVirome/DNAvirus/4_annotation/2_function/AMGsearch/8_uniprotKB/tmp -s 7




#### Gene expression
sed 's/#.*//' /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.protein.meta.trueVirus.fna > /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.gene.trueVirus.cleanheader.fna

conda activate seqkit
seqkit grep -f /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/genuine.phage.list /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.gene.trueVirus.cleanheader.fna \
            -o /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.gene.trueVirus.finalCleaned.fna

for f in $(cat /data/SongyiMetatranscriptome/SampleList2023) 
do
    n=${f}  
    mkdir -p /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}
done



bwa index /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.gene.trueVirus.finalCleaned.fna

for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)  
do
    n=${f}
    bwa mem -t 48 -B 0 /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.gene.trueVirus.finalCleaned.fna \
                  /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_1_trimmed.fastq.gz /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_2_trimmed.fastq.gz > \
                  /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sam
    samtools view -bS /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.bam
    samtools sort /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.bam -o /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sorted.bam
    samtools index /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sorted.bam
    samtools flagstat /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sorted.bam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.stat.txt
done



for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)  
do
    n=${f}
    bwa mem -t 48 /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.gene.trueVirus.finalCleaned.fna \
                  /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_1_trimmed.fastq.gz /data/SongyiMetatranscriptome/data2023/1_adapter/${n}_2_trimmed.fastq.gz > \
                  /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sam
    samtools view -bS /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.bam
    samtools sort /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.bam -o /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sorted.bam
    samtools index /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sorted.bam
    samtools flagstat /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sorted.bam > /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.stat.txt
done

#### Using salmon

for f in $(cat /data/SongyiMetatranscriptome/SampleList2023) 
do
    n=${f}  
    mkdir -p /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/strict/salmon
done

for f in $(cat /data/SongyiMetatranscriptome/SampleList2023)
do
    n=${f}
    salmon quant -t /data/ForestVirome/DNAvirus/4_annotation/2_function/0_gene/genes/repre.votu.gene.trueVirus.finalCleaned.fna -l ISR \
                 -a /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/bwa.strict.sorted.bam \
                 -o /data/ForestVirome/DNAvirus/4_annotation/2_function/7_abundance/RNA/${n}/strict/salmon \
                 -p 48 --meta
done


####### Prediction of life style of phages
mamba activate bacphlip



