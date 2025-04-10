## import libraries : phyloseq and microbiome

### Install required packages
##Phyloseq
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")
# BiocManager::install("metagenomeSeq") 
BiocManager::install("microbiome")
# 
# install.packages("devtools")
# library(devtools)  
# install_github("hallucigenia-sparsa/seqtime") 
# 
# remotes::install_github("vmikk/metagMisc")
#install.packages("writexl")

### Load required packages
library(dplyr)
library(forcats) 
library(metagenomeSeq)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(scales)
library(grid)
library(reshape2)
#library(seqtime)
library(agricolae)
library(RColorBrewer)
library(magrittr)
library(indicspecies)
library(Hmisc)
library(igraph)
library(qgraph)
library(randomForest)
library(multifunc)
library(FSA)
library(rcompanion)
library(seqinr)
library(metagMisc)
library(readxl)
library(writexl)
library(ggrepel)

# Set plotting theme
theme_set(theme_bw())


######## For the bacterial community ##########
### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
bac_phylo=import_biom("./Bacteria/asv_table_final.biom")
### merge with metadata
# Import sample metadata

## in metadata erase # (This step is essential)
map <- read.table(file = '../SampleMetadata/sample_metadata_for bacteria.tsv', sep = '\t', header = TRUE)
map <- sample_data(map)

head(map)
dim(map)
summary(map)
str(map)

summary(map)
colnames(map)
rownames(map)
nrow(map)

# Assign rownames to be Sample ID's
map$sampleID
rownames(map) <- map$sampleID
rownames(map)
dim(map)

# Merge biom data object with sample metadata + tree data(this is rooted!)
phy_tree = read_tree("./Bacteria/tree.nwk")
phy <- merge_phyloseq(bac_phylo, map, phy_tree)

class(phy)
phy   ## 7967 ASVs


## changing rank names
colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
phy  ## 7967 ASVs
sort(colSums(otu_table(phy)))



######## For the fungal community ##########
### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
fun_phylo=import_biom("./Fungi/asv_table_final.biom")
### merge with metadata
# Import sample metadata

## in metadata erase # (This step is essential)
f.map <- read.table(file = '../SampleMetadata/sample_metadata_for fungi.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$sampleID
rownames(f.map) <- f.map$sampleID
rownames(f.map)
dim(f.map)

# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("./Fungi/tree.nwk")
fun <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun)
fun   ## 2082 ASVs


## changing rank names
colnames(tax_table(fun)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
fun  ## 2082 ASVs
sort(colSums(otu_table(fun)))
