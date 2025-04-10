### Normalization
library(metagenomeSeq)
## Let's do normalization with CSS
## phyloseq to metagenomeSeq


#phy.clean? or phy.clean --> let's start from phy.clean


## Remove residual taxa that do not have any sequences
#Bacteria
sum(taxa_sums(bac.clean.ss.AR.AP) == 0)
taxa_sums(bac.clean.ss.AR.AP)

bac.clean.ss.AR.AP <- phyloseq::filter_taxa(bac.clean.ss.AR.AP, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(bac.clean.ss.AR.AP) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(bac.clean.ss.AR.AP))


## only keep samples over 0 read
(filt.sample <- sample_sums(bac.clean.ss.AR.AP) > 0)
sum(sample_sums(bac.clean.ss.AR.AP) <= 0)  ## 0 sample discarded
bac.clean.ss.AR.AP.f <- prune_samples(sample_sums(bac.clean.ss.AR.AP) > 0, bac.clean.ss.AR.AP)
bac.clean.ss.AR.AP.f  

colSums(otu_table(bac.clean.ss.AR.AP.f))
sort(sample_sums(bac.clean.ss.AR.AP.f))

## CODE for CSS normalization using preloaded data
bac.clean.filt <- bac.clean.ss.AR.AP.f
bac.clean.filt <- phyloseq::filter_taxa(bac.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples

met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)
# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean) 
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog.AR.AP <- bac.clean.filt
otu_table(bac.clean.nolog.AR.AP) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log.AR.AP <- bac.clean.filt
otu_table(bac.clean.log.AR.AP) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


### Fungi
fun.clean.ss2.AR.AP

## Remove residual taxa that do not have any sequences
#Bacteria
sum(taxa_sums(fun.clean.ss2.AR.AP) == 0)
taxa_sums(fun.clean.ss2.AR.AP)

fun.clean.ss2.AR.AP <- phyloseq::filter_taxa(fun.clean.ss2.AR.AP, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(fun.clean.ss2.AR.AP) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(fun.clean.ss2.AR.AP))


## only keep samples over 0 read
(filt.sample <- sample_sums(fun.clean.ss2.AR.AP) > 0)
sum(sample_sums(fun.clean.ss2.AR.AP) <= 0)  ## 0 sample discarded
fun.clean.ss2.AR.AP.f <- prune_samples(sample_sums(fun.clean.ss2.AR.AP) > 0, fun.clean.ss2.AR.AP)
fun.clean.ss2.AR.AP.f  

colSums(otu_table(fun.clean.ss2.AR.AP.f))
sort(sample_sums(fun.clean.ss2.AR.AP.f))

## CODE for CSS normalization using preloaded data
fun.clean.filt <- fun.clean.ss2.AR.AP.f
fun.clean.filt <- phyloseq::filter_taxa(fun.clean.filt, function(x) sum(x) != 0, TRUE)
# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
fun.clean.filt    ## use all samples

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)



# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.fun.clean) 
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## funk to the phyloseq file
fun.clean.nolog.AR.AP <- fun.clean.filt
otu_table(fun.clean.nolog.AR.AP) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log.AR.AP <- fun.clean.filt
otu_table(fun.clean.log.AR.AP) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)



#### Rarefying (input 1 for calculating alpha diversity)
set.seed(8046)

bac.rarefied<-rarefy_even_depth(bac.clean.ss.dna.f, sample.size = min(sample_sums(bac.clean.ss.dna.f)),
                  rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

fun.rarefied<-rarefy_even_depth(fun.clean.ss.f, sample.size = min(sample_sums(fun.clean.ss.f)),
                                rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


#### hellinger transformation (input 2 for calculating alpha diversity)
fun.hellinger.tab <- decostand(otu_table(fun.clean.ss2.AR.AP.f), method = "hellinger")
fun.hellinger <- fun.clean.ss2.AR.AP.f

otu_table(fun.hellinger)<- otu_table(data.frame(fun.hellinger.tab), taxa_are_rows = T)
