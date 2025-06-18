#### 1. Functional annotation results
phrog.res<-read.table("~/Desktop/SNU/KNU_Prof.HangilKim/ForestVirome/DNAvirus/FunctionalAnnotation/phrog/RevisedRes/phrogs.results.tsv", sep = '\t', header =F)
head(phrog.res)
nrow(phrog.res)
names(phrog.res) <- c("query","target","align.score","seq.identity","e.value","query.start","query.end","query.length","target.start","target.end","target.length")

phrog.res$align.length.q <- phrog.res$query.end-phrog.res$query.start +1
phrog.res$align.length.t <- phrog.res$target.end-phrog.res$target.start +1
phrog.res$align.portion.q <- phrog.res$align.length.q/phrog.res$query.length
phrog.res$align.portion.t <- phrog.res$align.length.t/phrog.res$target.length
matched.phrogs<-unique(phrog.res$query)

###### Load a reference table
phrog.annot.tab <- read.table("~/Desktop/SNU/KNU_Prof.HangilKim/ForestVirome/DNAvirus/FunctionalAnnotation/phrog/phrog_annot_v4.tsv", sep = '\t', header =T, quote="")
phrog.annot.tab$query <- paste0("phrog_",phrog.annot.tab$phrog)
head(phrog.annot.tab)
nrow(phrog.annot.tab)

nrow(subset(phrog.annot.tab, query %in% matched.phrogs))
length(unique(phrog.res$target)) #21077

initial.phrog.res<-merge(phrog.res,phrog.annot.tab, by ="query")
nrow(phrog.res) #5,780,199
nrow(initial.phrog.res) #5,780,199


#### Exclude results from the excluded vMAGs
#### Load GFF file to get gene list which needs to be removed from the phrog result
vMAG.gene.annot<-read.delim("~/Desktop/SNU/KNU_Prof.HangilKim/ForestVirome/DNAvirus/FunctionalAnnotation/viral_proteins/repre.votu.protein.meta.trueVirus.gff", sep ="\t", header = F)
head(vMAG.gene.annot)
vMAGs.exclude.rename <- gsub("||full","",vMAGs.exclude)
vMAGs.exclude.rename <- gsub("\\||","",vMAGs.exclude.rename)

vMAGs.exclude.rename<-c("TmD1IndAssemContig0359700","TmDContig00072623")

initial.phrog.res.wo.unknowns <- initial.phrog.res

for (i in as.character(vMAGs.exclude.rename)){
  initial.phrog.res.wo.unknowns <- initial.phrog.res.wo.unknowns[!(grepl(i,initial.phrog.res.wo.unknowns$target)),]
}

nrow(initial.phrog.res)
names(initial.phrog.res)
nrow(initial.phrog.res.wo.unknowns)
unique(initial.phrog.res.wo.unknowns$target) #61299
head(initial.phrog.res.wo.unknowns)
##### Get significant matching results
#### Threshold: alignment rate >= 70%, e.value < 0.001, seq.identity > 30%
sig.phrog.res<-subset(initial.phrog.res.wo.unknowns, align.portion.t >= 0.7 & e.value < 1e-05 & seq.identity >= 0.3 & align.score >=50)
nrow(sig.phrog.res) #2518544
length(unique(sig.phrog.res$target)) #61,285 proteins
write.table(sig.phrog.res,"./sig.PHROG.res.txt",sep= "\t",col.names = T,row.names = F)

##### Get the best results - more strict condition
library(dplyr)
sig.phrog.res_sorted <- sig.phrog.res[order(sig.phrog.res$e.value, -sig.phrog.res$align.score), ]

sig.phrog.res.best.hit <- sig.phrog.res_sorted %>%
  group_by(target) %>%
  slice(1) %>%
  ungroup()

print(sig.phrog.res.best.hit)
writexl::write_xlsx(sig.phrog.res.best.hit,"./PHROG_besthit.xlsx")

####### Add KofamScan, VOGs, CheckV, pVOGs results
kofam.AMG.candidates<-read.table("./AMGannotation/kofam/kofamscan_results_allNonviralORFs.txt",sep="\t",header=F, stringsAsFactors = F,quote = "")
head(kofam.AMG.candidates)
nrow(kofam.AMG.candidates)
colnames(kofam.AMG.candidates) <- c("gene", "KO_ID", "threshold", "Bit_Score", "E_value", "KO_Description")

nonviralORF.list <- read.table("./AMGannotation/nonviral_orfs.txt",sep="\t",quote="",header =F)
head(nonviralORF.list)
names(nonviralORF.list)[1]<-"ORF"

kofam.AMG.candidates.nonviral <-subset(kofam.AMG.candidates , gene %in% nonviralORF.list$ORF)
kofam.AMG.candidates.sig <- kofam.AMG.candidates.nonviral[kofam.AMG.candidates.nonviral$Bit_Score>=50&kofam.AMG.candidates.nonviral$E_value<1e-05,]
nrow(kofam.AMG.candidates.sig) #1,932 hits
unique(kofam.AMG.candidates.sig$gene) #626 genes
head(kofam.AMG.candidates.sig)
kofam.AMG.candidates.sig$Bit_Score <- as.numeric(kofam.AMG.candidates.sig$Bit_Score)
kofam.AMG.candidates.sig_sorted <- kofam.AMG.candidates.sig[order(kofam.AMG.candidates.sig$E_value, -kofam.AMG.candidates.sig$Bit_Score), ]
head(kofam.AMG.candidates.sig_sorted)


kofam.AMG.candidates.sig.best <- kofam.AMG.candidates.sig_sorted  %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup()
writexl::write_xlsx(kofam.AMG.candidates.sig,"./sig.kofam_hit.xlsx")
writexl::write_xlsx(kofam.AMG.candidates.sig.best,"./sig.kofam_besthit.xlsx")

### VOGs
VOG.AMG.candidates<-read.table("./AMGannotation/VOG/hmm_search_res_filtered.txt",sep="\t",header=F, stringsAsFactors = F,quote = "")
head(VOG.AMG.candidates)
nrow(VOG.AMG.candidates)
colnames(VOG.AMG.candidates) <-c("Target_Name", "Query_Name", "Full_Evalue", "Full_Score", "Full_Bias",
                                 "Domain_Evalue", "Domain_Score", "Domain_Bias", "Best_HMM_Domain",
                                 "Total_Domains", "Domain_Start", "Domain_End", "Target_Start",
                                 "Target_End", "Envelope_Start", "Envelope_End", "Accuracy")

length(unique(VOG.AMG.candidates$Target_Name)) #12,196

VOG.ref.tab <- read.delim(gzfile("./AMGannotation/VOG/vog.annotations.tsv.gz"),sep="\t",quote="")
VOG.ref.tab<- VOG.ref.tab[c(1,4,5)]
names(VOG.ref.tab)[1] <-"Query_Name"
VOG.AMG.candidates <- merge(VOG.ref.tab,VOG.AMG.candidates,by="Query_Name")


VOG.AMG.candidates.nonviral <-subset(VOG.AMG.candidates , Target_Name %in% nonviralORF.list$ORF)
unique(VOG.AMG.candidates.nonviral$ConsensusFunctionalDescription)
length(unique(VOG.AMG.candidates.nonviral$Target_Name)) #2,289 genes
head(VOG.AMG.candidates.nonviral)

VOG.AMG.candidates.nonviral.sig <- VOG.AMG.candidates.nonviral[VOG.AMG.candidates.nonviral$Full_Score>=50&VOG.AMG.candidates.nonviral$Full_Evalue<1e-05,]
VOG.AMG.candidates.nonviral_sorted <- VOG.AMG.candidates.nonviral.sig [order(VOG.AMG.candidates.nonviral.sig$Full_Evalue, -VOG.AMG.candidates.nonviral.sig$Full_Score), ]
head(VOG.AMG.candidates.nonviral_sorted)

VOG.AMG.candidates.nonviral.best <- VOG.AMG.candidates.nonviral_sorted  %>%
  group_by(Target_Name) %>%
  slice(1) %>%
  ungroup()

writexl::write_xlsx(VOG.AMG.candidates.nonviral.sig,"./sig.VOG_hit.xlsx")
writexl::write_xlsx(VOG.AMG.candidates.nonviral.best,"./sig.VOG_besthit.xlsx")


#### Comparison of annotation results of three databases
library(dplyr)
kofam.AMG.candidates.sig.best
VOG.AMG.candidates.nonviral.best

head(kofam.AMG.candidates.sig)
kofam.AMG.candidates.sig.best$KO_Description <- gsub("\"","",kofam.AMG.candidates.sig.best$KO_Description)
KOfam.summary<-kofam.AMG.candidates.sig.best%>% group_by(gene) %>% summarize(KOfam = paste0(KO_Description, collapse = ";"))

head(VOG.AMG.candidates.nonviral.best)
VOG.summary<-VOG.AMG.candidates.nonviral.best%>% group_by(Target_Name) %>% summarize(VOG = paste0(ConsensusFunctionalDescription, collapse = ";"))
names(VOG.summary)[1] <- "gene"

head(sig.phrog.res)
sig.phrog.res$annot[is.na(sig.phrog.res$annot)] <- "unknown"
sig.phrog.res.nonviral <- subset(sig.phrog.res, target %in% nonviralORF.list$ORF)
sig.phrog.res.nonviral <- sig.phrog.res.nonviral[order(sig.phrog.res.nonviral$e.value, -sig.phrog.res.nonviral$align.score), ]
sig.phrog.res.nonviral.best <- sig.phrog.res.nonviral %>%
  group_by(target) %>%
  slice(1) %>%
  ungroup()

PHROGs.summary <- VOG.AMG.candidates.nonviral.best %>% group_by(target) %>% summarize(PHROGs = paste0(annot, collapse = ";"))
names(PHROGs.summary)[1] <- "gene"
nrow(PHROGs.summary)
nonviral.annot.res<-merge(KOfam.summary,VOG.summary, by = "gene",all=T)
nrow(nonviral.annot.res)
nonviral.annot.res<-merge(nonviral.annot.res,PHROGs.summary, by = "gene",all=T)
nrow(nonviral.annot.res) #49207

###### Other DBs
###uniprotKB###
initial.uniprot.nonviral <- read.delim("./AMGannotation/uniprotKB/uniprotKB.res",
                                       sep="\t",quote="",header =F)
nrow(initial.uniprot.nonviral) # 1,643,867
head(initial.uniprot.nonviral)
names(initial.uniprot.nonviral) <- c("query","target","seqIden","align.length","mismatch",
                                     "gapOpening","qStart","qEnd","hStart","hEnd","e.val","bitScore")
#initial.uniprot.nonviral <- merge(initial.uniprot.nonviral,nonviral.AA.length.tab,by="query")
#nrow(initial.uniprot.nonviral) #1643867
#initial.uniprot.nonviral$align.portion.q <- initial.uniprot.nonviral$align.length/initial.uniprot.nonviral$qTotalLength
#head(initial.uniprot.nonviral)

sig.uniprot.nonviral.res<-subset(initial.uniprot.nonviral, e.val < 1e-05 & seqIden >= 0.3 & bitScore >=50)
nrow(sig.uniprot.nonviral.res) # 1,152,590
length(unique(sig.uniprot.nonviral.res$query)) #18,269 non-viral proteins
write.table(sig.uniprot.nonviral.res,"./sig.uniprot.annotation.res.txt",
            sep="\t",row.names =F,col.names = T)

sig.uniprot.nonviral.res_sorted <- sig.uniprot.nonviral.res[order(sig.uniprot.nonviral.res$e.val, -sig.uniprot.nonviral.res$bitScore), ]
sig.uniprot.nonviral.res.best.hit <- sig.uniprot.nonviral.res_sorted %>%
  group_by(query) %>%
  slice(1) %>%
  ungroup()

length(unique(sig.uniprot.nonviral.res.best.hit$target)) #13848

head(sig.uniprot.nonviral.res.best.hit)

uniprotKB.annot.tab<-read.delim("./AMGannotation/uniprotKB/idmapping_2025_02_10.tsv",sep="\t", quote="", fill=T,header =T)
head(uniprotKB.annot.tab)
uniprotKB.annot.tab <-uniprotKB.annot.tab[-c(2)]
names(uniprotKB.annot.tab)[1] <- "target"
sig.uniprotKB.nonviral.res.annot <- merge(sig.uniprot.nonviral.res.best.hit,uniprotKB.annot.tab, by ="target")
nrow(sig.uniprotKB.nonviral.res.annot)
nrow(sig.uniprot.nonviral.res.best.hit)
write.table(sig.uniprot.nonviral.res.best.hit,"./sig.uniprot.annotation.res.besthit.txt",
            sep="\t",row.names =F,col.names = T)

##### pfamA
initial.pfamA.nonviral <- read.delim("./AMGannotation/pfamA/pfamA.res",
                                     sep="\t",quote="",header =F)
nrow(initial.pfamA.nonviral) # 3670
head(initial.pfamA.nonviral)
names(initial.pfamA.nonviral) <- c("query","target","seqIden","align.length","mismatch",
                                   "gapOpening","qStart","qEnd","hStart","hEnd","e.val","bitScore")

sig.pfamA.nonviral.res<-subset(initial.pfamA.nonviral, e.val < 1e-05 & seqIden >= 0.3 & bitScore >=50)
nrow(sig.pfamA.nonviral.res) # 254
length(unique(sig.pfamA.nonviral.res$query)) #249 non-viral proteins

pfamA.annot.tab<-read.delim(gzfile("./AMGannotation/pfamA/Pfam-A.clans.tsv.gz"),sep="\t", quote="", fill=T,header =F)
head(pfamA.annot.tab)
names(pfamA.annot.tab) <- c("pfamID","ClanID","Clan","Short.Name","Description")

head(sig.pfamA.nonviral.res)
sig.pfamA.nonviral.res$pfamID <- gsub("\\.+\\w+","",sig.pfamA.nonviral.res$target)
sig.pfamA.nonviral.res$target
sig.pfamA.nonviral.res.annot <- merge(sig.pfamA.nonviral.res,pfamA.annot.tab, by ="pfamID")
nrow(sig.pfamA.nonviral.res.annot)
nrow(sig.pfamA.nonviral.res)

write.table(sig.pfamA.nonviral.res.annot,"./sig.pfamA.annotation.res.txt",
            sep="\t",row.names =F,col.names = T)

sig.pfamA.nonviral.res_sorted <- sig.pfamA.nonviral.res[order(sig.pfamA.nonviral.res$e.val, -sig.pfamA.nonviral.res$bitScore), ]
sig.pfamA.nonviral.res.best.hit <- sig.pfamA.nonviral.res_sorted %>%
  group_by(query) %>%
  slice(1) %>%
  ungroup()
sig.pfamA.nonviral.res.annot.best <- merge(sig.pfamA.nonviral.res.best.hit,pfamA.annot.tab, by ="pfamID")
write.table(sig.pfamA.nonviral.res.annot.best,"./sig.pfamA.annotation.res.best.txt",
            sep="\t",row.names =F,col.names = T)

##### dbcan (CAZyme)
overview.dbcan.nonviral <- read.delim("./AMGannotation/dbcan/nonviral_dbcanoverview.txt",
                                      sep="\t",quote="",header =T)
nrow(overview.dbcan.nonviral) # 784

initial.dbcandiamond.nonviral <- read.delim("./AMGannotation/dbcan/nonviral_dbcandiamond.out",
                                            sep="\t",quote="",header =T)

initial.hmmer.nonviral <- read.delim("./AMGannotation/dbcan/nonviral_dbcanhmmer.out",
                                     sep="\t",quote="",header =T)
initial.dbcansub.nonviral <- read.delim("./AMGannotation/dbcan/nonviral_dbcandbcan-sub.hmm.out",
                                        sep="\t",quote="",header =F)

initial.dbcansub.nonviral<- initial.dbcansub.nonviral[-c(14)]
colnames(initial.dbcansub.nonviral) <- initial.dbcansub.nonviral[1,]
head(initial.dbcansub.nonviral)
initial.dbcansub.nonviral <-initial.dbcansub.nonviral[-1,]
nrow(initial.dbcansub.nonviral) #217
head(initial.dbcansub.nonviral)

initial.dbcandiamond.nonviral$seqIden <- initial.dbcandiamond.nonviral$X..Identical/100

sig.dbcandiamond.nonviral.res<-subset(initial.dbcandiamond.nonviral, E.Value < 1e-05 & seqIden >= 0.3 & Bit.Score >=50)

head(overview.dbcan.nonviral)
names(overview.dbcan.nonviral)[1] <-"gene"
nonviral.annot.res.2<-merge(nonviral.annot.res,overview.dbcan.nonviral,by="gene",all=T)
nrow(nonviral.annot.res.2) #49207


#### Eggnog
initial.eggnog.nonviral <- read.delim("./AMGannotation/eggnog/EggNogRes.emapper.annotations",
                                      sep="\t",quote="",header =T)
head(initial.eggnog.nonviral)
nrow(initial.eggnog.nonviral)
sig.eggnog.nonviral.res <- subset(initial.eggnog.nonviral,score>=50)
nrow(sig.eggnog.nonviral.res) #8,635
unique(sig.eggnog.nonviral.res$query)
names(sig.eggnog.nonviral.res)[1] <- "gene"

head(sig.eggnog.nonviral.res)
write.table(sig.eggnog.nonviral.res,"./sig.eggnog.annotation.res.txt",
            sep="\t",row.names =F,col.names = T)

simple.eggnog.tab <-sig.eggnog.nonviral.res[c(1,2,5,8,9)]
names(simple.eggnog.tab)
simple.eggnog.tab$Description[which(simple.eggnog.tab$Description == "-" & simple.eggnog.tab$eggNOG_OGs != "-")] <- "unknown"



sig.eggnog.nonviral.res
names(sig.eggnog.nonviral.res)[1] <- "gene"


nonviral.annot.res.3 <- nonviral.annot.res.2[-c(5,9)] 
nonviral.annot.res.3<-merge(nonviral.annot.res.3,simple.eggnog.tab,by ="gene",all=T)
nrow(nonviral.annot.res.3) #49220

sig.pfamA.nonviral.res.annot.best
names(sig.pfamA.nonviral.res.annot.best)
simple.pfamA.tab <- sig.pfamA.nonviral.res.annot.best[c(2,17)]
names(simple.pfamA.tab) <- c("gene","pfamA")

nonviral.annot.res.3<-merge(nonviral.annot.res.3,simple.pfamA.tab,by ="gene",all=T)

sig.uniprotKB.nonviral.res.annot
names(sig.uniprotKB.nonviral.res.annot)
simple.uniprotKB.tab <- sig.uniprotKB.nonviral.res.annot[c(2,15)]
names(simple.uniprotKB.tab) <- c("gene","uniprotKB")

nonviral.annot.res.3<-merge(nonviral.annot.res.3,simple.uniprotKB.tab,by ="gene",all=T)
nrow(nonviral.annot.res.3) #49235
write.table(nonviral.annot.res.3,"./merge_annotations_reanalyzed.txt",col.names = T,row.names = F,quote = F,sep="\t")
writexl::write_xlsx(nonviral.annot.res.3,"./merge_annotations.xlsx")



######## 
initial.merged.annotation <- readxl::read_xlsx("./merge_annotations_edited.xlsx")
nrow(initial.merged.annotation) #49,235
head(initial.merged.annotation)
initial.merged.annotation <-data.frame(initial.merged.annotation)

initial.merged.annotation$uniprotKB[grepl("Uncharacterized protein",initial.merged.annotation$uniprotKB)] <-"unknown"
initial.merged.annotation$VOG[grepl("hypothetical protein",initial.merged.annotation$VOG)] <-"unknown"


nrow(initial.merged.annotation)

initial.merged.annotation[is.na(initial.merged.annotation)] <- "NotAnnotated"
head(initial.merged.annotation)
merged.annotation.melt<-reshape2::melt(initial.merged.annotation,id.var="gene")
head(merged.annotation.melt)
names(merged.annotation.melt)[2] <- "DB"
names(merged.annotation.melt)[3] <- "Description"

merged.annotation.melt.annotated <- subset(merged.annotation.melt, !(Description%in%c("unknown","NotAnnotated")))
nrow(merged.annotation.melt.annotated)
unique(merged.annotation.melt.annotated$gene) #38,891

merged.annotation.melt.unknown <- subset(merged.annotation.melt, Description%in%c("unknown","NotAnnotated"))
unknown.gene.list<-unique(merged.annotation.melt.unknown$gene)
true.unknown.genes<-unknown.gene.list[-which(unknown.gene.list%in%unique(merged.annotation.melt.annotated$gene))] #10,344


#### Find sole annotations
annotated.DB.count <- merged.annotation.melt.annotated
annotated.DB.count$Count <- 1
annotated.DB.count <- annotated.DB.count %>% group_by(gene) %>% reframe(Count = sum(Count))
sole.annotated.gene <- annotated.DB.count$gene[which(annotated.DB.count$Count == 1)]
sole.annotated.gene.tab<-merged.annotation.melt.annotated[merged.annotation.melt.annotated$gene%in% sole.annotated.gene ,]
nrow(sole.annotated.gene.tab) #29429
head(sole.annotated.gene.tab)

unknown.gene.tab <- data.frame(gene=true.unknown.genes,DB="NotAvailable",Description="Hypothetical protein")
consensus.annotation.tab<-rbind(sole.annotated.gene.tab,unknown.gene.tab)
unique(consensus.annotation.tab$DB)
nrow(consensus.annotation.tab) # 39773
nrow(initial.merged.annotation) #49235
find.consensus.tab<-initial.merged.annotation[!(initial.merged.annotation$gene%in%consensus.annotation.tab$gene), ]
writexl::write_xlsx(find.consensus.tab,"./find.consensus.tab.xlsx")


merged.annotation.melt.annotated


##### COG
init.COG.res<-read.delim("./AMGannotation/COG/COGres",sep="\t",header =F,fill=T, quote="")
head(init.COG.res)
length(unique(init.COG.res$V1)) #11,389

nonviralORF.list
init.COG.res.nonviral <- subset(init.COG.res,V1%in%nonviralORF.list$ORF)
length(unique(init.COG.res.nonviral$V1))

COG.res.nonviral.sig <- init.COG.res.nonviral[init.COG.res.nonviral$V11<1e-05&init.COG.res.nonviral$V12>=50,]
length(unique(COG.res.nonviral.sig$V1)) #9572

COG.res.nonviral.sig.sorted <- COG.res.nonviral.sig[order(COG.res.nonviral.sig$V11, -COG.res.nonviral.sig$V12),]
COG.res.nonviral.sig.sorted.tophit<- COG.res.nonviral.sig.sorted%>%
  group_by(V1) %>%
  slice(1) %>%
  ungroup()

COG.res.nonviral.sig.sorted.tophit2 <- COG.res.nonviral.sig.sorted.tophit[c(1,2)]
names(COG.res.nonviral.sig.sorted.tophit2) <- c("gene","COGtarget")
COG.res.nonviral.sig.sorted.tophit2$COGtarget
COG.gene.annotation<-read.csv("./AMGannotation/COG/cog-24.cog.csv",header=F)
head(COG.gene.annotation)
COG.gene.annotation.sub <- unique(COG.gene.annotation[c(3,7,8)])
nrow(COG.gene.annotation.sub)
head(COG.gene.annotation.sub)
names(COG.gene.annotation.sub) <- c("COGtarget","COG_1","COG_2")

unique(COG.gene.annotation.sub$COGtarget)
length(COG.gene.annotation$V3)
intersect(unique(COG.res.nonviral.sig.sorted.tophit2$COGtarget),COG.gene.annotation.sub$COGtarget) #1648
COG.res.nonviral.sig.sorted.tophit3 <- merge(COG.res.nonviral.sig.sorted.tophit2,COG.gene.annotation.sub,by="COGtarget")
nrow(COG.res.nonviral.sig.sorted.tophit3)


###### Obtain a table with identity, bit score, and e values, coverage
sig.pfamA.nonviral.res.best.hit
sig.uniprot.nonviral.res.best.hit
VOG.AMG.candidates.nonviral.best
sig.eggnog.nonviral.res
VOG.AMG.candidates.nonviral.best
kofam.AMG.candidates.sig.best

######### check if the data have coverage and aligned length info
head(sig.pfamA.nonviral.res.best.hit)
sig.pfamA.nonviral.res.best.hit$align.length.q <- sig.pfamA.nonviral.res.best.hit$qEnd-sig.pfamA.nonviral.res.best.hit$qStart +1
head(sig.uniprot.nonviral.res.best.hit)
sig.uniprot.nonviral.res.best.hit$align.length.q <- sig.uniprot.nonviral.res.best.hit$qEnd-sig.uniprot.nonviral.res.best.hit$qStart +1
head(VOG.AMG.candidates.nonviral.best)
VOG.AMG.candidates.nonviral.best$align.length.q <- VOG.AMG.candidates.nonviral.best$Target_End-VOG.AMG.candidates.nonviral.best$Target_Start +1

head(kofam.AMG.candidates.sig.best)
kofam.AMG.candidates.sig.best$align.length.q <- "-"

initial.eggnog.nonviral.hit <- read.delim("./AMGannotation/eggnog/EggNogRes.emapper.hits",
                                          sep="\t",quote="",header =F)
head(initial.eggnog.nonviral.hit)
nrow(initial.eggnog.nonviral.hit)
sig.eggnog.nonviral.hit.res <- subset(initial.eggnog.nonviral.hit,V1 %in% sig.eggnog.nonviral.res$query & V2 %in% sig.eggnog.nonviral.res$seed_ortholog & V8 >= 50 & V11 < 1e-05)
sig.eggnog.nonviral.hit.res <- subset(initial.eggnog.nonviral.hit,V1 %in% sig.eggnog.nonviral.res$gene)
head(sig.eggnog.nonviral.hit.res)
colnames(sig.eggnog.nonviral.hit.res ) <- c("qseqid",	"sseqid",	"pident",	"length",	"mismatch",	"gapopen",	"qstart",	"qend",	
                                            "sstart",	"send",	"evalue",	"bitscore",	"qcovhsp",	"scovhsp")

sig.eggnog.res.best.hit.info <-sig.eggnog.nonviral.hit.res  %>%
  group_by(qseqid) %>% slice(which.min(evalue))

write.table(sig.eggnog.res.best.hit.info,"./sig.eggnog.res.best.hit.info.txt",sep="\t",col.names = T,row.names = F)

sig.eggnog.res.best.hit.info$align.length.q <- sig.eggnog.res.best.hit.info$qend-sig.eggnog.res.best.hit.info$qstart+1

head(sig.eggnog.nonviral.res)


head(sig.phrog.res.nonviral.best.hit)
VOG.AMG.candidates.nonviral.best$Domain_Start
VOG.AMG.candidates.nonviral.best$Domain_End
cand.consensus.pfamA <- sig.pfamA.nonviral.res.best.hit[c(1,3,12,11)]
cand.consensus.uniprot <- sig.uniprot.nonviral.res.best.hit[c(1,3,12,11)]
cand.consensus.phrog <- sig.phrog.res.nonviral.best.hit[c(2,4,3,5)]
cand.consensus.eggnog <- sig.eggnog.res.best.hit.info[c(1,3,12,11)]
cand.consensus.VOG <- VOG.AMG.candidates.nonviral.best[c(4,6,5)]
cand.consensus.VOG$seqIden <- "-"
cand.consensus.VOG <- cand.consensus.VOG[c(1,4,2,3)]
cand.consensus.kofam <- kofam.AMG.candidates.sig.best[c(1,4,5)]
cand.consensus.kofam$seqIden <- "-"
cand.consensus.kofam <- cand.consensus.kofam[c(1,4,2,3)]

head(sig.eggnog.res.best.hit.info)


####### From dbcan results
initial.dbcandiamond.nonviral 
initial.hmmer.nonviral
initial.dbcansub.nonviral

head(initial.dbcansub.nonviral)
initial.dbcandiamond.nonviral.targetgene<-subset(initial.dbcandiamond.nonviral,Gene.ID %in%find.consensus.tab$gene )
initial.hmmer.nonviral.targetgene<-subset(initial.hmmer.nonviral,Gene.ID %in%find.consensus.tab$gene )
initial.dbcansub.nonviral.targetgene<-subset(initial.dbcansub.nonviral,`Gene ID` %in%find.consensus.tab$gene )

head(initial.dbcandiamond.nonviral.targetgene)
cand.consensus.dbcan.diamond <- initial.dbcandiamond.nonviral.targetgene[c(1,13,12,11)]
cand.consensus.dbcan.diamond$DB <- "dbcan_DIAMOND"
head(cand.consensus.dbcan.diamond)

head(initial.hmmer.nonviral.targetgene)
cand.consensus.dbcan.hmmer <- initial.hmmer.nonviral.targetgene[c(3,5)]
cand.consensus.dbcan.hmmer$DB <- "dbcan_HMMER"
cand.consensus.dbcan.hmmer$seqIden <- "-"
cand.consensus.dbcan.hmmer$Bit.Score <- "-"
head(cand.consensus.dbcan.hmmer)
cand.consensus.dbcan.hmmer<-cand.consensus.dbcan.hmmer[c(1,5,4,2,3)]

head(initial.dbcansub.nonviral.targetgene)
cand.consensus.dbcan.dbcansub <- initial.dbcansub.nonviral.targetgene[c(6,8)]
cand.consensus.dbcan.dbcansub$DB <- "dbcan_sub"
cand.consensus.dbcan.dbcansub$seqIden <- "-"
cand.consensus.dbcan.dbcansub$Bit.Score <- "-"
head(cand.consensus.dbcan.dbcansub)
cand.consensus.dbcan.dbcansub<-cand.consensus.dbcan.dbcansub[c(1,5,4,2,3)]

names(cand.consensus.pfamA) <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.uniprot) <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.phrog) <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.eggnog) <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.VOG) <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.kofam) <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.dbcan.diamond)[1:4] <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.dbcan.hmmer)[1:4] <-c("query" ,   "seqIden" , "bitScore" ,"e.val")
names(cand.consensus.dbcan.dbcansub)[1:4] <-c("query" ,   "seqIden" , "bitScore" ,"e.val")

cand.consensus.pfamA$DB <- "pfamA"
cand.consensus.uniprot$DB <- "uniprotKB"
cand.consensus.phrog$DB <- "PHROGs"
cand.consensus.eggnog$DB <- "Eggnog"
cand.consensus.VOG$DB <- "VOG"
cand.consensus.kofam$DB <- "KOfam"
names(cand.consensus.phrog)
cand.consensus.eggnog$seqIden <-cand.consensus.eggnog$seqIden/100

min(cand.consensus.eggnog$bitScore)
cand.consensus.eggnog <- subset(cand.consensus.eggnog, bitScore>= 50)

cand.consensus.info.tab<-rbind(cand.consensus.pfamA,cand.consensus.uniprot,cand.consensus.phrog,
      cand.consensus.eggnog,cand.consensus.VOG,cand.consensus.kofam,cand.consensus.dbcan.diamond,
      cand.consensus.dbcan.hmmer,cand.consensus.dbcan.dbcansub)

cand.consensus.info.tab$seqIden <- as.numeric(cand.consensus.info.tab$seqIden)
cand.consensus.info.tab$bitScore<- as.numeric(as.character(cand.consensus.info.tab$bitScore))

cand.consensus.info.tab.sub <- cand.consensus.info.tab[cand.consensus.info.tab$query %in% find.consensus.tab$gene,]
nrow(cand.consensus.info.tab)
nrow(cand.consensus.info.tab.sub)
##### Remove unknowns
head(cand.consensus.phrog)
phrog.unknown<-find.consensus.tab[find.consensus.tab$PHROGs%in%c("unknown","unknown function"),]$gene
VOG.unknown<-find.consensus.tab[find.consensus.tab$VOG=="unknown",]$gene
dbcan.hmmer.unknown<-find.consensus.tab[find.consensus.tab$dbcan_HMMER=="-",]$gene
dbcan.diamond.unknown<-find.consensus.tab[find.consensus.tab$dbcan_DIAMOND=="-",]$gene
dbcan.sub.unknown<-find.consensus.tab[find.consensus.tab$dbcan_sub=="-",]$gene
eggnog.unknown<-find.consensus.tab[find.consensus.tab$Eggnog=="unknown",]$gene
pfamA.unknown<-find.consensus.tab[find.consensus.tab$pfamA=="unknown",]$gene
uniprot.unknown<-find.consensus.tab[find.consensus.tab$uniprotKB=="unknown",]$gene

min(cand.consensus.info.tab.sub[cand.consensus.info.tab.sub$DB=="Eggnog",]$bitScore)
cand.consensus.info.tab.sub2 <- subset(cand.consensus.info.tab.sub, !(query%in% phrog.unknown & DB=="PHROGs"))
cand.consensus.info.tab.sub2 <- subset(cand.consensus.info.tab.sub2, !(query%in% VOG.unknown & DB=="VOG"))
cand.consensus.info.tab.sub2 <- subset(cand.consensus.info.tab.sub2, !(query%in% dbcan.hmmer.unknown & DB=="dbcan_HMMER"))
cand.consensus.info.tab.sub2 <- subset(cand.consensus.info.tab.sub2, !(query%in% dbcan.sub.unknown & DB=="dbcan_sub"))
cand.consensus.info.tab.sub2 <- subset(cand.consensus.info.tab.sub2, !(query%in% dbcan.diamond.unknown & DB=="dbcan_DIAMOND"))
cand.consensus.info.tab.sub2 <- subset(cand.consensus.info.tab.sub2, !(query%in% eggnog.unknown & DB=="Eggnog"))
cand.consensus.info.tab.sub2 <- subset(cand.consensus.info.tab.sub2, !(query%in% uniprot.unknown & DB=="uniprotKB"))



cand.consensus.info.tab.sub.sorted <- cand.consensus.info.tab.sub2[order(-cand.consensus.info.tab.sub2$seqIden, -cand.consensus.info.tab.sub2$bitScore),]
cand.consensus.info.tab.sorted.target <- cand.consensus.info.tab.sub.sorted%>%
  group_by(query) %>%
  slice(1) %>%
  ungroup()
nrow(cand.consensus.info.tab.sorted.target) #9462
nrow(find.consensus.tab)
cand.consensus.info.tab.sorted.target[cand.consensus.info.tab.sorted.target$DB=="Eggnog",]$bitScore
unique(cand.consensus.info.tab.sorted.target$DB)

writexl::write_xlsx(cand.consensus.info.tab.sorted.target,"./target_sorted_by_seqiden_bitscore_re4.xlsx")

##### Add descriptions
names(cand.consensus.info.tab.sorted.target)[1] <- "gene"
cand.consensus.info.tab.sorted.target.info<-merge(cand.consensus.info.tab.sorted.target,merged.annotation.melt.annotated.sub, by =c("gene"="gene","DB"="DB"))
head(merged.annotation.melt.annotated.sub)
nrow(merge(cand.consensus.info.tab.sorted.target,merged.annotation.melt.annotated.sub, by =c("gene"="gene","DB"="DB")))


cand.consensus.info.tab.sub2[cand.consensus.info.tab.sub2$query == "TmD1IndAssemContig0568181||full|provirus_1_33291_27",]
sig.eggnog.nonviral.hit.res[sig.eggnog.nonviral.hit.res$qseqid== "TmDContig01126438||full_5",]
sig.uniprot.nonviral.res.best.hit[sig.uniprot.nonviral.res.best.hit$query == "TmDContig01126438||full_5",]
#####
names(cand.consensus.info.tab.sorted.target.info)
consensus.tab.info <-cand.consensus.info.tab.sorted.target.info[c(1,4,5,6,2)]
names(consensus.tab.info)[4] <- "ConsensusAnnotation"

find.consensus.tab2<-merge(find.consensus.tab,consensus.tab.info, by = "gene")
writexl::write_xlsx(find.consensus.tab2,"./consensus_nonviral_annotation_check2.xlsx")

find.consensus.tab2[find.consensus.tab2$DB == "Eggnog",]$bitScore
TmD3IndAssemContig0606142||full_6

sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="TmD3IndAssemContig0606142||full_6",]

sig.phrog.res.nonviral.best.hit[sig.phrog.res.nonviral.best.hit$target=="TmD2IndAssemContig0679135||full_10",]
sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="TmD2IndAssemContig0679135||full_10",]
sig.pfamA.nonviral.res.annot.best[sig.pfamA.nonviral.res.annot.best$query=="TmD1IndAssemContig0568181||full|provirus_1_33291_36",]
VOG.AMG.candidates.nonviral.best[VOG.AMG.candidates.nonviral.best$Target_Name=="TmD1IndAssemContig0568181||full|provirus_1_33291_36",]
sig.eggnog.nonviral.hit.res[sig.eggnog.nonviral.hit.res$qseqid== "TmD1IndAssemContig0568181||full|provirus_1_33291_36",]
sig.dbcandiamond.nonviral.res[sig.dbcandiamond.nonviral.res$Gene.ID== "TmD1IndAssemContig0568181||full|provirus_1_33291_36",]

initial.hmmer.nonviral[initial.hmmer.nonviral$Gene.ID=="TmMContig00748229||full_2",]


initial.uniprot.nonviral[initial.uniprot.nonviral$query=="Y23TmD5IndAssemContig00121940||full_42",]
sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="Y23TmD3IndAssemContig00444329||full_33",]
sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="TmDContig00499598||full_4",]
sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="TmDContig01120711||full_3",]
sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="TmDContig01174769||full_13",]
sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="TmMContig01060093||full_3",]
sig.uniprotKB.nonviral.res.annot[sig.uniprotKB.nonviral.res.annot$query=="TmMContig01060093||full_8",]

nrow(consensus.annotation.tab)
head(consensus.annotation.tab)

writexl::write_xlsx(consensus.annotation.tab,"./sole and unknowns.xlsx")

####### Manual categorization
multi.annotated.tab <- readxl::read_xlsx("./consensus_nonviral_annotation_edited.xlsx")
nrow(multi.annotated.tab )
multi.annotated.tab[is.na(multi.annotated.tab$ConsensusAnnotation),]$gene
names(multi.annotated.tab)
multi.annotated.tab.ref <- multi.annotated.tab[c(11,12,13)]
multi.annotated.tab.ref  <- unique(multi.annotated.tab.ref )
nrow(multi.annotated.tab.ref)
names(multi.annotated.tab.ref) <- c("Description","DB","Category")
multi.annotated.tab.ref[is.na(multi.annotated.tab.ref$Description),]
head(consensus.annotation.tab)
consensus.annotation.tab2 <- consensus.annotation.tab
consensus.annotation.tab2$Category <- "0"

kofam.manual.annot<-multi.annotated.tab.ref$Description[which(multi.annotated.tab.ref$DB == "KOfam")]
for (i in kofam.manual.annot){
  consensus.annotation.tab2$Category[which(consensus.annotation.tab2$DB == "KOfam"&
                                    consensus.annotation.tab2$Description == i)] <- multi.annotated.tab.ref$Category[which(multi.annotated.tab.ref$DB=="KOfam"&
                                                                                                                             multi.annotated.tab.ref$Description == i)]
  
  }

VOG.manual.annot<-multi.annotated.tab.ref$Description[which(multi.annotated.tab.ref$DB == "VOG")]
for (i in VOG.manual.annot){
  consensus.annotation.tab2$Category[which(consensus.annotation.tab2$DB == "VOG"&
                                    consensus.annotation.tab2$Description == i)] <- multi.annotated.tab.ref$Category[which(multi.annotated.tab.ref$DB=="VOG"&
                                                                                                                             multi.annotated.tab.ref$Description == i)]
  
}

PHROGs.manual.annot<-multi.annotated.tab.ref$Description[which(multi.annotated.tab.ref$DB == "PHROGs")]
for (i in PHROGs.manual.annot){
  consensus.annotation.tab2$Category[which(consensus.annotation.tab2$DB == "PHROGs"&
                                    consensus.annotation.tab2$Description == i)] <- multi.annotated.tab.ref$Category[which(multi.annotated.tab.ref$DB=="PHROGs"&
                                                                                                                             multi.annotated.tab.ref$Description == i)]
  
}

uniprotKB.manual.annot<-multi.annotated.tab.ref$Description[which(multi.annotated.tab.ref$DB == "uniprotKB")]
for (i in uniprotKB.manual.annot){
  consensus.annotation.tab2$Category[which(consensus.annotation.tab2$DB == "uniprotKB"&
                                    consensus.annotation.tab2$Description == i)] <- multi.annotated.tab.ref$Category[which(multi.annotated.tab.ref$DB=="uniprotKB"&
                                                                                                                             multi.annotated.tab.ref$Description == i)]
  
}


pfamA.manual.annot<-multi.annotated.tab.ref$Description[which(multi.annotated.tab.ref$DB == "pfamA")]
for (i in pfamA.manual.annot){
  consensus.annotation.tab2$Category[which(consensus.annotation.tab2$DB == "pfamA"&
                                    consensus.annotation.tab2$Description == i)] <- multi.annotated.tab.ref$Category[which(multi.annotated.tab.ref$DB=="pfamA"&
                                                                                                                             multi.annotated.tab.ref$Description == i)]
  
}

Eggnog.manual.annot<-multi.annotated.tab.ref$Description[which(multi.annotated.tab.ref$DB == "Eggnog")]
for (i in Eggnog.manual.annot){
  consensus.annotation.tab2$Category[which(consensus.annotation.tab2$DB == "Eggnog"&
                                    consensus.annotation.tab2$Description == i)] <- multi.annotated.tab.ref$Category[which(multi.annotated.tab.ref$DB=="Eggnog"&
                                                                                                                             multi.annotated.tab.ref$Description == i)]
  
}

writexl::write_xlsx(consensus.annotation.tab2,"./sole.unknowns.annotated.xlsx")

phage.gene.list <- readxl::read_xlsx("./consensus_annotation_final_2.xlsx")
unique.phage.genes<-unique(phage.gene.list$Description)
write.table(unique.phage.genes, "./unique_phage_gene_list.txt",sep="\t",quote=F, row.names = F, col.names = F)

DUF.genes<-unique.phage.genes[grepl("DUF",unique.phage.genes)]
write.table(DUF.genes, "./unique_phage_DUF_gene_list.txt",sep="\t",quote=F, row.names = F, col.names = F)





##### 2. Load prediction results
phatyp.res<-read.table("../DNAvirus/reprePhages/reanalyzed/phatyp_prediction.tsv", sep ="\t", header =T)
names(phatyp.res)[1] <- "Genome2"
phatyp.res <- phatyp.res[c(1,3)]
nrow(phatyp.res) #3,177
active.genes.info$Genome2 <- active.genes.info$Genome
active.genes.info$Genome2<-gsub("\\|","_",active.genes.info$Genome2)
nrow(active.genes.info)

active.phage.list<-unique(active.genes.info$Genome2)
active.phage.list[active.phage.list%in%phatyp.res$Genome2] #3170
active.phage.list[!(active.phage.list%in%phatyp.res$Genome2)] #71

phatyp.res.add <- data.frame(Genome2=active.phage.list[!(active.phage.list%in%phatyp.res$Genome2)],
                             TYPE = "-")

phatyp.res<-rbind(phatyp.res,phatyp.res.add)

active.genes.info.Lcres<-merge(active.genes.info,phatyp.res, by ="Genome2")
unique(active.genes.info.Lcres[active.genes.info.Lcres$TYPE == "temperate" & active.genes.info.Lcres$category =="integration and excision",]$Genome)




#### add lifestyle prediction result
sig.phrog.res.best.hit2$Genome <- paste0(sig.phrog.res.best.hit2$query,".S")
sig.phrog.res.best.hit2$Genome <- gsub("_+\\d+.S","", sig.phrog.res.best.hit2$Genome)
unique(sig.phrog.res.best.hit2$Genome) #3,246 vMAGs
sig.phrog.res.best.hit2$Genome2 <- gsub("\\|","_",sig.phrog.res.best.hit2$Genome)
phage.all.list<-unique(sig.phrog.res.best.hit2$Genome2)
phage.all.list[!(phage.all.list%in%phatyp.res$Genome2)]
phatyp.res2<-subset(phatyp.res,Genome2 %in% phage.all.list)

sig.phrog.res.best.hit2.LC<-merge(sig.phrog.res.best.hit2,phatyp.res2, by ="Genome2")
nrow(sig.phrog.res.best.hit2.LC)

##### 3. Add additional information - Gene abundance and expression data
###### 3-1. gene abundance
gene.abundance<-function(SampleName){
  count.tab <- read.delim(paste0("../DNAvirus/FunctionalAnnotation/GeneAbundance/DNA/",SampleName,"/quant.sf"), sep ="\t", header = T)
  count.tab <- data.frame(count.tab)
  count.tab$Sample <- SampleName
  return(count.tab)
}

gene.abundance.Y22TmD1 <- gene.abundance("9D1PG")
gene.abundance.Y22TmD2 <- gene.abundance("9D2PG")
gene.abundance.Y22TmD3 <- gene.abundance("9D3PH")
gene.abundance.Y22TmD4 <- gene.abundance("9D4PH")
gene.abundance.Y22TmD5 <- gene.abundance("9D5PG")

gene.abundance.Y22TmM1 <- gene.abundance("9D1NH")
gene.abundance.Y22TmM2 <- gene.abundance("9D2NG")
gene.abundance.Y22TmM3 <- gene.abundance("9D3NG")
gene.abundance.Y22TmM4 <- gene.abundance("9D4NG")
gene.abundance.Y22TmM5 <- gene.abundance("9D5NG")

gene.abundance.Y23TmD1 <- gene.abundance("23D1PA")
gene.abundance.Y23TmD2 <- gene.abundance("23D2PB")
gene.abundance.Y23TmD3 <- gene.abundance("23D3PA")
gene.abundance.Y23TmD4 <- gene.abundance("23D4PA")
gene.abundance.Y23TmD5 <- gene.abundance("23D5PB")

gene.abundance.Y23TmM1 <- gene.abundance("23D1NA")
gene.abundance.Y23TmM2 <- gene.abundance("23D2NA")
gene.abundance.Y23TmM3 <- gene.abundance("23D3NA")
gene.abundance.Y23TmM4 <- gene.abundance("23D4NA")
gene.abundance.Y23TmM5 <- gene.abundance("23D5NA")

head(gene.abundance.Y22TmD1)

gene.abundance.tabs<-rbind(gene.abundance.Y22TmD1,gene.abundance.Y22TmD2,gene.abundance.Y22TmD3,gene.abundance.Y22TmD4,gene.abundance.Y22TmD5,
                           gene.abundance.Y22TmM1,gene.abundance.Y22TmM2,gene.abundance.Y22TmM3,gene.abundance.Y22TmM4,gene.abundance.Y22TmM5,
                           gene.abundance.Y23TmD1,gene.abundance.Y23TmD2,gene.abundance.Y23TmD3,gene.abundance.Y23TmD4,gene.abundance.Y23TmD5,
                           gene.abundance.Y23TmM1,gene.abundance.Y23TmM2,gene.abundance.Y23TmM3,gene.abundance.Y23TmM4,gene.abundance.Y23TmM5)

head(gene.abundance.tabs)
### RPKM (remove genes from Tricholoma matsutake contigs)
gene.abundance.tabs.true<-subset(gene.abundance.tabs, Name %in% unique(vMAG.incidence.with.genes$query))
gene.abundance.tabs.true <- gene.abundance.tabs.true %>% group_by(Sample) %>% mutate(TotalReads = sum(NumReads))

gene.abundance.tabs.true <- gene.abundance.tabs.true %>%
  mutate(RPKM = (NumReads * 1e9) / (TotalReads * Length))


gene.abundance.tabs.wide.count <- reshape2::dcast(gene.abundance.tabs.true, Name ~ Sample, value.var = "NumReads")
gene.abundance.tabs.wide.RPKM <- reshape2::dcast(gene.abundance.tabs.true, Name ~ Sample, value.var = "RPKM")
gene.abundance.tabs.wide.TPM <- reshape2::dcast(gene.abundance.tabs.true, Name ~ Sample, value.var = "TPM")



rownames(gene.abundance.tabs.wide.count) <- gene.abundance.tabs.wide.count$Name
gene.abundance.tabs.wide.count <- gene.abundance.tabs.wide.count[-c(1)]
colnames(gene.abundance.tabs.wide.count)

gene.abundance.tabs.wide.count <- gene.abundance.tabs.wide.count[c(12,14,16,18,20,11,13,15,17,19,2,4,6,8,10,1,3,5,7,9)]
colnames(gene.abundance.tabs.wide.count) <- c("Y22TmD1","Y22TmD2","Y22TmD3","Y22TmD4","Y22TmD5",
                                              "Y22TmM1","Y22TmM2","Y22TmM3","Y22TmM4","Y22TmM5",
                                              "Y23TmD1","Y23TmD2","Y23TmD3","Y23TmD4","Y23TmD5",
                                              "Y23TmM1","Y23TmM2","Y23TmM3","Y23TmM4","Y23TmM5")

rownames(gene.abundance.tabs.wide.RPKM) <- gene.abundance.tabs.wide.RPKM$Name
gene.abundance.tabs.wide.RPKM <- gene.abundance.tabs.wide.RPKM[-c(1)]
colnames(gene.abundance.tabs.wide.RPKM)

gene.abundance.tabs.wide.RPKM <- gene.abundance.tabs.wide.RPKM[c(12,14,16,18,20,11,13,15,17,19,2,4,6,8,10,1,3,5,7,9)]
colnames(gene.abundance.tabs.wide.RPKM) <- c("Y22TmD1","Y22TmD2","Y22TmD3","Y22TmD4","Y22TmD5",
                                              "Y22TmM1","Y22TmM2","Y22TmM3","Y22TmM4","Y22TmM5",
                                              "Y23TmD1","Y23TmD2","Y23TmD3","Y23TmD4","Y23TmD5",
                                              "Y23TmM1","Y23TmM2","Y23TmM3","Y23TmM4","Y23TmM5")


rownames(gene.abundance.tabs.wide.TPM) <- gene.abundance.tabs.wide.TPM$Name
gene.abundance.tabs.wide.TPM <- gene.abundance.tabs.wide.TPM[-c(1)]
colnames(gene.abundance.tabs.wide.TPM)

gene.abundance.tabs.wide.TPM <- gene.abundance.tabs.wide.TPM[c(12,14,16,18,20,11,13,15,17,19,2,4,6,8,10,1,3,5,7,9)]
colnames(gene.abundance.tabs.wide.TPM) <- c("Y22TmD1","Y22TmD2","Y22TmD3","Y22TmD4","Y22TmD5",
                                              "Y22TmM1","Y22TmM2","Y22TmM3","Y22TmM4","Y22TmM5",
                                              "Y23TmD1","Y23TmD2","Y23TmD3","Y23TmD4","Y23TmD5",
                                              "Y23TmM1","Y23TmM2","Y23TmM3","Y23TmM4","Y23TmM5")

#### Gene abundance - count
gene.abundance.tabs.wide.count$query <- rownames(gene.abundance.tabs.wide.count)
names(sig.phrog.res.best.hit2.LC)
vMAG.gene.abundance<-merge(sig.phrog.res.best.hit2.LC, gene.abundance.tabs.wide.count, by ="query",all=T)

#### Gene abundance - RPKM
gene.abundance.tabs.wide.RPKM$query <- rownames(gene.abundance.tabs.wide.RPKM)
vMAG.gene.abundance.RPKM<-merge(sig.phrog.res.best.hit2.LC, gene.abundance.tabs.wide.RPKM, by ="query",all=T)

#### Gene abundance - TPM
gene.abundance.tabs.wide.TPM$query <- rownames(gene.abundance.tabs.wide.TPM)
vMAG.gene.abundance.TPM<-merge(sig.phrog.res.best.hit2.LC, gene.abundance.tabs.wide.TPM, by ="query",all=T)


## Add host info
host.info.table <- php.res.tab.true.tax
names(host.info.table)[3] <- "Genome"
host.info.table.best.hit <- host.info.table %>%group_by(Genome) %>% slice(which.max(score))
host.info.table.best.hit$Genome <- gsub("full+\\|+\\|","full\\|",host.info.table.best.hit$Genome)
host.info.table.best.hit$Genome <- gsub("0_partial+\\|+\\|","0_partial\\|",host.info.table.best.hit$Genome)
#writexl::write_xlsx(host.info.table.best.hit,"../DNAvirus/Host prediction/best.hit.info.xlsx")
vMAG.abundance.count.with.host <- merge(vMAG.gene.abundance, host.info.table.best.hit, by ="Genome")
writexl::write_xlsx(vMAG.abundance.count.with.host, "../DNAvirus/FunctionalAnnotation/phrog/viral_gene_abundance_by_Sample_with_host_info_phpRes.xlsx")

vMAG.abundance.RPKM.with.host <- merge(vMAG.gene.abundance.RPKM,host.info.table.best.hit, by ="Genome")
head(vMAG.abundance.RPKM.with.host)

writexl::write_xlsx(vMAG.abundance.RPKM.with.host, "../DNAvirus/FunctionalAnnotation/phrog/viral_gene_abundance_RPKM_by_Sample_with_host_info_phpRes.xlsx")

vMAG.abundance.TPM.with.host <- merge(vMAG.gene.abundance.TPM, host.info.table.best.hit, by ="Genome")
head(vMAG.abundance.TPM.with.host)

writexl::write_xlsx(vMAG.abundance.TPM.with.host, "../DNAvirus/FunctionalAnnotation/phrog/viral_gene_abundance_TPM_by_Sample_with_host_info_phpRes.xlsx")



###### 3-2. gene expression
gene.expression<-function(SampleName){
  count.tab <- read.delim(paste0("../DNAvirus/FunctionalAnnotation/GeneAbundance/RNA/",SampleName,"/quant.sf"), sep ="\t", header = T)
  count.tab <- data.frame(count.tab)
  count.tab$Sample <- SampleName
  return(count.tab)
}

gene.expression.Y23TmD1 <- gene.expression("23R1PB")
gene.expression.Y23TmD4 <- gene.expression("23R4PE")
gene.expression.Y23TmD5 <- gene.expression("23R5PE")

gene.expression.Y23TmM1 <- gene.expression("23R1NA")
gene.expression.Y23TmM4 <- gene.expression("23R4NE")
gene.expression.Y23TmM5 <- gene.expression("23R5NE")

gene.expression.tabs<-rbind(gene.expression.Y23TmD1,gene.expression.Y23TmD4,gene.expression.Y23TmD5,
                            gene.expression.Y23TmM1,gene.expression.Y23TmM4,gene.expression.Y23TmM5)

head(gene.expression.tabs)
### RPKM (remove genes from Tricholoma matsutake contigs)
gene.expression.tabs.true<-subset(gene.expression.tabs, Name %in% unique(vMAG.incidence.with.genes$query))
gene.expression.tabs.true <- gene.expression.tabs.true %>% group_by(Sample) %>% mutate(TotalReads = sum(NumReads))

gene.expression.tabs.true <- gene.expression.tabs.true %>%
  mutate(RPKM = (NumReads * 1e9) / (TotalReads * Length))

gene.expression.tabs.wide.count <- reshape2::dcast(gene.expression.tabs.true, Name ~ Sample, value.var = "NumReads")
gene.expression.tabs.wide.RPKM <- reshape2::dcast(gene.expression.tabs.true, Name ~ Sample, value.var = "RPKM")
gene.expression.tabs.wide.TPM <- reshape2::dcast(gene.expression.tabs.true, Name ~ Sample, value.var = "TPM")


rownames(gene.expression.tabs.wide.count) <- gene.expression.tabs.wide.count$Name
gene.expression.tabs.wide.count <- gene.expression.tabs.wide.count[-c(1)]
colnames(gene.expression.tabs.wide.count)

gene.expression.tabs.wide.count <- gene.expression.tabs.wide.count[c(2,4,6,1,3,5)]
colnames(gene.expression.tabs.wide.count) <- c("Y23TmD1_R","Y23TmD4_R","Y23TmD5_R",
                                               "Y23TmM1_R","Y23TmM4_R","Y23TmM5_R")

rownames(gene.expression.tabs.wide.RPKM) <- gene.expression.tabs.wide.RPKM$Name
gene.expression.tabs.wide.RPKM <- gene.expression.tabs.wide.RPKM[-c(1)]
colnames(gene.expression.tabs.wide.RPKM)

gene.expression.tabs.wide.RPKM <- gene.expression.tabs.wide.RPKM[c(2,4,6,1,3,5)]
colnames(gene.expression.tabs.wide.RPKM) <- c("Y23TmD1_R","Y23TmD4_R","Y23TmD5_R",
                                              "Y23TmM1_R","Y23TmM4_R","Y23TmM5_R")


rownames(gene.expression.tabs.wide.TPM) <- gene.expression.tabs.wide.TPM$Name
gene.expression.tabs.wide.TPM <- gene.expression.tabs.wide.TPM[-c(1)]
colnames(gene.expression.tabs.wide.TPM)

gene.expression.tabs.wide.TPM <- gene.expression.tabs.wide.TPM[c(2,4,6,1,3,5)]
colnames(gene.expression.tabs.wide.TPM) <- c("Y23TmD1_R","Y23TmD4_R","Y23TmD5_R",
                                              "Y23TmM1_R","Y23TmM4_R","Y23TmM5_R")

#### Gene expression - count
gene.expression.tabs.wide.count$query <- rownames(gene.expression.tabs.wide.count)
names(sig.phrog.res.best.hit2.LC)
vMAG.gene.expression<-merge(sig.phrog.res.best.hit2.LC, gene.expression.tabs.wide.count, by ="query",all=T)

#### Gene expression - RPKM
gene.expression.tabs.wide.RPKM$query <- rownames(gene.expression.tabs.wide.RPKM)
vMAG.gene.expression.RPKM<-merge(sig.phrog.res.best.hit2.LC, gene.expression.tabs.wide.RPKM, by ="query",all=T)

#### Gene expression - TPM
gene.expression.tabs.wide.TPM$query <- rownames(gene.expression.tabs.wide.TPM)
vMAG.gene.expression.TPM<-merge(sig.phrog.res.best.hit2.LC, gene.expression.tabs.wide.TPM, by ="query",all=T)


## Add host info
vMAG.expression.count.with.host <- merge(vMAG.gene.expression, host.info.table.best.hit, by ="Genome")
head(vMAG.expression.count.with.host)

writexl::write_xlsx(vMAG.expression.count.with.host, "../DNAvirus/FunctionalAnnotation/phrog/viral_gene_expression_by_Sample_with_host_info_phpRes.xlsx")

vMAG.expression.RPKM.with.host <- merge(vMAG.gene.expression.RPKM,  host.info.table.best.hit, by ="Genome")
head(vMAG.expression.RPKM.with.host)

writexl::write_xlsx(vMAG.expression.RPKM.with.host, "../DNAvirus/FunctionalAnnotation/phrog/viral_gene_expression_RPKM_by_Sample_with_host_info_phpRes.xlsx")

vMAG.expression.TPM.with.host <- merge(vMAG.gene.expression.TPM,  host.info.table.best.hit, by ="Genome")
head(vMAG.expression.TPM.with.host)

writexl::write_xlsx(vMAG.expression.TPM.with.host, "../DNAvirus/FunctionalAnnotation/phrog/viral_gene_expression_TPM_by_Sample_with_host_info_phpRes.xlsx")


########## Assign activity information in each sample 
vMAG.expression.TPM.with.host
#### Agglomerate gene expression data into each phage
names(vMAG.expression.TPM.with.host)
phage.expr.tab <- vMAG.expression.TPM.with.host[c(1,2,4,17,18,19:25,30:36)]
head(phage.expr.tab)
phage.expr.tab.melt <- reshape2::melt(phage.expr.tab)
phage.expr.tab.melt.summary <- phage.expr.tab.melt %>% group_by(Genome,superkingdom,phylum,class,order,family,genus,species,TYPE,variable) %>% reframe(TPM = sum(value))

####### active and inactive phages in each sample
length(unique(phage.expr.tab.melt.summary$Genome)) #3,246 DNA viruses are active in at least one samples
phage.expr.tab.melt.summary$TYPE[which(phage.expr.tab.melt.summary$TYPE=="-")] <- "unknown"
phage.expr.tab.melt.summary$activity <- ifelse(phage.expr.tab.melt.summary$TPM > 0,paste0("active_",phage.expr.tab.melt.summary$TYPE),paste0("inactive_",phage.expr.tab.melt.summary$TYPE))
unique(phage.expr.tab.melt.summary$activity)

active.lytic.Y23TmM5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_virulent" & phage.expr.tab.melt.summary$variable == "Y23TmM5_R")]
active.lytic.Y23TmM4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_virulent" & phage.expr.tab.melt.summary$variable == "Y23TmM4_R")]
active.lytic.Y23TmM1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_virulent" & phage.expr.tab.melt.summary$variable == "Y23TmM1_R")]
active.lytic.Y23TmD1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_virulent" & phage.expr.tab.melt.summary$variable == "Y23TmD1_R")]
active.lytic.Y23TmD4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_virulent" & phage.expr.tab.melt.summary$variable == "Y23TmD4_R")]
active.lytic.Y23TmD5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_virulent" & phage.expr.tab.melt.summary$variable == "Y23TmD5_R")]

active.temper.Y23TmM5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_temperate" & phage.expr.tab.melt.summary$variable == "Y23TmM5_R")]
active.temper.Y23TmM4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_temperate" & phage.expr.tab.melt.summary$variable == "Y23TmM4_R")]
active.temper.Y23TmM1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_temperate" & phage.expr.tab.melt.summary$variable == "Y23TmM1_R")]
active.temper.Y23TmD1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_temperate" & phage.expr.tab.melt.summary$variable == "Y23TmD1_R")]
active.temper.Y23TmD4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_temperate" & phage.expr.tab.melt.summary$variable == "Y23TmD4_R")]
active.temper.Y23TmD5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_temperate" & phage.expr.tab.melt.summary$variable == "Y23TmD5_R")]

active.unknown.Y23TmM5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_unknown" & phage.expr.tab.melt.summary$variable == "Y23TmM5_R")]
active.unknown.Y23TmM4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_unknown" & phage.expr.tab.melt.summary$variable == "Y23TmM4_R")]
active.unknown.Y23TmM1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_unknown" & phage.expr.tab.melt.summary$variable == "Y23TmM1_R")]
active.unknown.Y23TmD1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_unknown" & phage.expr.tab.melt.summary$variable == "Y23TmD1_R")]
active.unknown.Y23TmD4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_unknown" & phage.expr.tab.melt.summary$variable == "Y23TmD4_R")]
active.unknown.Y23TmD5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity=="active_unknown" & phage.expr.tab.melt.summary$variable == "Y23TmD5_R")]

inactive.Y23TmM5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity%in%c("inactive_virulent","inactive_temperate","inactive_unknown") & phage.expr.tab.melt.summary$variable == "Y23TmM5_R")]
inactive.Y23TmM4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity%in%c("inactive_virulent","inactive_temperate","inactive_unknown") & phage.expr.tab.melt.summary$variable == "Y23TmM4_R")]
inactive.Y23TmM1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity%in%c("inactive_virulent","inactive_temperate","inactive_unknown") & phage.expr.tab.melt.summary$variable == "Y23TmM1_R")]
inactive.Y23TmD1 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity%in%c("inactive_virulent","inactive_temperate","inactive_unknown") & phage.expr.tab.melt.summary$variable == "Y23TmD1_R")]
inactive.Y23TmD4 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity%in%c("inactive_virulent","inactive_temperate","inactive_unknown") & phage.expr.tab.melt.summary$variable == "Y23TmD4_R")]
inactive.Y23TmD5 <- phage.expr.tab.melt.summary$Genome[which(phage.expr.tab.melt.summary$activity%in%c("inactive_virulent","inactive_temperate","inactive_unknown") & phage.expr.tab.melt.summary$variable == "Y23TmD5_R")]



######## Bacterial host of active-lytic, active-other, and inactive
host.info.table.best.hit2 <-host.info.table.best.hit[c(3,6:12)]
names(host.info.table.best.hit2)[2:8] <- paste0("host.",names(host.info.table.best.hit2)[2:8])
names(host.info.table.best.hit2)[1] <- "OTU"

##### rename unknowns
host.info.table.best.hit2$host.superkingdom2 <- host.info.table.best.hit2$host.superkingdom
host.info.table.best.hit2$host.phylum2 <- host.info.table.best.hit2$host.phylum
host.info.table.best.hit2$host.class2 <- host.info.table.best.hit2$host.class
host.info.table.best.hit2$host.order2 <- host.info.table.best.hit2$host.order
host.info.table.best.hit2$host.family2 <- host.info.table.best.hit2$host.family
host.info.table.best.hit2$host.genus2 <- host.info.table.best.hit2$host.genus
host.info.table.best.hit2$host.species2 <- host.info.table.best.hit2$host.species

unique(host.info.table.best.hit2$host.superkingdom2)
host.info.table.best.hit2$host.phylum2[which(host.info.table.best.hit2$host.superkingdom != "Unknown"&
                                                  host.info.table.best.hit2$host.phylum == "Unknown")] <- paste0(host.info.table.best.hit2$host.superkingdom[which(host.info.table.best.hit2$host.superkingdom != "Unknown"&
                                                                                                                                                                           host.info.table.best.hit2$host.phylum == "Unknown")],"_unknown")
sort(unique(host.info.table.best.hit2$host.phylum2))

unique(host.info.table.best.hit2$host.class2)
host.info.table.best.hit2$host.class2[which(host.info.table.best.hit2$host.superkingdom != "Unknown" & 
                                                 host.info.table.best.hit2$host.phylum == "Unknown"& 
                                                 host.info.table.best.hit2$host.class == "Unknown")] <- paste0(host.info.table.best.hit2$host.superkingdom2[which(host.info.table.best.hit2$host.superkingdom != "Unknown" & 
                                                                                                                                                                          host.info.table.best.hit2$host.phylum == "Unknown"& 
                                                                                                                                                                          host.info.table.best.hit2$host.class == "Unknown")],"_unknown")
host.info.table.best.hit2$host.class2[which(host.info.table.best.hit2$host.superkingdom != "Unknown" & 
                                                 host.info.table.best.hit2$host.phylum != "Unknown" & 
                                                 host.info.table.best.hit2$host.class == "Unknown")] <- paste0(host.info.table.best.hit2$host.phylum2[which(host.info.table.best.hit2$host.superkingdom != "Unknown" & 
                                                                                                                                                                    host.info.table.best.hit2$host.phylum != "Unknown" & 
                                                                                                                                                                    host.info.table.best.hit2$host.class == "Unknown")],"_unknown")
unique(host.info.table.best.hit2$host.class2)



unique(host.info.table.best.hit2$host.order2)
host.info.table.best.hit2$host.order2[which(host.info.table.best.hit2$host.order =="Unknown"&
                                                 grepl("_unknown",host.info.table.best.hit2$host.class2))] <-host.info.table.best.hit2$host.class2[which(host.info.table.best.hit2$host.order =="Unknown"&
                                                                                                                                                                 grepl("_unknown",host.info.table.best.hit2$host.class2))] 
host.info.table.best.hit2$host.order2[which(host.info.table.best.hit2$host.order =="Unknown"&
                                                 !(grepl("_unknown",host.info.table.best.hit2$host.class2)))] <-paste0(host.info.table.best.hit2$host.class2[which(host.info.table.best.hit2$host.order =="Unknown"&
                                                                                                                                                                           !(grepl("_unknown",host.info.table.best.hit2$host.class2)))],"_unknown")
host.info.table.best.hit2[host.info.table.best.hit2$host.order2 =="Unknown",]
unique(host.info.table.best.hit2$host.order2)

unique(host.info.table.best.hit2$host.family2)
host.info.table.best.hit2$host.family2[which(host.info.table.best.hit2$host.family =="Unknown"&
                                                  grepl("_unknown",host.info.table.best.hit2$host.order2))] <-host.info.table.best.hit2$host.order2[which(host.info.table.best.hit2$host.family =="Unknown"&
                                                                                                                                                                  grepl("_unknown",host.info.table.best.hit2$host.order2))]

host.info.table.best.hit2$host.family2[which(host.info.table.best.hit2$host.family =="Unknown"&
                                                  !(grepl("_unknown",host.info.table.best.hit2$host.order2)))] <-paste0(host.info.table.best.hit2$host.order2[which(host.info.table.best.hit2$host.family =="Unknown"&
                                                                                                                                                                            !(grepl("_unknown",host.info.table.best.hit2$host.order2)))],"_unknown")

sort(unique(host.info.table.best.hit2$host.family2))

host.info.table.best.hit2[host.info.table.best.hit2$host.family2 =="Unknown",]


unique(host.info.table.best.hit2$host.genus2)
host.info.table.best.hit2$host.genus2[which(host.info.table.best.hit2$host.genus =="Unknown"&
                                                 grepl("_unknown",host.info.table.best.hit2$host.family2))] <-host.info.table.best.hit2$host.family2[which(host.info.table.best.hit2$host.genus =="Unknown"&
                                                                                                                                                                   grepl("_unknown",host.info.table.best.hit2$host.family2))]

host.info.table.best.hit2$host.genus2[which(host.info.table.best.hit2$host.genus =="Unknown"&
                                                 !(grepl("_unknown",host.info.table.best.hit2$host.family2)))] <-paste0(host.info.table.best.hit2$host.family2[which(host.info.table.best.hit2$host.genus =="Unknown"&
                                                                                                                                                                             !(grepl("_unknown",host.info.table.best.hit2$host.family2)))],"_unknown")

sort(unique(host.info.table.best.hit2$host.genus2))

host.info.table.best.hit2[host.info.table.best.hit2$host.genus2 =="Unknown",]


unique(host.info.table.best.hit2$host.species2)
host.info.table.best.hit2$host.species2[which(host.info.table.best.hit2$host.species =="Unknown"&
                                              !(grepl("_unknown",host.info.table.best.hit2$host.genus2)))] <-paste0(host.info.table.best.hit2$host.genus2[which(host.info.table.best.hit2$host.species =="Unknown"&
                                                                                                                                                                    !(grepl("_unknown",host.info.table.best.hit2$host.genus2)))],"_unknown")

sort(unique(host.info.table.best.hit2$host.species2))

host.info.table.best.hit2[host.info.table.best.hit2$host.species2 =="Unknown",]


writexl::write_xlsx(host.info.table.best.hit2,"../Manuscript/Supplementary Table/Table S_Host prediction result using PHP_species level_reanalyzed.xlsx")


host.info.table.best.hit2
host.info.table.best.hit2
writexl::write_xlsx(host.info.table.best.hit2, "../../../SNU_IHBae/Manuscript/Comprehensive/Supplementary Files/Data/DataS15_hostprediction.xlsx")


###### Viral abundance (RNA reads mapped to genes)
head(vMAG.gene.expression)
vMAG.gene.expression.reads <- vMAG.gene.expression[c(1,17:23)]
vMAG.gene.expression.melt <- reshape2::melt(vMAG.gene.expression.reads)
head(vMAG.gene.expression.melt )
vMAG.rna.abundance.from.genes <- vMAG.gene.expression.melt %>% group_by(Genome,variable) %>% reframe(reads = sum(value))
max(vMAG.rna.abundance.from.genes$reads)

vMAG.rna.abundance.from.genes.wide <- reshape2::dcast(vMAG.rna.abundance.from.genes,Genome~variable, value.var = "reads")
colSums(vMAG.rna.abundance.from.genes.wide[-c(1)])
sum(colSums(vMAG.rna.abundance.from.genes.wide[-c(1)]))

nrow(vMAG.rna.abundance.from.genes.wide)

### Relative abundance
vMAG.rna.abundance.from.genes <- vMAG.rna.abundance.from.genes%>% group_by(variable) %>% mutate(total = sum(reads), ratio = reads/total)
sum(vMAG.rna.abundance.from.genes$reads)
writexl::write_xlsx(vMAG.rna.abundance.from.genes,"../DNAvirus/Host prediction/viral abundance_MetaT_from_mapped_reads.xlsx")

### add host information
host.info.table.best.hit3 <- host.info.table.best.hit2
names(host.info.table.best.hit3)[1] <- "Genome"
vMAG.rna.abundance.from.genes.host.php<-merge(vMAG.rna.abundance.from.genes,host.info.table.best.hit3, by="Genome")

vMAG.rna.abundance.from.genes.host.php.summary <- vMAG.rna.abundance.from.genes.host.php %>% group_by(host.genus2) %>% reframe(sumReads =sum(reads)) %>% arrange(desc(sumReads))
head(vMAG.rna.abundance.from.genes.host.php.summary)
print(vMAG.rna.abundance.from.genes.host.php.summary, n=700)
vMAG.rna.abundance.from.genes.host.php.summary$ratio <- vMAG.rna.abundance.from.genes.host.php.summary$sumReads/sum(vMAG.rna.abundance.from.genes.host.php.summary$sumReads)
vMAG.rna.abundance.from.genes.host.php.summary
vMAG.rna.abundance.from.genes.host.php.summary$host.genus3 <- vMAG.rna.abundance.from.genes.host.php.summary$host.genus2
vMAG.rna.abundance.from.genes.host.php.summary$host.genus3[which(vMAG.rna.abundance.from.genes.host.php.summary$ratio < 0.005)] <- "Others"

phage.order.rna.re<-unique(vMAG.rna.abundance.from.genes.host.php.summary$host.genus3)

vMAG.host.RelAbund.tab.php 
vMAG.rna.abundance.from.genes.host.php.summary$host.genus2 <- as.character(vMAG.rna.abundance.from.genes.host.php.summary$host.genus2)
vMAG.rna.abundance.from.genes.host.php.summary$host.genus3<-vMAG.rna.abundance.from.genes.host.php.summary$host.genus2
max(vMAG.rna.abundance.from.genes.host.php.summary$ratio)
vMAG.rna.abundance.from.genes.host.php.summary$host.genus3[which(vMAG.rna.abundance.from.genes.host.php.summary$ratio < 0.005)] <- "Others"
vMAG.rna.abundance.from.genes.host.php.summary.2 <- vMAG.rna.abundance.from.genes.host.php.summary[c(1,4)]
head(vMAG.rna.abundance.from.genes.host.php)
vMAG.rna.abundance.from.genes.host.php.2<-merge(vMAG.rna.abundance.from.genes.host.php,vMAG.rna.abundance.from.genes.host.php.summary.2, by = "host.genus2")

head(vMAG.rna.abundance.from.genes.host.php)
head(vMAG.rna.abundance.from.genes.host.php)
vMAG.rna.abundance.from.genes.host.php.2$host.genus3 <- factor(vMAG.rna.abundance.from.genes.host.php.2$host.genus3, levels = rev(phage.order.rna.re))
vMAG.rna.abundance.from.genes.host.php.2$variable <- factor(vMAG.rna.abundance.from.genes.host.php.2$variable, levels = order.rna.sample)


host.genus_colors <- c(
  "Bradyrhizobium" = "#ff7f0e",               # Orange (Consistent)
  "Nocardia" = "#c71585",                     # Medium Violet Red
  "Mycobacterium" = "#2ca02c",                # Green (Consistent)
  "Alphaproteobacteria_unknown" = "#d62728",  # Red (Consistent)
  "Arthrobacter" = "#e377c2",                 # Pink (Consistent)
  "Betaproteobacteria_unknown" = "#98df8a",   # Light Green (Consistent)
  "Streptomyces" = "#bcbd22",                 # Olive (Consistent)
  "Acidobacteriota_unknown" = "#7f7f7f",      # Gray (Consistent)
  "Amycolatopsis" = "#bcbd22",                # Olive (Consistent)
  "Dietzia" = "#4682b4",                      # Steel Blue
  "Chloroflexota_unknown" = "#ffbb78",        # Light Orange (Consistent)
  "Frankia" = "#c5b0d5",                      # Light Purple (Consistent)
  "Actinomadura" = "#c7c7c7",                 # Light Gray (Consistent)
  "Corynebacterium" = "#c49c94",              # Light Brown (Consistent)
  "Actinoplanes" = "#637939",                 # Dark Green (Consistent)
  "Bordetella" = "#ffd700",                   # Gold
  "Gemmatimonadota_unknown" = "#9edae5",      # Soft Cyan (Consistent)
  "Pseudonocardiales_unknown" = "#843c39",    # Dark Red (Consistent)
  "Afipia" = "#ff9896",                       # Light Red (Consistent)
  "Thermocrispum" = "#bd9e39",                # Dark Yellow (Consistent)
  "Rhodospirillales_unknown" = "#7b4173",     # Dark Purple (Consistent)
  "Bacteria_unknown" = "#e7969c",             # Light Magenta (Consistent)
  "Deltaproteobacteria_unknown" = "#aec7e8",  # Light Blue (Consistent)
  "Mycolicibacterium" = "#8c6d31",            # Dark Brown (Consistent)
  "Saccharomonospora" = "#ce6dbd",            # Magenta (Consistent)
  "Sphingomonas" = "#393b79",                 # Dark Blue (Consistent)
  "Thermostaphylospora" = "#ba55d3",          # Medium Orchid
  "Actinomycetota_unknown" = "#66c2a5",       # Aqua
  "Rhodococcus" = "#ff6347",                  # Tomato
  "Actinosynnema" = "#00bfff",                # Deep Sky Blue
  "Ktedonobacter" = "#1f77b4",                # Blue (Consistent)
  "Beijerinckiaceae_unknown" = "#ff33cc",     # Bright Pink (Consistent)
  "Gammaproteobacteria_unknown" = "#8c564b",  # Brown (Consistent)
  "Demequinaceae_unknown" = "#db7093",        # Pale Violet Red (Consistent)
  "Thermomonospora" = "#8b0000",              # Dark Red
  "Phenylobacterium" = "#9acd32",             # Yellow Green (Consistent)
  "Candidatus Bathyarchaeota_unknown" = "#f08080",  # Light Coral
  "Anaeromyxobacter" = "#32cd32",             # Lime Green
  "Dactylosporangium" = "#dda0dd",            # Plum (Consistent)
  "Actinopolymorpha" = "#ffb6c1",             # Light Pink (Consistent)
  "Others" = "#999999"                        # Soft Pink (Consistent)
)


ggplot(vMAG.rna.abundance.from.genes.host.php.2, aes(x=variable, y = ratio, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 3,reverse = T))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


#

ggplot(vMAG.rna.abundance.from.genes.host.php.2, aes(x=variable, y = ratio, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 3,reverse = T))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

### Divide by life styles
##### active-lytic
######## Split by the activity info
### Active-lytic phages
head(vMAG.rna.abundance.from.genes.host.php.2)
Y23TmM5.active.lytic.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.lytic.Y23TmM5 &variable == "Y23TmM5_R")
Y23TmM4.active.lytic.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.lytic.Y23TmM4 & variable == "Y23TmM4_R")
Y23TmM1.active.lytic.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.lytic.Y23TmM1 & variable == "Y23TmM1_R")
Y23TmD1.active.lytic.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.lytic.Y23TmD1 &variable == "Y23TmD1_R")
Y23TmD4.active.lytic.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.lytic.Y23TmD4 & variable == "Y23TmD4_R")
Y23TmD5.active.lytic.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.lytic.Y23TmD5 & variable == "Y23TmD5_R")

active.lytic.host<-rbind(Y23TmM5.active.lytic.host,Y23TmM4.active.lytic.host,Y23TmM1.active.lytic.host,
                         Y23TmD1.active.lytic.host,Y23TmD4.active.lytic.host,Y23TmD5.active.lytic.host)
active.lytic.host$variable <- factor(active.lytic.host$variable, levels = order.rna.sample)

ggplot(active.lytic.host, aes(x=variable, y = ratio, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())



### Active-other and inactive (dormant) phages
Y23TmM5.active.other.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.temper.Y23TmM5 &variable == "Y23TmM5_R")
Y23TmM4.active.other.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.temper.Y23TmM4 & variable == "Y23TmM4_R")
Y23TmM1.active.other.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.temper.Y23TmM1 & variable == "Y23TmM1_R")
Y23TmD1.active.other.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.temper.Y23TmD1 &variable == "Y23TmD1_R")
Y23TmD4.active.other.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.temper.Y23TmD4 & variable == "Y23TmD4_R")
Y23TmD5.active.other.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.temper.Y23TmD5 & variable == "Y23TmD5_R")

active.other.host<-rbind(Y23TmM5.active.other.host,Y23TmM4.active.other.host,Y23TmM1.active.other.host,
                         Y23TmD1.active.other.host,Y23TmD4.active.other.host,Y23TmD5.active.other.host)
active.other.host$variable <- factor(active.other.host$variable, levels = order.rna.sample)

ggplot(active.other.host, aes(x=variable, y = ratio, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


### Active-unknown and inactive (dormant) phages
Y23TmM5.active.unknown.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.unknown.Y23TmM5 &variable == "Y23TmM5_R")
Y23TmM4.active.unknown.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.unknown.Y23TmM4 & variable == "Y23TmM4_R")
Y23TmM1.active.unknown.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.unknown.Y23TmM1 & variable == "Y23TmM1_R")
Y23TmD1.active.unknown.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.unknown.Y23TmD1 &variable == "Y23TmD1_R")
Y23TmD4.active.unknown.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% active.unknown.Y23TmD4 & variable == "Y23TmD4_R")
Y23TmD5.active.unknown.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% active.unknown.Y23TmD5 & variable == "Y23TmD5_R")

active.unknown.host<-rbind(Y23TmM5.active.unknown.host,Y23TmM4.active.unknown.host,Y23TmM1.active.unknown.host,
                         Y23TmD1.active.unknown.host,Y23TmD4.active.unknown.host,Y23TmD5.active.unknown.host)
active.unknown.host$variable <- factor(active.unknown.host$variable, levels = order.rna.sample)

ggplot(active.unknown.host, aes(x=variable, y = ratio, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

### Inactive phages
Y23TmM5.inactive.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% inactive.Y23TmM5 &variable == "Y23TmM5_R")
Y23TmM4.inactive.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% inactive.Y23TmM4 & variable == "Y23TmM4_R")
Y23TmM1.inactive.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% inactive.Y23TmM1 & variable == "Y23TmM1_R")
Y23TmD1.inactive.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% inactive.Y23TmD1 &variable == "Y23TmD1_R")
Y23TmD4.inactive.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome %in% inactive.Y23TmD4 & variable == "Y23TmD4_R")
Y23TmD5.inactive.host<- subset(vMAG.rna.abundance.from.genes.host.php.2, Genome%in% inactive.Y23TmD5 & variable == "Y23TmD5_R")

inactive.host<-rbind(Y23TmM5.inactive.host,Y23TmM4.inactive.host,Y23TmM1.inactive.host,
                     Y23TmD1.inactive.host,Y23TmD4.inactive.host,Y23TmD5.inactive.host)
inactive.host$variable <- factor(inactive.host$variable, levels = order.rna.sample)

ggplot(inactive.host, aes(x=variable, y = ratio, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


##### Supplementary Tables
writexl::write_xlsx(active.lytic.host.2, "../Manuscript/Supplementary Table/Table S_host distribution_metaT_active_lytic_phages.xlsx")
writexl::write_xlsx(active.other.host.2, "../Manuscript/Supplementary Table/Table S_host distribution_metaT_active_lysogenic_phages.xlsx")
writexl::write_xlsx(inactive.host.2, "../Manuscript/Supplementary Table/Table S_host distribution_metaT_inactive_phages.xlsx")
writexl::write_xlsx(vMAG.host.RelAbund.tab.php, "../Manuscript/Supplementary Table/Table S_host distribution_metaG_all_phages.xlsx")






df.vMAG.woNC.phylo.dna.rel.host.2<-merge(df.vMAG.woNC.phylo.dna.rel,host.info.table.best.hit2, by="OTU")
df.vMAG.woNC.phylo.dna.rel.host.2


### Active-lytic phages
Y23TmM5.active.lytic.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.lytic.Y23TmM5 & Sample == "Y23TmM5")
Y23TmM4.active.lytic.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.lytic.Y23TmM4 & Sample == "Y23TmM4")
Y23TmM1.active.lytic.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.lytic.Y23TmM1 & Sample == "Y23TmM1")
Y23TmD1.active.lytic.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.lytic.Y23TmD1 & Sample == "Y23TmD1")
Y23TmD4.active.lytic.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.lytic.Y23TmD4 & Sample == "Y23TmD4")
Y23TmD5.active.lytic.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.lytic.Y23TmD5 & Sample == "Y23TmD5")
head(Y23TmM4.active.lytic.host)
active.lytic.host<-rbind(Y23TmM5.active.lytic.host,Y23TmM4.active.lytic.host,Y23TmM1.active.lytic.host,
                         Y23TmD1.active.lytic.host,Y23TmD4.active.lytic.host,Y23TmD5.active.lytic.host)
active.lytic.host$Sample <- factor(active.lytic.host$Sample, levels = c("Y23TmM5","Y23TmM4","Y23TmM1",
                                                                        "Y23TmD1","Y23TmD4","Y23TmD5"))

active.lytic.host$host.genus3 <- active.lytic.host$host.genus2
active.lytic.host$host.genus3[-which(active.lytic.host$host.genus2 %in% phage.order.rna.re)] <-"Others"
active.lytic.host$host.genus3 <- factor(active.lytic.host$host.genus3, levels = rev(phage.order.rna.re))
ggplot(active.lytic.host, aes(x=Sample, y = Abundance, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

### Active-other phages
Y23TmM5.active.other.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.temper.Y23TmM5 & Sample == "Y23TmM5")
Y23TmM4.active.other.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.temper.Y23TmM4 & Sample == "Y23TmM4")
Y23TmM1.active.other.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.temper.Y23TmM1 & Sample == "Y23TmM1")
Y23TmD1.active.other.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.temper.Y23TmD1 & Sample == "Y23TmD1")
Y23TmD4.active.other.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.temper.Y23TmD4 & Sample == "Y23TmD4")
Y23TmD5.active.other.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.temper.Y23TmD5 & Sample == "Y23TmD5")
head(Y23TmM4.active.other.host)
active.other.host<-rbind(Y23TmM5.active.other.host,Y23TmM4.active.other.host,Y23TmM1.active.other.host,
                         Y23TmD1.active.other.host,Y23TmD4.active.other.host,Y23TmD5.active.other.host)
active.other.host$Sample <- factor(active.other.host$Sample, levels = c("Y23TmM5","Y23TmM4","Y23TmM1",
                                                                        "Y23TmD1","Y23TmD4","Y23TmD5"))

active.other.host$host.genus3 <- active.other.host$host.genus2
active.other.host$host.genus3[-which(active.other.host$host.genus2 %in% phage.order.rna.re)] <-"Others"
active.other.host$host.genus3 <- factor(active.other.host$host.genus3, levels = rev(phage.order.rna.re))

active.other.host$host.genus3 <- factor(active.other.host$host.genus3, levels = rev(phage.order.rna.re))
ggplot(active.other.host, aes(x=Sample, y = Abundance, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


### Active-unknown phages
Y23TmM5.active.unknown.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.unknown.Y23TmM5 & Sample == "Y23TmM5")
Y23TmM4.active.unknown.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.unknown.Y23TmM4 & Sample == "Y23TmM4")
Y23TmM1.active.unknown.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.unknown.Y23TmM1 & Sample == "Y23TmM1")
Y23TmD1.active.unknown.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.unknown.Y23TmD1 & Sample == "Y23TmD1")
Y23TmD4.active.unknown.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.unknown.Y23TmD4 & Sample == "Y23TmD4")
Y23TmD5.active.unknown.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% active.unknown.Y23TmD5 & Sample == "Y23TmD5")
head(Y23TmM4.active.unknown.host)
active.unknown.host<-rbind(Y23TmM5.active.unknown.host,Y23TmM4.active.unknown.host,Y23TmM1.active.unknown.host,
                         Y23TmD1.active.unknown.host,Y23TmD4.active.unknown.host,Y23TmD5.active.unknown.host)
active.unknown.host$Sample <- factor(active.unknown.host$Sample, levels = c("Y23TmM5","Y23TmM4","Y23TmM1",
                                                                        "Y23TmD1","Y23TmD4","Y23TmD5"))


active.unknown.host$host.genus3 <- active.unknown.host$host.genus2
active.unknown.host$host.genus3[-which(active.unknown.host$host.genus2 %in% phage.order.rna.re)] <-"Others"
active.unknown.host$host.genus3 <- factor(active.unknown.host$host.genus3, levels = rev(phage.order.rna.re))
ggplot(active.unknown.host, aes(x=Sample, y = Abundance, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())



### Inactive phages
Y23TmM5.inactive.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% inactive.Y23TmM5 & Sample == "Y23TmM5")
Y23TmM4.inactive.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% inactive.Y23TmM4 & Sample == "Y23TmM4")
Y23TmM1.inactive.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% inactive.Y23TmM1 & Sample == "Y23TmM1")
Y23TmD1.inactive.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% inactive.Y23TmD1 & Sample == "Y23TmD1")
Y23TmD4.inactive.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% inactive.Y23TmD4 & Sample == "Y23TmD4")
Y23TmD5.inactive.host<- subset(df.vMAG.woNC.phylo.dna.rel.host.2, OTU %in% inactive.Y23TmD5 & Sample == "Y23TmD5")
head(Y23TmM4.inactive.host)
inactive.host<-rbind(Y23TmM5.inactive.host,Y23TmM4.inactive.host,Y23TmM1.inactive.host,
                     Y23TmD1.inactive.host,Y23TmD4.inactive.host,Y23TmD5.inactive.host)
inactive.host$Sample <- factor(inactive.host$Sample, levels = c("Y23TmM5","Y23TmM4","Y23TmM1",
                                                                "Y23TmD1","Y23TmD4","Y23TmD5"))

inactive.host$host.genus3 <- inactive.host$host.genus2
inactive.host$host.genus3[-which(inactive.host$host.genus2 %in% phage.order.rna.re)] <-"Others"
inactive.host$host.genus3 <- factor(inactive.host$host.genus3, levels = rev(phage.order.rna.re))
ggplot(inactive.host, aes(x=Sample, y = Abundance, fill = host.genus3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.genus_colors) +
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+ylim(0,1)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())



############ Relative abundance of phages indexed by host taxonomy #######
salmon.wo.NC.count.tabs.2
unique(salmon.wo.NC.count.tabs.2$Name)
salmon.wo.NC.count.tabs.3 <- salmon.wo.NC.count.tabs.2
salmon.wo.NC.count.tabs.3$ratio <- salmon.wo.NC.count.tabs.3$NumReads/salmon.wo.NC.count.tabs.3$totalReads
names(salmon.wo.NC.count.tabs.3)[1] <- "Genome"
Sample.list <- data.frame(Sample = unique(salmon.wo.NC.count.tabs.3$Sample),SampleRename = c("Y22TmD1","Y22TmD2","Y22TmD3","Y22TmD4","Y22TmD5",
                                                                                             "Y22TmM1","Y22TmM2","Y22TmM3","Y22TmM4","Y22TmM5",
                                                                                             "Y23TmD1","Y23TmD2","Y23TmD3","Y23TmD4","Y23TmD5",
                                                                                             "Y23TmM1","Y23TmM2","Y23TmM3","Y23TmM4","Y23TmM5"))

salmon.wo.NC.count.tabs.3<-merge(salmon.wo.NC.count.tabs.3,Sample.list,by="Sample")

phage.metaG.abund.tab.host<-merge(salmon.wo.NC.count.tabs.3,host.info.table.best.hit3,by ="Genome")
unique(phage.metaG.abund.tab.host$Genome)


phage.metaG.abund.tab.host.phylum <- phage.metaG.abund.tab.host %>% group_by(host.phylum2, SampleRename) %>% reframe(ratio = sum(ratio))
phage.metaG.abund.tab.host.phylum$host.phylum3 <- phage.metaG.abund.tab.host.phylum $host.phylum2
phage.metaG.abund.tab.host.phylum.mean <- phage.metaG.abund.tab.host.phylum %>% group_by(host.phylum2) %>% reframe(meanRatio =mean(ratio)) %>% arrange(desc(meanRatio))

dorm.phyla.host<-phage.metaG.abund.tab.host.phylum.mean$host.phylum2[which(phage.metaG.abund.tab.host.phylum.mean $meanRatio > 0.005)]
phage.metaG.abund.tab.host.phylum$host.phylum3 <-ifelse(phage.metaG.abund.tab.host.phylum$host.phylum3%in%dorm.phyla.host,
                                                        phage.metaG.abund.tab.host.phylum$host.phylum3,"Others")

order.phyla.host <- append(dorm.phyla.host,"Others")
order.phyla.host<-order.phyla.host[c(1,2,3,4,6:12,5,13)]
phage.metaG.abund.tab.host.phylum$host.phylum3 <- factor(phage.metaG.abund.tab.host.phylum$host.phylum3, levels=rev(order.phyla.host))
phage.metaG.abund.tab.host.phylum$SampleRename <- factor(phage.metaG.abund.tab.host.phylum$SampleRename, levels=order.all.sample)
ggplot(phage.metaG.abund.tab.host.phylum, aes(x=SampleRename, y = ratio, fill = host.phylum3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + 
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Pseudomonadota"="#99cc66",
                               "Actinomycetota"="#ffcc33",
                               "Acidobacteriota"="#99ccff",
                               "Planctomycetota"="#ff9933",
                               "Chloroflexota"="#00c19f",
                               "Verrucomicrobiota"="#cc99ff",
                               "Myxococcota"="#619cff",
                               "RCP2-54"= "#B6B0FF",
                               "Bacteria_unknown"="#474747",
                               "Cyanobacteriota"="#00A5CF",
                               "Gemmatimonadota"="#FF7F50",
                               "Thermomicrobiota"="#E69F00",
                               "Bacteroidota"="indianred2",
                               "Candidatus_Saccharibacteria"="#AF916D",
                               "Vulcanimicrobiota" ="#3257A8",
                               "Others"="#999999")) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


###### Visualization at the phylum level- host composition of active phages ##### (discard inactive or dormant phages)
vMAG.rna.abundance.from.genes.host.php.2$host.phylum3 <- vMAG.rna.abundance.from.genes.host.php.2$host.phylum2
active.phage.phylum.summary <- vMAG.rna.abundance.from.genes.host.php.2 %>% group_by(host.phylum2,variable) %>% reframe(ratio = sum(ratio))
active.phage.phylum.summary.wide <- reshape2::dcast(active.phage.phylum.summary, host.phylum2~variable, value.var = "ratio")
active.phage.phylum.summary.all <- active.phage.phylum.summary%>% group_by(host.phylum2) %>% reframe(ratio = mean(ratio)) %>% arrange(desc(ratio))
active.dorm.host<-active.phage.phylum.summary.all[active.phage.phylum.summary.all$ratio > 0.005,]$host.phylum2

vMAG.rna.abundance.from.genes.host.php.2$host.phylum3[-which(vMAG.rna.abundance.from.genes.host.php.2$host.phylum3 %in% active.dorm.host)] <- "Others"
ord.host.phyla.active <- append(active.dorm.host,"Others")
ord.host.phyla.active <- ord.host.phyla.active[c(1:4,6:13,5,14)]

vMAG.rna.abundance.from.genes.host.php.2$host.phylum3 <- factor(vMAG.rna.abundance.from.genes.host.php.2$host.phylum3,
                                                                levels=rev(ord.host.phyla.active))
vMAG.rna.abundance.from.genes.host.php.2$variable <- factor(vMAG.rna.abundance.from.genes.host.php.2$variable,
                                                                levels=c("Y23TmM5_R","Y23TmM4_R","Y23TmM1_R",
                                                                         "Y23TmD1_R","Y23TmD4_R","Y23TmD5_R"))

ggplot(vMAG.rna.abundance.from.genes.host.php.2, aes(x=variable, y = ratio, fill = host.phylum3)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + 
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Pseudomonadota"="#99cc66",
                               "Actinomycetota"="#ffcc33",
                               "Acidobacteriota"="#99ccff",
                               "Planctomycetota"="#ff9933",
                               "Chloroflexota"="#00c19f",
                               "Verrucomicrobiota"="#cc99ff",
                               "Myxococcota"="#619cff",
                               "RCP2-54"= "#B6B0FF",
                               "Bacteria_unknown"="#474747",
                               "Cyanobacteriota"="#00A5CF",
                               "Gemmatimonadota"="#FF7F50",
                               "Thermomicrobiota"="#E69F00",
                               "Bacteroidota"="indianred2",
                               "Candidatus_Saccharibacteria"="#AF916D",
                               "Vulcanimicrobiota" ="#3257A8",
                               "Candidatus Bathyarchaeota"="#8B5A2B",
                               "Others"="#999999")) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

######## Other taxonomic levels
##### Metagenomic abundance
phage.metaG.abund.tab.host.phylum <- phage.metaG.abund.tab.host %>% group_by(host.phylum2, SampleRename) %>% reframe(ratio = sum(ratio))
phage.metaG.abund.tab.host.class <- phage.metaG.abund.tab.host %>% group_by(host.phylum2,host.class2, SampleRename) %>% reframe(ratio = sum(ratio))
phage.metaG.abund.tab.host.order <- phage.metaG.abund.tab.host %>% group_by(host.phylum2,host.class2,host.order2, SampleRename) %>% reframe(ratio = sum(ratio))
phage.metaG.abund.tab.host.family <- phage.metaG.abund.tab.host %>% group_by(host.phylum2,host.class2,host.order2,host.family2,SampleRename) %>% reframe(ratio = sum(ratio))
phage.metaG.abund.tab.host.genus <- phage.metaG.abund.tab.host %>% group_by(host.phylum2,host.class2,host.order2,host.family2,host.genus2, SampleRename) %>% reframe(ratio = sum(ratio))
phage.metaG.abund.tab.host.species <- phage.metaG.abund.tab.host %>% group_by(host.phylum2,host.class2,host.order2,host.family2,host.genus2,host.species2,SampleRename) %>% reframe(ratio = sum(ratio))


phage.metaG.abund.tab.host.phylum.wide <-reshape2::dcast(phage.metaG.abund.tab.host.phylum,host.phylum2~SampleRename, value.var = "ratio")
phage.metaG.abund.tab.host.class.wide <-reshape2::dcast(phage.metaG.abund.tab.host.class,host.phylum2+host.class2~SampleRename, value.var = "ratio")
phage.metaG.abund.tab.host.order.wide <-reshape2::dcast(phage.metaG.abund.tab.host.order,host.phylum2+host.class2+host.order2~SampleRename, value.var = "ratio")
phage.metaG.abund.tab.host.family.wide <-reshape2::dcast(phage.metaG.abund.tab.host.family,host.phylum2+host.class2+host.order2+host.family2~SampleRename, value.var = "ratio")
phage.metaG.abund.tab.host.genus.wide <-reshape2::dcast(phage.metaG.abund.tab.host.genus,host.phylum2+host.class2+host.order2+host.family2+host.genus2~SampleRename, value.var = "ratio")
phage.metaG.abund.tab.host.species.wide <-reshape2::dcast(phage.metaG.abund.tab.host.species,host.phylum2+host.class2+host.order2+host.family2+host.genus2+host.species2~SampleRename, value.var = "ratio")


##### Metatranscriptomic abundance
vMAG.rna.abundance.from.genes.host.php.2.phylum <- vMAG.rna.abundance.from.genes.host.php.2 %>% group_by(host.phylum2, variable ) %>% reframe(ratio = sum(ratio))
vMAG.rna.abundance.from.genes.host.php.2.class <- vMAG.rna.abundance.from.genes.host.php.2 %>% group_by(host.phylum2,host.class2, variable ) %>% reframe(ratio = sum(ratio))
vMAG.rna.abundance.from.genes.host.php.2.order <- vMAG.rna.abundance.from.genes.host.php.2 %>% group_by(host.phylum2,host.class2,host.order2, variable ) %>% reframe(ratio = sum(ratio))
vMAG.rna.abundance.from.genes.host.php.2.family <- vMAG.rna.abundance.from.genes.host.php.2 %>% group_by(host.phylum2,host.class2,host.order2,host.family2,variable ) %>% reframe(ratio = sum(ratio))
vMAG.rna.abundance.from.genes.host.php.2.genus <- vMAG.rna.abundance.from.genes.host.php.2 %>% group_by(host.phylum2,host.class2,host.order2,host.family2,host.genus2, variable ) %>% reframe(ratio = sum(ratio))
vMAG.rna.abundance.from.genes.host.php.2.species <- vMAG.rna.abundance.from.genes.host.php.2 %>% group_by(host.phylum2,host.class2,host.order2,host.family2,host.genus2,host.species2,variable ) %>% reframe(ratio = sum(ratio))


vMAG.rna.abundance.from.genes.host.php.2.phylum.wide <-reshape2::dcast(vMAG.rna.abundance.from.genes.host.php.2.phylum,host.phylum2~variable, value.var = "ratio")
vMAG.rna.abundance.from.genes.host.php.2.class.wide <-reshape2::dcast(vMAG.rna.abundance.from.genes.host.php.2.class,host.phylum2+host.class2~variable, value.var = "ratio")
vMAG.rna.abundance.from.genes.host.php.2.order.wide <-reshape2::dcast(vMAG.rna.abundance.from.genes.host.php.2.order,host.phylum2+host.class2+host.order2~variable, value.var = "ratio")
vMAG.rna.abundance.from.genes.host.php.2.family.wide <-reshape2::dcast(vMAG.rna.abundance.from.genes.host.php.2.family,host.phylum2+host.class2+host.order2+host.family2~variable, value.var = "ratio")
vMAG.rna.abundance.from.genes.host.php.2.genus.wide <-reshape2::dcast(vMAG.rna.abundance.from.genes.host.php.2.genus,host.phylum2+host.class2+host.order2+host.family2+host.genus2~variable, value.var = "ratio")
vMAG.rna.abundance.from.genes.host.php.2.species.wide <-reshape2::dcast(vMAG.rna.abundance.from.genes.host.php.2.species,host.phylum2+host.class2+host.order2+host.family2+host.genus2+host.species2~variable, value.var = "ratio")

##### Ratio of active and dormant phages in metatranscriptome samples #####
gene.expression.tabs.wide.RPKM
sig.phrog.res.best.hit2.LC
gene.expression.tabs.wide.count

gene.expression.tabs.wide.count$query <- rownames(gene.expression.tabs.wide.count)

phage.gene.expr.count<-merge(sig.phrog.res.best.hit2.LC,gene.expression.tabs.wide.count,by ="query")
names(phage.gene.expr.count)
phage.gene.expr.count.sub <- phage.gene.expr.count[c(1,16,17,18,19:25)]
phage.gene.expr.count.sub.melt<-reshape2::melt(phage.gene.expr.count.sub)
phage.gene.expr.count.sub.melt.summary <- phage.gene.expr.count.sub.melt %>% group_by(Genome,variable,TYPE) %>% reframe(Count = sum(value))
phage.gene.expr.count.sub.melt.summary.wide <- reshape2::dcast(phage.gene.expr.count.sub.melt.summary,
                                                               Genome+TYPE~variable, value.var = "Count")

nrow(phage.gene.expr.count.sub.melt.summary.wide)


phage.gene.expr.count.sub.melt.summary$TYPE[which(phage.gene.expr.count.sub.melt.summary$TYPE == "-")] <- "unknown"
phage.gene.expr.count.sub.melt.summary$activity <- ifelse(phage.gene.expr.count.sub.melt.summary$Count > 0,paste0(phage.gene.expr.count.sub.melt.summary$TYPE,"_active"),
                                                          paste0(phage.gene.expr.count.sub.melt.summary$TYPE,"_inactive"))

phage.gene.expr.count.sub.melt.summary$Number <- 1
phage.gene.expr.count.sub.melt.summary2 <- phage.gene.expr.count.sub.melt.summary %>% group_by(activity, variable) %>% reframe(Count = sum(Number))
phage.gene.expr.count.sub.melt.summary2
phage.gene.expr.count.sub.melt.summary2.wide <- reshape2::dcast(phage.gene.expr.count.sub.melt.summary2,activity~variable, value.var = "Count")

phage.gene.expr.count.sub.melt.summary2 <- phage.gene.expr.count.sub.melt.summary2 %>% group_by(variable) %>% mutate(Total = sum(Count))
phage.gene.expr.count.sub.melt.summary2$Ratio <- phage.gene.expr.count.sub.melt.summary2$Count/phage.gene.expr.count.sub.melt.summary2$Total


phage.gene.expr.count.sub.melt.summary2$activity <- factor(phage.gene.expr.count.sub.melt.summary2$activity,
                                                                levels=rev(c("virulent_active","temperate_active","unknown_active",
                                                                         "virulent_inactive","temperate_inactive","unknown_inactive")))
phage.gene.expr.count.sub.melt.summary2$variable <- factor(phage.gene.expr.count.sub.melt.summary2$variable,
                                                            levels=c("Y23TmM5_R","Y23TmM4_R","Y23TmM1_R",
                                                                     "Y23TmD1_R","Y23TmD4_R","Y23TmD5_R"))

ggplot(phage.gene.expr.count.sub.melt.summary2, aes(x=variable, y = Ratio, fill = activity)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + 
  #scale_fill_discrete() +
  scale_fill_manual(values = c("virulent_active"="#99cc66",
                               "temperate_active"="#ffcc33",
                               "unknown_active"="#99ccff",
                               "virulent_inactive"="#ff9933",
                               "temperate_inactive"="#00c19f",
                               "unknown_inactive"="#cc99ff"
                               )) +
  
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 1,reverse = F))+
  theme(legend.position="right") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


head(initial.phrog.res)
sig.phrog.res.re<-subset(initial.phrog.res.wo.unknowns, align.portion >= 0.7 & e.value < 1e-05 & seq.identity >= 0.3 & align.score >= 50)
nrow(sig.phrog.res.re)
nrow(sig.phrog.res)
length(unique(sig.phrog.res.re$query)) #61285
length(unique(sig.phrog.res$query)) #61285

initial.phrog.res$align.portion <- initial.phrog.res$align.length/initial.phrog.res$target.length
head(initial.phrog.res)
sig.phrog.res.re[sig.phrog.res.re$query =="Y23TmM1IndAssemContig00848982||full_1",]

