######## DEG analysis using NOIseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NOISeq")

library(NOISeq)
library(phyloseq)

data.bac.pre.active <- readData(data = data.frame(otu_table(expr.gene.bac.phy.active.pre.f)), factors = data.frame(sample_data(expr.gene.bac.phy.active.pre.f)))
noiseq_result.bac.pre.active <- noiseqbio(data.bac.pre.active, k = 0.5, norm = "n", factor = "TmDominanceStatus")
head(noiseq_result.bac.pre.active@results)
noiseq.res.bac.pre.active<-data.frame(noiseq_result.bac.pre.active@results[[1]])
names(noiseq.res.bac.pre.active)
subset(noiseq.res.bac.pre.active, prob > 0.95 & log2FC > 2)
subset(noiseq.res.bac.pre.active, prob > 0.95 & log2FC < -2)

data.bac.post.active <- readData(data = data.frame(otu_table(expr.gene.bac.phy.active.post.f)), factors = data.frame(sample_data(expr.gene.bac.phy.active.post.f)))
noiseq_result.bac.post.active <- noiseqbio(data.bac.post.active, k = 0.5, norm = "n", factor = "TmDominanceStatus")
head(noiseq_result.bac.post.active@results)
noiseq.res.bac.post.active<-data.frame(noiseq_result.bac.post.active@results[[1]])
names(noiseq.res.bac.post.active)
subset(noiseq.res.bac.post.active, prob > 0.95 & log2FC > 2)
subset(noiseq.res.bac.post.active, prob > 0.95 & log2FC < -2)

data.bac.pre.post <- readData(data = data.frame(otu_table(expr.gene.bac.phy.pre.post.f)), factors = data.frame(sample_data(expr.gene.bac.phy.pre.post.f)))
noiseq_result.bac.pre.post <- noiseqbio(data.bac.pre.post, k = 0.5, norm = "n", factor = "TmDominanceStatus")
head(noiseq_result.bac.pre.post@results)
noiseq.res.bac.pre.post<-data.frame(noiseq_result.bac.pre.post@results[[1]])
names(noiseq.res.bac.pre.post)
subset(noiseq.res.bac.pre.post, prob > 0.95 & log2FC > 2)
subset(noiseq.res.bac.pre.post, prob > 0.95 & log2FC < -2)

data.fun.pre.active <- readData(data = data.frame(otu_table(expr.gene.fun.phy.active.pre.f)), factors = data.frame(sample_data(expr.gene.fun.phy.active.pre.f)))
noiseq_result.fun.pre.active <- noiseqbio(data.fun.pre.active, k = 0.5, norm = "n", factor = "TmDominanceStatus")
head(noiseq_result.fun.pre.active@results)
noiseq.res.fun.pre.active<-data.frame(noiseq_result.fun.pre.active@results[[1]])
names(noiseq.res.fun.pre.active)
subset(noiseq.res.fun.pre.active, prob > 0.95 & log2FC > 2)
subset(noiseq.res.fun.pre.active, prob > 0.95 & log2FC < -2)

data.fun.post.active <- readData(data = data.frame(otu_table(expr.gene.fun.phy.active.post.f)), factors = data.frame(sample_data(expr.gene.fun.phy.active.post.f)))
noiseq_result.fun.post.active <- noiseqbio(data.fun.post.active, k = 0.5, norm = "n", factor = "TmDominanceStatus")
head(noiseq_result.fun.post.active@results)
noiseq.res.fun.post.active<-data.frame(noiseq_result.fun.post.active@results[[1]])
names(noiseq.res.fun.post.active)
subset(noiseq.res.fun.post.active, prob > 0.95 & log2FC > 2)
subset(noiseq.res.fun.post.active, prob > 0.95 & log2FC < -2)

data.fun.pre.post <- readData(data = data.frame(otu_table(expr.gene.fun.phy.pre.post.f)), factors = data.frame(sample_data(expr.gene.fun.phy.pre.post.f)))
noiseq_result.fun.pre.post <- noiseqbio(data.fun.pre.post, k = 0.5, norm = "n", factor = "TmDominanceStatus")
head(noiseq_result.fun.pre.post@results)
noiseq.res.fun.pre.post<-data.frame(noiseq_result.fun.pre.post@results[[1]])
names(noiseq.res.fun.pre.post)
subset(noiseq.res.fun.pre.post, prob > 0.95 & log2FC > 2)
subset(noiseq.res.fun.pre.post, prob > 0.95 & log2FC < -2)

common.enriched.active.bac<-intersect(rownames(subset(noiseq.res.bac.pre.active, prob > 0.95 & log2FC >= 2)),
          rownames(subset(noiseq.res.bac.post.active, prob > 0.95 & log2FC >= 2)))
save(common.enriched.active.bac, file = "./GeneExpression/Expression/common.enriched.active.bac.RData")

common.enriched.pre.bac<-intersect(rownames(subset(noiseq.res.bac.pre.active, prob > 0.95 & log2FC <= -2)),
          rownames(subset(noiseq.res.bac.pre.post, prob > 0.95 & log2FC <= -2)))
save(common.enriched.pre.bac, file = "./GeneExpression/Expression/common.enriched.pre.bac.RData")

common.enriched.post.bac<-intersect(rownames(subset(noiseq.res.bac.post.active, prob > 0.95 & log2FC <= -2)),
                                   rownames(subset(noiseq.res.bac.pre.post, prob > 0.95 & log2FC >= 2)))
save(common.enriched.post.bac, file = "./GeneExpression/Expression/common.enriched.post.bac.RData")


common.enriched.active.fun<-intersect(rownames(subset(noiseq.res.fun.pre.active, prob > 0.95 & log2FC >= 2)),
                                      rownames(subset(noiseq.res.fun.post.active, prob > 0.95 & log2FC >= 2)))
save(common.enriched.active.fun, file = "./GeneExpression/Expression/common.enriched.active.fun.RData")

common.enriched.pre.fun<-intersect(rownames(subset(noiseq.res.fun.pre.active, prob > 0.95 & log2FC <= -2)),
                                   rownames(subset(noiseq.res.fun.pre.post, prob > 0.95 & log2FC <= -2)))
save(common.enriched.pre.fun, file = "./GeneExpression/Expression/common.enriched.pre.fun.RData")

common.enriched.post.fun<-intersect(rownames(subset(noiseq.res.fun.post.active, prob > 0.95 & log2FC <= -2)),
                                    rownames(subset(noiseq.res.fun.pre.post, prob > 0.95 & log2FC >= 2)))
save(common.enriched.post.fun, file = "./GeneExpression/Expression/common.enriched.post.fun.RData")


###### Assign enrichment groups
head(noiseq.res.bac.pre.active)
remove.list.bac.pre.active<-rownames(noiseq.res.bac.pre.active)[is.na(noiseq.res.bac.pre.active$Active_mean)]
noiseq.res.bac.pre.active.sub<-subset(noiseq.res.bac.pre.active, !(rownames(noiseq.res.bac.pre.active) %in% remove.list.bac.pre.active))
nrow(noiseq.res.bac.pre.active.sub)
noiseq.res.bac.pre.active.sub$Enrichment <- ifelse(noiseq.res.bac.pre.active.sub$log2FC >= 2 & noiseq.res.bac.pre.active.sub$prob > 0.95,"Active",
                                               ifelse(noiseq.res.bac.pre.active.sub$log2FC <= -2 & noiseq.res.bac.pre.active.sub$prob > 0.95,"Pre","Non-differential"))

remove.list.bac.post.active<-rownames(noiseq.res.bac.post.active)[is.na(noiseq.res.bac.post.active$Active_mean)]
noiseq.res.bac.post.active.sub<-subset(noiseq.res.bac.post.active, !(rownames(noiseq.res.bac.post.active) %in% remove.list.bac.post.active))
nrow(noiseq.res.bac.post.active.sub)
noiseq.res.bac.post.active.sub$Enrichment <- ifelse(noiseq.res.bac.post.active.sub$log2FC >= 2 & noiseq.res.bac.post.active.sub$prob > 0.95,"Active",
                                                   ifelse(noiseq.res.bac.post.active.sub$log2FC <= -2 & noiseq.res.bac.post.active.sub$prob > 0.95,"Post","Non-differential"))

remove.list.bac.pre.post<-rownames(noiseq.res.bac.pre.post)[is.na(noiseq.res.bac.pre.post$Pre_mean)]
noiseq.res.bac.pre.post.sub<-subset(noiseq.res.bac.pre.post, !(rownames(noiseq.res.bac.pre.post) %in% remove.list.bac.pre.post))
nrow(noiseq.res.bac.pre.post.sub)
noiseq.res.bac.pre.post.sub$Enrichment <- ifelse(noiseq.res.bac.pre.post.sub$log2FC >= 2 & noiseq.res.bac.pre.post.sub$prob > 0.95,"Post",
                                                   ifelse(noiseq.res.bac.pre.post.sub$log2FC <= -2 & noiseq.res.bac.pre.post.sub$prob > 0.95,"Pre","Non-differential"))


head(noiseq.res.fun.pre.active)
remove.list.fun.pre.active<-rownames(noiseq.res.fun.pre.active)[is.na(noiseq.res.fun.pre.active$Active_mean)]
noiseq.res.fun.pre.active.sub<-subset(noiseq.res.fun.pre.active, !(rownames(noiseq.res.fun.pre.active) %in% remove.list.fun.pre.active))
nrow(noiseq.res.fun.pre.active.sub)
noiseq.res.fun.pre.active.sub$Enrichment <- ifelse(noiseq.res.fun.pre.active.sub$log2FC >= 2 & noiseq.res.fun.pre.active.sub$prob > 0.95,"Active",
                                                   ifelse(noiseq.res.fun.pre.active.sub$log2FC <= -2 & noiseq.res.fun.pre.active.sub$prob > 0.95,"Pre","Non-differential"))

remove.list.fun.post.active<-rownames(noiseq.res.fun.post.active)[is.na(noiseq.res.fun.post.active$Active_mean)]
noiseq.res.fun.post.active.sub<-subset(noiseq.res.fun.post.active, !(rownames(noiseq.res.fun.post.active) %in% remove.list.fun.post.active))
nrow(noiseq.res.fun.post.active.sub)
noiseq.res.fun.post.active.sub$Enrichment <- ifelse(noiseq.res.fun.post.active.sub$log2FC >= 2 & noiseq.res.fun.post.active.sub$prob > 0.95,"Active",
                                                    ifelse(noiseq.res.fun.post.active.sub$log2FC <= -2 & noiseq.res.fun.post.active.sub$prob > 0.95,"Post","Non-differential"))

remove.list.fun.pre.post<-rownames(noiseq.res.fun.pre.post)[is.na(noiseq.res.fun.pre.post$Pre_mean)]
noiseq.res.fun.pre.post.sub<-subset(noiseq.res.fun.pre.post, !(rownames(noiseq.res.fun.pre.post) %in% remove.list.fun.pre.post))
nrow(noiseq.res.fun.pre.post.sub)
noiseq.res.fun.pre.post.sub$Enrichment <- ifelse(noiseq.res.fun.pre.post.sub$log2FC >= 2 & noiseq.res.fun.pre.post.sub$prob > 0.95,"Post",
                                                 ifelse(noiseq.res.fun.pre.post.sub$log2FC <= -2 & noiseq.res.fun.pre.post.sub$prob > 0.95,"Pre","Non-differential"))

write.csv(noiseq.res.bac.pre.active.sub,"./GeneExpression/Expression/noiseq.res.bac.pre.active.sub.csv")
write.csv(noiseq.res.bac.post.active.sub,"./GeneExpression/Expression/noiseq.res.bac.post.active.sub.csv")
write.csv(noiseq.res.bac.pre.post.sub,"./GeneExpression/Expression/noiseq.res.bac.pre.post.sub.csv")

write.csv(noiseq.res.fun.pre.active.sub,"./GeneExpression/Expression/noiseq.res.fun.pre.active.sub.csv")
write.csv(noiseq.res.fun.post.active.sub,"./GeneExpression/Expression/noiseq.res.fun.post.active.sub.csv")
write.csv(noiseq.res.fun.pre.post.sub,"./GeneExpression/Expression/noiseq.res.fun.pre.post.sub.csv")


###### GO enrichment test using topGO
library(topGO)

#### Bacteria-active
####### TopGO package (better than clusterProfiler)
library(topGO)
head(GO.gene.coexpr)

GO.gene.expr <- subset(eggnog.rep.bac.tax.express.RPKM,GeneID %in% rownames(gene.express.bac.RPKM.input))[c(1,11)]
GO.gene.expr$GOs[is.na(GO.gene.expr$GOs)] <-"-"
head(GO.gene.expr)

GO.gene.expr$Active <- ifelse(GO.gene.expr$GeneID %in% common.enriched.active.bac,1,0)
GO.gene.expr$Pre <- ifelse(GO.gene.expr$GeneID %in% common.enriched.pre.bac,1,0)
GO.gene.expr$Post <- ifelse(GO.gene.expr$GeneID %in% common.enriched.post.bac,1,0)
head(GO.gene.expr)

bac.active.geneGO.list<-GO.gene.expr$Active
names(bac.active.geneGO.list) <- GO.gene.expr$GeneID

expr.gene2go<-GO.gene.expr[c(1,2)]
rownames(expr.gene2go) <- expr.gene2go$GeneID
expr.gene2go <- expr.gene2go[-c(1)]
head(expr.gene2go)
bac.expr.gene2go<-strsplit(expr.gene2go$GOs,split = ",")
names(bac.expr.gene2go) <- rownames(expr.gene2go)

bac.active.GOdata <- new(
  "topGOdata",
  ontology = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  allGenes = bac.active.geneGO.list,
  geneSelectionFun = function(x) x == 1,  # 선택 기준: 1이면 유의미한 유전자
  annot = annFUN.gene2GO,  # Annotation 함수
  gene2GO = bac.expr.gene2go       # Annotation 데이터
)

bac.active.resultFisher <- runTest(bac.active.GOdata, algorithm = "classic", statistic = "fisher")

bac.active.allResults <- GenTable(
  bac.active.GOdata,
  classicFisher = bac.active.resultFisher,
  topNodes = 5293  # 상위 10개 GO term만 표시
)
tail(allResults)
bac.active.allResults$classicFisher <- as.numeric(bac.active.allResults$classicFisher)
sig.bac.active.allResults <- subset(bac.active.allResults, classicFisher < 0.05)
active.bac.GO.nitrogen<-sig.bac.active.allResults[grepl("nitrogen", sig.bac.active.allResults$Term),]
active.bac.GO.iron<-sig.bac.active.allResults[grepl("iron", sig.bac.active.allResults$Term),]
active.bac.GO.redox<-sig.bac.active.allResults[grepl("redox", sig.bac.active.allResults$Term),]
active.bac.GO.proteo<-sig.bac.active.allResults[grepl("proteo", sig.bac.active.allResults$Term),]

active.bac.GOs<-rbind(active.bac.GO.nitrogen,active.bac.GO.iron,active.bac.GO.redox,active.bac.GO.proteo)
nrow(active.bac.GOs)
active.bac.GOs <- active.bac.GOs %>% arrange(classicFisher)
writexl::write_xlsx(sig.bac.active.allResults,"../../../Manuscript/Comprehensive/Supplementary Files/Data/DataS_GOenrich_bacDEG_active_significant.xlsx")

#### Plotting


#### Bacteria-Pre
bac.pre.geneGO.list<-GO.gene.expr$Pre
names(bac.pre.geneGO.list) <- GO.gene.expr$GeneID

bac.pre.GOdata <- new(
  "topGOdata",
  ontology = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  allGenes = bac.pre.geneGO.list,
  geneSelectionFun = function(x) x == 1,  # 선택 기준: 1이면 유의미한 유전자
  annot = annFUN.gene2GO,  # Annotation 함수
  gene2GO = bac.expr.gene2go       # Annotation 데이터
)

bac.pre.resultFisher <- runTest(bac.pre.GOdata, algorithm = "classic", statistic = "fisher")
whichAlgorithms()
bac.pre.allResults <- GenTable(
  bac.pre.GOdata,
  classicFisher = bac.pre.resultFisher,
  topNodes = 4000  # 상위 10개 GO term만 표시
)

bac.pre.allResults$classicFisher <- as.numeric(bac.pre.allResults$classicFisher)
sig.bac.pre.allResults <- subset(bac.pre.allResults, classicFisher < 0.05)
writexl::write_xlsx(sig.bac.pre.allResults,"../../../Manuscript/Comprehensive/Supplementary Files/Data/DataS_GOenrich_bacDEG_pre_significant.xlsx")



#### Bacteria -Post
bac.post.geneGO.list<-GO.gene.expr$Post
names(bac.post.geneGO.list) <- GO.gene.expr$GeneID

bac.post.GOdata <- new(
  "topGOdata",
  ontology = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  allGenes = bac.post.geneGO.list,
  geneSelectionFun = function(x) x == 1,  # 선택 기준: 1이면 유의미한 유전자
  annot = annFUN.gene2GO,  # Annotation 함수
  gene2GO = bac.expr.gene2go       # Annotation 데이터
)

bac.post.resultFisher <- runTest(bac.post.GOdata, algorithm = "classic", statistic = "fisher")

bac.post.allResults <- GenTable(
  bac.post.GOdata,
  classicFisher = bac.post.resultFisher,
  topNodes = 4000  # 상위 10개 GO term만 표시
)
bac.post.allResults
bac.post.allResults[grepl("nitrogen",bac.post.allResults$Term),]

bac.post.allResults$classicFisher <- as.numeric(bac.post.allResults$classicFisher)
sig.bac.post.allResults <- subset(bac.post.allResults, classicFisher < 0.05)
writexl::write_xlsx(sig.bac.post.allResults,"../../../Manuscript/Comprehensive/Supplementary Files/Data/DataS_GOenrich_bacDEG_post_significant.xlsx")


##### Fungi
fun.GO.gene.expr <- subset(eggnog.rep.fun.tax.express.RPKM,GeneID %in% rownames(gene.express.fun.RPKM.input))[c(1,11)]
fun.GO.gene.expr$GOs[is.na(fun.GO.gene.expr$GOs)] <-"-"
head(fun.GO.gene.expr)

fun.GO.gene.expr$Active <- ifelse(fun.GO.gene.expr$GeneID %in% common.enriched.active.fun,1,0)
fun.GO.gene.expr$Pre <- ifelse(fun.GO.gene.expr$GeneID %in% common.enriched.pre.fun,1,0)
fun.GO.gene.expr$Post <- ifelse(fun.GO.gene.expr$GeneID %in% common.enriched.post.fun,1,0)
head(fun.GO.gene.expr)

fun.active.geneGO.list<-fun.GO.gene.expr$Active
names(fun.active.geneGO.list) <- fun.GO.gene.expr$GeneID

fun.expr.gene2go<-fun.GO.gene.expr[c(1,2)]
rownames(fun.expr.gene2go) <- fun.expr.gene2go$GeneID
fun.expr.gene2go <- fun.expr.gene2go[-c(1)]
head(fun.expr.gene2go)
fun.expr.gene2go.list<-strsplit(fun.expr.gene2go$GOs,split = ",")
names(fun.expr.gene2go.list) <- rownames(fun.expr.gene2go)

fun.active.GOdata <- new(
  "topGOdata",
  ontology = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  allGenes = fun.active.geneGO.list,
  geneSelectionFun = function(x) x == 1,  # 선택 기준: 1이면 유의미한 유전자
  annot = annFUN.gene2GO,  # Annotation 함수
  gene2GO = fun.expr.gene2go.list       # Annotation 데이터
)

fun.active.resultFisher <- runTest(fun.active.GOdata, algorithm = "classic", statistic = "fisher")

fun.active.allResults <- GenTable(
  fun.active.GOdata,
  classicFisher = fun.active.resultFisher,
  topNodes = 4000  # 상위 10개 GO term만 표시
)


fun.active.allResults$classicFisher <- as.numeric(fun.active.allResults$classicFisher)
sig.fun.active.allResults <- subset(fun.active.allResults, classicFisher < 0.05)
active.fun.GO.nitrogen<-sig.fun.active.allResults[grepl("nitrogen", sig.fun.active.allResults$Term),]
active.fun.GO.iron<-sig.fun.active.allResults[grepl("iron", sig.fun.active.allResults$Term),]
active.fun.GO.pheromone<-sig.fun.active.allResults[grepl("pheromone", sig.fun.active.allResults$Term),]
active.fun.GO.proteo<-sig.fun.active.allResults[grepl("proteo", sig.fun.active.allResults$Term),]

active.fun.GOs<-rbind(active.fun.GO.nitrogen,active.fun.GO.iron,active.fun.GO.pheromone,active.fun.GO.proteo)
nrow(active.fun.GOs)
active.fun.GOs <- active.fun.GOs %>% arrange(classicFisher)

writexl::write_xlsx(sig.fun.active.allResults,"../../../Manuscript/Comprehensive/Supplementary Files/Data/DataS_GOenrich_funDEG_active_significant.xlsx")

#### Bacteria-Pre
fun.pre.geneGO.list<-fun.GO.gene.expr$Pre
names(fun.pre.geneGO.list) <- fun.GO.gene.expr$GeneID

fun.pre.GOdata <- new(
  "topGOdata",
  ontology = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  allGenes = fun.pre.geneGO.list,
  geneSelectionFun = function(x) x == 1,  # 선택 기준: 1이면 유의미한 유전자
  annot = annFUN.gene2GO,  # Annotation 함수
  gene2GO = fun.expr.gene2go.list       # Annotation 데이터
)

fun.pre.resultFisher <- runTest(fun.pre.GOdata, algorithm = "classic", statistic = "fisher")

fun.pre.allResults <- GenTable(
  fun.pre.GOdata,
  classicFisher = fun.pre.resultFisher,
  topNodes = 4000  # 상위 10개 GO term만 표시
)
tail(allResults)

fun.pre.allResults$classicFisher <- as.numeric(fun.pre.allResults$classicFisher)
sig.fun.pre.allResults <- subset(fun.pre.allResults, classicFisher < 0.05)
writexl::write_xlsx(sig.fun.pre.allResults,"../../../Manuscript/Comprehensive/Supplementary Files/Data/DataS_GOenrich_funDEG_pre_significant.xlsx")



#### Bacteria -Post
fun.post.geneGO.list<-fun.GO.gene.expr$Post
names(fun.post.geneGO.list) <- fun.GO.gene.expr$GeneID

fun.post.GOdata <- new(
  "topGOdata",
  ontology = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  allGenes = fun.post.geneGO.list,
  geneSelectionFun = function(x) x == 1,  # 선택 기준: 1이면 유의미한 유전자
  annot = annFUN.gene2GO,  # Annotation 함수
  gene2GO = fun.expr.gene2go.list       # Annotation 데이터
)

fun.post.resultFisher <- runTest(fun.post.GOdata, algorithm = "classic", statistic = "fisher")

fun.post.allResults <- GenTable(
  fun.post.GOdata,
  classicFisher = fun.post.resultFisher,
  topNodes = 4000  # 상위 10개 GO term만 표시
)
tail(allResults)
fun.post.allResults$classicFisher <- as.numeric(fun.post.allResults$classicFisher)
sig.fun.post.allResults <- subset(fun.post.allResults, classicFisher < 0.05)
writexl::write_xlsx(sig.fun.post.allResults,"../../../Manuscript/Comprehensive/Supplementary Files/Data/DataS_GOenrich_funDEG_post_significant.xlsx")





##### Plotting GO overrepresentation test (active site)
active.bac.GOs
active.fun.GOs <- active.fun.GOs[!(active.fun.GOs$GO.ID%in%c("GO:0097428","GO:1901564","GO:1903050","GO:1903052","GO:0042810","GO:0042811")),]
nrow(active.fun.GOs)

active.bac.GOs$Domain <-"Bacteria"
active.fun.GOs$Domain <-"Fungi"

active.bac.GOs <- active.bac.GOs%>% arrange(desc(classicFisher))
active.bac.GOs$GO.ID <- factor(active.bac.GOs$GO.ID, levels =active.bac.GOs$GO.ID )
active.fun.GOs <- active.fun.GOs%>% arrange(desc(classicFisher))
active.fun.GOs$GO.ID <- factor(active.fun.GOs$GO.ID, levels =active.fun.GOs$GO.ID )
levels(unique(c(active.bac.GOs$GO.ID,active.fun.GOs$GO.ID)))
active.GOs.res <- rbind(active.bac.GOs,active.fun.GOs)
active.GOs.res$GO.ID  <- factor(active.GOs.res$GO.ID, levels = levels(unique(c(active.bac.GOs$GO.ID,active.fun.GOs$GO.ID))))

ggplot(active.GOs.res, aes(x=-log(classicFisher), y=GO.ID, color = Domain))+
  geom_point(size = 3,shape=16, alpha = 0.8)+ facet_wrap(~Domain,scales = 'free', ncol=1)+
  #coord_flip()+ 
  # theme(strip.placement = "outside",
  #       strip.background = element_rect(fill=NA,colour=NA),
  #       panel.spacing=unit(0,"cm"), axis.title.y = element_blank()) +
  # annotation_custom(grob = linesGrob(), xmin = -0.75, xmax = -0.75, ymin = -3.25, ymax = -0.75) +
  # coord_cartesian(clip="off")+
  guides(color = guide_legend(reverse =F, ncol = 1))+
  theme(legend.position="right") + #theme(aspect.ratio = 3)+
  #scale_color_manual(values = kegg.path.color)+
  #theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) +
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme_classic()+theme(aspect.ratio = 2)




##### Quinolinate-related genes
head(eggnog.rep.bac.tax.abund.RPKM)
nrow(eggnog.rep.bac.tax.express.RPKM.refine.sub)
eggnog.rep.bac.tax.express.RPKM.refine.sub[grepl("GO:0004514",eggnog.rep.bac.tax.express.RPKM.refine.sub$GOs),]
eggnog.rep.fun.tax.express.RPKM.refine.sub[grepl("GO:0004514",eggnog.rep.fun.tax.express.RPKM.refine.sub$GOs),]
