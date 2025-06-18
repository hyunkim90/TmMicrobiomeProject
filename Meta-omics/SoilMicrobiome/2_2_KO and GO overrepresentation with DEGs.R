###### 1. KO and GO enrichment overrepresentation test of DEGs (on the terminal)######
load(file="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/common.enriched.active.bac.RData")
load(file="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/common.enriched.pre.bac.RData")
load(file="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/common.enriched.post.bac.RData")

KOgene_list<-common.enriched.active.bac
kegg_result.bac<-enricher(KOgene_list,
                          TERM2GENE=KOannotation.bac[c(2,1)],
                          TERM2NAME=KOannotation.bac[c(2,4)])
df.kegg_result.bac <- as.data.frame(kegg_result.bac@result)
write.csv(df.kegg_result.bac,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_deg_active_kegg_result.csv")

KOgene_list<-common.enriched.pre.bac
kegg_result.bac<-enricher(KOgene_list,
                          TERM2GENE=KOannotation.bac[c(2,1)],
                          TERM2NAME=KOannotation.bac[c(2,4)])
df.kegg_result.bac <- as.data.frame(kegg_result.bac@result)
write.csv(df.kegg_result.bac,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_deg_pre_kegg_result.csv")

KOgene_list<-common.enriched.post.bac
kegg_result.bac<-enricher(KOgene_list,
                          TERM2GENE=KOannotation.bac[c(2,1)],
                          TERM2NAME=KOannotation.bac[c(2,4)])
df.kegg_result.bac <- as.data.frame(kegg_result.bac@result)
write.csv(df.kegg_result.bac,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_deg_post_kegg_result.csv")


KOgene_list<-common.enriched.active.bac
keggpath_result.bac.active<-enricher(KOgene_list,
                                 TERM2GENE=KOannotation.bac[c(3,1)],
                                 TERM2NAME=KOannotation.bac[c(3,4)])
df.keggpath_result.bac.active <- as.data.frame(keggpath_result.bac.active@result)
write.csv(df.keggpath_result.bac.active,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_active_kegg_path_result.csv")

KOgene_list<-common.enriched.pre.bac
keggpath_result.bac.pre<-enricher(KOgene_list,
                                     TERM2GENE=KOannotation.bac[c(3,1)],
                                     TERM2NAME=KOannotation.bac[c(3,4)])
df.keggpath_result.bac.pre <- as.data.frame(keggpath_result.bac.pre@result)
write.csv(df.keggpath_result.bac.pre,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_pre_kegg_path_result.csv")

KOgene_list<-common.enriched.post.bac
keggpath_result.bac.post<-enricher(KOgene_list,
                                     TERM2GENE=KOannotation.bac[c(3,1)],
                                     TERM2NAME=KOannotation.bac[c(3,4)])
df.keggpath_result.bac.post <- as.data.frame(keggpath_result.bac.post@result)
write.csv(df.keggpath_result.bac.post,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_post_kegg_path_result.csv")


GOgene_list<-common.enriched.active.bac
gene_go.bac.active<-enricher(GOgene_list,
                           TERM2GENE=GOannotation.bac[c(2,1)],
                           TERM2NAME=GOinfo[1:2])
df.gene_go.bac.active <- as.data.frame(gene_go.bac.active@result)
write.csv(df.gene_go.bac.active,file = "/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_deg_active_GO_result.csv", row.names=F)

GOgene_list<-common.enriched.pre.bac
gene_go.bac.pre<-enricher(GOgene_list,
                             TERM2GENE=GOannotation.bac[c(2,1)],
                             TERM2NAME=GOinfo[1:2])
df.gene_go.bac.pre <- as.data.frame(gene_go.bac.pre@result)
write.csv(df.gene_go.bac.pre,file = "/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_deg_pre_GO_result.csv", row.names=F)


GOgene_list<-common.enriched.post.bac
gene_go.bac.post<-enricher(GOgene_list,
                             TERM2GENE=GOannotation.bac[c(2,1)],
                             TERM2NAME=GOinfo[1:2])
df.gene_go.bac.post <- as.data.frame(gene_go.bac.post@result)
write.csv(df.gene_go.bac.post,file = "/data/SongyiMetagenome/7_annotation/ProteinClustering/bacteria/7_KOGO/DEG/Bacteria_deg_post_GO_result.csv", row.names=F)



load(file="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/common.enriched.active.fun.RData")
load(file="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/common.enriched.pre.fun.RData")
load(file="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/common.enriched.post.fun.RData")

KOgene_list<-common.enriched.active.fun
kegg_result.fun<-enricher(KOgene_list,
                          TERM2GENE=KOannotation.fun[c(2,1)],
                          TERM2NAME=KOannotation.fun[c(2,4)])
df.kegg_result.fun <- as.data.frame(kegg_result.fun@result)
write.csv(df.kegg_result.fun,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_deg_active_kegg_result.csv")

KOgene_list<-common.enriched.pre.fun
kegg_result.fun<-enricher(KOgene_list,
                          TERM2GENE=KOannotation.fun[c(2,1)],
                          TERM2NAME=KOannotation.fun[c(2,4)])
df.kegg_result.fun <- as.data.frame(kegg_result.fun@result)
write.csv(df.kegg_result.fun,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_deg_pre_kegg_result.csv")

KOgene_list<-common.enriched.post.fun
kegg_result.fun<-enricher(KOgene_list,
                          TERM2GENE=KOannotation.fun[c(2,1)],
                          TERM2NAME=KOannotation.fun[c(2,4)])
df.kegg_result.fun <- as.data.frame(kegg_result.fun@result)
write.csv(df.kegg_result.fun,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_deg_post_kegg_result.csv")


KOgene_list<-common.enriched.active.fun
keggpath_result.fun.active<-enricher(KOgene_list,
                                     TERM2GENE=KOannotation.fun[c(3,1)],
                                     TERM2NAME=KOannotation.fun[c(3,4)])
df.keggpath_result.fun.active <- as.data.frame(keggpath_result.fun.active@result)
write.csv(df.keggpath_result.fun.active,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_active_kegg_path_result.csv")

KOgene_list<-common.enriched.pre.fun
keggpath_result.fun.pre<-enricher(KOgene_list,
                                  TERM2GENE=KOannotation.fun[c(3,1)],
                                  TERM2NAME=KOannotation.fun[c(3,4)])
df.keggpath_result.fun.pre <- as.data.frame(keggpath_result.fun.pre@result)
write.csv(df.keggpath_result.fun.pre,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_pre_kegg_path_result.csv")

KOgene_list<-common.enriched.post.fun
keggpath_result.fun.post<-enricher(KOgene_list,
                                   TERM2GENE=KOannotation.fun[c(3,1)],
                                   TERM2NAME=KOannotation.fun[c(3,4)])
df.keggpath_result.fun.post <- as.data.frame(keggpath_result.fun.post@result)
write.csv(df.keggpath_result.fun.post,file ="/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_post_kegg_path_result.csv")



GOgene_list<-common.enriched.active.fun
gene_go.fun.active<-enricher(GOgene_list,
                             TERM2GENE=GOannotation.fun[c(2,1)],
                             TERM2NAME=GOinfo[1:2])
df.gene_go.fun.active <- as.data.frame(gene_go.fun.active@result)
write.csv(df.gene_go.fun.active,file = "/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_deg_active_GO_result.csv", row.names=F)

GOgene_list<-common.enriched.pre.fun
gene_go.fun.pre<-enricher(GOgene_list,
                          TERM2GENE=GOannotation.fun[c(2,1)],
                          TERM2NAME=GOinfo[1:2])
df.gene_go.fun.pre <- as.data.frame(gene_go.fun.pre@result)
write.csv(df.gene_go.fun.pre,file = "/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_deg_pre_GO_result.csv", row.names=F)


GOgene_list<-common.enriched.post.fun
gene_go.fun.post<-enricher(GOgene_list,
                           TERM2GENE=GOannotation.fun[c(2,1)],
                           TERM2NAME=GOinfo[1:2])
df.gene_go.fun.post <- as.data.frame(gene_go.fun.post@result)
write.csv(df.gene_go.fun.post,file = "/data/SongyiMetagenome/7_annotation/ProteinClustering/fungi/7_KOGO/DEG/Fungi_deg_post_GO_result.csv", row.names=F)


####### 2. Extract significant GO terms belonging to biological process #######
##### 2-1. load GO info table
Go.tab<-read.delim("~/Desktop/SNU/SNU_BaeInHyup/Metagenome/231026_2023Data/FunctionalAnnotation/EggNogRes/GeneOntology/go.tb",sep="\t",header = T)
head(Go.tab)

##### 2-2. Extract significant GO terms
df.gene_go.bac.active<-read.csv("./GeneEnrichment/Bacteria/DEG/Bacteria_deg_active_GO_result.csv")
df.gene_go.bac.active<-df.gene_go.bac.active[-c(8)]
names(df.gene_go.bac.active)[1] <- "GO"
df.gene_go.bac.active.sig <-subset(df.gene_go.bac.active, qvalue < 0.05)
df.gene_go.bac.active.sig.merged<-merge(Go.tab,df.gene_go.bac.active.sig,by=c("GO"="GO","Description"="Description"))
nrow(df.gene_go.bac.active.sig)
df.gene_go.bac.active.sig[!(df.gene_go.bac.active.sig$GO%in%df.gene_go.bac.active.sig.merged$GO),]
df.gene_go.bac.active.sig.BP <- subset(df.gene_go.bac.active.sig.merged, level == "BP")
writexl::write_xlsx(df.gene_go.bac.active.sig,"./GeneEnrichment/Bacteria/DEG/Bacteria_deg_active_GO_result_sig.xlsx")


df.gene_go.fun.active<-read.csv("./GeneEnrichment/Fungi/DEG/Fungi_deg_active_GO_result.csv")
df.gene_go.fun.active<-df.gene_go.fun.active[-c(8)]
names(df.gene_go.fun.active)[1] <- "GO"
df.gene_go.fun.active.sig <-subset(df.gene_go.fun.active, qvalue < 0.05)
df.gene_go.fun.active.sig.merged<-merge(Go.tab,df.gene_go.fun.active.sig,by=c("GO"="GO","Description"="Description"))
nrow(df.gene_go.fun.active.sig)
df.gene_go.fun.active.sig[!(df.gene_go.fun.active.sig$GO%in%df.gene_go.fun.active.sig.merged$GO),]
df.gene_go.fun.active.sig.woMF<- subset(df.gene_go.fun.active.sig.merged, level != "MF")
writexl::write_xlsx(df.gene_go.fun.active.sig,"./GeneEnrichment/Fungi/DEG/Fungi_deg_active_GO_result_sig.xlsx")



df.gene_go.bac.pre<-read.csv("./GeneEnrichment/Bacteria/DEG/Bacteria_deg_pre_GO_result.csv")
df.gene_go.bac.pre<-df.gene_go.bac.pre[-c(8)]
names(df.gene_go.bac.pre)[1] <- "GO"
df.gene_go.bac.pre.sig <-subset(df.gene_go.bac.pre, qvalue < 0.05)
df.gene_go.bac.pre.sig.merged<-merge(Go.tab,df.gene_go.bac.pre.sig,by=c("GO"="GO","Description"="Description"))
nrow(df.gene_go.bac.pre.sig)
df.gene_go.bac.pre.sig[!(df.gene_go.bac.pre.sig$GO%in%df.gene_go.bac.pre.sig.merged$GO),]
df.gene_go.bac.pre.sig.BP <- subset(df.gene_go.bac.pre.sig.merged, level == "BP")



ggplot(df.gene_go.bac.active.sig.BP, aes(x=Count, y=GO, color = Description))+
  geom_point(aes(size = -log(qvalue)),shape=16, alpha = 0.8)+ 
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
  theme_classic()




####### 3. Extract significant KEGG orthologs #######
##### 3-1. Load a KEGG table
KOlist <- readxl::read_xlsx("./GeneAbundance/KOlist.xlsx",1)
##### 3-2. Load KEGG ortholog over-representation results
df.gene_kegg.bac.active <- read.csv("./GeneEnrichment/Bacteria/DEG/Bacteria_deg_active_kegg_result.csv")
df.gene_kegg.bac.active <- df.gene_kegg.bac.active[-c(1,9)]
head(df.gene_kegg.bac.active)

df.gene_kegg.bac.pre <- read.csv("./GeneEnrichment/Bacteria/DEG/Bacteria_deg_pre_kegg_result.csv")
df.gene_kegg.bac.pre <- df.gene_kegg.bac.pre[-c(1,9)]
head(df.gene_kegg.bac.pre)

df.gene_kegg.bac.post <- read.csv("./GeneEnrichment/Bacteria/DEG/Bacteria_deg_post_kegg_result.csv")
df.gene_kegg.bac.post <- df.gene_kegg.bac.post[-c(1,9)]
head(df.gene_kegg.bac.post)

df.gene_kegg.fun.active <- read.csv("./GeneEnrichment/Fungi/DEG/Fungi_deg_active_kegg_result.csv")
df.gene_kegg.fun.active <- df.gene_kegg.fun.active[-c(1,9)]
head(df.gene_kegg.fun.active)

df.gene_kegg.fun.pre <- read.csv("./GeneEnrichment/Fungi/DEG/Fungi_deg_pre_kegg_result.csv")
df.gene_kegg.fun.pre <- df.gene_kegg.fun.pre[-c(1,9)]
head(df.gene_kegg.fun.pre)

df.gene_kegg.fun.post <- read.csv("./GeneEnrichment/Fungi/DEG/Fungi_deg_post_kegg_result.csv")
df.gene_kegg.fun.post <- df.gene_kegg.fun.post[-c(1,9)]
head(df.gene_kegg.fun.post)


##### 3-3. Assign KO descriptions
names(df.gene_kegg.bac.active)[1] <- "KO"
names(df.gene_kegg.bac.pre)[1] <- "KO"
names(df.gene_kegg.bac.post)[1] <- "KO"

df.gene_kegg.bac.active.info<-merge(KOlist, df.gene_kegg.bac.active, by = "KO")
df.gene_kegg.bac.active.info.sig <- df.gene_kegg.bac.active.info[df.gene_kegg.bac.active.info$qvalue <0.05,]
df.gene_kegg.bac.active.info.sig <- df.gene_kegg.bac.active.info.sig %>% group_by(KO)%>% 
                                          arrange(desc(Description),Count)

unique(df.gene_kegg.bac.active.info.sig$Description)

df.gene_kegg.bac.pre.info<-merge(KOlist, df.gene_kegg.bac.pre, by = "KO")
df.gene_kegg.bac.pre.info.sig <- df.gene_kegg.bac.pre.info[df.gene_kegg.bac.pre.info$qvalue <0.05,]
df.gene_kegg.bac.pre.info.sig <- df.gene_kegg.bac.pre.info.sig %>% group_by(KO)%>% 
  arrange(desc(Description),Count)

unique(df.gene_kegg.bac.pre.info.sig$Description)

df.gene_kegg.bac.post.info<-merge(KOlist, df.gene_kegg.bac.post, by = "KO")
df.gene_kegg.bac.post.info.sig <- df.gene_kegg.bac.post.info[df.gene_kegg.bac.post.info$qvalue <0.05,]
df.gene_kegg.bac.post.info.sig <- df.gene_kegg.bac.post.info.sig %>% group_by(KO)%>% 
  arrange(desc(Description),Count)

unique(df.gene_kegg.bac.post.info.sig$Description)


names(df.gene_kegg.fun.active)[1] <- "KO"
names(df.gene_kegg.fun.pre)[1] <- "KO"
names(df.gene_kegg.fun.post)[1] <- "KO"

df.gene_kegg.fun.active.info<-merge(KOlist, df.gene_kegg.fun.active, by = "KO")
df.gene_kegg.fun.active.info.sig <- df.gene_kegg.fun.active.info[df.gene_kegg.fun.active.info$qvalue <0.05,]
df.gene_kegg.fun.active.info.sig <- df.gene_kegg.fun.active.info.sig %>% group_by(KO)%>% 
  arrange(desc(Description),Count)

unique(df.gene_kegg.fun.active.info.sig$Description)

df.gene_kegg.fun.pre.info<-merge(KOlist, df.gene_kegg.fun.pre, by = "KO")
df.gene_kegg.fun.pre.info.sig <- df.gene_kegg.fun.pre.info[df.gene_kegg.fun.pre.info$qvalue <0.05,]
df.gene_kegg.fun.pre.info.sig <- df.gene_kegg.fun.pre.info.sig %>% group_by(KO)%>% 
  arrange(desc(Description),Count)

unique(df.gene_kegg.fun.pre.info.sig$Description)

df.gene_kegg.fun.post.info<-merge(KOlist, df.gene_kegg.fun.post, by = "KO")
df.gene_kegg.fun.post.info.sig <- df.gene_kegg.fun.post.info[df.gene_kegg.fun.post.info$qvalue <0.05,]
df.gene_kegg.fun.post.info.sig <- df.gene_kegg.fun.post.info.sig %>% group_by(KO)%>% 
  arrange(desc(Description),Count)

unique(df.gene_kegg.fun.post.info.sig$Description)

####### 3-4. Select KOs involved in metabolism and nutrient cycling
kopath.list.bac.active<-c("Butanoate metabolism","Carotenoid biosynthesis","Glyoxylate and dicarboxylate metabolism",
  "Cysteine and methionine metabolism","Quorum sensing","Arginine and proline metabolism",
  "Glycolysis / Gluconeogenesis","Methane metabolism","ABC transporters",
  "Phenylalanine metabolism","Porphyrin metabolism","Inositol phosphate metabolism",
  "Polycyclic aromatic hydrocarbon degradation","Sulfur metabolism","Alanine, aspartate and glutamate metabolism",
  "Mineral absorption")

kopath.list.bac.pre<-c("Butanoate metabolism","Carotenoid biosynthesis","Glyoxylate and dicarboxylate metabolism",
                          "Cysteine and methionine metabolism","Quorum sensing","Arginine and proline metabolism",
                          "Glycolysis / Gluconeogenesis","Methane metabolism","ABC transporters",
                          "Phenylalanine metabolism","Porphyrin metabolism","Inositol phosphate metabolism",
                          "Polycyclic aromatic hydrocarbon degradation","Sulfur metabolism","Alanine, aspartate and glutamate metabolism",
                          "Mineral absorption","Polyketide sugar unit biosynthesis","Starch and sucrose metabolism",
                       "Terpenoid backbone biosynthesis")

kopath.list.bac.post<-c("Butanoate metabolism","Carotenoid biosynthesis","Glyoxylate and dicarboxylate metabolism",
                       "Cysteine and methionine metabolism","Quorum sensing","Arginine and proline metabolism",
                       "Glycolysis / Gluconeogenesis","Methane metabolism","ABC transporters",
                       "Phenylalanine metabolism","Porphyrin metabolism","Inositol phosphate metabolism",
                       "Polycyclic aromatic hydrocarbon degradation","Sulfur metabolism","Alanine, aspartate and glutamate metabolism",
                       "Mineral absorption","Polyketide sugar unit biosynthesis","Starch and sucrose metabolism",
                       "Terpenoid backbone biosynthesis","Riboflavin metabolism","Fructose and mannose metabolism",
                       "Biosynthesis of various other secondary metabolites")


############# Make a plot (KO level) - Fig. 4b, Fig. S5 ##########

kegg.path.color <-c(
  "Butanoate metabolism" = "#ADD8E6",
  "Carotenoid biosynthesis" = "#FFA500",
  "Glyoxylate and dicarboxylate metabolism" = "#90EE90",
  "Cysteine and methionine metabolism" = "#FFFFE0",
  "Quorum sensing" = "#800080",
  "Arginine and proline metabolism" = "#00FFFF",
  "Glycolysis / Gluconeogenesis" = "#8B0000",
  "Methane metabolism" = "#808080",
  "ABC transporters" = "#006400",
  "Phenylalanine metabolism" = "#FFB6C1",
  "Porphyrin metabolism" = "#008080",
  "Inositol phosphate metabolism" = "#D8BFD8",
  "Polycyclic aromatic hydrocarbon degradation" = "#A52A2A",
  "Sulfur metabolism" = "#FFD700",
  "Alanine, aspartate and glutamate metabolism" = "#FA8072",
  "Mineral absorption" = "#00008B",
  "Polyketide sugar unit biosynthesis" = "#EE82EE",
  "Starch and sucrose metabolism" = "#32CD32",
  "Terpenoid backbone biosynthesis" = "#FFFF00",
  "Riboflavin metabolism" = "#FFA07A",
  "Fructose and mannose metabolism" = "#800020",
  "Biosynthesis of various other secondary metabolites" = "#40E0D0"
)

df.gene_kegg.bac.active.sig.plot<-subset(df.gene_kegg.bac.active.info.sig, Description %in% kopath.list.bac.active)
df.gene_kegg.bac.active.sig.plot$KO <- factor(df.gene_kegg.bac.active.sig.plot$KO, levels = df.gene_kegg.bac.active.sig.plot$KO)

ggplot(df.gene_kegg.bac.active.sig.plot, aes(x=Count, y=KO, color = Description))+
  geom_point(aes(size = -log(qvalue)),shape=16, alpha = 0.8)+ 
  #coord_flip()+ 
  # theme(strip.placement = "outside",
  #       strip.background = element_rect(fill=NA,colour=NA),
  #       panel.spacing=unit(0,"cm"), axis.title.y = element_blank()) +
  # annotation_custom(grob = linesGrob(), xmin = -0.75, xmax = -0.75, ymin = -3.25, ymax = -0.75) +
  # coord_cartesian(clip="off")+
  guides(color = guide_legend(reverse =F, ncol = 1))+
  theme(legend.position="right") + #theme(aspect.ratio = 3)+
scale_color_manual(values = kegg.path.color)+
#theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) +
theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme_classic()

######### Composition of KO and GOs of differentially expressed genes ########


KOannotation.bac <- read.table("../ProteinClustering/EggNog/Bacteria/KOannotation.tsv", sep ="\t", header =T)
KOannotation.fun <- read.table("../ProteinClustering/EggNog/Fungi/KOannotation.tsv", sep ="\t", header =T)

GOannotation.bac <- read.table("../ProteinClustering/EggNog/Bacteria/GOannotation.tsv", sep ="\t", header =T)
GOannotation.fun <- read.table("../ProteinClustering/EggNog/Fungi/GOannotation.tsv", sep ="\t", header =T)

KO.deg.active.bac<-subset(KOannotation.bac, gene %in% common.enriched.active.bac)
KO.deg.active.bac$count <- 1
KO.deg.active.bac.ko.summary <- KO.deg.active.bac %>% group_by(KO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
KO.deg.active.bac.kopath.summary <- KO.deg.active.bac %>% group_by(description) %>% reframe(count = sum(count)) %>% arrange(desc(count))

KO.deg.pre.bac<-subset(KOannotation.bac, gene %in% common.enriched.pre.bac)
KO.deg.pre.bac$count <- 1
KO.deg.pre.bac.ko.summary <- KO.deg.pre.bac %>% group_by(KO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
KO.deg.pre.bac.kopath.summary <- KO.deg.pre.bac %>% group_by(description) %>% reframe(count = sum(count)) %>% arrange(desc(count))

KO.deg.post.bac<-subset(KOannotation.bac, gene %in% common.enriched.post.bac)
KO.deg.post.bac$count <- 1
KO.deg.post.bac.ko.summary <- KO.deg.post.bac %>% group_by(KO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
KO.deg.post.bac.kopath.summary <- KO.deg.post.bac %>% group_by(description) %>% reframe(count = sum(count)) %>% arrange(desc(count))

GO.deg.active.bac<-subset(GOannotation.bac, gene %in% common.enriched.active.bac)
GO.deg.active.bac$count <- 1
GO.deg.active.bac.GO.summary <- GO.deg.active.bac %>% group_by(GO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
GO.deg.active.bac.GO.summary<-merge(Go.tab, GO.deg.active.bac.GO.summary, by ="GO")
GO.deg.active.bac.GO.summary[GO.deg.active.bac.GO.summary$level =="BP",]$Description
GO.deg.active.bac.GO.summary$TmStatus <- "Active"

GO.deg.pre.bac<-subset(GOannotation.bac, gene %in% common.enriched.pre.bac)
GO.deg.pre.bac$count <- 1
GO.deg.pre.bac.GO.summary <- GO.deg.pre.bac %>% group_by(GO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
GO.deg.pre.bac.GO.summary<-merge(Go.tab, GO.deg.pre.bac.GO.summary, by ="GO")
GO.deg.pre.bac.GO.summary[GO.deg.pre.bac.GO.summary$level =="BP",]$Description
GO.deg.pre.bac.GO.summary$TmStatus <- "Pre"

GO.deg.post.bac<-subset(GOannotation.bac, gene %in% common.enriched.post.bac)
GO.deg.post.bac$count <- 1
GO.deg.post.bac.GO.summary <- GO.deg.post.bac %>% group_by(GO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
GO.deg.post.bac.GO.summary<-merge(Go.tab, GO.deg.post.bac.GO.summary, by ="GO")
GO.deg.post.bac.GO.summary[GO.deg.post.bac.GO.summary$level =="BP",]$Description
GO.deg.post.bac.GO.summary$TmStatus <- "Post"


KO.deg.bac.kopath <- rbind(KO.deg.active.bac.kopath.summary,KO.deg.pre.bac.kopath.summary,KO.deg.post.bac.kopath.summary)
KO.deg.bac.kopath.wide <- reshape2::dcast(KO.deg.bac.kopath, description~TmStatus, value.var = "count")
KO.deg.bac.kopath.wide[is.na(KO.deg.bac.kopath.wide)] <- 0
head(KO.deg.bac.kopath.wide)


KO.deg.active.fun<-subset(KOannotation.fun, gene %in% common.enriched.active.fun)
KO.deg.active.fun$count <- 1
KO.deg.active.fun.ko.summary <- KO.deg.active.fun %>% group_by(KO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
KO.deg.active.fun.kopath.summary <- KO.deg.active.fun %>% group_by(description) %>% reframe(count = sum(count)) %>% arrange(desc(count))

KO.deg.pre.fun<-subset(KOannotation.fun, gene %in% common.enriched.pre.fun)
KO.deg.pre.fun$count <- 1
KO.deg.pre.fun.ko.summary <- KO.deg.pre.fun %>% group_by(KO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
KO.deg.pre.fun.kopath.summary <- KO.deg.pre.fun %>% group_by(description) %>% reframe(count = sum(count)) %>% arrange(desc(count))

KO.deg.post.fun<-subset(KOannotation.fun, gene %in% common.enriched.post.fun)
KO.deg.post.fun$count <- 1
KO.deg.post.fun.ko.summary <- KO.deg.post.fun %>% group_by(KO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
KO.deg.post.fun.kopath.summary <- KO.deg.post.fun %>% group_by(description) %>% reframe(count = sum(count)) %>% arrange(desc(count))


GO.deg.active.fun<-subset(GOannotation.fun, gene %in% common.enriched.active.fun)
GO.deg.active.fun$count <- 1
GO.deg.active.fun.GO.summary <- GO.deg.active.fun %>% group_by(GO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
GO.deg.active.fun.GO.summary<-merge(Go.tab, GO.deg.active.fun.GO.summary, by ="GO")
GO.deg.active.fun.GO.summary[GO.deg.active.fun.GO.summary$level =="BP",]$Description
GO.deg.active.fun.GO.summary$TmStatus <- "Active"

GO.deg.pre.fun<-subset(GOannotation.fun, gene %in% common.enriched.pre.fun)
GO.deg.pre.fun$count <- 1
GO.deg.pre.fun.GO.summary <- GO.deg.pre.fun %>% group_by(GO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
GO.deg.pre.fun.GO.summary<-merge(Go.tab, GO.deg.pre.fun.GO.summary, by ="GO")
GO.deg.pre.fun.GO.summary[GO.deg.pre.fun.GO.summary$level =="BP",]$Description
GO.deg.pre.fun.GO.summary$TmStatus <- "Pre"

GO.deg.post.fun<-subset(GOannotation.fun, gene %in% common.enriched.post.fun)
GO.deg.post.fun$count <- 1
GO.deg.post.fun.GO.summary <- GO.deg.post.fun %>% group_by(GO) %>% reframe(count = sum(count)) %>% arrange(desc(count))
GO.deg.post.fun.GO.summary<-merge(Go.tab, GO.deg.post.fun.GO.summary, by ="GO")
GO.deg.post.fun.GO.summary[GO.deg.post.fun.GO.summary$level =="BP",]$Description
GO.deg.post.fun.GO.summary$TmStatus <- "Post"



KO.deg.active.bac.kopath.summary$TmStatus <- "Active"
KO.deg.pre.bac.kopath.summary$TmStatus <- "Pre"
KO.deg.post.bac.kopath.summary$TmStatus <- "Post"
KO.deg.bac.kopath <- rbind(KO.deg.active.bac.kopath.summary,KO.deg.pre.bac.kopath.summary,KO.deg.post.bac.kopath.summary)
KO.deg.bac.kopath.wide <- reshape2::dcast(KO.deg.bac.kopath, description~TmStatus, value.var = "count")
KO.deg.bac.kopath.wide[is.na(KO.deg.bac.kopath.wide)] <- 0
head(KO.deg.bac.kopath.wide)

KO.deg.active.fun.kopath.summary$TmStatus <- "Active"
KO.deg.pre.fun.kopath.summary$TmStatus <- "Pre"
KO.deg.post.fun.kopath.summary$TmStatus <- "Post"
KO.deg.fun.kopath <- rbind(KO.deg.active.fun.kopath.summary,KO.deg.pre.fun.kopath.summary,KO.deg.post.fun.kopath.summary)
KO.deg.fun.kopath.wide <- reshape2::dcast(KO.deg.fun.kopath, description~TmStatus, value.var = "count")
KO.deg.fun.kopath.wide[is.na(KO.deg.fun.kopath.wide)] <- 0
head(KO.deg.fun.kopath.wide)

KO.deg.bac.kopath.wide <- KO.deg.bac.kopath.wide[c(1,4,2,3)]
KO.deg.fun.kopath.wide <- KO.deg.fun.kopath.wide[c(1,4,2,3)]

KO.deg.bac.kopath.wide[(grepl("metabolism",KO.deg.bac.kopath.wide$description)),]
KO.deg.fun.kopath.wide[(grepl("metabolism",KO.deg.fun.kopath.wide$description)),]

KO.deg.bac.kopath.wide[(grepl("MAPK signaling pathway",KO.deg.bac.kopath.wide$description)),]
KO.deg.fun.kopath.wide[(grepl("MAPK signaling pathway",KO.deg.fun.kopath.wide$description)),]

GO.deg.bac.GO <- rbind(GO.deg.active.bac.GO.summary,GO.deg.pre.bac.GO.summary,GO.deg.post.bac.GO.summary)
GO.deg.bac.GO.wide <- reshape2::dcast(GO.deg.bac.GO, GO+Description~TmStatus, value.var = "count")
GO.deg.bac.GO.wide[is.na(GO.deg.bac.GO.wide)] <- 0
head(GO.deg.bac.GO.wide)

GO.deg.fun.GO <- rbind(GO.deg.active.fun.GO.summary,GO.deg.pre.fun.GO.summary,GO.deg.post.fun.GO.summary)
GO.deg.fun.GO.wide <- reshape2::dcast(GO.deg.fun.GO, GO+Description~TmStatus, value.var = "count")
GO.deg.fun.GO.wide[is.na(GO.deg.fun.GO.wide)] <- 0
head(GO.deg.fun.GO.wide)

GO.deg.bac.GO.wide <- GO.deg.bac.GO.wide[c(1,2,5,3,4)]
GO.deg.fun.GO.wide <- GO.deg.fun.GO.wide[c(1,2,5,3,4)]

GO.deg.bac.GO.wide[(grepl("nitrogen",GO.deg.bac.GO.wide$Description)),]
GO.deg.fun.GO.wide[(grepl("interaction",GO.deg.fun.GO.wide$Description)),]

GO.deg.bac.GO.wide[(grepl("MAPK signaling pathway",GO.deg.bac.GO.wide$description)),]
GO.deg.fun.GO.wide[(grepl("MAPK signaling pathway",GO.deg.fun.GO.wide$description)),]

writexl::write_xlsx(GO.deg.bac.GO.wide,"./GeneExpression/Expression/GOs of Bacterial DEGs.xlsx")
writexl::write_xlsx(GO.deg.fun.GO.wide,"./GeneExpression/Expression/GOs of Fungal DEGs.xlsx")

writexl::write_xlsx(KO.deg.bac.kopath.wide,"./GeneExpression/Expression/KEGG pathways of Bacterial DEGs.xlsx")
writexl::write_xlsx(KO.deg.fun.kopath.wide,"./GeneExpression/Expression/KEGG pathways of Fungal DEGs.xlsx")

######### Fig. 4a and 4b - Bubble plots for up-regulated genes ######
###### Fig. 4a. Top 15 KEGG pathway
KO.deg.bac.kopath.wide$Domain <- "Bacteria"
KO.deg.fun.kopath.wide$Domain <- "Fungi"
kopath.degs<-rbind(KO.deg.bac.kopath.wide,KO.deg.fun.kopath.wide)

writexl::write_xlsx(kopath.degs,"./GeneExpression/Expression/DEG_KEGG_path_affiliation.xlsx")


kopath.degs.plot <- reshape2::melt(kopath.degs)
head(kopath.degs.plot)
names(kopath.degs.plot)[c(3,4)] <- c("TmStatus","Count")

kopath.list<-unique(kopath.degs.plot$description)
kopath.list.target1 <- kopath.list[grepl("metabolism|Metabolism", kopath.list)] 
kopath.list.target2 <- kopath.list[grepl("transport|Transport", kopath.list)] 
kopath.list.target3 <- kopath.list[grepl("Two-component|Type I|Type II|aromatic compounds|Carbon fixation|Biosynthesis of ", kopath.list)] 
kopath.list.targets<-c(kopath.list.target1,kopath.list.target2,kopath.list.target3)

kopath.degs.plot.sub <- subset(kopath.degs.plot, description %in%kopath.list.targets)
kopath.degs.plot.sub <- subset(kopath.degs.plot.sub, !(grepl("diabetes", kopath.degs.plot.sub$description)))
kopath.degs.plot.sub.summary <- kopath.degs.plot.sub %>% group_by(description,Domain) %>% reframe(Count =sum(Count))
kopath.degs.plot.sub.summary.bac.top20 <- kopath.degs.plot.sub.summary %>% filter(Domain == "Bacteria") %>% arrange(desc(Count)) %>% head(20)
kopath.degs.plot.sub.summary.fun.top20 <- kopath.degs.plot.sub.summary %>% filter(Domain == "Fungi") %>% arrange(desc(Count)) %>% head(20)

kopath.degs.plot.sub.top20.bac <- subset(kopath.degs.plot.sub, description %in%kopath.degs.plot.sub.summary.bac.top20$description & Domain =="Bacteria")
kopath.degs.plot.sub.top20.fun <- subset(kopath.degs.plot.sub, description %in%kopath.degs.plot.sub.summary.fun.top20$description & Domain =="Fungi")

ggplot(data.frame(kopath.degs.plot.sub.top20.bac), aes(x=TmStatus, y=description, color = TmStatus))+
  geom_point(aes(size = Count),shape=16, alpha = 0.8)+ facet_wrap(~Domain, scales = "free_y")+
  #coord_flip()+ 
  # theme(strip.placement = "outside",
  #       strip.background = element_rect(fill=NA,colour=NA),
  #       panel.spacing=unit(0,"cm"), axis.title.y = element_blank()) +
  # annotation_custom(grob = linesGrob(), xmin = -0.75, xmax = -0.75, ymin = -3.25, ymax = -0.75) +
  # coord_cartesian(clip="off")+
  guides(color = guide_legend(reverse =F, ncol = 1))+
  theme(legend.position="right") + #theme(aspect.ratio = 3)+
  scale_color_manual(values = c("Pre" = "khaki3","Active" = "dodgerblue4","Post"="salmon3"))+
  #theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) +
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme_classic()+theme(aspect.ratio = 2)

ggplot(data.frame(kopath.degs.plot.sub.top20.fun), aes(x=TmStatus, y=description, color = TmStatus))+
  geom_point(aes(size = Count),shape=16, alpha = 0.8)+ facet_wrap(~Domain, scales = "free_y")+
  #coord_flip()+ 
  # theme(strip.placement = "outside",
  #       strip.background = element_rect(fill=NA,colour=NA),
  #       panel.spacing=unit(0,"cm"), axis.title.y = element_blank()) +
  # annotation_custom(grob = linesGrob(), xmin = -0.75, xmax = -0.75, ymin = -3.25, ymax = -0.75) +
  # coord_cartesian(clip="off")+
  guides(color = guide_legend(reverse =F, ncol = 1))+
  theme(legend.position="right") + #theme(aspect.ratio = 3)+
  scale_color_manual(values = c("Pre" = "khaki3","Active" = "dodgerblue4","Post"="salmon3"))+
  #theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) +
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme_classic()+theme(aspect.ratio = 2)


###### Data S8 GO term part
GO.deg.bac.GO.wide$Domain <- "Bacteria"
GO.deg.fun.GO.wide$Domain <- "Fungi"
GO.degs<-rbind(GO.deg.bac.GO.wide,GO.deg.fun.GO.wide)

GOs.BP<-Go.tab$Description[which(Go.tab$level == "BP")]
GO.degs.BP<-subset(GO.degs, Description %in%GOs.BP)
GO.degs.plot <- reshape2::melt(GO.degs.BP)
GO.degs.BP.sub<- GO.degs.BP[!(grepl("obsolete", GO.degs.BP$Description)),]

writexl::write_xlsx(GO.degs.BP.sub,"./GeneExpression/Expression/DEG_GO_Biological_process_affiliation.xlsx")

#### Microbes possessing DEGs
#### 
head(diamond.rep.gene.bac.contig)
diamond.rep.gene.bac.contig$Phylum2 <- diamond.rep.gene.bac.contig$Phylum
diamond.rep.gene.bac.contig$Class2 <- diamond.rep.gene.bac.contig$Class
diamond.rep.gene.bac.contig$Order2 <- diamond.rep.gene.bac.contig$Order
diamond.rep.gene.bac.contig$Family2 <- diamond.rep.gene.bac.contig$Family
diamond.rep.gene.bac.contig$Genus2 <- diamond.rep.gene.bac.contig$Genus
diamond.rep.gene.bac.contig$Species2 <- diamond.rep.gene.bac.contig$Species

diamond.rep.gene.bac.contig$Phylum2[which(diamond.rep.gene.bac.contig$Phylum == "Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Domain[which(diamond.rep.gene.bac.contig$Phylum == "Unidentified")],"_unknown")
diamond.rep.gene.bac.contig$Class2[which(diamond.rep.gene.bac.contig$Class == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Phylum!="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Phylum[which(diamond.rep.gene.bac.contig$Class == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.bac.contig$Phylum!="Unidentified")],"_unknown")
diamond.rep.gene.bac.contig$Class2[which(diamond.rep.gene.bac.contig$Class == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Phylum2=="Bacteria_unknown")] <- "Bacteria_unknown"

diamond.rep.gene.bac.contig$Order2[which(diamond.rep.gene.bac.contig$Order == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Class!="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Class[which(diamond.rep.gene.bac.contig$Order == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.bac.contig$Class!="Unidentified")],"_unknown")

diamond.rep.gene.bac.contig$Order2[which(diamond.rep.gene.bac.contig$Order == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Class=="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Class2[which(diamond.rep.gene.bac.contig$Order == "Unidentified"& 
                                                                                                                                                   diamond.rep.gene.bac.contig$Class=="Unidentified")],"_unknown")


diamond.rep.gene.bac.contig$Order2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.bac.contig$Order2)

diamond.rep.gene.bac.contig$Family2[which(diamond.rep.gene.bac.contig$Family == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Order!="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Order[which(diamond.rep.gene.bac.contig$Family == "Unidentified"& 
                                                                                                                                                   diamond.rep.gene.bac.contig$Order!="Unidentified")],"_unknown")


diamond.rep.gene.bac.contig$Family2[which(diamond.rep.gene.bac.contig$Family == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Order=="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Order2[which(diamond.rep.gene.bac.contig$Family == "Unidentified"& 
                                                                                                                                                    diamond.rep.gene.bac.contig$Order=="Unidentified")],"_unknown")


diamond.rep.gene.bac.contig$Family2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.bac.contig$Family2)


diamond.rep.gene.bac.contig$Genus2[which(diamond.rep.gene.bac.contig$Genus == "Unidentified" & 
                                            diamond.rep.gene.bac.contig$Family!="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Family[which(diamond.rep.gene.bac.contig$Genus == "Unidentified"& 
                                                                                                                                                    diamond.rep.gene.bac.contig$Family!="Unidentified")],"_unknown")

diamond.rep.gene.bac.contig$Genus2[which(diamond.rep.gene.bac.contig$Genus == "Unidentified" & 
                                            diamond.rep.gene.bac.contig$Family=="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Family2[which(diamond.rep.gene.bac.contig$Genus == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.bac.contig$Family=="Unidentified")],"_unknown")


diamond.rep.gene.bac.contig$Genus2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.bac.contig$Genus2)

diamond.rep.gene.bac.contig$Species2[which(diamond.rep.gene.bac.contig$Species == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Genus!="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Genus[which(diamond.rep.gene.bac.contig$Species == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.bac.contig$Genus!="Unidentified")],"_unknown")

diamond.rep.gene.bac.contig$Species2[which(diamond.rep.gene.bac.contig$Species == "Unidentified" & 
                                           diamond.rep.gene.bac.contig$Genus=="Unidentified")] <- paste0(diamond.rep.gene.bac.contig$Genus2[which(diamond.rep.gene.bac.contig$Species == "Unidentified"& 
                                                                                                                                                      diamond.rep.gene.bac.contig$Genus=="Unidentified")],"_unknown")


diamond.rep.gene.bac.contig$Species2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.bac.contig$Species2)

### Fungi
diamond.rep.gene.fun.contig$Phylum2 <- diamond.rep.gene.fun.contig$Phylum
diamond.rep.gene.fun.contig$Class2 <- diamond.rep.gene.fun.contig$Class
diamond.rep.gene.fun.contig$Order2 <- diamond.rep.gene.fun.contig$Order
diamond.rep.gene.fun.contig$Family2 <- diamond.rep.gene.fun.contig$Family
diamond.rep.gene.fun.contig$Genus2 <- diamond.rep.gene.fun.contig$Genus
diamond.rep.gene.fun.contig$Species2 <- diamond.rep.gene.fun.contig$Species

diamond.rep.gene.fun.contig$Phylum2[which(diamond.rep.gene.fun.contig$Phylum == "Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Kingdom[which(diamond.rep.gene.fun.contig$Phylum == "Unidentified")],"_unknown")
diamond.rep.gene.fun.contig$Class2[which(diamond.rep.gene.fun.contig$Class == "Unidentified" & 
                                           diamond.rep.gene.fun.contig$Phylum!="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Phylum[which(diamond.rep.gene.fun.contig$Class == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.fun.contig$Phylum!="Unidentified")],"_unknown")

diamond.rep.gene.fun.contig$Class2[which(diamond.rep.gene.fun.contig$Class == "Unidentified" & 
                                           diamond.rep.gene.fun.contig$Phylum2=="Fungi_unknown")] <- "Fungi_unknown"

diamond.rep.gene.fun.contig$Order2[which(diamond.rep.gene.fun.contig$Order == "Unidentified" & 
                                           diamond.rep.gene.fun.contig$Class!="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Class[which(diamond.rep.gene.fun.contig$Order == "Unidentified"& 
                                                                                                                                                   diamond.rep.gene.fun.contig$Class!="Unidentified")],"_unknown")

diamond.rep.gene.fun.contig$Order2[which(diamond.rep.gene.fun.contig$Order == "Unidentified" & 
                                           diamond.rep.gene.fun.contig$Class=="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Class2[which(diamond.rep.gene.fun.contig$Order == "Unidentified"& 
                                                                                                                                                    diamond.rep.gene.fun.contig$Class=="Unidentified")],"_unknown")


diamond.rep.gene.fun.contig$Order2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.fun.contig$Order2)

diamond.rep.gene.fun.contig$Family2[which(diamond.rep.gene.fun.contig$Family == "Unidentified" & 
                                            diamond.rep.gene.fun.contig$Order!="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Order[which(diamond.rep.gene.fun.contig$Family == "Unidentified"& 
                                                                                                                                                    diamond.rep.gene.fun.contig$Order!="Unidentified")],"_unknown")


diamond.rep.gene.fun.contig$Family2[which(diamond.rep.gene.fun.contig$Family == "Unidentified" & 
                                            diamond.rep.gene.fun.contig$Order=="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Order2[which(diamond.rep.gene.fun.contig$Family == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.fun.contig$Order=="Unidentified")],"_unknown")


diamond.rep.gene.fun.contig$Family2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.fun.contig$Family2)


diamond.rep.gene.fun.contig$Genus2[which(diamond.rep.gene.fun.contig$Genus == "Unidentified" & 
                                           diamond.rep.gene.fun.contig$Family!="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Family[which(diamond.rep.gene.fun.contig$Genus == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.fun.contig$Family!="Unidentified")],"_unknown")

diamond.rep.gene.fun.contig$Genus2[which(diamond.rep.gene.fun.contig$Genus == "Unidentified" & 
                                           diamond.rep.gene.fun.contig$Family=="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Family2[which(diamond.rep.gene.fun.contig$Genus == "Unidentified"& 
                                                                                                                                                      diamond.rep.gene.fun.contig$Family=="Unidentified")],"_unknown")


diamond.rep.gene.fun.contig$Genus2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.fun.contig$Genus2)

diamond.rep.gene.fun.contig$Species2[which(diamond.rep.gene.fun.contig$Species == "Unidentified" & 
                                             diamond.rep.gene.fun.contig$Genus!="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Genus[which(diamond.rep.gene.fun.contig$Species == "Unidentified"& 
                                                                                                                                                     diamond.rep.gene.fun.contig$Genus!="Unidentified")],"_unknown")

diamond.rep.gene.fun.contig$Species2[which(diamond.rep.gene.fun.contig$Species == "Unidentified" & 
                                             diamond.rep.gene.fun.contig$Genus=="Unidentified")] <- paste0(diamond.rep.gene.fun.contig$Genus2[which(diamond.rep.gene.fun.contig$Species == "Unidentified"& 
                                                                                                                                                      diamond.rep.gene.fun.contig$Genus=="Unidentified")],"_unknown")


diamond.rep.gene.fun.contig$Species2 <- gsub("unknown_unknown","unknown",diamond.rep.gene.fun.contig$Species2)

names(diamond.rep.gene.bac.contig)
diamond.rep.gene.bac.contig.re <- diamond.rep.gene.bac.contig[c(1:6,13:18)]
diamond.rep.gene.fun.contig.re <- diamond.rep.gene.fun.contig[c(1:6,13:18)]
names(diamond.rep.gene.bac.contig.re)[c(7:12)] <-gsub("2","", names(diamond.rep.gene.bac.contig.re)[c(7:12)])
names(diamond.rep.gene.fun.contig.re)[c(7:12)] <-gsub("2","", names(diamond.rep.gene.fun.contig.re)[c(7:12)])



tax.gene.bac<-merge(rep.genes.bac.gtf.trim,diamond.rep.gene.bac.contig.re, by = c("contigName"="contigName"))
tax.gene.fun<-merge(rep.genes.fun.gtf.trim,diamond.rep.gene.fun.contig.re, by = c("contigName"="contigName"))

##### Overall DEGs
rep.bac.tax.active <- subset(tax.gene.bac, GeneID %in% common.enriched.active.bac)
rep.bac.tax.active
rep.bac.tax.active$count <- 1
rep.bac.tax.active.species <- rep.bac.tax.active%>%group_by(Genus)%>% reframe(count=sum(count))%>% arrange(desc(count))
rep.bac.tax.active.phylum <- rep.bac.tax.active%>%group_by(Phylum)%>% reframe(count=sum(count))%>% arrange(desc(count))

rep.bac.tax.pre <- subset(tax.gene.bac, GeneID %in% common.enriched.pre.bac)
rep.bac.tax.pre
rep.bac.tax.pre$count <- 1
rep.bac.tax.pre.species <- rep.bac.tax.pre%>%group_by(Phylum)%>% reframe(count=sum(count))%>% arrange(desc(count))

rep.bac.tax.post <- subset(tax.gene.bac, GeneID %in% common.enriched.post.bac)
rep.bac.tax.post
rep.bac.tax.post$count <- 1
rep.bac.tax.post.species <- rep.bac.tax.post%>%group_by(Phylum)%>% reframe(count=sum(count))%>% arrange(desc(count))

rep.fun.tax.active <- subset(tax.gene.fun, GeneID %in% common.enriched.active.fun)
rep.fun.tax.active
rep.fun.tax.active$count <- 1
rep.fun.tax.active.species <- rep.fun.tax.active%>%group_by(Species)%>% reframe(count=sum(count))%>% arrange(desc(count))

rep.fun.tax.pre <- subset(tax.gene.fun, GeneID %in% common.enriched.pre.fun)
rep.fun.tax.pre
rep.fun.tax.pre$count <- 1
rep.fun.tax.pre.species <- rep.fun.tax.pre%>%group_by(Species)%>% reframe(count=sum(count))%>% arrange(desc(count))

rep.fun.tax.post <- subset(tax.gene.fun, GeneID %in% common.enriched.post.fun)
rep.fun.tax.post
rep.fun.tax.post$count <- 1
rep.fun.tax.post.species <- rep.fun.tax.post%>%group_by(Genus)%>% reframe(count=sum(count))%>% arrange(desc(count))


######### Supplementary Figure - Taxonomic composition of DEGs (pie chart) ######
### Phylum level
rep.bac.tax.pre$TmStatus <- "Pre"
rep.bac.tax.active$TmStatus <- "Active"
rep.bac.tax.post$TmStatus <- "Post"

rep.bac.upregul.tax <- rbind(rep.bac.tax.pre,rep.bac.tax.active,rep.bac.tax.post)
rep.bac.upregul.tax.phylum <- rep.bac.upregul.tax%>%group_by(Phylum,TmStatus)%>% reframe(count=sum(count))
rep.bac.upregul.tax.phylum <- rep.bac.upregul.tax.phylum%>%group_by(TmStatus)%>% mutate(total=sum(count))
rep.bac.upregul.tax.phylum$ratio <- rep.bac.upregul.tax.phylum$count/rep.bac.upregul.tax.phylum$total


rep.bac.upregul.tax.phylum.summary <- rep.bac.upregul.tax.phylum%>%group_by(Phylum) %>% reframe(count = sum(count)) %>% arrange(desc(count))
rep.bac.upregul.tax.phylum.summary$ratio <- rep.bac.upregul.tax.phylum.summary$count/sum(rep.bac.upregul.tax.phylum.summary$count)
other.phylum.bac<-rep.bac.upregul.tax.phylum.summary$Phylum[which(rep.bac.upregul.tax.phylum.summary$ratio < 0.005)]

rep.bac.upregul.tax.phylum$Phylum2 <- rep.bac.upregul.tax.phylum$Phylum
rep.bac.upregul.tax.phylum$Phylum2[which(rep.bac.upregul.tax.phylum$Phylum %in% other.phylum.bac)] <- "Others"
unique(rep.bac.upregul.tax.phylum$Phylum2)


rep.bac.upregul.tax.phylum$Phylum2 <- factor(rep.bac.upregul.tax.phylum$Phylum2,levels=c("Pseudomonadota", "Actinomycetota",
                                                                                   "Acidobacteriota","Planctomycetota",
                                                                                   "Chloroflexota","Verrucomicrobiota",
                                                                                   "Myxococcota",
                                                                                   "Bacteria_unknown",
                                                                                   "Others"))

rep.bac.upregul.tax.phylum$TmStatus <- factor(rep.bac.upregul.tax.phylum$TmStatus,levels = c("Pre","Active","Post"))
rep.bac.upregul.tax.phylum.plot <- rep.bac.upregul.tax.phylum %>% group_by(Phylum2,TmStatus) %>% reframe(ratio = sum(ratio))

ggplot(rep.bac.upregul.tax.phylum.plot, aes(x='', y=ratio, fill=Phylum2))+
  geom_bar(stat='identity')+facet_wrap(~TmStatus)+
  theme_void()+
  coord_polar('y', start=0)+
  ggrepel::geom_text_repel(aes(label=paste0(round(ratio*100,2), '%')),
            position=position_stack(vjust=0.5), max.overlaps = 100)


rep.fun.tax.pre$TmStatus <- "Pre"
rep.fun.tax.active$TmStatus <- "Active"
rep.fun.tax.post$TmStatus <- "Post"

rep.fun.upregul.tax <- rbind(rep.fun.tax.pre,rep.fun.tax.active,rep.fun.tax.post)
rep.fun.upregul.tax.class <- rep.fun.upregul.tax%>%group_by(Phylum,Class,TmStatus)%>% reframe(count=sum(count))
rep.fun.upregul.tax.class <- rep.fun.upregul.tax.class%>%group_by(TmStatus)%>% mutate(total=sum(count))
rep.fun.upregul.tax.class$ratio <- rep.fun.upregul.tax.class$count/rep.fun.upregul.tax.class$total


rep.fun.upregul.tax.class.summary <- rep.fun.upregul.tax.class%>%group_by(Phylum,Class) %>% reframe(count = sum(count)) %>% arrange(desc(count))
rep.fun.upregul.tax.class.summary$ratio <- rep.fun.upregul.tax.class.summary$count/sum(rep.fun.upregul.tax.class.summary$count)
other.class.fun<-rep.fun.upregul.tax.class.summary$Class[which(rep.fun.upregul.tax.class.summary$ratio < 0.01)]

rep.fun.upregul.tax.class$Class2 <- rep.fun.upregul.tax.class$Class
rep.fun.upregul.tax.class$Class2[which(rep.fun.upregul.tax.class$Class %in% other.class.fun)] <- "Others"
unique(rep.fun.upregul.tax.class$Class2)


rep.fun.upregul.tax.class$Class2 <- factor(rep.fun.upregul.tax.class$Class2,levels=c("Agaricomycetes", "Dothideomycetes",
                                                                                         "Saccharomycetes","Leotiomycetes",
                                                                                         "Eurotiomycetes","Chytridiomycetes","Sordariomycetes",
                                                                                         "Mucoromycetes",
                                                                                         "Ascomycota_unknown",
                                                                                         "Fungi_unknown",
                                                                                         "Others"))

rep.fun.upregul.tax.class$TmStatus <- factor(rep.fun.upregul.tax.class$TmStatus,levels = c("Pre","Active","Post"))
rep.fun.upregul.tax.class.plot <- rep.fun.upregul.tax.class %>% group_by(Class2,TmStatus) %>% reframe(ratio = sum(ratio))
ggplot(rep.fun.upregul.tax.class.plot, aes(x='', y=ratio, fill=Class2))+
  geom_bar(stat='identity')+facet_wrap(~TmStatus)+
  theme_void()+
  coord_polar('y', start=0)+
  ggrepel::geom_text_repel(aes(label=paste0(round(ratio*100,2), '%')),
                           position=position_stack(vjust=0.5),max.overlaps=100)


rep.bac.upregul.tax.species <- rep.bac.upregul.tax%>%group_by(Phylum,Class,Order,Family,Genus,Species,TmStatus)%>% reframe(count=sum(count))
rep.bac.upregul.tax.species <- rep.bac.upregul.tax.species%>%group_by(TmStatus)%>% mutate(total=sum(count))
rep.bac.upregul.tax.species$ratio <- rep.bac.upregul.tax.species$count/rep.bac.upregul.tax.species$total

rep.fun.upregul.tax.species <- rep.fun.upregul.tax%>%group_by(Phylum,Class,Order,Family,Genus,Species,TmStatus)%>% reframe(count=sum(count))
rep.fun.upregul.tax.species <- rep.fun.upregul.tax.species%>%group_by(TmStatus)%>% mutate(total=sum(count))
rep.fun.upregul.tax.species$ratio <- rep.fun.upregul.tax.species$count/rep.fun.upregul.tax.species$total

writexl::write_xlsx(rep.bac.upregul.tax.species,"./GeneExpression/Expression/Taxonomic composition of upregulated genes in each condition_bacteria.xlsx")
writexl::write_xlsx(rep.fun.upregul.tax.species,"./GeneExpression/Expression/Taxonomic composition of upregulated genes in each condition_fungi.xlsx")

########  
####### Taxa contribution to DEGs #######
names(KO.deg.active.bac)[1] <- "GeneID"
KO.deg.active.bac<-merge(KO.deg.active.bac, rep.bac.tax.active, by ="GeneID")

names(KO.deg.pre.bac)[1] <- "GeneID"
KO.deg.pre.bac<-merge(KO.deg.pre.bac, rep.bac.tax.pre, by ="GeneID")

names(KO.deg.post.bac)[1] <- "GeneID"
KO.deg.post.bac<-merge(KO.deg.post.bac, rep.bac.tax.post, by ="GeneID")

KO.deg.bac <- rbind(KO.deg.active.bac,KO.deg.pre.bac,KO.deg.post.bac)
KO.deg.bac.top20 <- subset(KO.deg.bac, description %in% examine.pathway.list)

names(target.path.gene.bac.active)
names(KO.deg.active.fun)[1] <- "GeneID"
KO.deg.active.fun<-merge(KO.deg.active.fun, rep.fun.tax.active, by ="GeneID")

names(KO.deg.pre.fun)[1] <- "GeneID"
KO.deg.pre.fun<-merge(KO.deg.pre.fun, rep.fun.tax.pre, by ="GeneID")

names(KO.deg.post.fun)[1] <- "GeneID"
KO.deg.post.fun<-merge(KO.deg.post.fun, rep.fun.tax.post, by ="GeneID")


kegg.bac.deg.top20<-unique(kopath.degs.plot.sub.top20.bac$description)
kegg.fun.deg.top20<-unique(kopath.degs.plot.sub.top20.fun$description)

KO.deg.bac <- rbind(KO.deg.active.bac,KO.deg.pre.bac,KO.deg.post.bac)
KO.deg.fun <- rbind(KO.deg.active.fun,KO.deg.pre.fun,KO.deg.post.fun)


##### All pathways
KO.deg.bac.species <- KO.deg.bac %>% group_by(Phylum, Class, Order, Family, 
                                              Genus, Species, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.bac.species,"./GeneExpression/Expression/DEG_KEGGpath_species level contribution.xlsx")

KO.deg.fun.species <- KO.deg.fun %>% group_by(Phylum, Class, Order, Family, 
                                              Genus, Species, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.fun.species,"./GeneExpression/Expression/DEG_KEGGpath_species level contribution_fungi.xlsx")

#### 
KO.deg.bac.genus <- KO.deg.bac %>% group_by(Phylum, Class, Order, Family, 
                                            Genus, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.bac.genus,"./GeneExpression/Expression/DEG_KEGGpath_genus level contribution_bacteria.xlsx")

KO.deg.fun.genus <- KO.deg.fun %>% group_by(Phylum, Class, Order, Family, 
                                            Genus, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.fun.genus,"./GeneExpression/Expression/DEG_KEGGpath_genus level contribution_fungi.xlsx")

#### 
KO.deg.bac.family <- KO.deg.bac %>% group_by(Phylum, Class, Order, Family, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.bac.family,"./GeneExpression/Expression/DEG_KEGGpath_family level contribution_bacteria.xlsx")

KO.deg.fun.family <- KO.deg.fun %>% group_by(Phylum, Class, Order, Family, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.fun.family,"./GeneExpression/Expression/DEG_KEGGpath_family level contribution_fungi.xlsx")

#### 
KO.deg.bac.class <- KO.deg.bac %>% group_by(Phylum, Class,TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.bac.class,"./GeneExpression/Expression/DEG_KEGGpath_class level contribution_bacteria.xlsx")

KO.deg.fun.class <- KO.deg.fun %>% group_by(Phylum, Class, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.fun.class,"./GeneExpression/Expression/DEG_KEGGpath_class level contribution_fungi.xlsx")

#### 
KO.deg.bac.phylum <- KO.deg.bac %>% group_by(Phylum, TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.bac.phylum,"./GeneExpression/Expression/DEG_KEGGpath_phylum level contribution_bacteria.xlsx")

KO.deg.fun.phylum <- KO.deg.fun %>% group_by(Phylum,  TmStatus, pathway, description) %>% reframe(count = sum(count))

writexl::write_xlsx(KO.deg.fun.phylum,"./GeneExpression/Expression/DEG_KEGGpath_phylum level contribution_fungi.xlsx")


########## Fig. 4b Taxa contribution to each pathway
KO.deg.bac.examine <- subset(KO.deg.bac, description %in% kegg.bac.deg.top20)
KO.deg.fun.examine <- subset(KO.deg.fun, description %in% kegg.fun.deg.top20)

KO.deg.bac.examine.phylum <- KO.deg.bac.examine %>% 
  group_by(Phylum, TmStatus, pathway, description) %>% reframe(count = sum(count))
KO.deg.bac.examine.phylum<-KO.deg.bac.examine.phylum %>% group_by(TmStatus, pathway, description) %>% mutate(total = sum(count))
KO.deg.bac.examine.phylum$ratio <- KO.deg.bac.examine.phylum$count/KO.deg.bac.examine.phylum$total
KO.deg.bac.examine.phylum$Phylum2 <- KO.deg.bac.examine.phylum$Phylum
KO.deg.bac.examine.phylum$Phylum2[-(which(KO.deg.bac.examine.phylum$Phylum%in%c("Pseudomonadota","Actinomycetota","Acidobacteriota","Planctomycetota","Chloroflexota","Verrucomicrobiota","Myxococcota","Bacteria_unknown")))] <- "Others"
KO.deg.bac.examine.phylum.plot <- KO.deg.bac.examine.phylum %>% group_by(Phylum2,TmStatus, pathway, description) %>% reframe(count=sum(count),ratio=sum(ratio))
KO.deg.bac.examine.phylum.plot <- KO.deg.bac.examine.phylum.plot %>% arrange(description)
writexl::write_xlsx(KO.deg.bac.examine.phylum.plot,"./GeneExpression/Expression/KO.DEG_taxonomic data_bacteria_for plot.xlsx")

KO.deg.bac.examine.phylum.plot$TmStatus <- factor(KO.deg.bac.examine.phylum.plot$TmStatus,levels=c("Pre","Active","Post"))
KO.deg.bac.examine.phylum.plot$description <- factor(KO.deg.bac.examine.phylum.plot$description, levels =rev(unique(KO.deg.bac.examine.phylum.plot$description)))
KO.deg.bac.examine.phylum.plot$Phylum2 <- factor(KO.deg.bac.examine.phylum.plot$Phylum2, levels =c("Pseudomonadota","Actinomycetota","Acidobacteriota","Planctomycetota","Chloroflexota","Verrucomicrobiota","Myxococcota","Bacteria_unknown","Others"))


ggplot(KO.deg.bac.examine.phylum.plot, aes(x="", y=ratio, fill=Phylum2))+
  geom_bar(stat='identity')+facet_grid(rows = vars(description), cols = vars(TmStatus))+
  theme_void()+scale_fill_manual(values = c("Pseudomonadota"="#99cc66",
                                            "Actinomycetota"="#ffcc33",
                                            "Acidobacteriota"="#99ccff",
                                            "Planctomycetota"="#ff9933",
                                            "Chloroflexota"="#00c19f",
                                            "Verrucomicrobiota"="#cc99ff",
                                            "Myxococcota"="#619cff",
                                            "Bacteria_unknown"="#474747",
                                            "Others"="#999999"))+
  coord_polar('y', start=0)





KO.deg.fun.examine.class <- KO.deg.fun.examine %>% 
  group_by(Class, TmStatus, pathway, description) %>% reframe(count = sum(count))
KO.deg.fun.examine.class<-KO.deg.fun.examine.class %>% group_by(TmStatus, pathway, description) %>% mutate(total = sum(count))
KO.deg.fun.examine.class$ratio <- KO.deg.fun.examine.class$count/KO.deg.fun.examine.class$total
KO.deg.fun.examine.class$Class2 <- KO.deg.fun.examine.class$Class
KO.deg.fun.examine.class$Class2[-(which(KO.deg.fun.examine.class$Class%in%c("Agaricomycetes",
                                                                          "Dothideomycetes",
                                                                          "Saccharomycetes",
                                                                          "Leotiomycetes",
                                                                          "Eurotiomycetes",
                                                                          "Chytridiomycetes",
                                                                          "Sordariomycetes",
                                                                          "Mucoromycetes",
                                                                          "Ascomycota_unknown",
                                                                          "Fungi_unknown")))] <- "Others"
KO.deg.fun.examine.class.plot <- KO.deg.fun.examine.class %>% group_by(Class2,TmStatus, pathway, description) %>% reframe(count=sum(count),ratio=sum(ratio))
KO.deg.fun.examine.class.plot <- KO.deg.fun.examine.class.plot %>% arrange(description)
writexl::write_xlsx(KO.deg.fun.examine.class.plot,"./GeneExpression/Expression/KO.DEG_taxonomic data_fungi_for plot.xlsx")

KO.deg.fun.examine.class.plot$TmStatus <- factor(KO.deg.fun.examine.class.plot$TmStatus,levels=c("Pre","Active","Post"))
KO.deg.fun.examine.class.plot$description <- factor(KO.deg.fun.examine.class.plot$description, levels =rev(unique(KO.deg.fun.examine.class.plot$description)))
KO.deg.fun.examine.class.plot$Class2 <- factor(KO.deg.fun.examine.class.plot$Class2, levels =c("Agaricomycetes",
                                                                                               "Dothideomycetes",
                                                                                               "Saccharomycetes",
                                                                                               "Leotiomycetes",
                                                                                               "Eurotiomycetes",
                                                                                               "Chytridiomycetes",
                                                                                               "Sordariomycetes",
                                                                                               "Mucoromycetes",
                                                                                               "Ascomycota_unknown",
                                                                                               "Fungi_unknown","Others"))


ggplot(KO.deg.fun.examine.class.plot, aes(x="", y=ratio, fill=Class2))+
  geom_bar(stat='identity')+facet_grid(rows = vars(description), cols = vars(TmStatus))+
  theme_void()+scale_fill_manual(values = c("Agaricomycetes"="#ff6a72",
                                            "Dothideomycetes"="#33ccff",
                                            "Saccharomycetes"="#ffa633",
                                            "Leotiomycetes"="#64b200",
                                            "Eurotiomycetes"="#006699",
                                            "Chytridiomycetes"="#00cc99",
                                            "Sordariomycetes"="#cccc99",
                                            "Mucoromycetes"="#00a6ff",
                                            "Ascomycota_unknown"="#b385ff",
                                            "Fungi_unknown"="#474747",
                                            "Others"="#999999"))+
  coord_polar('y', start=0)

head(eggnog.rep.fun.tax)
eggnog.rep.fun.tax.express.RPKM.refine.sub[grepl("4.1.1.74",eggnog.rep.fun.tax.express.RPKM.refine.sub$EC),]$Species
eggnog.rep.fun.tax.express.RPKM.refine.sub[grepl("4.1.1.74",eggnog.rep.fun.tax.express.RPKM.refine.sub$EC),]$Species

eggnog.rep.bac.tax.express.RPKM.refine.sub[grepl("4.1.1.74",eggnog.rep.bac.tax.express.RPKM.refine.sub$EC),]$Species
eggnog.rep.bac.tax.express.RPKM.refine.sub[grepl("4.1.1.74",eggnog.rep.bac.tax.express.RPKM.refine.sub$EC),]$Species


######## GO overrepresentation plot
imp20.GOs <- readxl::read_xlsx("./20GOsfor.xlsx")
head(imp20.GOs)
imp20.GOs <- imp20.GOs %>% group_by(GO) %>% arrange(Count)
imp20.GOs$GO <- factor(imp20.GOs$GO, levels = imp20.GOs$GO)
ggplot(imp20.GOs, aes(x=Count, y=GO, color = Domain))+
  geom_point(aes(size = -log(qvalue)),shape=16, alpha = 0.8)+ facet_wrap(~Domain, scales = "free", ncol=1)+
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
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+theme(aspect.ratio = 2)
dev.off()
