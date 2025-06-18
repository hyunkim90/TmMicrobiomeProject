library(phyloseq)
library(dplyr)
library(ggplot2)
######### Find AMGs involved in biogeochemical processes
gene.function.tab <- readxl::read_xlsx("./consensus_annotation_finalized.xlsx")
head(gene.function.tab)
cand.AMG.tab<-gene.function.tab[grepl("AMG_",gene.function.tab$FunctionalCategory),]

##### Find what I intended to search
nrow(cand.AMG.tab) #1,680
writexl::write_xlsx(cand.AMG.tab,"./putative pahge AMGs.xlsx")

###### Abundance and expression data
##### phage abundance (metaG-based and metaT-based)
load(file="./vMAG.abund.act.read.RData")
load(file="./vMAG.abund.RPKM.RData")
load(file="./vMAG.abund.TPM.RData")

load(file="./vMAG.expr.act.read.RData")
load(file="./vMAG.expr.RPKM.RData")
load(file="./vMAG.expr.TPM.RData")

###### Gene abundance and expression
load(file="./gene.abundance.tabs.true.RData")
load(file="./gene.expression.tabs.true.RData")

###### Host prediction results (php)
load(file="./host.prediction.res.RData")
php.res.tab.true.tax
host.info.table <- php.res.tab.true.tax
names(host.info.table)[3] <- "Genome"
host.info.table.best.hit <- host.info.table %>%group_by(Genome) %>% slice(which.max(score))
host.info.table.best.hit$Genome <- gsub("full+\\|+\\|","full\\|",host.info.table.best.hit$Genome)
host.info.table.best.hit$Genome <- gsub("0_partial+\\|+\\|","0_partial\\|",host.info.table.best.hit$Genome)
###### Life style prediction result
phatyp.res<-read.table("../../../../../KNU_Prof.HangilKim/ForestVirome/DNAvirus/reprePhages/reanalyzed/phatyp_prediction.tsv", sep ="\t", header =T)
names(phatyp.res)[1] <- "Genome"
phatyp.res <- phatyp.res[c(1,3)]
nrow(phatyp.res) #3,246
unique(phatyp.res$TYPE)
phage.list<-taxa_names(vMAG.abund.act.read)
phage.list2 <-  gsub("\\|",".",phage.list)
length(phage.list2[which(phage.list2%in%phatyp.res$Genome)])

phatyp.res$TYPE2 <- ifelse(phatyp.res$TYPE == "virulent","lytic",ifelse(phatyp.res$TYPE == "temperate","temperate","unknown"))
unique(phatyp.res$TYPE2)
phage.list2 <-  gsub("\\|",".",phage.list)
phatyp.res$Genome2 <- gsub("\\.","\\|",phatyp.res$Genome)
names(phatyp.res)
phatyp.res2 <- phatyp.res[c(4,3)]
names(phatyp.res2) <- c("Genome","type")

##### manual life style annotation

####### Manual classification ######
#### find lysogeny-related genes
lysogeny_keywords <- c(  "integrase","Integrase", "excisionase", "Excisionase","repressor", "Repressor",
                         "lysogeny", "Lysogeny","prophage", "Prophage","attachment site", "Attachment site","antirepressor","Antirepressor",
                         "recombinase","Recombinase","ParA","ParB","cro/C1","recombination","Recombination","excision","Excision")
lysogeny_keywords <- paste0(lysogeny_keywords,collapse = "|")
gene.function.tab.final3 <- gene.function.tab.final2
putative.temperate.phages<-unique(gene.function.tab.final3[grepl(lysogeny_keywords,gene.function.tab.final3$Description2),]$genome)
gene.function.tab.final3$Type[which(gene.function.tab.final3$genome %in% putative.temperate.phages)] <- "temperate"
length(unique(gene.function.tab.final3$genome[which(gene.function.tab.final3$Type == "temperate")])) #1670
length(unique(gene.function.tab.final3$genome[which(gene.function.tab.final3$Type == "lytic")])) #1576

lytic.phage.list <- unique(gene.function.tab.final3$genome[which(gene.function.tab.final3$Type == "lytic")])
temp.phage.list <- unique(gene.function.tab.final3$genome[which(gene.function.tab.final3$Type == "temperate")])


names(gene.function.tab.final3)
phatyp.res3 <-unique(gene.function.tab.final3[c(1,16)])
nrow(phatyp.res3)
names(phatyp.res3)[1] <- "Genome"
######### 
head(gene.abundance.tabs.true)
nrow(gene.abundance.tabs.true)

head(gene.expression.tabs.true)
nrow(gene.expression.tabs.true)

######### prevalence table ######
vMAG.abund.TPM
df.phage.metaG.abund <- data.frame(otu_table(vMAG.abund.TPM))
df.phage.metaG.prev <- df.phage.metaG.abund
df.phage.metaG.prev[df.phage.metaG.prev>0] <- 1
colSums(df.phage.metaG.prev)

df.phage.metaT.abund <- data.frame(otu_table(vMAG.expr.TPM))
df.phage.metaT.prev <- df.phage.metaT.abund
df.phage.metaT.prev[df.phage.metaT.prev>0] <- 1
colSums(df.phage.metaT.prev)


####### 1. proportion of lytic and temperate phages in the metagenomic samples
phage.metaG.rel <- microbiome::transform(vMAG.abund.act.read, transform = "compositional")
phage.metaG.rel.melt <- psmelt(phage.metaG.rel)
names(phage.metaG.rel.melt)[1] <-"Genome"
phage.metaG.rel.melt <- merge(phage.metaG.rel.melt,phatyp.res3, by = "Genome")
nrow(phage.metaG.rel.melt) #64,920
unique(phage.metaG.rel.melt$Genome[which(phage.metaG.rel.melt$Type == "lytic")])
unique(phage.metaG.rel.melt$Genome[which(phage.metaG.rel.melt$Type == "temperate")])
phage.metaG.rel.melt.type <- phage.metaG.rel.melt %>% group_by(Type,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
phage.metaG.rel.melt.type$Type <- factor(phage.metaG.rel.melt.type$Type, levels=rev(c("lytic","temperate")))
order.metaG.sample <-c("Y22TmM5","Y22TmM4","Y22TmM3","Y22TmM2","Y22TmM1","Y22TmD1","Y22TmD2","Y22TmD3","Y22TmD4","Y22TmD5",
                       "Y23TmM5","Y23TmM4","Y23TmM3","Y23TmM2","Y23TmM1","Y23TmD1","Y23TmD2","Y23TmD3","Y23TmD4","Y23TmD5")
phage.metaG.rel.melt.type$Sample <- factor(phage.metaG.rel.melt.type$Sample, levels=order.metaG.sample)

##### Bar plot
ggplot(phage.metaG.rel.melt.type, aes(x=Sample, y = Ratio, fill = Type)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  #scale_fill_manual(values = host.genus_colors) +
  
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

##### box plot (grouped by the Tm colonization status)
phage.metaG.rel.melt.type$TmColonization <- factor(phage.metaG.rel.melt.type$TmColonization, levels = c("Pre","Active","Post"))
ggplot(data=phage.metaG.rel.melt.type, aes(x=TmColonization, y=Ratio, fill=Type)) + geom_boxplot(width = 0.8) +
  stat_summary(fun=mean, colour="darkred", geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = .75, linetype = "dashed") + facet_wrap(~Type)+
  geom_point(position='jitter',shape=1, alpha=.3, size = 2)+
  #xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Ratio \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) +
  #scale_fill_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

##### Statistical analysis
##### Normality test
shapiro.test(phage.metaG.rel.melt.type[phage.metaG.rel.melt.type$Type=="lytic",]$Ratio) # not follow the normal distribution since p-value = 0.02371
shapiro.test(phage.metaG.rel.melt.type[phage.metaG.rel.melt.type$Type=="temperate",]$Ratio) # not follow the normal distribution since p-value = 0.007196

#### On the basis of the normality test, appropriate statistical analyses should be applied.
#### Kruskal-wallis test
kw.lytic.ratio<-kruskal.test(Ratio ~ TmColonization, data = phage.metaG.rel.melt.type[phage.metaG.rel.melt.type$Type=="lytic",])
kw.lytic.ratio$p.value #0.001524466

#library(FSA)
DT.lytic.ratio= FSA::dunnTest(Ratio ~ TmColonization, data = phage.metaG.rel.melt.type[phage.metaG.rel.melt.type$Type=="lytic",],
                               method="bh")
PT.lytic.ratio = DT.lytic.ratio$res
writexl::write_xlsx(PT.fun.asv.rich,"Richness_KW_Dunn_statistics_ASV_level.xlsx")
#library(rcompanion)

dunn.lytic.ratio<-rcompanion::cldList(P.adj ~ Comparison,
                          data = PT.lytic.ratio,
                          threshold = 0.05)

kw.temperate.ratio<-kruskal.test(Ratio ~ TmColonization, data = phage.metaG.rel.melt.type[phage.metaG.rel.melt.type$Type=="temperate",])
kw.temperate.ratio$p.value #0.8852026

#library(FSA)

DT.temperate.ratio= FSA::dunnTest(Ratio ~ TmColonization, data = phage.metaG.rel.melt.type[phage.metaG.rel.melt.type$Type=="temperate",],
                              method="bh")
PT.temperate.ratio = DT.temperate.ratio$res
writexl::write_xlsx(PT.fun.asv.rich,"Richness_KW_Dunn_statistics_ASV_level.xlsx")
#library(rcompanion)
dunn.temperate.ratio<-rcompanion::cldList(P.adj ~ Comparison,
                                      data = PT.temperate.ratio,
                                      threshold = 0.05)

###### 2. phage alpha diversity
set.seed(1234)
######## 2-1. all phages
vMAG.abund.act.read.rarefy <- phyloseq::rarefy_even_depth(vMAG.abund.act.read, min(colSums(phyloseq::otu_table(vMAG.abund.act.read))))
colSums(phyloseq::otu_table(vMAG.abund.act.read.rarefy ))### all phages were used
all.alpha.index<-microbiome::alpha(vMAG.abund.act.read.rarefy, index = "all")
head(all.alpha.index)
phage.metadata <- data.frame(sample_data(vMAG.abund.act.read.rarefy ))
phage.meta.diversity <- merge(phage.metadata, all.alpha.index, by = "row.names")
rownames(phage.meta.diversity) <- phage.meta.diversity$Row.names
phage.meta.diversity <- phage.meta.diversity[-c(1)]

######## 2-2. lytic phages
lytic.phage.metaG <- phyloseq::prune_taxa(lytic.phage.list, vMAG.abund.act.read)
temp.phage.metaG <- phyloseq::prune_taxa(temp.phage.list, vMAG.abund.act.read)

lytic.phage.metaG.rarefy <- phyloseq::rarefy_even_depth(lytic.phage.metaG, min(colSums(phyloseq::otu_table(lytic.phage.metaG))))
colSums(phyloseq::otu_table(lytic.phage.metaG.rarefy ))### all phages were used
lytic.alpha.index<-microbiome::alpha(lytic.phage.metaG.rarefy, index = "all")
head(lytic.alpha.index)


lytic.phage.meta.diversity <- merge(phage.metadata, lytic.alpha.index, by = "row.names")
rownames(lytic.phage.meta.diversity) <- lytic.phage.meta.diversity$Row.names
lytic.phage.meta.diversity <- lytic.phage.meta.diversity[-c(1)]
lytic.phage.meta.diversity$Type <- "lytic"

temp.phage.metaG.rarefy <- phyloseq::rarefy_even_depth(temp.phage.metaG, min(colSums(phyloseq::otu_table(temp.phage.metaG))))
colSums(phyloseq::otu_table(temp.phage.metaG.rarefy ))### all phages were used
temp.alpha.index<-microbiome::alpha(temp.phage.metaG.rarefy, index = "all")
head(temp.alpha.index)

temp.phage.meta.diversity <- merge(phage.metadata, temp.alpha.index, by = "row.names")
rownames(temp.phage.meta.diversity) <- temp.phage.meta.diversity$Row.names
temp.phage.meta.diversity <- temp.phage.meta.diversity[-c(1)]
temp.phage.meta.diversity$Type <- "temperate"


type.alpha.meta <- rbind(lytic.phage.meta.diversity,temp.phage.meta.diversity)
names(type.alpha.meta)


type.alpha.meta$TmColonization <- factor(type.alpha.meta$TmColonization,levels=c("Pre","Active","Post"))

ggplot(type.alpha.meta, aes(x=TmColonization, y=chao1, fill=Type)) +
  geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.8), dotsize = 0.8)+facet_wrap(~Type)+
  #scale_fill_manual(values = color_treatment)+
  # stat_summary(fun.data=mean_sdl, mult=1, 
  #              geom="pointrange", color="red")+
  theme(legend.position = "top")+
  # geom_errorbar(aes(ymin=meanChao1-sdChao1, ymax=meanChao1+sdChao1), width=.2,
  #               position=position_dodge(.8))+
  #geom_point(aes(y=chao1), size = 2.5,position=position_dodge(.8))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


ggplot(type.alpha.meta, aes(x=TmColonization, y=diversity_shannon, fill=Type)) +
  geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.8), dotsize = 0.8)+facet_wrap(~Type)+
  #scale_fill_manual(values = color_treatment)+
  # stat_summary(fun.data=mean_sdl, mult=1, 
  #              geom="pointrange", color="red")+
  theme(legend.position = "top")+
  # geom_errorbar(aes(ymin=meanChao1-sdChao1, ymax=meanChao1+sdChao1), width=.2,
  #               position=position_dodge(.8))+
  #geom_point(aes(y=chao1), size = 2.5,position=position_dodge(.8))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())




########## 3. Host prediction results ##########
########## 
load(file="./host.prediction.res.RData")
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


host.info.table.best.hit2

names(host.info.table.best.hit2)
host.info.table.best.hit3 <- host.info.table.best.hit2[c(1,9:15)]
colnames(host.info.table.best.hit3) <- c("Genome","Host.domain","Host.phylum","Host.class","Host.order",
                                         "Host.family","Host.genus","Host.species")
######## bacterial host of lytic and temperate phages 
names(phage.metaG.rel.melt)
names(host.info.table.best.hit3)
host.info.table.best.hit4 <- host.info.table.best.hit3
names(host.info.table.best.hit4)[1]<- "Genome"
phage.metaG.rel.melt <- merge(phage.metaG.rel.melt,host.info.table.best.hit4, by = "Genome")
nrow(phage.metaG.rel.melt) #64,920

phage.metaG.rel.melt.type.host <- phage.metaG.rel.melt %>% group_by(Type,Host.phylum,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
phage.metaG.rel.melt.type.host$Type <- factor(phage.metaG.rel.melt.type.host$Type, levels=c("lytic","temperate"))
phage.metaG.rel.melt.type.host$Sample <- factor(phage.metaG.rel.melt.type.host$Sample, levels=order.metaG.sample)

##### Bar plot
ggplot(phage.metaG.rel.melt.type.host, aes(x=Sample, y = Ratio, fill = Host.phylum)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+facet_wrap(~Type)+
  #scale_fill_discrete() +
  #scale_fill_manual(values = host.genus_colors) +
  
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

  
  ######## simplified the phyla
  phage.metaG.rel.melt2 <- phage.metaG.rel.melt
  phage.metaG.rel.melt2$Host.phylum2 <- phage.metaG.rel.melt2$Host.phylum
  phage.metaG.rel.melt2$Host.phylum2[which(phage.metaG.rel.melt2$Host.class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
  phage.metaG.rel.melt2$Host.phylum2[which(phage.metaG.rel.melt2$Host.class == "Betaproteobacteria")] <- "Betaproteobacteria"
  phage.metaG.rel.melt2$Host.phylum2[which(phage.metaG.rel.melt2$Host.class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
  
  ord.host.phyla <- phage.metaG.rel.melt2 %>% group_by(Type,Host.phylum2,Sample) %>% reframe(Ratio = sum(Abundance))
  ord.host.phyla <- ord.host.phyla %>% group_by(Host.phylum2) %>% reframe(Ratio = mean(Ratio)) %>% arrange(desc(Ratio))
  phage.metaG.rel.melt.type.host$Host.phylum2 <- phage.metaG.rel.melt.type.host$Host.phylum
  print(  ord.host.phyla, n=80)
  dorm.host.phyla<-ord.host.phyla$Host.phylum2[which(ord.host.phyla$Ratio > 0.003)]
  min.host.phyla<-ord.host.phyla$Host.phylum2[which(ord.host.phyla$Ratio <= 0.003)]
  phage.metaG.rel.melt2$Host.phylum2[which(phage.metaG.rel.melt2$Host.phylum2 %in% min.host.phyla)] <- "Others"
  
  
  ord.vec.host.phyla<-append(dorm.host.phyla[c(1:6,8:14,7)],"Others")
  
  phage.metaG.rel.melt.type.host <- phage.metaG.rel.melt2 %>% group_by(Type,Host.phylum2,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
  phage.metaG.rel.melt.type.host$Type <- factor(phage.metaG.rel.melt.type.host$Type, levels=c("lytic","temperate"))
  phage.metaG.rel.melt.type.host$Sample <- factor(phage.metaG.rel.melt.type.host$Sample, levels=order.metaG.sample)
  phage.metaG.rel.melt.type.host$Host.phylum2 <- factor(phage.metaG.rel.melt.type.host$Host.phylum2,levels =rev(ord.vec.host.phyla))
  unique(  phage.metaG.rel.melt.type.host$Host.phylum2)
  
  
  library(circlize)
  col_fun.function = colorRamp2(c(0,3), c("white", "#99cc66"))
  col_fun.function(seq(1,3))
  
  host.phyla_colors<- c("Pseudomonadota"="#99cc66",
    "Actinomycetota"="#ffcc33",
    "Acidobacteriota"="#99ccff",
    "Planctomycetota"="#ff9933",
    "Chloroflexota"="#00c19f",
    "Verrucomicrobiota"="#cc99ff",
    "Myxococcota"="#619cff",
    "Gemmatimonadota"="#FF7F50",
    "Cyanobacteria"="#00A5CF",
    "Bacteria_unknown"="#474747",
    "Bacteroidota"="indianred2",
    "Thermomicrobiota"="#8B5A2B",
    "Alphaproteobacteria"="#99cc66",
    "Betaproteobacteria"="#BDDD99",
    "Gammaproteobacteria"="#DEEECB",
    "Bacillota"="#CC7722",
    "Others"="#999999")
  
  ##### Bar plot
  ggplot(phage.metaG.rel.melt.type.host, aes(x=Sample, y = Ratio, fill = Host.phylum2)) + 
    geom_bar(stat="identity", width = 0.8, position = 'stack')+facet_wrap(~Type)+
  #scale_fill_discrete() +
  scale_fill_manual(values = host.phyla_colors) +
  
  xlab('')+
    ylab("Ratio \n") +
    #ggtitle("genus Community Composition by Sample \n") +
    ## adjust positions
    guides(fill = guide_legend(ncol = 1,reverse = F))+
    theme(legend.position="right") +
    theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    scale_y_continuous(breaks=seq(0,1,.2))+
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
  
  writexl::write_xlsx(phage.metaG.rel.melt.type.host,"./Source_host_phage_life_strategy_all_rep_250228_final.xlsx")
  
  ##### Agglomerated at the Tm colonization level
  phage.metaG.rel.melt.type.host.mean <- phage.metaG.rel.melt.type.host %>% group_by(TmColonization,Host.phylum2,Type) %>% reframe(Ratio = mean(Ratio))
  phage.metaG.rel.melt.type.host.mean$TmColonization <- factor(phage.metaG.rel.melt.type.host.mean$TmColonization,
                                                               levels=c("Pre","Active","Post"))
  ggplot(phage.metaG.rel.melt.type.host.mean, aes(x=TmColonization, y = Ratio, fill = Host.phylum2)) + 
    geom_bar(stat="identity", width = 0.8, position = 'stack')+facet_wrap(~Type)+
    #scale_fill_discrete() +
    scale_fill_manual(values = host.phyla_colors) +
    
    xlab('')+
    ylab("Ratio \n") +
    #ggtitle("genus Community Composition by Sample \n") +
    ## adjust positions
    guides(fill = guide_legend(ncol = 1,reverse = F))+
    theme(legend.position="right") +
    theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    scale_y_continuous(breaks=seq(0,1,.2))+
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
  
###### ratio at the other taxonomic levels
  #### Phylum - original
  phage.metaG.rel.melt.type.host.phyla <- phage.metaG.rel.melt %>% group_by(type,Host.phylum,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
  #### Class
  phage.metaG.rel.melt.type.host.class <- phage.metaG.rel.melt %>% group_by(type,Host.phylum,Host.class,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
  #### Order
  phage.metaG.rel.melt.type.host.order <- phage.metaG.rel.melt %>% group_by(type,Host.phylum,Host.class,Host.order,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
  #### Family
  phage.metaG.rel.melt.type.host.family <- phage.metaG.rel.melt %>% group_by(type,Host.phylum,Host.class,Host.order,Host.family,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
  #### Genus
  phage.metaG.rel.melt.type.host.genus <- phage.metaG.rel.melt %>% group_by(type,Host.phylum,Host.class,Host.order,Host.family,Host.genus,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
  #### Species
  phage.metaG.rel.melt.type.host.species <- phage.metaG.rel.melt %>% group_by(type,Host.phylum,Host.class,Host.order,Host.family,Host.genus,Host.species,Sample,TmColonization) %>% reframe(Ratio = sum(Abundance))
  
  ######## Consider genome completeness
  host.info.table.best.hit3.qualified <- host.info.table.best.hit3[host.info.table.best.hit3$genome %in% qualified.phages,]
  nrow(host.info.table.best.hit3.qualified)
  unique(host.info.table.best.hit3.qualified$Host.order)
  
  host.info.table.qualified.lytic <- host.info.table.best.hit3.qualified[host.info.table.best.hit3.qualified$genome %in% lytic.phage.list,]
  unique(host.info.table.qualified.lytic$Host.genus)
  unique(host.info.table.qualified.lytic$Host.family)
  
  host.info.table.qualified.lytic$count <- 1
  host.info.table.qualified.lytic.summary <- host.info.table.qualified.lytic %>% group_by(Host.phylum,Host.class,Host.order,Host.family,Host.genus) %>% reframe(count =sum(count)) %>% arrange(desc(count))
  
  host.info.table.qualified.temp <- host.info.table.best.hit3.qualified[host.info.table.best.hit3.qualified$genome %in% temp.phage.list,]
  unique(host.info.table.qualified.temp$Host.family)
  
  host.info.table.qualified.temp$count <- 1
  host.info.table.qualified.temp.summary <- host.info.table.qualified.temp %>% group_by(Host.phylum,Host.class,Host.order,Host.family,Host.genus) %>% reframe(count =sum(count)) %>% arrange(desc(count))
  
  
  ###### 4. AMGs in lytic and temperate phages ########
  ######## Use genomes showing over or equal to 30% completeness
  quality.tab<- read.table("~/Desktop/SNU/KNU_Prof.HangilKim/ForestVirome/DNAvirus/BasicStatistics/GenomeQuality/quality_summary.tsv",
                           sep = "\t", header=T)
  head(quality.tab)
  quality.tab.2<- subset(quality.tab, contig_id %in% taxa_names(vMAG.abund.act.read))
  qualified.phages<-quality.tab.2$contig_id[which(quality.tab.2$completeness >= 30)] #396 phages
  
  phatyp.res.qualified<-phatyp.res3[phatyp.res3$Genome %in% qualified.phages,]
  unique(phatyp.res.qualified$Type)
  
  gene.function.tab
  head(gene.function.tab)
  gene.function.tab$genome <- paste0(gene.function.tab$gene,".S")
  gene.function.tab$genome <- gsub("_+\\d+.S+","",gene.function.tab$genome)
  gene.function.tab.qualified <- subset(gene.function.tab, genome %in% qualified.phages)
  gene.function.tab.qualified.AMG <- gene.function.tab.qualified[grepl("AMG_",gene.function.tab.qualified$FunctionalCategory),] #613 genes
  
  unique(gene.function.tab.qualified.AMG$FunctionalCategory)
  gene.function.tab.qualified.AMG$gene
  
  ##### Add phrog - Kegg link
  ref.phrog.ko.tab<-read.table("./AMGannotation/phrog/Unique_KO_Accessions.csv", sep =",", header =T,fill=T,quote="")
  head(ref.phrog.ko.tab)
  ref.phrog.ko.tab$KO_accessions[which(ref.phrog.ko.tab$KO_accessions=="")] <- "-"
  names(ref.phrog.ko.tab) <- c("accession","KO")
  gene.function.tab2.phrog <- gene.function.tab2[gene.function.tab2$AnntotedDB == "PHROGs",]
  names(gene.function.tab2.phrog)
  gene.function.tab2.phrog2<-merge(gene.function.tab2.phrog,ref.phrog.ko.tab,by="accession")
  nrow(gene.function.tab2.phrog2)
  nrow(gene.function.tab2.phrog)
  
  gene.function.tab2.others <- gene.function.tab2[gene.function.tab2$AnntotedDB != "PHROGs",]
  gene.function.tab2.others$KO <- "-"
  names(gene.function.tab2.phrog2)
  names(gene.function.tab2.others)
  gene.function.tab2.phrog2 <- gene.function.tab2.phrog2[c(2,1,3:8)]
  gene.function.tab3<-rbind(gene.function.tab2.phrog2,gene.function.tab2.others)
  
  gene.function.tab3.AMG<-gene.function.tab3[grepl("AMG_", gene.function.tab3$FunctionalCategory),]
  writexl::write_xlsx( gene.function.tab3.AMG,"./gene.function.tab3.AMG.xlsx")
  
  
  ####### Host-specific AMGs #####
  ##### Assign host information to the gene annotation result
  ###### without KO information
  gene.function.tab2
  names(gene.function.tab2)
  host.info.table.best.hit3
  names(host.info.table.best.hit3)[1] <- "genome"
  
  
  gene.function.tab3<-gene.function.tab2[!(grepl("TmD1IndAssemContig0359700|TmDContig00072623",gene.function.tab2$genome)),]
  length(unique(gene.function.tab3$genome))
  gene.function.tab4<-merge(gene.function.tab3,host.info.table.best.hit3, by ="genome")
nrow(gene.function.tab4)

gene.function.tab4$count <- 1
gene.function.host.summary <- gene.function.tab4 %>% group_by(Host.genus,Description2,FunctionalCategory) %>%reframe(count = sum(count))
gene.function.host.summary <- gene.function.tab4 %>% group_by(Host.family,FunctionalCategory) %>%reframe(count = sum(count))

gene.function.host.summary <- gene.function.tab4 %>% group_by(Host.phylum,Description2,FunctionalCategory) %>%reframe(count = sum(count))
writexl::write_xlsx(gene.function.host.summary,"./gene.function.host.summary_family_@.xlsx")

################## Appendix #####################
####### 1. Active and dormant phages (metaG-based approach showed gradual increase in the relative abundance of lytic phages from pre to post)
###### proportion of lytic and temperate phages in the metatranscriptomic samples ####
####### Separate active and dormant phages first
gene.expression.tabs.true
gene.expression.tabs.true2 <- gene.expression.tabs.true
gene.expression.tabs.true2$Genome <- paste0(gene.expression.tabs.true2$Name,".S")
gene.expression.tabs.true2$Genome  <- gsub("_+\\d+.S","", gene.expression.tabs.true2$Genome )
head(gene.expression.tabs.true2)
phage.metaT.abund.tab <- gene.expression.tabs.true2 %>% group_by(Genome,Sample) %>% reframe(TPM = sum(TPM),RPKM=sum(RPKM),read = sum(NumReads))

active.phage <- unique(phage.metaT.abund.tab[phage.metaT.abund.tab$TPM > 0,]$Genome) #3,241 active phages (at least one metaT sample)


##### Metagenomic abundance of active and dormant phages
phage.metaG.rel.melt.sub <-phage.metaG.rel.melt[phage.metaG.rel.melt$Sample %in% order.metaT.sample,]

unique(phage.metaT.abund.tab$Sample)
active.phage.Y23TmM5 <- unique(phage.metaT.abund.tab[phage.metaT.abund.tab$TPM > 0&phage.metaT.abund.tab$Sample=="23R5NE",]$Genome)
active.phage.Y23TmM4 <- unique(phage.metaT.abund.tab[phage.metaT.abund.tab$TPM > 0&phage.metaT.abund.tab$Sample=="23R4NE",]$Genome)
active.phage.Y23TmM1 <- unique(phage.metaT.abund.tab[phage.metaT.abund.tab$TPM > 0&phage.metaT.abund.tab$Sample=="23R1NA",]$Genome)
active.phage.Y23TmD1 <- unique(phage.metaT.abund.tab[phage.metaT.abund.tab$TPM > 0&phage.metaT.abund.tab$Sample=="23R1PB",]$Genome)
active.phage.Y23TmD4 <- unique(phage.metaT.abund.tab[phage.metaT.abund.tab$TPM > 0&phage.metaT.abund.tab$Sample=="23R4PE",]$Genome)
active.phage.Y23TmD5 <- unique(phage.metaT.abund.tab[phage.metaT.abund.tab$TPM > 0&phage.metaT.abund.tab$Sample=="23R5PE",]$Genome)

phage.metaG.rel.melt.sub$activity <- "0"
phage.metaG.rel.melt.sub$activity[which(phage.metaG.rel.melt.sub$Sample=="Y23TmM5" & phage.metaG.rel.melt.sub$Genome%in%active.phage.Y23TmM5)] <- "active"
phage.metaG.rel.melt.sub$activity[which(phage.metaG.rel.melt.sub$Sample=="Y23TmM4" & phage.metaG.rel.melt.sub$Genome%in%active.phage.Y23TmM4)] <- "active"
phage.metaG.rel.melt.sub$activity[which(phage.metaG.rel.melt.sub$Sample=="Y23TmM1" & phage.metaG.rel.melt.sub$Genome%in%active.phage.Y23TmM1)] <- "active"
phage.metaG.rel.melt.sub$activity[which(phage.metaG.rel.melt.sub$Sample=="Y23TmD1" & phage.metaG.rel.melt.sub$Genome%in%active.phage.Y23TmD1)] <- "active"
phage.metaG.rel.melt.sub$activity[which(phage.metaG.rel.melt.sub$Sample=="Y23TmD4" & phage.metaG.rel.melt.sub$Genome%in%active.phage.Y23TmD4)] <- "active"
phage.metaG.rel.melt.sub$activity[which(phage.metaG.rel.melt.sub$Sample=="Y23TmD5" & phage.metaG.rel.melt.sub$Genome%in%active.phage.Y23TmD5)] <- "active"
phage.metaG.rel.melt.sub$activity[which(phage.metaG.rel.melt.sub$activity=="0")] <- "dormant"

phage.metaG.rel.melt.sub.active <- phage.metaG.rel.melt.sub[phage.metaG.rel.melt.sub$activity == "active",]
phage.metaG.rel.melt.sub.active
head(phage.metaG.rel.melt.sub.active)


phage.metaG.rel.melt.sub.active.type<- phage.metaG.rel.melt.sub.active %>% group_by(type,Sample,TmColonization,activity) %>% reframe(Ratio = sum(Abundance))
phage.metaG.rel.melt.sub.active.type$Type <- factor(phage.metaG.rel.melt.sub.active.type$Type, levels=rev(c("lytic","temperate","unknown")))
order.metaG.sample <-c("Y22TmM5","Y22TmM4","Y22TmM3","Y22TmM2","Y22TmM1","Y22TmD1","Y22TmD2","Y22TmD3","Y22TmD4","Y22TmD5",
                       "Y23TmM5","Y23TmM4","Y23TmM3","Y23TmM2","Y23TmM1","Y23TmD1","Y23TmD2","Y23TmD3","Y23TmD4","Y23TmD5")
phage.metaG.rel.melt.sub.active.type$Sample <- factor(phage.metaG.rel.melt.sub.active.type$Sample, levels=order.metaT.sample)

##### Bar plot
ggplot(phage.metaG.rel.melt.sub.active.type, aes(x=Sample, y = Ratio, fill = type)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  #scale_fill_manual(values = host.genus_colors) +
  
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









