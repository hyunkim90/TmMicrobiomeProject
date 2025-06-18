####### Fig. 4a. Annotated and expressed gene count at the taxa level #######
diamond.rep.gene.bac.contig.re
nrow(diamond.rep.gene.bac.contig.re)
nrow(eggnog.rep.bac.tax.abund.read)
head(eggnog.rep.bac.tax.abund.read)

load("./GeneExpression/Abundance/eggnog.rep.bac.tax.abund.read.RData")
load("./GeneExpression/Abundance/eggnog.rep.bac.tax.abund.RPKM.RData")
load("./GeneExpression/Abundance/eggnog.rep.bac.tax.abund.TPM.RData")
load("./GeneExpression/Abundance/eggnog.rep.fun.tax.abund.read.RData")
load("./GeneExpression/Abundance/eggnog.rep.fun.tax.abund.RPKM.RData")
load("./GeneExpression/Abundance/eggnog.rep.fun.tax.abund.TPM.RData")


###### 1-0. Refine taxonomy information
eggnog.rep.bac.tax.abund.read.refine <-eggnog.rep.bac.tax.abund.read
eggnog.rep.bac.tax.abund.read.refine$Phylum[which(eggnog.rep.bac.tax.abund.read.refine$Phylum == "Unidentified")] <-
  paste0(eggnog.rep.bac.tax.abund.read.refine$Domain[which(eggnog.rep.bac.tax.abund.read.refine$Phylum == "Unidentified")],"_unknown")
eggnog.rep.bac.tax.abund.read.refine$Class[which(eggnog.rep.bac.tax.abund.read.refine$Class == "Unidentified")] <-
  paste0(eggnog.rep.bac.tax.abund.read.refine$Phylum[which(eggnog.rep.bac.tax.abund.read.refine$Class == "Unidentified")],"_unknown")
eggnog.rep.bac.tax.abund.read.refine$Class<-gsub("unknown_unknown","unknown",eggnog.rep.bac.tax.abund.read.refine$Class)
eggnog.rep.bac.tax.abund.read.refine$Order[which(eggnog.rep.bac.tax.abund.read.refine$Order == "Unidentified")] <-
  paste0(eggnog.rep.bac.tax.abund.read.refine$Class[which(eggnog.rep.bac.tax.abund.read.refine$Order == "Unidentified")],"_unknown")
unique(eggnog.rep.bac.tax.abund.read.refine$Order)
eggnog.rep.bac.tax.abund.read.refine$Order<-gsub("unknown_unknown","unknown",eggnog.rep.bac.tax.abund.read.refine$Order)
eggnog.rep.bac.tax.abund.read.refine$Family[which(eggnog.rep.bac.tax.abund.read.refine$Family == "Unidentified")] <-
  paste0(eggnog.rep.bac.tax.abund.read.refine$Order[which(eggnog.rep.bac.tax.abund.read.refine$Family == "Unidentified")],"_unknown")
unique(eggnog.rep.bac.tax.abund.read.refine$Family)
eggnog.rep.bac.tax.abund.read.refine$Family<-gsub("unknown_unknown","unknown",eggnog.rep.bac.tax.abund.read.refine$Family)
eggnog.rep.bac.tax.abund.read.refine$Genus[which(eggnog.rep.bac.tax.abund.read.refine$Genus == "Unidentified")] <-
  paste0(eggnog.rep.bac.tax.abund.read.refine$Family[which(eggnog.rep.bac.tax.abund.read.refine$Genus == "Unidentified")],"_unknown")
unique(eggnog.rep.bac.tax.abund.read.refine$Genus)
eggnog.rep.bac.tax.abund.read.refine$Genus<-gsub("unknown_unknown","unknown",eggnog.rep.bac.tax.abund.read.refine$Genus)
eggnog.rep.bac.tax.abund.read.refine$Species[which(eggnog.rep.bac.tax.abund.read.refine$Species == "Unidentified")] <-
  paste0(eggnog.rep.bac.tax.abund.read.refine$Genus[which(eggnog.rep.bac.tax.abund.read.refine$Species == "Unidentified")],"_unknown")
unique(eggnog.rep.bac.tax.abund.read.refine$Species)
eggnog.rep.bac.tax.abund.read.refine$Species<-gsub("unknown_unknown","unknown",eggnog.rep.bac.tax.abund.read.refine$Species)


eggnog.rep.fun.tax.abund.read.refine <- eggnog.rep.fun.tax.abund.read
eggnog.rep.fun.tax.abund.read.refine$Phylum[which(eggnog.rep.fun.tax.abund.read.refine$Phylum == "Unidentified")] <-
  paste0(eggnog.rep.fun.tax.abund.read.refine$Kingdom[which(eggnog.rep.fun.tax.abund.read.refine$Phylum == "Unidentified")],"_unknown")
eggnog.rep.fun.tax.abund.read.refine$Class[which(eggnog.rep.fun.tax.abund.read.refine$Class == "Unidentified")] <-
  paste0(eggnog.rep.fun.tax.abund.read.refine$Phylum[which(eggnog.rep.fun.tax.abund.read.refine$Class == "Unidentified")],"_unknown")
eggnog.rep.fun.tax.abund.read.refine$Class<-gsub("unknown_unknown","unknown",eggnog.rep.fun.tax.abund.read.refine$Class)
eggnog.rep.fun.tax.abund.read.refine$Order[which(eggnog.rep.fun.tax.abund.read.refine$Order == "Unidentified")] <-
  paste0(eggnog.rep.fun.tax.abund.read.refine$Class[which(eggnog.rep.fun.tax.abund.read.refine$Order == "Unidentified")],"_unknown")
unique(eggnog.rep.fun.tax.abund.read.refine$Order)
eggnog.rep.fun.tax.abund.read.refine$Order<-gsub("unknown_unknown","unknown",eggnog.rep.fun.tax.abund.read.refine$Order)
eggnog.rep.fun.tax.abund.read.refine$Family[which(eggnog.rep.fun.tax.abund.read.refine$Family == "Unidentified")] <-
  paste0(eggnog.rep.fun.tax.abund.read.refine$Order[which(eggnog.rep.fun.tax.abund.read.refine$Family == "Unidentified")],"_unknown")
unique(eggnog.rep.fun.tax.abund.read.refine$Family)
eggnog.rep.fun.tax.abund.read.refine$Family<-gsub("unknown_unknown","unknown",eggnog.rep.fun.tax.abund.read.refine$Family)
eggnog.rep.fun.tax.abund.read.refine$Genus[which(eggnog.rep.fun.tax.abund.read.refine$Genus == "Unidentified")] <-
  paste0(eggnog.rep.fun.tax.abund.read.refine$Family[which(eggnog.rep.fun.tax.abund.read.refine$Genus == "Unidentified")],"_unknown")
unique(eggnog.rep.fun.tax.abund.read.refine$Genus)
eggnog.rep.fun.tax.abund.read.refine$Genus<-gsub("unknown_unknown","unknown",eggnog.rep.fun.tax.abund.read.refine$Genus)
eggnog.rep.fun.tax.abund.read.refine$Species[which(eggnog.rep.fun.tax.abund.read.refine$Species == "Unidentified")] <-
  paste0(eggnog.rep.fun.tax.abund.read.refine$Genus[which(eggnog.rep.fun.tax.abund.read.refine$Species == "Unidentified")],"_unknown")
unique(eggnog.rep.fun.tax.abund.read.refine$Species)
eggnog.rep.fun.tax.abund.read.refine$Species<-gsub("unknown_unknown","unknown",eggnog.rep.fun.tax.abund.read.refine$Species)


names(eggnog.rep.fun.tax.abund.read.refine)
eggnog.rep.bac.tax.abund.read.refine$Count <- 1
rep.gene.bac.annotated.phylum <- eggnog.rep.bac.tax.abund.read.refine %>% group_by(Phylum) %>% reframe(Count=sum(Count)) %>% arrange(desc(Count))
rep.gene.bac.annotated.phylum$Group <- "Annotated"
print(rep.gene.bac.annotated.phylum, n=300)
eggnog.rep.fun.tax.abund.read.refine$Count <- 1
rep.gene.fun.annotated.phylum <- eggnog.rep.fun.tax.abund.read.refine %>% group_by(Phylum) %>% reframe(Count=sum(Count)) %>% arrange(desc(Count))
rep.gene.fun.annotated.phylum$Group <- "Annotated"

##### Count expressed genes
rep.gene.bac.express<-eggnog.rep.bac.tax.express.read[rowSums(eggnog.rep.bac.tax.express.read[c(34:39)])>0,]
rep.gene.fun.express<-eggnog.rep.fun.tax.express.read[rowSums(eggnog.rep.fun.tax.express.read[c(34:39)])>0,]

rep.gene.bac.express$Phylum[which(rep.gene.bac.express$Phylum == "Unidentified")] <-
  paste0(rep.gene.bac.express$Domain[which(rep.gene.bac.express$Phylum == "Unidentified")],"_unknown")
rep.gene.bac.express$Class[which(rep.gene.bac.express$Class == "Unidentified")] <-
  paste0(rep.gene.bac.express$Phylum[which(rep.gene.bac.express$Class == "Unidentified")],"_unknown")
rep.gene.bac.express$Class<-gsub("unknown_unknown","unknown",rep.gene.bac.express$Class)
rep.gene.bac.express$Order[which(rep.gene.bac.express$Order == "Unidentified")] <-
  paste0(rep.gene.bac.express$Class[which(rep.gene.bac.express$Order == "Unidentified")],"_unknown")
unique(rep.gene.bac.express$Order)
rep.gene.bac.express$Order<-gsub("unknown_unknown","unknown",rep.gene.bac.express$Order)
rep.gene.bac.express$Family[which(rep.gene.bac.express$Family == "Unidentified")] <-
  paste0(rep.gene.bac.express$Order[which(rep.gene.bac.express$Family == "Unidentified")],"_unknown")
unique(rep.gene.bac.express$Family)
rep.gene.bac.express$Family<-gsub("unknown_unknown","unknown",rep.gene.bac.express$Family)
rep.gene.bac.express$Genus[which(rep.gene.bac.express$Genus == "Unidentified")] <-
  paste0(rep.gene.bac.express$Family[which(rep.gene.bac.express$Genus == "Unidentified")],"_unknown")
unique(rep.gene.bac.express$Genus)
rep.gene.bac.express$Genus<-gsub("unknown_unknown","unknown",rep.gene.bac.express$Genus)
rep.gene.bac.express$Species[which(rep.gene.bac.express$Species == "Unidentified")] <-
  paste0(rep.gene.bac.express$Genus[which(rep.gene.bac.express$Species == "Unidentified")],"_unknown")
unique(rep.gene.bac.express$Species)
rep.gene.bac.express$Species<-gsub("unknown_unknown","unknown",rep.gene.bac.express$Species)


rep.gene.fun.express $Phylum[which(rep.gene.fun.express$Phylum == "Unidentified")] <-
  paste0(rep.gene.fun.express$Kingdom[which(rep.gene.fun.express$Phylum == "Unidentified")],"_unknown")
rep.gene.fun.express$Class[which(rep.gene.fun.express$Class == "Unidentified")] <-
  paste0(rep.gene.fun.express$Phylum[which(rep.gene.fun.express$Class == "Unidentified")],"_unknown")
rep.gene.fun.express$Class<-gsub("unknown_unknown","unknown",rep.gene.fun.express$Class)
rep.gene.fun.express$Order[which(rep.gene.fun.express$Order == "Unidentified")] <-
  paste0(rep.gene.fun.express$Class[which(rep.gene.fun.express$Order == "Unidentified")],"_unknown")
unique(rep.gene.fun.express$Order)
rep.gene.fun.express$Order<-gsub("unknown_unknown","unknown",rep.gene.fun.express$Order)
rep.gene.fun.express$Family[which(rep.gene.fun.express$Family == "Unidentified")] <-
  paste0(rep.gene.fun.express$Order[which(rep.gene.fun.express$Family == "Unidentified")],"_unknown")
unique(rep.gene.fun.express$Family)
rep.gene.fun.express$Family<-gsub("unknown_unknown","unknown",rep.gene.fun.express$Family)
rep.gene.fun.express$Genus[which(rep.gene.fun.express$Genus == "Unidentified")] <-
  paste0(rep.gene.fun.express$Family[which(rep.gene.fun.express$Genus == "Unidentified")],"_unknown")
unique(rep.gene.fun.express$Genus)
rep.gene.fun.express$Genus<-gsub("unknown_unknown","unknown",rep.gene.fun.express$Genus)
rep.gene.fun.express$Species[which(rep.gene.fun.express$Species == "Unidentified")] <-
  paste0(rep.gene.fun.express$Genus[which(rep.gene.fun.express$Species == "Unidentified")],"_unknown")
unique(rep.gene.fun.express$Species)
rep.gene.fun.express$Species<-gsub("unknown_unknown","unknown",rep.gene.fun.express$Species)

rep.gene.bac.express$Count <- 1
rep.gene.bac.express.phylum <- rep.gene.bac.express %>% group_by(Phylum) %>% reframe(Count=sum(Count)) %>% arrange(desc(Count))
rep.gene.bac.express.phylum$Group <- "Expressed"
print(rep.gene.bac.express.phylum, n=300)
rep.gene.fun.express$Count <- 1
rep.gene.fun.express.phylum <- rep.gene.fun.express %>% group_by(Phylum) %>% reframe(Count=sum(Count)) %>% arrange(desc(Count))
rep.gene.fun.express.phylum$Group <- "Expressed"



##### Plotting
rep.gene.bac.count<-rbind(rep.gene.bac.annotated.phylum,rep.gene.bac.express.phylum)
rep.gene.bac.count.summary <- rep.gene.bac.count%>% group_by(Phylum) %>% reframe(Count=sum(Count)) %>% arrange(desc(Count))
dominant.bac.phyla <- head(rep.gene.bac.count.summary$Phylum,10)

rep.gene.bac.count$Phylum2 <- rep.gene.bac.count$Phylum
rep.gene.bac.count$Phylum2[-(which(rep.gene.bac.count$Phylum%in%dominant.bac.phyla))] <- "Others" 
rep.gene.bac.count.plot <- rep.gene.bac.count %>% group_by(Phylum2,Group) %>% reframe(Count=sum(Count)) %>% arrange(desc(Count))

rep.gene.bac.count.plot$Group <- factor(rep.gene.bac.count.plot$Group, levels=rev(c("Annotated","Expressed")))
bac.phyla.order<-c(dominant.bac.phyla[c(1,3,4,5,6,7,8,9,10,2)],"Others")
rep.gene.bac.count.plot$Phylum2 <- factor(rep.gene.bac.count.plot$Phylum2, levels=bac.phyla.order)

ggplot(rep.gene.bac.count.plot, aes(x=Phylum2, y=Group, color = Group))+geom_point(aes(size = Count))+
  #scale_color_manual(value = c())+
  theme(aspect.ratio=0.5)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position="right")+geom_text(data=rep.gene.bac.count.plot, aes(label = Count))

#### Stacked area plot
rep.gene.bac.count.plot$Group <-factor(rep.gene.bac.count.plot$Group, levels =c("Annotated","Expressed"))
rep.gene.bac.count.plot$Group2 <- as.numeric(as.factor(rep.gene.bac.count.plot$Group))

ggplot(rep.gene.bac.count.plot, aes(x=Group2, y=Count, fill=Phylum2)) + 
  geom_area()+scale_fill_manual(values=c("Pseudomonadota"="#99cc66",
                                         "Actinomycetota"="#ffcc33",
                                         "Acidobacteriota"="#99ccff",
                                         "Planctomycetota"="#ff9933",
                                         "Chloroflexota"="#00c19f",
                                         "Verrucomicrobiota"="#cc99ff",
                                         "Myxococcota"="#619cff",
                                         "Bacteroidota"="tomato3",
                                         "Gemmatimonadota"="lightpink2",
                                         "Bacteria_unknown"="#474747",
                                         "Others"="#999999"))+theme_classic()



rep.gene.fun.count<-rbind(rep.gene.fun.annotated.phylum,rep.gene.fun.express.phylum)
rep.gene.fun.count.summary <- rep.gene.fun.count%>% group_by(Phylum) %>% reframe(Count=sum(Count)) %>% arrange(desc(Count))

rep.gene.fun.count$Group <- factor(rep.gene.fun.count$Group, levels=rev(c("Annotated","Expressed")))
fun.phyla.order<-c(rep.gene.fun.count.summary$Phylum[c(1,2,3,5,6,7,8,9,10,4)])
rep.gene.fun.count$Phylum <- factor(rep.gene.fun.count$Phylum, levels=fun.phyla.order)

ggplot(rep.gene.fun.count, aes(x=Phylum, y=Group, color = Group))+geom_point(aes(size = Count))+
  #scale_color_manual(value = c())+
  theme(aspect.ratio=0.5)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=10, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  theme(legend.position="right")+geom_text(data=rep.gene.fun.count, aes(label = Count))


#### Stacked area plot
rep.gene.fun.count$Group <-factor(rep.gene.fun.count$Group, levels =c("Annotated","Expressed"))
rep.gene.fun.count$Group2 <- as.numeric(as.factor(rep.gene.fun.count$Group))
unique(rep.gene.fun.count$Phylum)
ggplot(rep.gene.fun.count, aes(x=Group2, y=Count, fill=Phylum)) + 
  geom_area()+scale_fill_manual(values=c("Basidiomycota"="tan",
                                         "Ascomycota"="steelblue",
                                         "Mucoromycota"="#ffa633",
                                         "Chytridiomycota"="indianred2",
                                         "Zoopagomycota"="darkolivegreen3",
                                         "Olpidiomycota"="plum",
                                         "Blastocladiomycota"="bisque2",
                                         "Microsporidia"="honeydew3",
                                         "Cryptomycota"="navajowhite4",
                                         "Fungi_unknown"="#474747",
                                         "Others"="#999999"))+theme_classic()

writexl::write_xlsx(rep.gene.fun.count,"../../../Manuscript/Comprehensive/Supplementary Files/Data/Finalized/Source_Fig3a_fungi.xlsx")
writexl::write_xlsx(rep.gene.bac.count,"../../../Manuscript/Comprehensive/Supplementary Files/Data/Finalized/Source_Fig3a_bacteria.xlsx")




