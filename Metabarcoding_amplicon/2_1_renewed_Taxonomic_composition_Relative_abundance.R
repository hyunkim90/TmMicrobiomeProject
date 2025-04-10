###### Update bacterial taxonomy to sync the amplicon info to metagenome data ###
#### load refined taxonomy information 
renewed.tax.tab <- readxl::read_xlsx("~/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/renewed_tax_table.xlsx")
nrow(renewed.tax.tab)
head(renewed.tax.tab)

origin.tax.tab<-data.frame(phyloseq::tax_table(bac.clean.ss.dna.f))
origin.tax.tab
renewed.tax.tab.sub <- data.frame(subset(renewed.tax.tab, OTU %in% rownames(origin.tax.tab)))
nrow(renewed.tax.tab.sub)
names(renewed.tax.tab.sub) <- c("OTU","Phylum","Class","Order","Family","Genus")
rownames(renewed.tax.tab.sub) <- renewed.tax.tab.sub$OTU
renewed.tax.tab.sub <-renewed.tax.tab.sub[-c(1)]
bac.clean.ss.renew<-bac.clean.ss.dna.f
phyloseq::tax_table(bac.clean.ss.renew) <- phyloseq::tax_table(as.matrix(renewed.tax.tab.sub))


############### Taxonomic composition with renewed taxonomy information ###########
bac.renew.rel <- microbiome::transform(bac.clean.ss.renew, "compositional")
bac.renew.rel.tab<-phyloseq::psmelt(bac.renew.rel)
names(bac.renew.rel.tab)[7] <- "TmStatus"
bac.renew.rel.tab$TmStatus <- ifelse(bac.renew.rel.tab$TmStatus == "yes","TmD","TmM")
head(bac.renew.rel.tab)
bac.renew.rel.tab$Replicate <- paste0(bac.renew.rel.tab$TmStatus,bac.renew.rel.tab$compartment)

######## Agglomerated at the phylum level
######### 1. Each sample level
unique(bac.renew.rel.tab$Phylum)
library(dplyr)
bac.renew.rel.phylum.each.sample <- bac.renew.rel.tab %>% group_by(Replicate,month,Phylum) %>% reframe(Ratio = sum(Abundance))
head(bac.renew.rel.phylum.each.sample)

######### 2. TmStatus level
bac.renew.rel.phylum.TmStatus <- bac.renew.rel.tab %>% group_by(TmStatus,Replicate,month,Phylum) %>% reframe(Ratio = sum(Abundance))
head(bac.renew.rel.phylum.TmStatus)
bac.renew.rel.phylum.TmStatus <- bac.renew.rel.phylum.TmStatus %>% group_by(TmStatus,month,Phylum) %>% reframe(meanRatio = mean(Ratio), sdRatio = sd(Ratio))

####### Add color index into data frames
bac.renew.rel.abundant <- bac.renew.rel.phylum.TmStatus %>% group_by(Phylum) %>% reframe(Ratio = mean(meanRatio)) %>% arrange(desc(Ratio)) %>% filter(Ratio >= 0.005)

bac.renew.rel.phylum.each.sample$Phylum2 <- bac.renew.rel.phylum.each.sample$Phylum
bac.renew.rel.phylum.TmStatus$Phylum2 <- bac.renew.rel.phylum.TmStatus$Phylum

bac.renew.rel.phylum.each.sample$Phylum2[!(bac.renew.rel.phylum.each.sample$Phylum%in% bac.renew.rel.abundant$Phylum)] <- "Others"
bac.renew.rel.phylum.TmStatus$Phylum2[!(bac.renew.rel.phylum.TmStatus$Phylum%in% bac.renew.rel.abundant$Phylum)] <- "Others"

ord.bac.renew.phyla <- c(bac.renew.rel.abundant$Phylum,"Others")


order.rep <- c("TmM5","TmM4","TmM3","TmM2","TmM1","TmD1","TmD2","TmD3","TmD4","TmD5")
bac.renew.rel.phylum.each.sample$month <- factor(bac.renew.rel.phylum.each.sample$month, levels = order.month)
bac.renew.rel.phylum.each.sample$Replicate <- factor(bac.renew.rel.phylum.each.sample$Replicate, levels = order.rep)
bac.renew.rel.phylum.each.sample$Phylum2 <- factor(bac.renew.rel.phylum.each.sample$Phylum2, levels = rev(ord.bac.renew.phyla))

bac.renew.rel.phylum.TmStatus$month <- factor(bac.renew.rel.phylum.TmStatus$month, levels = order.month)
bac.renew.rel.phylum.TmStatus$TmStatus <- factor(bac.renew.rel.phylum.TmStatus$TmStatus, levels = c("TmM","TmD"))
bac.renew.rel.phylum.TmStatus$Phylum2 <- factor(bac.renew.rel.phylum.TmStatus$Phylum2, levels = rev(ord.bac.renew.phyla))

####### plotting
##### each replicate
library(ggplot2)
ggplot(bac.renew.rel.phylum.each.sample, aes(x=Replicate, y = Ratio, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~month, ncol=1)+
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


ggplot(bac.renew.rel.phylum.TmStatus, aes(x=month, y = meanRatio, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~TmStatus, ncol=1)+
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



######## Agglomerated at the genus level
######### 1. Each sample level
unique(bac.renew.rel.tab$Genus)
library(dplyr)
bac.renew.rel.Genus.each.sample <- bac.renew.rel.tab %>% group_by(Replicate,month,Phylum,Genus) %>% reframe(Ratio = sum(Abundance))
head(bac.renew.rel.Genus.each.sample)

######### 2. TmStatus level
bac.renew.rel.Genus.TmStatus <- bac.renew.rel.tab %>% group_by(TmStatus,Replicate,month,Phylum,Genus) %>% reframe(Ratio = sum(Abundance))
head(bac.renew.rel.Genus.TmStatus)
bac.renew.rel.Genus.TmStatus <- bac.renew.rel.Genus.TmStatus %>% group_by(TmStatus,month,Phylum,Genus) %>% reframe(meanRatio = mean(Ratio), sdRatio = sd(Ratio))

####### Add color index into data frames
bac.renew.rel.abundant.genera <- bac.renew.rel.Genus.TmStatus %>% group_by(Phylum,Genus) %>% reframe(Ratio = mean(meanRatio)) %>% arrange(desc(Ratio)) %>% filter(Ratio >= 0.005)
print(bac.renew.rel.abundant.genera, n=50)
bac.renew.rel.Genus.each.sample$Genus2 <- bac.renew.rel.Genus.each.sample$Genus
bac.renew.rel.Genus.TmStatus$Genus2 <- bac.renew.rel.Genus.TmStatus$Genus

bac.renew.rel.Genus.each.sample$Genus2[!(bac.renew.rel.Genus.each.sample$Genus%in% bac.renew.rel.abundant.genera$Genus)] <- "Others"
bac.renew.rel.Genus.TmStatus$Genus2[!(bac.renew.rel.Genus.TmStatus$Genus%in% bac.renew.rel.abundant.genera$Genus)] <- "Others"

ord.bac.renew.genus <- c(bac.renew.rel.abundant.genera$Genus,"Others")


order.rep <- c("TmM5","TmM4","TmM3","TmM2","TmM1","TmD1","TmD2","TmD3","TmD4","TmD5")
bac.renew.rel.Genus.each.sample$month <- factor(bac.renew.rel.Genus.each.sample$month, levels = order.month)
bac.renew.rel.Genus.each.sample$Replicate <- factor(bac.renew.rel.Genus.each.sample$Replicate, levels = order.rep)
bac.renew.rel.Genus.each.sample$Genus2 <- factor(bac.renew.rel.Genus.each.sample$Genus2, levels = rev(ord.bac.renew.genus))

bac.renew.rel.Genus.TmStatus$month <- factor(bac.renew.rel.Genus.TmStatus$month, levels = order.month)
bac.renew.rel.Genus.TmStatus$TmStatus <- factor(bac.renew.rel.Genus.TmStatus$TmStatus, levels = c("TmM","TmD"))
bac.renew.rel.Genus.TmStatus$Genus2 <- factor(bac.renew.rel.Genus.TmStatus$Genus2, levels = rev(ord.bac.renew.genus))
length(unique(bac.renew.rel.Genus.TmStatus$Genus2))
####### plotting
##### each replicate
colors_40 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1f78b4", "#fb9a99",
  "#33a02c", "#ff7f00", "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a",
  "#cab2d6", "#fdbf6f", "#ff7f7f", "#6b6ecf", "#c49c94", "#9edae5",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
  "#98df8a", "#ff9896", "#9467bd", "#c7c7c7", "#999999"
)
library(ggplot2)
ggplot(bac.renew.rel.Genus.each.sample, aes(x=Replicate, y = Ratio, fill = Genus2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~month, ncol=1)+
  #scale_fill_discrete() +
  scale_fill_manual(values = rev(colors_40)) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 2,reverse = T))+
  theme(legend.position="right") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


ggplot(bac.renew.rel.Genus.TmStatus, aes(x=month, y = meanRatio, fill = Genus2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~TmStatus, ncol=1)+
  #scale_fill_discrete() +
  scale_fill_manual(values = rev(colors_40)) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 2,reverse = T))+
  theme(legend.position="right") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

writexl::write_xlsx(bac.renew.rel.Genus.TmStatus,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_baceria.genus.xlsx")
writexl::write_xlsx(bac.renew.rel.Genus.each.sample,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_baceria.genus_allreps.xlsx")
writexl::write_xlsx(bac.renew.rel.phylum.TmStatus,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_baceria.phylum.xlsx")
writexl::write_xlsx(bac.renew.rel.phylum.each.sample,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_baceria.phylum_allreps.xlsx")


bac.renew.rel.supple.tab <- bac.renew.rel.tab %>% group_by(Replicate,TmStatus,month,Phylum,Class,Order,Family,Genus) %>% reframe(Ratio = sum(Abundance))
bac.renew.rel.supple.tab.wide <- reshape2::dcast(bac.renew.rel.supple.tab,Phylum+Class+Order+Family+Genus~month+Replicate, value.var = "Ratio")

writexl::write_xlsx(bac.renew.rel.supple.tab.wide,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/TableS_relAbund_baceria_allreps.xlsx" )


####### Fungi ####
tax.fun.tab<-data.frame(phyloseq::tax_table(fun.clean.ss.f))


tax.fun.tab$Phylum[is.na(tax.fun.tab$Phylum)] <- "Fungi_unknown"
unique(tax.fun.tab$Phylum)

tax.fun.tab$Class[is.na(tax.fun.tab$Class)] <- paste0(tax.fun.tab$Phylum[is.na(tax.fun.tab$Class)],"_unknown")
tax.fun.tab$Class <- gsub("unknown_unknown","unknown",tax.fun.tab$Class)
unique(tax.fun.tab$Class)

tax.fun.tab$Order[is.na(tax.fun.tab$Order)] <- paste0(tax.fun.tab$Class[is.na(tax.fun.tab$Order)],"_unknown")
tax.fun.tab$Order <- gsub("unknown_unknown","unknown",tax.fun.tab$Order)
unique(tax.fun.tab$Order)

tax.fun.tab$Family[is.na(tax.fun.tab$Family)] <- paste0(tax.fun.tab$Order[is.na(tax.fun.tab$Family)],"_unknown")
tax.fun.tab$Family <- gsub("unknown_unknown","unknown",tax.fun.tab$Family)
unique(tax.fun.tab$Family)

tax.fun.tab$Genus[is.na(tax.fun.tab$Genus)] <- paste0(tax.fun.tab$Family[is.na(tax.fun.tab$Genus)],"_unknown")
tax.fun.tab$Genus <- gsub("unknown_unknown","unknown",tax.fun.tab$Genus)
unique(tax.fun.tab$Genus)

fun.renew<-fun.clean.ss.f
phyloseq::tax_table(fun.renew) <- phyloseq::tax_table(as.matrix(tax.fun.tab))
fun.renew.rel <- microbiome::transform(fun.renew, "compositional")
fun.renew.rel.tab<-phyloseq::psmelt(fun.renew.rel)
names(fun.renew.rel.tab)[7] <- "TmStatus"
fun.renew.rel.tab$TmStatus <- ifelse(fun.renew.rel.tab$TmStatus == "yes","TmD","TmM")
head(fun.renew.rel.tab)
fun.renew.rel.tab$Replicate <- paste0(fun.renew.rel.tab$TmStatus,fun.renew.rel.tab$compartment)

######## Agglomerated at the phylum level
######### 1. Each sample level
unique(fun.renew.rel.tab$Phylum)
fun.renew.rel.phylum.each.sample <- fun.renew.rel.tab %>% group_by(Replicate,month,Phylum,Class) %>% reframe(Ratio = sum(Abundance))
head(fun.renew.rel.phylum.each.sample)

######### 2. TmStatus level
fun.renew.rel.phylum.TmStatus <- fun.renew.rel.tab %>% group_by(TmStatus,Replicate,month,Phylum,Class) %>% reframe(Ratio = sum(Abundance))
head(fun.renew.rel.phylum.TmStatus)
fun.renew.rel.phylum.TmStatus <- fun.renew.rel.phylum.TmStatus %>% group_by(TmStatus,month,Phylum,Class) %>% reframe(meanRatio = mean(Ratio), sdRatio = sd(Ratio))

####### Add color index into data frames
fun.renew.rel.abundant <- fun.renew.rel.phylum.TmStatus %>% group_by(Phylum,Class) %>% reframe(Ratio = mean(meanRatio)) %>% arrange(desc(Ratio)) %>% filter(Ratio >= 0.005)

fun.renew.rel.phylum.each.sample$Class2 <- fun.renew.rel.phylum.each.sample$Class
fun.renew.rel.phylum.TmStatus$Class2 <- fun.renew.rel.phylum.TmStatus$Class

fun.renew.rel.phylum.each.sample$Class2[!(fun.renew.rel.phylum.each.sample$Class%in% fun.renew.rel.abundant$Class)] <- "Others"
fun.renew.rel.phylum.TmStatus$Class2[!(fun.renew.rel.phylum.TmStatus$Class%in% fun.renew.rel.abundant$Class)] <- "Others"

ord.fun.renew.phyla <- c(fun.renew.rel.abundant$Class[c(1:7,9:12,8)],"Others")


fun.renew.rel.phylum.each.sample$month <- factor(fun.renew.rel.phylum.each.sample$month, levels = order.month)
fun.renew.rel.phylum.each.sample$Replicate <- factor(fun.renew.rel.phylum.each.sample$Replicate, levels = order.rep)
fun.renew.rel.phylum.each.sample$Class2 <- factor(fun.renew.rel.phylum.each.sample$Class2, levels = rev(ord.fun.renew.phyla))

fun.renew.rel.phylum.TmStatus$month <- factor(fun.renew.rel.phylum.TmStatus$month, levels = order.month)
fun.renew.rel.phylum.TmStatus$TmStatus <- factor(fun.renew.rel.phylum.TmStatus$TmStatus, levels = c("TmM","TmD"))
fun.renew.rel.phylum.TmStatus$Class2 <- factor(fun.renew.rel.phylum.TmStatus$Class2, levels = rev(ord.fun.renew.phyla))

####### plotting
##### each replicate
library(ggplot2)

ggplot(fun.renew.rel.phylum.each.sample, aes(x=Replicate, y = Ratio, fill = Class2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~month, ncol=1)+
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Agaricomycetes"="#ff6a72",
                               "Dothideomycetes"="#33ccff",
                               "Saccharomycetes"="#ffa633",
                               "Leotiomycetes"="#64b200",
                               "Eurotiomycetes"="#006699",
                               "Chytridiomycetes"="#00cc99",
                               "Sordariomycetes"="#cccc99",
                               "Mucoromycetes"="#00a6ff",
                               "Ascomycota_unknown"="#b385ff",
                               "Fungi_unknown"="#474747",
                               "Geminibasidiomycetes"="#54B5FB",
                               "Mortierellomycetes"="#F3C911",
                               "Pezizomycetes" ="#652926",
                               "Rozellomycotina_cls_Incertae_sedis"="#E5D1D0",
                               "Rozellomycota_cls_Incertae_sedis"="#F8BCBD",
                               "Umbelopsidomycetes"="#8B3D88" ,
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


ggplot(fun.renew.rel.phylum.TmStatus, aes(x=month, y = meanRatio, fill = Class2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~TmStatus, ncol=1)+
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Agaricomycetes"="#ff6a72",
                               "Dothideomycetes"="#33ccff",
                               "Saccharomycetes"="#ffa633",
                               "Leotiomycetes"="#64b200",
                               "Eurotiomycetes"="#006699",
                               "Chytridiomycetes"="#00cc99",
                               "Sordariomycetes"="#cccc99",
                               "Mucoromycetes"="#00a6ff",
                               "Ascomycota_unknown"="#b385ff",
                               "Fungi_unknown"="#474747",
                               "Geminibasidiomycetes"="#54B5FB",
                               "Mortierellomycetes"="#F3C911",
                               "Pezizomycetes" ="#652926",
                               "Rozellomycotina_cls_Incertae_sedis"="#E5D1D0",
                               "Rozellomycota_cls_Incertae_sedis"="#F8BCBD",
                               "Umbelopsidomycetes"="#8B3D88" ,
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



######## Agglomerated at the genus level
######### 1. Each sample level
unique(fun.renew.rel.tab$Genus)
library(dplyr)
fun.renew.rel.Genus.each.sample <- fun.renew.rel.tab %>% group_by(Replicate,month,Phylum,Genus) %>% reframe(Ratio = sum(Abundance))
head(fun.renew.rel.Genus.each.sample)

######### 2. TmStatus level
fun.renew.rel.Genus.TmStatus <- fun.renew.rel.tab %>% group_by(TmStatus,Replicate,month,Phylum,Genus) %>% reframe(Ratio = sum(Abundance))
head(fun.renew.rel.Genus.TmStatus)
fun.renew.rel.Genus.TmStatus <- fun.renew.rel.Genus.TmStatus %>% group_by(TmStatus,month,Phylum,Genus) %>% reframe(meanRatio = mean(Ratio), sdRatio = sd(Ratio))

####### Add color index into data frames
fun.renew.rel.abundant.genera <- fun.renew.rel.Genus.TmStatus %>% group_by(Phylum,Genus) %>% reframe(Ratio = mean(meanRatio)) %>% arrange(desc(Ratio)) %>% filter(Ratio >= 0.005)
print(fun.renew.rel.abundant.genera, n=50)
fun.renew.rel.Genus.each.sample$Genus2 <- fun.renew.rel.Genus.each.sample$Genus
fun.renew.rel.Genus.TmStatus$Genus2 <- fun.renew.rel.Genus.TmStatus$Genus

fun.renew.rel.Genus.each.sample$Genus2[!(fun.renew.rel.Genus.each.sample$Genus%in% fun.renew.rel.abundant.genera$Genus)] <- "Others"
fun.renew.rel.Genus.TmStatus$Genus2[!(fun.renew.rel.Genus.TmStatus$Genus%in% fun.renew.rel.abundant.genera$Genus)] <- "Others"


ord.fun.renew.genus <- c(fun.renew.rel.abundant.genera$Genus[c(1:7,9:36,8)],"Others")


fun.renew.rel.Genus.each.sample$month <- factor(fun.renew.rel.Genus.each.sample$month, levels = order.month)
fun.renew.rel.Genus.each.sample$Replicate <- factor(fun.renew.rel.Genus.each.sample$Replicate, levels = order.rep)
fun.renew.rel.Genus.each.sample$Genus2 <- factor(fun.renew.rel.Genus.each.sample$Genus2, levels = rev(ord.fun.renew.genus))

fun.renew.rel.Genus.TmStatus$month <- factor(fun.renew.rel.Genus.TmStatus$month, levels = order.month)
fun.renew.rel.Genus.TmStatus$TmStatus <- factor(fun.renew.rel.Genus.TmStatus$TmStatus, levels = c("TmM","TmD"))
fun.renew.rel.Genus.TmStatus$Genus2 <- factor(fun.renew.rel.Genus.TmStatus$Genus2, levels = rev(ord.fun.renew.genus))
length(unique(fun.renew.rel.Genus.TmStatus$Genus2))
####### plotting
##### each replicate

colors_fun.gen<- c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33",
  "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62", "#8da0cb",
  "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#3288bd",
  "#66c2a5", "#abdda4", "#e6f598", "#fee08b", "#fdae61", "#f46d43",
  "#d53e4f", "#9e0142", "#5e4fa2", "#3288bd", "#66c2a5", "#abdda4",
  "#e6f598", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#474747","#999999"
)

colors_40 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
  "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1f78b4", "#fb9a99",
  "#33a02c", "#ff7f00", "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a",
  "#cab2d6", "#fdbf6f", "#ff7f7f", "#6b6ecf", "#c49c94", "#9edae5",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
  "#98df8a", "#ff9896", "#9467bd", "#c7c7c7", "#999999"
)
library(ggplot2)
ggplot(fun.renew.rel.Genus.each.sample, aes(x=Replicate, y = Ratio, fill = Genus2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~month, ncol=1)+
  #scale_fill_discrete() +
  scale_fill_manual(values = rev(colors_fun.gen)) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 2,reverse = T))+
  theme(legend.position="right") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())


ggplot(fun.renew.rel.Genus.TmStatus, aes(x=month, y = meanRatio, fill = Genus2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~TmStatus, ncol=1)+
  #scale_fill_discrete() +
  scale_fill_manual(values = rev(colors_fun.gen)) +
  
  xlab('')+
  ylab("Relative abundance\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 2,reverse = T))+
  theme(legend.position="right") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

writexl::write_xlsx(fun.renew.rel.Genus.TmStatus,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_fungi.genus.xlsx")
writexl::write_xlsx(fun.renew.rel.Genus.each.sample,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_fungi.genus_allreps.xlsx")
writexl::write_xlsx(fun.renew.rel.phylum.TmStatus,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_fungi.class.xlsx")
writexl::write_xlsx(fun.renew.rel.phylum.each.sample,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/sourcedata_relAbund_fungi.class_allreps.xlsx")


fun.renew.rel.supple.tab <- fun.renew.rel.tab %>% group_by(Replicate,TmStatus,month,Phylum,Class,Order,Family,Genus) %>% reframe(Ratio = sum(Abundance))
fun.renew.rel.supple.tab.wide <- reshape2::dcast(fun.renew.rel.supple.tab,Phylum+Class+Order+Family+Genus~month+Replicate, value.var = "Ratio")

writexl::write_xlsx(fun.renew.rel.supple.tab.wide,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/TableS_relAbund_fungi_allreps.xlsx" )




#### Relative abundance table
melt.bac.rel <- phyloseq::psmelt(bac.renew.rel)
melt.fun.rel <- phyloseq::psmelt(fun.renew.rel)

names(melt.bac.rel)[7] <- "TmStatus"
melt.bac.rel$TmStatus <- ifelse(melt.bac.rel$TmStatus == "yes","TmD","TmM")
head(melt.bac.rel)
melt.bac.rel$Replicate <- paste0(melt.bac.rel$TmStatus,melt.bac.rel$compartment)

names(melt.fun.rel)[7] <- "TmStatus"
melt.fun.rel$TmStatus <- ifelse(melt.fun.rel$TmStatus == "yes","TmD","TmM")
head(melt.fun.rel)
melt.fun.rel$Replicate <- paste0(melt.fun.rel$TmStatus,melt.fun.rel$compartment)


melt.bac.rel$month <- factor(melt.bac.rel$month, levels = order.month)
melt.bac.rel$Replicate <- factor(melt.bac.rel$Replicate, levels = order.rep)
melt.fun.rel$month <- factor(melt.fun.rel$month, levels = order.month)
melt.fun.rel$Replicate <- factor(melt.fun.rel$Replicate, levels = order.rep)

###Bacteria
bac.rel.phylum <- melt.bac.rel %>% group_by(Replicate,month,TmStatus,Phylum) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.phylum,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial phyla and class_relative abundance_replicates.csv")

bac.rel.class <- melt.bac.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.class,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial class_relative abundance_replicates.csv")

bac.rel.order <- melt.bac.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class,Order) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.order,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial order_relative abundance_replicates.csv")

bac.rel.family <- melt.bac.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class,Order,Family) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.family,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial family_relative abundance_replicates.csv")

bac.rel.genus <- melt.bac.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class,Order,Family,Genus) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.genus,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial genus_relative abundance_replicates.csv")


### Abundance table for month data (mean and standard deviation)
bac.rel.phylum.month<-bac.rel.phylum %>% group_by(month, TmStatus, Phylum)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.phylum.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial phyla_relative abundance_month.csv")

bac.rel.class.month<-bac.rel.class %>% group_by(month, TmStatus, Phylum, Class)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.class.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial class_relative abundance_month.csv")

bac.rel.order.month<-bac.rel.order %>% group_by(month, TmStatus, Phylum, Class, Order)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.order.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial order_relative abundance_month.csv")

bac.rel.family.month<-bac.rel.family %>% group_by(month, TmStatus, Phylum, Class, Order, Family)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.family.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial family_relative abundance_month.csv")

bac.rel.genus.month<-bac.rel.genus %>% group_by(month, TmStatus, Phylum, Class, Order, Family, Genus)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.genus.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/bacterial genus_relative abundance_month.csv")

### fungi
fun.rel.phylum <- melt.fun.rel %>% group_by(Replicate,month,TmStatus,Phylum) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.phylum,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal phyla and class_relative abundance_replicates.csv")

fun.rel.class <- melt.fun.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.class,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal class_relative abundance_replicates.csv")

fun.rel.order <- melt.fun.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class,Order) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.order,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal order_relative abundance_replicates.csv")

fun.rel.family <- melt.fun.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class,Order,Family) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.family,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal family_relative abundance_replicates.csv")

fun.rel.genus <- melt.fun.rel %>% group_by(Replicate,month,TmStatus,Phylum,Class,Order,Family,Genus) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.genus,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal genus_relative abundance_replicates.csv")


### Abundance table for month data (mean and standard deviation)
fun.rel.phylum.month<-fun.rel.phylum %>% group_by(month, TmStatus, Phylum)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.phylum.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal phyla_relative abundance_month.csv")

fun.rel.class.month<-fun.rel.class %>% group_by(month, TmStatus, Phylum, Class)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.class.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal class_relative abundance_month.csv")

fun.rel.order.month<-fun.rel.order %>% group_by(month, TmStatus, Phylum, Class, Order)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.order.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal order_relative abundance_month.csv")

fun.rel.family.month<-fun.rel.family %>% group_by(month, TmStatus, Phylum, Class, Order, Family)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.family.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal family_relative abundance_month.csv")

fun.rel.genus.month<-fun.rel.genus %>% group_by(month, TmStatus, Phylum, Class, Order, Family, Genus)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.genus.month,"/Users/testaccount/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SNU/SNU_IHBae/Manuscript/Comprehensive/fungal genus_relative abundance_month.csv")

