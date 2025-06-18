########## Prove "Kill the Winner" hypothesis #########
load("./RData/bac.metaG.RPKM.RData")
load("./RData/phage.metaG.RPKM.RData")
########### Relative abundance normalized analysis ########
###### Bacterial contig abundance (transform to relative abundance)
load("../../TaxonomicAssignment/tax.reassign.bac.RData")
head(tax.reassign.bac)
bac.contig.act.read <- tax.reassign.bac[c(1,2,17,18,19:24)]

phage.list.all<-rownames(df.phage.metaG.RPKM)
phage.contig.list <- unique(gsub("\\|+\\w+","",phage.list.all))
bac.contig.act.read<-subset(bac.contig.act.read, !(contigName %in% phage.contig.list))
bac.contig.act.read$Year <- gsub("20","Y",bac.contig.act.read$Year)
bac.contig.act.read$Sample <- paste0(bac.contig.act.read$Year,bac.contig.act.read$Group2)
bac.contig.genus <- bac.contig.act.read %>%group_by(Phylum2,Class2,Order2,Family2,Genus2,Species2,Sample) %>% reframe(Count = sum(Count))
head(bac.contig.genus)

bac.contig.genus <- subset(bac.contig.genus, Sample!="Y23NC")
bac.contig.genus <- bac.contig.genus %>% group_by(Sample) %>% mutate(total = sum(Count))
head(bac.contig.genus)
bac.contig.genus$ratio <- bac.contig.genus$Count/bac.contig.genus$total

bac.contig.genus.tab <- bac.contig.genus %>% group_by(Phylum2,Genus2,Sample) %>% reframe(ratio = sum(ratio))
nrow(bac.contig.genus.tab)
head(bac.contig.genus.tab)

names(bac.contig.genus.tab) <- c("phylum","genus","variable","ratio")


bac.contig.genus.tab$phylum <- gsub("_Unidentified","_unknown",bac.contig.genus.tab$phylum)
bac.contig.genus.tab$genus <- gsub("_Unidentified","_unknown",bac.contig.genus.tab$genus)



##### 
vMAG.abund.act.read
vMAG.abund.act.read
df.phage.metaG.act.read <- data.frame(otu_table(vMAG.abund.act.read))
df.phage.metaG.act.read
df.phage.metaG.act.read$genome <- rownames(df.phage.metaG.act.read)
head(df.phage.metaG.act.read)
host.info.table.best.hit3

df.phage.metaG.act.read.host.info <- merge(df.phage.metaG.act.read,host.info.table.best.hit3,by = "genome")
nrow(df.phage.metaG.act.read.host.info)

##### Lytic and temperate
df.phage.metaG.act.read.host.info$type <- ifelse(df.phage.metaG.act.read.host.info$genome%in%lytic.phage.list,"lytic","temperate")

##### agglomerate phage abundance at the predicted host genus level
df.phage.metaG.act.read.host.info.melt <- reshape2::melt(df.phage.metaG.act.read.host.info)

df.phage.metaG.act.read.host.genus <- df.phage.metaG.act.read.host.info.melt %>% group_by(Host.phylum,Host.genus, type,variable) %>% reframe(act.read = sum(value))
head(df.phage.metaG.act.read.host.genus)
unique(df.phage.metaG.act.read.host.genus$Host.phylum)

df.phage.metaG.rel.host.genus <- df.phage.metaG.act.read.host.genus %>% group_by(variable) %>% mutate(total = sum(act.read))
df.phage.metaG.rel.host.genus$ratio <- df.phage.metaG.rel.host.genus$act.read/df.phage.metaG.rel.host.genus$total

df.phage.metaG.rel.host.genus.lytic <- df.phage.metaG.rel.host.genus[df.phage.metaG.rel.host.genus$type == "lytic",]
df.phage.metaG.rel.host.genus.temperate <- df.phage.metaG.rel.host.genus[df.phage.metaG.rel.host.genus$type == "temperate",]


##### Find consistent genera between bacterial abundance and predicted host abundance
matched.bac.genera.rel<-intersect(unique(bac.contig.genus.tab$genus), unique(df.phage.metaG.rel.host.genus$Host.genus))

##### Lytic
phage.metaG.rel.matched.lytic<-df.phage.metaG.rel.host.genus.lytic[df.phage.metaG.rel.host.genus.lytic$Host.genus%in%matched.bac.genera.rel,]
bac.metaG.rel.matched<-bac.contig.genus.tab[bac.contig.genus.tab$genus%in%matched.bac.genera.rel,]

phage.metaG.rel.matched.lytic.wide <- reshape2::dcast(phage.metaG.rel.matched.lytic,Host.genus~variable, value.var = "ratio")
rownames(phage.metaG.rel.matched.lytic.wide) <- phage.metaG.rel.matched.lytic.wide$Host.genus
phage.metaG.rel.matched.lytic.wide <- phage.metaG.rel.matched.lytic.wide[-c(1)]

bac.metaG.rel.matched.wide <- reshape2::dcast(bac.metaG.rel.matched,genus~variable, value.var = "ratio")
rownames(bac.metaG.rel.matched.wide) <- bac.metaG.rel.matched.wide$genus
bac.metaG.rel.matched.wide <- bac.metaG.rel.matched.wide[-c(1)]
bac.metaG.rel.matched.wide[is.na(bac.metaG.rel.matched.wide)] <- 0
##### Temperate
phage.metaG.rel.matched.temperate<-df.phage.metaG.rel.host.genus.temperate[df.phage.metaG.rel.host.genus.temperate$Host.genus%in%matched.bac.genera,]
phage.metaG.rel.matched.temperate.wide <- reshape2::dcast(phage.metaG.rel.matched.temperate,Host.genus~variable, value.var = "ratio")
rownames(phage.metaG.rel.matched.temperate.wide) <- phage.metaG.rel.matched.temperate.wide$Host.genus
phage.metaG.rel.matched.temperate.wide <- phage.metaG.rel.matched.temperate.wide[-c(1)]

####### Prove 1 #######
###### Mantel's correlation test
###### Community distance at the genus level (considering only matched genera)
dist.phage.metaG.rel <- vegan::vegdist(t(phage.metaG.rel.matched.wide), method = "canberra") 
dist.bac.metaG.rel <- vegan::vegdist(t(bac.metaG.rel.matched.wide), method = "canberra") 
mantel.bray.metaG<-vegan::mantel(dist.phage.metaG.rel, dist.bac.metaG.rel, method = "pearson")
## Mantel's R = 0.8981, P = 0.001

dist.phage.metaG.rel.lytic <- vegan::vegdist(t(phage.metaG.rel.matched.lytic.wide), method = "bray") 
mantel.bray.metaG.rel.lytic<-vegan::mantel(dist.phage.metaG.rel.lytic, dist.bac.metaG.rel, method = "spearman")
## Mantel's R = 0.8247, P = 0.001

dist.phage.metaG.rel.temperate <- vegan::vegdist(t(phage.metaG.rel.matched.temperate.wide), method = "bray") 
mantel.bray.metaG.rel.temperate<-vegan::mantel(dist.phage.metaG.rel.temperate, dist.bac.metaG.rel, method = "spearman")
## Mantel's R = 0.8089, P = 0.001


####### Prove 2 #######
##### linear regression test
phage.metaG.rel.matched.lytic
lm.tab.phage.rel.lytic <- phage.metaG.rel.matched.lytic
lm.tab.phage.rel.lytic <- lm.tab.phage.rel.lytic[-c(3,5,6)]
names(lm.tab.phage.rel.lytic) <- c("phylum","genus","sample","phage")

lm.tab.phage.rel.temperate <- phage.metaG.rel.matched.temperate
lm.tab.phage.rel.temperate <- lm.tab.phage.rel.temperate[-c(3,5,6)]
names(lm.tab.phage.rel.temperate) <- c("phylum","genus","sample","phage")

lm.tab.bac.rel <- bac.metaG.rel.matched
names(lm.tab.bac.rel) <- c("phylum","genus","sample","bacteria")

######### lytic and temperate
lm.tab.phage.rel.lytic.bac.genus<-merge(lm.tab.phage.rel.lytic,lm.tab.bac.rel, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.lytic.rel <- data.frame(genus =unique(lm.tab.phage.rel.lytic.bac.genus$genus),group ="1",pval=1,tval=1)
for (i in unique(lm.tab.phage.rel.lytic.bac.genus$genus)){
  test.res<- lm.tab.phage.rel.lytic.bac.genus[lm.tab.phage.rel.lytic.bac.genus$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(phage ~ bacteria, .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.lytic.rel$group[which(lm.res.lytic.rel$genus == i)] <- test.res$group
  lm.res.lytic.rel$pval[which(lm.res.lytic.rel$genus == i)] <- test.res$pval
  lm.res.lytic.rel$tval[which(lm.res.lytic.rel$genus == i)] <- test.res$tval
}



lm.tab.phage.temperate.bac.genus<-merge(lm.tab.phage.temperate,lm.tab.bac, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.temperate <- data.frame(genus =unique(lm.tab.phage.temperate.bac.genus$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.temperate.bac.genus$genus)){
  test.res<- lm.tab.phage.temperate.bac.genus[lm.tab.phage.temperate.bac.genus$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(phage ~ bacteria, .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.temperate$group[which(lm.res.temperate$genus == i)] <- test.res$group
  lm.res.temperate$pval[which(lm.res.temperate$genus == i)] <- test.res$pval
  lm.res.temperate$tval[which(lm.res.temperate$genus == i)] <- test.res$tval
}


########### RPKM normalized analysis ########
##### Add host prediction results to phage abundance table
vMAG.abund.RPKM
df.phage.metaG.RPKM <- data.frame(otu_table(vMAG.abund.RPKM))
df.phage.metaG.RPKM
df.phage.metaG.RPKM$genome <- rownames(df.phage.metaG.RPKM)
head(df.phage.metaG.RPKM)
host.info.table.best.hit3

df.phage.metaG.RPKM.host.info <- merge(df.phage.metaG.RPKM,host.info.table.best.hit3,by = "genome")
nrow(df.phage.metaG.RPKM.host.info)

##### Lytic and temperate
lytic.phage.list
df.phage.metaG.RPKM.host.info$type <- ifelse(df.phage.metaG.RPKM.host.info$genome%in%lytic.phage.list,"lytic","temperate")


##### agglomerate phage abundance at the predicted host genus level
df.phage.metaG.RPKM.host.info.melt <- reshape2::melt(df.phage.metaG.RPKM.host.info)

df.phage.metaG.RPKM.host.genus <- df.phage.metaG.RPKM.host.info.melt %>% group_by(Host.phylum,Host.genus, type,variable) %>% reframe(RPKM = sum(value))
head(df.phage.metaG.RPKM.host.genus)
unique(df.phage.metaG.RPKM.host.genus$Host.phylum)

df.phage.metaG.RPKM.host.genus.lytic <- df.phage.metaG.RPKM.host.genus[df.phage.metaG.RPKM.host.genus$type == "lytic",]
df.phage.metaG.RPKM.host.genus.temperate <- df.phage.metaG.RPKM.host.genus[df.phage.metaG.RPKM.host.genus$type == "temperate",]

##### Agglomerate bacterial abundance at the genus level
head(bac.metaG.RPKM)
bac.metaG.RPKM.melt <- reshape2::melt(bac.metaG.RPKM)
head(bac.metaG.RPKM.melt)
bac.metaG.RPKM.genus <- bac.metaG.RPKM.melt %>% group_by(phylum, genus,variable) %>% reframe(RPKM =sum(value))
head(bac.metaG.RPKM.genus)
unique(bac.metaG.RPKM.genus$genus)
unique(bac.metaG.RPKM.genus$phylum)
bac.metaG.RPKM.genus$phylum <- gsub("_Unclassified","_unknown",bac.metaG.RPKM.genus$phylum)
bac.metaG.RPKM.genus$genus <- gsub("_Unclassified","_unknown",bac.metaG.RPKM.genus$genus)

##### Find consistent genera between bacterial abundance and predicted host abundance
matched.bac.genera<-intersect(unique(bac.metaG.RPKM.genus$genus), unique(df.phage.metaG.RPKM.host.genus$Host.genus))

##### Lytic
phage.metaG.RPKM.matched.lytic<-df.phage.metaG.RPKM.host.genus.lytic[df.phage.metaG.RPKM.host.genus.lytic$Host.genus%in%matched.bac.genera,]
bac.metaG.RPKM.matched<-bac.metaG.RPKM.genus[bac.metaG.RPKM.genus$genus%in%matched.bac.genera,]

phage.metaG.RPKM.matched.lytic.wide <- reshape2::dcast(phage.metaG.RPKM.matched.lytic,Host.genus~variable, value.var = "RPKM")
rownames(phage.metaG.RPKM.matched.lytic.wide) <- phage.metaG.RPKM.matched.lytic.wide$Host.genus
phage.metaG.RPKM.matched.lytic.wide <- phage.metaG.RPKM.matched.lytic.wide[-c(1)]

bac.metaG.RPKM.matched.wide <- reshape2::dcast(bac.metaG.RPKM.matched,genus~variable, value.var = "RPKM")
rownames(bac.metaG.RPKM.matched.wide) <- bac.metaG.RPKM.matched.wide$genus
bac.metaG.RPKM.matched.wide <- bac.metaG.RPKM.matched.wide[-c(1)]

##### Temperate
phage.metaG.RPKM.matched.temperate<-df.phage.metaG.RPKM.host.genus.temperate[df.phage.metaG.RPKM.host.genus.temperate$Host.genus%in%matched.bac.genera,]
phage.metaG.RPKM.matched.temperate.wide <- reshape2::dcast(phage.metaG.RPKM.matched.temperate,Host.genus~variable, value.var = "RPKM")
rownames(phage.metaG.RPKM.matched.temperate.wide) <- phage.metaG.RPKM.matched.temperate.wide$Host.genus
phage.metaG.RPKM.matched.temperate.wide <- phage.metaG.RPKM.matched.temperate.wide[-c(1)]

####### Prove 1 #######
###### Mantel's correlation test
###### Community distance at the genus level (considering only matched genera)
dist.phage.metaG.RPKM <- vegan::vegdist(t(phage.metaG.RPKM.matched.wide), method = "canberra") 
dist.bac.metaG.RPKM <- vegan::vegdist(t(bac.metaG.RPKM.matched.wide), method = "canberra") 
mantel.bray.metaG<-vegan::mantel(dist.phage.metaG.RPKM, dist.bac.metaG.RPKM, method = "pearson")
## Mantel's R = 0.9345, P = 0.001

dist.phage.metaG.RPKM.lytic <- vegan::vegdist(t(phage.metaG.RPKM.matched.lytic.wide), method = "canberra") 
mantel.bray.metaG.lytic<-vegan::mantel(dist.phage.metaG.RPKM.lytic, dist.bac.metaG.RPKM, method = "pearson")
## Mantel's R = 0.9391, P = 0.001

dist.phage.metaG.RPKM.temperate <- vegan::vegdist(t(phage.metaG.RPKM.matched.temperate.wide), method = "canberra") 
mantel.bray.metaG.temperate<-vegan::mantel(dist.phage.metaG.RPKM.temperate, dist.bac.metaG.RPKM, method = "pearson")
## Mantel's R = 0.9288, P = 0.001


####### Prove 2 #######
##### linear regression test
phage.metaG.RPKM.matched.lytic
lm.tab.phage.lytic <- phage.metaG.RPKM.matched.lytic
lm.tab.phage.lytic <- lm.tab.phage.lytic[-c(3)]
names(lm.tab.phage.lytic) <- c("phylum","genus","sample","phage")

lm.tab.phage.temperate <- phage.metaG.RPKM.matched.temperate
lm.tab.phage.temperate <- lm.tab.phage.temperate[-c(3)]
names(lm.tab.phage.temperate) <- c("phylum","genus","sample","phage")

lm.tab.bac <- bac.metaG.RPKM.matched
names(lm.tab.bac) <- c("phylum","genus","sample","bacteria")

lm.tab.phage.bac.genus<-merge(lm.tab.phage,lm.tab.bac, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res <- data.frame(genus =unique(lm.tab.phage.bac.genus$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.bac.genus$genus)){
  test.res<- lm.tab.phage.bac.genus[lm.tab.phage.bac.genus$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(phage ~ bacteria, .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res$group[which(lm.res$genus == i)] <- test.res$group
  lm.res$pval[which(lm.res$genus == i)] <- test.res$pval
  lm.res$tval[which(lm.res$genus == i)] <- test.res$tval
}
nrow(lm.res[lm.res$group == "Positive",])
lm.res[lm.res$genus == "Bradyrhizobium",]

######### lytic and temperate
lm.tab.phage.lytic.bac.genus<-merge(lm.tab.phage.lytic,lm.tab.bac, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.lytic <- data.frame(genus =unique(lm.tab.phage.lytic.bac.genus$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.lytic.bac.genus$genus)){
  test.res<- lm.tab.phage.lytic.bac.genus[lm.tab.phage.lytic.bac.genus$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(phage ~ bacteria, .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.lytic$group[which(lm.res.lytic$genus == i)] <- test.res$group
  lm.res.lytic$pval[which(lm.res.lytic$genus == i)] <- test.res$pval
  lm.res.lytic$tval[which(lm.res.lytic$genus == i)] <- test.res$tval
}



lm.tab.phage.temperate.bac.genus<-merge(lm.tab.phage.temperate,lm.tab.bac, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.temperate <- data.frame(genus =unique(lm.tab.phage.temperate.bac.genus$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.temperate.bac.genus$genus)){
  test.res<- lm.tab.phage.temperate.bac.genus[lm.tab.phage.temperate.bac.genus$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(phage ~ bacteria, .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.temperate$group[which(lm.res.temperate$genus == i)] <- test.res$group
  lm.res.temperate$pval[which(lm.res.temperate$genus == i)] <- test.res$pval
  lm.res.temperate$tval[which(lm.res.temperate$genus == i)] <- test.res$tval
}

####### log-transformed
######### lytic and temperate
lm.tab.phage.lytic.bac.genus<-merge(lm.tab.phage.lytic,lm.tab.bac, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.lytic.log <- data.frame(genus =unique(lm.tab.phage.lytic.bac.genus$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.lytic.bac.genus$genus)){
  test.res<- lm.tab.phage.lytic.bac.genus[lm.tab.phage.lytic.bac.genus$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(log(phage+1) ~ log(bacteria+1), .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.lytic.log$group[which(lm.res.lytic.log$genus == i)] <- test.res$group
  lm.res.lytic.log$pval[which(lm.res.lytic.log$genus == i)] <- test.res$pval
  lm.res.lytic.log$tval[which(lm.res.lytic.log$genus == i)] <- test.res$tval
}



lm.tab.phage.temperate.bac.genus<-merge(lm.tab.phage.temperate,lm.tab.bac, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.temperate.log <- data.frame(genus =unique(lm.tab.phage.temperate.bac.genus$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.temperate.bac.genus$genus)){
  test.res<- lm.tab.phage.temperate.bac.genus[lm.tab.phage.temperate.bac.genus$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(log(phage+1) ~ log(bacteria+1), .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.temperate.log$group[which(lm.res.temperate.log$genus == i)] <- test.res$group
  lm.res.temperate.log$pval[which(lm.res.temperate.log$genus == i)] <- test.res$pval
  lm.res.temperate.log$tval[which(lm.res.temperate.log$genus == i)] <- test.res$tval
}



##### patterns of abundant genera (winner) 
#### Abundant bacterial genus
bac.contig.genus.meanRatio <- bac.contig.genus.tab %>% group_by(phylum,genus) %>% reframe(ratio =mean(ratio))%>% arrange(desc(ratio)) %>% head(31)
print(bac.contig.genus.meanRatio,n=40)
dorminant.genus<-bac.contig.genus.meanRatio$genus
dorminant.genus <- dorminant.genus[-which(dorminant.genus == "Bacteria_unknown")]

abundant.genera.lytic<-lm.res.lytic[lm.res.lytic$genus%in%dorminant.genus,]
abundant.genera.temperate<-lm.res.temperate[lm.res.temperate$genus%in%dorminant.genus,]

abundant.genera.lytic.log<-lm.res.lytic.log[lm.res.lytic.log$genus%in%dorminant.genus,]
abundant.genera.temperate.log<-lm.res.temperate.log[lm.res.temperate.log$genus%in%dorminant.genus,]




lm.tab.phage.temperate.bac.genus$type <- "temperate"
lm.tab.phage.lytic.bac.genus$type <- "lytic"

dorm.lm.tab.phage.lytic.bac.genus<-lm.tab.phage.lytic.bac.genus[lm.tab.phage.lytic.bac.genus$genus %in% dorminant.genus,]
dorm.lm.tab.phage.temperate.bac.genus<-lm.tab.phage.temperate.bac.genus[lm.tab.phage.temperate.bac.genus$genus %in% dorminant.genus,]

abundant.genera.lytic.log$type <- "lytic"
abundant.genera.temperate.log$type <- "temperate"

linetype.tab <- rbind(abundant.genera.lytic.log,abundant.genera.temperate.log)
linetype.tab$shape <- ifelse(linetype.tab$group=="No trend","P > 0.05","P < 0.05")
linetype.tab <- linetype.tab[c(1,5,6)]


dorm.lm.res.tab<-rbind(dorm.lm.tab.phage.lytic.bac.genus,dorm.lm.tab.phage.temperate.bac.genus)
dorm.lm.res.tab$phage <- log10(dorm.lm.res.tab$phage+1)
dorm.lm.res.tab$bacteria <- log10(dorm.lm.res.tab$bacteria+1)

dorm.lm.res.tab.plot<-merge(dorm.lm.res.tab,linetype.tab, by= c("genus"="genus","type"="type"))
head(dorm.lm.res.tab.plot)

dorm.lm.res.tab.plot$TmColonizationStatus <- ifelse(dorm.lm.res.tab.plot$sample %in% c("Y22TmD1","Y22TmD2","Y23TmD1","Y23TmM1"),"Active",
                                                    ifelse(dorm.lm.res.tab.plot$sample%in%c("Y22TmD3","Y22TmD4","Y22TmD5","Y23TmD2","Y23TmD3","Y23TmD4","Y23TmD5"),"Post","Pre"))

dorm.lm.res.tab.plot$TmColonizationStatus <- factor(dorm.lm.res.tab.plot$TmColonizationStatus,levels =c("Pre","Active","Post"))
ggplot(dorm.lm.res.tab.plot,aes(x=bacteria,y=phage,color=TmColonizationStatus))+ geom_point(aes(shape=type))+
  geom_smooth(method = "lm",se=T,aes(color=type,linetype=shape))+
  scale_color_manual(values = c("Pre" = "#CDC673", "Active" = "#104E8B","Post" = "#CD7054",
                                "lytic"="#9B59B6","temperate"="#8FBC8F"))+
  scale_linetype_manual(values = c("P < 0.05" = "solid", "P > 0.05" = "dashed"))+
  scale_shape_manual(values = c("lytic" = 16, "temperate" = 1))+
  facet_wrap(~genus, scales="free")

#### Select top 16 genera for visualization
dorminant.genus

unique(abundant.genera.lytic.log$genus)
unique(abundant.genera.temperate.log$genus)

lytic.dorminant.genus<-dorminant.genus[which(dorminant.genus %in% unique(abundant.genera.lytic.log$genus))]
temperate.dorminant.genus<-dorminant.genus[which(dorminant.genus %in% unique(abundant.genera.temperate.log$genus))]
dorm.genus.for.plot<-head(unique(c(lytic.dorminant.genus,temperate.dorminant.genus)),16)

dorm.lm.res.tab.plot.final<-dorm.lm.res.tab.plot[dorm.lm.res.tab.plot$genus%in%dorm.genus.for.plot,]
dorm.lm.res.tab.plot.final$genus <- factor(dorm.lm.res.tab.plot.final$genus, levels = dorm.genus.for.plot)
ggplot(dorm.lm.res.tab.plot.final,aes(x=bacteria,y=phage,color=TmColonizationStatus))+ geom_point(aes(shape=type))+
  geom_smooth(method = "lm",se=T,aes(color=type,linetype=shape))+
  scale_color_manual(values = c("Pre" = "#CDC673", "Active" = "#104E8B","Post" = "#CD7054",
                                "lytic"="#9B59B6","temperate"="#8FBC8F"))+
  scale_linetype_manual(values = c("P < 0.05" = "solid", "P > 0.05" = "dashed"))+
  scale_shape_manual(values = c("lytic" = 16, "temperate" = 1))+
  facet_wrap(~genus, scales="free")+theme(aspect.ratio = 1)+
  #theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 8,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 8,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5,size=8, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=8, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

writexl::write_xlsx(dorm.lm.res.tab.plot.final,"./Source data_abundance correlation_Phage and host_250228.xlsx")

lm.res.lytic.log$type <- "lytic"
lm.res.temperate.log$type <- "temperate"

lm.res.tab.log<-rbind(lm.res.lytic.log,lm.res.temperate.log)
writexl::write_xlsx(lm.res.tab.log,"./Supple data_linear regression result_Phage and host_250228.xlsx")


unique(dorm.lm.res.tab.plot.final$genus)



############## metatranscriptome-based abundance #########
expr.bac.abund<-read.delim("./RNAabundance/metag.rpkm.prokaryote.rna.abundance_bwa.tsv",sep="\t",header =T)
names(expr.bac.abund) <- c("Contig","Y23TmD1","Y23TmD4","Y23TmD5","Y23TmM1","Y23TmM4","Y23TmM5")
head(expr.bac.abund)
nrow(expr.bac.abund)
expr.bac.abund <- expr.bac.abund[-c(8)]
expr.bac.abund.filtered<-expr.bac.abund[rowSums(expr.bac.abund[-c(1)]) > 0,]
nrow(expr.bac.abund.filtered) #5464027

load(file= "../RepProtein_taxonomy/diamond.res.all.contigs_final.RData")
diamond.res.all.bac <- diamond.res.all[diamond.res.all$Domain=="Bacteria",]
nrow(diamond.res.all.bac)

diamond.res.all.bac$Phylum[which(diamond.res.all.bac$Phylum == "Unidentified")] <- "Bacteria_unknown"
unique(diamond.res.all.bac$Phylum)

diamond.res.all.bac$Class[which(diamond.res.all.bac$Class == "Unidentified" &
                                                      diamond.res.all.bac$Phylum == "Bacteria_unknown")] <- "Bacteria_unknown"

diamond.res.all.bac$Class[which(diamond.res.all.bac$Class == "Unidentified" &
                                                      diamond.res.all.bac$Phylum != "Bacteria_unknown")] <- paste0(diamond.res.all.bac$Phylum[which(diamond.res.all.bac$Class == "Unidentified" &
                                                                                                                                                                                              diamond.res.all.bac$Phylum != "Bacteria_unknown")],"_unknown")
unique(diamond.res.all.bac$Class)


diamond.res.all.bac$Order[which(diamond.res.all.bac$Order == "Unidentified" &
                                                      diamond.res.all.bac$Class == "Bacteria_unknown"&
                                                      diamond.res.all.bac$Phylum == "Bacteria_unknown")] <- "Bacteria_unknown"

diamond.res.all.bac$Order[which(diamond.res.all.bac$Order == "Unidentified" &
                                                      diamond.res.all.bac$Class != "Bacteria_unknown"&
                                                      grepl("_unknown",diamond.res.all.bac$Class))] <- diamond.res.all.bac$Class[which(diamond.res.all.bac$Order == "Unidentified" &
                                                                                                                                                                                 diamond.res.all.bac$Class != "Bacteria_unknown"&
                                                                                                                                                                                 grepl("_unknown",diamond.res.all.bac$Class))]

diamond.res.all.bac$Order[which(diamond.res.all.bac$Order == "Unidentified")] <- paste0(diamond.res.all.bac$Class[which(diamond.res.all.bac$Order == "Unidentified")],"_unknown")


unique(diamond.res.all.bac$Order)



diamond.res.all.bac$Family[which(diamond.res.all.bac$Family == "Unidentified" &
                                                       diamond.res.all.bac$Order == "Bacteria_unknown" &
                                                       diamond.res.all.bac$Class == "Bacteria_unknown"&
                                                       diamond.res.all.bac$Phylum == "Bacteria_unknown")] <- "Bacteria_unknown"

diamond.res.all.bac$Family[which(diamond.res.all.bac$Family == "Unidentified" &
                                                       diamond.res.all.bac$Order != "Bacteria_unknown" )] <- paste0(diamond.res.all.bac$Order[which(diamond.res.all.bac$Family == "Unidentified" &
                                                                                                                                                                                              diamond.res.all.bac$Order != "Bacteria_unknown" )],"_unknown")

diamond.res.all.bac$Family <- gsub("_unknown_unknown","_unknown",diamond.res.all.bac$Family)



diamond.res.all.bac$Genus[which(diamond.res.all.bac$Genus == "Unidentified" &
                                                      diamond.res.all.bac$Family == "Bacteria_unknown" &
                                                      diamond.res.all.bac$Order == "Bacteria_unknown" &
                                                      diamond.res.all.bac$Class == "Bacteria_unknown"&
                                                      diamond.res.all.bac$Phylum == "Bacteria_unknown")] <- "Bacteria_unknown"

diamond.res.all.bac$Genus[which(diamond.res.all.bac$Genus == "Unidentified" &
                                                      diamond.res.all.bac$Family != "Bacteria_unknown" )] <- paste0(diamond.res.all.bac$Family[which(diamond.res.all.bac$Genus == "Unidentified" &
                                                                                                                                                                                               diamond.res.all.bac$Family != "Bacteria_unknown" )],"_unknown")

diamond.res.all.bac$Genus <- gsub("_unknown_unknown","_unknown",diamond.res.all.bac$Genus)
unique(diamond.res.all.bac$Genus)
names(diamond.res.all.bac)

bac.tax <-unique(diamond.res.all.bac[c(1,7:12)])
nrow(bac.tax)
names(bac.tax)[1] <- "Contig"
expr.bac.abund.filtered2 <- merge(bac.tax,expr.bac.abund.filtered, by ="Contig")
nrow(expr.bac.abund.filtered2)

expr.bac.abund.filtered2.melt <- reshape2::melt(expr.bac.abund.filtered2)
bac.metaT.RPKM.genus <- expr.bac.abund.filtered2.melt%>%group_by(Phylum,Class,Order,Family, Genus,variable) %>% reframe(RPKM=sum(value))
head(bac.metaT.RPKM.genus)


### bac
load(file = "~/Desktop/SNU/KNU_Prof.HangilKim/ForestVirome/DNAvirus/host_phageome/host.RPKM.rna.RData")

head(df.phage.metaT.RPKM.species.summary)


df.phage.metaT.RPKM<-data.frame(otu_table(vMAG.expr.RPKM))

df.phage.metaT.RPKM$genome <- rownames(df.phage.metaT.RPKM)
head(df.phage.metaT.RPKM)
host.info.table.best.hit3

df.phage.metaT.RPKM.host.info <- merge(df.phage.metaT.RPKM,host.info.table.best.hit3,by = "genome")
nrow(df.phage.metaT.RPKM.host.info)

##### Lytic and temperate
lytic.phage.list
df.phage.metaT.RPKM.host.info$type <- ifelse(df.phage.metaT.RPKM.host.info$genome%in%lytic.phage.list,"lytic","temperate")

##### agglomerate phage abundance at the predicted host genus level
df.phage.metaT.RPKM.host.info.melt <- reshape2::melt(df.phage.metaT.RPKM.host.info)
head(df.phage.metaT.RPKM.host.info.melt)
df.phage.metaT.RPKM.host.genus <- df.phage.metaT.RPKM.host.info.melt %>% group_by(Host.phylum,Host.genus, type,variable) %>% reframe(RPKM = sum(value))
head(df.phage.metaT.RPKM.host.genus)
unique(df.phage.metaT.RPKM.host.genus$Host.phylum)

df.phage.metaT.RPKM.host.genus.lytic <- df.phage.metaT.RPKM.host.genus[df.phage.metaT.RPKM.host.genus$type == "lytic",]
df.phage.metaT.RPKM.host.genus.temperate <- df.phage.metaT.RPKM.host.genus[df.phage.metaT.RPKM.host.genus$type == "temperate",]

##### Agglomerate bacterial abundance at the genus level
unique(bac.metaT.RPKM.genus$Genus)
unique(bac.metaT.RPKM.genus$Phylum)
##### Find consistent genera between bacterial abundance and predicted host abundance
matched.bac.genera.metaT<-intersect(unique(bac.metaT.RPKM.genus$Genus), unique(df.phage.metaT.RPKM.host.genus$Host.genus))

##### Lytic
phage.metaT.RPKM.matched.lytic<-df.phage.metaT.RPKM.host.genus.lytic[df.phage.metaT.RPKM.host.genus.lytic$Host.genus%in%matched.bac.genera.metaT,]
bac.metaT.RPKM.matched<-bac.metaT.RPKM.genus[bac.metaT.RPKM.genus$Genus%in%matched.bac.genera.metaT,]

phage.metaT.RPKM.matched.lytic.wide <- reshape2::dcast(phage.metaT.RPKM.matched.lytic,Host.genus~variable, value.var = "RPKM")
rownames(phage.metaT.RPKM.matched.lytic.wide) <- phage.metaT.RPKM.matched.lytic.wide$Host.genus
phage.metaT.RPKM.matched.lytic.wide <- phage.metaT.RPKM.matched.lytic.wide[-c(1)]

bac.metaT.RPKM.matched.wide <- reshape2::dcast(bac.metaT.RPKM.matched,Genus~variable, value.var = "RPKM")
rownames(bac.metaT.RPKM.matched.wide) <- bac.metaT.RPKM.matched.wide$Genus
bac.metaT.RPKM.matched.wide <- bac.metaT.RPKM.matched.wide[-c(1)]

##### Temperate
phage.metaT.RPKM.matched.temperate<-df.phage.metaT.RPKM.host.genus.temperate[df.phage.metaT.RPKM.host.genus.temperate$Host.genus%in%matched.bac.genera,]
phage.metaT.RPKM.matched.temperate.wide <- reshape2::dcast(phage.metaT.RPKM.matched.temperate,Host.genus~variable, value.var = "RPKM")
rownames(phage.metaT.RPKM.matched.temperate.wide) <- phage.metaT.RPKM.matched.temperate.wide$Host.genus
phage.metaT.RPKM.matched.temperate.wide <- phage.metaT.RPKM.matched.temperate.wide[-c(1)]

phage.metaT.RPKM.matched.wide <- rbind(phage.metaT.RPKM.matched.lytic.wide,phage.metaT.RPKM.matched.temperate.wide)
####### Prove 1 #######
###### Mantel's correlation test
###### Community distance at the genus level (considering only matched genera)
dist.phage.metaT.RPKM <- vegan::vegdist(t(phage.metaT.RPKM.matched.wide), method = "canberra") 
dist.bac.metaT.RPKM <- vegan::vegdist(t(bac.metaT.RPKM.matched.wide), method = "canberra") 
mantel.bray.metaT<-vegan::mantel(dist.phage.metaT.RPKM, dist.bac.metaT.RPKM, method = "pearson")
## Mantel's R = 0.478, P = 0.0013889 

dist.phage.metaT.RPKM.lytic <- vegan::vegdist(t(phage.metaT.RPKM.matched.lytic.wide), method = "canberra") 
mantel.bray.metaT.lytic<-vegan::mantel(dist.phage.metaT.RPKM.lytic, dist.bac.metaT.RPKM, method = "pearson")
# 0.4691, 0.0041667
dist.phage.metaT.RPKM.temperate <- vegan::vegdist(t(phage.metaT.RPKM.matched.temperate.wide), method = "canberra") 
mantel.bray.metaT.temperate<-vegan::mantel(dist.phage.metaT.RPKM.temperate, dist.bac.metaT.RPKM, method = "pearson")
# 0.4866 0.0027778

####### Prove 2 #######
##### linear regression test
phage.metaT.RPKM.matched.lytic
lm.tab.phage.lytic.metaT <- phage.metaT.RPKM.matched.lytic
lm.tab.phage.lytic.metaT <- lm.tab.phage.lytic.metaT[-c(3)]
names(lm.tab.phage.lytic.metaT) <- c("phylum","genus","sample","phage")

lm.tab.phage.temperate.metaT <- phage.metaT.RPKM.matched.temperate
lm.tab.phage.temperate.metaT <- lm.tab.phage.temperate.metaT[-c(3)]
names(lm.tab.phage.temperate.metaT) <- c("phylum","genus","sample","phage")

lm.tab.bac.metaT <- bac.metaT.RPKM.matched[c(1,5,6,7)]
names(lm.tab.bac.metaT) <- c("phylum","genus","sample","bacteria")

#lm.tab.phage.bac.genus<-merge(lm.tab.phage,lm.tab.bac, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

# lm.res <- data.frame(genus =unique(lm.tab.phage.bac.genus$genus,group ="1",pval=1,tval=1))
# for (i in unique(lm.tab.phage.bac.genus$genus)){
#   test.res<- lm.tab.phage.bac.genus[lm.tab.phage.bac.genus$genus==i,] %>% group_by(genus) %>% 
#     do(mod = summary(lm(phage ~ bacteria, .))) %>%
#     mutate(
#       pval = mod$coefficients[[8]],
#       tval = mod$coefficients[[6]],
#       group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
#       group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
#     select(genus, group, pval, tval)
#   
#   lm.res$group[which(lm.res$genus == i)] <- test.res$group
#   lm.res$pval[which(lm.res$genus == i)] <- test.res$pval
#   lm.res$tval[which(lm.res$genus == i)] <- test.res$tval
# }
# nrow(lm.res[lm.res$group == "Positive",])
# lm.res[lm.res$genus == "Bradyrhizobium",]

######### lytic and temperate
lm.tab.phage.lytic.bac.genus.metaT<-merge(lm.tab.phage.lytic.metaT,lm.tab.bac.metaT, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.lytic.metaT <- data.frame(genus =unique(lm.tab.phage.lytic.bac.genus.metaT$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.lytic.bac.genus.metaT$genus)){
  test.res<- lm.tab.phage.lytic.bac.genus.metaT[lm.tab.phage.lytic.bac.genus.metaT$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(phage ~ bacteria, .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.lytic.metaT$group[which(lm.res.lytic.metaT$genus == i)] <- test.res$group
  lm.res.lytic.metaT$pval[which(lm.res.lytic.metaT$genus == i)] <- test.res$pval
  lm.res.lytic.metaT$tval[which(lm.res.lytic.metaT$genus == i)] <- test.res$tval
}



lm.tab.phage.temperate.bac.genus.metaT<-merge(lm.tab.phage.temperate.metaT,lm.tab.bac.metaT, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.temperate.metaT <- data.frame(genus =unique(lm.tab.phage.temperate.bac.genus.metaT$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.temperate.bac.genus.metaT$genus)){
  test.res<- lm.tab.phage.temperate.bac.genus.metaT[lm.tab.phage.temperate.bac.genus.metaT$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(phage ~ bacteria, .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.temperate.metaT$group[which(lm.res.temperate.metaT$genus == i)] <- test.res$group
  lm.res.temperate.metaT$pval[which(lm.res.temperate.metaT$genus == i)] <- test.res$pval
  lm.res.temperate.metaT$tval[which(lm.res.temperate.metaT$genus == i)] <- test.res$tval
}

####### log-transformed
######### lytic and temperate
lm.tab.phage.lytic.bac.genus.metaT<-merge(lm.tab.phage.lytic.metaT,lm.tab.bac.metaT, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.lytic.log.metaT <- data.frame(genus =unique(lm.tab.phage.lytic.bac.genus.metaT$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.lytic.bac.genus.metaT$genus)){
  test.res<- lm.tab.phage.lytic.bac.genus.metaT[lm.tab.phage.lytic.bac.genus.metaT$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(log(phage+1) ~ log(bacteria+1), .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.lytic.log.metaT$group[which(lm.res.lytic.log.metaT$genus == i)] <- test.res$group
  lm.res.lytic.log.metaT$pval[which(lm.res.lytic.log.metaT$genus == i)] <- test.res$pval
  lm.res.lytic.log.metaT$tval[which(lm.res.lytic.log.metaT$genus == i)] <- test.res$tval
}



lm.tab.phage.temperate.bac.genus.metaT<-merge(lm.tab.phage.temperate.metaT,lm.tab.bac.metaT, by= c("phylum"="phylum","genus"="genus","sample"="sample"))

lm.res.temperate.log.metaT <- data.frame(genus =unique(lm.tab.phage.temperate.bac.genus.metaT$genus,group ="1",pval=1,tval=1))
for (i in unique(lm.tab.phage.temperate.bac.genus.metaT$genus)){
  test.res<- lm.tab.phage.temperate.bac.genus.metaT[lm.tab.phage.temperate.bac.genus.metaT$genus==i,] %>% group_by(genus) %>% 
    do(mod = summary(lm(log(phage+1) ~ log(bacteria+1), .))) %>%
    mutate(
      pval = mod$coefficients[[8]],
      tval = mod$coefficients[[6]],
      group = ifelse(pval < 0.05 & tval < 0, 'Negative', 'No trend'),
      group = ifelse(pval < 0.05 & tval > 0, 'Positive', group)) %>%
    select(genus, group, pval, tval)
  
  lm.res.temperate.log.metaT$group[which(lm.res.temperate.log.metaT$genus == i)] <- test.res$group
  lm.res.temperate.log.metaT$pval[which(lm.res.temperate.log.metaT$genus == i)] <- test.res$pval
  lm.res.temperate.log.metaT$tval[which(lm.res.temperate.log.metaT$genus == i)] <- test.res$tval
}



##### patterns of abundant genera (winner) 
#### Abundant bacterial genus
abundant.genera.lytic.log.metaT<-lm.res.lytic.log.metaT[lm.res.lytic.log.metaT$genus%in%dorminant.genus,]
abundant.genera.temperate.log.metaT<-lm.res.temperate.log.metaT[lm.res.temperate.log.metaT$genus%in%dorminant.genus,]

lm.tab.phage.temperate.bac.genus.metaT$type <- "temperate"
lm.tab.phage.lytic.bac.genus.metaT$type <- "lytic"

dorm.lm.tab.phage.lytic.bac.genus.metaT<-lm.tab.phage.lytic.bac.genus.metaT[lm.tab.phage.lytic.bac.genus.metaT$genus %in% dorminant.genus,]
dorm.lm.tab.phage.temperate.bac.genus.metaT<-lm.tab.phage.temperate.bac.genus.metaT[lm.tab.phage.temperate.bac.genus.metaT$genus %in% dorminant.genus,]

abundant.genera.lytic.log.metaT$type <- "lytic"
abundant.genera.temperate.log.metaT$type <- "temperate"

linetype.tab.metaT <- rbind(abundant.genera.lytic.log.metaT,abundant.genera.temperate.log.metaT)
linetype.tab.metaT$shape <- ifelse(linetype.tab.metaT$group=="No trend","P > 0.05","P < 0.05")
linetype.tab.metaT <- linetype.tab.metaT[c(1,5,6)]


dorm.lm.res.tab.metaT<-rbind(dorm.lm.tab.phage.lytic.bac.genus.metaT,dorm.lm.tab.phage.temperate.bac.genus.metaT)
dorm.lm.res.tab.metaT$phage <- log10(dorm.lm.res.tab.metaT$phage+1)
dorm.lm.res.tab.metaT$bacteria <- log10(dorm.lm.res.tab.metaT$bacteria+1)

dorm.lm.res.tab.plot.metaT<-merge(dorm.lm.res.tab.metaT,linetype.tab.metaT, by= c("genus"="genus","type"="type"))
head(dorm.lm.res.tab.plot.metaT)

dorm.lm.res.tab.plot.metaT$TmColonizationStatus <- ifelse(dorm.lm.res.tab.plot.metaT$sample %in% c("Y23TmD1","Y23TmM1"),"Active",
                                                    ifelse(dorm.lm.res.tab.plot.metaT$sample%in%c("Y23TmD4","Y23TmD5"),"Post","Pre"))

dorm.lm.res.tab.plot.metaT$TmColonizationStatus <- factor(dorm.lm.res.tab.plot.metaT$TmColonizationStatus,levels =c("Pre","Active","Post"))
ggplot(dorm.lm.res.tab.plot.metaT,aes(x=bacteria,y=phage,color=TmColonizationStatus))+ geom_point(aes(shape=type))+
  geom_smooth(method = "lm",se=T,aes(color=type,linetype=shape))+
  scale_color_manual(values = c("Pre" = "#CDC673", "Active" = "#104E8B","Post" = "#CD7054",
                                "lytic"="#9B59B6","temperate"="#8FBC8F"))+
  scale_linetype_manual(values = c("P < 0.05" = "solid", "P > 0.05" = "dashed"))+
  scale_shape_manual(values = c("lytic" = 16, "temperate" = 1))+
  facet_wrap(~genus, scales="free")

#### Select top 12 genera for visualization
dorminant.genus

unique(abundant.genera.lytic.log$genus)
unique(abundant.genera.temperate.log$genus)

lytic.dorminant.genus.metaT<-dorminant.genus[which(dorminant.genus %in% unique(abundant.genera.lytic.log.metaT$genus))]
temperate.dorminant.genus.metaT<-dorminant.genus[which(dorminant.genus %in% unique(abundant.genera.temperate.log.metaT$genus))]
dorm.genus.for.plot.metaT<-head(unique(c(lytic.dorminant.genus.metaT,temperate.dorminant.genus.metaT)),16)

dorm.lm.res.tab.plot.final.metaT<-dorm.lm.res.tab.plot.metaT[dorm.lm.res.tab.plot.metaT$genus%in%dorm.genus.for.plot.metaT,]
dorm.lm.res.tab.plot.final.metaT$genus <- factor(dorm.lm.res.tab.plot.final.metaT$genus, levels = dorm.genus.for.plot.metaT)
ggplot(dorm.lm.res.tab.plot.final.metaT,aes(x=bacteria,y=phage,color=TmColonizationStatus))+ geom_point(aes(shape=type))+
  geom_smooth(method = "lm",se=T,aes(color=type,linetype=shape))+
  scale_color_manual(values = c("Pre" = "#CDC673", "Active" = "#104E8B","Post" = "#CD7054",
                                "lytic"="#9B59B6","temperate"="#8FBC8F"))+
  scale_linetype_manual(values = c("P < 0.05" = "solid", "P > 0.05" = "dashed"))+
  scale_shape_manual(values = c("lytic" = 16, "temperate" = 1))+
  facet_wrap(~genus, scales="free")+theme(aspect.ratio = 1)+
  #theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 8,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 8,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5,size=8, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=8, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

writexl::write_xlsx(dorm.lm.res.tab.plot.final,"./Source data_abundance correlation_Phage and host_250228.xlsx")

lm.res.lytic.log$type <- "lytic"
lm.res.temperate.log$type <- "temperate"

lm.res.tab.log<-rbind(lm.res.lytic.log,lm.res.temperate.log)
writexl::write_xlsx(lm.res.tab.log,"./Supple data_linear regression result_Phage and host_250228.xlsx")


