##### Refining diamond taxonomic assignment results
#### For metagenome data
TaxIdLineageTab.nr <- read.csv("/data/SongyiMetagenome/ncbi_lineages_2023-11-06.csv", fill=T)
head(TaxIdLineageTab.nr)

TaxIdLineageTab.nr$no.rank[which(TaxIdLineageTab.nr$no.rank == "")] <- "Unidentified"
TaxIdLineageTab.nr$superkingdom[which(TaxIdLineageTab.nr$superkingdom == "")] <- "Unidentified"
TaxIdLineageTab.nr$kingdom[which(TaxIdLineageTab.nr$kingdom == "")] <- "Unidentified"
TaxIdLineageTab.nr$phylum[which(TaxIdLineageTab.nr$phylum == "")] <- "Unidentified"
TaxIdLineageTab.nr$class[which(TaxIdLineageTab.nr$class == "")] <- "Unidentified"
TaxIdLineageTab.nr$order[which(TaxIdLineageTab.nr$order == "")] <- "Unidentified"
TaxIdLineageTab.nr$family[which(TaxIdLineageTab.nr$family == "")] <- "Unidentified"
TaxIdLineageTab.nr$genus[which(TaxIdLineageTab.nr$genus == "")] <- "Unidentified"
TaxIdLineageTab.nr$species[which(TaxIdLineageTab.nr$species == "")] <- "Unidentified"
TaxIdLineageTab.nr$tax_id <- as.character(TaxIdLineageTab.nr$tax_id)
head(TaxIdLineageTab.nr)


#### for individually assembled samples
CallDiamondNRRes <- function(SampleName){
  diamond.nr<-read.table(paste0("./6_tax/Diamond/",SampleName,"/","contig_nr_diamond_tax.tsv"), sep ='\t', header = F)
  names(diamond.nr) <- c("contigName","TaxID", "Evalue")
  return(diamond.nr)
}

#CreateInputTaxID <- function(diamond.nr, SampleName){
#  diamond.taxid<-diamond.nr$TaxID[which(!(diamond.nr$TaxID %in% c("0","1")))]
#  diamond.taxid<- unique(diamond.taxid)
#  write.table(diamond.taxid,paste0("./6_tax/Diamond/",SampleName,"/input_primaryTaxID_NR.tsv"), sep = '\t',row.names = F, col.names = F)
#}

###Assign primary tax id
DiamondLineageInfoPrep <- function(diamond.nr,SampleName){
#  primary.tax.id <- read.table(paste0("./6_tax/Diamond/",SampleName,"/tax_report.txt"),sep="|", header = T, quote = "")
#  diamond.nr$Primary.TaxID <- 0
#  for (i in unique(diamond.nr$TaxID)){
#    diamond.nr$Primary.TaxID[which(diamond.nr$TaxID == i)]<-ifelse(diamond.nr$TaxID[which(diamond.nr$TaxID == i)]=="0","0", 
#                                                                   ifelse(diamond.nr$TaxID[which(diamond.nr$TaxID == i)]=="1", "1", primary.tax.id$primary.taxid[which(primary.tax.id$taxid == i)]))
#  }
  
  diamond.nr.classified <- subset(diamond.nr, TaxID != "0")
  diamond.nr.unclassified <- subset(diamond.nr, TaxID == "0")
  
  TaxIdLineageTab.nr.sub <- TaxIdLineageTab.nr[c(1,40,38,2,3,4,5,6,7,8)]
  print(colnames( TaxIdLineageTab.nr.sub))
  names(TaxIdLineageTab.nr.sub)[1] <- "TaxID"
  diamond.nr.classified.2 <-merge(diamond.nr.classified, TaxIdLineageTab.nr.sub, by = c("TaxID"="TaxID"))
  diamond.nr.classified.2 <- diamond.nr.classified.2 [c(2,1,3,4,5,6,7,8,9,10,11,12)]
  print(colnames(diamond.nr.classified.2))
  diamond.nr.unclassified$no.rank <- "Unidentified"
  diamond.nr.unclassified$superkingdom <- "Unidentified"
  diamond.nr.unclassified$kingdom <- "Unidentified"
  diamond.nr.unclassified$phylum <- "Unidentified"
  diamond.nr.unclassified$class <- "Unidentified"
  diamond.nr.unclassified$order <- "Unidentified"
  diamond.nr.unclassified$family <- "Unidentified"
  diamond.nr.unclassified$genus <- "Unidentified"
  diamond.nr.unclassified$species <- "Unidentified"
  print(colnames(diamond.nr.unclassified))
  diamond.nr.lineage<-rbind(diamond.nr.classified.2, diamond.nr.unclassified)
  print(nrow(diamond.nr))
  print(nrow(diamond.nr.lineage))
  write.csv(diamond.nr.lineage,paste0("./6_tax/Diamond/",SampleName,"/diamond.contig.lineage.csv"))
  return(diamond.nr.lineage)

}



#### Running
### loading data
##2022 data
diamond.nr.9D1PG<-CallDiamondNRRes("9D1PG")
diamond.nr.9D2PG<-CallDiamondNRRes("9D2PG")
diamond.nr.9D3PH<-CallDiamondNRRes("9D3PH")
diamond.nr.9D4PH<-CallDiamondNRRes("9D4PH")
diamond.nr.9D5PG<-CallDiamondNRRes("9D5PG")

diamond.nr.9D1NH<-CallDiamondNRRes("9D1NH")
diamond.nr.9D2NG<-CallDiamondNRRes("9D2NG")
diamond.nr.9D3NG<-CallDiamondNRRes("9D3NG")
diamond.nr.9D4NG<-CallDiamondNRRes("9D4NG")
diamond.nr.9D5NG<-CallDiamondNRRes("9D5NG")

##2023 data
diamond.nr.23D1PA<-CallDiamondNRRes("23D1PA")
diamond.nr.23D2PB<-CallDiamondNRRes("23D2PB")
diamond.nr.23D3PA<-CallDiamondNRRes("23D3PA")
diamond.nr.23D4PA<-CallDiamondNRRes("23D4PA")
diamond.nr.23D5PB<-CallDiamondNRRes("23D5PB")

diamond.nr.23D1NA<-CallDiamondNRRes("23D1NA")
diamond.nr.23D2NA<-CallDiamondNRRes("23D2NA")
diamond.nr.23D3NA<-CallDiamondNRRes("23D3NA")
diamond.nr.23D4NA<-CallDiamondNRRes("23D4NA")
diamond.nr.23D5NA<-CallDiamondNRRes("23D5NA")

diamond.nr.23DNCA<-CallDiamondNRRes("23DNCA")

### Coassembled samples
### 2022+2023
diamond.nr.Dominant<-CallDiamondNRRes("Dominant")
diamond.nr.Minor<-CallDiamondNRRes("Minor")


#### Assing lineage info to diamond res
DiamondLineageInfoPrep(diamond.nr.9D1PG,"9D1PG")
DiamondLineageInfoPrep(diamond.nr.9D2PG,"9D2PG")
DiamondLineageInfoPrep(diamond.nr.9D3PH,"9D3PH")
DiamondLineageInfoPrep(diamond.nr.9D4PH,"9D4PH")
DiamondLineageInfoPrep(diamond.nr.9D5PG,"9D5PG")

DiamondLineageInfoPrep(diamond.nr.9D1NH,"9D1NH")
DiamondLineageInfoPrep(diamond.nr.9D2NG,"9D2NG")
DiamondLineageInfoPrep(diamond.nr.9D3NG,"9D3NG")
DiamondLineageInfoPrep(diamond.nr.9D4NG,"9D4NG")
DiamondLineageInfoPrep(diamond.nr.9D5NG,"9D5NG")

DiamondLineageInfoPrep(diamond.nr.23D1PA,"23D1PA")
DiamondLineageInfoPrep(diamond.nr.23D2PB,"23D2PB")
DiamondLineageInfoPrep(diamond.nr.23D3PA,"23D3PA")
DiamondLineageInfoPrep(diamond.nr.23D4PA,"23D4PA")
DiamondLineageInfoPrep(diamond.nr.23D5PB,"23D5PB")
DiamondLineageInfoPrep(diamond.nr.23D1NA,"23D1NA")
DiamondLineageInfoPrep(diamond.nr.23D2NA,"23D2NA")
DiamondLineageInfoPrep(diamond.nr.23D3NA,"23D3NA")
DiamondLineageInfoPrep(diamond.nr.23D4NA,"23D4NA")
DiamondLineageInfoPrep(diamond.nr.23D5NA,"23D5NA")

DiamondLineageInfoPrep(diamond.nr.23DNCA,"23DNCA")

diamond.nr.lineage.Dominant <- DiamondLineageInfoPrep(diamond.nr.Dominant,"Dominant")
diamond.nr.lineage.Minor <- DiamondLineageInfoPrep(diamond.nr.Minor,"Minor")

diamond.nr<-diamond.nr.Dominant
primary.tax.id <- read.table(paste0("./6_tax/Diamond/","Dominant","/tax_report.txt"),sep="|", header = T, quote = "")
  diamond.nr$Primary.TaxID <- 0
  for (i in unique(diamond.nr$TaxID)){
    diamond.nr$Primary.TaxID[which(diamond.nr$TaxID == i)]<-ifelse(diamond.nr$TaxID[which(diamond.nr$TaxID == i)]=="0","0", 
                                                                   ifelse(diamond.nr$TaxID[which(diamond.nr$TaxID == i)]=="1", "1", primary.tax.id$primary.taxid[which(primary.tax.id$taxid == i)]))
  }
  
  diamond.nr$ExRank <- 0
  diamond.nr$Domain <- 0
  diamond.nr$Kingdom <- 0
  diamond.nr$Phylum <- 0
  diamond.nr$Class <- 0
  diamond.nr$Order <- 0
  diamond.nr$Family <- 0
  diamond.nr$Genus <- 0
  diamond.nr$Species <- 0
  
  diamond.nr.classified <- subset(diamond.nr, TaxID != "0")
  diamond.nr.unclassified <- subset(diamond.nr, TaxID == "0")
  
colnames(TaxIdLineageTab.nr)
TaxIdLineageTab.nr.sub <- c(1,40,38,2,3,4,5,6,7,8)
names(TaxIdLineageTab.nr.sub)[1] <- "TaxID"
test2<-merge(diamond.nr.classified, TaxIdLineageTab.nr.sub, by = c("TaxID"="TaxID"))
nrow(test)
for (i in unique(as.character(diamond.nr.classified$Primary.TaxID))){
    diamond.nr.classified$ExRank[which(diamond.nr.classified$Primary.TaxID== i)] <- TaxIdLineageTab.nr$no.rank[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Domain[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$superkingdom[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Kingdom[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$kingdom[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Phylum[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$phylum[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Class[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$class[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Order[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$order[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Family[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$family[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Genus[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$genus[which(TaxIdLineageTab.nr$tax_id == i)]
    diamond.nr.classified$Species[which(diamond.nr.classified$Primary.TaxID == i)] <- TaxIdLineageTab.nr$species[which(TaxIdLineageTab.nr$tax_id == i)]
  }
  
  diamond.nr.unclassified$Domain <- "Unidentified"
  diamond.nr.unclassified$Kingdom <- "Unidentified"
  diamond.nr.unclassified$Phylum <- "Unidentified"
  diamond.nr.unclassified$Class <- "Unidentified"
  diamond.nr.unclassified$Order <- "Unidentified"
  diamond.nr.unclassified$Family <- "Unidentified"
  diamond.nr.unclassified$Genus <- "Unidentified"
  diamond.nr.unclassified$Species <- "Unidentified"
  
  diamond.nr.lineage<-rbind(diamond.nr.classified, diamond.nr.unclassified)

##### separate bacteria and fungi (#### in the linux terminal)
##loading
diamond.nr.lineage.9D1PG <- read.csv("./6_tax/Diamond/9D1PG/diamond.contig.lineage.csv")
diamond.nr.lineage.9D2PG <- read.csv("./6_tax/Diamond/9D2PG/diamond.contig.lineage.csv")
diamond.nr.lineage.9D3PH <- read.csv("./6_tax/Diamond/9D3PH/diamond.contig.lineage.csv")
diamond.nr.lineage.9D4PH <- read.csv("./6_tax/Diamond/9D4PH/diamond.contig.lineage.csv")
diamond.nr.lineage.9D5PG <- read.csv("./6_tax/Diamond/9D5PG/diamond.contig.lineage.csv")

diamond.nr.lineage.9D1NH <- read.csv("./6_tax/Diamond/9D1NH/diamond.contig.lineage.csv")
diamond.nr.lineage.9D2NG <- read.csv("./6_tax/Diamond/9D2NG/diamond.contig.lineage.csv")
diamond.nr.lineage.9D3NG <- read.csv("./6_tax/Diamond/9D3NG/diamond.contig.lineage.csv")
diamond.nr.lineage.9D4NG <- read.csv("./6_tax/Diamond/9D4NG/diamond.contig.lineage.csv")
diamond.nr.lineage.9D5NG <- read.csv("./6_tax/Diamond/9D5NG/diamond.contig.lineage.csv")

diamond.nr.lineage.23D1PA <- read.csv("./6_tax/Diamond/23D1PA/diamond.contig.lineage.csv")
diamond.nr.lineage.23D2PB <- read.csv("./6_tax/Diamond/23D2PB/diamond.contig.lineage.csv")
diamond.nr.lineage.23D3PA <- read.csv("./6_tax/Diamond/23D3PA/diamond.contig.lineage.csv")
diamond.nr.lineage.23D4PA <- read.csv("./6_tax/Diamond/23D4PA/diamond.contig.lineage.csv")
diamond.nr.lineage.23D5PB <- read.csv("./6_tax/Diamond/23D5PB/diamond.contig.lineage.csv")

diamond.nr.lineage.23D1NA <- read.csv("./6_tax/Diamond/23D1NA/diamond.contig.lineage.csv")
diamond.nr.lineage.23D2NA <- read.csv("./6_tax/Diamond/23D2NA/diamond.contig.lineage.csv")
diamond.nr.lineage.23D3NA <- read.csv("./6_tax/Diamond/23D3NA/diamond.contig.lineage.csv")
diamond.nr.lineage.23D4NA <- read.csv("./6_tax/Diamond/23D4NA/diamond.contig.lineage.csv")
diamond.nr.lineage.23D5NA <- read.csv("./6_tax/Diamond/23D5NA/diamond.contig.lineage.csv")

diamond.nr.lineage.23DNCA <- read.csv("./6_tax/Diamond/23DNCA/diamond.contig.lineage.csv")
diamond.nr.lineage.dominant <- read.csv("./6_tax/Diamond/Dominant/diamond.contig.lineage.csv")
diamond.nr.lineage.minor <- read.csv("./6_tax/Diamond/Minor/diamond.contig.lineage.csv")


diamond.nr.lineage.9D1PG<-diamond.nr.lineage.9D1PG[-c(1,5)]
diamond.nr.lineage.9D2PG<-diamond.nr.lineage.9D2PG[-c(1,5)]
diamond.nr.lineage.9D3PH<-diamond.nr.lineage.9D3PH[-c(1,5)]
diamond.nr.lineage.9D4PH<-diamond.nr.lineage.9D4PH[-c(1,5)]
diamond.nr.lineage.9D5PG<-diamond.nr.lineage.9D5PG[-c(1,5)]
diamond.nr.lineage.9D1NH<-diamond.nr.lineage.9D1NH[-c(1,5)]
diamond.nr.lineage.9D2NG<-diamond.nr.lineage.9D2NG[-c(1,5)]
diamond.nr.lineage.9D3NG<-diamond.nr.lineage.9D3NG[-c(1,5)]
diamond.nr.lineage.9D4NG<-diamond.nr.lineage.9D4NG[-c(1,5)]
diamond.nr.lineage.9D5NG<-diamond.nr.lineage.9D5NG[-c(1,5)]
diamond.nr.lineage.23D1PA<-diamond.nr.lineage.23D1PA[-c(1,5)]
diamond.nr.lineage.23D2PB<-diamond.nr.lineage.23D2PB[-c(1,5)]
diamond.nr.lineage.23D3PA<-diamond.nr.lineage.23D3PA[-c(1,5)]
diamond.nr.lineage.23D4PA<-diamond.nr.lineage.23D4PA[-c(1,5)]
diamond.nr.lineage.23D5PB<-diamond.nr.lineage.23D5PB[-c(1,5)]
diamond.nr.lineage.23D1NA<-diamond.nr.lineage.23D1NA[-c(1,5)]
diamond.nr.lineage.23D2NA<-diamond.nr.lineage.23D2NA[-c(1,5)]
diamond.nr.lineage.23D3NA<-diamond.nr.lineage.23D3NA[-c(1,5)]
diamond.nr.lineage.23D4NA<-diamond.nr.lineage.23D4NA[-c(1,5)]
diamond.nr.lineage.23D5NA<-diamond.nr.lineage.23D5NA[-c(1,5)]
diamond.nr.lineage.23DNCA<-diamond.nr.lineage.23DNCA[-c(1,5)]
diamond.nr.lineage.dominant<-diamond.nr.lineage.dominant[-c(1)]
diamond.nr.lineage.minor<-diamond.nr.lineage.minor[-c(1)]


names(diamond.nr.lineage.dominant)[4:12] <- c("ExRank","Kingdom","Domain","Phylum","Class","Order","Family","Genus","Species")
names(diamond.nr.lineage.minor)[4:12] <- c("ExRank","Kingdom","Domain","Phylum","Class","Order","Family","Genus","Species")

diamond.nr.lineage.dominant <- diamond.nr.lineage.dominant[c(1,2,3,4,6,5,7:12)]
diamond.nr.lineage.minor <- diamond.nr.lineage.minor[c(1,2,3,4,6,5,7:12)]
diamond.res.all<-rbind(diamond.nr.lineage.9D1PG,diamond.nr.lineage.9D2PG,diamond.nr.lineage.9D3PH,diamond.nr.lineage.9D4PH,diamond.nr.lineage.9D5PG,
      diamond.nr.lineage.9D1NH,diamond.nr.lineage.9D2PG,diamond.nr.lineage.9D2NG,diamond.nr.lineage.9D4NG,diamond.nr.lineage.9D5NG,
      diamond.nr.lineage.23D1PA,diamond.nr.lineage.23D2PB,diamond.nr.lineage.23D3PA,diamond.nr.lineage.23D4PA,diamond.nr.lineage.23D5PB,
      diamond.nr.lineage.23D1NA,diamond.nr.lineage.23D2NA,diamond.nr.lineage.23D3NA,diamond.nr.lineage.23D4NA,diamond.nr.lineage.23D5NA,
      diamond.nr.lineage.23DNCA,diamond.nr.lineage.dominant,diamond.nr.lineage.minor)



##### Check the counts of contigs classified prokaryotes, eukaryotes, and unidentified
### make a summary table
Summ.tab <- data.frame(Sample=c("9D1PG","9D2PG","9D3PH","9D4PH","9D5PG","9D1NH","9D2NG","9D3NG","9D4NG","9D5NG",
                                "23D1PA","23D2PB","23D3PA","23D4PA","23D5PB","23D1NA","23D2NA","23D3NA","23D4NA","23D5NA","23DNCA"),Bacteria = NA, Archaea = NA, Fungi = NA, OtherEukaryotes = NA, Total = NA)

Summ.tab$Total[which(Summ.tab$Sample == "9D1PG")]<-nrow(diamond.nr.lineage.9D1PG)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D1PG")]<-nrow(diamond.nr.lineage.9D1PG[diamond.nr.lineage.9D1PG$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D1PG")]<-nrow(diamond.nr.lineage.9D1PG[diamond.nr.lineage.9D1PG$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D1PG")]<-nrow(diamond.nr.lineage.9D1PG[diamond.nr.lineage.9D1PG$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D1PG")]<-nrow(diamond.nr.lineage.9D1PG[diamond.nr.lineage.9D1PG$Domain == "Eukaryota" & diamond.nr.lineage.9D1PG$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D2PG")]<-nrow(diamond.nr.lineage.9D2PG)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D2PG")]<-nrow(diamond.nr.lineage.9D2PG[diamond.nr.lineage.9D2PG$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D2PG")]<-nrow(diamond.nr.lineage.9D2PG[diamond.nr.lineage.9D2PG$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D2PG")]<-nrow(diamond.nr.lineage.9D2PG[diamond.nr.lineage.9D2PG$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D2PG")]<-nrow(diamond.nr.lineage.9D2PG[diamond.nr.lineage.9D2PG$Domain == "Eukaryota" & diamond.nr.lineage.9D2PG$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D3PH")]<-nrow(diamond.nr.lineage.9D3PH)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D3PH")]<-nrow(diamond.nr.lineage.9D3PH[diamond.nr.lineage.9D3PH$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D3PH")]<-nrow(diamond.nr.lineage.9D3PH[diamond.nr.lineage.9D3PH$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D3PH")]<-nrow(diamond.nr.lineage.9D3PH[diamond.nr.lineage.9D3PH$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D3PH")]<-nrow(diamond.nr.lineage.9D3PH[diamond.nr.lineage.9D3PH$Domain == "Eukaryota" & diamond.nr.lineage.9D3PH$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D4PH")]<-nrow(diamond.nr.lineage.9D4PH)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D4PH")]<-nrow(diamond.nr.lineage.9D4PH[diamond.nr.lineage.9D4PH$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D4PH")]<-nrow(diamond.nr.lineage.9D4PH[diamond.nr.lineage.9D4PH$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D4PH")]<-nrow(diamond.nr.lineage.9D4PH[diamond.nr.lineage.9D4PH$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D4PH")]<-nrow(diamond.nr.lineage.9D4PH[diamond.nr.lineage.9D4PH$Domain == "Eukaryota" & diamond.nr.lineage.9D4PH$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D5PG")]<-nrow(diamond.nr.lineage.9D5PG)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D5PG")]<-nrow(diamond.nr.lineage.9D5PG[diamond.nr.lineage.9D5PG$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D5PG")]<-nrow(diamond.nr.lineage.9D5PG[diamond.nr.lineage.9D5PG$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D5PG")]<-nrow(diamond.nr.lineage.9D5PG[diamond.nr.lineage.9D5PG$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D5PG")]<-nrow(diamond.nr.lineage.9D5PG[diamond.nr.lineage.9D5PG$Domain == "Eukaryota" & diamond.nr.lineage.9D5PG$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D1NH")]<-nrow(diamond.nr.lineage.9D1NH)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D1NH")]<-nrow(diamond.nr.lineage.9D1NH[diamond.nr.lineage.9D1NH$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D1NH")]<-nrow(diamond.nr.lineage.9D1NH[diamond.nr.lineage.9D1NH$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D1NH")]<-nrow(diamond.nr.lineage.9D1NH[diamond.nr.lineage.9D1NH$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D1NH")]<-nrow(diamond.nr.lineage.9D1NH[diamond.nr.lineage.9D1NH$Domain == "Eukaryota" & diamond.nr.lineage.9D1NH$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D2NG")]<-nrow(diamond.nr.lineage.9D2NG)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D2NG")]<-nrow(diamond.nr.lineage.9D2NG[diamond.nr.lineage.9D2NG$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D2NG")]<-nrow(diamond.nr.lineage.9D2NG[diamond.nr.lineage.9D2NG$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D2NG")]<-nrow(diamond.nr.lineage.9D2NG[diamond.nr.lineage.9D2NG$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D2NG")]<-nrow(diamond.nr.lineage.9D2NG[diamond.nr.lineage.9D2NG$Domain == "Eukaryota" & diamond.nr.lineage.9D2NG$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D3NG")]<-nrow(diamond.nr.lineage.9D3NG)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D3NG")]<-nrow(diamond.nr.lineage.9D3NG[diamond.nr.lineage.9D3NG$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D3NG")]<-nrow(diamond.nr.lineage.9D3NG[diamond.nr.lineage.9D3NG$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D3NG")]<-nrow(diamond.nr.lineage.9D3NG[diamond.nr.lineage.9D3NG$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D3NG")]<-nrow(diamond.nr.lineage.9D3NG[diamond.nr.lineage.9D3NG$Domain == "Eukaryota" & diamond.nr.lineage.9D3NG$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D4NG")]<-nrow(diamond.nr.lineage.9D4NG)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D4NG")]<-nrow(diamond.nr.lineage.9D4NG[diamond.nr.lineage.9D4NG$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D4NG")]<-nrow(diamond.nr.lineage.9D4NG[diamond.nr.lineage.9D4NG$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D4NG")]<-nrow(diamond.nr.lineage.9D4NG[diamond.nr.lineage.9D4NG$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D4NG")]<-nrow(diamond.nr.lineage.9D4NG[diamond.nr.lineage.9D4NG$Domain == "Eukaryota" & diamond.nr.lineage.9D4NG$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "9D5NG")]<-nrow(diamond.nr.lineage.9D5NG)
Summ.tab$Bacteria[which(Summ.tab$Sample == "9D5NG")]<-nrow(diamond.nr.lineage.9D5NG[diamond.nr.lineage.9D5NG$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "9D5NG")]<-nrow(diamond.nr.lineage.9D5NG[diamond.nr.lineage.9D5NG$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "9D5NG")]<-nrow(diamond.nr.lineage.9D5NG[diamond.nr.lineage.9D5NG$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "9D5NG")]<-nrow(diamond.nr.lineage.9D5NG[diamond.nr.lineage.9D5NG$Domain == "Eukaryota" & diamond.nr.lineage.9D5NG$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D1PA")]<-nrow(diamond.nr.lineage.23D1PA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D1PA")]<-nrow(diamond.nr.lineage.23D1PA[diamond.nr.lineage.23D1PA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D1PA")]<-nrow(diamond.nr.lineage.23D1PA[diamond.nr.lineage.23D1PA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D1PA")]<-nrow(diamond.nr.lineage.23D1PA[diamond.nr.lineage.23D1PA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D1PA")]<-nrow(diamond.nr.lineage.23D1PA[diamond.nr.lineage.23D1PA$Domain == "Eukaryota" & diamond.nr.lineage.23D1PA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D2PB")]<-nrow(diamond.nr.lineage.23D2PB)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D2PB")]<-nrow(diamond.nr.lineage.23D2PB[diamond.nr.lineage.23D2PB$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D2PB")]<-nrow(diamond.nr.lineage.23D2PB[diamond.nr.lineage.23D2PB$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D2PB")]<-nrow(diamond.nr.lineage.23D2PB[diamond.nr.lineage.23D2PB$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D2PB")]<-nrow(diamond.nr.lineage.23D2PB[diamond.nr.lineage.23D2PB$Domain == "Eukaryota" & diamond.nr.lineage.23D2PB$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D3PA")]<-nrow(diamond.nr.lineage.23D3PA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D3PA")]<-nrow(diamond.nr.lineage.23D3PA[diamond.nr.lineage.23D3PA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D3PA")]<-nrow(diamond.nr.lineage.23D3PA[diamond.nr.lineage.23D3PA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D3PA")]<-nrow(diamond.nr.lineage.23D3PA[diamond.nr.lineage.23D3PA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D3PA")]<-nrow(diamond.nr.lineage.23D3PA[diamond.nr.lineage.23D3PA$Domain == "Eukaryota" & diamond.nr.lineage.23D3PA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D4PA")]<-nrow(diamond.nr.lineage.23D4PA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D4PA")]<-nrow(diamond.nr.lineage.23D4PA[diamond.nr.lineage.23D4PA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D4PA")]<-nrow(diamond.nr.lineage.23D4PA[diamond.nr.lineage.23D4PA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D4PA")]<-nrow(diamond.nr.lineage.23D4PA[diamond.nr.lineage.23D4PA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D4PA")]<-nrow(diamond.nr.lineage.23D4PA[diamond.nr.lineage.23D4PA$Domain == "Eukaryota" & diamond.nr.lineage.23D4PA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D5PB")]<-nrow(diamond.nr.lineage.23D5PB)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D5PB")]<-nrow(diamond.nr.lineage.23D5PB[diamond.nr.lineage.23D5PB$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D5PB")]<-nrow(diamond.nr.lineage.23D5PB[diamond.nr.lineage.23D5PB$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D5PB")]<-nrow(diamond.nr.lineage.23D5PB[diamond.nr.lineage.23D5PB$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D5PB")]<-nrow(diamond.nr.lineage.23D5PB[diamond.nr.lineage.23D5PB$Domain == "Eukaryota" & diamond.nr.lineage.23D5PB$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D1NA")]<-nrow(diamond.nr.lineage.23D1NA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D1NA")]<-nrow(diamond.nr.lineage.23D1NA[diamond.nr.lineage.23D1NA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D1NA")]<-nrow(diamond.nr.lineage.23D1NA[diamond.nr.lineage.23D1NA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D1NA")]<-nrow(diamond.nr.lineage.23D1NA[diamond.nr.lineage.23D1NA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D1NA")]<-nrow(diamond.nr.lineage.23D1NA[diamond.nr.lineage.23D1NA$Domain == "Eukaryota" & diamond.nr.lineage.23D1NA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D2NA")]<-nrow(diamond.nr.lineage.23D2NA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D2NA")]<-nrow(diamond.nr.lineage.23D2NA[diamond.nr.lineage.23D2NA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D2NA")]<-nrow(diamond.nr.lineage.23D2NA[diamond.nr.lineage.23D2NA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D2NA")]<-nrow(diamond.nr.lineage.23D2NA[diamond.nr.lineage.23D2NA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D2NA")]<-nrow(diamond.nr.lineage.23D2NA[diamond.nr.lineage.23D2NA$Domain == "Eukaryota" & diamond.nr.lineage.23D2NA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D3NA")]<-nrow(diamond.nr.lineage.23D3NA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D3NA")]<-nrow(diamond.nr.lineage.23D3NA[diamond.nr.lineage.23D3NA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D3NA")]<-nrow(diamond.nr.lineage.23D3NA[diamond.nr.lineage.23D3NA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D3NA")]<-nrow(diamond.nr.lineage.23D3NA[diamond.nr.lineage.23D3NA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D3NA")]<-nrow(diamond.nr.lineage.23D3NA[diamond.nr.lineage.23D3NA$Domain == "Eukaryota" & diamond.nr.lineage.23D3NA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D4NA")]<-nrow(diamond.nr.lineage.23D4NA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D4NA")]<-nrow(diamond.nr.lineage.23D4NA[diamond.nr.lineage.23D4NA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D4NA")]<-nrow(diamond.nr.lineage.23D4NA[diamond.nr.lineage.23D4NA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D4NA")]<-nrow(diamond.nr.lineage.23D4NA[diamond.nr.lineage.23D4NA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D4NA")]<-nrow(diamond.nr.lineage.23D4NA[diamond.nr.lineage.23D4NA$Domain == "Eukaryota" & diamond.nr.lineage.23D4NA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23D5NA")]<-nrow(diamond.nr.lineage.23D5NA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23D5NA")]<-nrow(diamond.nr.lineage.23D5NA[diamond.nr.lineage.23D5NA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23D5NA")]<-nrow(diamond.nr.lineage.23D5NA[diamond.nr.lineage.23D5NA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23D5NA")]<-nrow(diamond.nr.lineage.23D5NA[diamond.nr.lineage.23D5NA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23D5NA")]<-nrow(diamond.nr.lineage.23D5NA[diamond.nr.lineage.23D5NA$Domain == "Eukaryota" & diamond.nr.lineage.23D5NA$Kingdom != "Fungi",])

Summ.tab$Total[which(Summ.tab$Sample == "23DNCA")]<-nrow(diamond.nr.lineage.23DNCA)
Summ.tab$Bacteria[which(Summ.tab$Sample == "23DNCA")]<-nrow(diamond.nr.lineage.23DNCA[diamond.nr.lineage.23DNCA$Domain == "Bacteria",])
Summ.tab$Archaea[which(Summ.tab$Sample == "23DNCA")]<-nrow(diamond.nr.lineage.23DNCA[diamond.nr.lineage.23DNCA$Domain == "Archaea",])
Summ.tab$Fungi[which(Summ.tab$Sample == "23DNCA")]<-nrow(diamond.nr.lineage.23DNCA[diamond.nr.lineage.23DNCA$Kingdom == "Fungi",])
Summ.tab$OtherEukaryotes[which(Summ.tab$Sample == "23DNCA")]<-nrow(diamond.nr.lineage.23DNCA[diamond.nr.lineage.23DNCA$Domain == "Eukaryota" & diamond.nr.lineage.23DNCA$Kingdom != "Fungi",])

###### Divide contigs based on the classification results
###Loading seq files
library(seqinr)
contig.9D1PG <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D1PG/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D2PG <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D2PG/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D3PH <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D3PH/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D4PH <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D4PH/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D5PG <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D5PG/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")

contig.9D1NH <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D1NH/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D2NG <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D2NG/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D3NG <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D3NG/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D4NG <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D4NG/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.9D5NG <- read.fasta("./4_assembly/megahit_PE/MetaLarge/9D5NG/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")

contig.23D1PA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D1PA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D2PB <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D2PB/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D3PA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D3PA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D4PA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D4PA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D5PB <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D5PB/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")

contig.23D1NA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D1NA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D2NA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D2NA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D3NA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D3NA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D4NA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D4NA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.23D5NA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23D5NA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")

contig.23DNCA <- read.fasta("./4_assembly/megahit_PE/MetaLarge/23DNCA/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")

#### Coassembly
contig.Dominant <- read.fasta("/data/SongyiMetagenome/4_assembly/Coassembly/coassembled_contig_Dominant_rename.fa", as.string = TRUE, seqtype = "DNA")
contig.Minor <- read.fasta("/data/SongyiMetagenome/4_assembly/Coassembly/coassembled_contig_Minor_rename.fa", as.string = TRUE, seqtype = "DNA")

contig.Dominant.2022 <- read.fasta("/data/SongyiMetagenome/4_assembly/Coassembly/Y2022/Dominant/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")
contig.Minor.2022 <- read.fasta("/data/SongyiMetagenome/4_assembly/Coassembly/Y2022/Minor/final.contigs.rename.fa", as.string = TRUE, seqtype = "DNA")


#### Divide and save fasta files
###Define a function
GetProkEukSeqs<-function(diamond.res, contigFastaFile,SampleName){
#### Bacterial and archaeal contigs
prokaryote.contig<-diamond.res$contigName[which(diamond.res$Domain %in% c("Bacteria","Archaea"))]
prokaryote.seq<-contigFastaFile[names(contigFastaFile) %in% prokaryote.contig]
print(length(prokaryote.seq))
#### Fungi and other eukaryotic contigs
eukaryote.contig<-diamond.res$contigName[which(diamond.res$Domain == "Eukaryota")]
eukaryote.seq<-contigFastaFile[names(contigFastaFile) %in% eukaryote.contig]
print(length(eukaryote.seq))

write.fasta(prokaryote.seq, names(prokaryote.seq), paste0("./4_assembly/megahit_PE/MetaLarge/",SampleName,"/final.contigs.rename.prokaryote.fa"))
write.fasta(eukaryote.seq, names(eukaryote.seq), paste0("./4_assembly/megahit_PE/MetaLarge/",SampleName,"/final.contigs.rename.eukaryote.fa"))
}


GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D1PG,contigFastaFile = contig.9D1PG, SampleName = "9D1PG")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D2PG,contigFastaFile = contig.9D2PG, SampleName = "9D2PG")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D3PH,contigFastaFile = contig.9D3PH, SampleName = "9D3PH")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D4PH,contigFastaFile = contig.9D4PH, SampleName = "9D4PH")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D5PG,contigFastaFile = contig.9D5PG, SampleName = "9D5PG")

GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D1NH,contigFastaFile = contig.9D1NH, SampleName = "9D1NH")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D2NG,contigFastaFile = contig.9D2NG, SampleName = "9D2NG")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D3NG,contigFastaFile = contig.9D3NG, SampleName = "9D3NG")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D4NG,contigFastaFile = contig.9D4NG, SampleName = "9D4NG")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.9D5NG,contigFastaFile = contig.9D5NG, SampleName = "9D5NG")

GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D1PA,contigFastaFile = contig.23D1PA, SampleName = "23D1PA")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D2PB,contigFastaFile = contig.23D2PB, SampleName = "23D2PB")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D3PA,contigFastaFile = contig.23D3PA, SampleName = "23D3PA")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D4PA,contigFastaFile = contig.23D4PA, SampleName = "23D4PA")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D5PB,contigFastaFile = contig.23D5PB, SampleName = "23D5PB")

GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D1NA,contigFastaFile = contig.23D1NA, SampleName = "23D1NA")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D2NA,contigFastaFile = contig.23D2NA, SampleName = "23D2NA")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D3NA,contigFastaFile = contig.23D3NA, SampleName = "23D3NA")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D4NA,contigFastaFile = contig.23D4NA, SampleName = "23D4NA")
GetProkEukSeqs(diamond.res = diamond.nr.lineage.23D5NA,contigFastaFile = contig.23D5NA, SampleName = "23D5NA")

GetProkEukSeqs(diamond.res = diamond.nr.lineage.23DNCA,contigFastaFile = contig.23DNCA, SampleName = "23DNCA")


GetBacFunSeqs<-function(diamond.res, contigFastaFile,SampleName){
#### Bacterial and archaeal contigs
prokaryote.contig<-diamond.res$contigName[which(diamond.res$Domain =="Bacteria")]
prokaryote.seq<-contigFastaFile[names(contigFastaFile) %in% prokaryote.contig]
print(length(prokaryote.seq))
#### Fungi and other eukaryotic contigs
eukaryote.contig<-diamond.res$contigName[which(diamond.res$Kingdom == "Fungi")]
eukaryote.seq<-contigFastaFile[names(contigFastaFile) %in% eukaryote.contig]
print(length(eukaryote.seq))

write.fasta(prokaryote.seq, names(prokaryote.seq), paste0("./4_assembly/megahit_PE/MetaLarge/",SampleName,"/final.contigs.rename.bacteria.fa"))
write.fasta(eukaryote.seq, names(eukaryote.seq), paste0("./4_assembly/megahit_PE/MetaLarge/",SampleName,"/final.contigs.rename.fungi.fa"))
}

GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D1PG,contigFastaFile = contig.9D1PG, SampleName = "9D1PG")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D2PG,contigFastaFile = contig.9D2PG, SampleName = "9D2PG")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D3PH,contigFastaFile = contig.9D3PH, SampleName = "9D3PH")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D4PH,contigFastaFile = contig.9D4PH, SampleName = "9D4PH")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D5PG,contigFastaFile = contig.9D5PG, SampleName = "9D5PG")

GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D1NH,contigFastaFile = contig.9D1NH, SampleName = "9D1NH")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D2NG,contigFastaFile = contig.9D2NG, SampleName = "9D2NG")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D3NG,contigFastaFile = contig.9D3NG, SampleName = "9D3NG")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D4NG,contigFastaFile = contig.9D4NG, SampleName = "9D4NG")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.9D5NG,contigFastaFile = contig.9D5NG, SampleName = "9D5NG")

GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D1PA,contigFastaFile = contig.23D1PA, SampleName = "23D1PA")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D2PB,contigFastaFile = contig.23D2PB, SampleName = "23D2PB")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D3PA,contigFastaFile = contig.23D3PA, SampleName = "23D3PA")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D4PA,contigFastaFile = contig.23D4PA, SampleName = "23D4PA")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D5PB,contigFastaFile = contig.23D5PB, SampleName = "23D5PB")

GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D1NA,contigFastaFile = contig.23D1NA, SampleName = "23D1NA")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D2NA,contigFastaFile = contig.23D2NA, SampleName = "23D2NA")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D3NA,contigFastaFile = contig.23D3NA, SampleName = "23D3NA")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D4NA,contigFastaFile = contig.23D4NA, SampleName = "23D4NA")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.23D5NA,contigFastaFile = contig.23D5NA, SampleName = "23D5NA")

GetBacFunSeqs(diamond.res = diamond.nr.lineage.23DNCA,contigFastaFile = contig.23DNCA, SampleName = "23DNCA")



############ For Coassembled samples
#### 2022+2023
GetBacFunSeqs<-function(diamond.res, contigFastaFile,SampleName){
#### Bacterial and archaeal contigs
prokaryote.contig<-diamond.res$contigName[which(diamond.res$superkingdom =="Bacteria")]
prokaryote.seq<-contigFastaFile[names(contigFastaFile) %in% prokaryote.contig]
print(length(prokaryote.seq))
#### Fungi and other eukaryotic contigs
eukaryote.contig<-diamond.res$contigName[which(diamond.res$kingdom == "Fungi")]
eukaryote.seq<-contigFastaFile[names(contigFastaFile) %in% eukaryote.contig]
print(length(eukaryote.seq))

write.fasta(prokaryote.seq, names(prokaryote.seq), paste0("/data/SongyiMetagenome/4_assembly/Coassembly/",SampleName,"/final.contigs.rename.bacteria.fa"))
write.fasta(eukaryote.seq, names(eukaryote.seq), paste0("/data/SongyiMetagenome/4_assembly/Coassembly/",SampleName,"/final.contigs.rename.fungi.fa"))
}

GetBacFunSeqs(diamond.res = diamond.nr.lineage.Dominant,contigFastaFile = contig.Dominant, SampleName = "Dominant")
GetBacFunSeqs(diamond.res = diamond.nr.lineage.Minor,contigFastaFile = contig.Minor, SampleName = "Minor")