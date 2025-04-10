
###### 2. make barplot with phyloseq
### Phyloseqs
bac.clean.ss.f
fun.clean.ss.f

sum(otu_table(fun.clean.ss.f))


order.sample <- c("n6", "n7", "n8", "n9", "n10", "n1","n2","n3","n4","n5",
                  "p6", "p7", "p8", "p9", "p10", "p1","p2","p3","p4","p5",
                  "j6", "j7", "j8", "j9", "j10", "j1","j2","j3","j4","j5",
                  "k6", "k7", "k8", "k9", "k10", "k1","k2","k3","k4","k5",
                  "t6", "t7", "t8", "t9", "t10", "t1","t2","t3","t4","t5")

order.month <- c("March","May","July","September","November")


### Class and phylum level
df.phylum <- bac.clean.ss.dna %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
#df.phylum$Phylum2[which(df.phylum$Class=="Betaproteobacteria")] <- "Betaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
#df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

head(df.phylum)
df.phylum$sampleID <- factor(df.phylum$sampleID, levels = order.sample)

df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))
unique(df.phylum$Phylum2)

levels(df.phylum$Phylum2)
levels(df.phylum$Phylum2) = c(levels(df.phylum$Phylum2), 'Low abundance')

# we need to group by month and inoculation status
df.phylum.rel <- df.phylum %>%  
  group_by(month, inoculation_status) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 0.5,]$Phylum2 <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified", "Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.proteobacteria <- c("Gammaproteobacteria","Alphaproteobacteria")
vec.reorder <- append(vec.uniden.Low, vec.order)
vec.reorder <- append(vec.reorder, vec.proteobacteria)

df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 
df.phylum.rel$month <- factor(df.phylum.rel$month, levels = order.month)
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=month, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~inoculation_status)+
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkolivegreen","Gammaproteobacteria" = "darkolivegreen3",
                               "Verrucomicrobiota" = "#438FFF", "Chloroflexi"= "#ECC846", "Acidobacteriota"= "#AF916D", 
                               "Planctomycetota"= "#8D6FD1", 
                               "Actinobacteriota"="indianred2","WPS-2"= "#71C0A7", "Bacteroidota"="steelblue1",
                               "RCP2-54"= "#B6B0FF", "Firmicutes" ="tan1", "Patescibacteria"= "#3257A8", 
                               "Myxococcota"= "#DD915F","unidentified" = "black", "Low abundance" = "light grey")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1

pdf(file = "230220_Bacteria_Relative abundance_month.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
df.phylum.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()



#### Bar plot of biological replicates
df.phylum.rel <- df.phylum %>%  
  group_by(sampleID, inoculation_status) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance < 0.5,]$Phylum2 <- 'Low abundance'
ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified", "Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.proteobacteria <- c("Gammaproteobacteria","Alphaproteobacteria")
vec.reorder <- append(vec.uniden.Low, vec.order)
vec.reorder <- append(vec.reorder, vec.proteobacteria)

df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 
df.phylum.rel$month <- factor(df.phylum.rel$month, levels = order.month)
df.phylum.rel$sampleID <- factor(df.phylum.rel$sampleID, levels = order.sample)
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=sampleID, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkolivegreen","Gammaproteobacteria" = "darkolivegreen3",
                               "Verrucomicrobiota" = "#438FFF", "Chloroflexi"= "#ECC846", 
                               "Acidobacteriota"= "#AF916D", "Actinobacteriota"="indianred2", "Planctomycetota"= "#8D6FD1",  
                              "WPS-2"= "#71C0A7", "RCP2-54"= "#B6B0FF",
                               "Bacteroidota"="steelblue1","Firmicutes" ="tan1", "Patescibacteria"= "#3257A8", 
                              "Myxococcota"= "#DD915F", "Methylomirabilota"= "#C4B07B", 
                               "Gemmatimonadota"= "#54B5FB", "Cyanobacteria"= "#57B956",
                               "Armatimonadota"= "#F8BCBD", "Dependentiae" = "#F2D7EE","GAL15"= "#D8D7BF", 
                               "unidentified" = "black", "Low abundance" = "light grey")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.phylum.rel.p1

pdf(file = "230220_Bacteria_Relative abundance_each sample_phylum.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
df.phylum.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()


#### Genus level
df.genus <- bac.clean.ss.dna %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

head(df.genus)
df.genus$sampleID <- factor(df.genus$sampleID, levels = order.sample)

df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))
unique(df.genus$Genus)

levels(df.genus$Genus)
levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')

# we need to group by month and inoculation status
df.genus.rel <- df.genus %>%  
  group_by(month, inoculation_status) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance < 0.2,]$Genus <- 'Low abundance'
ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified", "Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.reorder <- append(vec.uniden.Low, vec.order)

df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 
df.genus.rel$month <- factor(df.genus.rel$month, levels = order.month)
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=month, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~inoculation_status)+
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_otu2) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 9,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.genus.rel.p1

pdf(file = "230220_Bacteria_Relative abundance_genus_month.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 9) # The height of the plot in inches

# Step 2: Create the plot with R code
df.genus.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()


### Genus level- replicates
df.genus.rel <- df.genus %>%  
  group_by(sampleID) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance < 0.2,]$Genus <- 'Low abundance'
ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac%in%c("Low abundance", "unidentified", "Alphaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria"))]
vec.uniden.Low <- c("Low abundance", "unidentified")
vec.reorder <- append(vec.uniden.Low, vec.order)

df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 
df.genus.rel$month <- factor(df.genus.rel$month, levels = order.month)
df.genus.rel$sampleID <- factor(df.genus.rel$sampleID, levels = order.sample)
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=sampleID, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack')+
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_otu2) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 15,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.genus.rel.p1

pdf(file = "230220_Bacteria_Relative abundance_genus_replicates.pdf",   # The directory you want to save the file in
    width = 25, # The width of the plot in inches
    height = 12) # The height of the plot in inches

# Step 2: Create the plot with R code
df.genus.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()


#### Fungi #####
## Phylum level
##Grouped by sample
df.phylum.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE) %>% 
  psmelt()

df.phylum.fun %<>% mutate(Phylum = fct_explicit_na(Phylum, na_level = "unidentified"))
levels(df.phylum.fun$Phylum) = c(levels(df.phylum.fun$Phylum), 'Low abundance')

# we need to group by samples
df.phylum.fun.rel <- df.phylum.fun %>%  
  group_by(month, inoculation_status) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.phylum.fun.rel[df.phylum.fun.rel$RelAbundance < 0.1,]$Phylum <- 'Low abundance'
unique(df.phylum.fun$Phylum)

ord.f <- df.phylum.fun.rel %>% group_by(Phylum) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Phylum
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.phylum.fun.rel$Phylum <- factor(df.phylum.fun.rel$Phylum, levels = vec.reorder.f) 
df.phylum.fun.rel$sampleID <-factor(df.phylum.fun.rel$sampleID, levels = order.sample)
df.phylum.fun.rel$month <-factor(df.phylum.fun.rel$month, levels = order.month)
## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.fun.rel.p1 <- ggplot(df.phylum.fun.rel, aes(x=month, y = RelAbundance, fill = Phylum)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~inoculation_status)+
  #scale_fill_discrete() +
  scale_fill_manual(values = c( "Basidiomycota"= "#BE4146","Ascomycota" = "#11335F",
                                "Mortierellomycota"= "#F3C911", "Rozellomycota"= "#F8BCBD",
                                "Mucoromycota"= "#478F48", "Basidiobolomycota"= "#8B3D88", "Glomeromycota" = "#652926",
                                "Calcarisporiellomycota"= "#54B5FB",  "Kickxellomycota" = "#E5D1D0",
                               "unidentified" ="#000000",
                               "Low abundance"="#BFBEBE")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
df.phylum.fun.rel.p1


pdf(file = "230220_Fungi_Relative abundance_Phylum_month.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
df.phylum.fun.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()


###class
df.class.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.rel <- df.class.fun %>%  
  group_by(month, inoculation_status) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.class.fun.rel[df.class.fun.rel$RelAbundance < 0.5,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.class.fun.rel$Class <- factor(df.class.fun.rel$Class, levels = vec.reorder.f) 
df.class.fun.rel$sampleID <- factor(df.class.fun.rel$sampleID, levels = order.sample) 
df.class.fun.rel$month <- factor(df.class.fun.rel$month, levels = order.month) 


## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.class.fun.rel.p1 <- ggplot(df.class.fun.rel, aes(x=month, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~inoculation_status)+
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Agaricomycetes" = "#CC6C71","Leotiomycetes"= "#2878BD",
                               "Geminibasidiomycetes"="#478F48" , "Mortierellomycetes"="#F3C911",  "Eurotiomycetes"= "#6DA9DC",
                               "Dothideomycetes" =  "#99cc99",
                               "Pezizomycetes"="#ff985E","Sordariomycetes"= "#1E63AF",
                               "Rozellomycotina_cls_Incertae_sedis" = "#F8BCBD", "Rozellomycota_cls_Incertae_sedis" = "#D3D0CB",
                               "Umbelopsidomycetes" ="#AF916D", "Archaeorhizomycetes" = "#A871AE",
                               "Tritirachiomycetes" = "#FFF275", "GS25" = "#FA7921", 
                               "unidentified" ="#000000", "Low abundance" = "light grey")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse =T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
df.class.fun.rel.p1


pdf(file = "230220_Fungi_Relative abundance_Class_month.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
df.class.fun.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()

#### Class - replicates
df.class.fun.rel <- df.class.fun %>%  
  group_by(sampleID) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.class.fun.rel[df.class.fun.rel$RelAbundance < 0.5,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.class.fun.rel$Class <- factor(df.class.fun.rel$Class, levels = vec.reorder.f) 
df.class.fun.rel$sampleID <- factor(df.class.fun.rel$sampleID, levels = order.sample) 
df.class.fun.rel$month <- factor(df.class.fun.rel$month, levels = order.month) 


## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.class.fun.rel.p1 <- ggplot(df.class.fun.rel, aes(x=sampleID, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Agaricomycetes" = "#CC6C71","Leotiomycetes"= "#2878BD",
                               "Geminibasidiomycetes"="#478F48" , "Mortierellomycetes"="#F3C911",  
                               "Eurotiomycetes"= "#6DA9DC","Dothideomycetes" =  "#99cc99","Sordariomycetes"= "#1E63AF",
                               "Pezizomycetes"="#ff985E",
                               "Rozellomycota_cls_Incertae_sedis" = "#D3D0CB",
                               "Rozellomycotina_cls_Incertae_sedis" = "#F8BCBD", 
                               "Umbelopsidomycetes" ="#AF916D",
                               "Mucoromycotina_cls_Incertae_sedis"="#AE3535",
                               "Tremellomycetes"="#B6B0FF",
                               "Tritirachiomycetes"="#4A8DDC",                             
                               "Archaeorhizomycetes" = "#A871AE",
                               "Basidiobolomycetes"="#C39B6A", 
                               "GS25" = "#FA7921", 
                               "Microbotryomycetes" ="#D8D7BF",
                               "Calcarisporiellomycetes"="#FF977E",
                               "Ramicandelaberomycetes"="#B66DFF",
                               "unidentified" ="#000000", "Low abundance" = "light grey")) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = F))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
df.class.fun.rel.p1


pdf(file = "230220_Fungi_Relative abundance_Class_replicates.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Step 2: Create the plot with R code
df.class.fun.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()


#### Genus level - month
df.genus.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

df.genus.fun %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))
levels(df.genus.fun$Genus) = c(levels(df.genus.fun$Genus), 'Low abundance')

# we need to group by samples
df.genus.fun.rel <- df.genus.fun %>%  
  group_by(month, inoculation_status) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.fun.rel[df.genus.fun.rel$RelAbundance < 0.2,]$Genus <- 'Low abundance'
unique(df.genus.fun$Genus)

ord.f <- df.genus.fun.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Genus
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.genus.fun.rel$Genus <- factor(df.genus.fun.rel$Genus, levels = vec.reorder.f) 
df.genus.fun.rel$sampleID <- factor(df.genus.fun.rel$sampleID, levels = order.sample) 
df.genus.fun.rel$month <- factor(df.genus.fun.rel$month, levels = order.month) 


## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.fun.rel.p1 <- ggplot(df.genus.fun.rel, aes(x=month, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + facet_wrap(~inoculation_status)+
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_otu2) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 6,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
df.genus.fun.rel.p1


pdf(file = "230220_Fungi_Relative abundance_Genus_month.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 9) # The height of the plot in inches

# Step 2: Create the plot with R code
df.genus.fun.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()


#### Genus - replicates
df.genus.fun.rel <- df.genus.fun %>%  
  group_by(sampleID) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*100/sum(Abundance))  # Transform to rel. abundance

df.genus.fun.rel[df.genus.fun.rel$RelAbundance < 0.2,]$Genus <- 'Low abundance'
unique(df.genus.fun$Genus)

ord.f <- df.genus.fun.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Genus
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.genus.fun.rel$Genus <- factor(df.genus.fun.rel$Genus, levels = vec.reorder.f) 
df.genus.fun.rel$sampleID <- factor(df.genus.fun.rel$sampleID, levels = order.sample) 
df.genus.fun.rel$month <- factor(df.genus.fun.rel$month, levels = order.month) 


## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.genus.fun.rel.p1 <- ggplot(df.genus.fun.rel, aes(x=sampleID, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') + 
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_otu2) +
  
  xlab('')+
  ylab("Relative abundance (%)\n") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 6,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,100,20))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())
df.genus.fun.rel.p1


pdf(file = "230220_Fungi_Relative abundance_Genus_replicate.pdf",   # The directory you want to save the file in
    width = 25, # The width of the plot in inches
    height = 12) # The height of the plot in inches

# Step 2: Create the plot with R code
df.genus.fun.rel.p1

# Step 3: Run dev.off() to create the file!
dev.off()


##### Relative abundance table
bac.clean.ss.rel <- transform(bac.clean.ss.dna, transform = "compositional")
otu.bac.rel<-otu_table(bac.clean.ss.rel)
otu.bac.rel <- data.frame(otu.bac.rel)
otu.bac.rel$Total <- rowSums(otu.bac.rel)
otu.bac.rel$OTU <- rownames(otu.bac.rel)

bac.list.sub <- bac.list %>% select(OTU, OTU_id, number, Phylum,Class,Order,Family,Genus)

otu.bac.rel.tab<- merge(otu.bac.rel,bac.list.sub, by ="OTU")
write.csv(otu.bac.rel.tab,"230220_Distribution of bacterial ASVs_dna samples.csv")

fun.clean.ss.rel <- transform(fun.clean.ss, transform = "compositional")
otu.fun.rel<-otu_table(fun.clean.ss.rel)
otu.fun.rel <- data.frame(otu.fun.rel)
otu.fun.rel$Total <- rowSums(otu.fun.rel)
otu.fun.rel$OTU <- rownames(otu.fun.rel)

fun.list.sub <- fun.list %>% select(OTU, OTU_id, number, Phylum,Class,Order,Family,Genus)

otu.fun.rel.tab<- merge(otu.fun.rel,fun.list.sub, by ="OTU")
write.csv(otu.fun.rel.tab,"230220_Distribution of fungal ASVs.csv")

#### Relative abundance table
melt.bac.rel <- psmelt(bac.clean.ss.rel)
melt.fun.rel <- psmelt(fun.clean.ss.rel)

melt.bac.rel$month <- factor(melt.bac.rel$month, levels = order.month)
melt.bac.rel$sampleID <- factor(melt.bac.rel$sampleID, levels = order.sample)
melt.fun.rel$month <- factor(melt.fun.rel$month, levels = order.month)
melt.fun.rel$sampleID <- factor(melt.fun.rel$sampleID, levels = order.sample)

###Bacteria
bac.rel.phylum <- melt.bac.rel %>% group_by(sampleID,month,inoculation_status,Phylum) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.phylum,"230220_bacterial phyla and class_relative abundance_replicates.csv")

bac.rel.class <- melt.bac.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.class,"230220_bacterial class_relative abundance_replicates.csv")

bac.rel.order <- melt.bac.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class,Order) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.order,"230220_bacterial order_relative abundance_replicates.csv")

bac.rel.family <- melt.bac.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class,Order,Family) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.family,"230220_bacterial family_relative abundance_replicates.csv")

bac.rel.genus <- melt.bac.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class,Order,Family,Genus) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(bac.rel.genus,"230220_bacterial genus_relative abundance_replicates.csv")


### Abundance table for month data (mean and standard deviation)
bac.rel.phylum.month<-bac.rel.phylum %>% group_by(month, inoculation_status, Phylum)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.phylum.month,"230220_bacterial phyla_relative abundance_month.csv")

bac.rel.class.month<-bac.rel.class %>% group_by(month, inoculation_status, Phylum, Class)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.class.month,"230220_bacterial class_relative abundance_month.csv")

bac.rel.order.month<-bac.rel.order %>% group_by(month, inoculation_status, Phylum, Class, Order)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.order.month,"230220_bacterial order_relative abundance_month.csv")

bac.rel.family.month<-bac.rel.family %>% group_by(month, inoculation_status, Phylum, Class, Order, Family)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.family.month,"230220_bacterial family_relative abundance_month.csv")

bac.rel.genus.month<-bac.rel.genus %>% group_by(month, inoculation_status, Phylum, Class, Order, Family, Genus)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(bac.rel.genus.month,"230220_bacterial genus_relative abundance_month.csv")

### fungi
fun.rel.phylum <- melt.fun.rel %>% group_by(sampleID,month,inoculation_status,Phylum) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.phylum,"230220_fungal phyla and class_relative abundance_replicates.csv")

fun.rel.class <- melt.fun.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.class,"230220_fungal class_relative abundance_replicates.csv")

fun.rel.order <- melt.fun.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class,Order) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.order,"230220_fungal order_relative abundance_replicates.csv")

fun.rel.family <- melt.fun.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class,Order,Family) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.family,"230220_fungal family_relative abundance_replicates.csv")

fun.rel.genus <- melt.fun.rel %>% group_by(sampleID,month,inoculation_status,Phylum,Class,Order,Family,Genus) %>% summarise(RelAbundance = sum(Abundance), Percent = RelAbundance*100)
write.csv(fun.rel.genus,"230220_fungal genus_relative abundance_replicates.csv")


### Abundance table for month data (mean and standard deviation)
fun.rel.phylum.month<-fun.rel.phylum %>% group_by(month, inoculation_status, Phylum)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.phylum.month,"230220_fungal phyla_relative abundance_month.csv")

fun.rel.class.month<-fun.rel.class %>% group_by(month, inoculation_status, Phylum, Class)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.class.month,"230220_fungal class_relative abundance_month.csv")

fun.rel.order.month<-fun.rel.order %>% group_by(month, inoculation_status, Phylum, Class, Order)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.order.month,"230220_fungal order_relative abundance_month.csv")

fun.rel.family.month<-fun.rel.family %>% group_by(month, inoculation_status, Phylum, Class, Order, Family)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.family.month,"230220_fungal family_relative abundance_month.csv")

fun.rel.genus.month<-fun.rel.genus %>% group_by(month, inoculation_status, Phylum, Class, Order, Family, Genus)%>% summarise(MeanRelAbundance = mean(RelAbundance), SdRelAbund = sd(RelAbundance), MeanPercent = mean(RelAbundance)*100, SdPercent=SdRelAbund*100)
write.csv(fun.rel.genus.month,"230220_fungal genus_relative abundance_month.csv")

