### Alpha diversity
bac.rarefied
fun.rarefied

### Metadata
b.meta <- data.frame(sample_data(bac.clean.ss.dna.f))
f.meta <- data.frame(sample_data(fun.clean.ss.f))

### Estimate alpha diversity indices
library(microbiome)

## Bacteria
all.index.bac<-microbiome::alpha(bac.rarefied, index = "all")
head(all.index.bac)
colnames(all.index.bac)
# [1] "observed"                   "chao1"                      "diversity_inverse_simpson"  "diversity_gini_simpson"    
# [5] "diversity_shannon"          "diversity_fisher"           "diversity_coverage"         "evenness_camargo"          
# [9] "evenness_pielou"            "evenness_simpson"           "evenness_evar"              "evenness_bulla"            
# [13] "dominance_dbp"              "dominance_dmn"              "dominance_absolute"         "dominance_relative"        
# [17] "dominance_simpson"          "dominance_core_abundance"   "dominance_gini"             "rarity_log_modulo_skewness"
# [21] "rarity_low_abundance"       "rarity_rare_abundance" 

##Fungi
all.index.fun<-microbiome::alpha(fun.rarefied, index = "all")
head(all.index.fun)


## Add alpha diversity values into the metadata
b.meta$Chao1 <- all.index.bac$chao1
b.meta$Shannon <- all.index.bac$diversity_shannon
b.meta$Simpson <- all.index.bac$dominance_simpson
#b.meta$Rarity<-all.index.bac$rarity_log_modulo_skewness

f.meta$Chao1 <- all.index.fun$chao1
f.meta$Shannon <- all.index.fun$diversity_shannon
f.meta$Simpson <- all.index.fun$dominance_simpson
#f.meta$Rarity<-all.index.fun$rarity_log_modulo_skewness


### Dot plot (the alternative of a box plot)
## Add a data summary
b.meta$Status <- ifelse(b.meta$Status == "TmD","TmD","TmM")
f.meta$Status <- ifelse(f.meta$Status == "TmD","TmD","TmM")


b.meta.data.summary <- b.meta %>% group_by(Status,month) %>% mutate(meanChao1 = mean(Chao1), sdChao1 = sd(Chao1), 
                                                                                    meanShannon = mean(Shannon), sdShannon = sd(Shannon),
                                                                                    meanSimpson = mean(Simpson), sdSimpson = sd(Simpson))

f.meta.data.summary <- f.meta %>% group_by(Status,month) %>% mutate(meanChao1 = mean(Chao1), sdChao1 = sd(Chao1), 
                                                                                 meanShannon = mean(Shannon), sdShannon = sd(Shannon),
                                                                                 meanSimpson = mean(Simpson), sdSimpson = sd(Simpson))
color_status <- c("TmM" = "khaki3","TmD" = "dodgerblue4")
b.meta.data.summary$month <- factor(b.meta.data.summary$month, levels =c("March","May","July","September","November"))
f.meta.data.summary$month <- factor(f.meta.data.summary$month, levels =c("March","May","July","September","November"))

### Bacteria - Chao1
dp <- ggplot(b.meta.data.summary, aes(x=month, y=Chao1, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.8))+
  scale_fill_manual(values = color_status)+
  # stat_summary(fun.data=mean_sdl, mult=1, 
  #              geom="pointrange", color="red")+
  theme(legend.position = "top")+
  geom_errorbar(aes(ymin=meanChao1-sdChao1, ymax=meanChao1+sdChao1), width=.2,
                position=position_dodge(.8))+
  geom_point(aes(y=meanChao1), size = 3,position=position_dodge(.8))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.1,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
dp

### Bacteria - Shannon
dp <- ggplot(b.meta.data.summary, aes(x=month, y=Shannon, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.8))+
  scale_fill_manual(values = color_status)+
  # stat_summary(fun.data=mean_sdl, mult=1, 
  #              geom="pointrange", color="red")+
  theme(legend.position = "top")+
  geom_errorbar(aes(ymin=meanShannon-sdShannon, ymax=meanShannon+sdShannon), width=.2,
                position=position_dodge(.8))+
  geom_point(aes(y=meanShannon), size = 3,position=position_dodge(.8))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.1,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
dp


### Fungi - Chao1
dp <- ggplot(f.meta.data.summary, aes(x=Condition, y=Chao1, fill=Treatment)) +
  geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.8))+
  scale_fill_manual(values = color_treatment)+
  # stat_summary(fun.data=mean_sdl, mult=1, 
  #              geom="pointrange", color="red")+
  theme(legend.position = "top")+
  geom_errorbar(aes(ymin=meanChao1-sdChao1, ymax=meanChao1+sdChao1), width=.2,
                position=position_dodge(.8))+
  geom_point(aes(y=meanChao1), size = 3,position=position_dodge(.8))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.1,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
dp

### Fungi - Shannon
dp <- ggplot(f.meta.data.summary, aes(x=month, y=Shannon, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.8))+
  scale_fill_manual(values = color_status)+
  # stat_summary(fun.data=mean_sdl, mult=1, 
  #              geom="pointrange", color="red")+
  theme(legend.position = "top")+
  geom_errorbar(aes(ymin=meanShannon-sdShannon, ymax=meanShannon+sdShannon), width=.2,
                position=position_dodge(.8))+
  geom_point(aes(y=meanShannon), size = 3,position=position_dodge(.8))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.1,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
dp



dev.off()


### Save raw data tables
write.csv(b.meta.data.summary, "~/Desktop/SNU/SNU_BaeInHyup/Manuscript/Comprehensive/241112_Metadata with alpha dversity indices_bacteria.csv")
write.csv(f.meta.data.summary, "~/Desktop/SNU/SNU_BaeInHyup/Manuscript/Comprehensive/241112_Metadata with alpha dversity indices_fungi.csv")

### Statistical analysis
## Bacteria - Shannon
## Make a data frame
stat.tab.bac.Shannon<-data.frame(Month = c("March", "May", "July", "September","November"), tValue = c(0,0,0,0,0), df = c(0,0,0,0,0), Pvalue = c(0,0,0,0,0),
                               Index = c("Shannon","Shannon","Shannon","Shannon","Shannon"))
## Normality test
for (i in c("March", "May", "July", "September","November")){
  normalityCheck<-shapiro.test(b.meta.data.summary$Shannon[which(b.meta.data.summary$month == i)]) 
  hist(b.meta.data.summary$Shannon[which(b.meta.data.summary$month == i)])
  print(normalityCheck)
  
}

## March W = 0.86753, p-value = 0.09353 # follow normal distribution
## May W = 0.90052, p-value = 0.222# follow normal distribution
## July W = 0.86782, p-value = 0.09429 # follow normal distribution
## Sep W = 0.86782, p-value = 0.09429 # follow normal distribution

## Variance
for (i in c("March", "May", "July", "September","November")){
  var.TmM<-var(b.meta.data.summary$Shannon[which(b.meta.data.summary$Status == "TmM" & b.meta.data.summary$month == i)]) 
  var.TmD<-var(b.meta.data.summary$Shannon[which(b.meta.data.summary$Status == "TmD" & b.meta.data.summary$month == i)]) 
  print(var.TmM)
  print(var.TmD)
}


## Welch t-test
for (i in c(c("March", "May", "July", "September","November"))){
  statResult<-t.test(b.meta.data.summary$Shannon[which(b.meta.data.summary$Status == "TmM" & b.meta.data.summary$month == i)],
                     b.meta.data.summary$Shannon[which(b.meta.data.summary$Status == "TmD" & b.meta.data.summary$month == i)],
                     var.equal= F)
  
  stat.tab.bac.Shannon$tValue[which(stat.tab.bac.Shannon$Month == i)]<-statResult$statistic
  stat.tab.bac.Shannon$Pvalue[which(stat.tab.bac.Shannon$Month == i)]<-statResult$p.value
  stat.tab.bac.Shannon$df[which(stat.tab.bac.Shannon$Month == i)]<-statResult$parameter
}


## Fungi - Shannon
## Make a data frame
stat.tab.fun.Shannon<-data.frame(Month = c("March", "May", "July", "September","November"), tValue = c(0,0,0,0,0), df = c(0,0,0,0,0), Pvalue = c(0,0,0,0,0),
                                 Index = c("Shannon","Shannon","Shannon","Shannon","Shannon"))
## Normality test
for (i in c("March", "May", "July", "September","November")){
  normalityCheck<-shapiro.test(f.meta.data.summary$Shannon[which(f.meta.data.summary$month == i)]) 
  hist(f.meta.data.summary$Shannon[which(f.meta.data.summary$month == i)])
  print(normalityCheck)
  
}

## March W = 0.91555, p-value = 0.3213 # follow normal distribution
## May W = 0.85669, p-value = 0.06975# follow normal distribution
## July W = 0.87132, p-value = 0.1036 # follow normal distribution
## Sep W = 0.91922, p-value = 0.3505 # follow normal distribution



## Variance
for (i in c("March", "May", "July", "September","November")){
  var.TmM<-var(f.meta.data.summary$Shannon[which(f.meta.data.summary$Status == "TmM" & f.meta.data.summary$month == i)]) 
  var.TmD<-var(f.meta.data.summary$Shannon[which(f.meta.data.summary$Status == "TmD" & f.meta.data.summary$month == i)]) 
  print(var.TmM)
  print(var.TmD)
}


## Welch t-test
for (i in c(c("March", "May", "July", "September","November"))){
  statResult<-t.test(f.meta.data.summary$Shannon[which(f.meta.data.summary$Status == "TmM" & f.meta.data.summary$month == i)],
                     f.meta.data.summary$Shannon[which(f.meta.data.summary$Status == "TmD" & f.meta.data.summary$month == i)],
                     var.equal= F)
  
  stat.tab.fun.Shannon$tValue[which(stat.tab.fun.Shannon$Month == i)]<-statResult$statistic
  stat.tab.fun.Shannon$Pvalue[which(stat.tab.fun.Shannon$Month == i)]<-statResult$p.value
  stat.tab.fun.Shannon$df[which(stat.tab.fun.Shannon$Month == i)]<-statResult$parameter
}


write.csv(stat.tab.fun.Shannon, "~/Desktop/SNU/SNU_BaeInHyup/Manuscript/Comprehensive/241112_Alpha diversity_stat_fungi.csv")
write.csv(stat.tab.bac.Shannon, "~/Desktop/SNU/SNU_BaeInHyup/Manuscript/Comprehensive/241112_Alpha diversity_stat_bacteria.csv")
