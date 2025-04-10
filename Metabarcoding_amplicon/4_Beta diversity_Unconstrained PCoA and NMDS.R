### For bacteria ###
##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
##Unconstrained PCoA
color_inocul <- c("no" = "khaki3","yes" = "dodgerblue4")

color_month <- c("March" = "darkslategray4",
                 "May" = "steelblue3","July" = "salmon3","September" = "lavenderblush4", "November" = "#652926")



bray1.bac <-  ordinate(bac.clean.log.dna, 'PCoA', 'bray')

# write.csv(bray1.bac$vectors, "Bacteria_PCoA_2017.csv")
# write.csv(bray2.bac$vectors, "Bacteria_PCoA_2018.csv")


pcoa.bac.inocul<-plot_ordination(bac.clean.log.dna, bray1.bac, type = "samples", color='inoculation_status', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_inocul)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
  #+stat_ellipse(type = "t", linetype = 2, level = 0.95)

pcoa.bac.month<-plot_ordination(bac.clean.log.dna, bray1.bac, type = "samples", color='month', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_month)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
#+stat_ellipse(type = "t", linetype = 2, level = 0.95)

pdf(file = "230220_bacteria_PCoA_inoculation.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
pcoa.bac.inocul

# Step 3: Run dev.off() to create the file!
dev.off()


pdf(file = "230220_bacteria_PCoA_month.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
pcoa.bac.month

# Step 3: Run dev.off() to create the file!
dev.off()


### Fungi
bray1.fun <-  ordinate(fun.clean.log, 'PCoA', 'bray')

# write.csv(bray1.fun$vectors, "Bacteria_PCoA_2017.csv")
# write.csv(bray2.bac$vectors, "Bacteria_PCoA_2018.csv")


pcoa.fun.inocul<-plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='inoculation_status', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_inocul)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
#+stat_ellipse(type = "t", linetype = 2, level = 0.95)

pcoa.fun.month<-plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='month', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_color_manual(values=color_month)+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
#+stat_ellipse(type = "t", linetype = 2, level = 0.95)

pdf(file = "230220_fungi_PCoA_inoculation.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
pcoa.fun.inocul

# Step 3: Run dev.off() to create the file!
dev.off()


pdf(file = "230220_fungi_PCoA_month.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
pcoa.fun.month

# Step 3: Run dev.off() to create the file!
dev.off()


##PERMANOVA
b.otu.all <- otu_table(bac.clean.nolog.dna)
b.meta.all <- data.frame(sample_data(bac.clean.log.dna))

b.permanova.all <- adonis2(formula = t(b.otu.all) ~ (inoculation_status+month), data = b.meta.all, permutations=9999, method = "bray")
b.permanova.all

# adonis2(formula = t(b.otu.all) ~ (inoculation_status + month), data = b.meta.all, permutations = 9999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# inoculation_status  1   2.7282 0.16255 9.6561 0.0001 ***
#   month               4   1.6243 0.09677 1.4372 0.0085 ** 
#   Residual           44  12.4317 0.74068                  
# Total              49  16.7841 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


f.otu.all <- otu_table(fun.clean.nolog)
f.meta.all <- data.frame(sample_data(fun.clean.log))

f.permanova.all <- adonis2(formula = t(f.otu.all) ~ (inoculation_status+month), data = f.meta.all, permutations=9999, method = "bray")
f.permanova.all

# adonis2(formula = t(f.otu.all) ~ (inoculation_status + month), data = f.meta.all, permutations = 9999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# inoculation_status  1   1.7845 0.17377 12.383  1e-04 ***
#   month               4   2.1438 0.20876  3.719  1e-04 ***
#   Residual           44   6.3410 0.61747                  
# Total              49  10.2694 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



### PERMANOVA with soil chemical properties
sampleList<-c("j1","j2","j3", "j8","j9","j10", "k1","k2","k3", "k8","k9","k10", "t1","t2","t3", "t8","t9","t10")

bac.clean.nolog.chemicals<-subset_samples(bac.clean.nolog.dna, sampleID %in% sampleList)
bac.clean.nolog.chemicals <- phyloseq::filter_taxa(bac.clean.nolog.chemicals, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.chemicals<-subset_samples(fun.clean.nolog, sampleID %in% sampleList)
fun.clean.nolog.chemicals <- phyloseq::filter_taxa(fun.clean.nolog.chemicals, function(x) sum(x) != 0, TRUE)

bac.clean.log.chemicals<-subset_samples(bac.clean.log.dna, sampleID %in% sampleList)
bac.clean.log.chemicals <- phyloseq::filter_taxa(bac.clean.log.chemicals, function(x) sum(x) != 0, TRUE)

fun.clean.log.chemicals<-subset_samples(fun.clean.log, sampleID %in% sampleList)
fun.clean.log.chemicals <- phyloseq::filter_taxa(fun.clean.log.chemicals, function(x) sum(x) != 0, TRUE)



b.otu.soil <- otu_table(bac.clean.log.chemicals)
b.meta.soil <- data.frame(sample_data(bac.clean.log.chemicals))

b.permanova.soil <- adonis2(formula = t(b.otu.soil) ~ (pH+SOM+TN+TP+P2O5+Na+K+Ca+Mg+Cu+Zn+Fe+Mn+Water), data = b.meta.soil, permutations=9999, method = "bray")
b.permanova.soil



f.otu.soil <- otu_table(fun.clean.log.chemicals)
f.meta.soil <- data.frame(sample_data(fun.clean.log.chemicals))

f.permanova.soil <- adonis2(formula = t(f.otu.soil) ~ (pH+SOM+TN+TP+P2O5+Na+K+Ca+Mg+Cu+Zn+Fe+Mn+Water), data = f.meta.soil, permutations=9999, method = "bray")
f.permanova.soil
