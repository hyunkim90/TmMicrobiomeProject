#### cor0.2
get_gephi_input_withp('IAA',2,0.2)
IAA.node2 <- read.table("./2_0.2_node_IAA.tsv",sep="\t",header = T)
nrow(IAA.node2) #3657 genes (out of 4,048)
IAA.edge2 <- read.table("./2_0.2_edge_IAA.tsv",sep="\t",header = T)
head(IAA.edge)
nrow(IAA.edge2) #64,516


#### Add additional node information from eggnog results
names(IAA.node2)[2] <- "NewID"
IAA.node2.detail.info<-merge(IAA.node2,rep.gene.info.TPM, by = "NewID")
nrow(IAA.node2.detail.info)
nrow(IAA.node2)
head(IAA.node2.detail.info)
#### Remove unnecessary information
names(IAA.node2.detail.info)
IAA.node2.detail.info.essence <- IAA.node2.detail.info[c(1,2,3,4,13,14,15,16,29:36)]
head(IAA.node2.detail.info.essence)

#### Refine taxonomy information (treat unidentified)
IAA.node2.detail.info.essence$Phylum[which(IAA.node2.detail.info.essence$Domain=="Bacteria" & IAA.node2.detail.info.essence$Phylum == "Unidentified")] <- "Bacteria_unknown"
IAA.node2.detail.info.essence$Phylum[which(IAA.node2.detail.info.essence$Kingdom=="Fungi" & IAA.node2.detail.info.essence$Phylum == "Unidentified")] <- "Fungi_unknown"
unique(IAA.node2.detail.info.essence$Phylum)

IAA.node2.detail.info.essence$Class[grepl("_unknown",IAA.node2.detail.info.essence$Phylum)&IAA.node2.detail.info.essence$Class == "Unidentified"] <- IAA.node2.detail.info.essence$Phylum[grepl("_unknown",IAA.node2.detail.info.essence$Phylum)&IAA.node2.detail.info.essence$Class == "Unidentified"]
IAA.node2.detail.info.essence$Class[!(grepl("_unknown",IAA.node2.detail.info.essence$Phylum))&IAA.node2.detail.info.essence$Class == "Unidentified"] <- paste0(IAA.node2.detail.info.essence$Phylum[!(grepl("_unknown",IAA.node2.detail.info.essence$Phylum))&IAA.node2.detail.info.essence$Class == "Unidentified"],
                                                                                                                                                            "_unknown")

unique(IAA.node2.detail.info.essence$Class)

IAA.node2.detail.info.essence$Order[grepl("_unknown",IAA.node2.detail.info.essence$Class)&IAA.node2.detail.info.essence$Order == "Unidentified"] <- IAA.node2.detail.info.essence$Class[grepl("_unknown",IAA.node2.detail.info.essence$Class)&IAA.node2.detail.info.essence$Order == "Unidentified"]
IAA.node2.detail.info.essence$Order[!(grepl("_unknown",IAA.node2.detail.info.essence$Class))&IAA.node2.detail.info.essence$Order == "Unidentified"] <- paste0(IAA.node2.detail.info.essence$Class[!(grepl("_unknown",IAA.node2.detail.info.essence$Class))&IAA.node2.detail.info.essence$Order == "Unidentified"],
                                                                                                                                                           "_unknown")

unique(IAA.node2.detail.info.essence$Order)


IAA.node2.detail.info.essence$Family[grepl("_unknown",IAA.node2.detail.info.essence$Order)&IAA.node2.detail.info.essence$Family == "Unidentified"] <- IAA.node2.detail.info.essence$Order[grepl("_unknown",IAA.node2.detail.info.essence$Order)&IAA.node2.detail.info.essence$Family == "Unidentified"]
IAA.node2.detail.info.essence$Family[!(grepl("_unknown",IAA.node2.detail.info.essence$Order))&IAA.node2.detail.info.essence$Family == "Unidentified"] <- paste0(IAA.node2.detail.info.essence$Order[!(grepl("_unknown",IAA.node2.detail.info.essence$Order))&IAA.node2.detail.info.essence$Family == "Unidentified"],
                                                                                                                                                             "_unknown")

unique(IAA.node2.detail.info.essence$Family)

IAA.node2.detail.info.essence$Genus[grepl("_unknown",IAA.node2.detail.info.essence$Family)&IAA.node2.detail.info.essence$Genus == "Unidentified"] <- IAA.node2.detail.info.essence$Family[grepl("_unknown",IAA.node2.detail.info.essence$Family)&IAA.node2.detail.info.essence$Genus == "Unidentified"]
IAA.node2.detail.info.essence$Genus[!(grepl("_unknown",IAA.node2.detail.info.essence$Family))&IAA.node2.detail.info.essence$Genus == "Unidentified"] <- paste0(IAA.node2.detail.info.essence$Family[!(grepl("_unknown",IAA.node2.detail.info.essence$Family))&IAA.node2.detail.info.essence$Genus == "Unidentified"],
                                                                                                                                                            "_unknown")

unique(IAA.node2.detail.info.essence$Genus)

IAA.node2.detail.info.essence$Species[grepl("_unknown",IAA.node2.detail.info.essence$Genus)&IAA.node2.detail.info.essence$Species == "Unidentified"] <- IAA.node2.detail.info.essence$Genus[grepl("_unknown",IAA.node2.detail.info.essence$Genus)&IAA.node2.detail.info.essence$Species == "Unidentified"]
IAA.node2.detail.info.essence$Species[!(grepl("_unknown",IAA.node2.detail.info.essence$Genus))&IAA.node2.detail.info.essence$Species == "Unidentified"] <- paste0(IAA.node2.detail.info.essence$Genus[!(grepl("_unknown",IAA.node2.detail.info.essence$Genus))&IAA.node2.detail.info.essence$Species == "Unidentified"],
                                                                                                                                                               "_unknown")

unique(IAA.node2.detail.info.essence$Species)


unique(NH3transport.genes$PFAMs)
#### Add functional category information for each gene
target.genes<-unique(c(nitrate.reduc.genes$GeneID,NH3transport.genes$GeneID,glutamate.produc.genes$GeneID,Tryptophan.produc.genes$GeneID,
                       IAA.produc.genes$GeneID)) #4,048 genes
IAA.node2.detail.info.essence$FuncCategory <- ifelse(IAA.node2.detail.info.essence$GeneID%in% nitrate.reduc.genes$GeneID,"Nitrate reduction",
                                                    ifelse(IAA.node2.detail.info.essence$GeneID%in%NH3transport.genes$GeneID,"Ammonium transporter",
                                                           ifelse(IAA.node2.detail.info.essence$GeneID%in%glutamate.produc.genes$GeneID,"Glutamate biosynthesis",
                                                                  ifelse(IAA.node2.detail.info.essence$GeneID%in%Tryptophan.produc.genes$GeneID,"Tryptophan biosynthesis","IAA biosynthesis"))))



IAA.node2.detail.info.essence$FuncCategory2 <- ifelse(IAA.node2.detail.info.essence$Domain == "Bacteria", paste0("Bac_",IAA.node2.detail.info.essence$FuncCategory),
                                                     paste0("Fun_",IAA.node2.detail.info.essence$FuncCategory))
write.table(IAA.node2.detail.info.essence,"./2_0.2_node_IAA_additionalInfo.tsv", sep = "\t",quote=F, row.names = F, col.names = T)
nrow(IAA.edge2)

######## Make subnetwork for main figure ########
##### Networks related to Tm
##### Select links associated with IAA production (EC 2.6.1.27, 4.1.1.28)
Tmgenes<-IAA.node2.detail.info.essence$NewID[which(IAA.node2.detail.info.essence$Species == "Tricholoma matsutake")]
net.centered.Tm<-IAA.edge2[IAA.edge2$Source%in%Tmgenes|IAA.edge2$Target%in%Tmgenes,]
unique(c(net.centered.Tm$Source,net.centered.Tm$Target))
write.table(net.centered.Tm,"./2_0.2_edge_centered_Tm.tsv", sep = "\t",quote=F, row.names = F, col.names = T)

IAA.node2.detail.info.essence.Tm<-subset(IAA.node2.detail.info.essence, NewID %in% unique(c(net.centered.Tm$Source,net.centered.Tm$Target)))
write.table(IAA.node2.detail.info.essence.Tm,"./2_0.2_node_centered_Tm.tsv", sep = "\t",quote=F, row.names = F, col.names = T)




###### Module-level analyses ######
library(igraph)

IAA.edge2 <- read.table("./2_0.2_edge_IAA.tsv",sep="\t",header = T)

nodeattribIAA_net <- data.frame(node=union(IAA.edge2$Source,IAA.edge2$Target))
nodeattribIAA_net$kingdom <- 0

for (i in as.character(nodeattribIAA_net$node))
{
  if (i%in%nodeattribIAA_net$node[(grep("^B", nodeattribIAA_net$node))] == TRUE)
  {nodeattribIAA_net[nodeattribIAA_net$node==i,"kingdom"] <- "Bacteria"}
  
  else
  {nodeattribIAA_net[nodeattribIAA_net$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattribIAA_net) <- as.character(nodeattribIAA_net$node)
names(IAA.node2.detail.info.essence)
names(IAA.edge2)
head(IAA.node2.detail.info.essence)

nodeattribIAA_net<-IAA.node2.detail.info.essence
names(nodeattribIAA_net)[1] <- "node"
all_IAA_net<- graph_from_data_frame(IAA.edge2,direct=F, vertices=nodeattribIAA_net)

## Number of nodes
length(V(all_IAA_net))  #4,048

## Number of bacteria and fungi nodes
length(grep("Bac",names(V(all_IAA_net)))) #3891
length(grep("Fun",names(V(all_IAA_net)))) #157


## Connections 
bb_occur_IAA <- droplevels(IAA.edge2[with(IAA.edge2, grepl("Bac",Source) & grepl("Bac",Target)),])
nrow(bb_occur_IAA) #57,981

ff_occur_IAA <- droplevels(IAA.edge2[with(IAA.edge2, grepl("Fun",Source) & grepl("Fun",Target)),])
nrow(ff_occur_IAA) #623

fb_occur_IAA <- droplevels(IAA.edge2[with(IAA.edge2, grepl("Fun",Source) & grepl("Bac",Target)),])
nrow(fb_occur_IAA) #5,912

##### Extract modules
cluster.louvain_IAA <- cluster_louvain(as.undirected(all_IAA_net))

print(membership(cluster.louvain_IAA))
unique(cluster.louvain_IAA$membership) ## Five modules were identified

# Get membership of each node
IAA.membership <- cluster.louvain_IAA$membership

# Loop through each module to extract subnetworks
subnetworks <- list()
for (mod in unique(IAA.membership)) {
  # Get nodes in the current module
  nodes_in_module <- which(IAA.membership == mod)
  
  # Create a subgraph for the module
  subgraph <- induced_subgraph(all_IAA_net, vids = nodes_in_module)
  
  # Store the subgraph
  subnetworks[[paste0("Module_", mod)]] <- subgraph
}

# Example: Access a specific subnetwork
module1_graph <- subnetworks[["Module_1"]]
module2_graph <- subnetworks[["Module_2"]]
module3_graph <- subnetworks[["Module_3"]]
module4_graph <- subnetworks[["Module_4"]]
module5_graph <- subnetworks[["Module_5"]]

##### Find modules containing Tm genes
names(V(module1_graph))
nrow(IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module1_graph))& 
                                       IAA.node2.detail.info.essence$Species == "Tricholoma matsutake",]) #46
IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module2_graph))& 
                                IAA.node2.detail.info.essence$Species == "Tricholoma matsutake",]
IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module3_graph))& 
                                IAA.node2.detail.info.essence$Species == "Tricholoma matsutake",] #3
IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module4_graph))& 
                                IAA.node2.detail.info.essence$Species == "Tricholoma matsutake",] #1
IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module5_graph))& 
                                IAA.node2.detail.info.essence$Species == "Tricholoma matsutake",] #5
#### Module 1 contained the highest number of Tm genes compared to other modules.
##### Details in module 1
node.module1<-IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module1_graph)),]
unique(node.module1$FuncCategory2)
writexl::write_xlsx(node.module1,"./Info_Module1.xlsx")
unique(node.module1$Species)
edge.module1<-IAA.edge2[IAA.edge2$Source%in% node.module1$NewID &IAA.edge2$Target%in% node.module1$NewID,]
writexl::write_xlsx(edge.module1,"./Info_Module1_edge.xlsx")



##### Assign module information in the node and edge tables
IAA.edge.info <- read.table("./2_0.2_edge_IAA.tsv",sep="\t",header = T)

IAA.node2.detail.info.essence$Module <- 0
IAA.node2.detail.info.essence$Module[which(IAA.node2.detail.info.essence$NewID %in% names(V(module1_graph)))] <- "Module1"
IAA.node2.detail.info.essence$Module[which(IAA.node2.detail.info.essence$NewID %in% names(V(module2_graph)))] <- "Module2"
IAA.node2.detail.info.essence$Module[which(IAA.node2.detail.info.essence$NewID %in% names(V(module3_graph)))] <- "Module3"
IAA.node2.detail.info.essence$Module[which(IAA.node2.detail.info.essence$NewID %in% names(V(module4_graph)))] <- "Module4"
IAA.node2.detail.info.essence$Module[which(IAA.node2.detail.info.essence$NewID %in% names(V(module5_graph)))] <- "Module5"


write.table(IAA.node2.detail.info.essence,"./2_0.2_node_IAA_additionalInfo_module.tsv", sep = "\t",quote=F, row.names = F, col.names = T)


node.module1<-IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module1_graph)),]
node.module2<-IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module2_graph)),]
node.module3<-IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module3_graph)),]
node.module4<-IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module4_graph)),]
node.module5<-IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module5_graph)),]

IAA.edge2$module <-0
IAA.edge2$module[which(IAA.edge2$Source%in% node.module1$NewID &IAA.edge2$Target%in% node.module1$NewID)] <-"module1"
IAA.edge2$module[which(IAA.edge2$Source%in% node.module2$NewID &IAA.edge2$Target%in% node.module2$NewID)] <-"module2"
IAA.edge2$module[which(IAA.edge2$Source%in% node.module3$NewID &IAA.edge2$Target%in% node.module3$NewID)] <-"module3"
IAA.edge2$module[which(IAA.edge2$Source%in% node.module4$NewID &IAA.edge2$Target%in% node.module4$NewID)] <-"module4"
IAA.edge2$module[which(IAA.edge2$Source%in% node.module5$NewID &IAA.edge2$Target%in% node.module5$NewID)] <-"module5"
IAA.edge2$module[which(IAA.edge2$module==0)] <-"Between-module association"

write.table(IAA.edge2,"./2_0.2_edge_module.tsv", sep = "\t",quote=F, row.names = F, col.names = T)

##### Module composition
module.comp.tab <- IAA.node2.detail.info.essence
module.comp.tab$Count <- 1

module.comp.tab.summary <- module.comp.tab %>% group_by(Module, Domain, Kingdom, Phylum, Class,
                                                        Order, Family, Genus, Species) %>% reframe(Count = sum(Count))

module.comp.tab.summary <- module.comp.tab %>% group_by(Module, Domain, Kingdom, Phylum, Class,
                                                        Order, Family, Genus) %>% reframe(Count = sum(Count))


module.comp.tab.summary <- module.comp.tab.summary %>%group_by(Module) %>%mutate(Total = sum(Count))

module.comp.tab.summary$Ratio <- module.comp.tab.summary$Count/module.comp.tab.summary$Total
writexl::write_xlsx(module.comp.tab.summary,"./IAAnet_Module composition_raw table.xlsx")


module1.comp <- module.comp.tab.summary[module.comp.tab.summary$Module == "Module1",]

#### Directly associated with Tm genes
module1.Tm.genes<-IAA.node2.detail.info.essence[IAA.node2.detail.info.essence$NewID %in% names(V(module1_graph))& 
                                                  IAA.node2.detail.info.essence$Species == "Tricholoma matsutake",]


asso.Tm.bac<-edge.module1[edge.module1$Source %in% module1.Tm.genes$NewID &grepl("Bac",edge.module1$Target),]
nrow(asso.Tm.bac) #1,942 associations

gene.list.Tm<-unique(c(asso.Tm.bac$Source,asso.Tm.bac$Target))
bac.gene.asso.Tm<-gene.list.Tm[-(which(gene.list.Tm%in%module1.Tm.genes$NewID))]#405 genes

module1.bac.genes <- node.module1[node.module1$NewID %in% bac.gene.asso.Tm,]


module1.bac.genes$Count <- 1
module1.bac.dir.asso.genes.summary <- module1.bac.genes %>% group_by(Domain, Kingdom, Phylum, Class,
                                                        Order, Family, Genus) %>% reframe(Count = sum(Count))

module1.bac.dir.asso.genes.summary$Ratio <- module1.bac.dir.asso.genes.summary$Count/sum(module1.bac.dir.asso.genes.summary$Count)
writexl::write_xlsx(module.comp.tab.summary,"./IAAnet_Module composition_raw table.xlsx")
print(module1.bac.dir.asso.genes.summary, n=60)

print(module1.comp, n=80)
names(module1.comp)

module1.comp <- module1.comp[-c(1,10)]


module1.comp$Group <- "All"
module1.bac.dir.asso.genes.summary$Group <- "Directly linked"
dir.bac.helper.list<-module1.bac.dir.asso.genes.summary$Genus
save(dir.bac.helper.list,file="../dir.bac.helper.list.RData")

tab.plot.module.comp<-rbind(module1.comp,module1.bac.dir.asso.genes.summary)
unique(tab.plot.module.comp$Genus)

#### Plotting ###
head(tab.plot.module.comp)
tab.plot.module.comp$Genus2 <- tab.plot.module.comp$Genus
tab.plot.module.comp$Genus2[tab.plot.module.comp$Count == 1] <- "Others"
tab.plot.module.comp$Genus2[which(tab.plot.module.comp$Domain == "Bacteria"&tab.plot.module.comp$Genus2=="Others")] <-"Bacteria_others"
tab.plot.module.comp$Genus2[which(tab.plot.module.comp$Kingdom == "Fungi"&tab.plot.module.comp$Genus2=="Others")] <-"Fungi_others"
unique(tab.plot.module.comp$Genus2)
writexl::write_xlsx(tab.plot.module.comp,"./module1 taxonomic composition.xlsx")

### Ordering genera
order.module1 <- tab.plot.module.comp %>% group_by(Domain, Kingdom,Genus,Genus2) %>% reframe(Ratio = mean(Ratio)) %>% arrange(desc(Ratio))
order.module1.bac <- order.module1[order.module1$Kingdom != "Fungi",]$Genus2
order.module1.fun <- order.module1[order.module1$Kingdom == "Fungi",]$Genus2
order.module1.list<- unique(c(order.module1.bac,order.module1.fun))
order.module1.list <- order.module1.list[c(2:33,1,34:38)]
tab.plot.module.comp$Genus2 <- factor(tab.plot.module.comp$Genus2, levels = rev(order.module1.list))

ggplot(tab.plot.module.comp, aes(x=Group, y = Ratio, fill = Genus2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Alphaproteobacteria_unknown" = "#1f77b4",
                               "Acidobacteriota_unknown" = "#ff7f0e",
                               "Planctomycetota_unknown" = "#2ca02c",
                               "Pseudomonadota_unknown" = "#d62728",
                               "Actinomycetota_unknown" = "#9467bd",
                               "Chloroflexota_unknown" = "#8c564b",
                               "Planctomycetia_unknown" = "#e377c2",
                               "Verrucomicrobiota_unknown" = "#7f7f7f",
                               "Burkholderiaceae_unknown" = "#bcbd22",
                               "Hyphomicrobiales_unknown" = "#17becf",
                               "Actinomycetes_unknown" = "#aec7e8",
                               "Bradyrhizobium" = "#ffbb78",
                               "Aliidongia" = "#98df8a",
                               "Nitrobacteraceae_unknown" = "#ff9896",
                               "Trebonia" = "#c5b0d5",
                               "Gammaproteobacteria_unknown" = "#c49c94",
                               "Edaphobacter" = "#f7b6d2",
                               "Streptosporangiales_unknown" = "#c7c7c7",
                               "Gemmatimonadota_unknown" = "#dbdb8d",
                               "Paraburkholderia" = "#9edae5",
                               "Deltaproteobacteria_unknown" = "#3182bd",
                               "Isosphaeraceae_unknown" = "#e6550d",
                               "Betaproteobacteria_unknown" = "#31a354",
                               "Singulisphaera" = "#756bb1",
                               "Caballeronia" = "#636363",
                               "Candidatus Rokubacteria_unknown" = "#fd8d3c",
                               "Fimbriiglobus" = "#6baed6",
                               "Gemmataceae_unknown" = "#bcbddc",
                               "Gemmatales_unknown" = "#d9d9d9",
                               "Mycobacteriaceae_unknown" = "#969696",
                               "Pseudonocardiaceae_unknown" = "#fdae6b",
                               "Bacteria_others" = "#6baed6",
                               "Bacteria_unknown" = "#474747",
                               "Tricholoma" = "#fdbe85",
                               "Agaricomycetes_unknown" = "#74c476",
                               "Agaricales_unknown" = "#e7298a",
                               "Ascomycota_unknown" = "#ffed6f",
                               "Fungi_others" = "#999999")) +
  
  xlab('')+
  ylab("Ratio \n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(ncol = 4,reverse = T))+
  theme(legend.position="bottom") + theme(aspect.ratio = 1.2)+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())



module1.bac.genes[module1.bac.genes$Genus == "Bradyrhizobium",]
module1.bac.genes[module1.bac.genes$Genus == "Caballeronia",]
module1.bac.genes[module1.bac.genes$Genus == "Trebonia",]
module1.bac.genes[module1.bac.genes$Genus == "Paraburkholderia",]
module1.bac.genes[module1.bac.genes$Family == "Nitrobacteraceae",]
module1.bac.genes[module1.bac.genes$Family == "Paraburkholderia",]


NH3.reduc.gene.bac.net<-module1.bac.genes[module1.bac.genes$FuncCategory2 == "Bac_Nitrate reduction",]$NewID
NH3.asso.TmGene<-unique(edge.module1[edge.module1$Source %in% module1.Tm.genes$NewID &edge.module1$Target %in% NH3.reduc.gene.bac.net,]$Source)
node.module1[node.module1$NewID %in% NH3.asso.TmGene,]

unique(module1.bac.genes$Genus)
unique(node.module1[node.module1$NewID %in% NH3.reduc.gene.bac.net,]$Genus)

glutamate.biosyn.gene.bac.net<-module1.bac.genes[module1.bac.genes$FuncCategory2 == "Bac_Glutamate biosynthesis",]$NewID
unique(node.module1[node.module1$NewID %in% glutamate.biosyn.gene.bac.net,]$Genus)
glutamate.biosyn.asso.TmGene<-unique(edge.module1[edge.module1$Source %in% module1.Tm.genes$NewID &edge.module1$Target %in% glutamate.biosyn.gene.bac.net,]$Source)
node.module1[node.module1$NewID %in%glutamate.biosyn.asso.TmGene,]

Tryptophan.biosyn.gene.bac.net<-module1.bac.genes[module1.bac.genes$FuncCategory2 == "Bac_Tryptophan biosynthesis",]$NewID
unique(node.module1[node.module1$NewID %in% Tryptophan.biosyn.gene.bac.net,]$Genus)


IAA.biosyn.gene.bac.net<-module1.bac.genes[module1.bac.genes$FuncCategory2 == "Bac_IAA biosynthesis",]$NewID
unique(node.module1[node.module1$NewID %in% IAA.biosyn.gene.bac.net,]$Genus)

IAA.biosyn.gene.bac.net.info<-module1.bac.genes[module1.bac.genes$FuncCategory2 == "Bac_IAA biosynthesis",]
IAA.biosyn.gene.bac.net.info[grepl("4.1.1.74|1.4.3.4|1.4.3.22",IAA.biosyn.gene.bac.net.info$EC),]
IAA.biosyn.gene.bac.net.info[grepl("1.13.12.3",IAA.biosyn.gene.bac.net.info$EC),]
IAA.biosyn.gene.bac.net.info[grepl("1.14.13.168",IAA.biosyn.gene.bac.net.info$EC),]


## Network properties
meta_degree <- sort(igraph::degree(all_IAA_net,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree))) #"max degree = " "69"     
print(c('mean degree = ', mean(meta_degree))) #"mean degree = "   "15.04"
meta_degree <- as.data.frame(meta_degree)
ggplot(meta_degree, aes(x=meta_degree)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_IAA_net, directed = FALSE)))

print(paste("mean clustering coefficient = ", igraph::transitivity(all_IAA_net, "global")))

print(paste("mean betweenness centrality = ", mean(betweenness(all_IAA_net, directed = FALSE, normalized = TRUE))))

print(paste("mean closeness centrality = ", mean(closeness(all_IAA_net, normalized = TRUE))))

print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_IAA_net, V(all_IAA_net))))))

##

net <- all_IAA_net
IAA_all_deg <- igraph::degree(net,mode="all")
IAA_all_betweenness <- betweenness(net, normalized = TRUE)
IAA_all_closeness <- closeness(net, normalized = TRUE)
IAA_all_transitivity <- igraph::transitivity(net, "local", vids = V(net))
names(IAA_all_transitivity)<- V(net)$name
IAA_all_transitivity[is.na(IAA_all_transitivity)] <- 0

### network properties of bacteria and fungi
df.IAA.degree<-data.frame(IAA_all_deg)
head(df.IAA.degree)
df.IAA.degree$Group <- "IAA"
names(df.IAA.degree)[1] <- c("Degree")

df.IAA.closeness<-data.frame(IAA_all_closeness)
head(df.IAA.closeness)
df.IAA.closeness$Group <- "IAA"
names(df.IAA.closeness)[1] <- c("Closeness")

df.IAA.betweenness<-data.frame(IAA_all_betweenness)
head(df.IAA.betweenness)
df.IAA.betweenness$Group <- "IAA"
names(df.IAA.betweenness)[1] <- c("Betweenness")

df.IAA.degree$kingdom <- ifelse(grepl("^B",rownames(df.IAA.degree)),'Bacteria', 'Fungi')
df.IAA.closeness$kingdom <- ifelse(grepl("^B",rownames(df.IAA.closeness)),'Bacteria', 'Fungi')
df.IAA.betweenness$kingdom <- ifelse(grepl("^B",rownames(df.IAA.betweenness)),'Bacteria', 'Fungi')

shapiro.test(df.IAA.degree$Degree) #p-value < 2.2e-16, not normal
shapiro.test(df.IAA.betweenness$Betweenness) #p-value < 2.2e-16, not normal
shapiro.test(df.IAA.closeness$Closeness) #p-value = 4.498e-07, not normal

##Degree
p1 <- print_degree(df.IAA.degree)
p1

## Betweenness
p2 <- print_betweenness(df.IAA.betweenness)
p2

## Closeness
p3 <- print_closeness(df.IAA.closeness)
p3

gridExtra::grid.arrange(p1, p2, p3, ncol=3)

dev.off()

### Within and among-connectivity hub estimation
library(brainGraph)
library(igraph)

network_hub_z_p<-function(net.soil.spar, keyword){####Assigning community membership (module)
  ## Perform cluster analysis using greedy clustering algorithm 
  cfg_soil <- cluster_fast_greedy(as.undirected(net.soil.spar))
  
  ### Among-module connectivity
  print(membership(cfg_soil))
  Pi<-part_coeff(net.soil.spar, membership(cfg_soil))
  z_p_table<-data.frame(Pi)
  head(data.frame(Pi))
  ### within-module connectivity
  Zi<-within_module_deg_z_score(net.soil.spar, membership(cfg_soil))
  z_p_table$zi <- Zi
  
  z_p_table$zi[is.na(z_p_table$zi)] <- 0
  z_p_table$Pi[is.na(z_p_table$Pi)] <- 0
  ###Assigning Kingdom
  z_p_table$kingdom <- 0
  z_p_table$kingdom[grep("^B",rownames(z_p_table))] <- "Bacteria"
  z_p_table$kingdom[grep("^F",rownames(z_p_table))] <- "Fungi"
  
  z_p_table$Role <- 0
  z_p_table$Role[z_p_table$zi >= 2.5 & z_p_table$Pi >= 0.62] <- "Network hub"
  z_p_table$Role[z_p_table$zi < 2.5 & z_p_table$Pi < 0.62] <- "Peripheral"
  z_p_table$Role[z_p_table$zi < 2.5 & z_p_table$Pi >= 0.62] <- "Connector"
  z_p_table$Role[z_p_table$zi >= 2.5 & z_p_table$Pi < 0.62] <- "Module hub"
  
  
  ### Assigning Module numbers
  cluster.raw<- membership(cfg_soil)
  
  df.cluster <- data.frame(matrix(ncol = 1, nrow = 0))
  names(df.cluster)[1] <- "Module"
  
  for(i in 1:length(cluster.raw)){
    df.add<-data.frame(cluster.raw[i])
    names(df.add)[1] <- "Module"
    df.cluster <- rbind(df.cluster, df.add)
  }
  
  z_p_table <- merge(z_p_table, df.cluster, by = "row.names")
  names(z_p_table)[1] <- "Node"
  write.csv(z_p_table, paste0("z_p_table_igraph_",keyword,".csv"))
  
  ##plotting
  p<-ggplot(z_p_table, aes(x=Pi, y=zi, color=kingdom)) +
    xlab('\n Among-module connectivity (Pi)')+
    ylab("Within-module connectivity (Zi) \n") +
    geom_point(size=3, alpha=0.7) +
    theme(aspect.ratio = 1)+
    # ggtitle("Volcano Plot \n") +
    theme(legend.text=element_text(size=13)) + 
    theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
    guides(size="none") +
    #scale_x_continuous(breaks=seq(0,1,0.2))+
    #scale_y_continuous(breaks=seq(-20,0,-5))+
    scale_color_manual(values=c("Bacteria"="#EE7600","Fungi"="#9A32CD")) + 
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_rect(colour = "black", fill = "white"), plot.background=element_blank())+
    geom_hline(yintercept=2.5, linetype='dashed', color='black', linewidth = 0.75)+geom_vline(xintercept=0.62, linetype='dashed', color='black', linewidth = 0.75)
  
  return(p)
}


network_hub_z_p(all_IAA_net, "IAA")
network_hub_z_p(all_TmD_net, "TmD")

dev.off()






