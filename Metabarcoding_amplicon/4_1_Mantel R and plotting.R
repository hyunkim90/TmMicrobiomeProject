############# CCA and Mantel's R ############
########## 1. CCA ###########
bac.clean.log.chemicals
fun.clean.log.chemicals

taxa_data.b <- data.frame(otu_table(bac.clean.log.chemicals))
taxa_data.f <- data.frame(otu_table(fun.clean.log.chemicals))
### Excellent reference: https://rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html#google_vignette
set.seed(123)
cca_result <- cca(formula = t(taxa_data.f) ~., data = env_data)
finalmodel<- ordistep(cca_result, scope=formula(cca_result))
#### Check variacne inflation factor (vif values should be lower than 10)
vif.cca(finalmodel)
#### Check model significance (p value should be lower than 0.05)
anova.cca(finalmodel)
anova.cca(finalmodel, by="terms")
#### Plotting
# Plot the basic CCA results without sites
plot(finalmodel, display = c("species", "bp"), type = "n") # type = "n" creates an empty plot
# Add species points with specified shapes
species_scores <- scores(finalmodel, display = "species")
points(species_scores, pch = 19, col = "#A63446") # pch = 19 is for filled circles
# Extract environmental variable scores
env_scores <- scores(finalmodel, display = "bp")
# Add arrows for the environmental variables
arrows(0, 0, env_scores[, 1], env_scores[, 2], col = "#0C6291", length = 0.1)
# Add labels for environmental variables
text(env_scores[, 1], env_scores[, 2], labels = rownames(env_scores), col = "#0C6291", pos = 3, cex = 0.8)
# pch = 19 is for filled circles


########## 2. Mantel's R ###########
# prepare data
devtools::install_github("Github-Yilei/ggcor")
install.packages("file2meco")

library(microeco)
library(magrittr)
library(file2meco)

dataset.b.chem<-phyloseq2meco(bac.clean.log.chemicals)
dataset.f.chem<-phyloseq2meco(fun.clean.log.chemicals)
# extract two phyla to show the steps
d1 <- clone(dataset.b.chem)
d1$tidy_dataset()
d1$cal_betadiv()
d2 <- clone(dataset.f.chem)
d2$tidy_dataset()
d2$cal_betadiv()
d1$sample_table
# first perform mantel test
t1 <- trans_env$new(dataset = d1, env_cols = 6:19)
t1$cal_mantel(use_measure = "bray", partial_mantel = TRUE, method = "spearman")
t1$res_mantel
t2 <- trans_env$new(dataset = d2, env_cols = 5:18)
t2$cal_mantel(use_measure = "bray", partial_mantel = TRUE, method = "spearman")
t2$res_mantel
# extract a part of the results 
x1 <- data.frame(spec = "Bacteria", t1$res_mantel) %>% .[, c(1, 3, 6, 8)]
x2 <- data.frame(spec = "Fungi", t2$res_mantel) %>% .[, c(1, 3, 6, 8)]
# rename columns
colnames(x1) <- colnames(x2) <- c("spec", "env", "r", "p.value")
# generate interval data
x1 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
x2 %<>% dplyr::mutate(rd = cut(r, breaks = c(-Inf, 0.3, 0.6, Inf), labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
                      pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

# cobine two tables
plot_table <- rbind(x1, x2)
# install ggcor following the steps (https://chiliubio.github.io/microeco_tutorial/intro.html#github-packages)
library(ggplot2)
library(ggcor)
set_scale()

g1 <- quickcor(t1$data_env, type = "upper", cor.test = TRUE, show.diag = FALSE) +
  geom_square() +
  geom_mark(sig.thres = 0.05, markonly = TRUE, color = "black", size = 6) +
  anno_link(aes(colour = pd, size = rd), data = plot_table) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3))

g1


