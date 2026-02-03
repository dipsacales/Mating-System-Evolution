setwd("~/Desktop/Projects/Valeriana_Mating_System")

#install.packages("phyloglm")
#install.packages("phytools")
#library(phyloglm)
library(ape)
library(phytools)
library(geiger)
library(phylolm)

###LOAD TREE
time_tree <- read.tree(file= "time_tree_valeriana.nwk")
tip <- c("Valeriana_magellanica", "Valeriana_retrorsa", "Valeriana_bulbosa","Valeriana_lancifolia")
tree <-drop.tip(time_tree, tip)

### LOAD DATA
df1 <- read.csv(file = "~/Desktop/Projects/Valeriana_Mating_System/traits.csv") #mating system gyno or NOT 
df2 <- read.csv(file = "~/Desktop/Projects/Valerianaceae_biogeography/species_environmental_traits2.csv")

### Merge data.frames
traits <- merge(df1, df2)

### Match tree & data
tree$tip.label <- gsub("_", " ", tree$tip.label)
traits$species <- gsub("_", " ", traits$species)

shared_species <- intersect(tree$tip.label, traits$species)
length(shared_species)

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, shared_species))

library(tidyr)
library(dplyr)

traits_pruned <- traits %>%
  filter(species %in% shared_species) %>%
  arrange(match(species, tree_pruned$tip.label))

### Does tree mach dataset??? TRUE??
all(tree_pruned$tip.label == traits_pruned$species)


### ADD, , mean_pwm, mean_ps
BIO1  <- setNames(traits_pruned$mean_bio1, traits_pruned$species)
BIO2  <- setNames(traits_pruned$mean_bio2, traits_pruned$species)
BIO3  <- setNames(traits_pruned$mean_bio3, traits_pruned$species)
BIO4  <- setNames(traits_pruned$mean_bio4, traits_pruned$species)
BIO8  <- setNames(traits_pruned$mean_bio8, traits_pruned$species)
BIO9  <- setNames(traits_pruned$mean_bio9, traits_pruned$species)
BIO13  <- setNames(traits_pruned$mean_bio13, traits_pruned$species)
BIO14  <- setNames(traits_pruned$mean_bio14, traits_pruned$species)
BIO15  <- setNames(traits_pruned$mean_bio15, traits_pruned$species)
BIO19  <- setNames(traits_pruned$mean_bio19, traits_pruned$species)
Elev <- setNames(traits_pruned$mean_elevation, traits_pruned$species)


tree_pruned$edge.length <- tree_pruned$edge.length +0.0000001

dev.off()
### Plot Characters
contMap(tree_pruned, Elev, plot = TRUE, legend = 5, lwd = 2,  fsize=c(0.5,0.5), leg.txt="" , sig = 0, res =100)
title("Mean Elevation (°meters)")

#### Calculate Phylogenetic Signal
### BIO1
phylosig(tree_pruned, BIO1, method = "K", test = TRUE)
phylosig(tree_pruned, BIO1, method = "lambda", test = TRUE)

### BIO2
phylosig(tree_pruned, BIO2, method = "K", test = TRUE)
phylosig(tree_pruned, BIO2, method = "lambda", test = TRUE)

### BIO3
phylosig(tree_pruned, BIO3, method = "K", test = TRUE)
phylosig(tree_pruned, BIO3, method = "lambda", test = TRUE)

### BIO4
phylosig(tree_pruned, BIO4, method = "K", test = TRUE)
phylosig(tree_pruned, BIO4, method = "lambda", test = TRUE)

###BIO8
phylosig(tree_pruned, BIO8, method = "K", test = TRUE)
phylosig(tree_pruned, BIO8, method = "lambda", test = TRUE)

###BIO9
phylosig(tree_pruned, BIO9, method = "K", test = TRUE)
phylosig(tree_pruned, BIO9, method = "lambda", test = TRUE)

###BIO13
phylosig(tree_pruned, BIO13, method = "K", test = TRUE)
phylosig(tree_pruned, BIO13, method = "lambda", test = TRUE)

###BIO14
phylosig(tree_pruned, BIO14, method = "K", test = TRUE)
phylosig(tree_pruned, BIO14, method = "lambda", test = TRUE)

###BIO15
phylosig(tree_pruned, BIO15, method = "K", test = TRUE)
phylosig(tree_pruned, BIO15, method = "lambda", test = TRUE)

###BIO19
phylosig(tree_pruned, BIO19, method = "K", test = TRUE)
phylosig(tree_pruned, BIO19, method = "lambda", test = TRUE)

### Elev
phylosig(tree_pruned, Elev, method = "K", test = TRUE)
phylosig(tree_pruned, Elev, method = "lambda", test = TRUE)


## K ≈ 1, λ ≈ 1	Brownian motion (strong phylogenetic niche conservatism)
## K < 1, λ < 1	Convergence / lability
## K > 1	Strong clade-level niche conservatism
## λ ≈ 0	Trait independent of phylogeny


### Fit Brownian Motion (BM) vs Ornstein–Uhlenbeck (OU)#######################
##############################################################################
###Make sure tree is ulrametric!!! ###########################################
##############################################################################
is.ultrametric(tree_pruned)

##########################################################
## fix ultrametric tree:                                ##
## http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html ##
##########################################################

force.ultrametric <- function(tree, method = c("nnls", "extend")) {
  method <- method[1]
  if (method == "nnls") 
    tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE, trace = 0)
  else if (method == "extend") {
    h <- diag(vcv(tree))
    d <- max(h) - h
    ii <- sapply(1:Ntip(tree), function(x, y) which(y == x),
                 y = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

### Add tree_pruned to function #####

tree_ult <- ladderize(force.ultrametric(tree_pruned, method = "extend"), F)
plot(tree_ult)
is.ultrametric(tree_ult)

tree_pruned <- tree_ult
is.ultrametric(tree_pruned)

##############################################################################
##############################################################################

### BIO1
bm_BIO1 <- fitContinuous(tree_pruned, BIO1, model = "BM")
ou_BIO1 <- fitContinuous(tree_pruned, BIO1, model = "OU")

### BIO2
bm_BIO2 <- fitContinuous(tree_pruned, BIO2, model = "BM")
ou_BIO2 <- fitContinuous(tree_pruned, BIO2, model = "OU")

### BIO3
bm_BIO3 <- fitContinuous(tree_pruned, BIO3, model = "BM")
ou_BIO3 <- fitContinuous(tree_pruned, BIO3, model = "OU")

### BIO4
bm_BIO4 <- fitContinuous(tree_pruned, BIO4, model = "BM")
ou_BIO4 <- fitContinuous(tree_pruned, BIO4, model = "OU")

### BIOS8
bm_BIO8 <- fitContinuous(tree_pruned, BIO8, model = "BM")
ou_BIO8 <- fitContinuous(tree_pruned, BIO8, model = "OU")

### BIOS9
bm_BIO9 <- fitContinuous(tree_pruned, BIO9, model = "BM")
ou_BIO9 <- fitContinuous(tree_pruned, BIO9, model = "OU")

### BIOS18
bm_BIO13 <- fitContinuous(tree_pruned, BIO13, model = "BM")
ou_BIO13 <- fitContinuous(tree_pruned, BIO13, model = "OU")

### BIO14
bm_BIO14 <- fitContinuous(tree_pruned, BIO14, model = "BM")
ou_BIO14 <- fitContinuous(tree_pruned, BIO14, model = "OU")

### BIOS15
bm_BIO15 <- fitContinuous(tree_pruned, BIO15, model = "BM")
ou_BIO15 <- fitContinuous(tree_pruned, BIO15, model = "OU")

### BIOS19
bm_BIO19 <- fitContinuous(tree_pruned, BIO19, model = "BM")
ou_BIO19 <- fitContinuous(tree_pruned, BIO19, model = "OU")

### Elev
bm_elev <- fitContinuous(tree_pruned, Elev, model = "BM")
ou_elev <- fitContinuous(tree_pruned, Elev, model = "OU")

aic_table <- data.frame(
  Trait = c("BIO1", "BIO2","BIO3","BIO4","BIO8","BIO9", "BIO13", "BIO14", "BIO15", "BIO19", "Elevation"),
  BM_AIC = c(bm_BIO1$opt$aic, bm_BIO2$opt$aic, bm_BIO3$opt$aic, bm_BIO4$opt$aic, bm_BIO8$opt$aic, bm_BIO9$opt$aic, bm_BIO13$opt$aic, bm_BIO14$opt$aic, bm_BIO15$opt$aic, bm_BIO19$opt$aic, bm_elev$opt$aic),
  OU_AIC = c(ou_BIO1$opt$aic, ou_BIO2$opt$aic, ou_BIO3$opt$aic, ou_BIO4$opt$aic, ou_BIO8$opt$aic, ou_BIO9$opt$aic, ou_BIO13$opt$aic, ou_BIO14$opt$aic, ou_BIO15$opt$aic, ou_BIO19$opt$aic, ou_elev$opt$aic)
) %>% print()


aic_table$Delta_AIC = aic_table$BM_AIC - aic_table$OU_AIC
aic_table

### ΔAIC (BM − OU)	Conclusion
### > 2	OU better → stabilizing selection / climatic optima
### < −2	BM better → random drift / niche diffusion
### −2 to 2	No strong support

##########################################################
## Run l1ou using separate variables   ###################
##########################################################

#install_github("khabbazian/l1ou")
library(l1ou)

tree <- tree_pruned

plot(tree)
is.ultrametric(tree)

niche_dat <- read.csv("~/Desktop/Projects/Valerianaceae_biogeography/species_environmental_traits2.csv", header = T, row.names = 1)

### NEED to delete n_records from csv table/df

niche_dat <- niche_dat[, c("mean_bio1", "mean_bio2", "mean_bio3", "mean_bio4", "mean_bio8", "mean_bio9", "mean_bio14","mean_bio13","mean_bio15","mean_bio19","mean_elevation")]
tree_ult <- tree

hil <- adjust_data(tree_ult, as.matrix(niche_dat))

colnames(hil$Y) <- c("mean_bio1", "mean_bio2", "mean_bio3", "mean_bio4", "mean_bio8", "mean_bio9", "mean_bio14","mean_bio13","mean_bio15","mean_bio19","mean_elevation")
eModel <- estimate_shift_configuration(hil$tree, hil$Y, criterion = "pBIC", alpha.upper = 8, quietly = F, max.nShifts = 5)
eModel
plot(eModel, edge.shift.ann = F, cex = 0.5)
edgelabels(seq_along(tree$edge[, 1]), adj = c(0.5, -0.5), frame = "n", cex = 0.5, col = "blue")


### , "green", "purple", "white", "pink", "orange", "black", "#88E", "#DDC", "green3", "#DDCC66"

png(filename = "~/Desktop/l1ou.png", width = 1800, height = 1000)
plot(eModel, palette = c("#88CCEE", "#DDCC77", "gray90", "blue", "pink"), 
     asterisk = T, edge.shift.ann = F, cex = 0.75, edge.width = 5.5, label.offset = 0.2, bar.axis = F, 
     x.lim = c(0, 1))

dev.off()

result_pBIC <- l1ou_bootstrap_support(eModel, nItrs = 100, multicore = TRUE, nCores = 4)
result_pBIC$detection.rate




#############################################################################################################################
#############################################################################################################################
### Fucking around with glm/phylolm  ################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

fit <- lm(traits_pruned$mean_bio1 ~ traits_pruned$mean_bio2 + traits_pruned$mean_bio19)
bn_fit <- glm(traits_pruned$gynodioecy ~ traits_pruned$mean_ts + traits_pruned$mean_mat + traits_pruned$mean_elevation, family = binomial)
summary(bn_fit)

#############################################################################################################################
################################### Fucking Around with MCMCglmm ############################################################
#############################################################################################################################
#install.packages("MCMCglmm")
library(MCMCglmm)

?MCMCglmm



############################################################
######## Match variables to tip labels on tree #############
############################################################

tree <- tree_pruned
tree$edge.length <- tree$edge.length + 0.00001

response <-  traits_pruned$gynodioecy
names(response) <- tree$tip.label # Ensure names match tip labels

predictor1 <- traits_pruned$mean_elevation
names(predictor1) <- tree$tip.label # Ensure names match tip labels

predictor2 <- traits_pruned$mean_iso
names(predictor2) <- tree$tip.label # Ensure names match tip labels

predictor3 <- traits_pruned$mean_ts
names(predictor3) <- tree$tip.label # Ensure names match tip labels


##############################################################
### Build data frame with response and predictor variables ###
##############################################################

dat = data.frame(response = response, mean_elev = predictor1, mean_iso = predictor2, mean_ts =predictor3)

mod1 = phyloglm(response ~ log(predictor1),phy=tree,data=dat,boot=100, btol=100, method = "logistic_IG10")
mod2 = phyloglm(response ~ log(predictor2),phy=tree,data=dat,boot=100, btol=500, method = "logistic_IG10")
mod3 = phyloglm(response ~ log(predictor3),phy=tree,data=dat,boot=100, btol=500, method = "logistic_IG10")

Model_aic_table <- data.frame(
  Trait = c("predictor2","predictor2","predictor3"),
  AIC = c(mod1$aic, mod2$aic, mod3$aic),
  Alpha = c(mod1$alpha, mod2$alpha, mod3$alpha)
) %>% print()

### Plot Discrete Traits ################################################
### install ggtree ######################################################
### if (!requireNamespace("BiocManager", quietly = TRUE))################
###  install.packages("BiocManager")#####################################
###cBiocManager::install("ggtree")#######################################
#########################################################################
library(ggtree)
library(ape)
library(phytools)

disc_trait <- traits_pruned$gynodioecy
names(disc_trait) <- traits_pruned$species
levels(disc_trait)

dev.off()
er_model <- fitMk(tree_pruned, disc_trait, model = "ER", pi = "fitzjohn")  # equal rates
gyno_simmap <- simmap(er_model)

ttttt <- summary(gyno_simmap)

plot(ttttt, type= "fan", ftype= "off") # Plot Circle Tree w/ Tip & Node bales from Simmap


#### Another try

tree <- tree_pruned
char_data <-  traits_pruned$gynodioecy
names(char_data) <- tree$tip.label # Ensure names match tip labels

dotTree(tree, char_data, legend = F, ftype = "i", show.tip.label = F)

# Make a copy of the tree and remove tip labels
tree_no_labels <- tree
tree_no_labels$tip.label <- rep("", length(tree_no_labels$tip.label))

# Plot using dotTree with the modified tree object
dotTree(tree_no_labels, char_data, standardize = TRUE, length = 10)



