#!/usr/bin/env Rscript

### Plot bootstrap distribution of TE association 
#############################################################################
# The raw data I collected from the output of the Snakemake pipelines TreeCertainty.smk and Phylogeny.smk
# =======================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020-01-13
# Version 1
# =======================================
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE)
# library(cowplot) 
# library(reshape2) # for melt

# ---------
# Data
# ---------

# P. anserina
observed_an <- read.table("/Users/lore/Dropbox/PhD_UU/Analyses/JohannessonsServer/3_SpokBlockPaper/4_BootstrapTEs/distributions/PaWa137m-TEs-1000bp_obsvalues.txt", header = TRUE)
randvalues_an <- read.table("/Users/lore/Dropbox/PhD_UU/Analyses/JohannessonsServer/3_SpokBlockPaper/4_BootstrapTEs/distributions/PaWa137m-TEs-1000bp_randvalues_10000.txt", header = TRUE)

# P. pauciseta
observed_pa <- read.table("/Users/lore/Dropbox/PhD_UU/Analyses/JohannessonsServer/3_SpokBlockPaper/4_BootstrapTEs/distributions/CBS237.71m-TEs-1000bp_obsvalues.txt", header = TRUE)
randvalues_pa <- read.table("/Users/lore/Dropbox/PhD_UU/Analyses/JohannessonsServer/3_SpokBlockPaper/4_BootstrapTEs/distributions/CBS237.71m-TEs-1000bp_randvalues_10000.txt", header = TRUE)

# P. comata
observed_co <- read.table("/Users/lore/Dropbox/PhD_UU/Analyses/JohannessonsServer/3_SpokBlockPaper/4_BootstrapTEs/distributions/PcWa139m-TEs-1000bp_obsvalues.txt", header = TRUE)
randvalues_co <- read.table("/Users/lore/Dropbox/PhD_UU/Analyses/JohannessonsServer/3_SpokBlockPaper/4_BootstrapTEs/distributions/PcWa139m-TEs-1000bp_randvalues_10000.txt", header = TRUE)
# ---------

# What is the mean of the observed values in each species?
mean_an <- observed_an$Coverage %>% mean()
mean_pa <- observed_pa$Coverage %>% mean()
mean_co <- observed_co$Coverage %>% mean()

# Make a big data frame with all species

randvalues <- rbind(cbind(Strain = "PaWa137m", randvalues_an),
                    cbind(Strain = "CBS237.71m", randvalues_pa),
                    cbind(Strain = "PcWa139m", randvalues_co))

# Make a dataframe to add each mean of observed values to each histogram
meansframe <- data.frame(Strain = c("PaWa137m", "CBS237.71m", "PcWa139m"), Mean = c(mean_an, mean_pa, mean_co))
quantile95 <- data.frame(Strain = c("PaWa137m", "CBS237.71m", "PcWa139m"), 
                         quantile95 = c(quantile(randvalues_an$Coverage, probs = 0.95)[[1]], 
                                        quantile(randvalues_pa$Coverage, probs = 0.95)[[1]], 
                                        quantile(randvalues_co$Coverage, probs = 0.95)[[1]]))
# Plot
ggplot(randvalues, aes(x = Coverage)) + geom_histogram() + 
  geom_vline(aes(xintercept = Mean), data = meansframe, colour = "red", linetype="dashed") +
  geom_vline(aes(xintercept = quantile95), data = quantile95, colour = "darkred") +
  theme_bw() + ylab("Frequency") + facet_grid(Strain ~ .) 



