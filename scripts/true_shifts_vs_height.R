#asd
library(treeio)
library(ape)
library(tibble)
library(ggplot2)

setwd("~/projects/pesto_empirical/")


fpaths <- Sys.glob("data/simulations/age_scaling_effect/*.tre")

sample_fpaths <- sample(fpaths, 100)

indices <- list()
trees <- list()


pb <- txtProgressBar(min = 1, max = length(fpaths), initial = 1)
for (i in seq_along(fpaths)){
  bn <- basename(fpaths[i])
  name1 <- gsub("\\.tre", "", bn)## remove file extension
  tree_index <- strsplit(name1, "_")[[1]][2]
    
  tree <- read.beast.newick(fpaths[i])
  trees[[i]] <- tree
  indices[[i]] <- as.numeric(tree_index)
  setTxtProgressBar(pb,i)
}


N_total <- sapply(trees, function(tree) sum(sapply(tree@data$N, sum)))
tree_length <- sapply(trees, function(tree) sum(tree@phylo$edge.length))
heights <- sapply(trees, function(tree) max(node.depth.edgelength(tree@phylo)))

df <- tibble(
  "N" = N_total,
  "tree_length" = tree_length,
  "height" = heights,
  "N_per_time" = N_total / tree_length,
  "tree_index" = unlist(indices)
)

library(scales)

p <- ggplot(df, aes(x = height, y = N_per_time)) +
  geom_jitter(height = 0, width = 2, size = 0.1) + 
  theme_classic() +
  xlab("tree height (Ma)") +
  ylab("N per time (true)") +
  scale_x_continuous(breaks = seq(30, 100, length.out = 8)) +
  geom_hline(color = "#fc2861",yintercept = 0.0008, linetype = "dashed") +
  scale_y_log10(labels = scales::comma)

ggsave("figures/true_shift_vs_height.pdf", p, width = 120, height = 70, units = "mm")

#write.csv(df, "output/age_scaling_effect_true_nshift.csv")


