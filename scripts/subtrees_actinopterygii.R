library(ape)

fishtree <- read.tree("data/empirical/Actinopterygii_Rabosky2018.tree")


fishtree

ntip <- length(fishtree$tip.label)
internal_node_indices <- (ntip+2):max(fishtree$edge)


set.seed(123)
random_nodes <- sample(internal_node_indices, 1000)
random_nodes <- internal_node_indices

#n_subtrees <- 10000




trees <- list()
i <- 1L
used_nodes <- numeric()
for (node in random_nodes){
  tree <- extract.clade(fishtree, node)
  if (length(tree$tip.label) > 30){
    trees[[i]] <- tree
    used_nodes <- append(used_nodes, node)
    i <- i + 1L
  }
  if (i > 12000){
    break
  }
}


for (j in seq_along(trees)){
  fpath <- paste0("data/empirical_fish_subtrees/Actinopterygii_node_", used_nodes[j], ".tre")
  write.tree(trees[[j]], fpath)
}


heights <- sapply(trees, function(tree) max(node.depth.edgelength(tree)))
ntips <- sapply(trees, function(tree) length(tree$tip.label))

hist(heights)
hist(log(ntips))
hist(ntips)

ntips

plot(heights, ntips)

plot.phylo(fishtree, show.tip.label = FALSE);nodelabels()


extrema <- function(x){
  mi <- min(x)
  ma <- max(x)
  return(c(mi,ma))
}

library(treeio)
library(ggtree)
library(ggplot2)
#tr <- read.beast.newick("output/empirical_fish_subtrees/newick/Actinopterygii_node_22326.tre")
tr <- read.beast.newick("output/empirical_fish_subtrees/newick/Actinopterygii_node_12596.tre")

th <- max(node.depth.edgelength(tr@phylo))
p1 <- ggtree(tr, aes(color = mean_netdiv)) +
  theme(legend.position = "top") +
  ggtitle("net diversification (range 0.54 to 2.51)")
  #scale_color_continuous(low='darkgreen', high='red') +
  #xlim(c(0.0, 10.0))
p2 <- ggtree(tr, aes(color = nshift)) +
  theme(legend.position = "top") +
  scale_color_continuous(low='black', high='red') +
  geom_tiplab(size = 2) +
  xlim(c(0.0, th+2)) +
  ggtitle("number of rate shifts (range 0.002 to 12.5)")


p <- p1 | p2
ggsave("figures/weird_cichlid_clade.pdf", p, height = 200, width = 300, units = "mm")

















