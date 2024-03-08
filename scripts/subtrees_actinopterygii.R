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

max(node.depth.edgelength(fishtree))


extract.clade









