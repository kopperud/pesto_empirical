library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)
library(readr)


setwd("~/projects/pesto_empirical/")

#phy <- read.beast.newick("")





readNumberOfShifts <- function(name, subdir = "empirical"){
  fpath <- paste0("output/", subdir, "/jld2/", name, ".jld2")
  mu <- h5read(fpath, "mu")
  lambda <- h5read(fpath, "lambda")
  lambdaml <- h5read(fpath, "lambdaml")
  muml <- h5read(fpath, "muml")
  etaml <- h5read(fpath, "etaml")
  Ns <- h5read(fpath, "Nsum")
  N <- h5read(fpath, "N")
  
  
  Ntotal <- sum(N)
  
  fpath <- paste0("output/", subdir, "/rates/", name, ".csv")
  df10 <- read.csv(fpath) |> as_tibble()
  df11 <- df10 %>%
    filter(shift_bf > 10, nshift > 0.5)
  N_supported <- N[df11$edge,,,drop=FALSE]
  
  phypath <- paste0("output/", subdir, "/newick/", name, ".tre")
  tree <- read.beast.newick(phypath)
  edge_df <- tree@data
  edge_df <- arrange(edge_df, edge)
  
  bls <- tree@phylo$edge.length
  tree_netdiv <- sum(edge_df$mean_netdiv * bls) / sum(bls)
  tree_lambda <- sum(edge_df$mean_lambda * bls) / sum(bls)
  tree_mu <- sum(edge_df$mean_mu * bls) / sum(bls)
  tree_relext <- sum(edge_df$mean_relext * bls) / sum(bls)
  
  #phy <- read.tree(phypath)
  phy <- tree@phylo
  height <- max(node.depth.edgelength(phy))
  tl <- sum(phy$edge.length)
  
  df1 <- tibble(
    "name" = name,
    "height" = height,
    "N_total" = Ntotal,
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "tree_netdiv" = tree_netdiv,
    "tree_lambda" = tree_lambda,
    "tree_mu" = tree_mu,
    "tree_relext" = tree_relext,
    "inference" = subdir,
    "type" = "pooled",
    "how_many_supported" = nrow(df11)
  )
  
  df2 <- tibble(
    "name" = name,
    "height" = height,
    "N_total" = sum(N_supported),
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "tree_netdiv" = tree_netdiv,
    "tree_lambda" = tree_lambda,
    "tree_mu" = tree_mu,
    "tree_relext" = tree_relext,
    "inference" = subdir,
    "type" = "strong support",
    "how_many_supported" = nrow(df11)
  )
  
  df <- bind_rows(df1, df2)
  return(df)
}




fpaths1 <- Sys.glob("output/empirical/jld2/*.jld2")
#fpath <- fpaths[1]

dfs1 <- list()
ix <- length(fpaths1)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths1[i]
  name <- strsplit(basename(fpath), "\\.")[[1]][[1]]
  dfs1[[i]] <- readNumberOfShifts(name, "empirical")
  setTxtProgressBar(pb,i)
};close(pb)

fpaths2 <- Sys.glob("output/empirical_fixedprior/jld2/*.jld2")

dfs2 <- list()
ix <- length(fpaths2)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths2[i]
  name <- strsplit(basename(fpath), "\\.")[[1]][[1]]
  dfs2[[i]] <- readNumberOfShifts(name, "empirical_fixedprior")
  setTxtProgressBar(pb,i)
};close(pb)

fpaths3 <- Sys.glob("output/empirical_joint/jld2/*.jld2")

dfs3 <- list()
ix <- length(fpaths3)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths3[i]
  name <- strsplit(basename(fpath), "\\.")[[1]][[1]]
  dfs3[[i]] <- readNumberOfShifts(name, "empirical_joint")
  setTxtProgressBar(pb,i)
};close(pb)

df <- bind_rows(
  bind_rows(dfs1),
  bind_rows(dfs2),
  bind_rows(dfs3)
)

###############
##
## some quick plots
##
#############

#df <- bind_rows(dfs)
df[["N_per_time"]] <- df$N_total / df$treelength
df[["support_per_time"]] <- df$how_many_supported / df$treelength
write.csv(df, "output/empirical_munged.csv")

###

df_empirical_bayes <- df %>% 
  filter(inference == "empirical_fixedprior") %>%
  #filter(inference == "empirical") %>%
  filter(type == "pooled")

df_empirical_bayes %>%
  select(c(name, tree_lambda, tree_mu)) %>%
  mutate(r = tree_lambda - tree_mu) %>%
  arrange(-tree_mu) %>%
  print(n = 44)

df_empirical_bayes$lambaml

df_empirical_bayes %>%
  select(c(name, lambdaml, muml)) %>%
  mutate(r = lambdaml - muml) %>%
  print(n = 44)

#df %>%
#  filter(inference == "empirical_fixedprior") %>%
  

######################
##
## calculating which branches had the largest shifts (in netdiv)
##
#######################

fpaths <- Sys.glob("output/empirical_fixedprior/rates/*.csv")
names1 <- gsub("\\.csv", "", basename(fpaths))

rates1 <- lapply(fpaths, read.csv)
for (i in seq_along(rates1)){
  rates1[[i]][["name"]] <- names1[[i]]
}

rates <- bind_rows(rates1)

rates %>% 
  as_tibble() %>%
  filter(shift_bf > 10) %>%
  filter(nshift > 0.5) %>%
  dplyr::arrange(-delta_netdiv) %>%
  print(n = 50)

#tr <- read.beast.newick("output/empirical/newick/Asteraceae_Palazzesi2022.tre")


#keep.tip()

treeheight <- function(phy) max(node.depth.edgelength(phy))

extract.clade(tr@phylo, 3071)

#treeheight(tr@phylo)




library(ggplot2)

p1 <- df %>%
  filter(type == "pooled") %>% 
  ggplot(aes(x = height, y = N_per_time)) +
  geom_point() +
  labs(y = "Shifts per time (N/t)") +
  ggtitle("all branches pooled") +
  theme_classic()

p2 <- df %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = height, y = N_per_time)) +
  geom_point() +
  labs(y = "Shifts per time (N/t)") +
  ggtitle("strongly supported branches") +
  theme_classic()

p3 <- df %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = height, y = how_many_supported)) +
  geom_point() +
  labs(y = "no. supported branches") +
  ggtitle("strongly supported branches") +
  theme_classic()

p4 <- df %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = height, y = support_per_time)) +
  geom_point() +
  labs(y = "no. supported branches per time") +
  ggtitle("strongly supported branches") +
  theme_classic()

p5 <- df %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = height, y = log(support_per_time))) +
  geom_point() +
  labs(y = "log(no. supported branches per time)") +
  ggtitle("strongly supported branches") +
  theme_classic()


px <- (p1 | p2 | p3 | p4 | p5) &
  xlab("tree height (Ma)")

ggsave("figures/123.pdf", px, width = 250, height = 100, units = "mm")


df2 = tibble(
  "type" = df$type,
  "log_height" = log(df$height),
  "log_n_per_time" = log(df$support_per_time)
)

df2 %>%
  dplyr::filter(log_n_per_time > -100) %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = log_height, y = log_n_per_time)) +
  geom_point() + 
  stat_smooth(method = "lm")



df$how_many_supported




###############
## plot some trees


#td <- read.beast.newick("output/empirical_fixedprior/newick/Sigmodontinae_VallejosGarrido2023.tre")
td <- read.beast.newick("output/empirical/newick/Sigmodontinae_VallejosGarrido2023.tre")

td <- read.beast.newick("output/empirical_fixedprior/newick/Mimosa_Vasconcelos2020.tre")
#td <- read.beast.newick("output/empirical/newick/Mimosa_Vasconcelos2020.tre")


ggtree(td, aes(color = mean_netdiv))
ggtree(td, aes(color = delta_netdiv))

td@data$delta_netdiv |> sum()

















