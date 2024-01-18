library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)
library(readr)


setwd("~/projects/pesto_empirical/")


heights <- seq(30, 100, length.out = 8)
heights


fpaths <- Sys.glob("output/simulations/age_scaling_effect/jld2/*.jld2")

fpath <- fpaths[1]

readNumberOfShifts <- function(fpath){
  hi <- strsplit(
    strsplit(fpath, "/")[[1]][5],
    "_")[[1]] |>
    readr::parse_number()
  height <- hi[1]
  tree_index <- hi[2]
  
  
  fpath <- paste0("output/simulations/age_scaling_effect/jld2/h", height, "_", tree_index, ".jld2")
  mu <- h5read(fpath, "mu")
  lambda <- h5read(fpath, "lambda")
  lambdaml <- h5read(fpath, "lambdaml")
  muml <- h5read(fpath, "muml")
  etaml <- h5read(fpath, "etaml")
  Ns <- h5read(fpath, "Nsum")
  N <- h5read(fpath, "N")
  
  
  Ntotal <- sum(N)
  
  fpath <- paste0("output/simulations/age_scaling_effect/rates/h", height, "_", tree_index, ".csv")
  df10 <- read.csv(fpath) |> as_tibble()
  df11 <- df10 %>%
    filter(shift_bf > 10, nshift > 0.5)
  N_supported <- N[df11$edge,,,drop=FALSE]
  
  phypath <- paste0("output/simulations/age_scaling_effect/newick/h", height, "_", tree_index, ".tre")
  phy <- read.tree(phypath)
  tl <- sum(phy$edge.length)
  
  df1 <- tibble(
    "tree_index" = tree_index,
    "height" = height,
    "N_total" = Ntotal,
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "type" = "pooled",
    "how_many_supported" = nrow(df11)
  )
  
  df2 <- tibble(
    "tree_index" = tree_index,
    "height" = height,
    "N_total" = sum(N_supported),
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "type" = "strong support",
    "how_many_supported" = nrow(df11)
  )
  
  df <- bind_rows(df1, df2)
  return(df)
}



dfs <- list()
ix <- length(fpaths)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths[i]
  dfs[[i]] <- readNumberOfShifts(fpath)
  setTxtProgressBar(pb,i)
};close(pb)


df <- bind_rows(dfs)
df[["N_per_time"]] <- df$N_total / df$treelength
df[["support_per_time"]] <- df$how_many_supported / df$treelength

write.csv(df, "output/age_scaling_effect_munged.csv")

library(ggplot2)

p1 <- df %>%
  filter(type == "pooled") %>% 
  ggplot(aes(x = height, y = N_per_time)) +
  geom_point() +
  xlim(c(0.0, 110)) +
  labs(y = "Shifts per time (N/t)") +
  ggtitle("all branches pooled") +
  theme_classic()

p2 <- df %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = height, y = N_per_time)) +
  geom_point() +
  xlim(c(0.0, 110)) +
  labs(y = "Shifts per time (N/t)") +
  ggtitle("strongly supported branches") +
  theme_classic()

p3 <- df %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = height, y = how_many_supported)) +
  geom_point() +
  xlim(c(0.0, 110)) +
  labs(y = "no. supported branches") +
  ggtitle("strongly supported branches") +
  theme_classic()

p4 <- df %>%
  filter(type == "strong support") %>% 
  ggplot(aes(x = height, y = support_per_time)) +
  geom_point() +
  xlim(c(0.0, 110)) +
  labs(y = "no. supported branches per time") +
  ggtitle("strongly supported branches") +
  theme_classic()


(p1 | p2 | p3 | p4) &
  xlab("tree height (Ma)")


df$how_many_supported














  




