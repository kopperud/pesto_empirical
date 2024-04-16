library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)
library(readr)


setwd("~/projects/pesto_empirical/")


readNumberOfShifts <- function(dir = "age_scaling_effect", fpath){
  hi <- strsplit(
    strsplit(fpath, "/")[[1]][5],
    "_")[[1]] |>
    readr::parse_number()
  height <- hi[1]
  tree_index <- hi[2]
  
  
  fpath <- paste0("output/simulations/", dir, "/jld2/h", height, "_", tree_index, ".jld2")
  mu <- h5read(fpath, "mu")
  lambda <- h5read(fpath, "lambda")
  lambdaml <- h5read(fpath, "lambdaml")
  muml <- h5read(fpath, "muml")
  etaml <- h5read(fpath, "etaml")
  #Ns <- h5read(fpath, "Nsum")
  #N <- h5read(fpath, "N")
  
  
  
  
  fpath <- paste0("output/simulations/", dir, "/rates/h", height, "_", tree_index, ".csv")
  df10 <- read.csv(fpath) |> as_tibble()
  df11 <- df10 %>%
    #filter(shift_bf > 10, nshift > 0.5)
    filter(shift_bf > 10.0)
  #N_supported <- N[df11$edge,,,drop=FALSE]
  Ntotal <- sum(df10$nshift)
  N_supported <- df11$nshift
  
  phypath <- paste0("output/simulations/", dir, "/newick/h", height, "_", tree_index, ".tre")
  phy <- read.tree(phypath)
  tl <- sum(phy$edge.length)
  ntip <- length(phy$tip.label)
  
  df1 <- tibble(
    "tree_index" = tree_index,
    "ntip" = ntip,
    "height" = height,
    "N_total" = Ntotal,
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "type" = "pooled",
    "inference" = dir,
    "how_many_supported" = nrow(df11)
  )
  
  df2 <- tibble(
    "tree_index" = tree_index,
    "ntip" = ntip,
    "height" = height,
    "N_total" = sum(N_supported),
    "treelength" = tl,
    "muml" = muml,
    "etaml" = etaml,
    "lambdaml" = lambdaml,
    "type" = "strong support",
    "inference" = dir,
    "how_many_supported" = nrow(df11)
  )
  
  df <- bind_rows(df1, df2)
  return(df)
}

heights <- seq(30, 100, length.out = 8)
heights



fpaths<- Sys.glob("output/simulations/age_scaling_effect/jld2/*.jld2")

dfs1 <- list()
#dirs <- c("age_scaling_effect", "age_scaling_effect_fixedprior")
ix <- length(fpaths)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths[i]
  dfs1[[i]] <- readNumberOfShifts(dir = "age_scaling_effect", fpath)
  setTxtProgressBar(pb,i)
};close(pb)
df1 <- bind_rows(dfs1)

# ###########
# ## load the trees with fixed prior
# fpaths <- Sys.glob("output/simulations/age_scaling_effect_fixedprior/jld2/*.jld2")
# 
# dfs2 <- list()
# ix2 <- length(fpaths)
# pb <- txtProgressBar(min = 1, max = ix2, initial = 1) 
# for (i in 1:ix2){
#   fpath <- fpaths[i]
#   dfs2[[i]] <- readNumberOfShifts(dir = "age_scaling_effect_fixedprior", fpath)
#   setTxtProgressBar(pb,i)
# };close(pb)
# df2 <- bind_rows(dfs2)
# 
# df <- bind_rows(
#   df1, df2
# )
df <- df1
df[["N_per_time"]] <- df$N_total / df$treelength
df[["support_per_time"]] <- df$how_many_supported / df$treelength

write.csv(df, "output/age_scaling_effect_munged.csv")


##############
##
##  some quick plots
##
################

## subsample the empirical hyperprior analyses to match visually
number_of_samples <- nrow(df1[df1$height == 30 & df1$type == "pooled",])
dfx <- df %>%
  filter(type == "pooled") %>%
  filter(inference == "age_scaling_effect") %>%
  #group_by(height) %>%
  #slice_head(n = number_of_samples) %>%
  filter(tree_index < 51)
  #ungroup()

df3 <- bind_rows(
  df1, dfx
)
df3[["N_per_time"]] <- df3$N_total / df3$treelength

p_empirical_vs_fixedprior <- df3 %>%
  filter(type == "pooled") %>% 
  #filter(how_many_supported > 0) %>%
  ggplot(aes(x = height, y = N_per_time, color = inference)) +
  theme_classic() +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  facet_grid(cols = vars(inference))

ggsave("figures/empirical_vs_fixedrates.pdf", p_empirical_vs_fixedprior)

df4 <- df3 %>%
  filter(type == "pooled")

x <- filter(df4, inference == "age_scaling_effect")$etaml
y <- filter(df4, inference == "age_scaling_effect_fixedprior")$etaml

df5 <- tibble(
  "hyperprior" = x,
  "fixedprior" = y
)
ggplot(df5, aes( x = hyperprior,y = fixedprior)) +
  geom_point() +
  theme_classic() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0)
#plot(log(x), log(y))






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














  




