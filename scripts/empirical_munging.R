library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)
library(readr)


setwd("~/projects/pesto_empirical/")


readNumberOfShifts <- function(name){
  fpath <- paste0("output/empirical/jld2/", name, ".jld2")
  mu <- h5read(fpath, "mu")
  lambda <- h5read(fpath, "lambda")
  lambdaml <- h5read(fpath, "lambdaml")
  muml <- h5read(fpath, "muml")
  etaml <- h5read(fpath, "etaml")
  Ns <- h5read(fpath, "Nsum")
  N <- h5read(fpath, "N")
  
  
  Ntotal <- sum(N)
  
  fpath <- paste0("output/empirical/rates/", name, ".csv")
  df10 <- read.csv(fpath) |> as_tibble()
  df11 <- df10 %>%
    filter(shift_bf > 10, nshift > 0.5)
  N_supported <- N[df11$edge,,,drop=FALSE]
  
  phypath <- paste0("output/empirical/newick/", name, ".tre")
  phy <- read.tree(phypath)
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
    "type" = "strong support",
    "how_many_supported" = nrow(df11)
  )
  
  df <- bind_rows(df1, df2)
  return(df)
}



fpaths <- Sys.glob("output/empirical/jld2/*.jld2")
#fpath <- fpaths[1]

dfs <- list()
ix <- length(fpaths)
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  fpath <- fpaths[i]
  name <- strsplit(basename(fpath), "\\.")[[1]][[1]]
  dfs[[i]] <- readNumberOfShifts(name)
  setTxtProgressBar(pb,i)
};close(pb)


df <- bind_rows(dfs)
df[["N_per_time"]] <- df$N_total / df$treelength
df[["support_per_time"]] <- df$how_many_supported / df$treelength

write.csv(df, "output/empirical_munged.csv")

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



















