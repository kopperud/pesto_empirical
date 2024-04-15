#07_rate_shift_type_true_shifts.R

library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)

setwd("~/projects/pesto_empirical/")


shift_fractions <- function(N, lambda, mu, tree_index){
  n <- 5
  K <- n^2
  # if (length(dim(N)) != 3) {
  #   stop("N must have 3 dimensions")
  # }
  
  delta_mu <- matrix(mu, K, K) - t(matrix(mu, K, K))
  delta_lambda <- matrix(lambda, K, K) - t(matrix(lambda, K, K))
  
  N_delta_mu <- sum(delta_mu * N)
  N_delta_lambda <- sum(delta_lambda * N)
  
  eps <- 0.00001
  is_mu <- abs(delta_mu) > eps
  is_lambda <- abs(delta_lambda) > eps
  
  is_only_mu <- is_mu & !is_lambda
  is_only_lambda <- is_lambda & !is_mu
  is_both <- is_mu & is_lambda
  
  #N1 <- colSums(N, dims = 1) ## it's not really a column but I don't know why R names it this way
  N1 <- N
  
  N_sum <- sum(N)
  
  N_mu <- (N[is_only_mu] |> sum())
  N_lambda <- (N[is_only_lambda] |> sum())
  N_both <- (N[is_both] |> sum())
  
  
  N_mu_ratio <- N_mu / N_sum
  N_lambda_ratio <- N_lambda / N_sum
  N_both_ratio <- N_both / N_sum
  
  res <- tibble(
    "N_lambda" = N_lambda,
    "N_mu" = N_mu,
    "N_both" = N_both,
    "N_lambda_ratio" = N_lambda_ratio,
    "N_mu_ratio" = N_mu_ratio,
    "N_both_ratio" = N_both_ratio,
    "N_delta_mu" = N_delta_mu,
    "N_delta_lambda" = N_delta_lambda,
    #"N_delta_both" = N_delta_both,
    "N_sum" = N_sum,
    "tree_index" = tree_index
  )
  return(res)
}



plot.phylo(tree@phylo)

library(microbenchmark)


lambda_true <- c(0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767)
mu_true <- c(0.08781391971600606, 0.08781391971600606, 0.08781391971600606, 0.08781391971600606, 0.08781391971600606, 0.12822615661070985, 0.12822615661070985, 0.12822615661070985, 0.12822615661070985, 0.12822615661070985, 0.16666666666666669, 0.16666666666666669, 0.16666666666666669, 0.16666666666666669, 0.16666666666666669, 0.216631134489277, 0.216631134489277, 0.216631134489277, 0.216631134489277, 0.216631134489277, 0.31632545122245176, 0.31632545122245176, 0.31632545122245176, 0.31632545122245176, 0.31632545122245176)

ix <- 1000
dfs <- list()
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (tree_index in 1:ix){
  tree <- read.beast.newick(paste0("data/simulations/rate_shift_type/", tree_index, ".tre"))
  N <- do.call(cbind, tree@data$N) |> rowSums()
  df <- shift_fractions(N, lambda_true, mu_true, tree_index)
  dfs[[tree_index]] <- df
  setTxtProgressBar(pb, tree_index)
};close(pb)
dfs

df <- bind_rows(dfs)
write.csv(df, "output/rate_shift_type_true_shifts.csv")


plot(1:ix, df$N_lambda, xlim = c(0, 600))
points(201:(200+ix), df$N_mu, add = T)
points(401:(400+ix), df$N_both, add = T)


mean(df$N_lambda)
mean(df$N_mu)
mean(df$N_both)

par(mfrow = c(3,1))
hist(df$N_lambda, xlab = "", main = "Nlambda/Nall", col = "blue")
hist(df$N_mu, xlab = "", main = "Nmu/Nall", col = "red")
hist(df$N_both, xlab = "", main = "Nboth/Nall", col = "green")



sqrt(625)











