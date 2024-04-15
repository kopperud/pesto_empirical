library(rhdf5)
library(ape)
library(treeio)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)


setwd("~/projects/pesto_empirical/")

list_to_array <- function(l){
  n_branches <- length(l)
  n <- sqrt(length(l[[1]]))
  x <- do.call(c, l)
  array(x, dim = c(n_branches, n, n))
}

shift_fractions <- function(N, lambda, mu, tree_index, label){
  n <- 6
  K <- n^2
  
  #if(length(dim(N)) != 3){
  #  stop("N must have 3 dimensions")
  #}
  
  delta_mu <- matrix(mu, K, K) - t(matrix(mu, K, K))
  delta_lambda <- matrix(lambda, K, K) - t(matrix(lambda, K, K))
  
  eps <- 0.00001
  is_mu <- abs(delta_mu) > eps
  is_lambda <- abs(delta_lambda) > eps
  
  is_only_mu <- is_mu & !is_lambda
  is_only_lambda <- is_lambda & !is_mu
  is_both <- is_mu & is_lambda
  
  #N1 <- colSums(N, dims = 1) ## it's not really a column but I don't know why R names it this way
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
    "N_sum" = N_sum,
    "type" = label,
    "tree_index" = tree_index
  )
  return(res)
}


foo <- function(tree_index){
  #fpath <- paste0("data/simulations/rate_shift_type/", tree_index, ".tre")
  #phy_true <- read.beast.newick(fpath)
  #fpath <- paste0("output/simulations/rate_shift_type/newick/", tree_index, ".tre")
  #phy <- read.beast.newick(fpath)
  
  fpath <- paste0("output/simulations/rate_shift_type/jld2/", tree_index, ".jld2")
  mu <- h5read(fpath, "mu")
  lambda <- h5read(fpath, "lambda")
  lambdaml <- h5read(fpath, "lambdaml")
  muml <- h5read(fpath, "muml")
  etaml <- h5read(fpath, "etaml")
  Ns <- h5read(fpath, "Nsum")
  ntip <- h5read(fpath, "ntip")
  #N <- h5read(fpath, "N")
  
  fpath <- paste0("output/simulations/rate_shift_type/rates/", tree_index, ".csv")
  df <- read.csv(fpath) |> as_tibble()
  
  
  ## true
  #N_true <- list_to_array(phy_true@data$N)
  #df1 <- shift_fractions(N_true, lambda, mu, tree_index, "true")
  
  ## estimated, all
  df2 <- shift_fractions(Ns, lambda, mu, tree_index, "estimate, all pooled")
  df2$ntip <- ntip
  
  ## estimated, only supported branches
  #df5 <- phy@data[order(phy@data$edge),]
  #df5 <- df[order(df$edge),]
  #df5 <- df5 %>%
  #  dplyr::filter(shift_bf > 10)#, nshift > 0.5)
    
  #N_supported <- N[df5$edge,,,drop=FALSE]
  #df3 <- shift_fractions(N_supported, lambda, mu, tree_index, "estimate, strong support")
  
  #df4 <- bind_rows(df1, df2, df3)
  #df4 <- bind_rows(df2, df3)
  return(df2)
}

fpaths <- Sys.glob("~/projects/pesto_empirical/output/simulations/rate_shift_type/newick/*.tre")
valid_tree_indices <- as.numeric(gsub(".tre", "", basename(fpaths)))

dfs <- list()
ix <- 1000
pb <- txtProgressBar(min = 1, max = ix, initial = 1) 
for (i in 1:ix){
  if (i %in% valid_tree_indices){
    dfs[[i]] <- foo(i)  
  }
  setTxtProgressBar(pb,i)
};close(pb)

df <- bind_rows(dfs)

#shift_fractions(N, lambda, mu, tree_index, "mu")
#tree_index <- 3

x <- df %>%
  dplyr::filter(type == "estimate, strong support")
  #filter(N_sum > 0)

x$N_lambda |> mean()
x$N_mu |> mean()
x$N_both |> mean()
plot(x$N_sum, x$N_lambda)

library(tidyr)
df_long <- df %>% pivot_longer(
  cols = c("N_lambda", "N_mu", "N_both", "N_sum", "N_lambda_ratio", "N_mu_ratio", "N_both_ratio"),
  names_to = "Ntype",
  values_to = "N"
)

write.csv(df, "output/rate_shift_type_munged.csv")
write.csv(df_long, "output/rate_shift_type_munged_long.csv")


p1 <- df_long %>%
  filter(type == "estimate, strong support") %>%
  filter(Ntype != "N_sum") %>%
  ggplot(aes(color = Ntype, y = N)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("estimate, strong support")

p2 <- df_long %>%
  filter(type == "estimate, all pooled") %>%
  filter(Ntype != "N_sum") %>%
  ggplot(aes(color = Ntype, y = N)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("estimate, all pooled")

p1 | p2 + plot_layout(guides = "collect")


n <- 6
K <- n^2

prior_ratio <- (n-1)/(K-1)
prior_both_ratio <- ((n-1)^2/(K-1))

p1 <- ggplot(df, aes(color = type, y = N_lambda, x = 1)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = prior_ratio, linetype = "dashed") +
  ggtitle("shift in speciation rate")
  
p2 <- ggplot(df, aes(color = type, y = N_both, x = 1)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = prior_both_ratio, linetype = "dashed") +
  ggtitle("shift in both rates")

p3 <- ggplot(df, aes(color = type, y = N_mu, x = 1)) +
  geom_boxplot() +
  theme_classic() +
  geom_hline(yintercept = prior_ratio, linetype = "dashed") +
  ggtitle("shift in extinction rate")

(p1 | p2 | p3) &
  ylim(c(0.0, 1.0))


p3a <- ggplot(df, aes(y = N_mu)) +#, x = type)) +
  geom_histogram() +
  theme_classic() +
  geom_hline(yintercept = 4/24, linetype = "dashed") +
  ggtitle("shift in extinction rate")
p3a

ggplot(df, aes(x = N_both, fill = type)) +
  geom_histogram() +
  facet_wrap(~ type)



df[df$type=="true","N_lambda"]$N_lambda |> mean()
df[df$type=="true","N_mu"]$N_mu |> mean()
df[df$type=="true","N_both"]$N_both |> mean()







## expectation:
print(paste("single rate: ", (n-1)/(n^2-1)))
print(paste("joint shift: ", ((n-1)^2) / (n^2-1)))






























