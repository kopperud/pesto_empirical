library(tibble)
library(dplyr)
library(tibble)


df <- read.csv("output/empirical_munged.csv") %>%
    as_tibble()

df <- df %>%
    filter(type == "pooled") %>%
    filter(inference == "empirical_joint") 

m0 <- lm(log(N_per_time) ~ log(tree_netdiv), data = df)
m1 <- lm(log(N_per_time) ~ log(height), data = df)
m2 <- lm(log(N_per_time) ~ log(height) + log(tree_netdiv), data = df)

models <- list(m0, m1, m2)
## print R squared 

for (model in models){
    cat("\n\n")
    s <- summary(model)
    rsq <- s$r.squared
    cat("MODEL CALL: ")
    print(model$call)
    cat(paste("R squared: ", format(rsq, digits = 3)))
    aic1 <- AIC(model)
    cat(paste(",   AIC: ", format(aic1, digits = 3)))
}
cat("\n")
## print AIC

## print summaries

