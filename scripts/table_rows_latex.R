library(tibble)
library(dplyr)

df_analyses <- read.csv("output/empirical_munged.csv") |> 
    as_tibble() |>
    filter(inference == "empirical_joint") |>
    filter(type == "pooled")

df_metadata <- read.csv("data/empirical/Phylogenies for DeepILS project.csv") |> 
    as_tibble()

names1 <- gsub("\\.tree", "", df_metadata$Filename)
names1 <- gsub("\\.tre", "", names1)
df_metadata$name <- names1

df <- inner_join(df_analyses, df_metadata, by = "name") |>
    arrange(Clade)

for (i in 1:nrow(df)){
    cat(df$Clade[i])
    cat("  \t &  ")
    cat(df$Taxonomic.rank[i])
    cat("  \t &  ")
    cat(format(df$Root.age[i], digits = 1))
    cat("  \t &  ")
    cat(format(df$NTips[i], digits = 1))
    cat("  \t &  ")
    cat(format(df$P.extant.sampling[i], digits = 2, nsmall = 2))
    cat("  \t &  ")
    cat(format(df$tree_netdiv[i], digits = 3, nsmall = 3))
    cat("  \t &  ")
    cat(format(df$N_per_time[i], digits = 3, nsmall = 3))
    cat("  \t &  ")
    cite_label <- paste0("\\cite{", df$bibtex_key[i], "}")
    cat(cite_label)
    cat("  \t \\\\ ")

    cat("\n")
}
