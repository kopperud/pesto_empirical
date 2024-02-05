fpaths <- Sys.glob("data/simulations/rate_shift_type/*.tre")

library(treeio)


lambda <- c(0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767, 0.13172087957400913, 0.19233923491606472, 0.25, 0.3249467017339156, 0.47448817683367767)
mu <- c(0.08781391971600606, 0.08781391971600606, 0.08781391971600606, 0.08781391971600606, 0.08781391971600606, 0.12822615661070985, 0.12822615661070985, 0.12822615661070985, 0.12822615661070985, 0.12822615661070985, 0.16666666666666669, 0.16666666666666669, 0.16666666666666669, 0.16666666666666669, 0.16666666666666669, 0.216631134489277, 0.216631134489277, 0.216631134489277, 0.216631134489277, 0.216631134489277, 0.31632545122245176, 0.31632545122245176, 0.31632545122245176, 0.31632545122245176, 0.31632545122245176)

calculate_magnitude <- function(N, lambda, mu){
  r <- lambda - mu
  mag <- 0.0
  for (i in 1:nrow(N)){
    for (j in 1:ncol(N)){
      if (Ns[i,j] != 0){
        m <- Ns[i,j]
        mag <- mag + m * Ns[i,j] * (r[i] - r[j])
      }
    }
  }
  return(mag) 
}

mags <- numeric(length(fpaths))

n <- length(fpaths)
j <- 0L
pb <- txtProgressBar(min = 1, max = n, initial = 1) 
for (i in 1:n){
  j <- j + 1L
  fpath <- fpaths[i]
  tr <- read.beast.newick(fpath)
  N <- do.call(rbind, tr@data$N)
  Ns <- matrix(colSums(N), nrow = 25) ## it should be by column
  
  Ntotal <- sum(Ns)
  
  if (Ntotal > 0){
    mag <- calculate_magnitude(Ns, lambda, mu)
    mags[j] <- mag / Ntotal 
  }
  
  setTxtProgressBar(pb,i)
}
mags <- mags[1:j]


hist(mags)

#hist(r, breaks = 25)

rm <- matrix(r, 25, 25)
delta_r <- t(rm) - rm

dim(rm)
rm

delta_r

hist(unique(-delta_r))

hist(delta_r)
#delta_r[abs(delta_r)]
vdelta_r <- c(delta_r)


vdelta_r_nontrivial <- vdelta_r[!abs(vdelta_r) < 0.00001]



hist(r, breaks =25)
hist(vdelta_r_nontrivial, breaks = 16)
hist(vdelta_r_nontrivial, breaks = seq(-0.7,0.7, length.out = 15))

hist(r)



lambda




Ns



