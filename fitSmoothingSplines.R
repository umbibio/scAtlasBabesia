library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


## Read in the data.

L1 <- readRDS('./rds/all_pstime_fits_L1.rds')
L2 <- readRDS('./rds/all_pstime_fits_L2.rds')

## Run once or read from disk
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

sync.tc.df <- L2$bdiv_human$sync.tc.df
sc.tc.df.adj <- list(bbig = L2$bbig$sc.tc.df.adj, bbov = L2$bbov$sc.tc.df.adj,
                     bdiv_cow = L2$bdiv_cow$sc.tc.df.adj, bdiv_human = L2$bdiv_human$sc.tc.df.adj)


sc.spline.fits <- lapply(sc.tc.df.adj, function(df){
  genes <- unique(df$variable)
  
  spline.fits <- mclapply(1:length(genes), function(i){
    tmp <- df %>% dplyr::filter(variable == genes[i]) %>%
      transmute(GeneID = variable, x = tme, y = y)
    y <- tmp$y
    t <- tmp$x
    
    w <- rep(1, length(y))
    w[which(y == 0)] <- 1/2
    rna.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
    
    rna.sp <- predict(rna.sp, seq(0, 12, by = 1/3)) 
    plot(tmp$x, tmp$y)
    points(rna.sp$x, rna.sp$y, type = 'l', col = 'red')
    #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'blue')
    #points(sc.rna.pp$x, sc.rna.sp$y, type = 'l', col = 'green')
    mu <- data.frame(x = rna.sp$x, y = rna.sp$y) 
    mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
    return(mu)
  }, mc.cores = num.cores)
  
  return(spline.fits)
})

saveRDS(sc.spline.fits, './Input/bd_spline_fits_sc_tc_20min.rds')

## Sychrnonized bulk reads
genes <- unique(sync.tc.df$variable)
sync.spline.fits <- mclapply(1:length(genes), function(i){
  tmp <- sync.tc.df %>% dplyr::filter(variable == genes[i]) %>%
    transmute(GeneID = variable, x = tme, y = y)
  y <- tmp$y
  t <- tmp$x
  
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/2
  rna.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  
  rna.sp <- predict(rna.sp, seq(0, 12, by = 1/3)) 
  plot(tmp$x, tmp$y)
  points(rna.sp$x, rna.sp$y, type = 'l', col = 'red')
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'blue')
  #points(sc.rna.pp$x, sc.rna.sp$y, type = 'l', col = 'green')
  mu <- data.frame(x = rna.sp$x, y = rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

saveRDS(sync,spline.fits, './Input/all_spline_fits_sc_tc_20min.rds')


