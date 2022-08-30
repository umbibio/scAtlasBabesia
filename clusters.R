library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(dtwclust)

source('./util_funcs.R')

getCurvePeakLoc <- function(t, y){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.9))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}

# read fitted splines data with their peak order 

sc.tc.df.adj <- readRDS('../Input/compScBabesia/RData_new/all_sc_tc_df_adj.RData')
sc.tc.fits <- readRDS('../Input/compScBabesia/RData_new/all_sme_fits_sc_tc_20min.RData')

# keeping the genes that are cell cycle markers
markers.sig <- readRDS('../Input/compScBabesia/RData/all.markers.sig.cell.cycle.phase.RData')

sc.tc.mus <- lapply(1:length(sc.tc.fits), function(i){
  sc.tc.mu <- smoothSplineSmeFits(sc.tc.fits[[i]], unique(sc.tc.df.adj[[i]]$variable), extend = F)
  colnames(sc.tc.mu) <- c('GeneID', 't', 'y')
  sc.tc.mu <- sc.tc.mu %>% dplyr::filter(GeneID %in% markers.sig[[i]]$markers.sig$GeneID) %>%
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
    as.data.frame()
  
  return(sc.tc.mu)
})

names(sc.tc.mus) <- names(sc.tc.fits)


# read time course Bdiv Bulk rna from Brendan et al.
sync.tc.df <- readRDS('../Input/compScBabesia/RData_new/bd_sync_tc_df.RData')
sync.tc.fit <- readRDS('../Input/compScBabesia/RData_new/bd_sme_fits_sync_tc_20min.RData')


sync.tc.mu <- smoothSplineSmeFits(sync.tc.fit, unique(sync.tc.df$variable), extend = F)
colnames(sync.tc.mu) <- c('GeneID', 't', 'y')
sync.tc.mu <- sync.tc.mu %>% dplyr::filter(GeneID %in% markers.sig[[4]]$markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
  as.data.frame()


## Generate the clusters with dtw
num.clust <- 5L

sc.hc_dtws <- lapply(sc.tc.mus, function(x){
  sc.hc_dtw <- dtwClustCurves(x[,2:ncol(x)], nclust = num.clust)
  return(sc.hc_dtw)
})


# toxo top marker info 

toxo.markers <- readRDS('../Input/compScBabesia/RData/TG.markers.sig.RData')
toxo.bdiv.orth <- read.xlsx('../Input/compScBabesia/Orthologs/GT1_BDiv.xlsx')
toxo.markers <- inner_join(toxo.markers, toxo.bdiv.orth, by = c('GeneID' = 'GeneID_Toxo'))

# look at the top 20 toxo markers
top.toxo.markers <- toxo.markers %>% dplyr::filter(avg_log2FC > 1) %>% 
  group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 20)

### Extract the top markers in Each babesia spp. Filter for ones that appear in Toxo markers

top.markers <- lapply(1:length(markers.sig), function(i){
  tmp <- markers.sig[[i]]$markers.sig %>% dplyr::filter(avg_log2FC > 1) %>% 
    group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 100) %>%
    ungroup() %>% transmute(GeneID = GeneID, spp = names(markers.sig)[i])
  tmp <- inner_join(tmp, top.toxo.markers, by = c('GeneID' = 'GeneID.y'))
  tmp <- tmp %>% transmute(GeneID = GeneID, spp = spp, Phase = cluster, ProductDescription = ProductDescriptionToxo)
})
names(top.markers) <- names(markers.sig)


## transition phases identified by vision

B.big.ph.bnd <- c(0, 2.5, 6.75, 8.25, 11.35, 12)
B.bov.ph.bond <- c(0, 2.25, 7.5 , 8.25, 11.35, 12)
B.div.cow.ph.bond <- c(0, 2.75, 7.30, 8.75, 11.35, 12)
B.div.hum.ph.bond <- c(0, 2.25, 6.75, 8.75, 11.35 ,12)

ph.trans <- list(B.big.ph.bnd, B.bov.ph.bond, B.div.cow.ph.bond, B.div.hum.ph.bond)
names(ph.trans) <- names(sc.tc.fits)

#  reorder and scale genes for heatmap

sc.tc.mus.scale <- lapply(1:length(sc.tc.mus), function(i){
  
  sc.tc.mus.scale <- sc.tc.mus[[i]]
  
  sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                     center = T,scale = T)
  sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
    pivot_longer(cols = -t, names_to = 'GeneID', values_to = 'y')
  
  
  ## Add curve cluster info
  sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                             order = as.numeric(sc.hc_dtws[[i]]$order),
                             cluster = cutree(sc.hc_dtws[[i]],k = num.clust))
  
  sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')
  
  ## Reorder the genes within each cluster.
  sc.peak.order <- sc.tc.mus.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(t, y))
  sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.peak.order, by = 'GeneID')
  
  
  
  tmp <- sc.tc.mus.scale %>% ungroup() %>% group_by(cluster) %>% mutate(mean.peak = mean(peak.ord))
  tmp <- data.frame(cluster = unique(tmp$cluster[order(tmp$mean.peak)]), cluster.order = 1:length(unique(tmp$cluster)))
  tmp <- sc.tc.mus.scale %>% ungroup() %>% group_by(cluster) %>% mutate(mean.peak = mean(peak.ord)) %>%
    transmute(mean.peak = mean.peak) %>% distinct()
  tmp2 <- data.frame(cluster = unique(tmp$cluster[order(tmp$mean.peak)]), cluster.order = 1:length(unique(tmp$cluster)))
  tmp <- left_join(tmp, tmp2, by = 'cluster')
  
  
  sc.tc.mus.scale <- left_join(sc.tc.mus.scale, tmp, by = 'cluster')
  sc.tc.mus.scale$stage <- paste('dtw', sc.tc.mus.scale$cluster.order, sep = '')
  
  
  return(sc.tc.mus.scale)
})

names(sc.tc.mus.scale) <- names(sc.tc.fits)


# add  toxo phases 
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")

sc.tc.mus.scale.phase <- lapply(1:length(titles), function(i){
  
  sc.tc.mus.scale <- sc.tc.mus.scale[[i]]
  ph <- ph.trans[[i]]

  # defining cell cycle phase based on transition pooints
  sc.tc.mus.scale <- sc.tc.mus.scale %>%
    dplyr::mutate(cell.cycle.phase = ifelse((t >=  ph[1] & t < ph[2]), 'G1.L',
                                            ifelse((t >= ph[2] & t < ph[3]), 'SM',
                                                   ifelse(t >= ph[3] & t < ph[4] , 'MC',
                                                          ifelse(t >= ph[4] & t < ph[5], 'C', 'G1.E')))))
  
  # add transition point to data according to the peak time 
  sc.tc.mus.scale <- sc.tc.mus.scale %>%
    mutate(trans.time = ifelse((t >=  ph[1] & t < ph[2]), ph[2],
                               ifelse((t >= ph[2] & t < ph[3]), ph[3],
                                      ifelse(t >= ph[3] & t < ph[4] , ph[4],
                                             ifelse(t >= ph[4] & t < ph[5], ph[5], ph[6])))))

  
  sc.tc.mus.scale$cell.cycle.phase <- factor(sc.tc.mus.scale$cell.cycle.phase, levels = c('G1.L',  'SM', 'MC', 'C','G1.E'))
  
  return(sc.tc.mus.scale)
})

names(sc.tc.mus.scale.phase) <- names(sc.tc.mus.scale)

## add Marker phase info 
sc.tc.mus.scale.df <- lapply(1:length(sc.tc.mus.scale.phase), function(i){
  
  markers <- markers.sig[[i]]$markers.sig
  names(markers) <- gsub('cluster', 'marker.phase', names(markers))
  markers <- markers %>% group_by(marker.phase) %>% mutate(num.degs  = n())
  
  sc.tc.mus.tab <- sc.tc.mus.scale.phase[[i]]
  
  sc.tc.mus.scale.markers <- left_join(sc.tc.mus.tab, markers, by = 'GeneID')
  sc.tc.mus.scale.markers$cell.cycle.phase <-  factor( sc.tc.mus.scale.markers$cell.cycle.phase ,
                                                       levels = c('G1.L', 'SM', 'MC', 'C', 'G1.E'))
  sc.tc.mus.scale.markers$marker.phase <- factor( sc.tc.mus.scale.markers$marker.phase, 
                                                  levels = c('G', 'SM', 'MC', 'C'))
  
  sc.tc.mus.scale.markers$num.deg.phase <- paste(sc.tc.mus.scale.markers$num.degs, sc.tc.mus.scale.markers$marker.phase, sep = " ")
  
  sc.tc.mus.scale.markers <- sc.tc.mus.scale.markers %>%
    mutate(num.deg.phase = factor(num.deg.phase, unique(arrange(sc.tc.mus.scale.markers, marker.phase)$num.deg.phase)))
  
  return(sc.tc.mus.scale.markers)
  
})

names(sc.tc.mus.scale.df) <- names(sc.tc.mus.scale)



saveRDS(sc.tc.mus.scale.df,'./Input/all_sc_tc_mus_scale_toxo_inferred_cell_cycle_phases_Marker_phases_Progression_heatmap.RData')
sc.tc.mus.scale.df <- readRDS('./Input/all_sc_tc_mus_scale_toxo_inferred_cell_cycle_phases_Marker_phases_Progression_heatmap.RData')

# plot

titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")
title.phases <- c("G1.L", "SM", "MC", "C", "G1.E")
names(title.phases) <- c('G1.L',  'SM', 'MC', 'C','G1.E')

ps <- lapply(1:length(titles), function(i){
  
  p <- ggplot(sc.tc.mus.scale.df[[i]], 
              aes(x = t, reorder_within(GeneID, -peak.ord, stage ), fill = y)) +
    
    geom_tile() +
    scale_x_discrete(expand=c(0,0)) +
    ylab("Genes") + xlab("time/cells") +
    scale_fill_gradientn(colours = viridis::inferno(10)) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold"),
      #axis.title.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(),
      legend.position = "none") +
    facet_grid(num.deg.phase ~ cell.cycle.phase  , scales = 'free', space = 'free',
               labeller = labeller(cell.cycle.phase = title.phases)) +
    theme(strip.background=element_rect(fill='white', color = 'black'))+
    theme(panel.spacing = unit(0.1, "lines")) +
    theme(strip.text = element_text(size = 14, face = 'bold'))+
    ggtitle(titles[i])+
    theme(
      plot.title = element_text(size=18, face = "bold.italic", color = 'black', hjust = 0.5),
      axis.title.x = element_text(size=14, face="bold"),
    )
  
})
p <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol = 4)

