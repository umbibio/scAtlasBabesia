library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)


source('./util_funcs.R')

S.O.list <- readRDS('../Input/compScBabesia/rds_rev2/S.O.list.ortholog.rds')
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

phenos <- lapply(S.O.list, `[[`, 1)
S.Os <-  lapply(S.O.list, `[[`, 2)

spps <- unlist(lapply(phenos, function(x) x$spp[1]))
print(spps)

## Integrate Samples
ref.ind <- 4
data.ind <- c(1,2,3,4) ## Exclude Microti
ref.ind <- 4 ## Use B. Div human as reference

all.samples.integrated <- processeMergedS.O(S.O.list, file.info, data.ind, ref.in = ref.ind, res = 0.2, SC = F)

all.samples.integrated@meta.data$spp <- factor(all.samples.integrated@meta.data$spp, 
                                               levels = unique(all.samples.integrated@meta.data$spp))

saveRDS(all.samples.integrated, '../Input/compScBabesia/rds_rev2/all.samples.integrated.rds')

#Idents(all.samples.integrated) <- "phase.cond"

p <- DimPlot(all.samples.integrated, reduction = "pca", 
             #group.by = "cells", 
             split.by = 'spp',
             pt.size = 1,
             shape.by='spp',
             label = TRUE, label.size = 4) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)


flip.pc1 <- F
flip.pc2 <- T

all.samples.integrated.flip <- all.samples.integrated
all.samples.integrated.flip[['pca']]@cell.embeddings[,2] <- -all.samples.integrated.flip[['pca']]@cell.embeddings[,2]


p <- DimPlot(all.samples.integrated.flip, reduction = "pca", 
             #group.by = "cells", 
             split.by = 'spp',
             pt.size = 1,
             shape.by='spp',
             label = TRUE, label.size = 4) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)

## Split the integrated S.Os
DefaultAssay(all.samples.integrated.flip) <- "RNA" ## switch back to original RNA assay

S.O.integrated.list <- SplitObject(all.samples.integrated.flip, split.by = 'spp')


## Fit a pseudo-time to each using and align with Bdiv bulk data.
bd.tc.logCPM <- readRDS('../Input/compScBabesia/rds_rev2/bd_sync_tc_logCPM.rds')

plotUpdatedPstimeS.O <- function(L){
  Idents(L$S.O.bd.update) <- 'adj.time.idx'
  p <- DimPlot(L$S.O.bd.update, reduction = "pca", 
               #group.by = "cells", 
               #split.by = 'spp',
               pt.size = 1,
               shape.by='spp',
               label = TRUE, label.size = 6) + NoLegend() + 
    theme(panel.spacing = unit(0.5, "lines")) + 
    theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
    theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
    theme(
      axis.title.x = element_text(size=14, face="bold"),
      axis.title.y = element_text(size=14, face="bold")
    )
  
  return(p)
  
}

getMaximaLocs <- function(t, y){
  
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
    #maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.8))
    
    ## Peaks for entities of interest
    entity.x = maxval$x
    entity.y = maxval$y
    
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  L <- list(entity.x = entity.x, entity.y = entity.y)
  
  return(L)
}

### Fit pseudo-time, align with B div bulk, and manually adjust based on distribution of lag times
print(spps)

L1 <- list()
L2 <- list()


## BBIG
L1$bbig <- fitPseudoTime(S.O.integrated.list[[1]], reverset.time = T)
L2$bbig <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[1]], L1$bbig$cell.cycle.genes.df, 
                                           L1$bbig$sds.data, bd.tc.logCPM)

plot(L2$bbig$den)
print(L2$bbig$lag.time)
p <- plotUpdatedPstimeS.O(L2$bbig)
plot(p)


## Adjusting the start
getMaximaLocs(L2$bbig$den$x, L2$bbig$den$y) 

L2$bbig <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[1]], L1$bbig$cell.cycle.genes.df, 
                                           L1$bbig$sds.data, bd.tc.logCPM, lag.time = 30)

p <- plotUpdatedPstimeS.O(L2$bbig)
plot(p)


## BBOV
L1$bbov <- fitPseudoTime(S.O.integrated.list[[2]], reverset.time = T)
L2$bbov <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[2]], L1$bbov$cell.cycle.genes.df, 
                                           L1$bbov$sds.data, bd.tc.logCPM)

plot(L2$bbov$den)
print(L2$bbov$lag.time)
p <- plotUpdatedPstimeS.O(L2$bbov)
plot(p)

getMaximaLocs(L2$bbov$den$x, L2$bbov$den$y) ## take the first one
L2$bbov <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[2]], L1$bbov$cell.cycle.genes.df, 
                                           L1$bbov$sds.data, bd.tc.logCPM, lag.time = 29)


p <- plotUpdatedPstimeS.O(L2$bbov)

plot(p)

## BDIV_C
L1$bdiv_cow <- fitPseudoTime(S.O.integrated.list[[3]], reverset.time = T)
L2$bdiv_cow <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[3]], L1$bdiv_cow$cell.cycle.genes.df, 
                                               L1$bdiv_cow$sds.data, bd.tc.logCPM)

plot(L2$bdiv_cow$den)
print(L2$bdiv_cow$lag.time)
p <- plotUpdatedPstimeS.O(L2$bdiv_cow)
plot(p)

## Adjusting the start
getMaximaLocs(L2$bdiv_cow$den$x, L2$bdiv_cow$den$y) ## take the first one

L2$bdiv_cow <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[3]], L1$bdiv_cow$cell.cycle.genes.df, 
                                               L1$bdiv_cow$sds.data, bd.tc.logCPM, lag.time = 29)

p <- plotUpdatedPstimeS.O(L2$bdiv_cow)

plot(p)


## BDIV_H
L1$bdiv_human <- fitPseudoTime(S.O.integrated.list[[4]], reverset.time = T)
L2$bdiv_human <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[4]], L1$bdiv_human$cell.cycle.genes.df, 
                                                 L1$bdiv_human$sds.data, bd.tc.logCPM)

plot(L2$bdiv_human$den)
print(L2$bdiv_human$lag.time)

p <- plotUpdatedPstimeS.O(L2$bdiv_human)
plot(p)

## Adjusting the start
getMaximaLocs(L2$bdiv_human$den$x, L2$bdiv_human$den$y) ## take the first one

L2$bdiv_human <- alignBdivPseudoTimeWithBulkSync(S.O.integrated.list[[4]], L1$bdiv_human$cell.cycle.genes.df, 
                                                 L1$bdiv_human$sds.data, bd.tc.logCPM, lag.time = 28)

p <- plotUpdatedPstimeS.O(L2$bdiv_human)
plot(p)




saveRDS(L1, '../Input/compScBabesia/rds_rev2/all_pstime_fits_L1.rds')
saveRDS(L2, '../Input/compScBabesia/rds_rev2/all_pstime_fits_L2.rds')


cell_cycle_genes_df <- lapply(L1, `[[`, 1)
saveRDS(cell_cycle_genes_df, '../Input/compScBabesia/rds_rev2/all_cell_cycle_genes_df.rds')

cell_cycle_genes_df_adj <- lapply(L2, `[[`, 8)
saveRDS(cell_cycle_genes_df_adj, '../Input/compScBabesia/rds_rev2/all_cell_cycle_genes_df_adj.rds')

sc.tc.df.adj <- lapply(L2, `[[`, 10)
saveRDS(sc.tc.df.adj, '../Input/compScBabesia/rds_rev2/all_sc_tc_df_adj.rds')


#### No Rerun above

titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")
ps <- lapply(1:length(titles), function(i){
  p <- ggplot(L2[[i]]$pc.sds.adj, aes(x=PC_1,y=PC_2)) +
    geom_point(aes(
      fill = cluster.y
    ), shape=21, size = 1.5)+ 
    geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) +
    geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'black', shape=21, size = 5, stroke = 1.2)+
    geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'blue', shape=8, size = 4, stroke = 1.1)+
    theme_bw(base_size = 14) +
    theme(legend.position = "none") +
    ylab('PC2') + xlab('PC1') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(titles[i]) +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) +
    guides(color = 'none')
  
})


p <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol = 4)

plot(p)

ggsave(filename="../Output/compScBabsSpecies/figs/pstime_fits.png",
       plot=p,
       width = 12, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


#LabelClusters(plot = p1, id = 'adj.time.idx', size = 7, repel = T, color = 'black')


## Generate gam genes S.O update an apply network smoothing on common genes
gam.genes <- lapply(L1, function(L){
  gam.genes <- L$gam.genes
})

gam.common <- Reduce(intersect, gam.genes)

S.O.integrated.list.update.orth <- lapply(1:length(L2), function(i){
  S.O <- subset(L2[[i]]$S.O.bd.update, features = gam.common)
  
  return(subset(L2[[i]]$S.O.bd.update, features = gam.common))
})

names(S.O.integrated.list.update.orth) <- names(S.O.list)[1:4]
saveRDS(S.O.integrated.list.update.orth, '../Input/compScBabesia/rds_rev2/S.O.integrated.list.pstime.GAM.intersect.rds')

S.O.integrated.list.update <- lapply(L2, `[[`, 1)

names(S.O.integrated.list.update) <- names(S.O.list)[1:4]
saveRDS(S.O.integrated.list.update, '../Input/compScBabesia/rds_rev2/S.O.integrated.list.pstime.rds')


S.O.integrated.list.update.orth.indiv <- lapply(1:length(L2), function(i){
  S.O <- subset(L2[[i]]$S.O.bd.update, features = gam.genes[[i]])
  
  return(S.O)
})

names(S.O.integrated.list.update.orth.indiv) <- names(S.O.list)[1:4]
saveRDS(S.O.integrated.list.update.orth.indiv, '../Input/compScBabesia/rds_rev2/S.O.integrated.list.pstime.GAM.indiv.rds')

# S.O.gams <- lapply(S.O.integrated.list.update, function(S.O){
#   S.O.gam <- subset(S.O, features = gam.common)
#   #S.O.gam <- prep_S.O(S.O.gam)
#   #S.O.gam.smooth <- smooth.S.O(S.O.gam)
#   return(S.O.gam)
# })

#S.O.gam = lapply(S.O.gams, `[[`, 1)
#S.O.gam.smooth = lapply(S.O.gams, `[[`, 2)

#saveRDS(S.O.gams, '../Input/compScBabesia/RData/S.O.list.ortholog.gam.Rdata')
#saveRDS(S.O.gam.smooth, '../Input/compScBabesia/RData/S.O.list.ortholog.gam.smooth.Rdata')


######
Assays(S.O.integrated.list.update.orth.indiv[[1]])


my.gene <- 'Bdiv-015780c' ## ASF

DefaultAssay(S.O.integrated.list.update.orth.indiv[[2]]) <- 'RNA'
p2 <- FeaturePlot(object = S.O.integrated.list.update.orth.indiv[[2]], 
                  #shape.by = 'spp',
                  #split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = my.gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)

