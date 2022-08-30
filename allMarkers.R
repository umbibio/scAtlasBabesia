
## Marker analysis

library(Seurat)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ggVennDiagram)

source('./util_funcs.R')

## helper functions 

# generate matched contrasts between phases (ex. C vs C)
makeMatchedContrasts <- function(S.O.integrated){
  
  objs <- as.character(unique(S.O.integrated@meta.data$phase.spp))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.spp = gsub(":.*", "", ref), ref.phase = gsub(".*:", "", ref), 
           query.spp = gsub(":.*", "", query), query.phase = gsub(".*:", "", query))
  my.contrasts <- contrasts %>% dplyr::filter(ref.phase == query.phase & ref.spp != query.spp)
  
  return(my.contrasts)
  
}

# Find markers between desired contrasts
getSigGlobalMarkers <- function(S.O.obj){
  
  markers <- FindAllMarkers(object = S.O.obj, only.pos = TRUE, min.pct = 0) 
  markers$GeneID <- gsub('-', '_', markers$gene)
  markers.sig <- markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
  
  M <- list(markers.sig  = markers.sig)
  
  return(M)
}

# Find markers between cell cycle phases 
getCellCyclePhaseMarkers <- function(S.O.integrated){
  S.O.integrated.merge.list <- SplitObject(S.O.integrated, split.by = "spp")
  spps <- names(S.O.integrated.merge.list)
  all.markers.list <- mclapply(1:length(spps), function(i){
    S.O <- S.O.integrated.merge.list[[i]]
    DefaultAssay(S.O) <- "RNA"
    Idents(S.O) <- 'cell.cycle.phase'
    markers <- FindAllMarkers(object = S.O, only.pos = TRUE)
    markers$GeneID = gsub('-', '_', markers$gene)
    colnames(markers)[colnames(markers) == 'cluster'] <- 'phase'
    markers$spp <- spps[i]
    return(markers)
  }, mc.cores = num.cores) 
  all.markers <- bind_rows(all.markers.list)
  
  all.markers.sig <- all.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
  
  return(all.markers.sig)
}


# extract expression data for plotting 
exprExprData <- function(S.O, genes){
  expr <- data.frame(as.matrix(S.O@assays$RNA@data))
  expr$GeneID <- gsub('-', '_', rownames(as.matrix(S.O@assays$RNA@data)))
  expr <- expr %>% pivot_longer(-GeneID, names_to = 'Sample', values_to = 'expr')
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), Phase = S.O@meta.data$cell.cycle.phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  
  expr.filt <- expr %>% dplyr::filter(GeneID %in% genes)
  
  expr.filt <- inner_join(meta.data, expr.filt, by = 'Sample')
  return(expr.filt)  
}


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


## read in data 

prod.desc <- read.xlsx('./rds/BDiv_Prod_desc.xlsx')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, ProductDescriptionBdiv = Product.Description) %>% distinct()


toxo.bdiv.orth <- read.xlsx("./rds/GT1_BDiv.xlsx")
colnames(toxo.bdiv.orth) <- c('GeneID_Toxo', 'GeneID', 'ProductDescriptionToxo')

num.cores <- detectCores()


## cell cycle marker analysis done per spp (no common genes)

S.O.integrated.list <- readRDS("./rds/S.O.integrated.list.pstime.GAM.indiv.cell.cycle.phase.rds")
markers.sig.list <- lapply(S.O.integrated.list, function(S.O) {
  
  spp <- unique(S.O@meta.data$spp)
  Idents(S.O) <- "cell.cycle.phase"
  cat(paste('Finding sig global markers of', spp))
  cat('\n')
  
  getSigGlobalMarkers(S.O)
  
})

clusters <-  c('G', 'SM', 'MC', 'C')
spps <- c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human')
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")

for(i in 1:length(titles)){
  markers.sig.list[[i]]$markers.sig$spp <- spps[i]
}

all.markers.sig <- do.call(Map, c(f = rbind, markers.sig.list))
all.markers.sig <- all.markers.sig$markers.sig
all.markers.sig$spp <-  factor(all.markers.sig$spp, levels = c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human'))
names(all.markers.sig) <- gsub("cluster", "phase", names(all.markers.sig))
all.markers.sig$phase <- factor(all.markers.sig$phase, levels = c('G', 'SM', 'MC', 'C'))
all.markers.sig.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(genes = list(GeneID), num.deg = n()) 


#saveRDS(all.markers.sig.stat, './Input/all.markers.sig.cell.cycle.phase.RData')



## Doing diff exp using anchored datasets
## Orthologs S.Os after  integration


S.O.integrated <- readRDS('./rds/S.O.integrated.list.pstime.GAM.indiv.cell.cycle.phase.rds')
spps <- names(S.O.integrated)

## Differential expression must be done on genes present in all datasets
comm.genes <- lapply(S.O.integrated, function(S.O){
  genes <- rownames(S.O@assays$RNA@data)
})

comm.genes <- Reduce(intersect, comm.genes)

S.O.integrated <- lapply(S.O.integrated, function(S.O){
  S.O <- subset(S.O, features = comm.genes)
})
names(S.O.integrated) <- spps

## We nned to re-merge and integrate the data for cross spp marker analysis
S.O.integrated.merge <- merge(S.O.integrated[[1]], S.O.integrated[2:4], 
                              add.cell.ids = names(S.O.integrated))
S.O.integrated.merge.list <- SplitObject(S.O.integrated.merge, split.by = "spp")
S.O.integrated.merge.list <- lapply(X = S.O.integrated.merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = S.O.integrated.merge.list)
anchors <- FindIntegrationAnchors(object.list = S.O.integrated.merge.list, 
                                  anchor.features = features, reference = 4)

S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay for PCA. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
S.O.integrated <- ScaleData(S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, npcs = 30, verbose = FALSE)
S.O.integrated <- RunUMAP(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindNeighbors(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.12)

S.O.integrated@meta.data$spp <- factor(S.O.integrated.merge@meta.data$spp, 
                                       levels = unique(S.O.integrated.merge@meta.data$spp))

S.O.integrated@meta.data$phase.spp <- paste(S.O.integrated@meta.data$spp, S.O.integrated@meta.data$cell.cycle.phase, sep = ':')



S.O.integrated <- readRDS('./rds/S.O.integrated.list.pstime.GAM.indiv.cell.cycle.phase_sinlge_obj.rds')


## Cell cycle marker: Independently done per spp (common genes)
cell.cycle.markers.sig <- getCellCyclePhaseMarkers(S.O.integrated)
cell.cycle.markers.sig <- left_join(cell.cycle.markers.sig, prod.desc, by = 'GeneID') 
cell.cycle.markers.sig <- left_join(cell.cycle.markers.sig, toxo.bdiv.orth, by='GeneID') %>% arrange(spp, desc(avg_log2FC))

# we will use this markers for filtering genes on the following markeer analysis


## Cross spp differential expression: Phase-based unique markers

## spp specific cell cycle markers 
## Identity based comparison. EX: compare C to C accross all spps
## ident.1 case, ident.2 is control

DefaultAssay(S.O.integrated) <- "RNA"
Idents(S.O.integrated) <- "phase.spp"


contrasts <- makeMatchedContrasts(S.O.integrated)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.spp<- gsub(':.*', '', contrasts.groups$ref)
matched.DEGs <- mclapply(1:nrow(contrasts.groups), function(i){
  tmp <- FindMarkers(S.O.integrated, ident.1 = contrasts.groups$ref[i],
                     ident.2 = c(unlist(contrasts.groups$query[i])), 
                     only.pos = T, verbose = T)
  tmp$ref <- contrasts.groups$ref[i]
  tmp$ref.spp <- contrasts.groups$ref.spp[i]
  tmp$phase <- contrasts.groups$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})
markers.int.local <- bind_rows(matched.DEGs)
markers.int.local.sig <- markers.int.local %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
markers.int.local.sig <- markers.int.local.sig  %>% mutate(pct.ratio = (pct.1 / pct.2))

## Now consider markers that are phase specific
markers.int.local.sig.cc <- inner_join(markers.int.local.sig, cell.cycle.markers.sig, 
                                       by = c('GeneID', 'ref.spp' = 'spp', 'phase'))

markers.int.local.sig.specific <- markers.int.local.sig.cc %>% 
  arrange(GeneID, ref.spp, desc(pct.ratio)) %>% group_by(GeneID) %>% slice_max(n = 1, order_by = pct.ratio)
markers.int.local.sig.specific.df <- markers.int.local.sig.specific %>% dplyr::select(!contains(".y"))




## Find conserved markers in matched inferred phases


DefaultAssay(S.O.integrated) <- 'RNA'
Idents(S.O.integrated) <- 'cell.cycle.phase'
phases  <- c('G', 'SM', 'MC', 'C')

cons.markers <- lapply(phases, function(cc){
  cons.marker <- FindConservedMarkers(S.O.integrated, ident.1 = cc, grouping.var = "spp", verbose = T)
  cons.marker$gene <- rownames(cons.marker)
  cons.marker$GeneID <- gsub('-', '_', cons.marker$gene)
  cons.marker$phase <- cc
  return(cons.marker)
})


cons.markers.all.phases <- do.call("rbind", cons.markers)

cons.FC <- cons.markers.all.phases[grepl('GeneID|log2FC|phase', names(cons.markers.all.phases))]
colnames(cons.FC) <- gsub("_avg_log2FC", "", colnames(cons.FC))
cons.FC.lng <- cons.FC %>% pivot_longer(!c(GeneID , phase),names_to = "spp" , values_to = 'max_avg_log2FC')

cons.FC.max <- cons.FC.lng %>%  arrange(GeneID,  desc(max_avg_log2FC)) %>% group_by(GeneID) %>% slice_max(n= 1, order_by= max_avg_log2FC)
cons.FC.max <- cons.FC.max %>% dplyr::select(GeneID, phase,spp, max_avg_log2FC)

cons.markers.phases <- left_join(cons.FC.max, cons.markers.all.phases, by = c('GeneID',"phase")) 
cons.markers.phases.sig <- cons.markers.phases %>% dplyr::filter(max_pval < 0.01 & max_avg_log2FC > 1)
cons.markers.phases.sig.df <- cons.markers.phases.sig %>%
  dplyr::select(GeneID, phase, spp, gene, max_avg_log2FC, minimump_p_val, max_pval )


# ## Now consider markers that are phase specific
cons.markers.phases.sig.cc <- inner_join(cons.markers.phases.sig.df, cell.cycle.markers.sig, 
                                         by = c("gene" ,"GeneID", "phase", "spp"))


##################################################################
# second approach of finding conserved and spp specific markers


## looking at the intersection and differences in each phase across species

markers.sig <- readRDS('./rds/all.markers.sig.cell.cycle.phase.list.rds')
clusters <-  c('G', 'SM', 'MC', 'C')
spps <- c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human')
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")

for(i in 1:length(titles)){
  markers.sig[[i]]$markers.sig$spp <- spps[i]
}

all.markers.sig <- do.call(Map, c(f = rbind, markers.sig))
all.markers.sig <- all.markers.sig$markers.sig
all.markers.sig$spp <-  factor(all.markers.sig$spp, levels = c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human'))
names(all.markers.sig) <- gsub("cluster", "phase", names(all.markers.sig))
all.markers.sig$phase <- factor(all.markers.sig$phase, levels = c('G', 'SM', 'MC', 'C'))
all.markers.sig.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(num.deg = n()) 

all.markers <- all.markers.sig %>% group_by(spp, phase) %>% summarise(genes = list(GeneID), total = n()) %>%
  group_by(phase) %>% 
  mutate(intersect = list(reduce(genes, intersect))) %>% ungroup() %>% mutate(intersect.length = lengths(intersect))

## Identify shared markers in phases across all spp

all.markers.shared.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(genes = list(GeneID), total = n()) %>%
  group_by(phase) %>% 
  summarise(intersect = list(reduce(genes, intersect)), intersect.length = lengths(intersect)) 

## shared markers
all.markers.shared.stat$phase <- factor(all.markers.shared.stat$phase, levels = c('C', 'MC', 'SM', 'G'))


## Venn diagram to show the intersection and differences
markers.sig <- readRDS('./input//all.markers.sig.cell.cycle.phase.list.rds')
clusters <-  c('G', 'SM', 'MC', 'C')
spps <- c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human')
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")

for(i in 1:length(titles)){
  markers.sig[[i]]$markers.sig$spp <- spps[i]
}

getVenn <- function(markers.list, cluster){
  
  Venn.list <- list(
    Bbig = markers.list[[1]]$markers.sig$GeneID[markers.list[[1]]$markers.sig$cluster == cluster] ,
    Bbov = markers.list[[2]]$markers.sig$GeneID[markers.list[[2]]$markers.sig$cluster == cluster] ,
    Bdiv_cow = markers.list[[3]]$markers.sig$GeneID[markers.list[[3]]$markers.sig$cluster == cluster],
    Bdiv_human = markers.list[[4]]$markers.sig$GeneID[markers.list[[4]]$markers.sig$cluster == cluster]
  )
  
}


listforVenn <- lapply(1:length(clusters), function(i){
  
  getVenn(markers.sig, clusters[i])
  
})


saveRDS(listforVenn, "./rds/cell_cycle_markers_shared_across_speciies_Venn_diag.rds")

ps <- lapply(1:length(clusters), function(i){
  
  venn <- Venn(listforVenn[[i]])
  data <- process_data(venn)
  
  p <- ggplot() +
    # 1. region count layer
    geom_sf(aes(fill = count), data = venn_region(data)) +
    # 2. set edge layer
    geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
    # 4. region label layer
    geom_sf_label(aes(label = count), data = venn_region(data),size =10,  alpha = 0.5) +
    geom_sf_text(aes(label = name), color=c("bbig" = "firebrick","bbov" ="darkorchid3",
                                            'bdiv_cow' = 'darkslateblue', 'bdiv_human' = 'darkolivegreen4'),
                 fontface = "bold", size= 6,nudge_y = 0,nudge_x=0,
                 data = venn_setlabel(data))+
    geom_sf(color=c("bbig" = "firebrick","bbov" ="darkorchid3",
                    'bdiv_cow' = 'darkslateblue', 'bdiv_human' = 'darkolivegreen4'), 
            size = 1, data = venn_setedge(data), show.legend = FALSE) + 
    #geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), data = venn_region(data))+
    theme_void()+ 
    theme(legend.position = "none")+
    theme(plot.margin=unit(c(0,0,0,0),"cm"))+
    ggtitle(clusters[i])+
    theme(plot.title = element_text(size = "22", face = "bold.italic"))
  
  
})

p <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol = 4)


####### species specific markers in each phase subtracting the intersect 

all.markers.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(genes = list(GeneID), total = n())

XX <- full_join(all.markers.stat, all.markers.shared.stat, by = "phase")
XX.diff <- XX %>% rowwise() %>% 
  mutate(spp.genes = list(setdiff(genes, intersect)), num.spp.genes = length(unlist(setdiff(genes, intersect))))

spp.specific.markers <- XX.diff %>% select(spp, phase, spp.genes,num.spp.genes)

spp.specific.markers$phase <- factor(spp.specific.markers$phase, levels = c('G', 'SM', 'MC', 'C'))
saveRDS(spp.specific.markers, "./Input/spp.specific.markers_second_approach.rds")


#####
## Identifying host specific markers
#####

##  Bdiv_human vs Bdiv_cow: Phase_based

DefaultAssay(S.O.integrated) <- 'RNA'
Idents(S.O.integrated) <- 'spp'

S.O.int.bdivs <- subset(S.O.integrated, ident = c('Bdiv_human', 'Bdiv_cow'))
DefaultAssay(S.O.int.bdivs) <- 'RNA'
Idents(S.O.int.bdivs) <- 'phase.spp'

contrasts <- makeMatchedContrasts(S.O.integrated)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.spp<- gsub(':.*', '', contrasts.groups$ref)

contrasts.bdivs <- contrasts %>% dplyr::filter(grepl('Bdiv', ref.spp) & grepl('Bdiv', query.spp))
contrasts.bdivs$phase <- gsub('.*:', '', contrasts.bdivs$ref)
contrasts.bdivs$ref.spp<- gsub(':.*', '', contrasts.bdivs$ref)

# comparison of matched phases 
matched.DEGs.bdivs <- mclapply(1:nrow(contrasts.bdivs), function(i){
  tmp <- FindMarkers(S.O.int.bdivs, ident.1 = contrasts.bdivs$ref[i],
                     ident.2 = contrasts.bdivs$query[i], 
                     only.pos = T, verbose = T)
  tmp$ref <- contrasts.bdivs$ref[i]
  tmp$ref.spp <- contrasts.bdivs$ref.spp[i]
  tmp$phase <- contrasts.bdivs$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})
markers.int.local.bdivs <- bind_rows(matched.DEGs.bdivs)
markers.int.local.bdivs.sig <- markers.int.local.bdivs %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.01)
markers.int.local.bdivs.sig <- markers.int.local.bdivs.sig  %>% mutate(pct.ratio = (pct.1 / pct.2))

## Now consider markers that are phase specific
markers.int.local.bdivs.sig.cc <- inner_join(markers.int.local.bdivs.sig, cell.cycle.markers.sig, 
                                             by = c('GeneID', 'ref.spp' = 'spp', 'phase'))
markers.int.local.bdivs.sig.cc <- markers.int.local.bdivs.sig.cc %>% dplyr::select(!contains(".y"))

markers.int.local.bdivs.sig.df <- markers.int.local.bdivs.sig.cc %>% 
  arrange(GeneID, desc(pct.ratio)) %>% group_by(GeneID, ref.spp) %>% slice_max(n = 1, order_by = pct.ratio)


####
###  human vs cow: phase based (Bdiv_human vs merged species in cow host)
####

DefaultAssay(S.O.integrated) <- "RNA"
S.O.integrated@meta.data <- S.O.integrated@meta.data %>% mutate(host = ifelse( spp == "Bdiv_human", "human", "cow"))
S.O.integrated@meta.data$host.phase <- paste(S.O.integrated@meta.data$host, S.O.integrated@meta.data$cell.cycle.phase, sep = ":")

Idents(S.O.integrated) <- 'host.phase'

# comparison of matched phases between Bdiv_human and rest of species
makeHostMatchedContrasts <- function(S.O.integrated){
  
  objs <- as.character(unique(S.O.integrated@meta.data$host.phase))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.host = gsub(":.*", "", ref), ref.phase = gsub(".*:", "", ref), 
           query.host = gsub(":.*", "", query), query.phase = gsub(".*:", "", query))
  my.contrasts <- contrasts %>% dplyr::filter(ref.phase == query.phase & ref.host != query.host)
  
  return(my.contrasts)
  
}

contrasts <- makeHostMatchedContrasts(S.O.integrated)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.host <- gsub(':.*', '', contrasts.groups$ref)

matched.DEGs <- mclapply(1:nrow(contrasts.groups), function(i){
  tmp <- FindMarkers(S.O.integrated, ident.1 = contrasts.groups$ref[i],
                     ident.2 = c(unlist(contrasts.groups$query[i])), 
                     only.pos = T, verbose = T)
  tmp$ref <- contrasts.groups$ref[i]
  tmp$ref.host <- contrasts.groups$ref.host[i]
  tmp$phase <- contrasts.groups$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})


host.markers.int.local <- bind_rows(matched.DEGs)
host.markers.int.local.sig <- host.markers.int.local %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.01)
host.markers.int.local.sig  <- host.markers.int.local.sig %>% mutate(pct.ratio = (pct.1 / pct.2))


cell.cycle.markers.sig <- cell.cycle.markers.sig %>% mutate(host = ifelse( spp == "Bdiv_human", "human", "cow"))
host.markers.int.local.sig.cc <- inner_join(host.markers.int.local.sig, cell.cycle.markers.sig, 
                                            by = c('GeneID', 'phase', 'ref.host' = 'host'))
host.markers.int.local.sig.cc <- host.markers.int.local.sig.cc %>% dplyr::select(!contains(".y"))

# unique to spp and phase
host.markers.int.local.sig.specific <- host.markers.int.local.sig.cc  %>%  distinct(GeneID, .keep_all = T) %>%
  arrange(GeneID, desc(pct.ratio)) %>% group_by(GeneID) %>% slice_max(n = 1, order_by = pct.ratio)

## write.xlsx(host.markers.int.local.sig.specific, ./Input/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig.xlsx')

## this data will be filtered at the next step to exclude the markers that are specific to host cell 
host.markers.int.local.sig.specific <- read.xlsx('./input/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig.xlsx')
host.markers.int.local.sig.stats <- host.markers.int.local.sig.specific %>% group_by(phase, ref.host) %>% 
  summarise(genes = list(GeneID), num.deg = n())

host.markers.int.local.sig.stats$phase <- factor(host.markers.int.local.sig.stats$phase, levels = c('G', 'SM', 'MC', 'C'))
host.markers.int.local.sig.stats$ref.host <- factor(host.markers.int.local.sig.stats$ref.host, levels = c('cow', 'human'))
colnames(host.markers.int.local.sig.stats) <- gsub("ref.host", "host", colnames(host.markers.int.local.sig.stats))



# Perform phase based DE btw Bdiv_cow vs Bbig , Bdiv_cow vs Bbov upregulated & downregulated

DefaultAssay(S.O.integrated) <- "RNA"
S.O.integrated@meta.data <- S.O.integrated@meta.data %>% mutate(host = ifelse( spp == "Bdiv_human", "human", "cow"))
Idents(S.O.integrated) <- 'spp'

S.O.spp.cow <- subset(S.O.integrated, ident = c('Bbig', 'Bbov', 'Bdiv_cow'))
unique(S.O.spp.cow$phase.spp)
Idents(S.O.spp.cow) <- 'phase.spp'


makeMatchedContrastsCowSpp <- function(S.O.integrated){
  
  objs <- as.character(unique(S.O.integrated@meta.data$phase.spp))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.spp = gsub(":.*", "", ref), ref.phase = gsub(".*:", "", ref), 
           query.spp = gsub(":.*", "", query), query.phase = gsub(".*:", "", query))
  my.contrasts <- contrasts %>% dplyr::filter(ref.phase == query.phase & ref.spp != query.spp)
  my.contrasts <- my.contrasts %>% filter((str_detect( "Bdiv_cow", ref.spp)))
  
  
  return(my.contrasts)
  
}

contrasts <- makeMatchedContrastsCowSpp(S.O.spp.cow)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.host <- gsub(':.*', '', contrasts.groups$ref)



matched.DEGs.Cow.Spp <- mclapply(1:nrow(contrasts.groups), function(i){
  tmp <- FindMarkers(S.O.spp.cow, ident.1 = contrasts.groups$ref[i],
                     ident.2 = c(unlist(contrasts.groups$query[i])), 
                     only.pos = F, verbose = T)
  tmp$ref <- contrasts.groups$ref[i]
  tmp$ref.host <- contrasts.groups$ref.host[i]
  tmp$phase <- contrasts.groups$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})


cow.spp.phase.based.markers <- bind_rows(matched.DEGs.Cow.Spp)
cow.spp.phase.based.markers.sig <- cow.spp.phase.based.markers %>% dplyr::filter(abs(avg_log2FC) > 1 & p_val_adj < 0.01)
cow.spp.phase.based.markers.sig  <- cow.spp.phase.based.markers.sig %>% mutate(pct.ratio = (pct.1 / pct.2))
cow.spp.phase.based.markers.sig  <- cow.spp.phase.based.markers.sig %>% mutate(direction = ifelse(avg_log2FC > 0, "up", "down")) 

## Now consider markers that are phase specific
cow.spp.phase.based.markers.sig.cc <- inner_join(cow.spp.phase.based.markers.sig, cell.cycle.markers.sig, 
                                                 by = c('GeneID', 'ref.host' = 'spp', 'phase'))
cow.spp.phase.based.markers.sig.cc <- cow.spp.phase.based.markers.sig.cc %>% dplyr::select(!contains(".y"))

cow.spp.phase.based.markers.sig.cc.df <- cow.spp.phase.based.markers.sig.cc %>% 
  arrange(GeneID, desc(pct.ratio)) %>% group_by(GeneID, ref.host,) %>% slice_max(n = 1, order_by = pct.ratio)

cow.spp.stat <- cow.spp.phase.based.markers.sig.cc.df %>% 
  group_by(ref.host, phase, direction) %>% summarise(genes= list(GeneID), total = n())



host.markers.int.local.sig.specific <- read.xlsx('./Input/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig.xlsx')
host.stat <- host.markers.int.local.sig.specific %>% group_by(ref.host, phase) %>% summarise(genes= list(GeneID), total = n())
host.stat$direction <- "up"


## removing DE genes in Bdiv_cow vs Bbig , Bdiv_cow vs Bbov 

XX <- full_join(host.stat, cow.spp.stat, by = "phase") %>% rowwise() %>%
  mutate(overlap = length(intersect(unlist(genes.x), unlist(genes.y))), 
         overlap.genes = list(intersect(unlist(genes.x), unlist(genes.y))),
         num.diff.genes = length(setdiff(unlist(genes.x), unlist(genes.y))),
         filt.genes = list(setdiff(unlist(genes.x), unlist(genes.y))))

XX.human.up <- XX %>% filter(ref.host.x == "human" & direction.y == "up")
XX.cow.up <- XX  %>% filter(ref.host.x == "cow", direction.y == "down")


XX.human.up.lng <- XX.human.up  %>% group_by(phase) %>% unnest(filt.genes) 
XX.human.up.lng.filt <- XX.human.up.lng %>% dplyr::select(ref.host.x, phase, filt.genes)

XX.cow.up.lng <- XX.cow.up %>% group_by(phase) %>% unnest(filt.genes)
XX.cow.up.lng.filt <- XX.cow.up.lng %>% dplyr::select(ref.host.x, phase, filt.genes)

XX.human.cow.lng.filt <- rbind(XX.human.up.lng.filt, XX.cow.up.lng.filt)



XX.human.cow.up <-rbind(XX.human.up, XX.cow.up)
XX.human.cow.up$phase <- factor(XX.human.cow.up$phase, levels = c('G', 'SM', 'MC', 'C'))
XX.human.cow.up$ref.host.x <- factor(XX.human.cow.up$ref.host.x, levels = c('cow', 'human'))
colnames(XX.human.cow.up) <- gsub("ref.host", "host", colnames(XX.human.cow.up))

saveRDS(XX.human.cow.up, "./Input/human_vs_cow.rds")

# plot
p <-  ggplot(XX.human.cow.up, aes(x=num.diff.genes, y=host.x, fill = host.x,  color = host.x)) +
  geom_bar(stat="identity", size = 2, width = 0.75, alpha = 0.4, aes(fill=host.x))+
  scale_fill_manual(values = c('cow' = 'darkgoldenrod', 'human' = 'darkolivegreen4'))+
  scale_color_manual(values = c('cow' = 'darkgoldenrod', 'human' = 'darkolivegreen4'))+
  geom_text(aes(label=num.diff.genes), hjust= 1.15, color="black", size=10, fontface = 'bold')+
  theme_bw()+
  facet_grid(phase~., scales = "free", space='free',switch = "y",labeller=label_wrap_gen(multi_line = TRUE))+
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(face="bold", size = 20, color = "black", angle=0, hjust = 1),
    axis.title = element_blank(),
    #axis.title.y = element_text(size=24, face="bold"),
    axis.title.x = element_text(size=24, face="bold"),
    #axis.text.y = element_text(face="bold", size=12, angle=0),
    axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 20),
        legend.title = element_text(face = "bold", size = 22))+
  theme(strip.text.y.left = element_text(angle = 0))

p


XX.human.cow.lng.filt <- rbind(XX.human.up.lng.filt, XX.cow.up.lng.filt)
host.markers.int.local.sig.specific.filt <- 
  host.markers.int.local.sig.specific[host.markers.int.local.sig.specific$GeneID %in% XX.human.cow.lng.filt$filt.genes, ]

# .y correspond to info from cell cycle marker
final.tab <- host.markers.int.local.sig.specific.filt %>% select(GeneID, !contains(c("spp", ".y")))
colnames(final.tab) <- gsub("\\.x", "", colnames(final.tab))
View(final.tab) # for SI tab


#write.xlsx(host.markers.int.local.sig.specific.filt, './rds/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig_filtered_cow_spp_markers_rev.xlsx')

host.markers.int.local.sig.specific.filt <- read.xlsx('./rds/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig_filtered_cow_spp_markers_rev.xlsx')
saveRDS(host.markers.int.local.sig.specific.filt, "./Input/human_vs_cow.rds")

host.markers.int.local.sig.stats.filt <- host.markers.int.local.sig.specific.filt %>% group_by(phase, ref.host) %>% 
  summarise(genes = list(GeneID), num.deg = n())

host.markers.int.local.sig.stats.filt$phase <- factor(host.markers.int.local.sig.stats.filt$phase, levels = c('G', 'SM', 'MC', 'C'))
host.markers.int.local.sig.stats.filt$ref.host <- factor(host.markers.int.local.sig.stats.filt$ref.host, levels = c('cow', 'human'))
colnames(host.markers.int.local.sig.stats.filt) <- gsub("ref.host", "host", colnames(host.markers.int.local.sig.stats.filt))

host.markers.int.local.filt.top <- host.markers.int.local.sig.specific.filt %>% group_by(ref.host) %>% 
  slice_max(n = 1, order_by = avg_log2FC.x)

# feature plot of top markers

S.O.integrated[['pca']]@cell.embeddings[,2] = -S.O.integrated[['pca']]@cell.embeddings[,2]
S.O.integrated[['pca']]@cell.embeddings[,1] = -S.O.integrated[['pca']]@cell.embeddings[,1]
p <- FeaturePlot(object = S.O.integrated, features = rev.default(host.markers.int.local.filt.top$gene.x), label = F,repel = F,
                 #shape.by  = 'spp',
                 split.by = 'host',
                 cols = c("grey", "blue"), reduction = "pca")

plot(p)
num.plt <- 4
p2 <- lapply(1:num.plt, function(i){
  
  plt <- p[[i]] + 
    theme(axis.text.x  = element_text(face = "bold", size = 16, angle = 0, hjust = 0.5),
          axis.text.y  = element_text(face = "bold", size = 18),
          axis.title = element_blank(),
          #axis.title.x = element_text(face = "bold", size = 20), 
          #axis.title.y = element_text(face = "bold", size = 20),
          plot.title = element_text(face = "bold", size = 20)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_rect(color = "black", size = 0.5),
          axis.line.x = element_line(colour = "black", size = 0.5),
          axis.line.y = element_line(colour = "black", size = 0.5)
          
    ) 
  
  
})


DefaultAssay(S.O.integrated) <- "RNA"
S.O.integrated@meta.data <- S.O.integrated@meta.data %>% mutate(host = ifelse( spp == "Bdiv_human", "human", "cow"))
S.O.integrated@meta.data$host.phase <- paste(S.O.integrated@meta.data$host, S.O.integrated@meta.data$cell.cycle.phase, sep = ":")
unique(S.O.integrated@meta.data$host.phase)
Idents(S.O.integrated) <- 'host.phase'

top.markers <- host.markers.int.local.filt.top$GeneID

# violin plot of top markers
p <- VlnPlot(subset(S.O.integrated, idents = c('cow:C', 'human:C')), features = gsub('_', '-', top.markers))

pp <- lapply(1:length(unique(top.markers)), function(i){
  
  plt <- p[[i]] + 
    scale_fill_manual(values = c('cow:C' = 'darkgoldenrod', 'human:C' = 'darkolivegreen4'))+
    theme(axis.text.x  = element_text(face = "bold", size = 16, angle = 0, hjust = 0.5),
          axis.text.y  = element_text(face = "bold", size = 18),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(face = "bold", size = 20),
          plot.title = element_text(face = "bold", size = 20)) + 
    ylab("expr")
  
  
})

p2 <- grid.arrange(pp[[2]], pp[[1]], nrow = 2)

