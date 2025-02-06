#### Data processing
## qc
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(cowplot)
library(viridis)
library(future)
library(qs)

seurat.data <- qread(file = "./data/seurat_decontX.qs")
seurat.data[["percent.mt"]] <- PercentageFeatureSet(seurat.data,pattern = "^mt-")
head(seurat.data@meta.data, 5)

options(repr.plot.width=10, repr.plot.height=5)
VlnPlot(seurat.data,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        group.by = "group")

seurat.data <- subset(seurat.data,subset =  nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                          nCount_RNA > 500 & nCount_RNA < 50000  & percent.mt < 15)


## Standardization, Dimensionality Reduction, Clustering
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F)

seurat.data <- seurat.data %>% 
  RunPCA(verbose = F) %>% 
  RunUMAP(reduction = "pca", dims = 1:30, n.neighbors = 10, verbose = F)
ElbowPlot(seurat.data,ndims = 50)

options(repr.plot.width = 10, repr.plot.height = 4.5)
p1.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "pca", group.by = "sampleID")+NoAxes()+ggtitle("Before_PCA"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("Before_UMAP"),
                      guides = "collect"
)
p1.compare

seurat.data <- FindNeighbors(seurat.data, k.param = 20,reduction = "pca", dims = 1:30)
for (res in c(0.05,0.1,0.2,0.3,0.5)){
  print(res)
  seurat.data <- FindClusters(seurat.data,resolution = res, algorithm = 1)
}

cluster_umap <- wrap_plots(ncol = 2,
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.05", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T) & NoAxes()
                           )
cluster_umap

## annotation
seurat.data@meta.data$seurat_clusters <- seurat.data@meta.data$RNA_snn_res.0.1
Idents(seurat.data) <- "seurat_clusters"

check_genes <- c("Hbb-bt","Hba-a2","Hba-a1", ## erythrocyte
                 "Ppbp","Gp1bb", ## platelet
                  "Pdgfra", ## Fibroblasts
                 "Myh11","Acta2", ## SMC
                 "Ager","Rtkn2", ## AT1
                 "Lamp3","Sftpc", ## AT2
                 "Tppp3","Foxj1", ## Ciliated
                 "Scgb3a2","Muc5b","Scgb1a1", ## Club
                 "Chil3","Krt79", ## AM AM_proliferating
                 "C3ar1", ## IM
                 "Pecam1", ## EC
                 "Cd68","Plac8", ## Monocytes
                 "Cd79b","Ms4a1", ## B_cells
                 "Clec10a","Ccl17", ## DC DC_proliferating
                 "Cd3g","Cd3e","Trbc2", ## T_cells T_cells_proliferating
                 "S100a8","S100a9", ## Neutrophils
                 "Nkg7","Gzma", ## NK_cells
                 "Upk1b","Msln", ## Mesotheliocytes
                 "Mcpt8","Cpa3","Fcer1a", ## Basophils
                 "Mki67","Top2a") ## _proliferating

options(repr.plot.width = 7.5, repr.plot.height = 10)
DotPlot(object = seurat.data, features = check_genes, 
        assay = "RNA",scale = T) + coord_flip()

head(Idents(seurat.data))
Idents(seurat.data) <- "seurat_clusters"
DimPlot(seurat.data, reduction = "umap", label = T)
seurat.data <- RenameIdents(seurat.data,
                      "0"="T_cells",
                      "1"="B_cells", 
                      "2"="AM", 
                      "3"= "IM", 
                      "4"= "AT2", 
                      "5"= "EC",
                      "6"= "Neutrophils", 
                      "7"= "Fibroblasts", 
                      "8"= "DC",
                      "9"= "Monocytes", 
                      "10"= "NK_cells", 
                      "11"= "T_cells",
                      "12"= "Ciliated", 
                      "13"= "SMC",
                      "14"= "T_cells_proliferating",
                      "15"= "Club",
                      "16"= "Basophils",
                      "17"= "AT1",
                      "18"= "Unknown",
                      "19"= "Unknown",
                      "20"= "IM",
                      "21"= "Mesotheliocytes"
                
)

#### Visualization
## cell annotaion
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(qs)
library(cowplot)
source("./code/custom_seurat_functions.R")

options(repr.plot.width = 5, repr.plot.height = 6)
p1 <- TSNE.UMAP.Plot(seurat.data,
                     "RNA_snn_res.0.1",
                     Style=1,
                     plot.title = NA,
                     legend.position = "none",
                     #legend.point.size  = 3,
                     reduction = "umap",
                     label = T,
                     label.size = 5,
                     point.size = 1e-1)+labs(title = "Seurat_clusters")
p1


options(repr.plot.width = 5, repr.plot.height = 6)
p2 <- TSNE.UMAP.Plot(seurat.data,
                     "celltype_sub",
                     Style=1,
                     plot.title = NA,
                     legend.position = "none",
                     #legend.point.size = 3,
                     reduction = "umap",
                     label = T,
                     label.size = 5,
                     point.size = 1e-1)+labs(title = "Celltype")
p2

p1 | p2

## cell markers
seurat.data$celltype_sub <- factor(seurat.data$celltype_sub,levels = c("Fibroblasts","SMC","AT1","AT2","Ciliated","Club","AM","AM_proliferating",
                                                                       "IM","EC","Monocytes","B_cells","DC","T_cells","T_cells_proliferating",
                                                                       "Neutrophils","NK_cells","Mesotheliocytes","Basophils"))
Idents(seurat.data) <- seurat.data$celltype_sub
levels(seurat.data)

options(repr.plot.width = 6, repr.plot.height = 2.5)
check_genes <- c("Pdgfra", ## Fibroblasts
                 "Myh11","Acta2", ## SMC
                 "Ager","Rtkn2", ## AT1
                 "Lamp3","Sftpc", ## AT2
                 "Tppp3","Foxj1", ## Ciliated
                 "Scgb3a2","Muc5b","Scgb1a1", ## Club
                 "Chil3","Krt79", ## AM AM_proliferating
                 "C3ar1", ## IM
                 "Pecam1", ## EC
                 "Cd68","Plac8", ## Monocytes
                 "Cd79b","Ms4a1", ## B_cells
                 "Clec10a","Ccl17", ## DC DC_proliferating
                 "Cd3g","Cd3e","Trbc2", ## T_cells T_cells_proliferating
                 "S100a8","S100a9", ## Neutrophils
                 "Nkg7","Gzma", ## NK_cells
                 #"Upk1b","Msln", ## Mesotheliocytes
                 "Mcpt8","Cpa3","Fcer1a", ## Basophils
                 "Mki67","Top2a") ## _proliferating
p.dot = DotPlot_2(object = seurat.data,
                  Combine = F, dot.range.max = 3,
                  dot.range.min = 0, label.size = 4,
                  features = check_genes, legend.key.size = 0.4)&
  scale_color_distiller(palette = 'RdYlBu')
p.dot

## Percentage of cell number
if(T){
  text.size = 12
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), 
                   panel.grid=element_blank(), 
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
}

table(seurat.data$celltype_sub2)
color.use=c('Epithelial'='#1f77b4', 'T_cells'='#ff7f0e','AM'= '#279e68', 'DC'='#d62728', 
            'IM'='#aa40fc', 'B_cells'='#8c564b', 'EC'='#e377c2', 
            'Monocytes'= '#b5bd61', 'Fibroblasts'='#17becf','Neutrophils'='#aec7e8',
            'NK_cells'= '#98df8a','SMC' = '#FF1493','Basophils' = '#ffa500')

options(repr.plot.width = 5, repr.plot.height = 5)
p1 = plot.clusters.group(data = seurat.data,
                         Group = "group",
                         legend.position = "none",
                         celltype.id = "celltype_sub2",
                         widths = c(2,1),
                         log = F,
                         text.size = 12,
                         order = F,
                         legend.title = "Cell type",
                         color = 1,
                         xlab = "",
                         cell.counts = F)&mytheme&
  scale_fill_manual(values = color.use);p1

meta.data  <- seurat.data@meta.data
head(meta.data)
table(seurat.data$group)
table(seurat.data$celltype_sub2)

## OR
meta.data  <- seurat.data@meta.data
head(meta.data)
table(seurat.data$group)
table(seurat.data$celltype_sub2)
group.color = c("IR_17Gy" = "#f39b7f","Control" = "#00a087")
dat.plot = meta.data
dat.plot$group = factor(dat.plot$group,levels = c("Control","IR_17Gy"))
options(repr.plot.width = 11, repr.plot.height = 3.8)
p2 = gg.tissueDist(cellInfo.tb = dat.plot,
                   meta.cluster = dat.plot$celltype_sub2,
                   colname.patient = "sampleID",
                   loc = dat.plot$group,
                   cuts =c(0.0000000, 0.0000001, 0.8, 1.5, 2, Inf),
                   bin.label = c("-", "+/-", "+", "++", "+++"),
                   verbose = 1,
                   z.hi.OR = 4,
                   z.hi.ROE = 2,
                   text.size = 10,
                   text_color = paletteer::paletteer_d('ggsci::category20c_d3')
)

## NETosis_Score
netosis_features <- list(c("Csf3r","Cxcr2","Il1b","Selplg","Mmp9","Fpr2","Itgam","Padi4",
                      "Fpr1","Itgb2","Fcgr4","Cybb","Bst1","Mgam","Mpo"))

seurat.data <- AddModuleScore(seurat.data,
                          features = netosis_features,
                          ctrl = 100,
                          name = "NETosis_Score")
VlnPlot(seurat.data,features = 'NETosis_Score', 
        pt.size = 0, adjust = 2,group.by = "celltype_sub2")


#### Macrophages subclustering
Mac <- subset(seurat.data,celltype_sub %in% c("AM","IM"))
Mac <- Mac %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F) %>% 
  RunPCA(verbose = F) %>% 
  RunUMAP(reduction = "pca", dims = 1:30,verbose = F) %>% 
  FindNeighbors(k.param = 20, dims = 1:30) %>% 
  FindClusters(resolution = c(0.01,0.03,0.05,0.1,0.2,0.3,0.5), algorithm = 1)

Mac@meta.data$seurat_clusters <- Mac@meta.data$RNA_snn_res.0.05
Idents(Mac) <- "seurat_clusters"

Mac <- RenameIdents(Mac,
                    "0"="IM1",
                    "1"="AM1", 
                    "2"="AM2", 
                    "3"= "IM2", 
                    "4"= "AM3", 
                    "5"= "IM3"
)

DimPlot(Mac, reduction = "umap", label = T)
Mac@meta.data$celltype_sub_Mac = as.character(Idents(Mac))

## M1/M2 scoring
M1_Polarization <- unique(c("IL12", "IL23", "IL12", "TNF", "IL6", "CD86", "IL1B", "MARCO", "NOS2", "IL12", "CD64", "CD80", "CXCR10", "IL23", "CXCL9", "CXCL10", "CXCL11", "CD86", "IL1A", "IL1B", "IL6", "TNF", "CCL5", "IRF5", "IRF1", "CD40", "IDO1", "KYNU", "CCR7"))
M2_Polarization <- unique(c("ARG1", "ARG2", "IL10", "CD32", "CD163", "CD23", "CD200R1", "PDCD1LG2", "CD274", "MARCO", "CSF1R", "CD206", "IL1RN", "IL1R2", "IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", "CCL24", "LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", "CTSA", "CTSB", "CSTC", "CTSD", "TGFB1", "TGFB2", "TGFB3", "MMP14", "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASL", "TNFSF12", "TNFSF8", "CD276", "VTCN1", "MSR1", "FN1", "IRF4"))

library(homologene)
M1_Polarization_mouse <- homologene(M1_Polarization, inTax = 9606, outTax = 10090)
colnames(M1_Polarization_mouse) <- c("gene_human","gene_mouse","ID_human","ID_mouse")

M2_Polarization_mouse <- homologene(M2_Polarization, inTax = 9606, outTax = 10090)
colnames(M2_Polarization_mouse) <- c("gene_human","gene_mouse","ID_human","ID_mouse")

sce <- Mac
sce <- AddModuleScore(sce,
                      features = list(M1_Polarization_mouse$gene_mouse),
                      ctrl = 100,
                      name = "M1_Polarization")
colnames(sce@meta.data)
colnames(sce@meta.data)[19] <- "M1_Polarization"

library(Nebulosa)
plot_density(sce,"M1_Polarization")


sce <- AddModuleScore(sce,
                      features = list(M2_Polarization_mouse$gene_mouse),
                      ctrl = 100,
                      name = "M2_Polarization")
colnames(sce@meta.data)
colnames(sce@meta.data)[20] <- "M2_Polarization"

library(Nebulosa)
plot_density(sce,"M2_Polarization")


#### cellchat
## circle plot
gg11 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2),sources.use = 10,targets.use = c(1:9,11:13),arrow.size = 0.1,title.name = "")
gg22 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2),sources.use = c(1:9,11:13),targets.use = 10,arrow.size = 0.1,title.name = "")

##  L-R Communication probability
levels(cellchat@meta$celltype_sub2)
options(repr.plot.width = 7, repr.plot.height = 4)
netVisual_bubble(cellchat, sources.use = 10, targets.use = c(3,5,6),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(2,3,5,6), targets.use = 10,  comparison = c(1, 2), angle.x = 45)

i=1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

gg1 <- netVisual_bubble(cellchat, sources.use = 10, targets.use = c(3,5,6),  comparison = c(1, 2), max.dataset = 2, signaling = pathway.union,title.name = "Increased signaling in IR",angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 10, targets.use = c(3,5,6),  comparison = c(1, 2), max.dataset = 1, signaling = pathway.union,title.name = "Decreased signaling in IR",angle.x = 45, remove.isolate = T)
options(repr.plot.width = 10, repr.plot.height = 4)
gg1 + gg2

gg3 <- netVisual_bubble(cellchat, sources.use = c(2,3,5,6), targets.use = 10,  comparison = c(1, 2), signaling = pathway.union,max.dataset = 2, title.name = "Increased signaling in IR",angle.x = 45, remove.isolate = T,font.size = 8)
gg4 <- netVisual_bubble(cellchat, sources.use = c(2,3,5,6), targets.use = 10,  comparison = c(1, 2), signaling = pathway.union,max.dataset = 1, title.name = "Decreased signaling in IR",angle.x = 45, remove.isolate = T,font.size = 8)
options(repr.plot.width = 10, repr.plot.height = 6)
gg3 + gg4


## Transcriptional regulation
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(qs)
source("R/compute_module_score.R")
source("source_R/custom_seurat_functions.R")

seu <- qread("Mac_sub.qs")
regulons <- clusterProfiler::read.gmt("./data/02-pyscenic.regulons.gmt")
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
print(regulon.list[1])

## RAS matrix
seu <- ComputeModuleScore(seu, gene.sets = regulon.list, min.size = 10, cores = 10)

## regulon module
rasMat <- seu[["AUCell"]]@data
rasMat <- t(rasMat)
pccMat <- cor(rasMat)

# Connection Specificity Index (CSI)
CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}

csiMat <- pbapply::pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)

## h for cut tree
library(dendextend)
library(ggsci)
h = 7
row_dend = as.dendrogram(hclust(dist(pccMat), method = "complete"))
clusters <- dendextend::cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))
plot(row_dend)

library(ComplexHeatmap)
library(circlize)

col_range = c(0, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))

for(i in nrow(pccMat)) {
  pccMat[i, i] <- 0
}

ht <- Heatmap(
  matrix = pccMat,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)

lgd <- Legend(
  col_fun = col_fun,
  title = "",
  at = col_range,
  labels = c("low", "high"),
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
)

{
  draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
  decorate_heatmap_body("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    x1 = x1/length(ind)
    x2 = x2/length(ind)
    grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2),
              hjust = 0, vjust = 0, default.units = "npc",
              gp = gpar(fill=NA, col="#FCB800", lwd=3))
    grid.text(label = paste0("M",clusters),
              x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
              default.units = "npc",
              hjust = 1, vjust = 0.5,
              gp = gpar(fontsize=12, fontface="bold"))
  })
  decorate_column_dend("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
              default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
  })
}

tree <- column_dend(ht)
ind <- cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon=names(ind), cluster=paste0("M",ind))
table(regulon.clusters$cluster)

k = length(clusters)
cell.info <- seu@meta.data
moduleRasMat <- lapply(paste0("M",1:k), function(x){
  regulon.use <- subset(regulon.clusters, cluster == x)$regulon
  rowMeans(rasMat[, regulon.use, drop=FALSE])
})
names(moduleRasMat) <- paste0("M",1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
cell.info <- cbind(cell.info, moduleRasMat[rownames(cell.info), ])

cell.info <- cbind(cell.info, FetchData(seu, vars = paste0("umap_", 1:2)))

p.list <- lapply(paste0("M",1:k), function(module){
  data.use <- cell.info
  expression.color <- c("darkblue", "lightblue", "green", "yellow", "red")
  max.val <- quantile(data.use[, module], 0.99)
  low.val <- quantile(data.use[, module], 0.1)
  data.use[, module] <- ifelse(data.use[, module] > max.val, max.val, data.use[, module])
  ggplot(data.use, aes(umap_1, umap_2, color=get(module))) +
    geom_point(size=0.05) +
    theme_bw(base_size = 11,base_rect_size = 1) +
    ggtitle(module) +
    facet_wrap(~group) +
    scale_color_gradientn(name = NULL, colors = expression.color) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.title.x=element_text(size=15,face="bold"),
          axis.title.y=element_text(size=15,face="bold"),
          axis.text.y=element_text(size=12,face="bold"),
          axis.text.x=element_text(size=12,face="bold")
          
    )
})

cowplot::plot_grid(plotlist = p.list, ncol = 3)

## M8 module
cp.list[[8]]








































