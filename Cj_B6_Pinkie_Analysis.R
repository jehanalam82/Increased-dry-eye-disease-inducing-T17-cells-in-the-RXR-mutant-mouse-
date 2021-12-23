install.packages('Seurat', version = "3.X.X")
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(cowplot)
B.data <- Read10X("H:/Documents/Analysis/B6")
P.data <- Read10X("H:/Documents/Analysis/Pinkie")
B <- CreateSeuratObject(counts = B.data, project = "C57BL_6", min.cells = 3, min.features = 200)
P <- CreateSeuratObject(counts = P.data, project = "Pinkie", min.cells = 3, min.features = 200)
B<- NormalizeData(B, normalization.method = "LogNormalize", scale.factor = 10000)
P<- NormalizeData(P, normalization.method = "LogNormalize", scale.factor = 10000)
B<- FindVariableFeatures(B, selection.method = "vst", nfeatures = 2000)
P <- FindVariableFeatures(P, selection.method = "vst", nfeatures = 2000)
BP.list <- list(B, P)
## PERFORM INTEGRATION
Common.anchors <- FindIntegrationAnchors(object.list = BP.list, dims = 1:30)
Common.combined <- IntegrateData(anchorset = Common.anchors, dims = 1:30)
DefaultAssay(Common.combined) <- "integrated"
Common.combined <- ScaleData(Common.combined, verbose = FALSE)
Common.combined <- RunPCA(Common.combined, npcs = 30, verbose = FALSE)
Common.combined <- RunUMAP(Common.combined, reduction = "pca", dims = 1:30)
Common.combined <- FindNeighbors(Common.combined, reduction = "pca", dims = 1:30)
Common.combined <- FindClusters(Common.combined, resolution = 0.5)
p1 <- DimPlot(Common.combined, reduction = "umap" , group.by = "orig.ident")
p2 <- DimPlot(Common.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(Common.combined, reduction = "umap", split.by = "orig.ident")
new.cluster.ids <- c("M??MHCIIlow", "Neutrophils-LCN2low", "????-T", "Monocytes", "M??", "ILC2","NK cells", "CD4", "Neutrophils-LCN2high", "cDC2-Retnlahigh", "CD8", "B ", "cDC1", "Proliferating cell", "Mast", "Treg", "Naive CD4+", "mDC", "pDC")
names(new.cluster.ids) <- levels(Common.combined)
Common.combined <- RenameIdents(Common.combined, new.cluster.ids)
DimPlot(Common.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() 
DefaultAssay(Common.combined) <- "RNA"

########cells cluster count#########
Common.combined[["my.clusters"]] <- Idents(object = Common.combined)
table(object@meta.data$my.clusters, object@meta.data$orig.ident)
#########DEG-across the clusters###########
BP.list <-SplitObject(Common.combined, split.by = "orig.ident")
samples <- names(BP.list )
master.markers <- NULL
for (i in seq_along(BP.list)) {
  tmp.markers <- FindAllMarkers(BP.list[[i]], only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  tmp.markers$orig.ident <- samples[i]
  master.markers <- rbind(master.markers, tmp.markers)
}

master.markers %>%
  group_by(cluster, orig.ident) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(tmp.markers, file = "C:/Users/jalam/Documents/BP_DEG_IPA_tmp" )
write.csv(master.markers, file = "C:/Users/jalam/Documents/BP_IPA_master.markers" )

###########feature-plot#######

FeaturePlot(Common.combined, features = c("Serpinb2", "Lcn2", "Trdc", "Apoe","Cd209a","Gata3","Gzma","Cd28", "S100a9", "Retnla","Xcl1","Igkc","Naaa","Stmn1", "Mcpt4", "Foxp3", "Slpr1", "Fscn1","Siglec"))

VlnPlot(Common.combined, features = c("Il17a", "Il17f","Ltb","Cxcr6", "Rorc","Il1r1" ), group.by = "orig.ident",
        
        pt.size = 1.5, combine = TRUE)
########use "split.by" option to compare between differnt custers######
###########Violen-plot with significance #######

vp_case1 <- function(gene_signature, file_name, test_sign){
plot_case1 <- function(signature, y_max = NULL){
VlnPlot(Common.combined, features = signature,
pt.size = 0.1,
group.by = "orig.ident",
y.max = y_max 
) + stat_compare_means(comparisons = test_sign, label = "p.signif")
}
plot_list <- list()
y_max_list <- list()
for (gene in gene_signature) {
plot_list[[gene]] <- plot_case1(gene)
y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) 
plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
}
cowplot::plot_grid(plotlist = plot_list)
file_name <- paste0(file_name, ".pdf")
ggsave(file_name, width = 4, height = 4) 
}
gene_sig <- c("Il17a", "Il17f","Ltb","Cxcr6", "Rorc","Il1r1")
comparisons <- list(c("C57BL_6", "Pinkie"))
vp_case1(gene_signature = gene_sig, file_name = "signatureGenes2", test_sign = comparisons)

saveRDS(Common.combined, file='C:/Users/jalam/Documents/Jahan Paper Writting/Single Cell RNA Seq Paper/B6_Pinkie_Frontiers.rds')


################################
#Top 100 markers between strains
################################

Common.combined <- ScaleData(Common.combined)

#Get the Differential Genes
Idents(Common.combined) <- "orig.ident"
strain.markers <- FindMarkers(Common.combined, ident.1 = "C57BL_6", 
                              only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

#Selecting Top 40 Genes 
library(dplyr)
top40 <- strain.markers %>%
  filter(p_val_adj < 0.05) %>% #significant filter
  slice_max(abs(avg_log2FC), n = 40) %>% #Using absolute for high or low genes
  arrange(avg_log2FC) #Order by value
top40genes <- rownames(top40)


#Heatmap Across All Cells
DoHeatmap(Common.combined, features = top100genes, draw.lines = FALSE)

#Average-based Heatmap
Avg <- AverageExpression(Common.combined, return.seurat = TRUE)
DoHeatmap(Avg, features = top100genes, draw.lines = FALSE, slot = "data")

VlnPlot(Common.combined, features = c("Tnf"), group.by = "orig.ident",
        
        pt.size = 1.5, combine = FALSE)
