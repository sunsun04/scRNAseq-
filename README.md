< Data preprocessing > 
#SVF(sample)
SVF.matrix <- Read10X("/Users/sunyoung/Desktop/SVF sc data/") # matrix file은 barcode, matrix,features file 세 개를 모두 포함하는 파일
SVF.data <- Read10X("/Users/sunyoung/Desktop/SVF sc data/")
SVF <- CreateSeuratObject(counts = SVF.data, project = "SVFsc",min.cells = 3,min.features = 200)
SVF[["percent.mt"]] <- PercentageFeatureSet(SVF,pattern = "^mt-")
# mitocontrial gene pattern of mouse sample : "^mt-"(mouse sample은 HUMNA과 다르게 반드시 소문자 mt로 입력해야함!!!)
VlnPlot(SVF,features = c("nFeature_RNA","nCount_RNA", "percent.mt"),ncol = 3)
plot1 <- FeatureScatter(SVF,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(SVF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
SVF <- subset(SVF, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
save(SVF,file = "/Users/sunyoung/Desktop/SVF sc data/SVF.rda")

# LNSC(control)
LNSC.data <- Read10X("/Users/sunyoung/Desktop/single cell transcriptomics/LNSC sc data/GSM6093723/")
LNSC <- CreateSeuratObject(counts = LNSC.data,project = "LNSCsc", min.cells = 3,min.features = 200)
LNSC[["percent.mt"]] <- PercentageFeatureSet(LNSC,pattern = "^mt-")
VlnPlot(LNSC,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
plot1 <- FeatureScatter(LNSC,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(LNSC,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1 + plot2

write.csv(SVF@meta.data,"/Users/sunyoung/Documents/transcriptome /SVF transcriptome/SVF.csv")
write.csv(LNSC@meta.data, "/Users/sunyoung/Documents/transcriptome /SVF transcriptome/LNSC.csv")

# excel에 sample lane 추가 
SVF$sample <- "SVF"
LNSC$sample <- "LNSC"

write.csv(SVF@meta.data,"/Users/sunyoung/Documents/transcriptome /SVF transcriptome/SVF_QC.csv")
write.csv(LNSC@meta.data,"/Users/sunyoung/Documents/transcriptome /SVF transcriptome/LNSC_QC.csv")

#merge : SVF + LNSC 
publishedSVF <- merge( x= SVF ,y = LNSC)
write.csv(publishedSVF@meta.data,"/Users/sunyoung/Documents/transcriptome /SVF transcriptome/publishedSVF_total.csv")

# Normalization
publishedSVF.list <- SplitObject(publishedSVF,split.by = "sample")
publishedSVF.list <- lapply(X = new_SVFLNSC.list,FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,selection.method = "vst",nfeatures = 2000) 
  })
Features <- SelectIntegrationFeatures(object.list = publishedSVF.list)

#integration 
#integration 한 후에는 QC필요 없음 

immune.anchors <- FindIntegrationAnchors(object.list = publishedSVF.list,anchor.features = Features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(immune.combined, reduction = "umap",split.by = "sample")
saveRDS(immune.combined, file = "/Users/sunyoung/Documents/transcriptome /SVF transcriptome/intergration_SVF.rds")


#for finding DEG
DefaultAssay(immune.combined) <- "RNA"
s.markers <- FindAllMarkers(immune.combined,only.pos = TRUE ,min.pct =0.25, logfc.threshold = 0.25)
head(s.markers)
write.csv(s.markers,"/Users/sunyoung/Documents/transcriptome /SVF transcriptome/s_Markers.csv")

immune.combined <- RenameIdents(immune.combined,`0` = "MSC1", `1` = "FDC", `2` = "FRC",`3` = "Mono", `4` = "DC", `5` = "TRC", `6` = "MSC2", `7` = "BEC", `8` = "pDC", `9` = "CD8T",`10` = "DNC", `11` = "Mac1&2", `12` = "CD4T", `13` = "B cell", `14` = "LEC",'15'= "PvC" , '16' = "neuronal")
DimPlot(immune.combined, reduction = "umap",split
        
# remove clusters that have CD45 gene
immune.combined_rm_ptprc <- subset(immune.combined2, idents = c("Mono","DC","pDC","CD8T","Mac1&2","CD4T","B cell"), invert = TRUE) #invert =TURE : 언급된 cluster 제거, invert = TRUE 안하면 언급된 cluster가 남고 나ㅏ머지는 표시 안됨. 

# re-clustering code by 경석쌤 
immune.combined_rm_ptprc_final2 <- FindVariableFeatures(immune.combined_rm_ptprc_final2,selection.method = "vst", nfeatures = 2000)
DefaultAssay(immune.combined_rm_ptprc) <- "integrated"
all.genes <-row.names(immune.combined_rm_ptprc_final2)
immune.combined_rm_ptprc_final2 <- ScaleData(immune.combined_rm_ptprc_final2,features = all.genes)
immune.combined_rm_ptprc_final2 <- RunPCA(immune.combined_rm_ptprc_final2,features = VariableFeatures(object = immune.combined_rm_ptprc_final2))
ElbowPlot(immune.combined_rm_ptprc_final2)
immune.combined_rm_ptprc_final2 <-FindNeighbors(immune.combined_rm_ptprc_final2,dims = 1:30)
immune.combined_rm_ptprc_final2 <- FindClusters(immune.combined_rm_ptprc_final2,resolution = 0.5)
immune.combined_rm_ptprc_final2 <-RunUMAP(immune.combined_rm_ptprc_final2,dims = 1:30)
# dims ratio와 resolution value에 따라 cluster가 달라진다.왜지..? 적절한 값을 찾아야한다..


DimPlot(immune.combined_rm_ptprc_final2, reduction = "umap",split.by = "sample")
SVFLNSCmarkers <- FindAllMarkers(immune.combined_rm_ptprc_final2,only.pos = TRUE ,min.pct =0.25, logfc.threshold = 0.25)
write.csv(SVFLNSCmarkers,"C:/Users/HI/Documents/R data/SVFLNSCmarkers.csv")

# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# singleR labeling 
BiocManager::install("SingleCellExperiment")
ref <- celldex::MouseRNAseqData()
results3 <- SingleR(test = as.SingleCellExperiment(immune.combined_rm_ptprc_final2),ref = ref, labels = ref$label.main)
results3 # for results datasheet confirmation
immune.combined_rm_ptprc_final2$singlr_labels <- results3$labels
immune.combined_rm_ptprc_final2[[]] # for  labelging confirmation
DimPlot(immune.combined_rm_ptprc_final2, reduction = 'umap', group.by = 'singlr_labels',label = TRUE)

