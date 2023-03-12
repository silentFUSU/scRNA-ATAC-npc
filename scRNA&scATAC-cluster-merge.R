set.seed(1) 
myPaths <- .libPaths()
new <- c('/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library')
myPaths <- c(myPaths, new) 
.libPaths(myPaths)
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/")
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# Load the PBMC dataset
filenames=c("24h","44h","46h","48h","50h","52h","54h","72h","D12","D18","D25","D35","D6","WT")

#############scRNA_quality_control###################
for(i in 1:14){
  filename <- filenames[i]
  pbmc.data <- Read10X(data.dir = paste0("data/merge_data/cellranger_arc/NPC_",filename,"/outs/filtered_feature_bc_matrix/"))
  # Initialize the Seurat object with the raw (non-normalized data).
  pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = filename, min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  # Visualize QC metrics as a violin plot
  VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(filename = paste0("result/merge_data/quality_control/",filename,"/raw_feature_count_mt.png"),width = 15,height = 10)
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  # 
  # plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  # plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  # plot1 + plot2
  # ggsave(filename = paste0("result/merge_data/quality_control/",filename,"/feature_count_and_mt_count.png"),width = 15,height = 10)
  pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 20)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- ScaleData(pbmc)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- RunUMAP(pbmc, dims = 1:30)
  pbmc <- FindNeighbors(pbmc, dims = 1:30)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  DimPlot(pbmc,label = T,repel = T,pt.size = 1.5)
  ggsave(paste0("result/merge_data/cluster/",filename,"_0.5.png"),width = 15,height = 10)
  VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
  ggsave(filename = paste0("result/merge_data/quality_control/",filename,"/feature_count_and_mt_count_after_filter.png"),width = 15,height = 10)
  saveRDS(pbmc,paste0("data/merge_data/scRNA/cluster_0.5/",filename,".rds"))
}

#############scATAC quality control & jiont data quality control ###################
filenames=c("24h","44h","46h","48h","50h","52h","54h","72h","D12","D18","D25","D35","D6","WT")
for(i in 1:14){
  filename <- filenames[i]
  counts <- Read10X_h5(filename = paste0("data/merge_data/cellranger_arc/NPC_",filename,"/outs/filtered_feature_bc_matrix.h5"))
  metadata <- read.csv(
    file = paste0("data/merge_data/cellranger_arc/NPC_",filename,"/outs/per_barcode_metrics.csv"),
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0("data/merge_data/cellranger_arc/NPC_",filename,"/outs/atac_fragments.tsv.gz"),
    min.cells = 10,
    min.features = 200
  )
  
  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  Annotation(pbmc) <- annotations
  # compute nucleosome signal score per cell
  pbmc <- NucleosomeSignal(object = pbmc)
  # compute TSS enrichment score per cell
  pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
  pbmc$pct_reads_in_peaks <- pbmc$atac_peak_region_fragments / pbmc$atac_fragments * 100
  print(quantile(pbmc$pct_reads_in_peaks,c(0.01,0.99)))
  print(quantile(pbmc$atac_peak_region_fragments,c(0.01,0.95)))
  print(quantile(pbmc$TSS.enrichment,c(0.05,0.99)))
  print(quantile(pbmc$nucleosome_signal,c(0.01,0.99)))
  VlnPlot(
     object = pbmc,
     features = c('pct_reads_in_peaks', 'atac_peak_region_fragments',
                  'TSS.enrichment','nucleosome_signal'),
     pt.size = 0.1,
     ncol = 4
   )
  ggsave(paste0("result/merge_data/quality_control/",filename,"/raw_atac_reads_in_peaks.png"),width = 20,height = 10)
  pbmc <- subset(
    x = pbmc,
    subset = atac_peak_region_fragments > 300 &
      atac_peak_region_fragments < 10000 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 1.35 &
      nucleosome_signal > 0.55 &
      TSS.enrichment > 3
  )
  pbmc.rna <- readRDS(paste0("data/merge_data/scRNA/cluster_0.5/",filename,".rds"))
  rna_barcode <- rownames(pbmc.rna@meta.data)
  pbmc <- subset(x = pbmc,cells = rownames(pbmc@meta.data)[which(rownames(pbmc@meta.data)%in%rna_barcode)] )
  pbmc <- RunTFIDF(pbmc)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  pbmc <- RunSVD(pbmc)
  pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
  pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
  pbmc <- FindClusters(object = pbmc, verbose = TRUE, algorithm = 3,resolution = 0.5)
  DimPlot(pbmc,label = T,repel = T,pt.size = 1.5)
  ggsave(paste0("result/merge_data/cluster/",filename,"_atac_0.5.png"),width = 15,height = 10)
  saveRDS(pbmc,paste0("data/merge_data/scATAC/cluster_0.5/",filename,"_atac.rds"))
  atac_barcode <- rownames(pbmc@meta.data)
  pbmc.rna <- subset(pbmc.rna, cells = rownames(pbmc.rna@meta.data)[which(rownames(pbmc.rna@meta.data)%in%atac_barcode)])
  pbmc.rna <- NormalizeData(pbmc.rna, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc.rna <- FindVariableFeatures(pbmc.rna, selection.method = "vst", nfeatures = 2000)
  pbmc.rna <- ScaleData(pbmc.rna)
  pbmc.rna <- RunPCA(pbmc.rna, features = VariableFeatures(object = pbmc.rna))
  pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)
  pbmc.rna <- FindNeighbors(pbmc.rna, dims = 1:30)
  pbmc.rna <- FindClusters(pbmc.rna, resolution = 0.5)
  DimPlot(pbmc.rna,label = T,repel = T,pt.size = 1.5)
  ggsave(paste0("result/merge_data/cluster/",filename,"_0.5_final.png"),width = 15,height = 10)
  saveRDS(pbmc.rna,paste0("data/merge_data/scRNA/cluster_0.5/",filename,".rds"))
}
#######################scRNA-seq merge####################################
npc44h <- readRDS("data/merge_data/scRNA/cluster_0.5/44h.rds")
npc46h <- readRDS("data/merge_data/scRNA/cluster_0.5/46h.rds")
npc48h <- readRDS("data/merge_data/scRNA/cluster_0.5/48h.rds")
npc50h <- readRDS("data/merge_data/scRNA/cluster_0.5/50h.rds")
npc52h <- readRDS("data/merge_data/scRNA/cluster_0.5/52h.rds")
npc54h <- readRDS("data/merge_data/scRNA/cluster_0.5/54h.rds")
npc24h <- readRDS("data/merge_data/scRNA/cluster_0.5/24h.rds")
npc72h <- readRDS("data/merge_data/scRNA/cluster_0.5/72h.rds")
WT <- readRDS("data/merge_data/scRNA/cluster_0.5/WT.rds")
D6 <- readRDS("data/merge_data/scRNA/cluster_0.5/D6.rds")
D12 <- readRDS("data/merge_data/scRNA/cluster_0.5/D12.rds")
D18 <- readRDS("data/merge_data/scRNA/cluster_0.5/D18.rds")
D25 <- readRDS("data/merge_data/scRNA/cluster_0.5/D25.rds")
D35 <- readRDS("data/merge_data/scRNA/cluster_0.5/D35.rds")
npc <- merge(npc24h,y=c(npc44h,npc46h,npc48h,npc50h,npc52h,npc54h,npc72h,WT,D6,D12,D18,D25,D35),
             add.cell.ids = c("npc24h","npc44h", "npc46h", "npc48h","npc50h","npc52h","npc54h","npc72h","WT","D6","D12","D18","D25","D35"),
             project = "NPC_RNA")
npc <- NormalizeData(npc, normalization.method = "LogNormalize", scale.factor = 10000)
npc <- FindVariableFeatures(npc, selection.method = "vst", nfeatures = 2000)
npc <- ScaleData(npc)
npc <- RunPCA(npc, features = VariableFeatures(object = npc))
npc <- RunUMAP(npc, dims = 1:30)
npc <- FindNeighbors(npc, dims = 1:30)
npc <- FindClusters(npc, resolution = 0.5)
saveRDS(npc,"data/merge_data/scRNA/cluster_0.5/npc.rds")
DimPlot(npc,label = T,repel = T,pt.size = 1.5, group.by = "orig.ident",raster = F)
ggsave("result/merge_data/annotation/merge_orig.ident_final.png",width = 15,height = 10)

#######################scATAC-seq merge####################################
npc44h <- readRDS("data/merge_data/scATAC/cluster_0.5/44h_atac.rds")
npc46h <- readRDS("data/merge_data/scATAC/cluster_0.5/46h_atac.rds")
npc48h <- readRDS("data/merge_data/scATAC/cluster_0.5/48h_atac.rds")
npc50h <- readRDS("data/merge_data/scATAC/cluster_0.5/50h_atac.rds")
npc52h <- readRDS("data/merge_data/scATAC/cluster_0.5/52h_atac.rds")
npc54h <- readRDS("data/merge_data/scATAC/cluster_0.5/54h_atac.rds")
npc24h <- readRDS("data/merge_data/scATAC/cluster_0.5/24h_atac.rds")
npc72h <- readRDS("data/merge_data/scATAC/cluster_0.5/72h_atac.rds")
WT <- readRDS("data/merge_data/scATAC/cluster_0.5/WT_atac.rds")
D6 <- readRDS("data/merge_data/scATAC/cluster_0.5/D6_atac.rds")
D12 <- readRDS("data/merge_data/scATAC/cluster_0.5/D12_atac.rds")
D18 <- readRDS("data/merge_data/scATAC/cluster_0.5/D18_atac.rds")
D25 <- readRDS("data/merge_data/scATAC/cluster_0.5/D25_atac.rds")
D35 <- readRDS("data/merge_data/scATAC/cluster_0.5/D35_atac.rds")
npc_atac <- merge(npc24h,y=c(npc44h,npc46h,npc48h,npc50h,npc52h,npc54h,npc72h,WT,D6,D12,D18,D25,D35),
             add.cell.ids = c("npc24h","npc44h", "npc46h", "npc48h","npc50h","npc52h","npc54h","npc72h","WT","D6","D12","D18","D25","D35"),
             project = "NPC_ATAC")
npc_atac <- RunTFIDF(npc_atac)
npc_atac <- FindTopFeatures(npc_atac, min.cutoff = 'q0')
npc_atac <- RunSVD(npc_atac)
#npc <- RunUMAP(object = npc, reduction = 'lsi', dims = 2:30)
npc_atac <- FindNeighbors(object = npc_atac, reduction = 'lsi', dims = 2:30)
npc_atac <- FindClusters(object = npc_atac, verbose = TRUE, algorithm = 3,resolution = 0.5)
DimPlot(npc_atac,label = T,repel = T,pt.size = 1.5, group.by = "orig.ident",raster =F,label.size = 7)
ggsave("result/merge_data/annotation/merge_orig_atac_WT_D6_before_harmony.png",width = 15,height = 10)
npc_atac.harmony <- RunHarmony(object = npc_atac, 
                            group.by.vars = 'orig.ident', 
                            reduction = 'lsi', 
                            assay.use = 'peaks', 
                            project.dim = FALSE,
                            )
npc_atac.harmony <- RunUMAP(object = npc_atac.harmony, reduction = 'harmony', dims = 2:30)
DimPlot(npc_atac.harmony,label = T,repel = T,pt.size = 1.5, group.by = "orig.ident",raster = F,label.size = 7)
ggsave("result/merge_data/annotation/merge_orig_atac_WT_D6_after_harmony.png",width = 15,height = 10)
saveRDS(npc_atac,"data/merge_data/scATAC/cluster_0.5/npc_ATAC.rds")
saveRDS(npc_atac.harmony,"data/merge_data/scATAC/cluster_0.5/npc_ATAC_harmony.rds")
