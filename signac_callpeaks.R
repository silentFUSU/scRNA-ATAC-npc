set.seed(1) 
myPaths <- .libPaths()
new <- c('/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library')
myPaths <- c(myPaths, new) 
.libPaths(myPaths)
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(harmony)
filenames=c("WT","24h","44h","46h","48h","50h","52h","54h","72h","D6")
#rna seurat cellcyclescoring 计算出来的结果回溯到ATAC上，并callpeaks
for(i in 1:10){
  filename <- filenames[i]
  se <- readRDS(paste0("data/merge_data/scATAC/cluster_0.5/",filename,"_atac.rds"))
  rna_se <- readRDS(paste0("data/merge_data/scRNA/cluster_0.5/",filename,"_cellcycle.rds"))
  se$Phase <- rna_se$Phase
  DimPlot(se,group.by = "Phase")
  peaks <- CallPeaks(
    object = se,
    group.by = "Phase",
    macs2.path ="/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/bin/macs2",
    format = "BEDPE"
  )
  CoveragePlot(
    object = se,
    region = "PAX6",
    peaks.group.by="Phase",
    ranges = peaks,
    ranges.title = "MACS2"
  )
  ggsave(paste0("result/merge_data/signac/WT_D6_PAX6.png"),width = 5,height = 5)
  saveRDS(se,paste0("data/merge_data/scATAC/cluster_0.5/",filename,"_atac_cellcycle.rds"))
  saveRDS(peaks,"data/merge_data/scATAC/peaks/",filename,"_ATAC_cellcycle_peaks.rds")
}
#对WT-D6整体进行callpeaks
se <- readRDS("data/merge_data/scATAC/cluster_0.5/npc_ATAC.rds")
rna_se <- readRDS("data/merge_data/scRNA/cluster_0.5/npc_cellcycle.rds")
se$Phase <- rna_se$Phase
se <- subset(se,cells=rownames(se@meta.data)[which(se$orig.ident %in% c("WT","24h","44h","46h","48h","50h","52h","54h","72h","D6"))])
se <- RunTFIDF(se)
se <- FindTopFeatures(se, min.cutoff = 'q0')
se <- RunSVD(se)
se <- FindNeighbors(object = se, reduction = 'lsi', dims = 2:30)
se <- FindClusters(object = se, verbose = TRUE, algorithm = 3,resolution = 0.5)
DimPlot(se,group.by = "orig.ident",label=T)
saveRDS(se,"data/merge_data/scATAC/cluster_0.5/npc_WT_D6_ATAC_cellcycle.rds")
se<-readRDS("data/merge_data/scATAC/cluster_0.5/npc_WT_D6_ATAC_cellcycle.rds")
peaks_se <- CallPeaks(
  object = se,
  group.by = "orig.ident",
  macs2.path ="/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/bin/macs2",
  format = "BEDPE"
)
saveRDS(peaks_se,"data/merge_data/scATAC/peaks/npc_WT_D6_ATAC_cellcycle_peaks.rds")

#"chr11-31784779-31817961" PAX6范围
se$orig.ident<-factor(se$orig.ident,levels = c("WT","24h","44h","46h","48h","50h","52h","54h","72h","D6"))
ranges.show <- StringToGRanges(c("chr11-31795722-31796845","chr11-31804677-31805607","chr11-31809698-31810701",
        "chr11-31810835-31811854","chr11-31811872-31813306"))#pax6一些有明显变化的峰
CoveragePlot(
  object = se,
  region = c("chr11-31795700-31817961")#pax6的部分范围,
  group.by = "orig.ident",
  ranges = peaks_se,
  ranges.title = "MACS2",
  region.highlight=ranges.show
)

ranges<-as.data.frame(se@assays[["peaks"]]@ranges)
ranges<-ranges[which(ranges$seqnames=="chr11" & ranges$start>=31784779 & ranges$end<=31817961),]#选择得到pax6的范围
pax6range<-paste0(ranges$seqnames,"-",ranges$start,"-",ranges$end)
for(i in 2:10){
  filename1<-filenames[i]
  filename2<-filenames[i-1]
  da_peaks <- FindMarkers(
    object = se,
    ident.1 = rownames(se@meta.data)[which(se$orig.ident==filename1)],
    ident.2 = rownames(se@meta.data)[which(se$orig.ident==filename2)],
    test.use = 'LR',
    latent.vars = 'atac_peak_region_fragments',
    features = pax6range,min.pct = 0,logfc.threshold = 0)
  write.csv(da_peaks,paste0("result/merge_data/signac/",filename1,"_",filename2,"_PAX6_change.csv"))
  }#查看在PAX6上峰的差异
  

  
