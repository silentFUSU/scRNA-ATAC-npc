set.seed(1) 
myPaths <- .libPaths()
new <- c('/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library')
myPaths <- c(myPaths, new) 
.libPaths(myPaths)
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/")
library(Seurat)
library(ggplot2)
library(stringr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
filenames=c("24h","44h","46h","48h","50h","52h","54h","72h","D12","D18","D25","D35","D6","WT")
#所有样本的集合进行细胞周期预测
npc <- readRDS("data/merge_data/scRNA/cluster_0.5/npc.rds")
npc.cellcycle <- CellCycleScoring(npc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
npc.cellcycle$Phase<-factor(npc.cellcycle$Phase,levels = c("G1","S","G2M"))
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5,raster = F,group.by = "Phase")
ggsave("result/merge_data/Phase/RNA_merge_phase.png",width = 15,height = 10)
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5, cells.highlight = rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$Phase=="G1")],cols.highlight="#F8766D",raster = F)
ggsave("result/merge_data/Phase/RNA_merge_phase_G1.png",width = 15,height = 10)
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5,cells.highlight = rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$Phase=="S")],cols.highlight = "#0cb702",raster = F)
ggsave("result/merge_data/Phase/RNA_merge_phase_S.png",width = 15,height = 10)
DimPlot(npc.cellcycle,label = T,repel = T,pt.size = 1.5,cells.highlight = rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$Phase=="G2M")],cols.highlight="#00a9ff",raster = F)
ggsave("result/merge_data/Phase/RNA_merge_phase_G2M.png",width = 15,height = 10)
#所有样本的细胞周期占比
npc.cellcycle$orig.ident <- factor(npc.cellcycle$orig.ident,levels = c("WT","24h","44h","46h","48h","50h","52h","54h","72h","D6","D12","D18","D25","D35"))
ggplot(npc.cellcycle@meta.data, aes(x = orig.ident, fill = Phase)) +
  geom_bar(width = 0.3, position = "fill")+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("result/merge_data/Phase/merge_phase_percentage.png",width = 15,height = 10)
#所有样本进行细胞周期预测
for(i in 1:14){
  filename <- filenames[i]
  pbmc <- readRDS(paste0("data/merge_data/scRNA/cluster_0.5/",filename,".rds"))
  pbmc$Phase<-"S"
  barcode_G1 <- as.data.frame(rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$orig.ident==filename & npc.cellcycle$Phase=="G1")])
  colnames(barcode_G1)[1]<-"barcode"
  barcode_G1$barcode <- str_split_fixed(barcode_G1$barcode, "_", 2)
  barcode_G1 <- barcode_G1$barcode[,2]
  pbmc$Phase[which(rownames(pbmc@meta.data) %in% barcode_G1)]<-"G1"
  barcode_G2M <- as.data.frame(rownames(npc.cellcycle@meta.data)[which(npc.cellcycle$orig.ident==filename & npc.cellcycle$Phase=="G2M")])
  colnames(barcode_G2M)[1]<-"barcode"
  barcode_G2M$barcode <- str_split_fixed(barcode_G2M$barcode, "_", 2)
  barcode_G2M <- barcode_G2M$barcode[,2]
  pbmc$Phase[which(rownames(pbmc@meta.data) %in% barcode_G2M)]<-"G2M"
  pbmc$Phase<-factor(pbmc$Phase,levels = c("G1","S","G2M"))
  DimPlot(pbmc,label = T,repel = T,pt.size = 1.5,group.by = "Phase")
  ggsave(paste0("result/merge_data/Phase/",filename,"_phase.png"),width = 15,height = 10)
  saveRDS(pbmc,paste0("data/merge_data/scRNA/cluster_0.5/",filename,"_cellcycle.rds"))
}
