set.seed(1)  
myPaths <- .libPaths()
new <- c('/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library')
myPaths <- c(myPaths, new) 
.libPaths(myPaths)
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/")
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(tidyr)

npc <- readRDS("data/merge_data/scRNA/cluster_0.5/npc_cellcycle.rds")
samples <- c("WT","24h","44h","46h","48h","50h","52h","54h","72h","D6")
stage <- c("G1","S","G2M")
GO_database <- 'org.Hs.eg.db'
for(i in 2:10){
  sample1 <- samples[i]
  sample2 <- samples[i-1]
  for (j in 1:3){  
  marker <-  FindMarkers(npc, ident.1 = rownames(npc@meta.data)[which(npc$Phase==stage[j] & npc$orig.ident==sample1)],
                                          ident.2 = rownames(npc@meta.data)[which(npc$Phase==stage[j] & npc$orig.ident==sample2)])
  marker$Gene <- rownames(marker)
  marker$Significant <- ifelse(marker$p_val_adj < 0.05 & abs(marker$avg_log2FC) >= 0.7, 
                               ifelse(marker$avg_log2FC > 0.7, "Up", "Down"), "Stable")
  ggplot(
    # 数据、映射、颜色
    marker, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Significant), size=2) +
    scale_color_manual(values = c("blue","grey","red")) +
    # 注释
    geom_text_repel(
      data = subset(marker, p_val_adj < 0.05 & abs(marker$avg_log2FC) >= 0.7),
      aes(label = Gene),
      size = 5,max.overlaps = 100,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) +
    # 辅助线
    geom_vline(xintercept=c(-0.7,0.7),lty=4,col="black",lwd=0.8) +
    # geom_vline(xintercept=c(-1,1),lty=4,col="red",lwd=0.8)+
    geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
    # 坐标轴
    labs(x="log2(fold change)",
         y="-log10 (p-value)") +
    xlim(-3,3)+
    # 图例
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()

  ggsave(paste0("result/merge_data/DEA/",sample2,"_",sample1,"_",stage[j],"_volcano.png"),width = 13,height = 10)
  marker_gene <- bitr(marker$Gene[which(marker$avg_log2FC>0)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  marker_GO <- enrichGO( marker_gene$ENTREZID,#GO富集分析
                         OrgDb = GO_database,
                         keyType = "ENTREZID",#设定读取的gene ID类型
                         ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                         pvalueCutoff = 0.05,#设定p值阈值
                         qvalueCutoff = 0.05,#设定q值阈值
                         readable = T)
  if(nrow(marker_GO) > 0){ 
     dotplot(marker_GO)
     ggsave(paste0("result/merge_data/Go_KEGG/",sample2,"_",sample1,"_",stage[j],"_GO_up.png"),width = 10,height = 5)
     gene_list <- marker_GO@result
     gene_list <- gene_list[1:10,]
     gene_list<- separate(gene_list,geneID,paste0("gene",1:max(marker_GO$Count)),"/")
     write.csv(gene_list,paste0("result/merge_data/Go_KEGG/",sample2,"_",sample1,"_",stage[j],"_GO_up.csv"))
     }
  
  
  marker_gene <- bitr(marker$Gene[which(marker$avg_log2FC<0)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  marker_GO <- enrichGO( marker_gene$ENTREZID,#GO富集分析
                         OrgDb = GO_database,
                         keyType = "ENTREZID",#设定读取的gene ID类型
                         ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                         pvalueCutoff = 0.05,#设定p值阈值
                         qvalueCutoff = 0.05,#设定q值阈值
                         readable = T)
  if(nrow(marker_GO) > 0){ 
     dotplot(marker_GO)
     ggsave(paste0("result/merge_data/Go_KEGG/",sample2,"_",sample1,"_",stage[j],"_GO_down.png"),width = 10,height = 5)
     gene_list <- marker_GO@result
     gene_list <- gene_list[1:10,]
     gene_list<- separate(gene_list,geneID,paste0("gene",1:max(marker_GO$Count)),"/")
     write.csv(gene_list,paste0("result/merge_data/Go_KEGG/",sample2,"_",sample1,"_",stage[j],"_GO_down.csv"))
     }
  }
}
