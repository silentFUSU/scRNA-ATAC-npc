# npc-scRNA&ATAC-seq
一些在npc项目中，分析scRNA-seq与scATAC-seq数据时所用到的代码，以及相关注释

seurat 所计算出来的cellcycle，先对WT-D35的集合进行cell cycle的估计，然后在回溯到每个sample中。注意到这种方法，和单个sample进行的cell cycle预测会存在一定的出入，所以以这种方法作为标准。
