#!/bin/bash

lists=(H9_WT NPC_24h NPC_44h NPC_46h NPC_48h NPC_50h NPC_52h NPC_54h NPC_72h D6_NPC NPC_D12 NPC_D18 NPC_D25 neuron_D35)
for list in ${lists[*]}
do
    echo ${list} begin
    cellranger-arc count --id=${list} \
                        --reference=/storage/zhangyanxiaoLab/suzhuojie/ref_data/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                        --libraries=/storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/data/merge_data/${list}_library.csv \
                        --localcores=16 \
                        --localmem=64
    echo ${list} done
done
