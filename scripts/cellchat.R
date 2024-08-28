# cellchat
library(CellChat)
library(ComplexHeatmap)

cellchat <- createCellChat(object=ser, group.by='cell_type_2')
cellchat.ctrl.nafl <- createCellChat(object=subset(ser, disease_group2 %in% c("NAFL", "ctrl")), group.by='cell_type_2')
cellchat.nash <- createCellChat(object=subset(ser, disease_group2 %in% c("eNASH","aNASH")), group.by='cell_type_2')

cellchat.ctrl.nafl <- preprocess_cellchat(cellchat.ctrl.nafl)
cellchat.nash <- preprocess_cellchat(cellchat.nash)
cellchat <- preprocess_cellchat(cellchat)

object.list=list(ctrl_nafl = cellchat.ctrl.nafl, NASH = cellchat.nash)
cellchat.m <- mergeCellChat(object.list,
                            add.names = names(object.list))

nmf.out.ctrl.nafl <- selectK(object.list$ctrl_nafl, pattern='outgoing')
nmf.out.nash <- selectK(object.list$NASH, pattern='outgoing')

nmf.in.ctrl.nafl <- selectK(object.list$ctrl_nafl, pattern='incoming')
nmf.in.nash <- selectK(object.list$NASH, pattern='incoming')

groupSize <- as.numeric(table(cellchat@idents))
pdf(paste0(save.dir, "./all.ser.interaction.count.wieght.pdf"), width=15, height=8)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


df.net$pathway_name %>% unique()
pathways.show.list= c("TGFb", "BMP", "NRG", "PDGF", "VEGF", "IGF", "SPP1", "COLLAGEN", "LAMININ", "FN1", "CD45", "CD22", "VCAM", "SEMA3", "COMPLEMENT", "PARs", "VISFATIN", "ALCAM")
pathways.show.list %in% (df.net$pathway_name %>% unique())
for(i in 1:length(pathways.show.list)){
  pathway.show=pathways.show.list[i]
  pdf(paste0(save.dir, pathway.show, "_celltypeContribution.pdf"), width=4.5, height=2.5)
  netAnalysis_signalingRole_network(cellchat, signaling = pathway.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
}



cellchat.m <- computeNetSimilarityPairwise(cellchat.m, type = "structural")
cellchat.m <- netEmbedding(cellchat.m, type = "structural", umap.method = 'uwot')
cellchat.m <- netClustering(cellchat.m, type = "structural")

gg1 <- rankNet(cellchat.m, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.m, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
ggsave(paste0(save.dir, "ctrl.nash.info_flow_significant_DEP.pdf"), gg1+gg2, width=10, height=10)



pathways.show = c("NOTCH")
pdf(paste0(save.dir, pathways.show, "_NASH_celltypeContribution.pdf"), width=4.5, height=2.5)
netAnalysis_signalingRole_network(object.list[[2]], signaling = pathways.show, width = 8, height = 1.6, font.size = 10,
                                  color.heatmap="YlOrRd",color.use=full_celltype_color_list[levels(object.list[[1]]@idents)],
                                  measure.name=c("Sender", "Receiver", "Influencer"),
                                  measure=c("outdeg", 'indeg', 'info'))
dev.off()

pathways.show = c("CD6")
pdf(paste0(save.dir, pathways.show, "_NASH_celltypeContribution.pdf"), width=4.5, height=2.5)
netAnalysis_signalingRole_network(object.list[[2]], signaling = pathways.show, width = 8, height = 1.6, font.size = 10,
                                  color.heatmap="YlOrRd",color.use=full_celltype_color_list[levels(object.list[[1]]@idents)],
                                  measure.name=c("Sender", "Receiver", "Influencer"),
                                  measure=c("outdeg", 'indeg', 'info'))
dev.off()


pathways.show=c("IGF")
pdf(paste0(save.dir, pathways.show, "_ctrl_celltypeContribution.pdf"), width=4.5, height=2.5)
netAnalysis_signalingRole_network(object.list[[1]], signaling = pathways.show, width = 8, height = 1.6,
                                  font.size = 10, color.heatmap="YlGn",color.use=full_celltype_color_list[levels(object.list[[1]]@idents)],
                                  measure.name=c("Sender", "Receiver", "Influencer"),
                                  measure=c("outdeg", 'indeg', 'info'))

dev.off()
