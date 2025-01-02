## hepatocyte: 211005_n48_A03_4_5_13_Reanalyze_pseudotime_scWGCNA
## others: 211005_n48_A02_S13_AWS_scWGCNA_v2
library(ggplot2)
library(Seurat)
library(cowplot)
library(Matrix)
library(reshape2)
library(gridExtra)
library(grid)
library(SeuratWrappers)
library(tidyverse)
library("WGCNA")
library("hdWGCNA")
library(patchwork)
library(igraph)
theme_set(theme_cowplot())

###############
# hepatocytes #
###############
hep.ser <- celltype_ser_list[[1]]

set.seed(12345)
hep.ser <- SetupForWGCNA(hep.ser, gene_select="fraction", fraction=0.05, wgcna_name="hep_wgcna")
hep.ser <- MetacellsByGroups(hep.ser, group.by=c("orig.ident", "predicted.celltype.l2", "disease_group2"), k=25, ident.group="disease_group2") #
hep.ser <- NormalizeMetacells(hep.ser)

hep.ser <- SetDatExpr(hep.ser, group.by="predicted.celltype.l2", group_name = "hepatocyte", use_metacells=TRUE)
# select soft-power threshold
hep.ser <- TestSoftPowers(hep.ser, setDatExpr=FALSE)
plot_list <- PlotSoftPowers(hep.ser)
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(hep.ser)
head(power_table)


# construct co-expression network
hep.ser <- ConstructNetwork(hep.ser, setDatExpr=FALSE) # automatic detection of softpowers

png("./scWGCNA_results/Hepatocyte_dendogram.png", width=3000, height=3000, res=300)
PlotDendrogram(hep.ser, main="Hepatocyte scWGCNA Dendogram (v3)")
dev.off()

# Module eigengenes and intramodular connectivity
hep.ser <- Seurat::ScaleData(hep.ser, features=GetWGCNAGenes(hep.ser), vars.to.regress = c('nCount_RNA'))
hep.ser <- ModuleEigengenes(hep.ser, group.by.vars="orig.ident")
hMEs <- GetMEs(hep.ser)
MEs <- GetMEs(hep.ser, harmonized=FALSE)

# compute module connectivity
hep.ser <- ModuleConnectivity(hep.ser, group.by="predicted.celltype.l2", group_name = "hepatocyte")
hep.ser <- ResetModuleNames(hep.ser, new_name="Hep-M")
modules <- GetModules(hep.ser)
hep.ser <- ModuleExprScore(hep.ser, n_genes=25, method ="Seurat")

hep.ser@meta.data <- cbind(
  hep.ser@meta.data,
  GetMEs(hep.ser, harmonized=TRUE)
)

saveRDS(modules, "./scWGCNA_results/hep.ser_modules.rds")
# write.csv(modules, "./scWGCNA_results/hep.ser_modules.csv")
saveRDS(MEs, "./scWGCNA_results/hep.ser_MEs.rds")
# write.csv(MEs, "./scWGCNA_results/hep.ser_MEs.csv")
saveRDS(hMEs, "./scWGCNA_results/hep.ser_hMEs.rds")
# write.csv(hMEs, "./scWGCNA_results/hep.ser_hMEs.csv")
saveRDS(hep.ser, file="./scWGCNA_results/hep.ser_scWGCNA_obj.rds")
saveRDS(GetMetaCellObject(hep.ser), file='./scWGCNA_results/hep.ser_scWGCNA_MetaCell_obj.rds')
write.csv(hep.ser@meta.data, file='./scWGCNA_results/hep.ser_metaData.csv')



##################
# cholangiocytes #
##################
set.seed(12345)
chol.ser <- celltype_ser_list[[2]]
chol.ser <- SetupForWGCNA(chol.ser, gene_select="fraction", fraction=0.05, wgcna_name="chol_wgcna")
chol.ser <- MetacellsByGroups(chol.ser, group.by=c("orig.ident", "predicted.celltype.l2", "disease_group2"), k=25, ident.group="disease_group2") #
chol.ser <- NormalizeMetacells(chol.ser)

chol.ser <- SetDatExpr(chol.ser, group.by="predicted.celltype.l2", group_name = "cholangiocyte", use_metacells=TRUE)
# select soft-power threshold
chol.ser <- TestSoftPowers(chol.ser, setDatExpr=FALSE)
plot_list <- PlotSoftPowers(chol.ser)
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(chol.ser)

# construct co-expression network
chol.ser <- ConstructNetwork(chol.ser, setDatExpr=FALSE) # automatic detection of softpowers
png("./scWGCNA_results/Cholangiocyte_dendogram.png", width=2000, height=2000, res=300)
PlotDendrogram(chol.ser, main="Cholangiocyte scWGCNA Dendogram")
dev.off()

# Module eigengenes and intramodular connectivity
chol.ser <- Seurat::ScaleData(chol.ser, features=GetWGCNAGenes(chol.ser), vars.to.regress = c('nCount_RNA'))
chol.ser <- ModuleEigengenes(chol.ser, group.by.vars="orig.ident")
hMEs <- GetMEs(chol.ser)
MEs <- GetMEs(chol.ser, harmonized=FALSE)

# compute module connectivity
chol.ser <- ModuleConnectivity(chol.ser, group.by="predicted.celltype.l2", group_name = "cholangiocyte")
chol.ser <- ResetModuleNames(chol.ser, new_name="Chol-M")
modules <- GetModules(chol.ser)

# add hMEs to Seurat meta-data:
chol.ser@meta.data <- cbind(
  chol.ser@meta.data,
  GetMEs(chol.ser, harmonized=TRUE)
)


saveRDS(modules, "./scWGCNA_results/chol.ser_modules.rds")
# write.csv(modules, "./scWGCNA_results/chol.ser_modules.csv")
saveRDS(MEs, "./scWGCNA_results/chol.ser_MEs.rds")
# write.csv(MEs, "./scWGCNA_results/chol.ser_MEs.csv")
saveRDS(hMEs, "./scWGCNA_results/chol.ser_hMEs.rds")
# write.csv(hMEs, "./scWGCNA_results/chol.ser_hMEs.csv")
saveRDS(chol.ser, file="./scWGCNA_results/chol.ser_scWGCNA_obj.rds")
saveRDS(GetMetaCellObject(chol.ser), file='./scWGCNA_results/chol.ser_scWGCNA_MetaCell_obj.rds')
write.csv(chol.ser@meta.data, file='./scWGCNA_results/chol.ser_metaData.csv')



#################
# stellate cell #
#################
hsc.ser <- celltype_ser_list[[3]]

set.seed(12345)
hsc.ser <- SetupForWGCNA(hsc.ser, gene_select="fraction", fraction=0.05, wgcna_name="hsc_wgcna")
hsc.ser <- MetacellsByGroups(hsc.ser, group.by=c("orig.ident", "predicted.celltype.l2", "disease_group2"), k=25, ident.group="disease_group2") #
hsc.ser <- NormalizeMetacells(hsc.ser)
hsc.ser <- SetDatExpr(hsc.ser, group.by="predicted.celltype.l2", group_name = "stellate_cell", use_metacells=TRUE)
# select soft-power threshold
hsc.ser <- TestSoftPowers(hsc.ser, setDatExpr=FALSE)
plot_list <- PlotSoftPowers(hsc.ser)
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(hsc.ser)
# construct co-expression network
hsc.ser <- ConstructNetwork(hsc.ser, setDatExpr=FALSE) # automatic detection of softpowers
png("./scWGCNA_results/Stellate_cell_dendogram.png", width=2000, height=2000, res=300)
PlotDendrogram(hsc.ser, main="Stellate cell scWGCNA Dendogram")
dev.off()

# Module eigengenes and intramodular connectivity
hsc.ser <- Seurat::ScaleData(hsc.ser, features=GetWGCNAGenes(hsc.ser), vars.to.regress = c('nCount_RNA'))
hsc.ser <- ModuleEigengenes(hsc.ser, group.by.vars="orig.ident")
hMEs <- GetMEs(hsc.ser)
MEs <- GetMEs(hsc.ser, harmonized=FALSE)

# compute module connectivity
hsc.ser <- ModuleConnectivity(hsc.ser, group.by="predicted.celltype.l2", group_name = "stellate_cell")
hsc.ser <- ResetModuleNames(hsc.ser, new_name="HSC-M")
modules <- GetModules(hsc.ser)

# add hMEs to Seurat meta-data:
hsc.ser@meta.data <- cbind(
  hsc.ser@meta.data,
  GetMEs(hsc.ser, harmonized=TRUE)
)


saveRDS(modules, "./scWGCNA_results/hsc.ser_modules.rds")
# write.csv(modules, "./scWGCNA_results/hsc.ser_modules.csv")
saveRDS(MEs, "./scWGCNA_results/hsc.ser_MEs.rds")
# write.csv(MEs, "./scWGCNA_results/hsc.ser_MEs.csv")
saveRDS(hMEs, "./scWGCNA_results/hsc.ser_hMEs.rds")
# write.csv(hMEs, "./scWGCNA_results/hsc.ser_hMEs.csv")
saveRDS(hsc.ser, file="./scWGCNA_results/hsc.ser_scWGCNA_obj.rds")
saveRDS(GetMetaCellObject(hsc.ser), file='./scWGCNA_results/hsc.ser_scWGCNA_MetaCell_obj.rds')
# write.csv(hsc.ser@meta.data, file='./scWGCNA_results/hsc.ser_metaData.csv')


####################
# endothelial cell #
####################
endo.ser <- celltype_ser_list[[4]]

set.seed(12345)
endo.ser <- SetupForWGCNA(endo.ser, gene_select="fraction", fraction=0.05, wgcna_name="endo_wgcna")
endo.ser <- MetacellsByGroups(endo.ser, group.by=c("orig.ident", "predicted.celltype.l2", "disease_group2"), k=25, ident.group="disease_group2") #
endo.ser <- NormalizeMetacells(endo.ser)
endo.ser <- SetDatExpr(endo.ser, group.by="predicted.celltype.l2", group_name = "endothelial_cell", use_metacells=TRUE)
# select soft-power threshold
endo.ser <- TestSoftPowers(endo.ser, setDatExpr=FALSE)
plot_list <- PlotSoftPowers(endo.ser)
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(endo.ser)

# construct co-expression network
endo.ser <- ConstructNetwork(endo.ser, setDatExpr=FALSE) # automatic detection of softpowers
png("./scWGCNA_results/Endothelial_cell_dendogram.png", width=2000, height=2000, res=300)
PlotDendrogram(endo.ser, main="Endothelial cell scWGCNA Dendogram")
dev.off()

# Module eigengenes and intramodular connectivity
endo.ser <- Seurat::ScaleData(endo.ser, features=GetWGCNAGenes(endo.ser), vars.to.regress = c('nCount_RNA'))
endo.ser <- ModuleEigengenes(endo.ser, group.by.vars="orig.ident")
hMEs <- GetMEs(endo.ser)
MEs <- GetMEs(endo.ser, harmonized=FALSE)

# compute module connectivity
endo.ser <- ModuleConnectivity(endo.ser, group.by="predicted.celltype.l2", group_name = "endothelial_cell")
endo.ser <- ResetModuleNames(endo.ser, new_name="endo-M")
modules <- GetModules(endo.ser)

# add hMEs to Seurat meta-data:
endo.ser@meta.data <- cbind(
  endo.ser@meta.data,
  GetMEs(endo.ser, harmonized=TRUE)
)

saveRDS(modules, "./scWGCNA_results/endo.ser_modules.rds")
#write.csv(modules, "./scWGCNA_results/endo.ser_modules.csv")
saveRDS(MEs, "./scWGCNA_results/endo.ser_MEs.rds")
#write.csv(MEs, "./scWGCNA_results/endo.ser_MEs.csv")
saveRDS(hMEs, "./scWGCNA_results/endo.ser_hMEs.rds")
#write.csv(hMEs, "./scWGCNA_results/endo.ser_hMEs.csv")
saveRDS(endo.ser, file="./scWGCNA_results/endo.ser_scWGCNA_obj.rds")
saveRDS(GetMetaCellObject(endo.ser), file='./scWGCNA_results/endo.ser_scWGCNA_MetaCell_obj.rds')
write.csv(endo.ser@meta.data, file='./scWGCNA_results/endo.ser_metaData.csv')





# hepatocyte module projection
# 211005_n48_A02_S13_scWGCNA.R : hep21 (V4), babyhep (nat22),  plotting
# module 8/14
hep.21.orig <- readRDS("../Integrated_with_Subannotations.hep.comm.2021.rds")
hep.21.nuc = subset(hep.21.orig, assay_type == "single_nuc")
hep.21.v4 <- ProjectModules(
  seurat_obj = hep.21.v3,
  seurat_ref = hep.ser,
  wgcna_name="hep_wgcna",
  wgcna_name_proj = "hep_wgcna_projected",
  assay="RNA"
)

hep.21.v4 <- ModuleConnectivity(hep.21.v3, assay="RNA", slot='data')

hep.21.v4@meta.data$cell_type_sh <- "Hep"
hep.21.v4 <- SetDatExpr(hep.21.v4, group_name="Hep", group.by="cell_type_sh",use_metacells = F)
hep.21.v4 <- ModulePreservation(hep.21.v4, hep.ser, name="sn-hep", verbose=3, n_permutations=250)
mod_pres <- GetModulePreservation(hep.21.v4, "sn-hep")$Z
obs_df <- GetModulePreservation(hep.21.v4, "sn-hep")$obs

plot_list <- PlotModulePreservation(
  hep.21.v4,
  name="sn-hep",
  statistics = "summary"
)
pdf("hep.21.v4.plotmodule.pdf", width=12, height=6)
patchwork::wrap_plots(plot_list, ncol=2)
dev.off()

write.csv(mod_pres, "./hep21.v4.modulePreservation_z.csv", row.names = F)
write.csv(obs_df, "./hep21.v4.modulePreservation_obs.csv", row.names = F)




# module 3,4
tmp <- read_h5ad("../Hepatocytes_primary_ad.nat.cell.biol.2022.h5ad")
nat.22 <- CreateSeuratObject(count=t(tmp$X %>% as.matrix()), data=t(tmp$X %>% as.matrix()), meta.data=tmp$obs)
umap <- tmp$obsm$X_umap
rownames(umap) <- rownames(tmp$X)
pca <- tmp$obsm$X_pca
rownames(pca) <- rownames(tmp$X)
diffmap <- tmp$obsm$X_diffmap
rownames(diffmap) <- rownames(tmp$X)
tsne <- tmp$obsm$X_tsne
rownames(tsne) <- rownames(tmp$X)

nat.22.pca <- CreateDimReducObject(embeddings=pca, key="PC_", assay="RNA", global=T)
nat.22.umap <- CreateDimReducObject(embeddings=umap, key="UMAP_", assay="RNA", global=T)
nat.22.diffmap <- CreateDimReducObject(embeddings=diffmap, key="DIFFMAP_", assay="RNA", global=T)
nat.22.tsne <- CreateDimReducObject(embeddings=tsne, key="TSNE_", assay="RNA", global=T)

nat.22[['pca']] <- nat.22.pca
nat.22[['umap']] <- nat.22.umap
nat.22[['diffmap']] <- nat.22.diffmap
nat.22[['tsne']] <- nat.22.tsne

nat.22 <- FindVariableFeatures(nat.22)
nat.22 <- ScaleData(nat.22)
DimPlot(nat.22, reduction = 'umap', group.by="source")

nat.22 <- ProjectModules(
  seurat_obj = nat.22,
  seurat_ref = hep.ser,
  wgcna_name="hep_wgcna",
  wgcna_name_proj = "hep_wgcna_projected",
  assay="RNA"
)

nat.22 <- ModuleExprScore(
  nat.22, n_genes=25, method="Seurat", nbin=15)

nat.22 <- ModuleConnectivity(nat.22, assay="RNA", slot='data')
nat.22@meta.data <- cbind(nat.22@meta.data, GetMEs(nat.22, harmonized=TRUE))

nat.22 <- SetDatExpr(nat.22, group_name="Hep", group.by="cell_type_sh",use_metacells = F)
nat.22 <- ModulePreservation(nat.22, hep.ser, name="sn-hep", verbose=3, n_permutations=250)

mod_pres <- GetModulePreservation(nat.22, "sn-hep")$Z
obs_df <- GetModulePreservation(nat.22, "sn-hep")$obs

write.csv(mod_pres, "./nat22.modulePreservation_z.csv", row.names = F)
write.csv(obs_df, "./nat22.modulePreservation_obs.csv", row.names = F)

plot_list <- PlotModulePreservation(
  nat.22,
  name="sn-hep",
  statistics = "summary"
)


pdf("./nat.22.plotmodule.pdf", width=12, height=6)
wrap_plots(plot_list, ncol=2)
dev.off()



# modules 6,12
stm.human.samples = c("GSM6808760_Y31122A7", "GSM6808758_Y21822A2", "GSM6808755_Y21822A1",
                      "GSM6808756_Y31122A5", "GSM6808759_Y31122A6", "GSM6808757_Y31122A8")
stm.human.samples.short1=c("NASH.A7", "NASH.A2", "Normal.A1", "Normal.A5", "NASH.A6", "Normal.A8")
stm.human.samples.short2=c("NASH3", "NASH1", "Normal1", "Normal2", "NASH2", "Normal3")

stm.data.dir="./public_data/EphB2_NASH_2023_SciTransMed/"
stm.ser.list = vector('list', 6)
for(i in 1:6){
  expr_mtx <- ReadMtx(mtx= paste0(stm.data.dir, stm.human.samples[i], "_matrix.mtx.gz"),
                      features=paste0(stm.data.dir, stm.human.samples[i], "_features.tsv.gz"),
                      cells=paste0(stm.data.dir, stm.human.samples[i], "_barcodes.tsv.gz"))
  stm.ser.list[[i]] <- CreateSeuratObject(counts=expr_mtx, project=stm.human.samples.short1[i])
}

stm.ser <- merge(stm.ser.list[[1]], y=unlist(stm.ser.list)[2:length(stm.ser.list)],
                    add.cell.ids= stm.human.samples.short2,
                    merge.data=TRUE)


stm.ser[['percent.mt']] <- PercentageFeatureSet(stm.ser, pattern="^MT-")
stm.ser = stm.ser %>% subset(., subset=nFeature_RNA>500 & nFeature_RNA<8000 & percent.mt <10) %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(., vars.to.regress='nCount_RNA') %>% RunPCA(pc.genes=stm.ser@var.genes, npcs=40)
stm.ser <- stm.ser %>% RunHarmony("orig.ident", plot_converge=T)
stm.ser <- stm.ser %>% RunUMAP(reduction='harmony', dims=1:40) %>% FindNeighbors(reduction='harmony', dims=1:40) %>% FindClusters(resolution=1.1)


# metadata from original article
metadata.orig <- read.csv(paste0(stm.data.dir, "scitranslmed.adc9653_data_file_s2_humanSeuratMetadata.txt"),sep='\t') %>% column_to_rownames("X")
stm.ser <- subset(stm.ser, cells=rownames(metadata.orig))
stm.ser.cellorder = rownames(stm.ser@meta.data)
newmetadata = (merge(stm.ser@meta.data, metadata.orig[,c( "CellCluster", 'condition')], by=0) %>% column_to_rownames("Row.names"))[stm.ser.cellorder,]
stm.ser@meta.data <- newmetadata

stm.ser.hep <- ProjectModules(
  seurat_obj = stm.ser.hep,
  seurat_ref = hep.ser,
  #group.by.vars = "orig.ident",
  modules=NULL,
  wgcna_name="hep_wgcna",
  wgcna_name_proj = "hep_wgcna_projected",
  assay="RNA",
  vars.to.regress = "nCount_RNA")

stm.ser.hep <- ModuleExprScore(stm.ser.hep, n_genes=25, method="Seurat", nbin=15)
stm.ser.hep <- ModuleConnectivity(stm.ser.hep, assay="RNA", slot='data')
stm.ser.hep@meta.data <- cbind(stm.ser.hep@meta.data, GetMEs(stm.ser.hep, harmonized=TRUE))

stm.ser.hep <- ModulePreservation(stm.ser.hep, hep.ser, name="sn-hep", verbose=3, n_permutations=250, parallel=T)
mod_pres <- GetModulePreservation(stm.ser.hep, "sn-hep")$Z
obs_df <- GetModulePreservation(stm.ser.hep, "sn-hep")$obs

plot_list <- PlotModulePreservation(
  stm.ser.hep,
  name="sn-hep",
  statistics = "summary"
)


pdf(paste0(save.dir, "stm.ser.hep.Preservation.plots.pdf"), width=10, height=5)
wrap_plots(plot_list, ncol=2)
dev.off()

write.csv(mod_pres, paste0(save.dir, "stm.ser.hep.modulePreservation_z.csv"), row.names = F)
write.csv(obs_df, paste0(save.dir, "stm.ser.hep.modulePreservation_obs.csv"), row.names = F)


# module overlap heatmap
celltype='hepatocyte'
module.info = read.table("./public_data/EphB2_NASH_2023_SciTransMed/scitranslmed.adc9653_data_file_s2_humanModules.txt", sep='\t', header=1)
colnames(module.info) = c("Module1_2", "Module3", "Module4")

wgcna.module <- readRDS(paste0(input.dir, "hep.ser_modules.rds"))# [,2:3] #%>% filter(!duplicated(module,color)) %>% remove_rownames()

public.kappa = data.frame( matrix(rep(NA, 6), ncol=6))
colnames(public.kappa) = c("wgcna_gene_count", "public_module_count", "intersect_n", "union_n", 'wgcna_module', "public_module")
for(i in 1:18){
  for (j in 1:3){
    wgcna.genes = (wgcna.module %>% filter(module == !!paste0("Hep-M", i)))$gene_name
    public.genes= module.info[,j,drop=T]
    public.genes = public.genes[which(public.genes != "")]
    public.kappa[1,] = c(length(wgcna.genes), length(public.genes),
                     length(intersect(wgcna.genes, public.genes)),
                     length(union(wgcna.genes, public.genes)),
                     paste0("Hep-M", i),
                     colnames(module.info)[j])
    if(i==1 & j==1){ public.kappa.m = public.kappa}
    else{ public.kappa.m = rbind(public.kappa.m, public.kappa)}
  }
}

public.kappa.m[,c("intersect_n", "union_n")] = sapply(public.kappa.m[,c("intersect_n", "union_n")], as.numeric)
public.kappa.m$kappa = public.kappa.m$intersect_n / public.kappa.m$union_n

pdf(paste0(save.dir, "EphB2_module_WGCNA_kappa_index.pdf"), width=4.8, height=3.15)
(public.kappa.m[,c("wgcna_module", 'public_module', 'kappa')] %>% reshape2::dcast(.,wgcna_module~public_module) %>% column_to_rownames("wgcna_module") )[c("Hep-M6", "Hep-M12"),] %>% t() %>%
  pheatmap(., cluster_rows = F, cluster_cols = F, main="Kappa index (EphB2)",
           color = colorRampPalette(c(rev(brewer.pal(n= 7, name="PuBu")),
                                      brewer.pal(n = 7, name ="RdPu")))(100) ,
           fontsize=15)
dev.off()
