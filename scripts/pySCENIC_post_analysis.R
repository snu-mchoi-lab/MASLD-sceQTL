# post-analysis of pyscenic outputs
# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html#retrieve-auc-scores-per-cell
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(grid)
library(cowplot)


output_dir='path_to_output'
celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')
for(celltype_cur in celltypes){
  loom_tmp <- open_loom(paste0("./seurat_analysis/211005_n48_A02/SCENIC/v2/output/", celltype_cur, "_scenic_output.loom"))
  regulons_incidMat <- get_regulons(loom_tmp, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom_tmp, column.attr.name='RegulonsAUC')

saveRDS(list('regulons' = regulons, 'regulonAUC' = regulonAUC),
        file = paste0(output_dir, celltype_cur, '_pyscenic_output_AUC.rds'))
}



scenic_heatmap_preprocess <- function(celltype_cur, ser_cur, scenic_dir=output_dir){
  scenic_cur = readRDS(paste0(scenic_dir,celltype_cur, "_pyscenic_output_AUC.rds" ))
  regulons_cur = scenic_cur$regulons
  regulonAUC_cur <- scenic_cur$regulonAUC
  AUCmat_cur = AUCell::getAUC(regulonAUC_cur)
  rownames(AUCmat_cur) <- gsub('[(+)]', '', rownames(AUCmat_cur))

  ser_new = ser_cur
  ser_new[["AUC"]] <- CreateAssayObject(data=AUCmat_cur)

  ser_new <- ScaleData(ser_new, assay='AUC', features=rownames(AUCmat_cur))

  return(list('scenic' = scenic_cur,
              'regulons' = regulons_cur,
              'regulonAUC' =regulonAUC_cur,
              'AUCmat' = AUCmat_cur,
              "ser_tmp" = ser_new) )
}


scenic_heatmap_raw_data <- function(prepared_list, grouping_colname, grouping_order=NA, save_heatmap=F, write_deg=F, write_auc=F, save_dir=output_dir){
  DefaultAssay(prepared_list$ser_tmp  ) <- "AUC"
  Idents(prepared_list$ser_tmp) <- grouping_colname
  celltype_cur = prepared_list$ser_tmp@meta.data[1,"cell_type_1"]

  deg.ls <- FindAllMarkers(prepared_list$ser_tmp, logfc.threshold = 0.005)
  deg.ls = deg.ls %>% filter(p_val_adj < 0.05)

  if(nrow(deg.ls) <= 1){ stop( "Less than 2 DEGs!. Consider using different grouping column ")}
  deg.ls = deg.ls %>% filter(avg_log2FC > 0) %>% arrange(cluster) %>% arrange(-avg_log2FC)

  if( !all(colnames(prepared_list$ser_tmp) == colnames(prepared_list$ser_tmp@assays$AUC@data) )){stop("AUC matrix in seurat object does not match with metadata in seurat object")}

  deg_mean_auc = prepared_list$ser_tmp@assays$AUC@data[deg.ls$gene %>% unique,] %>% t() %>% data.frame() %>%
    mutate( !!grouping_colname := prepared_list$ser_tmp@meta.data[,grouping_colname] ) %>%
    group_by(!!sym(grouping_colname)) %>%
    summarise_all(mean) %>% column_to_rownames(grouping_colname) %>% as.matrix() %>% scale(., center=T, scale=T) %>% t()

  if( length(grouping_order) > 1){
    deg_mean_auc = deg_mean_auc[,grouping_order]
    }

  if(save_heatmap){
    grouping_col_n = Idents(prepared_list$ser_tmp) %>% unique() %>% length()
    pdf(paste0(save_dir, "heatmap_", celltype_cur,"_", grouping_colname,".pdf"), width=grouping_col_n/2+2, nrow(deg.ls)/8)
    draw(pheatmap(deg_mean_auc, cluster_cols=F))
    dev.off()
  }

  if(write_deg){
    write.csv(deg.ls, paste0(save_dir, "DE_TF_", celltype_cur,"_", grouping_colname, ".csv"))
  }
  if(write_auc){
    write.csv(prepared_list$AUCmat, paste0(save_dir, 'auc_',celltype_cur,"_", grouping_colname, ".csv" ))
  }
  return(list( "deg_mean_auc" = deg_mean_auc,
               "deg_list" = deg.ls))
}


tmp <- scenic_heatmap_preprocess('cholangiocyte', celltype_ser_list[[2]])
tmp2 <- scenic_heatmap_raw_data(tmp, "disease_group2", grouping_order=c('ctrl', "NAFL", "eNASH", "aNASH"), save_heatmap=T, write_deg=T, write_auc=T)

tmp <- scenic_heatmap_preprocess('hepatocyte', celltype_ser_list[[1]])
tmp2 <- scenic_heatmap_raw_data(tmp, "disease_group2", grouping_order=c('ctrl', "NAFL", "eNASH", "aNASH"), save_heatmap=T, write_deg=T, write_auc=T)

tmp <- scenic_heatmap_preprocess('stellate_cell', celltype_ser_list[[3]])
tmp2 <- scenic_heatmap_raw_data(tmp, "disease_group2", grouping_order=c('ctrl', "NAFL", "eNASH", "aNASH"), save_heatmap=T, write_deg=T, write_auc=T)

tmp <- scenic_heatmap_preprocess('endothelial_cell', celltype_ser_list[[4]])
tmp2 <- scenic_heatmap_raw_data(tmp, "disease_group2", grouping_order=c('ctrl', "NAFL", "eNASH", "aNASH"), save_heatmap=T, write_deg=T, write_auc=T)


# correlation with module expression
celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')
wgcna.module.name=list("hepatocyte" = 'Hep-M', 'cholangiocyte'= "Chol-M", "stellate_cell" = "HSC-M", 'endothelial_cell' = "endo-M")
wgcna.module.count = list('hepatocyte'=c(3,4,6,8,12,14), 'cholangiocyte'=c(1,2,3,4,6,7), 'stellate_cell' = 1:3, 'endothelial_cell' = c(1,4,7,9))
wgcna.filename.prefix= list('hepatocyte' = "hep", 'cholangiocyte' = 'chol', 'stellate_cell'='hsc','endothelial_cell'='endo' )
scwgcna.dir="path_to_WGCNA_module_data"

for(celltype_cur in celltypes){
  hME <- readRDS(paste0(scwgcna.dir,wgcna.filename.prefix[[celltype_cur]], '.ser_hMEs.rds'))
  modules <- readRDS(paste0(scwgcna.dir,wgcna.filename.prefix[[celltype_cur]], '.ser_modules.rds'))[,c('module','color')] %>% unique
  wgcna.name.color <- split(modules$module, modules$color)
  colnames(hME) <- lapply(colnames(hME), function(x){wgcna.name.color[[x]]}) %>% unlist()
  hME <- hME[,c('grey', paste0(wgcna.module.name[[celltype_cur]], wgcna.module.count[[celltype_cur]]) )]

  auc_mat = read.csv(paste0(output_dir, celltype_cur, "_disease_group2.csv"), check.names = F)
  colnames(auc_mat)[1]="X"
  auc_mat = auc_mat %>% column_to_rownames("X") %>% t() %>% data.frame(., check.names=F)
  if(all(rownames(auc_mat) == rownames(hME)) ){
    auc_modules_merged <- cbind(auc_mat, hME)
  }
  else{ stop("Cell order does not match between auc_mat and hME!")}

  # calculate correlation coefficients, p values
  # https://stackoverflow.com/questions/34326906/extracting-and-formatting-results-of-cor-test-on-multiple-pairs-of-columns
  corrFunc <- function(var1, var2, data) {
    result = cor.test(data[,var1], data[,var2])
    data.frame(var1, var2, result[c("estimate","p.value","statistic","method")],
               stringsAsFactors=FALSE)
  }

  vars = data.frame('var1' = rep( colnames(auc_mat), each=ncol(hME)),
                    'var2' = rep( colnames(hME), times=ncol(auc_mat)))
  corrs = do.call(rbind, mapply(corrFunc, vars[,1], vars[,2], MoreArgs=list(data=auc_modules_merged),
                                SIMPLIFY=FALSE))
  corrs$p_adj_bf = p.adjust(corrs$p.value, method='bonferroni')


  write.csv(corrs, paste0(output_dir, celltype_cur, "_module_scenic_corr.csv"))
  # fill 0 values with 1e-301
  corrs[which(corrs$p.value < 1e-300),"p.value"] = 1e-301

}
