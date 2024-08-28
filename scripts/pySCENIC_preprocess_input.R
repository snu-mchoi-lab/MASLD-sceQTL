devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
library(SCopeLoomR)

input_path="path_to_pyscenic_input"

hep.cols=c('orig.ident',   "disease_group", 'disease_group2', 'cell_type_1', 'cell_type_4')
chol.cols=c('orig.ident', "disease_group", 'disease_group2', 'cell_type_1', 'RNA_snn_res.0.05')
hsc.cols=c('orig.ident', "disease_group", 'disease_group2', 'cell_type_1', 'cell_type_4')
endo.cols=c('orig.ident', "disease_group", 'disease_group2', 'cell_type_1', 'RNA_snn_res.0.12')
cols_list = list(hep.cols, chol.cols, hsc.cols, endo.cols)

for(i in 1:4){
  ser.tmp= celltype_ser_list[[i]]
  celltype_cur= ser.tmp@meta.data[1,"cell_type_1"]
  cols_cur = cols_list[[i]]

  exprMat <- ser.tmp@assays$RNA@data
  exprMat_filter <- exprMat[which(rowSums(exprMat) > 1*.01*ncol(exprMat)), ]

  cellInfo <- ser.tmp@meta.data[,cols_cur]

  loom <- build_loom(paste0(input_path, celltype_cur, ".loom"), dgem=exprMat_filter)
  loom <- add_cell_annotation(loom, cellInfo, colnames(cellInfo))
  close_loom(loom)
}
