# cholangiocyte pseudotime analysis
library(slingshot)
library(tidyverse)
library(ggplot2)
library(Seurat)


chol_sling <- slingshot(Embeddings(celltype_ser_list[[2]], "umap"), clusterLabels=celltype_ser_list[[2]]$RNA_snn_res.0.05, start.clus=1, stretch=0) # start from bipotent progenitor

ksmooth_gene <- function(one_gene, pseudotime_mat, norm_mat){
  if( sum(rownames(pseudotime_mat)!= colnames(norm_mat)) !=0){
    print("ERROR! Cell names in pseudotime matrix and normalized matrix do not match!")
    return("error")
  }
  else{
    one_gene_idx_norm <- which(rownames(norm_mat) == one_gene)
    one_gene_df <- data.frame(gene_count = norm_mat[one_gene_idx_norm,],
                              pseudotime = pseudotime_mat[,"pseudotime"])
    one_gene_df <- one_gene_df[!is.na(one_gene_df$pseudotime),]
    bw <- bw.SJ(one_gene_df$pseudotime)
    ks <- as.data.frame(ksmooth(y=one_gene_df$gene_count,
                                x=one_gene_df$pseudotime,
                                bandwidth=bw*20))
    colnames(ks) <- c("x", one_gene)
    return(as.data.frame(ks))
  }
}


plot_ksmooth <- function(slingshot_res, ser_obj, gene_list, plot_nrow=2, lineage_n=1){
  pseudotime_mat = slingPseudotime(slingshot_res)
  colnames(pseudotime_mat)[lineage_n] <- "pseudotime"
  metadata = ser_obj@meta.data

  if( sum(rownames(pseudotime_mat)!= rownames(metadata)) !=0){
    print("ERROR! Cell names in pseudotime matrix and meta data do not match!")
    return("error")
  }

  for (one_gene in gene_list){
    if(one_gene == gene_list[1]){
      ksdf <- ksmooth_gene(one_gene, pseudotime_mat, GetAssayData(ser_obj, slot="data"))
    }
    else{
      ksdf2 <- ksmooth_gene(one_gene, pseudotime_mat, GetAssayData(ser_obj, slot="data"))
      ksdf <- merge(ksdf, ksdf2, by="x")
    }
  }
  ksdf.m <- melt(ksdf, id.vars="x")
  p <- ggplot(ksdf.m, aes(x=x, y=value)) + geom_line() + theme_bw() + xlab("pseudotime") + ylab("normalized_expression") + facet_wrap(~variable, nrow=plot_nrow, scales="free")
  return(p)

chol_gene_list = c("MUC6", "FGFR2", "FOXA1", "FGFR3", # progenitor
                   "HNF4A", "ASGR1","TF","AMBP", # hepatocyte
                   "KRT19", "KRT18", "MMP7", "CD24") # cholangiocyte
chol_p <- plot_ksmooth(chol_sling, celltype_ser_list[[2]], chol_gene_list, plot_nrow=3 )
ggsave(paste0(save.dir, "/cholangiocyte_lineage_progenitor_marker_dotplot.pdf"), chol_p, width=5.5, height=6)




# developmental pathway genes

chol.ser <- readRDS(paste0(save.dir, "chol.ser_scWGCNA_obj.rds"))

notch=c("DLL3","RBPJL","DTX2","CREBBP","CTBP1","CTBP2","DTX3L","PTCRA","JAG1","DTX1","DVL1","DVL2","DVL3","DTX3","EP300","SNW1","DTX4","NCSTN","KAT2A","DLL1","HDAC1","HDAC2","HES1","RBPJ","JAG2","HES5","LFNG","MFNG","NOTCH1","NOTCH2","NOTCH3","NOTCH4","APH1A","DLL4","MAML3","PSENEN","PSEN1","PSEN2","RFNG","ADAM17","MAML2","NUMB","KAT2B","NUMBL","CIR1","NCOR2","MAML1") # kegg
wnt=c("FRAT1","APC2","NFAT5","WIF1","FZD10","CHP1","CSNK1A1L","CREBBP","PRICKLE1","CSNK1A1","CSNK1E","CSNK2A1","CSNK2A2","CSNK2B","CTBP1","CTBP2","CTNNB1","PRICKLE2","DVL1","DVL2","DVL3","EP300","DKK1","DAAM1","PLCB1","FBXW11","FRAT2","DAAM2","FZD2","CACYBP","DKK4","DKK2","GSK3B","APC","JUN","RHOA","LRP6","LRP5","SMAD2","SMAD3","SMAD4","MMP7","MYC","NFATC1","NFATC2","NFATC3","NFATC4","LEF1","WNT16","NLK","PLCB2","PLCB3","PLCB4","WNT4","PPARD","PPP2CA","PPP2CB","PPP2R1A","PPP2R1B","PPP2R5A","PPP2R5B","PPP2R5C","PPP2R5D","PPP2R5E","PPP3CA","PPP3CB","PPP3CC","PPP3R1","PPP3R2","PRKACA","PRKACB","PRKACG","PRKCA","PRKCB","PRKCG","MAPK8","MAPK9","MAPK10","PRKX","PSEN1","CTNNBIP1","VANGL2","CHD8","RAC1","RAC2","RAC3","SENP2","CCND1","ROCK1","CHP2","SFRP1","SFRP2","SFRP4","SFRP5","SOX17","SIAH1","PORCN","SKP1","MAP3K7","TBL1X","TCF7","TCF7L2","TP53","SKP1P2","WNT1","WNT2","WNT3","WNT5A","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT10B","WNT11","WNT2B","WNT9A","WNT9B","FZD5","TBL1XR1","FZD3","CXXC4","WNT10A","FOSL1","WNT5B","CAMK2A","CAMK2B","CAMK2D","CAMK2G","VANGL1","AXIN1","AXIN2","FZD1","FZD4","FZD6","FZD7","FZD8","FZD9","TCF7L1","CUL1","NKD1","NKD2","RUVBL1","CCND2","BTRC","CCND3","WNT3A","TBL1Y","CER1","ROCK2","RBX1") # kegg
hedgehog = c("CSNK1A1L","CSNK1A1","CSNK1D","CSNK1E","CSNK1G2","CSNK1G3","FBXW11","GAS1","STK36","GLI1","GLI2","GLI3","GSK3B","BMP8A","IHH","LRP2","DHH","WNT16","SUFU","RAB23","CSNK1G1","WNT4","PRKACA","PRKACB","PRKACG","PRKX","PTCH1","HHIP","SHH","BMP2","BMP4","BMP5","BMP6","BMP7","BMP8B","SMO","WNT1","WNT2","WNT3","WNT5A","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT10B","WNT11","WNT2B","WNT9A","WNT9B","ZIC2","WNT10A","WNT5B","PTCH2","BTRC","WNT3A") #kegg
hippo=c("CRB2","CRB1","PARD3","PARD6A","PARD6G","PARD6B","PRKCZ","PRKCI","PATJ","PALS1","AMOT","YAP1","WWTR1","CDH1","LIMD1","AJUBA","WTIP","NF2","WWC1","FRMD1","FRMD6","SAV1","STK3","RASSF6","RASSF1","PPP2CA","PPP2CB","PPP2R1B","PPP2R1A","PPP2R2A","PPP2R2B","PPP2R2C","PPP2R2D","LATS2","LATS1","MOB1A","MOB1B","PPP1CA","PPP1CB","PPP1CC","TP53BP2","LLGL2","LLGL1","SCRIB","DLG1","DLG2","DLG3","DLG4","DLG5","CSNK1D","CSNK1E","TPTEP2-CSNK1E","BTRC","FBXW11","TP73","BBC3","TEAD1","TEAD4","TEAD3","TEAD2","CCN2","GLI2","AREG","BIRC5","AFP","ITGB2","FGF1","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD7","SMAD2","SMAD3","SMAD4","SERPINE1","BMP2","BMP4","BMP5","BMP6","BMP7","BMP8B","BMP8A","GDF5","GDF6","GDF7","AMH","BMPR1A","BMPR1B","BMPR2","SMAD1","ID1","ID2","WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","WNT10B","WNT10A","WNT11","WNT16","FZD1","FZD7","FZD2","FZD3","FZD4","FZD5","FZD8","FZD6","FZD10","FZD9","DVL3","DVL2","DVL1","YWHAZ","YWHAB","YWHAQ","YWHAE","YWHAH","YWHAG","GSK3B","CTNNB1","APC","APC2","AXIN1","AXIN2","NKD1","NKD2","TCF7","TCF7L1","TCF7L2","LEF1","MYC","CCND1","CCND2","CCND3","SOX2","SNAI2","BIRC2","BIRC3","ACTG1","ACTB","CTNNA3","CTNNA1","CTNNA2") # Kegg manual

chol.ser <- AddModuleScore(chol.ser, features=list(notch), name="NOTCH")
chol.ser <- AddModuleScore(chol.ser, features=list(wnt), name="WNT")
chol.ser <- AddModuleScore(chol.ser, features=list(hedgehog), name="HH")
chol.ser <- AddModuleScore(chol.ser, features=list(hippo), name="HIPPO")


module_develop_cor = (cor(chol.ser@meta.data[,c("HH1", "HIPPO1","NOTCH1", 'WNT1', paste0("Chol-M", 1:7))]))[paste0('Chol-M', c(1,4,5,3,6,2,7)), 1:4]
write.csv(module_develop_cor, paste0(save.dir, "cholangiocytes_developmental_pathway_cor.csv"))
colnames(module_develop_cor) = str_replace(colnames(module_develop_cor), "1","")
rownames(module_develop_cor) = str_replace(rownames(module_develop_cor), "Chol-","")

dev.off()
pdf(paste0(save.dir, 'cholangiocytes_developmental_pathway_cor.pdf'), width=2, height=2.6)
module_develop_cor %>% pheatmap(., cluster_rows = F, treeheight_col = 5,
                                color = colorRampPalette(c(rev(brewer.pal(n= 7, name="PuBu")),
                                                                                brewer.pal(n = 7, name ="RdPu")))(100) )
dev.off()
