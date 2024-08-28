
celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', "endothelial_cell")
max_chunks=c(19,17,14,15)
max_chunks = as.list(max_chunks)
names(max_chunks) = celltypes

input_results_dir="path_to_input_files"
#input_results_dir="./results/"
analysis_results_dir="path_to_output_files"
#analysis_results_dir="./results/"
input_pheno_dir = "path_to_other_input_files"

for(celltype_cur in c("endothelial_cell", "cholangiocyte", 'stellate_cell')){ #, "hepatocyte",
  for(i in 1:max_chunks[[celltype_cur]]){
    tmp <- fread(paste0(input_results_dir, "gt_liver_disease_", celltype_cur, "_chunk_", i, ".csv"))
    tmp = tmp %>% filter(term == "G")
    if(i==1) { pme.m = tmp}
    else{ pme.m = rbind(pme.m, tmp)}
  }

  input.genes = (fread(paste0(input_pheno_dir, "ge.df.", celltype_cur, ".txt"), nrows=2) %>% colnames())[-c(1:8)]
  tested.genes = pme.m$gene %>% unique
  missing.input.genes = input.genes[which(!input.genes %in% tested.genes)]
  print(paste0("Length of input genes: ", length(input.genes)))
  print(paste0("Length of result genes: ", length(tested.genes )))
  print(paste0("Following input genes are not included in result genes : ", paste0(missing.input.genes, collapse=",")))
  #print(paste0("cis-SNPs are: ", paste0(gene_snp_pairs[missing.input.genes], collapse=',')))
  pme.m$gene_snp = paste0(pme.m$gene, "_", pme.m$snp)
  pme.m = pme.m %>% group_by(gene) %>% mutate(padj_eSNP_null = p.adjust(pval_null, method="BH"),
                                              padj_eSNP_NAFLD = p.adjust(pval_NAFLD, method="BH"),
                                              padj_eSNP_ctrl = p.adjust(pval_ctrl, method='BH')) %>% ungroup() %>% data.frame()
  pme.m_bestP_null = pme.m %>% group_by(gene) %>% dplyr::slice(which.min(padj_eSNP_null))
  pme.m_bestP_NAFLD = pme.m %>% group_by(gene) %>% dplyr::slice(which.min(padj_eSNP_NAFLD))
  pme.m_bestP_ctrl = pme.m %>% group_by(gene) %>% dplyr::slice(which.min(padj_eSNP_ctrl))

  pme.m_bestP_null$padj_eSNP_eGene_null = p.adjust(pme.m_bestP_null$padj_eSNP_null, method="BH")
  pme.m_bestP_NAFLD$padj_eSNP_eGene_NAFLD = p.adjust(pme.m_bestP_NAFLD$padj_eSNP_NAFLD, method="BH")
  pme.m_bestP_ctrl$padj_eSNP_eGene_ctrl = p.adjust(pme.m_bestP_ctrl$padj_eSNP_ctrl, method="BH")

  max_p_0.05_thr_null = pme.m_bestP_null %>% filter(padj_eSNP_eGene_null < 0.05) %>% dplyr::slice(which.max(padj_eSNP_eGene_null)) %>% ungroup() %>% dplyr::select(padj_eSNP_null) %>% max()
  max_p_0.05_thr_NAFLD = pme.m_bestP_NAFLD %>% filter(padj_eSNP_eGene_NAFLD < 0.05) %>% dplyr::slice(which.max(padj_eSNP_eGene_NAFLD)) %>% ungroup() %>% dplyr::select(padj_eSNP_NAFLD) %>% max()
  max_p_0.05_thr_ctrl = pme.m_bestP_ctrl %>% filter(padj_eSNP_eGene_ctrl < 0.05) %>% dplyr::slice(which.max(padj_eSNP_eGene_ctrl)) %>% ungroup() %>% dplyr::select(padj_eSNP_ctrl) %>% max()

  pme.m.sig.null = pme.m %>% filter(padj_eSNP_null < max_p_0.05_thr_null)
  pme.m.sig.NAFLD = pme.m %>% filter(padj_eSNP_NAFLD < max_p_0.05_thr_NAFLD)
  pme.m.sig.ctrl = pme.m %>% filter(padj_eSNP_ctrl < max_p_0.05_thr_ctrl)

  pme.m.sig = pme.m %>% filter(gene_snp %in% c(pme.m.sig.null$gene_snp, pme.m.sig.NAFLD$gene_snp, pme.m.sig.ctrl$gene_snp))
  pme.m.sig = pme.m.sig %>% mutate("is_liver_eQTL" = (gene_snp %in% pme.m.sig.null$gene_snp),
                                   'is_nafld_eQTL' = (gene_snp %in% pme.m.sig.NAFLD$gene_snp),
                                   'is_ctrl_eQTL'=(gene_snp %in% pme.m.sig.ctrl$gene_snp))


  write.csv(pme.m.sig, paste0(analysis_results_dir, "pme.sig.0.05.", celltype_cur, '.csv'), row.names = F)
  fwrite(pme.m, paste0(analysis_results_dir, "pme.all.merged.", celltype_cur, '.csv.gz'), compress='gzip')
  saveRDS( list(max_p_005_null = max_p_0.05_thr_null, max_p_005_NAFLD = max_p_0.05_thr_NAFLD, max_p_005_ctrl = max_p_0.05_thr_ctrl,
                bestP_null = pme.m_bestP_null, bestP_NAFLD = pme.m_bestP_NAFLD, bestP_ctrl = pme.m_bestP_ctrl) , paste0(analysis_results_dir, "p_values_stats_", celltype_cur, ".rds"))
  saveRDS(tested.genes, paste0(analysis_results_dir, "eQTL_tested_genes_", celltype_cur, ".rds"))
}
