# integrative analysis
# parse meme-suite fimo results
fimores <- read.csv(paste0(output_path, "/pme.sig.eqtls/fimo.tsv"), sep='\t')
fimores = fimores %>% filter(start < 15) %>% filter(stop>15)
fimores = fimores %>% unique()
fimores = fimores %>% tidyr::separate(., col='sequence_name', sep="_", remove=F, into=c("snp", 'allele'))

# differential TF binding
## remove if same TF present in
fimores_ref = fimores %>% filter(allele == 'ref') %>% dplyr::select(-c("sequence_name", 'allele', 'motif_alt_id'))
fimores_alt= fimores %>% filter(allele=='alt') %>% dplyr::select(-c("sequence_name", 'allele', 'motif_alt_id'))
fimores_m = merge(fimores_ref, fimores_alt, by=c("motif_id", 'snp', 'start', 'strand'), suffixes=c("_ref", "_alt"), all=T)
fimores_m = fimores_m %>% mutate("ref_bind" = (!is.na(`q.value_ref`)) & (`p.value_ref` < 1e-4),
                                 "alt_bind" = (!is.na(`q.value_alt`)) & (`p.value_alt` < 1e-4),
                                 "differential_bind" = ((ref_bind + alt_bind) == 1) )
fimores_d = fimores_m %>% filter(differential_bind)
fimores_d$motif_id %>% table %>% sort()

# merge TF motif & module & ieQTL
wgcna.module.name=list("hepatocyte" = 'Hep-M', 'cholangiocyte'= "Chol-M", "stellate_cell" = "HSC-M", 'endothelial_cell' = "endo-M")
wgcna.module.count = list('hepatocyte'=c(3,4,6,8,12,14), 'cholangiocyte'=c(1,2,3,4,6,7), 'stellate_cell' = 1:3, 'endothelial_cell' = c(1,4,7,9))
wgcna.filename.prefix= list('hepatocyte' = "hep", 'cholangiocyte' = 'chol', 'stellate_cell'='hsc','endothelial_cell'='endo' )
ieQTL_results_path = 'path_to_ieQTL_summary_stats'

for(celltype_cur in celltypes){
  modules_c = paste0(wgcna.module.name[[celltype_cur]], wgcna.module.count[[celltype_cur]])
  corrs= read.csv(paste0(scenic_output_dir, celltype_cur, "_module_scenic_corr.csv"))
  colnames(corrs) <- c("X", "TF", "WGCNA_module", "TF_module_cor_estimate", "cor_p_value", "cor_statistic", "cor_method", "cor_p_adj_bf")

  for(module_cur in modules_c){
    ieqtl_stats=fread(paste0(ieQTL_results_path, "interaction.", celltype_cur, "_", module_cur, ".csv"))
    ieqtl_stats = ieqtl_stats %>% filter(term=="G") %>% as.data.frame()

    corrs_tmp = corrs %>% filter(cor_p_adj_bf < 0.05) %>% filter(TF_module_cor_estimate > 0) %>% filter(WGCNA_module == module_cur)
    fimores_tmp = merge(fimores_d, (pme.sig.anno %>% filter(!is.na(!!sym(module_cur))))[,c("snp", 'gene', module_cur)], by='snp', all.y=T)
    fimores_tmp$TF_short = lapply(fimores_tmp$motif_id, function(x){str_split(x, "_HUMAN")[[1]][1]}) %>% unlist

    fimores_tmp = merge(fimores_tmp, ieqtl_stats[,c('gene', 'snp', 'Estimate_full', "pval_full")],
                        by=c("gene",'snp'), all.x=T)
    motif_scenic = merge(corrs_tmp, fimores_tmp, by.x="TF", by.y="TF_short")
    motif_scenic = motif_scenic %>% arrange(-TF_module_cor_estimate)
    colnames(motif_scenic)[which(colnames(motif_scenic) == module_cur)] <- "Estimate_interaction"

    if(module_cur == modules_c[1]){ motif_scenic_merged = motif_scenic}
    else{motif_scenic_merged = rbind(motif_scenic_merged, motif_scenic)}
  }
  if(celltype_cur==celltypes[1]){
    motif_scenic_merged.m = motif_scenic_merged %>% mutate(celltype= celltype_cur)
  }
  else{
    motif_scenic_merged.m = rbind(motif_scenic_merged.m, motif_scenic_merged %>% mutate(celltype= celltype_cur))
  }

  write.csv(motif_scenic_merged, paste0(save.dir, "scenic_wgcna_meme_merged_", celltype_cur, ".csv"))
}

motif_scenic_merged.m = motif_scenic_merged.m %>% filter(TF_module_cor_estimate > 0.2)



# compare with GWAS snps
nafld.gwas ="path_to_MASLD_GWAS_summary_stats"
gwas_sumstats = fread(nafld.gwas) %>% filter(P < 5e-8)
gwas_sumstats = gwas_sumstats %>% select(ID, P, beta_h, REF_hg38, ALT_hg38) %>% mutate("cALT_MASLD"=T)
colnames(gwas_sumstats)[1:5] = c("snp", paste0(c("P", "beta_h"), "_cALT_MASLD"), "ref", "alt")
gwas_sumstats = gwas_sumstats %>% filter(!is.na(ref))

japan.data.dir = "path_to_BBJ_GWAS_summary_stats"
for( phenotype in c("BBJ.ALT", "BBJ.AST", "BBJ.GGT")){
  gwas_sumstats_tmp = fread(paste0(japan.data.dir, phenotype, "_GRCh38_harmonized.txt"))%>% filter(P < 5e-8) %>%
    select("ID_hg38", "P", "BETA", "REF", "ALT") %>% mutate(!!phenotype := T)
  colnames(gwas_sumstats_tmp)[1:5] =c("snp", paste0(c("P", "beta_h"), "_", phenotype), "ref", "alt")
  gwas_sumstats = merge(gwas_sumstats, gwas_sumstats_tmp, by=c("snp",'ref','alt'),all=T)
  print(phenotype)
}

ukbb.data.dir="path_to_UKBB_GWAS_summary_stats"
for (phenotype_long in c("GCST90019492_UKBB_ALT", "GCST90019497_UKBB_AST", "GCST90019507_UKBB_GGT")){
  phenotype=gsub(pattern="GCST[0-9]+_", "", phenotype_long)
  gwas_sumstats_tmp = fread(paste0(ukbb.data.dir, phenotype_long, "_GRCh38_harmonized.txt"))%>% filter(P < 5e-8) %>%
    select("ID", "P", "BETA", "REF", "ALT") %>% mutate(!!phenotype := T)
  colnames(gwas_sumstats_tmp)[1:5] =c("snp", paste0(c("P", "beta_h"), "_", phenotype), "ref", "alt")
  gwas_sumstats = merge(gwas_sumstats, gwas_sumstats_tmp, by=c("snp",'ref','alt'),all=T)
  print(phenotype)
}

write.table(gwas_sumstats, paste0(analysis_results_dir, "cALT_BBJ_UKBB_ALT_AST_GGT_gwas_5e-8_snps_union.txt"))

motif_scenic_merged.m %>% filter(snp %in% unique(gwas_sumstats$snp)) %>% select(TF, WGCNA_module, TF_module_cor_estimate, gene, snp, celltype) %>% unique() %>% dim() # 23
write.csv(motif_scenic_merged.m %>% filter(snp %in% unique(gwas_sumstats$snp)) , paste0(analysis_results_dir, "quartets_n23_gwas.csv"))
