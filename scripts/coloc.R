library(coloc)
analysis_results_dir="path_to_results_dir"
pme.sig.anno <-fread(paste0(analysis_results_dir, "pme.sig.anno.full.txt"), sep='\t') %>% data.frame(., check.names = F)


# use gw_merged data
calculate_coloc_pme <- function(eqtl_res_df, gwas_summary_df, category){
  if(category == "full"){
    eqtl_input = list ( "beta" = eqtl_res_df$Estimate_full,
                        "varbeta" = (eqtl_res_df$Std_Error_full)^2,
                        "snp" = eqtl_res_df$snp,
                        "position" = eqtl_res_df$position,
                        "type" = "quant",
                        "sdY" = 1)}
  if(category == "null"){
    eqtl_input = list ( "beta" = eqtl_res_df$Estimate_null,
                        "varbeta" = (eqtl_res_df$Std_Error_null)^2,
                        "snp" = eqtl_res_df$snp,
                        "position" = eqtl_res_df$position,
                        "type" = "quant",
                        "sdY" = 1)}
  if(category == "NAFLD"){
    eqtl_input = list ( "beta" = eqtl_res_df$Estimate_NAFLD,
                        "varbeta" = (eqtl_res_df$Std_Error_NAFLD)^2,
                        "snp" = eqtl_res_df$snp,
                        "position" = eqtl_res_df$position,
                        "type" = "quant",
                        "sdY" = 1)}
  if(category == "ctrl"){
    eqtl_input = list ( "beta" = eqtl_res_df$Estimate_ctrl,
                        "varbeta" = (eqtl_res_df$Std_Error_ctrl)^2,
                        "snp" = eqtl_res_df$snp,
                        "position" = eqtl_res_df$position,
                        "type" = "quant",
                        "sdY" = 1)}

  # set 2: GWAS
  gwas.input = list("beta"= gwas_summary_df$beta_h,
                    "varbeta" = (gwas_summary_df$SE)^2,
                    "snp" = gwas_summary_df$ID ,
                    "position"= gwas_summary_df$POS,
                    "type" = "cc")

  tryCatch({
    res <- coloc.abf(dataset1=eqtl_input,
                     dataset2=gwas.input)
    top.snp.info = (res$results %>% arrange(-SNP.PP.H4))[1,]
    res_df <- res$summary %>% t() %>% data.frame() %>% mutate("coloc_top_snp" = top.snp.info[1,'snp'],
                                                              'coloc_top_SNP.PP.H4' = top.snp.info[1,'SNP.PP.H4'])
  },error = function(cond){return(NA)})

  #return(res_df)
  return(res_df)
}

sig.egenes <- pme.merged.anno$gene %>% unique()
dataname_list=c("Nat22_cALT_NAFLD_GWAS", "GCST90019492_UKBB_ALT", "GCST90019496_UKBB_ApoB", "GCST90019494_UKBB_ALP", "GCST90019501_UKBB_TotalCholesterol", "GCST90019523_UKBB_Triglyceride", "GCST90019495_UKBB_ApoA1", "GCST90019497_UKBB_AST", "GCST90019507_UKBB_GGT", "GCST90019512_UKBB_LDL", "GCST90019521_UKBB_TotalBilirubin", "GCST90019499_UKBB_CRP")
dataname_list2=c("GCST90019500_UKBB_Ca", "GCST90019519_UKBB_UrineSodium", "GCST90019506_UKBB_eGFR", "GCST90019509_UKBB_HbA1C")
dataname_list3 = c("GCST90019517_UKBB_UrineK", "GCST90019524_UKBB_Urate", "GCST90019504_UKBB_CystatinC", "GCST90019510_UKBB_HDL", "GCST90019511_UKBB_IGF1", "GCST90019516_UKBB_Phosphate", "GCST90019526_UKBB_VitD")
dataname_list_m= c(dataname_list, dataname_list2, dataname_list3)
data.dir="path_to_gwas_sumstats"
coloc_results_dir="path_to_output"


celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')


for(celltype_cur in celltypes){
  gw_merged = fread(paste0(analysis_results_dir, "pme.all.merged.", celltype_cur, ".csv.gz"))
  gw_merged$position = lapply(gw_merged$snp, function(x){str_split(x, ":")[[1]][2]}) %>% unlist() %>% as.numeric()
  sig.egenes = (pme.sig.anno %>% filter(!!sym(paste0("is_liver_eQTL_", celltype_cur)) | !!sym(paste0("is_ieQTL_", celltype_cur))))$gene %>% unique

  for(m in 1:length(dataname_list_m)){
    data.name=dataname_list_m[m]
    gwas.data <- fread(paste0(data.dir, data.name, "_GRCh38_dataMerged.txt"))
    gwas.data = gwas.data %>% filter(!duplicated(ID))

    k=0
    for (egene in sig.egenes){
      print(paste0("Running coloc on.....", egene))
      k = k+1
      common.snps = intersect( (gw_merged %>% filter(gene == egene))$snp,
                               gwas.data$ID)
      eqtl_res_df_g = gw_merged %>% filter(gene ==egene) %>% filter(snp %in% common.snps) # CHEK2
      gwas_summary_df_g = gwas.data %>% filter(ID %in% common.snps)

      coloc_tmp <- data.frame("gene" = egene,
                              "nsnps" = NA,
                              "PP.H0.abf_null"=NA,
                              "PP.H1.abf_null"=NA,
                              "PP.H2.abf_null"=NA,
                              "PP.H3.abf_null"=NA,
                              "PP.H4.abf_null"=NA,
                              "coloc_top_snp_null"=NA,
                              "coloc_top_SNP.PP.H4_null"=NA,
                              "PP.H0.abf_ctrl"=NA,
                              "PP.H1.abf_ctrl"=NA,
                              "PP.H2.abf_ctrl"=NA,
                              "PP.H3.abf_ctrl"=NA,
                              "PP.H4.abf_ctrl"=NA,
                              "coloc_top_snp_ctrl"=NA,
                              "coloc_top_SNP.PP.H4_ctrl"=NA,
                              "PP.H0.abf_NAFLD"=NA,
                              "PP.H1.abf_NAFLD"=NA,
                              "PP.H2.abf_NAFLD"=NA,
                              "PP.H3.abf_NAFLD"=NA,
                              "PP.H4.abf_NAFLD"=NA,
                              "coloc_top_snp_NAFLD"=NA,
                              "coloc_top_SNP.PP.H4_NAFLD"=NA)

      if( min(gwas_summary_df_g$p_value) < 1e-6 ) {
        coloc_tmp[1,2:9] <- calculate_coloc_pme(eqtl_res_df_g, gwas_summary_df_g, category='null')
        coloc_tmp[1,10:16] <- calculate_coloc_pme(eqtl_res_df_g, gwas_summary_df_g, category='ctrl')[,2:8]
        coloc_tmp[1,17:23] <- calculate_coloc_pme(eqtl_res_df_g, gwas_summary_df_g, category='NAFLD')[,2:8]
      }

      if(k==1){
        coloc_res <- coloc_tmp
      }
      else{
        coloc_res = rbind(coloc_res, coloc_tmp)
      }
    }

    write.table(coloc_res, paste0(coloc_results_dir, data.name, "_coloc_res.",celltype_cur, ".txt"), sep='\t', row.names = F)
    gc()
  }
}


# filter, process coloc results
coloc_res_m_list = vector('list', 4)
celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')
names(coloc_res_m_list) <- celltypes

# should have calculated coloc on "GT_gd not null -> ref) liver eQTL calling method)
check_concordance = function(gene_cur, is_something_columns_celltype){
  result = (pme.sig.anno[,c('gene','snp',is_something_columns_celltype)] %>% filter(gene==gene_cur))[,paste0(c("is_nafld_eQTL", "is_ctrl_eQTL", "is_liver_eQTL"), "_", celltype_cur)] %>% colSums()>0
  return(result)
}

for(celltype_cur in celltypes){
  for(m in 1:length(dataname_list_m)){
    data.name=dataname_list_m[m]

    is_something_columns = colnames(pme.sig.anno)[str_starts(colnames(pme.sig.anno), "is_")]
    is_something_columns_celltype = is_something_columns[str_ends(is_something_columns, paste0("_", celltype_cur))]

    coloc_res = read.table(paste0(coloc_results_dir, data.name, "_coloc_res.", celltype_cur, ".txt"), header=1) %>% filter(!is.na(nsnps)) %>%
      dplyr::select(gene, nsnps,PP.H4.abf_null, PP.H4.abf_NAFLD, PP.H4.abf_ctrl) %>%  # PP.H4.abf_full, PP.H4.abf_null,
      filter( PP.H4.abf_null > 0.8 | PP.H4.abf_NAFLD > 0.8 | PP.H4.abf_ctrl > 0.8) #PP.H4.abf_null > 0.8 | PP.H4.abf_full > 0.8 |

    if(nrow(coloc_res) == 0 ){next}
    coloc_res = coloc_res %>% mutate(!!paste0("is_NAFLD_eGene_", celltype_cur) := NA,
                                     !!paste0("is_ctrl_eGene_", celltype_cur) := NA,
                                     !!paste0("is_liver_eGene_", celltype_cur) := NA)
    coloc_res[,paste0(c("is_NAFLD_eGene", "is_ctrl_eGene", "is_liver_eGene"), "_", celltype_cur)] <- lapply(coloc_res$gene, function(x){check_concordance(x, is_something_columns_celltype)}) %>% data.frame() %>% t()
    coloc_res = coloc_res %>% mutate(!!paste0("NAFLD_concordant_coloc_",data.name) := (PP.H4.abf_NAFLD > 0.8) & !!sym(paste0("is_NAFLD_eGene_", celltype_cur)),
                                     !!paste0("ctrl_concordant_coloc_", data.name) := (PP.H4.abf_ctrl > 0.8) & !!sym(paste0("is_ctrl_eGene_", celltype_cur)),
                                     !!paste0("liver_concordant_coloc_", data.name) := (PP.H4.abf_null > 0.8) & !!sym(paste0("is_liver_eGene_", celltype_cur)))
    coloc_res = coloc_res[,c("gene",c(paste0(c('is_NAFLD_eGene', 'is_ctrl_eGene', "is_liver_eGene"), "_", celltype_cur),
                                      paste0(c("NAFLD", "ctrl", "liver"), "_concordant_coloc", "_",data.name)))]
    if(m==1) { coloc_res_m = coloc_res}
    else{ coloc_res_m = merge(coloc_res_m, coloc_res, by=c("gene",paste0(c('is_NAFLD_eGene_', 'is_ctrl_eGene_', "is_liver_eGene_"), celltype_cur)), all=T) }
    }
 coloc_res_m_list[[celltype_cur]] <- coloc_res_m
 }


saveRDS(coloc_res_m_list, paste0(analysis_results_dir, "coloc_res_m.rds"))

coloc_counts = sapply(names(coloc_res_m_list),
                      function(x){(coloc_res_m_list[[x]] %>% column_to_rownames("gene"))[,4:(ncol(coloc_res_m_list[[x]])-1)] %>%
                          as.data.table()%>%
                          colSums(., na.rm=T) %>% data.frame('egene_n' = . ) %>%
                          dplyr::rename(.,!!sym(x):=egene_n) %>% mutate("gwas" = rownames(.))},
                      simplify = FALSE,USE.NAMES = TRUE)
coloc_counts_merged = Reduce(function(x, y) merge(x, y, all=TRUE, by='gwas'), coloc_counts)

coloc_counts_merged$gwas= lapply(coloc_counts_merged$gwas, function(x){str_replace(x, "concordant_coloc_", "")}) %>% unlist
coloc_counts_merged$gwas= lapply(coloc_counts_merged$gwas, function(x){str_replace(x, "_GCST[0-9]+", "")}) %>% unlist

gwas.levels=c("Nat22_cALT_NAFLD_GWAS", "UKBB_ALT", "UKBB_AST", "UKBB_GGT", "UKBB_ALP", "UKBB_TotalBilirubin",
              "UKBB_CRP", "UKBB_HbA1C", "UKBB_Triglyceride", "UKBB_TotalCholesterol","UKBB_LDL", "UKBB_HDL", "UKBB_ApoA1", "UKBB_ApoB",
              "UKBB_eGFR","UKBB_CystatinC", "UKBB_IGF1", "UKBB_Ca", "UKBB_Phosphate", "UKBB_Urate","UKBB_VitD",  "UKBB_UrineSodium", "UKBB_UrineK")
coloc_counts_merged$gwas = factor(coloc_counts_merged$gwas,
                                  levels= paste0(rep( c("liver_", "NAFLD_", "ctrl_"), length(gwas.levels)), rep(gwas.levels, each=3)))
coloc_counts_merged = coloc_counts_merged %>% arrange(gwas)
coloc_counts_merged = coloc_counts_merged %>% column_to_rownames("gwas")


pdf(paste0(analysis_results_dir, 'coloc_eGene_count_heatmap.pdf'), width=5.5, height=3)
pheatmap(coloc_counts_merged[str_starts(rownames(coloc_counts_merged), "liver_"),] %>% t(), cluster_cols=F, cluster_rows=F)
dev.off()
