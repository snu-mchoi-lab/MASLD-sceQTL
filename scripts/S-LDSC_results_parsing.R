japan.data.dir = "path_to_raw_BBJ_gwas_summary_stats"
japan.dataname_list = list.files(paste0(japan.data.dir, "raw_summary_stats/"), pattern=".txt.gz$")
japan.dataname_list = lapply(japan.dataname_list, function(x){str_replace(x, ".autosome.txt.gz", "") %>% str_replace(., "BBJ.", "")}) %>% unlist()

ldsc_categories=c("all_eqtl", "liver_eqtl", "ieqtl", lapply(celltypes, function(x){paste0(x, c("_i", "_liver"))}) %>% unlist(), "Hep_M3_4_6","Hep_M3_4_6_12",  "Hep_M8_14")

for(gwas_cur in japan.dataname_list){
  gwas_short=paste0("BBJ.", gwas_cur)

  print(gwas_short)
  for(lc in ldsc_categories){
    res_tmp <- read.table(paste0(ldsc_output_dir, 'output/BBJ.', gwas_cur, ".sumstats__", lc, ".results"), header=1) %>%
      mutate(ldsc_category=lc,
             gwas_phenotype=gwas_short)
    if(lc == ldsc_categories[1] & gwas_cur == japan.dataname_list[[1]]){ ldsc_res_merged = res_tmp}
    else{ ldsc_res_merged = rbind(ldsc_res_merged, res_tmp)}
  }
}

ldsc_res_merged.f = ldsc_res_merged %>% filter(Category == 'ANNOTL2_0')
ldsc_res_merged.f$ldsc_category = factor(ldsc_res_merged.f$ldsc_category, levels=ldsc_categories)

ldsc_res_merged.f %>%
  ggplot(., aes(y=ldsc_category, x=Enrichment, alpha=(Enrichment_p<0.05))) + geom_point(aes(size=-log10(Enrichment_p))) +
  geom_errorbar(aes(xmin=Enrichment-Enrichment_std_error*qnorm(0.975),
                    xmax=Enrichment+Enrichment_std_error*qnorm(0.975)), width=0.4) +
  facet_wrap(~gwas_phenotype, scales='free', nrow=4)  + # nrow=1, width=60, height=3.2
  scale_size(range=c(0.5,2))+
  geom_vline(xintercept=1, linetype='dashed') + theme_bw() + scale_alpha_discrete(range=c(0.3, 1)) #+
