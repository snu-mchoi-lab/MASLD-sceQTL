#####################################################
# S-LDSC permutation analysis: make input bed files
###################################################
# random sample SNPs from GT matrix x 100 times 
n_all = (pme.sig.anno)[,c("snp", "CHROM", "POS", "REF", "ALT")] %>%  distinct() %>% nrow()
n_liver = (pme.sig.anno[which(rowSums(pme.sig.anno[,paste0("is_liver_eQTL_", celltypes)]) > 0), ])[,c("snp", "CHROM", "POS", "REF", "ALT")] %>% distinct() %>% nrow()
n_i = (pme.sig.anno[which(rowSums(pme.sig.anno[,paste0("is_ieQTL_", celltypes)]) > 0), ])[,c("snp", "CHROM", "POS", "REF", "ALT")] %>% distinct() %>% nrow()
# liver-eQTLs
# ieQTLs

ldsc_dir_v4 = '...'

all.snps.data <- read.table("all.snps.tested.chrmState.nafld.atac.chipseq.gwasCatalog.parsed.txt", sep='\t') %>% select(ID)
all.anno.vcf <- read.vcfR("all.snps.tested.chrmState.nafld.atac.chipseq.gwasCatalog.vcf")@fix[,c('CHROM', "POS", "ID", "REF", "ALT")] %>% data.frame

all_snpset_list = vector('list', length=100)
liver_snpset_list = vector('list', length=100)
i_snpset_list = vector('list', length=100)

for(i in 1:100){
  all_snpset_list[[i]] = (all.snps.data %>% sample_n(n_all))$ID
  liver_snpset_list[[i]] = (all.snps.data %>% sample_n(n_liver))$ID
  i_snpset_list[[i]] = (all.snps.data %>% sample_n(n_i))$ID
}

saveRDS(list("all_snpsets" = all_snpset_list,
             'liver_snpsets' = liver_snpset_list,
             'i_snpsets' = i_snpset_list),
        paste0(ldsc_dir_v4, "input/snpsets_list.rds"))


make_eqtl_bed <- function(pme.sig.anno.f, prefix_cur){
  write.table(pme.sig.anno.f, paste0(prefix_cur, ".bed"), sep='\t', col.names=F, row.names=F, quote=F)
  for(chr_cur in 1:22){ write.table(pme.sig.anno.f %>% filter(chr==paste0("chr", chr_cur)),
                                    paste0(prefix_cur,".", chr_cur, ".bed"),
                                    sep='\t', col.names=F, row.names=F, quote=F) }
}

all.anno.vcf$POS = as.numeric(all.anno.vcf$POS)
for(i in 1:100){
  all.anno.vcf[,c("ID",  "CHROM", "POS", "REF", "ALT")] %>% 
    filter(ID %in% all_snpset_list[[i]]) %>% 
    mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>% 
    dplyr::select(chr, start, end) %>% arrange(chr, start) %>% 
    make_eqtl_bed(., paste0(ldsc_dir_v4, "input/all_snpset/all_snpset", i))
  
  all.anno.vcf[,c("ID",  "CHROM", "POS", "REF", "ALT")] %>% 
    filter(ID %in% liver_snpset_list[[i]]) %>% 
    mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>% 
    dplyr::select(chr, start, end) %>% arrange(chr, start) %>% 
    make_eqtl_bed(., paste0(ldsc_dir_v4, "input/liver_snpset/liver_snpset", i))
  
  all.anno.vcf[,c("ID",  "CHROM", "POS", "REF", "ALT")] %>% 
    filter(ID %in% i_snpset_list[[i]]) %>% 
    mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>% 
    dplyr::select(chr, start, end) %>% arrange(chr, start) %>% 
    make_eqtl_bed(., paste0(ldsc_dir_v4, "input/i_snpset/i_snpset", i))
}

##############################################
# merge & plot results 
##############################################
annot_list = c("all_snpset", 'liver_snpset', 'i_snpset')
gwas_list = c("BBJ.ALT", "BBJ.GGT", "Nat22_cALT_NAFLD")

ldsc_v4_output_dir="./eQTL_calling/211005_n48_A02/results/v15_DonorWise/PME/v3/analysis/S_LDSC/v4_permutation/output/ldsc_results/"
for(annot in annot_list){
  for(set_i in 1:100){
    for(gwas_cur in gwas_list){
      ldsc_res_tmp = read.table(paste0(ldsc_v4_output_dir, gwas_cur, ".sumstats_", annot, set_i, ".results"), header=1) %>% filter(Category == "ANNOTL2_0") %>%
        mutate('eQTL_category' = gsub("_snpset", "", annot), subset_i = set_i, gwas_phenotype=gwas_cur)
      if(annot == annot_list[1] & set_i==1 & gwas_cur == gwas_list[1]){ ldsc_res = ldsc_res_tmp}
      else{ ldsc_res = rbind(ldsc_res, ldsc_res_tmp)}
    }
  }
}


# for all
original_data = read.csv("./eQTL_calling/211005_n48_A02/results/v15_DonorWise/PME/v3/analysis/s_ldsc_enrichment_results.csv")
original_data = original_data %>% filter(gwas_phenotype %in% gwas_list) %>% filter(ldsc_category %in% paste0(c("all_", 'liver_', 'i'), 'eqtl'))
ldsc_eqtl_conv_df = data.frame(eQTL_category= gsub("_snpset", "", annot_list),
                               ldsc_category = c("all_eqtl", 'liver_eqtl', 'ieqtl'))
original_data = merge(original_data, ldsc_eqtl_conv_df, by='ldsc_category', all.x=T)

for(annot in annot_list){
  for(gwas_cur in gwas_list){
    annot_short = gsub("_snpset", '', annot)
    original_Enrichment_val = (original_data %>% filter(eQTL_category == annot_short) %>% filter(gwas_phenotype == gwas_cur))$Enrichment
    original_Enrichment_p = (original_data %>% filter(eQTL_category == annot_short) %>% filter(gwas_phenotype == gwas_cur))$Enrichment_p
    
    ldsc_tmp = ldsc_res %>% filter(eQTL_category == annot_short) %>% filter(gwas_phenotype == gwas_cur)
    enrich_val_higher = mean(ldsc_tmp$Enrichment > original_Enrichment_val)
    enrich_p_lower_fraction = mean(ldsc_tmp$Enrichment_p < original_Enrichment_p)
    
    res_tmp = data.frame("eQTL_category" = annot_short, gwas_phenotype=gwas_cur, Enrichment_higher=enrich_val_higher, Enrichment_P_lower = enrich_p_lower_fraction)
    
    if(annot == annot_list[1] & gwas_cur == gwas_list[1]) { res_merged = res_tmp }
    else{ res_merged = rbind(res_merged, res_tmp)}
  }
}

ldsc_res_for_plot = ldsc_res[,c("gwas_phenotype", 'eQTL_category', 'Enrichment', 'Enrichment_std_error', "Enrichment_p")] %>% mutate('origin' = 'permutation')
ldsc_res_for_plot = rbind(ldsc_res_for_plot, 
                          original_data[,c("gwas_phenotype", 'eQTL_category', 'Enrichment', 'Enrichment_std_error', "Enrichment_p")] %>% mutate('origin' = 'Real data'))

ggplot(ldsc_res_for_plot %>% filter(origin=='permutation'), aes(x=gwas_phenotype, y=Enrichment)) + 
  geom_jitter( aes(size = (-log10(Enrichment_p)), color=eQTL_category), position=position_jitterdodge(jitter.width=0.1, dodge.width=0.8), alpha=0.1) + 
  geom_boxplot(aes(fill=eQTL_category), alpha=0.8, position=position_dodge(width=0.8), width=0.5, outlier.shape=NA) +
  geom_point(data = ldsc_res_for_plot %>% filter(origin=='Real data'), aes(size = (-log10(Enrichment_p)), color=eQTL_category), position=position_dodge(width=0.8), shape=15)+
  theme_classic() + geom_hline(yintercept=1, linetype='dashed') + 
  scale_color_manual(values=c("black", 'slategrey', 'wheat4'), labels = c('all' = "all eQTLs", 'i' = "ieQTLs", "liver" = 'liver-eQTLs')) + scale_fill_manual(values=c("black", 'slategrey', 'wheat4'), labels = c('all' = "all eQTLs", 'i' = "ieQTLs", "liver" = 'liver-eQTLs'))  +
  xlab("GWAS phenotype") + scale_x_discrete(labels = c('BBJ.ALT' = "ALT", "BBJ.GGT" = "GGT", "Nat22_cALT_NAFLD"= "MASLD (cALT)")) 
ggsave("LDSC_permutation_Enrichment_results.pdf", width=6, height=3)

