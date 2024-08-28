# preprocess for LDSC

# eQTLs
make_eqtl_bed <- function(pme.sig.anno.f, prefix_cur){
  write.table(pme.sig.anno.f, paste0(ldsc_dir, "input/", prefix_cur, ".bed"), sep='\t', col.names=F, row.names=F, quote=F)
  for(chr_cur in 1:22){ write.table(pme.sig.anno.f %>% filter(chr==paste0("chr", chr_cur)),
                                    paste0(ldsc_dir_v2,"input/", prefix_cur,".", chr_cur, ".bed"),
                                    sep='\t', col.names=F, row.names=F, quote=F) }
}

(pme.sig.anno)[,c("snp", "CHROM", "POS", "REF", "ALT")] %>%
  distinct() %>%
  mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>%
  dplyr::select(chr, start, end) %>% arrange(chr, start) %>%
  make_eqtl_bed(., "all_eqtl")


(pme.sig.anno[which(rowSums(pme.sig.anno[,paste0("is_liver_eQTL_", celltypes)]) > 0), ])[,c("snp", "CHROM", "POS", "REF", "ALT")] %>%
  distinct() %>%
  mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>%
  dplyr::select(chr, start, end) %>% arrange(chr, start) %>%
  make_eqtl_bed(., "liver_eqtl")


(pme.sig.anno[which(rowSums(pme.sig.anno[,paste0("is_ieQTL_", celltypes)]) > 0), ])[,c("snp", "CHROM", "POS", "REF", "ALT")] %>%
  distinct() %>%
  mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>%
  dplyr::select(chr, start, end) %>% arrange(chr, start) %>%
  make_eqtl_bed(., "ieqtl")

for(celltype_cur in celltypes){
  (pme.sig.anno[which(rowSums(pme.sig.anno[,paste0(c("is_liver_eQTL_", "is_ieQTL_"), celltype_cur)]) > 0), ])[,c("snp", "CHROM", "POS", "REF", "ALT")] %>%
    distinct() %>%
    mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>%
    dplyr::select(chr, start, end) %>% arrange(chr, start) %>%
    make_eqtl_bed(.,paste0(celltype_cur,"_all"))
}

for(celltype_cur in celltypes){
  (pme.sig.anno[which(pme.sig.anno[,paste0("is_liver_eQTL_", celltype_cur)]), ])[,c("snp", "CHROM", "POS", "REF", "ALT")] %>%
    distinct() %>%
    mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>%
    dplyr::select(chr, start, end) %>% arrange(chr, start) %>%
    make_eqtl_bed(.,paste0(celltype_cur,"_liver"))
}

for(celltype_cur in celltypes){
  (pme.sig.anno[which(pme.sig.anno[,paste0("is_ieQTL_", celltype_cur)]), ])[,c("snp", "CHROM", "POS", "REF", "ALT")] %>%
    distinct() %>%
    mutate(chr=paste0('chr',CHROM), start= POS-1,end= POS) %>%
    dplyr::select(chr, start, end) %>% arrange(chr, start) %>%
    make_eqtl_bed(.,paste0(celltype_cur,"_i"))
}

# gwas summary stats
japan.data.dir = "path_to_raw_BBJ_gwas_summary_stats"
japan.dataname_list = list.files(paste0(japan.data.dir, "raw_summary_stats/"), pattern=".txt.gz$")
japan.dataname_list = lapply(japan.dataname_list, function(x){str_replace(x, ".autosome.txt.gz", "") %>% str_replace(., "BBJ.", "")}) %>% unlist()
japan.dataname_list
harmonized.dir="path_to_harmonized_BBJ_gwas_summary_stats"

output_dir= "path_to_output_directory"
for(data.name.cur in japan.dataname_list){
  gwas.data.m = fread(paste0(harmonized.dir, "BBJ.", data.name.cur,"_GRCh38_harmonized.txt"))
  gwas.data.m.tmp = gwas.data.m[,c("ID_hg38", "POS", "CHROM", "N", "BETA", "SE", "ALT", "REF")] %>%
    mutate("Z" = BETA/SE, SNP=paste0(str_replace(CHROM, "chr",""), "_", POS)) %>%
    dplyr::select("SNP", "N", "Z", 'ALT', "REF") %>%
    dplyr::rename(A1=ALT, A2=REF)
  write.table(gwas.data.m.tmp, paste0(output_dir, "BBJ.", data.name.cur, ".sumstats"), row.names=F, quote=F)
}
