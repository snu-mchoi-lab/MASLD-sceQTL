# meme suite analysis

half.window=15
for(snp_i in 1:nrow(pme.sig.motif)){
  snp_id = pme.sig.motif[snp_i, 'snp']
  snp_pos_cur = pme.sig.motif[snp_i,"POS"]
  snp_chrom_cur = pme.sig.motif[snp_i, "CHROM"]
  ref_allele = pme.sig.anno[snp_i,'REF']
  alt_allele = pme.sig.anno[snp_i,'ALT']
  ref_length = str_length(ref_allele)

  left_arm = getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0("chr", snp_chrom_cur), start=(snp_pos_cur - half.window), end=(snp_pos_cur-1))
  right_arm = getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0("chr", snp_chrom_cur), start=(snp_pos_cur + ref_length), end=(snp_pos_cur+ ref_length +half.window))

  ref_region = paste0(as.character(left_arm), ref_allele, as.character(right_arm))
  alt_region = paste0(as.character(left_arm), alt_allele, as.character(right_arm))
  meme_tmp = c(paste0(">", snp_id, "_ref"), ref_region,
               paste0(">", snp_id, "_alt"), alt_region)

  if(snp_i==1){ meme_list <- meme_tmp}
  else{meme_list <- c(meme_list, meme_tmp)}
}

fileConn <- file(paste0(analysis_results_dir, "MEME_suite/input_files/pme.sig.anno.fa"))
writeLines(meme_list, fileConn)
close(fileConn)


fimo_path = "./meme-5.5.1/src/fimo"
db_dir= "./motif_databases/HUMAN/"
output_path = paste0(analysis_results_dir, "MEME_suite/output_files")
system(paste0(fimo_path, " --oc ", output_path, " --verbosity 1 --thresh 1.0e-4 max-stored-scores 4000000 ",
              db_dir, "HOCOMOCOv11_core_HUMAN_mono_meme_format.meme ", paste0(analysis_results_dir, "MEME_suite/input_files/pme.sig.anno.fa")))
