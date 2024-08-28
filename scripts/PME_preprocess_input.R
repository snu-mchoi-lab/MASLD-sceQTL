# PME model

# preprocessing
excluded.samples = c("FL551", "FL568", "FL406", "FL597", "FL296", "FL605",
                    "FL284", "FL511", "FL531", "FL793") # too small cell count

included.samples = celltype_ser_list[[1]]@meta.data$orig.ident %>% unique()
included.samples = included.samples[which(!included.samples %in% excluded.samples)]

clininfo <- read.table(paste0(save.dir, "clininfo.txt"), sep='\t')
clininfo <- merge(clininfo,
                  data.frame("disease_group2" = c("ctrl", "NAFL", "eNASH", "aNASH"),
                             "disease_group3" = c("ctrl", "NAFLD", "NAFLD", "NAFLD"),
                             "disease_group_num" = c(0,1,1,1)), by.x="disease_group", by.y='disease_group2')

# write ge matrix, covariates matrix
for(i in 1:4){
  # remove low QC samples
  ser.tmp <- celltype_ser_list[[i]]
  ser.tmp = subset(ser.tmp, subset=(orig.ident %in% included.samples))
  if(i %in% c(2,4)){ ser.tmp@meta.data$cell_type_4 = ser.tmp@meta.data$cell_type_1 }
  metadata.tmp <- ser.tmp@meta.data[,c("orig.ident", 'nCount_RNA', "nFeature_RNA", "percent.mt", 'disease_group2', 'cell_type_1', "cell_type_4")] %>% mutate("cell_barcode" = rownames(.))

  # gene selection: expressed in >50% donors, expressed in > 5% of cells
  # expression per donor -> use matrix operation
  donorMatrix = reshape2::dcast(metadata.tmp[,c('orig.ident', 'cell_barcode')], formula = orig.ident~cell_barcode)
  donorMatrix[,2:ncol(donorMatrix)] = ifelse(is.na(donorMatrix[,2:ncol(donorMatrix)]), 0, 1 )
  geMatrix = GetAssayData(ser.tmp, slot = 'counts') %>% t()
  geMatrix = (geMatrix>0)
  if(!all(rownames(geMatrix) %in% colnames(donorMatrix))){ print("cell names do not match! Check metadata and counts slot"); break}
  if(!nrow(geMatrix)==(ncol(donorMatrix)-1)){ print("cell count does not match! Check metadata and counts slot"); break}

  donorMatrix = donorMatrix[,c("orig.ident", rownames(geMatrix))] %>% column_to_rownames("orig.ident") %>% as.matrix() %>% Matrix(., sparse=T)
  if(!all(rownames(geMatrix) %in% colnames(donorMatrix))){ print("Second check! cell names do not match! Check metadata and counts slot"); break}
  donorMatrix2 = donorMatrix %*% geMatrix
  selected.genes = colnames(donorMatrix2)[colSums(donorMatrix2 != 0) > (length(included.samples)/2)]
  # expressed in >5% of total cells
  common.genes = intersect(selected.genes,
                           colnames(geMatrix)[colSums(geMatrix) > (0.05*ncol(ser.tmp))])
  saveRDS(common.genes, paste0(save.dir, "/input_data/common_genes_", ser.tmp@meta.data$cell_type_1[1], ".rds"))

  cell.order = colnames(ser.tmp)
  saveRDS(cell.order, paste0(save.dir, "/input_data/cell_order_", ser.tmp@meta.data$cell_type_1[1], ".rds"))
  all(rownames(metadata.tmp) == colnames(GetAssayData(ser.tmp, slot='count')[common.genes,]))

  # Write GE df
  ge.df = cbind(metadata.tmp, GetAssayData(ser.tmp, slot='count')[common.genes,] %>% t()) # since same order use cbind rather than merge
  if(!all(rownames(ge.df) == cell.order)){print("ge df does not match cell order! check!"); break}
  fwrite(ge.df, paste0(save.dir, '/input_data/ge.df.', ser.tmp@meta.data$cell_type_1[1], '.txt'), sep='\t')

  # covariates matrix
  ser.tmp <- FindVariableFeatures(ser.tmp)
  ser.tmp <- ScaleData(ser.tmp)
  ser.tmp <- RunPCA(ser.tmp, npcs=20)

  pca.matrix = ser.tmp@reductions$pca@cell.embeddings[,1:15]

  all(rownames(pca.matrix) == cell.order)
  all(rownames(metadata.tmp) == cell.order)


  covariates = cbind(pca.matrix, metadata.tmp)
  covariates = merge(covariates, clininfo, by.x='orig.ident', by.y="FL_ID", all.x=T) %>% column_to_rownames("cell_barcode")
  covariates$biopsy_age = scale(covariates$biopsy_age)
  covariates$nCount_RNA = scale(log10(covariates$nCount_RNA))

  if(sum(!cell.order %in% rownames(covariates)) !=0 ){ print("Covariates df is missing some cells! check!"); break}
  if(sum(!rownames(covariates) %in% cell.order) !=0 ){ print("Covariates df has some spurious cells! check!"); break}
  covariates = covariates[cell.order, ]

  fwrite(covariates %>% mutate("cell_barcode" = rownames(.)),
         paste0(save.dir, 'input_data/covariates.', ser.tmp@meta.data$cell_type_1[1], '.txt'), sep='\t')

  if(i==1){common.genes.merged = common.genes}
  else {common.genes.merged = c(common.genes.merged, common.genes) %>% unique()}
  gc()
}


# filter SNPs
process_raw_gtf2 <- function(raw_gtf){
  processed = raw_gtf %>% dplyr::select(1,3,4,5,7, 9) %>% filter(V3 == "gene") %>%
    separate(., V9, into = c("pre", "after"), sep = 'gene_name ') %>%
    dplyr::select(-pre) %>%
    separate(., after, into = c("gene_name", "transcript_type"), sep = '; ', extra = "drop") %>%
    dplyr::select(-transcript_type, -V3)
  colnames(processed) = c("#Chr","old_start", "old_end", "strand", "ID")
  processed$start = processed$old_start
  processed$end = processed$old_end
  processed[which(processed$strand == "-"),"start"] = processed[which(processed$strand == "-"), "old_end"]
  processed[which(processed$strand == "-"),"end"] = processed[which(processed$strand == "-"), "old_start"]
  processed = processed[,c("#Chr", "start", "end", "ID")]
  head(processed)

  return(processed)
}



find_cis_snps_for_egene <- function(egene_input,
                                    gt_raw,
                                    cis_length=1e6,
                                    gtf_processed){
  egene_data = gtf_processed %>% filter(ID == egene_input)
  egene_chr = gsub("chr", "", egene_data[["#Chr"]])
  egene_tss = egene_data[["start"]]
  esnps= gt_raw$map %>% filter(chromosome == egene_chr) %>% filter( position > egene_tss - cis_length ) %>% filter( position < egene_tss + cis_length) %>% rownames(.)
  return(esnps)
}

plink_path=paste0(save.dir, "plink_prefix")
gt_raw <- read.plink(plink_path)
gt <- 2- data.frame(as(gt_raw$genotypes,'numeric'), check.names = F)
gt$orig.ident = rownames(gt)

cr38 <- read.table(paste0(save.dir, "cellranger_grch38_2020_A_annotations.gtf"), skip=5, sep='\t')
cr38_processed <- process_raw_gtf2(cr38)

all.cis.snps <- pbmclapply(1:length(common.genes.merged),
                           FUN = function(x){ find_cis_snps_for_egene(common.genes.merged[x],
                                                                      gt_raw, gtf_processed=cr38_processed) },
                           mc.cores=20, mc.preschedule=F)
# 3 genes with duplicate gene names in seurat object -> some are retained in common.genes.merged (MATR3 (MATR3.1 only), TBCE (TBCE.1 both), HSPA14 (HSPA14.1 only ))(two genes @ almost same position)
# those retained in common.genes.merged mostly do not have the same ID in gtf file causing error message in find_cis_snps_for_egene function -> error message is output, rather than list of eSNPs
# error message should be removed for cis.snps.merged

all.cis.snps.m <- Reduce(union, all.cis.snps)
length(all.cis.snps.m) # 589872

gene_snp_pairs = all.cis.snps
names(gene_snp_pairs) = common.genes.merged


saveRDS(all.cis.snps, paste0(save.dir, "/input_data/all.cis.snps.list.rds"))
saveRDS(gene_snp_pairs, paste0(save.dir, "/input_data/gene_snp_pairs.rds"))
saveRDS(all.cis.snps.m, paste0(save.dir, "./input_data/all.cis.snps.m.rds"))
saveRDS(cr38_processed, paste0(save.dir, "/input_data/cr38_processed.rds"))

# sample level gt matrix
# 3 genes with duplicate gene names in seurat object -> some are retained in common.genes.merged (MATR3 (MATR3.1 only), TBCE (TBCE.1 both), HSPA14 (HSPA14.1 only ))(two genes @ almost same position)
# those retained in common.genes.merged mostly do not have the same ID in gtf file causing error message in find_cis_snps_for_egene function -> error message is output, rather than list of eSNPs
# error message should be removed for all.cis.snps.m object 
all.cis.snps.m = all.cis.snps.m[-41468] #Error in filter(., chromosome == egene_chr) : \n  \033[1m\033[22mProblem while computing `..1 = chromosome == egene_chr`.\n\033[31m✖\033[39m Input `..1` must be of size 624678 or 1, not size 0.\n" (d/t duplicate gene names of 3 genes)
gt.f = gt[,all.cis.snps.m]
fwrite(gt.f, paste0(save.dir, "/input_data/gt.df.bulk.txt", sep='\t', row.names=T))
