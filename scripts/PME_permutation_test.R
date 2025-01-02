suppressMessages({
  library(lme4)
  library(Matrix)
  library(dplyr)
  library(lmtest)
  library(tidyverse)
  library(pbmcapply)
  library(data.table)})

merged_eqtl_res_dir = "..."
output_dir='./permutation/'
celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')


# obtain the list of significant eQTLs
for (celltype_cur in celltypes){
  sig.eqtls <- read.csv(paste0("../pme.sig.0.05.", celltype_cur, ".csv"))
  sig.eqtls = sig.eqtls %>% filter(is_liver_eQTL) 
  write.table(sig.eqtls, paste0(output_dir, 'sig_and_all_gene_snp_pair_', celltype_cur, ".txt"), row.names=F)
  if(celltype_cur==celltypes[[1]]) { top_snps_merged = sig.eqtls$snp}
  else{ top_snps_merged = c(top_snps_merged, sig.eqtls$snp) %>% unique}
}
saveRDS(top_snps_merged, paste0(output_dir, "sig_and_all_snps_merged.rds"))


# permutate 200 times per gene-snp pair, calculate its empirical p value
pme_gt_disease_permute <- function(ge, gt.b, covs, gene, snps){
  # check if rownames (=cell order) are the same
  print(paste0("Testing: ", gene, "; # of SNPs: ", length(snps)))
  if(!all(rownames(ge) == rownames(covs))){stop("GE vs covs rownames do not match")}
  
  # set cell order 
  cell.order.f= rownames(ge)
  gt <- merge(covs[,c("orig.ident"), drop=F]%>% mutate("cell_barcode"=rownames(.)),
              gt.b[,snps, drop=F], by.x='orig.ident',by.y=0, all.x=T) %>% column_to_rownames("cell_barcode")
  gt = gt[cell.order.f,]
  
  loopn=0
  for(snp in snps){
    # all LRT
    E <- ge[,gene]
    G <- gt[,snp]
    covariates <- covs[,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", "disease_group3", paste0("PC_", 1:5))] # "disease_group_num", scwgcna modules
    
    data.test <- data.frame(E,G,covariates)
    null_fm = paste0("E~ G + disease_group3 + (1|orig.ident) + biopsy_age + Sex + nCount_RNA +",paste0("PC_", 1:5,collapse="+ "))
    null_model <- lme4::glmer(formula=as.formula(null_fm),
                              family='poisson', nAGQ=0, data= data.test %>% filter(!is.na(G)) , control = lme4::glmerControl(optimizer='nloptwrap'))
    
    null_summary = summary(null_model)$coefficients
    colnames(null_summary) = paste0(c("Estimate", "Std_Error", "z_value", "pval"), "_null")
    all.output = data.frame(gene = gene, snp = snp,
                            term = rownames(null_summary),  null_summary) %>% remove_rownames()
  }
  return(all.output)
}

###########################
# permutation for liver-eQTLs
##########################
input_dir="./input_data/"
merged_eqtl_res_dir = "..."
output_dir='./permutation/'
celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')

n.cores=14
subset.gt.df <- readRDS(paste0(output_dir, "gt.df.bulk.rds"))
sample.names=rownames(subset.gt.df)
set.seed(123)

permuted.gt.df = apply(subset.gt.df, 2, function(column) sample(column)) %>% 
  as.data.frame(., check.names=F, row.names = sample.names) 


for(celltype in celltypes){
  ge.df <- fread(paste0(input_dir, '/ge.df.', celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
  covariates <- fread(paste0(input_dir, "covariates.", celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
  eqtls_of_interest = read.table(paste0(output_dir, 'sig_and_all_gene_snp_pair_', celltype, ".txt"), header=1)
  print(paste0("eqtls tested: ", nrow(eqtls_of_interest)))
  gc()
  set.seed(123)
  for(permut_i in 1:200){
    print(paste0('celltype: ', celltype, "/ permut_version: ", permut_i))
    permuted.gt.df = apply(permuted.gt.df, 2, function(column) sample(column)) %>% 
      as.data.frame(., check.names=F, row.names = sample.names) 
    
    print(paste0("permut_i: ", permut_i))
    start_time = Sys.time()
    res_tmp <- mclapply(1:nrow(eqtls_of_interest), 
                        FUN = function(x){
                          pme_gt_disease_permute(ge = ge.df, gt = permuted.gt.df, covs=covariates,
                                                 gene=eqtls_of_interest[x, 'gene'], snps=eqtls_of_interest[x, 'snp']) },
                        mc.cores=n.cores, mc.preschedule=F)
    res_tmp2 <- do.call('rbind', res_tmp)
    write.csv(res_tmp2, paste0(output_dir, "/sc_eQTL_permute_v", permut_i, "_", celltype, "_merged.csv"))
    
    end_time = Sys.time()
    print(end_time)
    time_elapsed = end_time -start_time
    print(time_elapsed)
    
  }
}



###########################
# permutation for ieQTLs
##########################
# to reduce memory, subset genotypes
input_dir="./input_data/"
gt.df <- fread(paste0(input_dir, "/gt.df.bulk.txt"), sep='\t') %>% column_to_rownames("V1")
ieqtls_snps = (pme.sig.anno %>% filter(is_ieQTL_cholangiocyte | is_ieQTL_hepatocyte | is_ieQTL_stellate_cell | is_ieQTL_endothelial_cell))$snp %>% unique()
gt.df.subset = gt.df[,ieqtls_snps]
saveRDS(gt.df.subset, paste0(output_dir, "gt.df.bulk.subset_ieQTLs.rds"))

input_dir="./input_data/"
output_dir='./permutation/'
celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')
pme_results_dir="..."

suppressMessages({
  library(lme4)
  library(Matrix)
  library(dplyr)
  library(lmtest)
  library(tidyverse)
  library(pbmcapply)
  library(data.table)})

start_time = Sys.time()
print(start_time)


pme_gt_state_v2 <- function(ge, gt, covariates_raw, gene, snp, phenotype){
  cell.order.f= rownames(ge)
  gt <- merge(covariates_raw[,c("orig.ident"), drop=F]%>% mutate("cell_barcode"=rownames(.)),
              gt[,snp, drop=F], by.x='orig.ident',by.y=0, all.x=T) %>% column_to_rownames("cell_barcode")
  G = gt[cell.order.f,snp]
  E <- ge[,gene]
  
  covariates <-covariates_raw[,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", paste0("PC_", 1:5))] # "disease_group_num", scwgcna modules
  pheno <- covariates_raw[,phenotype]
  
  data.test <- data.frame(E,G,covariates, pheno)
  null_fm = paste0("E~G + (1|orig.ident) + biopsy_age + Sex + nCount_RNA + ", paste0("PC_", 1:5, collapse=" +"), "+ pheno")
  full_fm = paste0("E~G + (1|orig.ident) + biopsy_age + Sex + nCount_RNA + ",paste0("PC_", 1:5, collapse=" +"), "+ pheno + G*pheno")
  
  full_model <- lme4::glmer(formula=as.formula(full_fm),
                            family='poisson', nAGQ=0, data= data.test %>% filter(!is.na(G)) , control = lme4::glmerControl(optimizer='nloptwrap'))
  null_model <- lme4::glmer(formula=as.formula(null_fm),
                            family='poisson', nAGQ=0, data= data.test %>% filter(!is.na(G)) , control = lme4::glmerControl(optimizer='nloptwrap'))
  model_lrt <- anova(null_model, full_model)
  model_lrt$`Pr(>Chisq)`[2]
  full_summary = summary(full_model)$coefficients
  colnames(full_summary) = paste0(c("Estimate", "Std_Error", "z_value", "pval"), "_full")
  full_output = data.frame(gene = gene, snp = snp,
                           term = rownames(full_summary),  full_summary) %>% remove_rownames()
  
  null_summary = summary(null_model)$coefficients
  colnames(null_summary) = paste0(c("Estimate", "Std_Error", "z_value", "pval"), "_null")
  null_output = data.frame(gene = gene, snp = snp,
                           term = rownames(null_summary),  null_summary) %>% remove_rownames()
  gc()
  output =  merge(full_output, null_output, by=c('gene', 'snp', 'term'), all=T) %>% mutate(lrt_pval=model_lrt$`Pr(>Chisq)`[2],
                                                                                           'Phenotype' = phenotype )
  return(output)
}



wgcna.filename.list = list("hepatocyte" = "hep", "cholangiocyte"='chol', 'stellate_cell' = "hsc", "endothelial_cell" = 'endo')
pme.sig.tmp = readRDS(paste0(analysis_results_dir, "pme.sig.anno.full.rds"))

for(celltype in celltypes){
  ge.df <- fread(paste0(pme_results_dir, 'input_data/ge.df.', celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
  gt.df <- readRDS(paste0(output_dir, "gt.df.bulk.subset_ieQTLs.rds"))
  covariates.df <- fread(paste0(pme_results_dir, "input_data/covariates.", celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
  
  wgcna.module.scores <- readRDS(paste0(pme_results_dir, "input_data/", wgcna.filename.list[[celltype]] %>% as.character, ".ser_hMEs.rds"))
  wgcna.module <- readRDS(paste0(pme_results_dir, "input_data/", wgcna.filename.list[[celltype]] %>% as.character,".ser_modules.rds"))[,2:3] %>% filter(!duplicated(module,color)) %>% remove_rownames()
  wgcna.module.list <- as.character(wgcna.module$module)
  names(wgcna.module.list) <- as.character(wgcna.module$color)
  colnames(wgcna.module.scores) <- wgcna.module.list[colnames(wgcna.module.scores)] %>% unname()
  wgcna.module.scores <- scale(wgcna.module.scores)
  wgcna.module.scores <- wgcna.module.scores[rownames(covariates),]
  wgcna.module.count = list('hepatocyte'=c(3,4,6,8,12,14), 'cholangiocyte'=c(1,2,3,4,6,7), 'stellate_cell' = 1:3, 'endothelial_cell' = c(1,4,7,9))
  selected.modules =paste0(wgcna.module.name[[celltype]], wgcna.module.count[[celltype]])
  
  wgcna.module.scores = wgcna.module.scores[,selected.modules]
  if(all(rownames(wgcna.module.scores) == rownames(covariates.df))){
    covariates.df <- cbind(covariates.df, wgcna.module.scores)
  } else{ 
    stop("ERROR! CHECK for cell order in covariates.df and wgcna.module.scores")}
  
  print(paste0("cell order (GE vs Cov): ", all(rownames(ge.df) == rownames(covariates.df))))
  
  for(phenotype_cur in c(selected.modules, 'disease_group3')){
    set.seed(123)
    if(phenotype_cur=='disease_group3') { eqtls_of_interest = pme.sig.tmp %>% filter(!is.na (!!sym(paste0('disease_interaction_', celltype)))) }
    else{ eqtls_of_interest = pme.sig.tmp %>% filter(!is.na (!!sym(phenotype_cur)) )}
    
    # random shuffling of phenotypes
    if(phenotype_cur=='disease_group3') { 
      donor_disease_matrix = covariates.df[,c('orig.ident', 'disease_group3')] %>% unique
      donor_disease_matrix$disease_group3 = sample(donor_disease_matrix$disease_group3)
      covariates.df = covariates.df %>% select(-disease_group3) %>% rownames_to_column("cell_barcode") %>% 
        merge(., donor_disease_matrix, by='orig.ident', all=T) %>% column_to_rownames('cell_barcode')
      covariates.df = covariates.df[rownames(wgcna.module.scores), ] }
    
    else{ 
      covariates.df[,phenotype_cur] = sample(covariates.df[,phenotype_cur])
    }
    
    for(permut_i in 1:200){
      print(paste0("permut_i: ", permut_i))
      start_time = Sys.time()
      res_tmp <- pbmclapply(1:nrow(eqtls_of_interest), 
                            FUN = function(x){ pme_gt_state_v2(ge=ge.df,
                                                               gt=gt.df,
                                                               covariates_raw=covariates.df,
                                                               gene=eqtls_of_interest[x,"gene"],
                                                               snp=eqtls_of_interest[x, "snp"],
                                                               phenotype=phenotype_cur)},
                            mc.cores=n.cores, mc.preschedule=F)
      res_tmp2 <- do.call('rbind', res_tmp)
      write.csv(res_tmp2, paste0(output_dir, "/ieQTLs/.", celltype, "_", phenotype_cur,"_", "permut_", permut_i,  ".csv"))
      
      end_time = Sys.time()
      print(end_time)
      time_elapsed = end_time -start_time
      print(time_elapsed)
      
    }
  }
}


