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

args=commandArgs(trailingOnly = T)
n.cores=as.numeric(args[1])
celltype=as.character(args[2])
phenotype_cur=as.character(args[3])

wgcna.filename.list = list("hepatocyte" = "hep", "cholangiocyte"='chol', 'stellate_cell' = "hsc", "endothelial_cell" = 'endo')
# celltype='cholangiocyte'


# call ieQTLs
ge.df <- fread(paste0(save.dir, './input_data/ge.df.', celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
gt.df <- fread(paste0(save.dir, "./input_data/gt.df.bulk.txt"), sep='\t') %>% column_to_rownames("V1")
covariates.df <- fread(paste0(save.dir, "./input_data/covariates.", celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
sig_gene_snp_pairs <- fread(paste0(save.dir, "./analysis/pme.sig.0.05.", celltype, ".csv"))[,c('gene', 'snp')] %>% data.frame()
wgcna.module.scores <- readRDS(paste0(save.dir, "./input_data/", wgcna.filename.list[[celltype]] %>% as.character, ".ser_hMEs.rds"))
wgcna.module <- readRDS(paste0(save.dir, "./input_data/", wgcna.filename.list[[celltype]] %>% as.character,".ser_modules.rds"))[,2:3] %>% filter(!duplicated(module,color)) %>% remove_rownames()
wgcna.module.list <- as.character(wgcna.module$module)
names(wgcna.module.list) <- as.character(wgcna.module$color)
colnames(wgcna.module.scores) <- wgcna.module.list[colnames(wgcna.module.scores)] %>% unname()

wgcna.module.scores <- scale(wgcna.module.scores)
wgcna.module.scores <- wgcna.module.scores[rownames(covariates.df),]

if(!(phenotype_cur %in% c(colnames(wgcna.module.scores), colnames(covariates.df)))){
  stop("Phenotype not present! please check!")}

if(all(rownames(wgcna.module.scores) == rownames(covariates.df))){
  covariates.df <- cbind(covariates.df, wgcna.module.scores)
} else{
  stop("ERROR! CHECK for cell order in covariates.df and wgcna.module.scores")}


ge.df = ge.df[,c(colnames(ge.df)[1:7], sig_gene_snp_pairs$gene %>% unique)]
gt.df = gt.df[,sig_gene_snp_pairs$snp %>% unique]
print(paste0("cell order (GE vs Cov): ", all(rownames(ge.df) == rownames(covariates.df))))

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

res_tmp <- pbmclapply(1:nrow(sig_gene_snp_pairs),
                      FUN = function(x){ pme_gt_state_v2(ge=ge.df,
                                                                gt=gt.df,
                                                                covariates_raw=covariates.df,
                                                                gene=sig_gene_snp_pairs[x,"gene"],
                                                                snp=sig_gene_snp_pairs[x, "snp"],
                                                                phenotype=phenotype_cur)},
                      mc.cores=n.cores, mc.preschedule=F)
res_tmp2 <- do.call('rbind', res_tmp)
write.csv(res_tmp2, paste0(save.dir, "./results_AWS/interaction.", celltype, "_", phenotype_cur, ".csv"))


end_time = Sys.time()
print(end_time)
time_elapsed = end_time -start_time
print(time_elapsed)


# run @ bash
# for i in $(seq 7)
#do
#  Rscript /scripts_dell/211005_n48_A02/PME_ieQTL_calling.R 10 "cholangiocyte" 'Chol-M'$i
#done
