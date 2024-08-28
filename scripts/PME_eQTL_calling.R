
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
chunk_i=as.numeric(args[2]) # starts at 1
celltype=args[3]


print(paste0("n_cores: ", n.cores, " / chunk: ", chunk_i, " / celltype: ", celltype))
print("Chunk should start from 1 !")



ge.df <- fread(paste0(save.dir, '/input_files/ge.df.', celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
gt.df <- fread(paste0(save.dir, "/input_files/gt.df.bulk.txt"), sep='\t') %>% column_to_rownames("V1")
covariates <- fread(paste0(save.dir, "/input_files/covariates.", celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")

common.genes = readRDS(paste0("../input_files/common_genes_", celltype, ".rds"))
gene_snp_pairs <- readRDS("../input_files/gene_snp_pairs.rds")

cis.snps.count = lapply(gene_snp_pairs, FUN=function(x){length(x)}) %>% unlist() %>% data.frame(., check.names = F)
colnames(cis.snps.count) = "eSNP_n"
cis.snps.count = cis.snps.count[common.genes,,drop=F] %>% filter(eSNP_n > 1) %>% arrange(eSNP_n) %>% filter(eSNP_n > 1)

common.genes_sorted = rownames(cis.snps.count)
ge.df = ge.df[,c(colnames(ge.df)[1:7], common.genes_sorted)]


chunk_size=500
max_chunk_i = (length(common.genes_sorted) %/% chunk_size) + 1

start.gene.index = chunk_size*(chunk_i -1 ) + 1
end.gene.index = chunk_size*(chunk_i)
if(chunk_i==max_chunk_i){ end.gene.index = length(common.genes_sorted )}

ge.df = ge.df[,c(1:7, (7+start.gene.index):(7+end.gene.index))]

pme_gt_disease <- function(ge, gt.b, covs, gene, snps){
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

    # NAFLD cells only
    nafld.cells = rownames(ge.df %>% filter(disease_group2 != 'ctrl'))
    E <- ge[nafld.cells,gene]
    G <- gt[nafld.cells, snp]
    covariates <- covs[nafld.cells,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", "disease_group3", paste0("PC_", 1:5))] # "disease_group_num", scwgcna modules

    data.test <- data.frame(E,G,covariates)
    # print("df ready")
    full_fm = paste0("E~ G + (1|orig.ident) + biopsy_age + Sex + nCount_RNA + ", paste0("PC_", 1:5,collapse="+ "))

    full_model <- lme4::glmer(formula=as.formula(full_fm),
                              family='poisson', nAGQ=0, data= data.test %>% filter(!is.na(G)) , control = lme4::glmerControl(optimizer='nloptwrap'))
    full_summary = summary(full_model)$coefficients
    colnames(full_summary) = paste0(c("Estimate", "Std_Error", "z_value", "pval"), '_NAFLD')

    nafld.output = data.frame(gene = gene, snp = snp,
                              term = rownames(full_summary),  full_summary) %>% remove_rownames()

    # ctrl cells only
    ctrl.cells = rownames(ge.df %>% filter(disease_group2 == 'ctrl'))
    E <- ge[ctrl.cells,gene]
    G <- gt[ctrl.cells,snp]
    covariates <- covs[ctrl.cells,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", "disease_group3", paste0("PC_", 1:5))] # "disease_group_num", scwgcna modules

    data.test <- data.frame(E,G,covariates)

    full_fm = paste0("E~ G + (1|orig.ident) + biopsy_age + Sex + nCount_RNA + ", paste0("PC_", 1:5,collapse="+ "))

    full_model <- lme4::glmer(formula=as.formula(full_fm),
                              family='poisson', nAGQ=0, data= data.test %>% filter(!is.na(G)) , control = lme4::glmerControl(optimizer='nloptwrap'))

    full_summary = summary(full_model)$coefficients
    colnames(full_summary) = paste0(c("Estimate", "Std_Error", "z_value", "pval"), "_ctrl")
    ctrl.output = data.frame(gene = gene, snp = snp,
                             term = rownames(full_summary),  full_summary) %>% remove_rownames()


    output.merged = merge(nafld.output, ctrl.output, by=c('gene', 'snp', 'term'), suffixes=c("_NAFLD", "_ctrl"), all=T) # all=T
    output.merged = merge(all.output, output.merged, by=c("gene", 'snp', 'term'), all=T) # all=T

    loopn= loopn+1
    if(loopn%%100 == 0) {print(paste0("gene: ", gene, " / processed snp #: ", loopn))}
    if(loopn%%200 == 0) {gc()}
    if(snp == snps[1]){ output_m = output.merged}
    else{ output_m = rbind(output_m, output.merged)}
  }
  return(output_m)
}


print(paste0("Start eQTL calling from ", start.gene.index, " - ", end.gene.index, " th gene"))
print(paste0("first gene : " ,common.genes_sorted[start.gene.index]))
print(paste0("Last gene : " ,common.genes_sorted[end.gene.index]))

res_tmp <- pbmclapply( start.gene.index:end.gene.index,
                       FUN = function(x){
                         pme_gt_disease(ge = ge.df, gt = gt.df, covs=covariates,
                                                 gene=common.genes_sorted[x], snps=gene_snp_pairs[[ common.genes_sorted[x] ]]) },
                       mc.cores=n.cores, mc.preschedule=F)

res_tmp2 <- do.call('rbind', res_tmp)
write.csv(res_tmp2, paste0(save.dir, "./results/", celltype, "_chunk_",chunk_i, ".csv"))

end_time = Sys.time()
print(end_time)
time_elapsed = end_time -start_time
print(time_elapsed)



# @ bash
# for i in {1..19}
#do
#  Rscript /scripts/PME_eQTL_calling.R 60 $i 'hepatocyte'
#done
