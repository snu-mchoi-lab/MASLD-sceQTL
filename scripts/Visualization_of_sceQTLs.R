# code for eQTL/ieQTL plots

library(scales)
library(tidyverse)
library(dplyr)
library(data.table)
library(reshape2)
library(patchwork)
library(grid)
library(gridExtra)

#############################################################################################
#############################################################################################
##########################         PLOT Liver-eQTLs       ###################################
#############################################################################################
#############################################################################################

# It is recommended to set ge.df , gt.df, covariates.df, eqtl.df variables to the corresponding data frames
# If it is not set, below function will read files from input.dir and analysis.dir

plot.liver.eqtls <- function(celltype, gene, snp, 
                             ge.df = NA, gt.df=NA, covariates.df = NA, eqtl.df=NA, 
                             custom_colors=c("burlywood1",'coral', 'coral3'),
                             input.dir=input_dir,
                             analysis.dir=input_dir.2){
  
  if(is.na(ge.df)){ge.df <- fread(paste0(input.dir, 'ge.df.', celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")}
  if(is.na(gt.df)){gt.df <- fread(paste0(input.dir, "gt.df.bulk.txt")) %>% column_to_rownames("V1")}
  if(is.na(covariates.df)){covariates.df <- fread(paste0(input.dir, "covariates.", celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")}
  if(is.na(eqtl.df)){eqtl.df <- fread(paste0(analysis.dir, "liver_eQTLs_summary_stats_", celltype, ".txt.gz"))}
  cell.order.f= rownames(ge.df)
  gt.df <- merge(covariates.df[,c("orig.ident"), drop=F]%>% mutate("cell_barcode"=rownames(.)),
                 gt.df[,snp, drop=F], by.x='orig.ident',by.y=0, all.x=T) %>% column_to_rownames("cell_barcode")
  G = gt.df[cell.order.f,snp]
  E <- ge.df[,gene]
  
  covariates <-covariates.df[,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", paste0("PC_", 1:5), "disease_group2", 'disease_group3')]
  
  data.test <- data.frame(E,G,covariates, check.names=F)
  
  tmp = data.test %>% group_by(orig.ident,disease_group3) %>% summarise(mean_gene = mean(E),
                                                                        mean_log2_1 = mean( log2(E + 1)))
  tmp = merge(tmp, data.test[,c("orig.ident", "G", "disease_group2")] %>% unique(),
              by='orig.ident', all=T)
  
  betaG = (eqtl.df %>% filter(gene == !!gene) %>% filter(snp == !!snp) )$Beta_all_cells
  P = (eqtl.df %>% filter(gene == !!gene) %>% filter(snp == !!snp) )$pval_all_cells
  
  tmp$disease_group2= factor(tmp$disease_group2, levels=c("ctrl", "NAFL", "eNASH", "aNASH"))
  
  p<- ggplot(tmp %>% filter(!is.na(G)), aes(y=mean_log2_1)) + 
    geom_jitter(height=0, width=0.2, aes(x=as.factor(G)), alpha=0.5)+ 
    theme_bw() + 
    scale_fill_manual(values=custom_colors, name='Genotype') +
    geom_boxplot(aes(x=as.factor(G), y=mean_log2_1, fill=as.factor(G)), alpha=0.7, outlier.shape=NA, width=0.6) +
    geom_smooth(aes(x=G+1, y=mean_log2_1), method = 'lm', se = F, color='grey50', size=0.5)+
    ylab(paste0("mean expression: ", gene)) + xlab(paste0("Genotype: ", snp))+
    ggtitle(paste0(gene, "_", snp)) + 
    labs(subtitle=paste0(celltype,  " /b=",round(betaG, 2)," /p=",scientific(P, 2)))
  return(list(p, tmp, betaG, P)) # returns a list of plot[[1]], and the data used for plotting (elements 2-4)
}

# EXAMPLE

input_dir = "..."
input_dir.2 = "./sceQTL_output/liver_eQTLs/"

celltype_cur='hepatocyte'
ge.df.raw <- fread(paste0(input_dir, 'ge.df.', celltype_cur, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
gt.df.raw <- fread(paste0(input_dir, "gt.df.bulk.txt")) %>% column_to_rownames("V1")
covariates.df.raw <- fread(paste0(input_dir, "covariates.", celltype_cur, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
eqtl.df.raw <- fread(paste0(input_dir.2, "liver_eQTLs_summary_stats_", celltype_cur, ".txt.gz"))

genes_cur=c("MLIP","AC091114.1","CYP2D6")
snps_cur=c("chr6:54092790:TC", "chr8:125479062:TC", "chr22:42129754:GA")
plot_list <- vector('list', length(genes_cur))

for(i in 1:length(genes_cur)){ #3 
  gene_cur=genes_cur[i]
  snp_cur=snps_cur[i]
  plot_list[[i]] <- plot.liver.eqtls(celltype=celltype_cur, gene=gene_cur, snp=snp_cur, 
                                     ge.df = ge.df.raw, 
                                     gt.df=gt.df.raw,
                                     covariates.df = covariates.df.raw, 
                                     eqtl.df=eqtl.df.raw)[[1]] 
}




#############################################################################################
#############################################################################################
##########################           PLOT ieQTLs          ###################################
#############################################################################################
#############################################################################################

# It is recommended to set wgcna.module.scores, ge.df , gt.df, covariates.df ,interaction.results, variables to the corresponding data frames
# If it is not set, below function will read files from input.dir and analysis.dir

plot.wgcna.ieqtls.pretty.v2 <- function(celltype, phenotype, gene, snp, 
                                        wgcna.module.scores=NA, ge.df = NA, gt.df=NA, covariates.df = NA, interaction.results=NA, # if these variables are set to corresponding dataframes, input.dir and analysis.dir is not required
                                        input.dir=input_dir, # should be set to directory with eQTL input files
                                        analysis.dir=input_dir.2, # should be set to directory with ieQTL summary statistics files
                                        bin_by_cell=T, # if True, cells are binned by module expression value regardless of donor / if FALSE, cells are binned by module expession within each donor
                                        nbin=3, # how many bins of module expression levels to split into
                                        sample_filter_column=NA, sample_filter_value=NA){ # if you want to plot with only selected samples, specify the column name and value for the samples. It should be included in the covariates matrix.
  if(is.na(ge.df)){ge.df <- fread(paste0(input.dir, 'ge.df.', celltype, ".txt"), sep='\t')}
  if(is.na(gt.df)){gt.df <- fread(paste0(input.dir, "gt.df.bulk.txt")) %>% column_to_rownames("V1")}
  if(is.na(covariates.df)){covariates.df <- fread(paste0(input.dir, "covariates.", celltype, ".txt"), sep='\t') }
  if(is.na(interaction.results)){ interaction.results <- fread(paste0(analysis.dir, "interaction.", celltype, "_", phenotype, ".csv"))}
  if(is.na(wgcna.module.scores)){
    wgcna.module.scores <- readRDS(paste0(input.dir, wgcna.filename.list[[celltype]] %>% as.character, ".ser_hMEs.rds"))
    wgcna.module <- readRDS(paste0(input.dir, wgcna.filename.list[[celltype]] %>% as.character,".ser_modules.rds"))[,2:3] %>% filter(!duplicated(module,color)) %>% remove_rownames()
    wgcna.module.list <- as.character(wgcna.module$module)
    names(wgcna.module.list) <- as.character(wgcna.module$color)
    colnames(wgcna.module.scores) <- wgcna.module.list[colnames(wgcna.module.scores)] %>% unname()
    wgcna.module.scores <- wgcna.module.scores[covariates.df$cell_barcode,] }
  
  wgcna.module.scores <- wgcna.module.scores[covariates.df$cell_barcode,]
  
  if(all(rownames(wgcna.module.scores) == covariates.df$cell_barcode)){
    covariates.df <- cbind(covariates.df, wgcna.module.scores) %>% data.frame(., check.names=F)
  } else{ stop("ERROR! CHECK for cell order in covariates.df and wgcna.module.scores") }
  
  cell.order.f= ge.df$cell_barcode
  gt <- merge(covariates.df[,c("orig.ident", "cell_barcode")],
              gt.df[,snp, drop=F], by.x='orig.ident',by.y=0, all.x=T) %>% column_to_rownames("cell_barcode")
  gt = gt[cell.order.f,]
  
  E <- (ge.df %>% data.frame())[,which(colnames(ge.df) == gene)]
  G <- gt[,which(colnames(gt) == snp)]
  
  if(!is.na(sample_filter_column) & !(sample_filter_column %in% colnames(covariates.df))){
    covariates <-covariates.df[,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", paste0("PC_", 1:5), 'disease_group2', "disease_group", 'disease_group3', sample_filter_column)] 
  }
  else{
    covariates <-covariates.df[,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", paste0("PC_", 1:5), 'disease_group2', "disease_group", 'disease_group3')] 
  }
  
  pheno <- covariates.df[,which(colnames(covariates.df)==phenotype)]
  pheno <- scale(pheno)
  
  data.test <- data.frame(E,G,covariates, pheno, check.names=F)
  data.test$pheno <- as.numeric(data.test$pheno)
  
  if(bin_by_cell){
    data.test <- data.test %>% mutate(module_bin = ntile(pheno, n=nbin)) 
    tmp = data.test %>% mutate(mean_gene = E,
                               mean_log2_1 = mean( log2(E + 1)))
  }
  else{
    tmp <- data.test %>% group_by(orig.ident) %>% mutate(module_bin = ntile(pheno, n=nbin)) %>% 
      ungroup  %>% group_by(orig.ident, module_bin) %>% summarise(mean_gene = mean(E),
                                                                  mean_log2_1 = mean( log2(E + 1)))
  }
  if(!is.na(sample_filter_column)){
    tmp = merge(tmp, data.test[,c("orig.ident", "G", sample_filter_column)] %>% unique(),
                by='orig.ident', all=T)
    tmp = tmp %>% filter(!!sym(sample_filter_column) %in% sample_filter_value)
  }
  else{
    tmp = merge(tmp, data.test[,c("orig.ident", "G")] %>% unique(),
                by='orig.ident', all=T)}
  
  betaG = (interaction.results %>% filter(gene == !!gene) %>% filter(snp == !!snp) %>% filter(term=="G"))$Estimate_full
  betaGxP=(interaction.results %>% filter(gene == !!gene) %>% filter(snp == !!snp) %>% filter(term=="G:pheno"))$Estimate_full
  anovaP = (interaction.results %>% filter(gene == !!gene) %>% filter(snp == !!snp) %>% filter(term=="G:pheno"))$lrt_pval
  
  p<- ggplot(tmp %>% filter(!is.na(G)), aes(y=mean_log2_1, fill=as.factor(G))) + 
    geom_jitter(aes(x=as.factor(G), group=as.factor(G)), height=0, width=0.3, size=0.5, alpha=0.5) + 
    geom_boxplot( aes(x=as.factor(G), group=as.factor(G), fill=as.factor(G)), alpha=0.7, outlier.shape=NA, width=0.6) +
    geom_smooth(aes(x=(G+1), y=mean_log2_1), method='lm', color="grey50", fill='grey80', se=F) +
    facet_wrap(~module_bin, nrow=1) + ylab(paste0("mean expression: ", gene)) + 
    xlab(paste0("Genotype: ", snp)) +
    scale_fill_manual(values=c("burlywood1",'coral', 'coral3')) + theme_test() +
    ggtitle(paste0(gene, "_", snp)) + 
    labs(subtitle=paste0("celltype: ", celltype, " / interaction: ", phenotype, " / ", sample_filter_value, 
                         "\n betaG=",round(betaG, 2)," / betaGxP=",round(betaGxP,2)," / LRT p=",scientific(anovaP, 2)))
  return(list(p, tmp)) }


# celltypes=c("hepatocyte", 'cholangiocyte', 'stellate_cell', 'endothelial_cell')
celltype_cur='hepatocyte'
phenotype_cur="Hep-M12"
gene_cur="EFHD1"
snp_cur="chr2:232655544:AT"

input_dir = "./sceQTL_input/"
input_dir.2 = "./sceQTL_output/ieQTLs/"
input_dir.3 = './hdWGCNA_results/'

wgcna.filename.list = list("hepatocyte" = "hep", "cholangiocyte"='chol', 'stellate_cell' = "hsc", "endothelial_cell" = 'endo')
wgcna.module.name=list("hepatocyte" = 'Hep-M', 'cholangiocyte'= "Chol-M", "stellate_cell" = "HSC-M", 'endothelial_cell' = "endo-M")
wgcna.module.count = list('hepatocyte'=c(3,4,6,8,12,14), 'cholangiocyte'=c(1,2,3,4,6,7), 'stellate_cell' = 1:3, 'endothelial_cell' = c(1,4,7,9))

wgcna.module.scores.raw <- readRDS(paste0(input_dir.3, wgcna.filename.list[[celltype_cur]] %>% as.character, ".ser_hMEs.rds"))
wgcna.module.raw <- readRDS(paste0(input_dir.3, wgcna.filename.list[[celltype_cur]] %>% as.character,".ser_modules.rds"))[,2:3] %>% filter(!duplicated(module,color)) %>% remove_rownames()
wgcna.module.list.raw <- as.character(wgcna.module.raw$module)
names(wgcna.module.list.raw) <- as.character(wgcna.module.raw$color)
colnames(wgcna.module.scores.raw) <- wgcna.module.list.raw[colnames(wgcna.module.scores.raw)] %>% unname()

gt.df.raw <- fread(paste0(input_dir, "gt.df.bulk.txt")) %>% column_to_rownames("V1")
ge.df.raw <- fread(paste0(input_dir, 'ge.df.', celltype_cur, ".txt"), sep='\t')
covariates.df.raw <- fread(paste0(input_dir, "covariates.", celltype_cur, ".txt"), sep='\t')

interaction.results_raw <- fread(paste0(input_dir.2, "interaction.", celltype_cur, "_",phenotype_cur, ".txt.gz"))

p1 <- plot.wgcna.ieqtls.pretty.v2(celltype=celltype_cur, 
                                  phenotype=phenotype_cur, 
                                  gene=gene_cur, snp=snp_cur, 
                                  wgcna.module.scores=wgcna.module.scores.raw, 
                                  ge.df = ge.df.raw, gt.df=gt.df.raw, covariates.df = covariates.df.raw, 
                                  interaction.results=interaction.results_raw,
                                  nbin=3, bin_by_cell=F)[[1]]
p2 <- plot.wgcna.ieqtls.pretty.v2(celltype=celltype_cur, 
                                  phenotype=phenotype_cur, 
                                  gene=gene_cur, snp=snp_cur, 
                                  wgcna.module.scores=wgcna.module.scores.raw, 
                                  ge.df = ge.df.raw, gt.df=gt.df.raw, covariates.df = covariates.df.raw, 
                                  interaction.results=interaction.results_raw,
                                  nbin=3, bin_by_cell=F, sample_filter_column='disease_group3', sample_filter_value=c("NAFLD"))[[1]]
p3<- plot.wgcna.ieqtls.pretty.v2(celltype=celltype_cur, 
                                 phenotype=phenotype_cur, 
                                 gene=gene_cur, snp=snp_cur, 
                                 wgcna.module.scores=wgcna.module.scores.raw, 
                                 ge.df = ge.df.raw, gt.df=gt.df.raw, covariates.df = covariates.df.raw, 
                                 interaction.results=interaction.results_raw,
                                 nbin=3, bin_by_cell=F, sample_filter_column='disease_group3', sample_filter_value=c("ctrl"))[[1]]
plot_grid(p1,p2,p3, nrow=3)
ggsave(paste0("./interaction_eQTL_plots/", phenotype_cur, "_",gene_cur,"_", snp_cur, "_pretty.v2.pdf"), 
       width=4.5, height=10.5)




#############################################################################################
#############################################################################################
##########################    PLOT disease-ieQTLs         ###################################
#############################################################################################
#############################################################################################

plot.disease.ieqtls <- function(celltype, gene, snp, ge.df = NA, gt.df=NA, covariates.df = NA, interaction.results=NA,
                                input.dir=input_dir, # should be set to directory with eQTL input files
                                analysis.dir=input_dir.2){ # should be set to directory with ieQTL summary statistics files
  
  if(is.na(ge.df)){ge.df <- fread(paste0(input.dir, 'ge.df.', celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")}
  if(is.na(gt.df)){gt.df <- fread(paste0(input.dir, "gt.df.bulk.txt")) %>% column_to_rownames("V1")}
  if(is.na(covariates.df)){covariates.df <- fread(paste0(input.dir, "covariates.", celltype, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")}
  if(is.na(interaction.results)){ interaction.results <- fread(paste0(analysis.dir, "interaction.", celltype, "_disease_group3.txt.gz"))}
  
  cell.order.f= rownames(ge.df)
  gt.df <- merge(covariates.df[,c("orig.ident"), drop=F]%>% mutate("cell_barcode"=rownames(.)),
                 gt.df[,snp, drop=F], by.x='orig.ident',by.y=0, all.x=T) %>% column_to_rownames("cell_barcode")
  G = gt.df[cell.order.f,snp]
  E <- ge.df[,gene]
  
  covariates <-covariates.df[,c("orig.ident", "biopsy_age", "Sex", "nCount_RNA", paste0("PC_", 1:5), "disease_group2", 'disease_group3')] # "disease_group_num", scwgcna modules
  
  data.test <- data.frame(E,G,covariates, check.names=F)
  
  tmp = data.test %>% group_by(orig.ident,disease_group3) %>% summarise(mean_gene = mean(E),
                                                                        mean_log2_1 = mean( log2(E + 1)))
  tmp = merge(tmp, data.test[,c("orig.ident", "G", "disease_group2")] %>% unique(),
              by='orig.ident', all=T)
  
  betaG = (interaction.results %>% filter(gene == !!gene) %>% filter(snp == !!snp) %>% filter(term=="G"))$Estimate_full
  betaGxP=(interaction.results %>% filter(gene == !!gene) %>% filter(snp == !!snp) %>% filter(term=="G:phenoNAFLD"))$Estimate_full
  anovaP = (interaction.results %>% filter(gene == !!gene) %>% filter(snp == !!snp) %>% filter(term=="G:phenoNAFLD"))$lrt_pval
  
  tmp$disease_group2= factor(tmp$disease_group2, levels=c("ctrl", "NAFL", "eNASH", "aNASH"))
  
  p<- ggplot(tmp %>% filter(!is.na(G)), aes(x=as.factor(G), y=mean_log2_1)) + 
    geom_jitter(aes(x=as.factor(G), group=as.factor(G)), height=0, width=0.3, size=0.5, alpha=0.5) + 
    geom_boxplot( aes(x=as.factor(G), group=as.factor(G), fill=as.factor(G)), alpha=0.7, outlier.shape=NA, width=0.6) +
    geom_smooth(aes(x=(G+1), y=mean_log2_1), method='lm', color="grey50", fill='grey80', se=F) +
    facet_wrap(~disease_group3) + ylab(paste0("mean expression: ", gene)) + 
    xlab(paste0("Genotype: ", snp)) +
    scale_fill_manual(values=c("burlywood1",'coral', 'coral3')) + theme_test() +
    ggtitle(paste0(gene, "_", snp)) + labs(subtitle=paste0("celltype: ", celltype, " / interaction: disease",
                                                           "\n betaG=",round(betaG, 2)," / betaGxP=",round(betaGxP,2)," / ANOVA p=",scientific(anovaP, 2)))
  
  
  return(p)
}


celltype_cur='cholangiocyte'
gene_cur="BLNK"
snp_cur="chr10:95544975:CA"

input_dir = "./sceQTL_input"
input_dir.2 = "./sceQTL_output/ieQTLs/"

gt.df.raw <- fread(paste0(input_dir, "gt.df.bulk.txt")) %>% column_to_rownames("V1")
ge.df.raw <- fread(paste0(input_dir, 'ge.df.', celltype_cur, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")
covariates.df.raw <- fread(paste0(input_dir, "covariates.", celltype_cur, ".txt"), sep='\t') %>% column_to_rownames("cell_barcode")

interaction.results_raw <- fread(paste0(input_dir.2, "interaction.", celltype_cur,"_disease_group3.txt.gz"))

p <- plot.disease.ieqtls(celltype=celltype_cur, gene=gene_cur, snp = snp_cur, 
                         ge.df=ge.df.raw, gt.df = gt.df.raw, covariates.df = covariates.df.raw, interaction.results = interaction.results_raw)

p