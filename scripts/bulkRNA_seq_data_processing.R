##################
##################
# EFHD1 siRNA treatment
##################
##################
counts_dir="./cowork_experiments/pf_SMJ/organoid_siRNA_RNA_seq/counts/"
analysis_dir="./cowork_experiments/pf_SMJ/organoid_siRNA_RNA_seq/analysis/"
genotype_df = data.frame('donor' = c("CW20009", "CW70221", "CRL2097", "CW10206", "CW10176", "CW60413"),
                         'genotype' = c(rep("AA", 3), rep("TT", 3)),
                         'gender' = c("F", "F", "M", "M", "F", "F"),
                         'age'= c(15, 42, 0, 32, 55, 30),
                         'ethnicity' = c(rep('Asian', 2), "Caucasian", "Asian", rep("Caucasian", 2)))

# merge counts 
i=1
for(counts_cur in list.files(counts_dir, '*.counts$') ){
  counts_tmp = fread(paste0(counts_dir, counts_cur))[,c(1,7)]
  colnames(counts_tmp)[2] = gsub("Aligned.out.bam", "", colnames(counts_tmp)[2])
  
  if(i==1) { counts_merged = counts_tmp}
  else{ counts_merged = merge(counts_merged, counts_tmp, by='Geneid')}
  i = i+1
}

fwrite(counts_merged,"raw_counts_merged.txt",sep='\t' )
counts_merged = counts_merged %>% column_to_rownames("Geneid")

# Retain samples with high-quality: n11
rnaseq.qc = read.csv(paste0(analysis_dir, "organoid_siEFHD1_FA_RNA_QC.csv"))
rnaseq.qc = rnaseq.qc %>% separate(col='cleanname', sep="_", into=c("donor", "rep", 'treatment'), remove=F)
rnaseq.qc = merge(rnaseq.qc, genotype_df, by='donor', all=T)
hiqc = rnaseq.qc %>% filter(quantity > 0.5) %>% filter(RIN > 5)


# make sample metadata
deseq_coldata = data.frame(sampleID = colnames(counts_merged)) %>% separate(col='sampleID', sep="_", into=c("donor", "rep", 'treatment'), remove=F)
deseq_coldata = merge(deseq_coldata, genotype_df, by='donor', all=T)
deseq_coldata = deseq_coldata %>% mutate(siRNA = case_when( treatment %in% c("siRNA", 'both') ~ T, TRUE ~ F))
deseq_coldata = deseq_coldata %>% column_to_rownames("sampleID")
deseq_coldata = deseq_coldata[colnames(counts_merged),]
deseq_coldata = deseq_coldata[hiqc$cleanname,]
deseq_coldata = deseq_coldata4 %>% mutate(FA = case_when( treatment %in% c("cont", 'siRNA') ~ F,treatment %in% c("FA", "both") ~ T,TRUE ~ NA))
deseq_coldata$FA = as.factor(deseq_coldata$FA)
deseq_coldata$gender = as.factor(deseq_coldata$gender)
deseq_coldata$genotype = as.factor(deseq_coldata$genotype)
deseq_coldata$siRNA = as.factor(deseq_coldata$siRNA)


# DESEQ normalization & DEG calculation
appropriate_gene_names = rownames(counts_merged)[!grepl("^ENSG", rownames(counts_merged))]
counts_merged = counts_merged[appropriate_gene_names ,hiqc$cleanname]

dds4 <- DESeqDataSetFromMatrix(countData = counts_merged, colData=deseq_coldata,design = ~ siRNA)
dds4 <- estimateSizeFactors(dds4)
keep4 <- rowSums(counts(dds4, normalized=TRUE) >= 10 ) >= 5.5
dds4 <- dds4[keep4,]
design(dds4) = ~ gender + age + FA + siRNA

dds4 <- DESeq(dds4)
dds4.deg.res = results(dds4, contrast = c('siRNA', TRUE, FALSE))
dds4.deg.res = data.frame(dds4.deg.res)


##################
##################
# FOXO1 inh in HepG2 
##################
##################
##############################################################
# ctrl vs (100nM samples + 1000nM samples)
##############################################################
# alignment: STAR / counting: featureCounts from subread
hepg2.analysis_dir="./hepg2_foxo1_inh/"
hepg2.counts.dir = paste0(hepg2.analysis_dir, "counts/") 

# merge counts 
i=1
for(counts_cur in list.files(hepg2.counts.dir, '*.counts$') ){
  counts_tmp = fread(paste0(hepg2.counts.dir, counts_cur))[,c(1,7)]
  colnames(counts_tmp)[2] = gsub("Aligned.out.bam", "", colnames(counts_tmp)[2])
  colnames(counts_tmp)[2] = gsub("../bam/", "", colnames(counts_tmp)[2])
  
  if(i==1) { hepg2.counts_merged = counts_tmp}
  else{ hepg2.counts_merged = merge(hepg2.counts_merged, counts_tmp, by='Geneid')}
  i = i+1
}

fwrite(hepg2.counts_merged,paste0(hepg2.analysis_dir, "analysis/hepg2_raw_counts_merged.txt"),sep='\t' )

# remove genes with no gene names (starts with ENSG)
hepg2.appropriate_gene_names = rownames(hepg2.counts_merged)[!grepl("^ENSG", rownames(hepg2.counts_merged))]
hepg2.counts_merged = hepg2.counts_merged[hepg2.appropriate_gene_names,]

# make sample metadata dataframe
hepg2.coldata = hepg2.summary_merged[,c("sampleID", "celltype", "treatment", 'rep')]
hepg2.coldata = hepg2.coldata %>% mutate( foxo1_inh = case_when( treatment == 'ctrl' ~ "control",
                                                                 TRUE ~ "FOXO1_inh"))
hepg2.coldata$foxo1_inh = factor(hepg2.coldata$foxo1_inh, levels=c('control', 'FOXO1_inh'))

hepg2.counts_merged =hepg2.counts_merged %>% column_to_rownames("Geneid")
hepg2.coldata = hepg2.coldata %>% column_to_rownames("sampleID")

hepg2.coldata = hepg2.coldata %>% mutate('batch' = ifelse( rownames(.) %in% paste0("HepG2_", c('ctrl_1', 'ctrl_2', "low_1", "high_1")), 1,2))
hepg2.coldata$treatment = factor(hepg2.coldata$treatment, levels=c('ctrl', 'low', 'high'))
hepg2.coldata$batch = factor(hepg2.coldata$batch)

write.table(hepg2.coldata %>% rownames_to_column('sample_name'), paste0(hepg2.analysis_dir, 'analysis/hepg2_sample_metadata.txt'), sep='\t', row.names = F)

hepg2.colors=c("grey", 'lightskyblue2', 'lightskyblue4')

# DESeq normalize
hepg2.dds <- DESeqDataSetFromMatrix(countData = hepg2.counts_merged,
                                    colData=hepg2.coldata,
                                    design = ~ factor(foxo1_inh))
hepg2.dds <- estimateSizeFactors(hepg2.dds)
hepg2.dds.normCounts <- counts(hepg2.dds, normalized=T)

# DEG
hepg2.dds <- hepg2.dds[hepg2.dds.features.selected, ]
design(hepg2.dds) = ~ batch + foxo1_inh

hepg2.dds <- DESeq(hepg2.dds)
hepg2.deg.res = results(hepg2.dds, contrast = c('foxo1_inh', "FOXO1_inh", 'control'))
hepg2.deg.res = data.frame(hepg2.deg.res)
hepg2.deg.res %>% rownames_to_column('gene_name') %>% fwrite(paste0(hepg2.analysis_dir, 'analysis/DEG_batch_foxo1_inh.txt'), sep='\t')



##############################################################
# ctrl vs 100nM
##############################################################
hepg2.dds.ctrl.low <- hepg2.dds[hepg2.dds.features.selected, (hepg2.coldata %>% filter(treatment %in% c('ctrl', 'low'))) %>% rownames ]
design(hepg2.dds.ctrl.low ) = ~ batch + foxo1_inh

hepg2.dds.ctrl.low  <- DESeq(hepg2.dds.ctrl.low )
hepg2.deg.res.ctrl.low  = results(hepg2.dds.ctrl.low , contrast = c('foxo1_inh', "FOXO1_inh", 'control'))
hepg2.deg.res.ctrl.low  = data.frame(hepg2.deg.res.ctrl.low )
hepg2.deg.res.ctrl.low  %>% rownames_to_column("gene_name") %>% fwrite(., paste0(hepg2.analysis_dir, 'analysis/DEG_batch_ctrl_vs_low.txt'), sep='\t')

##############################################################
# ctrl vs 1000nM
##############################################################
hepg2.dds.ctrl.high <- hepg2.dds[hepg2.dds.features.selected, (hepg2.coldata %>% filter(treatment %in% c('ctrl', 'high'))) %>% rownames ]
design(hepg2.dds.ctrl.high ) = ~ batch + foxo1_inh

hepg2.dds.ctrl.high  <- DESeq(hepg2.dds.ctrl.high )
hepg2.deg.res.ctrl.high  = results(hepg2.dds.ctrl.high , contrast = c('foxo1_inh', "FOXO1_inh", 'control'))
hepg2.deg.res.ctrl.high  = data.frame(hepg2.deg.res.ctrl.high )
hepg2.deg.res.ctrl.high  %>% rownames_to_column("gene_name") %>% fwrite(paste0(hepg2.analysis_dir, 'analysis/DEG_batch_ctrl_vs_high.txt'), sep='\t')



##############################################################
# Plot EFHD1 expression 
##############################################################
tmp = counts(hepg2.dds ,normalize=T)[c("EFHD1"),] %>% t() %>% as.data.frame() 
tmp = merge(tmp, hepg2.coldata, by=0, all=T)
tmp = tmp %>% mutate(treatment_group = case_when(treatment=='ctrl' ~ 1,
                                                 treatment=='low' ~ 2,
                                                 treatment=='high' ~ 3,
                                                 TRUE ~ NA))
tmp = tmp %>% mutate(drug_concentration = case_when(treatment=='ctrl' ~ "0",
                                                    treatment=='low' ~ "100nM",
                                                    treatment=='high' ~ "1000nM",
                                                    TRUE ~ NA))

ggplot(tmp, aes(y=EFHD1)) + geom_smooth(method='lm', alpha=0.3, color='grey80', fill="grey90", aes(x=treatment_group)) + 
  geom_jitter(width=0.2, height=0, aes(x=treatment, fill=drug_concentration), size=3, shape=21) + 
  scale_fill_manual(values=hepg2.colors) +theme_classic() + xlab("AS1842856") + scale_x_discrete(labels=c("0", "100nM","1000nM"))+ plot_layout(ncol=1)
ggsave(paste0(hepg2.analysis_dir, "analysis/FOXO1_EFHD1_expression.pdf"), width=4, height=2)







##################
##################
# FOXO1 inh in Organoid
##################
##################
org.analysis_dir="./organoid_foxo1_inh/"
org.counts.dir = paste0(org.analysis_dir, "counts/") 

# merge counts 
i=1
for(counts_cur in list.files(org.counts.dir, '*.counts$') ){
  counts_tmp = fread(paste0(org.counts.dir, counts_cur))[,c(1,7)]
  colnames(counts_tmp)[2] = gsub("Aligned.out.bam", "", colnames(counts_tmp)[2])
  colnames(counts_tmp)[2] = gsub("../bam/", "", colnames(counts_tmp)[2])
  
  if(i==1) { org.counts_merged = counts_tmp}
  else{ org.counts_merged = merge(org.counts_merged, counts_tmp, by='Geneid')}
  i = i+1
}

fwrite(org.counts_merged,paste0(org.analysis_dir, "analysis/organoid_raw_counts_merged.txt"),sep='\t' )

# remove genes without gene names (starts with ENSG)
org.counts_merged = org.counts_merged %>% column_to_rownames("Geneid")
org.appropriate_gene_names = rownames(org.counts_merged)[!grepl("^ENSG", rownames(org.counts_merged))]
org.counts_merged = org.counts_merged[org.appropriate_gene_names,]


# make sample metadata dataframe
org.coldata = org.summary_merged[,c("sampleID", "donor", "treatment", 'rep')]
org.coldata$treatment = factor(org.coldata$treatment, levels=c('ctrl', 'drug'))

org.counts_merged =org.counts_merged %>% column_to_rownames("Geneid")
org.coldata = org.coldata %>% column_to_rownames("sampleID")
org.coldata$donor = factor(org.coldata$donor)
org.coldata$treatment = factor(org.coldata$treatment)
org.coldata$rep = factor(org.coldata$rep)


# DESeq normalize
org.dds <- DESeqDataSetFromMatrix(countData = org.counts_merged,
                                  colData=org.coldata,
                                  design = ~ factor(treatment))
org.dds <- estimateSizeFactors(org.dds)
org.dds.normCounts <- counts(org.dds, normalized=T)


# DEG
org.dds <- org.dds[org.dds.features.selected, ]
design(org.dds) = ~ donor + treatment

org.dds <- DESeq(org.dds)
org.deg.res = results(org.dds, contrast = c('treatment', "drug", 'ctrl'))
org.deg.res = data.frame(org.deg.res)
org.deg.res %>% rownames_to_column("gene_name") %>% filter(padj < 0.05) %>% fwrite(paste0(org.analysis_dir, 'analysis/DEG_batch_foxo1_inh.txt'), sep='\t')


############################################################
### M12 vs FOXO1-inh-DEGs relationship: kME ~ DEG Log2FC ###
############################################################
tmp = modules_info %>% remove_rownames() %>% 
  merge(., org.deg.res %>% rownames_to_column("gene_name"), by='gene_name', all=F)

tmp$abs_log2FoldChange = abs(tmp$log2FoldChange)

for(module_cur in c(paste0("Hep-M", c(3,4,6,8,12,14)),'grey')){
  #if(module_cur == 'grey') {cortest_formula = as.formula( paste0("kME_",module_cur, "~ log2FoldChange" )) }
  cortest_formula = as.formula( paste0("~`kME_",module_cur, "`+ abs_log2FoldChange" )) 
  cortest_res_tmp = cor.test(formula=cortest_formula, data = tmp %>% filter(module == module_cur), method='spearman')
  cortest_res_tmp = data.frame(cortest_res_tmp$estimate[[1]], cortest_res_tmp$p.value, module=module_cur) 
  if(module_cur=='Hep-M3'){cortest_res = cortest_res_tmp}
  else{ cortest_res = rbind(cortest_res, cortest_res_tmp)}
}

colnames(cortest_res) = c("R", "P", 'module')

ggplot(cortest_res, aes(x=reorder(module, -R, mean), y=R, size=(-log10(P)), color=module)) + geom_point() + ylab( "R (kME vs -log2FC)") + theme_classic() + xlab("module") +
  scale_color_manual(limits=c("Hep-M12", 'grey', 'Hep-M3', 'Hep-M4', "Hep-M6", "Hep-M8", "Hep-M14"),
                     values=c("darkorchid",'grey', "#0066CC", "#007FFF", "#3399FF", "#6699CC", "#3366CC") ) + 
  geom_hline(yintercept=0, linetype='dashed') + scale_size(range=c(1,8)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + ylim(c(-0.15,0.35))
ggsave(paste0(org.analysis_dir, "analysis/kME_log2FC_module_genes_spearman_cor.pdf"), width=3.7, height=2.8)



for(module_cur in c(paste0("Hep-M", c(3,4,6,8,12,14)))){
  #if(module_cur == 'grey') {cortest_formula = as.formula( paste0("kME_",module_cur, "~ log2FoldChange" )) }
  egenes_cur =  (pme.sig.anno %>% filter(!is.na(!!sym( module_cur))))$gene %>% unique
  cortest_formula = as.formula( paste0("~`kME_",module_cur, "`+ abs_log2FoldChange" )) 
  cortest_res_tmp = cor.test(formula=cortest_formula, data = tmp %>% filter(gene_name %in% egenes_cur), method='spearman')
  cortest_res_tmp = data.frame(cortest_res_tmp$estimate[[1]], cortest_res_tmp$p.value, module=module_cur) 
  if(module_cur=='Hep-M3'){cortest_res_i = cortest_res_tmp}
  else{ cortest_res_i = rbind(cortest_res_i, cortest_res_tmp)}
}

colnames(cortest_res_i) = c("R", "P", 'module')

ggplot(cortest_res_i, aes(x=reorder(module, -R, mean), y=R, size=(-log10(P)), color=module)) + geom_point() + ylab( "R (kME vs -log2FC)") + theme_classic() + xlab("module") +
  scale_color_manual(limits=c("Hep-M12", 'grey', 'Hep-M3', 'Hep-M4', "Hep-M6", "Hep-M8", "Hep-M14"),
                     values=c("darkorchid",'grey', "#0066CC", "#007FFF", "#3399FF", "#6699CC", "#3366CC") ) + 
  geom_hline(yintercept=0, linetype='dashed') + scale_size(range=c(1,8)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + ylim(c(-0.15,0.2))
ggsave(paste0(org.analysis_dir, "analysis/kME_abs_log2FC_module_igenes_spearman_cor.pdf"), width=3.7, height=2.8)



############################################################
### M12 vs FOXO1-inh-DEGs relationship: DEG list overlap ###
############################################################
# fisher-test
all.degs= org.deg.res %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 1) %>% rownames()
all.genes=intersect(org.dds.features.selected, modules_info$gene_name )

for(module_cur in c(wgcna.module.count[["hepatocyte"]], 'grey')){
  if(module_cur != 'grey') { module_cur = paste0("Hep-M", module_cur) }
  
  egenes_cur =  (pme.sig.anno %>% filter(!is.na(!!sym( module_cur))))$gene %>% unique
  genes_cur =( modules_info %>% filter(module == module_cur))$gene_name
  
  fisher_res_tmp_egene = fisher_test_genes(egenes_cur, all.degs, all.genes) %>% 
    mutate('module'=module_cur, 'category' = "iegene", 'gene_count' = length(egenes_cur))
  fisher_res_tmp_gene = fisher_test_genes(genes_cur, all.degs, all.genes) %>% 
    mutate('module'=module_cur, 'category' = "module_gene", 'gene_count' = length(genes_cur))
  fisher_res_tmp = rbind(fisher_res_tmp_egene, fisher_res_tmp_gene)
  
  if(module_cur == 'grey') { fisher_res_tmp = fisher_res_tmp_gene}
  
  if(module_cur == paste0("Hep-M", wgcna.module.count[["hepatocyte"]][1])){ fisher_res_m = fisher_res_tmp}
  else{fisher_res_m = rbind(fisher_res_m, fisher_res_tmp)}
}
fisher_res_m %>% filter(category == 'iegene') %>% arrange(-OR)
fisher_res_m %>% filter(category == 'module_gene') %>% arrange(-OR)

fisher_res_m %>% remove_rownames() %>% fwrite(paste0(org.analysis_dir, 'analysis/org.deg.abs_log2FC_1.fisher_res_degs_n5281.txt'), sep='\t')


# plot fisher res : all DEG - order by OR
fisher_res_cur = fread(paste0(org.analysis_dir, 'analysis/org.deg.abs_log2FC_1.fisher_res_degs_n5281.txt'), sep='\t')
fisher_res_cur$module = factor(fisher_res_cur$module, levels=c("Hep-M12", 'grey', 'Hep-M3', 'Hep-M4', "Hep-M6", "Hep-M8", "Hep-M14"))


ggplot(fisher_res_cur %>% filter(category == 'iegene'), aes(x=reorder(module, -OR, FUN=mean), y=OR, color=module)) + geom_point(position=position_dodge(.9), aes(size=-log10(p_value))) + 
  geom_errorbar(aes(ymin=conf_int_low, ymax=conf_int_hi), width=.2,position=position_dodge(.9)) + 
  geom_hline(yintercept=1, linetype='dashed') + theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + 
  scale_color_manual(limits=c("Hep-M12", 'grey', 'Hep-M3', 'Hep-M4', "Hep-M6", "Hep-M8", "Hep-M14"),
                     values=c("darkorchid",'grey', "#0066CC", "#007FFF", "#3399FF", "#6699CC", "#3366CC") ) + ggtitle("DEG enrichment (ieGenes)") + xlab("Module") + 
  ggplot(fisher_res_cur %>% filter(category == 'module_gene'), aes(x=reorder(module, -OR, FUN=mean), y=OR, color=module)) + geom_point(position=position_dodge(.9), aes(size=-log10(p_value))) + 
  geom_errorbar(aes(ymin=conf_int_low, ymax=conf_int_hi), width=.2,position=position_dodge(.9)) + 
  geom_hline(yintercept=1, linetype='dashed') + theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + xlab("Module") +
  scale_color_manual(limits=c("Hep-M12", 'grey', 'Hep-M3', 'Hep-M4', "Hep-M6", "Hep-M8", "Hep-M14"),
                     values=c("darkorchid",'grey', "#0066CC", "#007FFF", "#3399FF", "#6699CC", "#3366CC") ) + ggtitle("DEG enrichment (module genes)") 
ggsave(paste0(org.analysis_dir, 'analysis/org.deg.abs_log2FC_1.fisher_res_degs_n5281.pdf'), width=8.5, height=3.5)

