
waterfall_df <- somatic_mutations %>% as_tibble() %>%
  mutate(variant_class=ExonicFunc_refGene)

top_genes <- somatic_mutations %>% arrange(desc(as.numeric(gt_AF))) %>% group_by(sample) %>% top_n(n=25, wt=gt_AF)

gene_list <- c(gene_list1, top_genes$Gene_refGene, levels(factor(cosmic.v89.cgc$`Gene Symbol`)))
rm(top_genes)

waterfall_df[waterfall_df$variant_class == '.',]$variant_class <- waterfall_df[waterfall_df$variant_class == '.',]$Func_refGene

waterfall_df <- waterfall_df %>% dplyr::select(CHROM, POS, REF, ALT, Gene_refGene, variant_class, sample)
colnames(waterfall_df) <- c('chr', 'pos', 'ref', 'alt', 'gene', 'variant_class', 'sample')
variant_priority <- levels(factor(waterfall_df$variant_class))[c(2,3,8,9,12,11,10,13,1,14:16,4,5:7)]
waterfall_df$translational_effect <- 'non_silent'
waterfall_df[waterfall_df$variant_class %in% variant_priority[8:16], ]$translational_effect <- 'silent'

tmb_df <- waterfall_df %>% group_by(sample, translational_effect) %>%
  summarise(count=n(),tmb=count/54)
tmb_df$sample <-  fct_relevel(tmb_df$sample, c('Primary', "Recurrent_A", 'Recurrent_B', 'Recurrent_C'))


tmb2 <- waterfall_df %>% group_by(sample) %>%
  summarise(count=n(),tmb=count/54)
tmb2$sample <-  fct_relevel(tmb2$sample, c('Primary', "Recurrent_A", 'Recurrent_B', 'Recurrent_C'))

mut_heatmap_df <- waterfall_df[waterfall_df$gene %in% gene_list1,]
mut_heatmap_df$sample <-  fct_relevel(mut_heatmap_df$sample, c('Primary', "Recurrent_A", 'Recurrent_B', 'Recurrent_C'))

#' 
#' 
#' ## set up cnv chart
## ------------------------------------------------------------------------
cnv_heatmap_df <- cns_table %>%
  filter(abs(log2) >= 0.25)

cnv_heatmap_df <- lapply(gene_list, function(x){
  y <- cnv_heatmap_df[grepl(x, cnv_heatmap_df$gene, ignore.case = TRUE), ]
  if (nrow(y)!= 0) {
    y$query_gene <- x
  }
  return(y)
})
cnv_heatmap_df <- do.call(rbind, cnv_heatmap_df)
cnv_heatmap_df <- unique(cnv_heatmap_df)
cnv_heatmap_df$query_gene <- as.character(cnv_heatmap_df$query_gene)

cnv_heatmap_df <- cnv_heatmap_df %>% group_by(query_gene, sample) %>% summarise(average_cnv=mean(log2))
levels(cnv_heatmap_df$sample) <- c('Recurrent_A', 'Recurrent_B', 'Recurrent_C', 'Primary')

cnv_list <- c(gene_list1, 'SND1', "MIER2", "VAV1")
cnv_heatmap_df <- cnv_heatmap_df %>% filter(query_gene %in% cnv_list)
cnv_heatmap_df$cnv_type <- cut(cnv_heatmap_df$average_cnv, include.lowest = TRUE, breaks = c(-10,-1, -0.25, 0.25, 1, 10))
levels(cnv_heatmap_df$cnv_type) <- c('deletion', 'mild deletion', 'neutral', 'mild amplification', 'amplification')


cnv_heatmap_df$sample <- fct_relevel(cnv_heatmap_df$sample, c('Primary', "Recurrent_A", 'Recurrent_B', 'Recurrent_C'))

cnv_heatmap_df <- cnv_heatmap_df %>% filter(cnv_type != 'neutral')


#' Calculate the genomic instability index for each sample using 54 Mb with a 0.5-fold change threshold
## ------------------------------------------------------------------------
case_cnv_granges <- lapply(cns_list, function(x){makeGRangesFromDataFrame(x[abs(x$log2) > 0.5 ,], keep.extra.columns = TRUE, seqnames.field = 'chromosome', start.field = 'start', end.field = 'end')})
case_cnv_granges <- GRangesList(case_cnv_granges)
case_total_cnv_lengths <- sum(width(reduce(case_cnv_granges, ignore.strand=TRUE)))
cnv_lenths_mb <- case_total_cnv_lengths / (1e7)
case_genomic_instability <- as.data.frame(abs(case_total_cnv_lengths))
case_genomic_instability$cnv_lenths_mb <- cnv_lenths_mb
case_genomic_instability$instability_idx <- (cnv_lenths_mb / 54) * 100
rm(case_total_cnv_lengths, case_cnv_granges, cnv_lenths_mb)

pheno_summary <- case_genomic_instability

# as calculated from TheTA2 (presented as a proportion)
pheno_summary$tumor_purity <- c(36.6, # A
    37.6, #B
    30.2, #C
    38.2 ) #PRIMARY

pheno_summary$sample <- as.factor(c("Recurrent_A", 'Recurrent_B', 'Recurrent_C', 'Primary'))
pheno_summary$sample <- fct_relevel(pheno_summary$sample, c('Primary', "Recurrent_A", 'Recurrent_B', 'Recurrent_C'))

pheno_summary <- pheno_summary %>% tidyr::gather(var, val, c('tumor_purity', 'instability_idx'))
pheno_summary$var <- as.factor(pheno_summary$var)
levels(pheno_summary$var) <- c('Genomic Instability Index', 'Tumor Purity')
pheno_summary$val <- round(pheno_summary$val, digits=1)

#' ## plot the mutation summary
#' 
## ----fig.height=16, fig.width=10-----------------------------------------
mut_burden_barplot <- ggplot(data=tmb_df) + 
  aes(x=sample,y=tmb, fill=translational_effect) + 
  geom_bar(stat='identity') +
  geom_text(data=tmb2, aes(label=paste('N =',count), fill=NULL),size=8.5, nudge_y = 0.5, show.legend=FALSE) +
  labs(x='', y='TMB (Muts/MB)', fill='Translational Effect') +
  theme_classic() +
  theme +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
   scale_fill_manual(values = c('grey30', "grey70"))

mut_heatmap <- ggplot(data=mut_heatmap_df) +
  aes(x=sample, y=gene, fill=variant_class) +
  geom_tile(color='white', size=1, alpha=0.9) +
  labs(x='',y='SNV', fill='SNV Type') +
  theme_classic() +
  theme +
  scale_fill_manual(values = rev(brewer.pal(7, 'Spectral')) ) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

cnv_heatmap <- ggplot(data=cnv_heatmap_df) +
  aes(x=sample, y=query_gene, fill=cnv_type) +
  geom_tile(color='white', size=1) +
  labs(x='',y='CNV', fill='CNV Type') +
  theme_classic() +
  theme +
  scale_fill_manual(values = rev(brewer.pal(5,"RdBu")),drop=FALSE, labels = c("deletion", "", "neutral", "", "amplification"))

phenotype <- ggplot(data=pheno_summary) +
  aes(x=sample, y=var, fill=val) +
  geom_tile(color='white', size=1) +
  geom_text(aes(label=paste0(val, '%')), show.legend=FALSE, size=8.5) +
  labs(x='',y='', fill = "Percent") +
  theme_classic() +
  theme +theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_gradientn(colors = brewer.pal(100,"PuBu"), limits = c(0,60), breaks = c(0,30,60)) 


oncoprint_chart <- ggarrange(mut_burden_barplot ,phenotype,mut_heatmap,cnv_heatmap, nrow=4, align = 'v', heights = c(1.2,1,4,5))