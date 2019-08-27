# prepare database for shiny app

library(tidyverse)
library(GenomicFeatures)
library(Gviz)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#------------
# dataset
#------------
dataset <- read.table('data/master_intergration_table_circRNA__bycirc_FACS__commonSNP_and_distalSNP_with_COLOC__DEG_biopsy_overlap_with_miRNA_binding_with_circbase.tab', header=T, check.names = F)

cid <- c('geneID', 'circID', "Exon1",
         "Exon2", "Strand", 'SpliceType',
         'AveExpr_log2cpm',  'AveExpr_log2cpm_Biopsy',
         'mean_circ_fraction', 'mean_circ_fraction_Biopsy',
         # DE information
         'logFC_CD_vs_Control', 'padj_CD_vs_Control',
         'logFC_UC_vs_Control','padj_UC_vs_Control',
         'logFC_vdjnorm', 'padj_vdjnorm',
         #!miRNA binding information
         'best_binding_miRNA', 'miRNA_bind_sites_per_kB_circRNA',
         #celltype enet model information
         'predictor', 'no_parent_regressed',
         'parent_regressed', 'delta_parent',
         'no_parent_regressed_rank', 'parent_regressed_rank',
         #circQTL
         'num_feature__cRNA', 'num_SNP_feature__cRNA',
         'num_SNP_feature__mRNA' , "middle_beta__cRNA",
         #eQTL parent
         'middle_beta__mRNA', 'cis_trans__mRNA',
         'cis_trans__cRNA','common_or_uncommon', 'opposite_QTL_effect',
         #coloc
         'trait', 'mRNA_mlog10_pval', 'cRNA_mlog10_pval',
         'method','is_cRNA_coloc_greater',
         # circDB
         'circbase_ID', 'genomic_length', 'spliced_seq_length',
         'repeats', 'annotation', 'circRNA_study')

# subset relevant columns for table display
dataset <- dataset[, cid]


# cid <- c('GeneSymbol', 'BackspliceLocation', "BackspliceExon1",
#          "BackspliceExon2", "Strand", 'SpliceType',
#          'AveExpr_log2cpm',  'AveExpr_log2cpm_Biopsy',
#          'MeanCircFraction', 'MeanCircFraction_Biopsy',
#          # DE information
#          'logFC_CD_vs_Control', 'padj_CD_vs_Control',
#          'logFC_UC_vs_Control','padj_UC_vs_Control',
#          'logFC_vdjnorm', 'padj_vdjnorm',
#          #!miRNA binding information
#          'best_binding_miRNA', 'miRNA_bind_sites_per_kB_circRNA',
#          #celltype enet model information
#          'celltype_enet_predictor', 'celltype_enet_no_parent_regressed',
#          'celltype_enet_parent_regressed', 'celltype_enet_delta_parent',
#          'celltype_enet_no_parent_regressed_rank', 'celltype_enet_parent_regressed_rank',
#          #circQTL
#          'circQTL_num_feature_cRNA', 'circQTL_num_SNP_feature_cRNA',
#          'circQTL_num_SNP_feature_mRNA' , "circQTL_middle_beta_cRNA",
#          #eQTL parent
#          'parent_eQTL_middle_beta_mRNA', 'parent_eQTL_cis_trans_mRNA',
#          'parent_eQTL_cis_trans_cRNA','parent_eQTL_common_or_uncommon', 'parent_eQTL_opposite_QTL_effect',
#          #coloc
#          'trait', 'mRNA_mlog10_pval', 'cRNA_mlog10_pval',
#          'method','is_cRNA_coloc_greater',
#          # circbase
#          'circbase_ID', 'genomic_length', 'spliced_seq_length',
#          'repeats', 'annotation', 'circRNA_study')


chars <- c('geneID', 'circID', 'predictor')
dataset[, chars] <- lapply(dataset[,chars], as.character)

# colnames(dataset) <- cid


# separate out numeric, factor, and string columns for their appropriate selections
nums <- sapply(dataset, is.numeric)
facs <- sapply(dataset, is.factor)
chars <- sapply(dataset, is.character)

exprs_cols <- c('logFC_CD_vs_Control', 'logFC_UC_vs_Control', 'logFC_vdjnorm')
pval_cols <- c('padj_CD_vs_Control', 'padj_UC_vs_Control', 'padj_vdjnorm')

# prep df for volcano plot
vdata <- dataset %>%
    mutate(padj_CD_vs_Control = -log(padj_CD_vs_Control),
           padj_UC_vs_Control = -log(padj_UC_vs_Control),
           padj_vdjnorm = -log(padj_vdjnorm))

ensids <- mapIds(org.Hs.eg.db, keys = as.character(dataset$geneID), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
entrezids <- mapIds(org.Hs.eg.db, keys = as.character(dataset$geneID), column="ENTREZID", keytype="SYMBOL", multiVals="first")
dataset$ensembl_id <- ensids
dataset$entrez_id <- entrezids

# define genomic browser tracks & info
gbrowsedf <- dataset %>% tidyr::separate(col = circID, into=c('chr', 'start', 'end'),'-|:', remove=FALSE)
gbrowsedf$gr_strand <- as.character(gbrowsedf$Strand)
gbrowsedf$gr_strand[gbrowsedf$gr_strand %in% c('-/+', '+/-')] <- '*'

axis_track <- GenomeAxisTrack()


browser_gr <- makeGRangesFromDataFrame(gbrowsedf, seqnames.field = 'chr', start.field = 'start', end.field = 'end', strand.field = 'gr_strand',keep.extra.columns = TRUE)

snps <- read_tsv('data/common_cRNA_mRNA_SNPS_by_gene.txt')


col_details <- read_csv('data/master_data_header.txt', col_names = c('colname', 'detail'))

sketch <- htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th('GeneID', title='GeneSymbol: geneic locus of transcription'),
      th('circID', title='Backsplice: start/end of backsplice mapping in hg19'),
      th('Exon1', title='Exon_1: annotated first exon within major isoform'),
      th('Exon2', title='Exon_2: annotated second exon within major isoform'),
      th('Strand', title='Strand: backsplice alignment strand'),
      th('SpliceType', title='SpliceType: exonic or intronic circRNA'),
      th('AveExpr_log2cpm', title='AverageExpression: log2cpm average MSCCR-blood expression'),
      th('AveExpr_log2cpm_Biopsy', title='AverageExpression: log2cpm average MSCCR-biopsy expression (if observed)'),
      th('mean_circ_fraction', title='MeanCircularization: backsplice support normalized by maximal forward splice (MSCCR-blood)'),
      th('mean_circ_fraction_Biopsy', title='MeanCircularization: backsplice support normalized by maximal forward splice (MSCCR-biopsy)'),
      th('logFC_CD_vs_Control', title='logFoldChange_CD_vs_Control: log fold change for CD vs Control (MSCCR-blood)'),
      th('padj_CD_vs_Control', title='p_adj_CD_vs_Control: BH-corrected p value for CD vs Control (MSCCR-blood)'),
      th('logFC_UC_vs_Control', title='logFoldChange_UC_vs_Control: log fold change for UC vs Control (MSCCR-blood)'),
      th('padj_UC_vs_Control', title='p_adj_UC_vs_Control: BH-corrected p value for UC vs Control (MSCCR-blood)'),
      th('logFC_vdjnorm', title='logFoldChange_VDJ_seq: log fold change for T/BCR seq normalized total expression (MSCCR-blood)'),
      th('padj_vdjnorm', title='p_adj_VDJ_seq: BH-corrected p value for for T/BCR seq normalized total expression (MSCCR-blood)'),
      th('best_binding_miRNA', title='Predicted_miRNA: best-3prime UTR binding miRNA'),
      th('miRNA_bind_sites_per_kB_circRNA', title='miRNA_Sponginess: number of miRNA binding sites per kB spliced circRNA'),
      th('predictor', title='CellType: FACS category used in elastic-net model for circRNA expression'),
      th('no_parent_regressed', title='CellType_enet_Coefficient_without_parent: value of elastic-net coefficient for predictor without parent mRNA regressed'),
      th('parent_regressed', title='CellType_enet_Coefficient_with_parent: value of elastic-net coefficient for predictor with parent mRNA regressed (<=CellType_enet_Coefficient)'),
      th('delta_parent', title='CellType_Parent_Effect: difference in elastic-net coefficients when regressing out parent'),
      th('no_parent_regressed_rank', title='CellType_enet_model_rank: term rank in model (lowest is more important; highest least important) without regressing parent mRNA'),
      th('parent_regressed_rank', title='CellType_enet_model_rank: term rank in model (lowest is more important; highest least important) with regressing parent mRNA'),
      th('num_feature__cRNA', title='Num_of_circRNA_isoforms: number of distinct circRNA transcribed from same locus controlled by SNP'),
      th('num_SNP_feature__cRNA', title='Num_SNPs: number of SNPs associated with circRNA'),
      th('num_SNP_feature__mRNA', title='Num_SNPs: number of SNPs associated with mRNA'),
      th('middle_beta__cRNA', title='AverageBeta_circRNA: average logFoldChange in circRNA for every copy of variant allele'),
      th('middle_beta__mRNA', title='AverageBeta_mRNA: average logFoldChange in parent mRNA for every copy of variant allele'),
      th('cis_trans__mRNA', title='Cis_Trans_mRNA: cis or trans SNPs for mRNA'),
      th('cis_trans__cRNA', title='Cis_Trans_mRNA: cis or trans SNPs for circRNA'),
      th('common_or_uncommon', title='Opposite_QTL_Effects: do AverageBeta_circRNA and AverageBeta_mRNA have same sign (common)'),
      th('opposite_QTL_effect', title='Opposite_QTL_Effects: do AverageBeta_circRNA and AverageBeta_mRNA have same sign (no)'),
      th('trait', title='coloc_Trait: phenotypic trait tested in de Lange GWAS (PMC5289481): IBD | CD | UC'),
      th('mRNA_mlog10_pval', title='mRNA_parent_cauaslity_pvalue: -log10(p_value) of causality test for risk mediation (pleiotropy null); SNP -> mRNA -> trait'),
      th('cRNA_mlog10_pval', title='circRNA_parent_cauaslity_pvalue: -log10(p_value) of causality test for risk mediation (pleiotropy null); SNP -> circRNA -> trait'),
      th('method', title='coloc_method: colcalization method used (SMR | MetaXCan | COLOC)'),
      th('is_cRNA_coloc_greater', title='circRNA_colocalize: do circRNA daughers mediate more risk for trait than mRNA parents'),
      th('circbase_ID', title='circbase_ID: circbase ID for this backsplice'),
      th('genomic_length', title='GenomicLength: length of unspliced circRNA'),
      th('spliced_seq_length', title='SplicedLength: length of spliced circRNA'),
      th('repeats', title='Repeats: overlapping known repeats'),
      th('annotation', title='Annotation: overlapping known hg19 annotation'),
      #th('in_other_samples', title='InOtherSamples: also found in these celltypes'),
      th('circRNA_study', title='InOther_circRNA_study: also observed in these studies'),
      th('ensembl_id', title='ensembl ID of the linear parent'),
      th('entrez_id', title='entrez ID of the linear parent')
    )
  )
))

save(nums, facs, chars, exprs_cols, pval_cols, dataset, gbrowsedf,browser_gr, vdata, axis_track, snps,sketch,col_details, file='data/_circdb.rdata')


