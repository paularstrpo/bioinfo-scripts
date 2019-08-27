# -------------------------------- #
# VCF parsing                      #
# author: Paula Restrepo           #
# email:  paularstrpo@gmail.com    #
# -------------------------------- #

# # # # # # # # # # #
# Set things up...  #
# # # # # # # # # # #
library(vcfR)
library(tidyverse)


extract.vaf <- function(var.table, AD.col='gt_AD') {
    var.table <- tidyr::separate(var.table, AD.col, c("AD_N", "AD_T"), sep=",")
    var.table[, c("AD_N", "AD_T")] <- lapply(var.table[, c("AD_N", "AD_T")], as.numeric)
    var.table$VAF_AD <- 100 * ( var.table$AD_T / (var.table$AD_N + var.table$AD_T) )
    return(var.table)
}

parse_extra <- function(df, annot){
    res <- substr(x=df$Extra, start=regexpr(pattern=paste(annot, '=', sep=''),
                                            text=df$Extra), stop=nchar(df$Extra) )
    res <- substr(x=res,  start=1, stop=regexpr(pattern=';', text=res))
    res <- gsub(x=gsub(pattern=paste(annot, '=', sep=''), replacement='', x=res),
                pattern=';', replacement='')
    res[!grepl(annot, df$Extra)] <- 'NA'
    return(res)
}



# # # # # # # # # # # #
# import & parse vcf  #
# # # # # # # # # # # #

# Note: This is originally written for use with annotated Mutect2 vcf files downloaded from TCGA.
#       Adjust for your purposes as necessary.

# if you have a sample sheet use this.
sample_sheet <- read_tsv('metadata/vcf_sample_sheet.tsv', 
                col_names=c('file_id', 'file_name', 'data_category', 'data_type',
                            'project_id', 'patient_id', 'sample_id', 'sample_type'),
                skip=1)
sample_sheet[, c('patient_id', 'sample_id', 'sample_type')] <- lapply(sample_sheet[, c('patient_id', 'sample_id', 'sample_type')], function(x){substr(x, 1, regexpr(',', x)-1)})

tcga_cesc_phenodata <- readRDS('metadata/cesc_parsed_clinical_phenodata.RDS')
sample_sheet <- sample_sheet %>% filter(sample_sheet$patient_id %in% tcga_cesc_phenodata$submitter_id)

path <- 'vcf/' 
# 1) import vcf files and use vcfR to read them in.
# use a sample sheet with the column containing the vcf file names of your samples so it knows which to read in.
raw_var<- lapply(paste0(path, sample_sheet$file_name), read.vcfR, verbose=FALSE, limit=1500)
raw_var <- lapply(raw_var, vcfR2tidy, single_frame=TRUE)
names(raw_var) <- sample_sheet$file_id

# comment these two lines out if no annotations
# this will grab the column name and format of the VEP annotations from the metadata part of the vcfR tidy object.
annot <- as.data.frame(raw_var[[1]]$meta) %>% filter(ID=='CSQ') %>% dplyr::select(Description) %>% deframe()
annot <- unlist(strsplit(gsub(pattern='Consequence annotations from Ensembl VEP. Format: ', replacement='', x=annot, fixed=TRUE), '|', fixed=TRUE))
#annot <- unlist(strsplit('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|ENTREZ|EVIDENCE', '|', fixed=TRUE))

# apply some basic filters. 
somatic_mutations <- bind_rows(lapply(names(raw_var), function(x){
       y <- raw_var[[x]]$dat %>%
              filter(FILTER == "PASS") %>% # keep variants that pass mutect2 filters
              filter(gt_AF > 0.05) %>% # keep only variants with at minimum 5% VAF. remove this if you feel it's too stringent.
               mutate(file_id=x,
                     file_name=sample_sheet$file_name[sample_sheet$file_id==x], # add sample specific info as you see fit.
                     patient_id=sample_sheet$patient_id[sample_sheet$file_id==x],
                     sample_id=sample_sheet$sample_id[sample_sheet$file_id==x],
                     sample_type=sample_sheet$sample_type[sample_sheet$file_id==x],
                     variant=paste0(CHROM, ":", POS, "_", REF, ">", ALT),
                     combo=paste0(sample_id, '_',variant)) %>%
                     tidyr::separate(col=CSQ, into=annot, sep='\\|') # comment this line and get rid of the %>% at the end of the above line if no annotation
       return(y)
}))
rm(annot, raw_var, path, sample_sheet)
# # # # SAVE DATA
saveRDS(somatic_mutations, file='metadata/tcga_somatic_mutations.RDS')

# # # # Gather some summary metrics (like mutation count, etcetera)
# # # #restart R session if too much memory footprint
somatic_mutations$variant_class <- somatic_mutations$Consequence
somatic_mutations$variant_class[grepl('stop_gained', x=somatic_mutations$variant_class)] <- 'nonsense'
somatic_mutations$variant_class[grepl('stop_lost', x=somatic_mutations$variant_class)] <- 'nonsense'
somatic_mutations$variant_class[grepl('start_lost', x=somatic_mutations$variant_class)] <- 'nonsense'
somatic_mutations$variant_class[grepl('NMD_transcript_variant', x=somatic_mutations$variant_class)] <- 'nonsense'

somatic_mutations$variant_class[grepl('frameshift', x=somatic_mutations$variant_class)] <- 'frameshift'
somatic_mutations$variant_class[grepl('deletion', x=somatic_mutations$variant_class)] <- 'frameshift'
somatic_mutations$variant_class[grepl('insertion', x=somatic_mutations$variant_class)] <- 'frameshift'

somatic_mutations$variant_class[grepl('splice', x=somatic_mutations$variant_class)] <- 'splice'
somatic_mutations$variant_class[grepl('missense', x=somatic_mutations$variant_class)] <- 'missense'
somatic_mutations$variant_class[!somatic_mutations$variant_class %in% c('frameshift', 'missense', 'splice')] <- 'silent'

somatic_mutations$variant_class <- as.factor(somatic_mutations$variant_class)
somatic_mutations$is_nonsilent <- somatic_mutations$variant_class!='silent'

mutation_summary<- somatic_mutations %>%
                        group_by(sample_id) %>%
                        summarise(total_mut_count=n_distinct(variant),
                        maxVAF=max(gt_AF),
                        minVAF=min(gt_AF),
                        medianVAF=median(gt_AF))

mutation_summary<- somatic_mutations %>% group_by(sample_id, is_nonsilent) %>%
 summarise(nonsilent_mutation_count=n_distinct(variant)) %>%
 filter(is_nonsilent==TRUE) %>%
 dplyr::select(-is_nonsilent) %>%
 inner_join(mutation_summary, by=c('sample_id'))
saveRDS(somatic_mutations, 'metadata/somatic_mutations.RDS')
saveRDS(mutation_summary, 'metadata/mutation_summary.RDS')
