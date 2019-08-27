library(tidyverse)
#library(TCGAutils)

setwd('metadata')
x=read_tsv('fastq_manifest.txt')
file_uuids=x$id
library(GenomicDataCommons)
library(magrittr)

TCGAtranslateID = function(file_ids, legacy = TRUE) {
    info = files(legacy = legacy) %>%
        filter( ~ file_id %in% file_ids) %>%
        select('cases.samples.submitter_id') %>%
        results_all()
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list = lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file = sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list)))
    }

res = TCGAtranslateID(file_uuids)
head(res)

new_manifest <- inner_join(x, res, by=c("id"="file_id"))

# write_tsv(new_manifest, path='fastq_file_sampleID_mapping.tsv')

#   Alternatively, you can use a compact string representation
#   where each character represents one column: c = character, i
#   = integer, n = number, d = double, l = logical, f = factor, D
#   = date, T = date time, t = time, ? = guess, or ‘_’/‘-’ to
#   skip the column.

#cesc_phenodata <- read_tsv('cesc_phenodata.tsv', col_types=c("ccffnfnffnnnnffffnfnfnffnfffffnnffffffffnf"), na =c("--","NA", 'not reported'))
#rm(x, res, file_uuids)

new_manifest$sample_barcode <- new_manifest$submitter_id
new_manifest$submitter_id <- substr(new_manifest$submitter_id, 1,12)

colnames(new_manifest)[1] <- 'fastq_uuid'
colnames(new_manifest)[2:4] <- paste0('fastq_', colnames(new_manifest)[2:4])

library(TCGAutils)
library(curatedTCGAData)
barcode_defs <- TCGAbiospec(as.character(new_manifest$sample_barcode))

cols <- getcesc_phenodataNames("CESC")
test2<-colData(curatedTCGAData(diseaseCode = "CESC",assays = c("RNAseq*"), dry.run = FALSE))

cesc_phenodata <- as.data.frame(test2[, cols])
cesc_phenodata$submitter_id <- rownames(cesc_phenodata)
cesc_phenodata$full_sample_barcode <- toupper(test2$patient.samples.sample.portions.portion.analytes.analyte.2.aliquots.aliquot.bcr_aliquot_barcode)
cesc_phenodata$batch_number <- as.factor(substr(cesc_phenodata$full_sample_barcode, nchar(cesc_phenodata$full_sample_barcode) -6, nchar(cesc_phenodata$full_sample_barcode)))
# cesc_phenodata <- test3
saveRDS(cesc_phenodata, file='cesc_parsed_clinical_phenodata.RDS')