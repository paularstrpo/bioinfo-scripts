getMatchedPheno <- function(project){
    # This function will get phenotype data and sample/file info for data download
    # for rna-seq counts information for those TCGA samples that have available matched
    # primary tumor and normal tissue data.

    # The function will save an RData file containing clinical, manifest, and
    # case id list information.

    # Dependencies include TCGAbiolinks which can be installed via dev_tools::install_github
    # > devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")

    library(TCGAbiolinks)
    query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",
                    sample.type = c("Primary solid Tumor", "Solid Tissue Normal"))

    manifest.table <- data.frame(query$results[[1]])
    manifest.table$case.id <- substr(manifest.table$cases,1,12)
    manifest.table$sample.id <- substr(manifest.table$cases,1,16)
    manifest.table$sample.type <- substr(manifest.table$cases, 14,16)
    #manifest.table <- manifest.table[! duplicated(manifest.table$sample.id), ]

    # find samples which have both recurrent and primary available
    primary <- manifest.table[(manifest.table$tissue.definition == "Primary solid Tumor") & (manifest.table$sample.type == '01A'), ]

    normal <- manifest.table[(manifest.table$tissue.definition == "Solid Tissue Normal") & (manifest.table$sample.type == '11A'), ]
    nrow(normal)


    matched.manifest <- rbind(primary, normal)

    case.ids <- merge(primary, normal, by=c('case.id'))
    case.ids <- levels(factor(case.ids$case.id))
    length(case.ids)

    manifest.table <- matched.manifest[which(matched.manifest$case.id %in% case.ids), ]
    nrow(manifest.table)

    clinical <- GDCquery_clinic(project = project, type = "clinical")
    clinical <- clinical[which(clinical$submitter_id %in% levels(factor(normal$case.id))) | which(clinical$submitter_id %in% levels(factor(primary$case.id))),]
    nrow(clinical)

    results <- list(clinical, manifest.table, normal, primary, project)
    names(results) <- c('clinical_data', 'manifest_table', "normal", 'primary_tumor', 'project')
    return(results)
}

TCGAtranslateID <- function(file_ids, legacy = TRUE) {
    library(GenomicDataCommons)
    library(tidyverse)
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