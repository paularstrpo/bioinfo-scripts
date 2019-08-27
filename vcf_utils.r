extract.vaf <- function(var.table, AD.col='gt_AD') {
    library(tidyverse)
    var.table <- tidyr::separate(var.table, AD.col, c("AD_N", "AD_T"), sep=",")
    var.table[, c("AD_N", "AD_T")] <- lapply(var.table[, c("AD_N", "AD_T")], as.numeric)
    var.table$VAF_AD <- 100 * ( var.table$AD_T / (var.table$AD_N + var.table$AD_T) )
    return(var.table)
}

parse.extra <- function(df, annot){
    library(tidyverse)
    res <- substr(x=df$Extra, start=regexpr(pattern=paste(annot, '=', sep=''),
                                            text=df$Extra), stop=nchar(df$Extra) )
    res <- substr(x=res,  start=1, stop=regexpr(pattern=';', text=res))
    res <- gsub(x=gsub(pattern=paste(annot, '=', sep=''), replacement='', x=res),
                pattern=';', replacement='')
    res[!grepl(annot, df$Extra)] <- 'NA'
    return(res)
}

parse.somatic.vcf <- function(path, sample_sheet){
            library(tidyverse)
            library(vcfR)

            # read in vcf files using vcfR
            somatic_mutations<- lapply(paste0(path, sample_sheet$file_name), read.vcfR, verbose=FALSE, limit=1500)
            somatic_mutations <- lapply(somatic_mutations, vcfR2tidy, single_frame=TRUE)
            names(somatic_mutations) <- sample_sheet$sample_id

            # comment these two lines out if no annotations
            # this will grab the column name and format of the VEP annotations from the
            # metadata part of the vcfR tidy object.
            annot <- as.data.frame(somatic_mutations[[1]]$meta) %>% filter(ID=='CSQ') %>% dplyr::select(Description) %>% deframe()
            annot <- unlist(strsplit(gsub(pattern='Consequence annotations from Ensembl VEP. Format: ', replacement='', x=annot, fixed=TRUE), '|', fixed=TRUE))
            
            # apply some basic filters for VAF and based on Mutect2's internal checks.
            somatic_mutations <- bind_rows(lapply(names(somatic_mutations), function(x){
                y <- somatic_mutations[[x]]$dat %>%
                        filter(FILTER == "PASS") %>% # keep variants that pass mutect2 filters
                        filter(gt_AF > 0.05) %>% # keep only variants with at minimum 5% VAF. remove this if you feel it's too stringent.
                        mutate(sample_id=x,
                                variant=paste0(CHROM, ":", POS, "_", REF, ">", ALT),
                                combo=paste0(sample_id, '_',variant)) %>%
                                tidyr::separate(col=CSQ, into=annot, sep='\\|') # comment this line and get rid of the %>% at the end of the above line if no annotation
                return(y)
            }))
            
            # join back to sample sheet to append some basic phenotype information.
            somatic_mutations <- somatic_mutations %>% inner_join(sample_sheet, by=c('sample_id'))

            # define variant class in a more generalized way and identify which mutations are nonsilent.
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

            return(somatic_mutations)
}

summarise.mutations <- function(somatic_mutations){
    library(tidyverse)
    
    # takes the direct output of parse.somatic.vcf() to get some basic summary metrics.
    mutation_summary <- somatic_mutations %>%
                            group_by(sample_id) %>%
                            summarise(total_mut_count=n_distinct(variant),
                            maxVAF=max(gt_AF),
                            minVAF=min(gt_AF),
                            medianVAF=median(gt_AF))

    mutation_summary <- somatic_mutations %>% group_by(sample_id, is_nonsilent) %>%
                            summarise(nonsilent_mutation_count=n_distinct(variant)) %>%
                            filter(is_nonsilent==TRUE) %>%
                            dplyr::select(-is_nonsilent) %>%
                            inner_join(mutation_summary, by=c('sample_id'))

    return(mutation_summary)
}

parse.topiary.neoantigens <- function(path, sample_sheet){
    # your sample sheet should have a file_name column with the name of the topiary output
    # filenames matched to sample ID.
    
    library(tidyverse)
    epitope_list <- lapply(sample_sheet$file_name, read_csv)
    names(epitope_list) <- sample_sheet$sample_id

    epitope_df <- bind_rows(lapply(names(epitope_list), function(x){
        epitope.df <- epitope_list[[x]] %>%
        mutate(sample=x,
               variant=factor(gsub(x=gsub(x=variant, pattern="chr", replacement=''), pattern=' g.', replacement=':', fixed=TRUE)),
               combo=as.factor(paste(sample, variant, sep='_'))) ## this is what we merge on!!
    })) %>% inner_join(sample_sheet, by=c('sample_id'))
    
    return(epitope_df)
}