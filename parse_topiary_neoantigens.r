# # # # # # # # # # # # # # # # # # # #
# load and parse epitope predictions  #
# # # # # # # # # # # # # # # # # # # #
library(tidyverse)
library(ggrepel)
load('_rdata/parsed.variants.dnaseq.rdata')

## copy csv files from original output to gbm-ith-analysis/raw_data/topiary/.
# 1) read in epitope predictions

# DNA
path <-"path/to/your/topiary_outputs.csv"

flist <-  list.files(path=path, pattern = "*.csv$", full.names = TRUE)

epitope.list <- lapply(flist, function(x){
    x %>% read_csv() %>% as_tibble() %>%
        mutate(vcf_file=gsub(x=gsub(x=x, pattern=path, replacement=''),
                             pattern='.csv', replacement='', fixed=TRUE)) %>%
        mutate(sample=gsub(x=`vcf_file`, pattern='_filtered.vcf', replacement=''))
})

epitope.df <- do.call(rbind, epitope.list) %>%
    as_tibble() %>%
    mutate(sample=as.factor(sample))

levels(epitope.df$sample) <- sample.names ## clean up sample names before anything else!
epitope.df <- epitope.df %>%
    mutate(variant=factor(gsub(x=gsub(x=variant, pattern="chr", replacement=''),
                               pattern=' g.', replacement=':', fixed=TRUE)),
           combo=as.factor(paste(`sample`, `variant`, sep='_'))) ## this is what we merge on!!


# merge with existing variant info into a new table



# Need to adjust the combo columns in each df to take on the same format:
deletions <- grep('del',epitope.df$variant)
inserts <- grep('ins',epitope.df$variant)
others <- grep('>', epitope.df$variant)
epitope.df$combo_pos <- NA
epitope.df$combo_chr <- gsub(":.*","",epitope.df$variant)
epitope.df[deletions,]$combo_pos <- as.numeric(gsub("_.*","",gsub(".*:","",epitope.df[deletions,]$variant)))
epitope.df[inserts,]$combo_pos <- as.numeric(gsub("_.*","",gsub(".*:","",epitope.df[inserts,]$variant)))
epitope.df[others,]$combo_pos <- as.numeric(gsub("[A-Z]>.*","",gsub(".*:","",epitope.df[others,]$variant)))
epitope.df[deletions,]$combo_pos <- epitope.df[deletions,]$combo_pos-1
epitope.df$combo <- paste0(epitope.df$sample,"_",epitope.df$combo_chr,":",epitope.df$combo_pos)
vcf_filtered_df$combo <- paste0(vcf_filtered_df$sample,'_',gsub(x=paste0(vcf_filtered_df$CHROM,':',vcf_filtered_df$POS),  pattern="chr", replacement=''))

