# global.R

libs <- c('shiny',
          'tidyverse',
          'DT',
          'shinydashboard',
          'ggrepel',
          'GenomicFeatures',
          'Gviz',
          'GenomicRanges',
          'TxDb.Hsapiens.UCSC.hg19.knownGene',
          'circlize'
          )

# install.packages(libs)
lapply(libs, library, character.only=TRUE)

load('data/_circdb.rdata')