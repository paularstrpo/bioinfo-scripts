
#----------------
# circlize
# ---------------
library(tidyverse)
library(circlize)
load('data/_circdb.rdata')

test_circle <- 'chr10:13233299-13234568'
snps <- read_tsv('data/common_cRNA_mRNA_SNPS_by_gene.txt')
snps_filtered <- snps %>%
    as_tibble() %>%
    dplyr::filter(trait_id == test_circle) %>%
    mutate(chr=paste0('chr', chr),
           pos=as.numeric(pos),
           pos2=as.numeric(pos))

snps_filtered <- as.data.frame(snps_filtered)
snps_filtered <- snps_filtered[, c(1,5,19, 2:4, 6:18)]


gene <- gbrowsedf[gbrowsedf$circID == test_circle,]$geneID[1]
chr <- gbrowsedf[gbrowsedf$circID == test_circle,]$chr[1]
junction <- browser_gr[browser_gr$circID == test_circle]




# get the exons with the gene coordinates
# to plot the backsplice
gquery <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)[which(genes(TxDb.Hsapiens.UCSC.hg19.knownGene)$gene_id == junction$entrez_id),]
tx <- subsetByOverlaps(transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene), gquery)
tx <- as.data.frame(tx)
exons <- subsetByOverlaps(exons(TxDb.Hsapiens.UCSC.hg19.knownGene), gquery)

exons <- as.data.frame(exons)
exons$gene <- as.character(gene)
exons$exon_id <- as.factor(exons$exon_id)
levels(exons$exon_id) <- as.character(1:length(levels(exons$exon_id)))



exons <- exons[, c('gene', 'start', 'end', 'seqnames', 'width', 'strand', 'exon_id')]

bs_exons <- exons[c(junction$Exon1,junction$Exon2),]
crna_exons <-  exons[c(junction$Exon1:junction$Exon2),]
bs_exon1 <- bs_exons[1, c(1:3, 5)]
bs_exon2 <- bs_exons[2, c(1:3, 5)]

circos.par(track.height = 0.1, gap.degree=180)
circos.genomicInitialize(exons)
circos.genomicTrack(snps_filtered, ylim=c(min(snps_filtered$beta) - 0.5, max(snps_filtered$beta) +0.5), panel.fun = function(region, beta, ...) {
    circos.genomicPoints(region, beta, col = 'blue', pch = 16, cex = 0.5, ...)
})

circos.genomicTrack(exons, ylim = c(-2.5, 2.5),
                    panel.fun = function(region, value, ...) {
                        # for each transcript
                        current_tx_start = min(region[, 1])
                        current_tx_end = max(region[, 2])
                        circos.lines(c(current_tx_start, current_tx_end),
                                     c(0, 0), col = "#CCCCCC")
                        circos.genomicRect(region, ytop = 1,
                                           ybottom = -1, col = "orange", border = NA)

                    }, bg.border = NA)


circos.genomicTrack(crna_exons, ylim = c(-2.5, 2.5),
                    panel.fun = function(region, value, ...) {
                        # for each transcript
                        current_tx_start = min(region[, 1])
                        current_tx_end = max(region[, 2])
                        circos.lines(c(current_tx_start, current_tx_end),
                                     c(0, 0), col = "#CCCCCC")
                        circos.genomicRect(region, ytop = 1,
                                           ybottom = -1, col = "red", border = NA)

                    }, bg.border = NA)

circos.genomicLink(bs_exon1, bs_exon2, col = 'red', border = NA, h=0.1, h2=0.1)

circos.clear()


