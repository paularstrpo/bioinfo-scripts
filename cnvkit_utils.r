parse.cnvkit <- function(path, ext, sample.names) {
    library(tidyverse)
    flist <-  list.files(path=path, pattern = ext, full.names = TRUE, recursive = TRUE)
    cns_list <- lapply(flist, read_tsv)
    names(cns_list) <- sample.names
    cns_list <- lapply(names(cns_list), function(sample){
        x <- cns_list[[sample]]
        x$sample <- sample
        x$segment_length = x$end - x$start
        return(x)})

    return(cns_list)
}

plot.segment.distribution <- function(cns_table){
    library(tidyverse)
    plot <- ggplot(data=cns_table) +
    aes(x=segment_length, fill=sample) + 
    geom_density(alpha=0.2) + 
    facet_wrap(~sample) + 
    theme_classic()
    return(plot)
}

plot.singlechr.circlize <- function(cn_dat, mut_dat, chrom, ref_genome){
    library(tidyverse)
    library(circlize)

    cn_dat <- cn_dat %>% 
            select(chromosome, start, end, log2) %>%
            filter(chromosome==chrom)

    mut_dat <- mut_dat %>%
            mutate(POS2=POS) %>% 
            select(CHROM, POS, POS2, VAF, gene_symbol) %>%
            filter(CHROM == chrom)
            
    return({
        circos.par("track.height" = 0.12, start.degree = 180,
                canvas.xlim = c(-1, 1), canvas.ylim = c(0, 1), gap.degree = 180,
                cell.padding = c(0, 0, 0, 0))

        circos.initializeWithIdeogram(species = ref_genome, chromosome.index = chrom)
        circos.genomicLabels(mut_dat, labels.column = 4, side = "outside", cex=0.5)

        circos.genomicTrack(cn_dat, ylim = c(-2, max(cn_dat$log2)+1), panel.fun = function(region, log2, ...) {
            circos.genomicLines(region, 0, col='black', pch = 16, cex = 0.5, ...)
            circos.genomicLines(region, log2, col="#f7a048", pch = 16, cex = 0.5, ...)
            circos.genomicPoints(region, log2, col="#f7a048", pch = 16, cex = 0.5, ...)
        })

        circos.clear()
    })
}


plot.wholesample.circlize <- function(cn_dat, mut_dat, ref_genome){
        library(tidyverse)
        library(circlize)

        cn_dat <- cn_dat %>% 
                select(chromosome, start, end, log2)

        mut_dat <- mut_dat %>%
                mutate(POS2=POS) %>% 
                select(CHROM, POS, POS2, VAF, gene_symbol)

        return({
            circos.initializeWithIdeogram(species = ref_genome, chromosome.index = paste0("chr", c(1:22, "X", "Y")))
            circos.par("track.height" = 0.12)
            circos.genomicLabels(mut_dat, labels.column = 4, side = "outside", cex=0.5)

            circos.genomicTrack(cn_dat, ylim = c(-2, max(cn_dat$log2)+1), panel.fun = function(region, log2, ...) {
                circos.genomicLines(region, 0, col='black', pch = 16, cex = 0.5, ...)
                circos.genomicLines(region, log2, col="#f7a048", pch = 16, cex = 0.5, ...)
                circos.genomicPoints(region, log2, col="#f7a048", pch = 16, cex = 0.5, ...)
            })

            circos.clear()
        })
}
