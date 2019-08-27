# takes mut_dat in the format outputted from the somatic_mutations dataframe from the vcf parser function.
oncoprinter <- function(mut_dat, cns_list, gene_list, total_seq_mb=54, cn_thresh=0.25, cn_levels=c(-2,-1, -0.25, 0.25, 1, 10), plot_heights=c(1.2,1,4,5)){
    library(GenomicRanges)
    library(tidyverse)
    library(ggpubr)

    # get TMB and VAF metrics
    summary_df <- mut_dat %>%
                        group_by(sample_id) %>%
                        summarise(total_mut_count=n_distinct(variant),
                        maxVAF=max(gt_AF),
                        minVAF=min(gt_AF),
                        medianVAF=median(gt_AF))

    summary_df <- mut_dat %>% group_by(sample_id, is_nonsilent) %>%
                        summarise(nonsilent_mutation_count=n_distinct(variant),
                                  tmb=nonsilent_mutation_count/total_seq_mb) %>%
                        filter(is_nonsilent==TRUE) %>%
                        dplyr::select(-is_nonsilent) %>%
                        inner_join(mut_summary, by=c('sample_id'))

    # select genes for the somatic mutation heatmap
    mut_dat <- mut_dat %>% filter(SYMBOL %in% gene_list)

    # select genes for cnv heatmap and prepare the dataframe
    cn_dat <- bind_rows(cns_list)
    cn_dat <- cn_dat %>%
    filter(abs(log2) >= cn_thresh)

    cn_dat <- unique(bind_rows(lapply(gene_list, function(x){
                        y <- cnv_heatmap_df[grepl(x, cn_dat$gene, ignore.case = TRUE), ]
                        if (nrow(y)!= 0) {
                            y$query_gene <- x
                        }
                        return(y)
                        })))

    cn_dat <- cn_dat %>%
              mutate(query_gene=as.character(query_gene)) %>%
              group_by(query_gene, sample) %>%
              summarise(average_cnv=mean(log2))

    cnv_dat$cnv_type <- cut(cnv_heatmap_df$average_cnv, include.lowest = TRUE, breaks = cn_levels)
    levels(cnv_dat$cnv_type) <- c('deletion', 'mild deletion', 'neutral', 'mild amplification', 'amplification')
    cnv_dat <- cnv_dat %>% filter(cnv_type != 'neutral')


    # calculate percentage of coding sequence altered by cnvs and append it to the summary table.
    cnv_granges <- lapply(cns_list, function(x){makeGRangesFromDataFrame(x[abs(x$log2) > cn_thresh ,], keep.extra.columns = TRUE, seqnames.field = 'chromosome', start.field = 'start', end.field = 'end')})
    cnv_granges <- GRangesList(case_cnv_granges)
    
    cnv_summary <- data.frame(
    total_cnv_lengths = sum(width(reduce(case_cnv_granges, ignore.strand=TRUE))),
    cnv_lenths_mb = (case_total_cnv_lengths / (1e7)),
    instability_idx = ((cnv_lenths_mb / 54) * 100)
    )
    cnv_summary$sample_id <- rownames(cnv_summary)

    summary_melt <- summary_df %>%
                  inner_join(cnv_summary, by=c('sample_id')) %>%
                  select(sample_id, maxVAF, medianVAF,instability_idx) %>%
                  tidyr::gather(var, val, c("maxVAF", "medianVAF","instability_idx")) 
    rm(cnv_granges, cnv_summary)

    ## begin plots
    
    # TMB barplot
    tmb_barplot <- ggplot(data=summary_df) + 
    aes(x=sample_id,y=tmb) + 
    geom_bar(stat='identity', fill='grey30') +
    geom_text(aes(label=paste('N =',nonsilent_mutation_count), fill=NULL),size=8.5, nudge_y = 0.5, show.legend=FALSE) +
    labs(x='', y='TMB (Muts/MB)', fill='Translational Effect') +
    theme_classic() +
    theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) 

    # Summary heatmap
    phenotype <- ggplot(data=summary_melt) +
    aes(x=sample, y=var, fill=val) +
    geom_tile(color='white', size=1) +
    geom_text(aes(label=paste0(val, '%')), show.legend=FALSE, size=8.5) +
    labs(x='',y='', fill = "Percent") +
    theme_classic() +
    theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

    # mutation heatmap
    mut_heatmap <- ggplot(data=mut_dat) +
    aes(x=sample_id, y=SYMBOL, fill=variant_class) +
    geom_tile(color='white', size=1, alpha=0.9) +
    labs(x='',y='Somatic Mutation', fill='Variant Class') +
    theme_classic() +
    theme +
    scale_fill_manual(values = rev(brewer.pal(length(levels(factor(variant_class))), 'Spectral')) ) +
    theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

    # cnv heatmap
    cnv_heatmap <- ggplot(data=cnv_dat) +
    aes(x=sample, y=query_gene, fill=cnv_type) +
    geom_tile(color='white', size=1) +
    labs(x='',y='CNV', fill='CNV Type') +
    theme_classic() +
    theme +
    scale_fill_manual(values = rev(brewer.pal(5,"RdBu")),drop=FALSE, labels = c("deletion", "", "neutral", "", "amplification"))

    # arrange using ggpubr
    oncoprint_chart <- ggarrange(mut_burden_barplot ,phenotype,mut_heatmap,cnv_heatmap, nrow=4, align = 'v', heights = plot_heights)
    return(oncoprint_chart)
}