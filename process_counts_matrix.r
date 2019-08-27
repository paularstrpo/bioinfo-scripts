setwd('~/Desktop/TCGA htseq counts')
cancer_type <- 'READ'
pheno <- '_rdata/TCGA-pan-master-phenoData.rdata'


load(pheno)
pheno <- phenoData[phenoData$project_id == paste('TCGA', cancer_type, sep='-') , ]

setwd(cancer_type)

read.gz <- function(x) { x <- read.delim(gzfile(x), header=FALSE, col.names=c('ensid', 'counts'))}

flist <- list.files()
flist <- flist[1:nrow(pheno)]

counts <- NULL
for (i in 1:length(flist)) {counts[[flist[i]]] <- read.gz(flist[i])}

counts <- do.call(cbind, counts)
rownames(counts) <- counts[, 1]

counts <- counts[, substr(colnames(counts), nchar(colnames(counts)) - 5 , nchar(colnames(counts))) == 'counts']
ncol(counts)

matchID <- data.frame(colnames(counts))
colnames(matchID) <- c('file_name')
matchID$file_name <- substr(as.character(matchID$file_name), 1, nchar(as.character(matchID$file_name)) - nchar(".counts"))
matchID <- merge(matchID, pheno[, c("file_name", "sample_id")], by=c('file_name'))
colnames(counts) <- matchID$sample_id

rm(matchID, flist, i, phenoData)

READ.counts <- counts
READ.phenoData <- pheno
setwd('~/Desktop/TCGA htseq counts/_rdata')
save(READ.counts, READ.phenoData, file='READ_raw-counts+phenoData.rdata')
