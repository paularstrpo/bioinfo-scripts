library(tidyverse)
library(fgsea)

gsea_lister <- function(sig_list, rank_list, sig_name){
  
  # perform the gsea, include the sample comparison
  y <- lapply(names(rank_list), function(i){
    
    x <- fgsea(pathways=sig_list, stats=rank_list[[i]], nperm=1000) %>%
    as_tibble() %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    mutate(sample_comparison=i, signature_name=sig_name)
  return(x)
})
  
  y <-  do.call(rbind, y) %>%
  as_tibble()
  
  return(y)
}