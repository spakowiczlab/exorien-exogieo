prepareVolLDA <- function(des.res, l2fc, nlda){
  
  # Format for combined volcanoplot. Swap l2fc direction for deseq2 results to put EO on the right of the plots.
  volplot.tcga <- des.res$res.TCGA %>%
    mutate(colgroup = ifelse(padj < 0.05, Taxa.Level, NA)) %>%
    mutate(collevel = fct_relevel(colgroup, c("kingdom", "phylum", "class", "order", "family", "genus")),
           dataset = "TCGA")
  volplot.orien<- des.res$res.ORIEN %>%
    mutate(colgroup = ifelse(padj < 0.05, Taxa.Level, NA)) %>%
    mutate(collevel = fct_relevel(colgroup, c("kingdom", "phylum", "class", "order", "family", "genus")),
           dataset = "ORIEN")
  volplot.16s <- des.res$res.16s %>%
    left_join(l2fc) %>%
    mutate(taxcode = str_sub(microbe,1,1),
           Taxa = microbe,
           Taxa.Level = case_when(taxcode == "k" ~ "kingdom",
                                  taxcode == "p" ~ "phylum",
                                  taxcode == "c" ~ "class",
                                  taxcode == "o" ~ "order",
                                  taxcode == "f" ~ "family",
                                  taxcode == "g" ~ "genus",
                                  taxcode == "s" ~ "species"),
           colgroup = ifelse(padj < 0.05, Taxa.Level, NA)) %>%
    # drop_na(colgroup) %>%
    mutate(
      # collevel = fct_relevel(colgroup, c("kingdom", "phylum", "class", "order", "family", "genus")),
      dataset = "16S",
      log2FoldChange = l2fc)
  
  volplot.in <- bind_rows(volplot.16s, volplot.orien, volplot.tcga) %>%
    mutate(collevel = fct_relevel(colgroup, c("kingdom", "phylum", "class", "order", "family", "genus")))
  
  
  # Format for combined LDA
  arrangeLDA <- function(ntax, plot.input, d){
    res.lda <- plot.input %>%
      arrange(padj) %>%
      mutate(sigrank = row_number(),
             taxname = gsub("^\\w__", "", Taxa))
    
    tmp <- res.lda %>%
      filter(sigrank <= ntax) %>%
      arrange(log2FoldChange) %>%
      mutate(dataset = d)
    tmp$tax.fact <- fct_relevel(tmp$taxname, tmp$taxname)
    
    return(tmp)
    
  }
  
  lda.in <- bind_rows(arrangeLDA(nlda, volplot.16s, "16S"),
                      arrangeLDA(nlda, volplot.orien, "ORIEN"),
                      arrangeLDA(nlda, volplot.tcga, "TCGA"))
  
  outs <- list(volplot.in, lda.in)
  names(outs) <- c("volplot", "LDA")
  return(outs)
}
