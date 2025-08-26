distanceBoxplot <- function(){
  samp_match <- read.csv("../data/16s_sample_matching.csv")
  krakenmet = read.delim("/fs/ess/PAS1695/exoticpipe/external-data/kraken2-metaphlan-noplants.txt",
                         header = F, stringsAsFactors = F) %>%
    dplyr::rename("Taxonomy" = "V1")
  exogieo <- read.csv("../data/16s_counts_long.csv") %>%
    dplyr::rename("microbe" = name,
                  "counts" = new_est_reads)
  exotic <- readRDS("/fs/ess/PAS1695/projects/exotic/data/drake-output/2022-03-08/tcc.counts.RDS")
  
  get_levels <- function(exora){
    exora <- exora %>%
      dplyr::select(sample, microbe, Taxonomy, counts) %>%
      mutate(domain = gsub("(d__\\w+).*", "\\1", Taxonomy),
             kingdom = ifelse(grepl("k__", Taxonomy),gsub(".*(k__\\w+).*", "\\1", Taxonomy), NA),
             phylum = ifelse(grepl("p__", Taxonomy),gsub(".*(p__\\w+).*", "\\1", Taxonomy), NA),
             order = ifelse(grepl("o__", Taxonomy),gsub(".*(o__\\w+).*", "\\1", Taxonomy), NA),
             class = ifelse(grepl("c__", Taxonomy),gsub(".*(c__\\w+).*", "\\1", Taxonomy), NA),
             family = ifelse(grepl("f__", Taxonomy),gsub(".*(f__\\w+).*", "\\1", Taxonomy), NA),
             genus = ifelse(grepl("g__", Taxonomy),gsub(".*(g__\\w+).*", "\\1", Taxonomy), NA),
             species = gsub(".*(s__.*)", "\\1", Taxonomy)
      ) %>%
      mutate(genus = ifelse(is.na(genus), paste0("g__unclassified-", species), genus),
             family = ifelse(is.na(family), paste0("f__unclassified-", gsub("g__unclassified-", "", genus)), family),
             order = ifelse(is.na(order), paste0("o__unclassified-", gsub("f__unclassified-", "", family)), order),
             class = ifelse(is.na(class), paste0("c__unclassified-", gsub("o__unclassified-", "", order)), class),
             phylum = ifelse(is.na(phylum), paste0("p__unclassified-", gsub("c__unclassified-", "", class)),phylum),
             kingdom = ifelse(is.na(kingdom), paste0("k__unclassified-", gsub("p__unclassified-", "", phylum)), kingdom),
             domain = ifelse(is.na(domain), paste0("d__unclassified-", gsub("k__unclassified-", "", kingdom)), domain)) %>%
      dplyr::select(sample, microbe, Taxonomy, domain, kingdom, phylum, class, order, family, genus, species, counts)
    return(exora)
  }
  
  taxkey <- krakenmet %>%
    dplyr::filter(grepl("s__", Taxonomy)) %>%
    mutate(microbe = gsub(".*s__(.*)", "\\1", Taxonomy),
           microbe = make.names(microbe))
  
  exotic.all <- exotic %>%
    pivot_longer(-sample, names_to = "microbe", values_to = "counts") %>%
    left_join(taxkey) %>%
    get_levels() %>%
    filter(domain == "d__Bacteria") # Only include bacteria to better match 16s
  
  exogieo.all <- exogieo %>%
    left_join(taxkey) %>%
    get_levels()
  
  exogieo.prev <- exogieo.all %>%
    select(sample, genus, counts) %>% 
    group_by(sample, genus) %>% 
    summarize(count = sum(counts)) %>%
    mutate(prev = ifelse(count > 0, 1, 0)) %>%
    select(-count)
  
  exotic.prev <- exotic.all %>%
    select(sample, genus, counts) %>% 
    group_by(sample, genus) %>% 
    summarize(count = sum(counts)) %>%
    mutate(prev = ifelse(count > 0, 1, 0)) %>%
    select(-count)
  
  samp_match <- samp_match %>%
    mutate(crop16s = substr(samp16s, 9, length(samp16s)))
  
  exotic.prev <- exotic.prev %>%
    filter(sample %in% samp_match$RNAseq)
  
  exogieo.prev$sample <- as.character(exogieo.prev$sample)
  
  combined.prev <- rbind(exotic.prev, exogieo.prev) %>%
    pivot_wider(names_from = genus, values_from = prev, values_fill = 0) %>%
    column_to_rownames(var = "sample") %>%
    as.matrix()
  
  set.seed(19971030)
  dist <- vegdist(combined.prev, method = "bray")
  
  dist.res <- as.matrix(dist) %>%
    as.data.frame() %>%
    rownames_to_column(var = "RNAseq") %>%
    pivot_longer(-RNAseq, names_to = "samp16s", values_to = "distance") %>%
    filter(RNAseq %in% samp_match$RNAseq,
           samp16s %in% samp_match$crop16s)
  
  pt.rna <- samp_match %>%
    select(patient.id, RNAseq) %>%
    dplyr::rename("pt.id.rna" = patient.id) %>%
    filter(!is.na(RNAseq)) %>%
    distinct()
  
  pt.16s <- samp_match %>%
    select(patient.id, crop16s, specimen.type) %>%
    dplyr::rename("pt.id.16s" = patient.id,
                  "samp16s" = crop16s)
  
  
  dist.res.lab <- dist.res %>%
    left_join(pt.rna) %>%
    left_join(pt.16s) %>%
    mutate(paired = (pt.id.16s == pt.id.rna)) %>%
    mutate(paired = ifelse(is.na(paired), "CONTROL", paired))
  
  return(dist.res.lab)
}
