validate16sStackedBar <- function(orderby){
  
  exorien.clin <- read.csv("../data/clinical_ORIEN_curated.csv")
  
  loadCounts <- function(){
    exo.ra <- read.csv("/fs/ess/PAS1695/projects/exorien/data/drake-output/2022-03-16/2022-03-16_unnormalized-microbes_humanRNAfilt_w_taxonomy.csv")
    
    goodsamps <- intersect(exo.ra$sample, exorien.clin$RNASeq)
    
    exo.ra.filt <- exo.ra %>%
      filter(sample %in% goodsamps)
    
    return(exo.ra.filt)
  }
  
  orien.counts <- loadCounts()
  tcga.counts <- read.csv("/fs/ess/PAS1695/projects/exogieo/data/ORIEN-processing/2022-03-16_unnorm-mics_filt.csv", stringsAsFactors = F) 
  
  
  samp_match <- read.csv("../data/16s_sample_matching.csv")
  krakenmet = read.delim("/fs/ess/PAS1695/exoticpipe/external-data/kraken2-metaphlan-noplants.txt",
                         header = F, stringsAsFactors = F) %>%
    rename("Taxonomy" = "V1")
  counts.16s <- read.csv("../data/16s_counts_long.csv") %>%
    rename("microbe" = name,
           "counts" = new_est_reads)
  
  taxkey <- krakenmet %>%
    dplyr::filter(grepl("s__", Taxonomy)) %>%
    mutate(microbe = gsub(".*s__(.*)", "\\1", Taxonomy),
           microbe = make.names(microbe))
  
  counts.16s.tax <- counts.16s %>%
    left_join(taxkey) %>%
    select(sample, microbe, Taxonomy, counts) %>%
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
  
  get_ra <- function(exora, level){
    exora <- exora %>%
      filter(microbe != "Homo.sapiens") %>%
      select(sample, !! sym(level), counts) %>%
      group_by(sample, !! sym(level)) %>%
      summarize(count = sum(counts)) %>%
      ungroup() %>%
      group_by(sample) %>%
      mutate(rel.abun = count / sum(count)) %>%
      ungroup()
    
    return(exora)
  }
  
  
  counts <- list(counts.16s.tax, tcga.counts, orien.counts)
  names(counts) <- c("16s", "TCGA RNASeq", "ORIEN RNASeq")
  
  ras <- lapply(names(counts), function(x) get_ra(counts[[x]], "phylum") %>% 
                  mutate(samptype = x,
                         sample = as.character(sample))) %>%
    bind_rows()
  
  large.phyl <- ras %>%
    group_by(phylum) %>%
    summarize(median.ra = median(rel.abun)) %>%
    arrange(desc(median.ra)) %>%
    mutate(x = row_number()) %>%
    dplyr::filter(x <= 8)
  
  sampord <- ras %>%
    filter(phylum == orderby) %>%
    arrange(desc(rel.abun))
  sampord <- sampord$sample
  
  stackedin <- ras %>%
    mutate(Phylum = ifelse(phylum %in% large.phyl$phylum, phylum, "Other"),
           Phylum = gsub("p__", "", Phylum),
           # Phylum = fct_relevel(Phylum, "Firmicutes"),
           Phylum = fct_relevel(Phylum, "Other", after = Inf),
           sample = fct_relevel(sample, sampord))
  
  return(stackedin)
  
}
