calcL2FC <- function(){
  subset_matched_clin <- function(){
    tmp1 <- read_excel("/fs/ess/PAS1695/projects/exogieo/221202_Spakowicz_GSL-PY-2916/Copy of DeIDd_Jin_GI_GreaterThan50_Frozen_Shipment1_ClinicalData_09.2020-updated in 2022.xlsx") %>%
      # select(`TCC ID`, `RNASeq-T SL ID`, `Specimen Type`) %>%
      rename("link" =  `TCC ID`,
             "sample" = `RNASeq-T SL ID`) %>%
      mutate(EO.stat = "lo")
    tmp2 <- read_excel("/fs/ess/PAS1695/projects/exogieo/221202_Spakowicz_GSL-PY-2916/Copy of DeIDd_Jin_GI_lessThan50_Frozen_Shipment1_ClinicalData_09.2020_lcm.xlsx") %>%
      # select(`TCC ID`, `RNASeq-T SL ID`, `Specimen Type`) %>%
      rename("link" =  `TCC ID`,
             "sample" = `RNASeq-T SL ID`) %>%
      mutate(EO.stat = "eo")
    tmp3 <- read_excel("/fs/ess/PAS1695/projects/exogieo/221202_Spakowicz_GSL-PY-2916/Spackowicz Sample DNA.xls") %>%
      # select(Sample, `Subject ID 1`) %>%
      rename("sample16s" = Sample,
             "link" = `Subject ID 1`) %>%
      drop_na(link) 
    
    meta <- bind_rows(tmp1,tmp2) %>%
      inner_join(tmp3) %>%
      select(sample16s, link, `Tumor/Normal`, EO.stat) %>%
      filter(duplicated(sample16s)) %>%
      filter(`Tumor/Normal` == "Tumor") %>%
      mutate(id = as.numeric(gsub(".*-0", "", sample16s)))
    
    return(meta)
  }
  
  exoRAtowide <- function(data, taxlev){
    tmp <- data
    tmp$Taxa <- tmp[[taxlev]]
    tmp.wide<- tmp %>%
      dplyr::filter(!is.na(Taxa)) %>%
      group_by(sample, Taxa) %>%
      summarise(ra = sum(counts, na.rm = T)) %>%
      pivot_wider(names_from = "Taxa", values_from = "ra")
    
    tmp.wide[is.na(tmp.wide)] <- 0
    return(tmp.wide)
  }
  
  exoToDF <- function(taxalevels = c("domain", "kingdom", "phylum", "class", 
                                     "order", "family", "genus", "species"), 
                      data){
    w.ls <- lapply(taxalevels, function(x) exoRAtowide(data, x))
    w.df <- reduce(w.ls, function(x,y) left_join(x,y)) 
    return(w.df)
  }
  
  exo.ra <- read.csv("../../../data/16s_taxa_allLevels_counts_wide.csv", stringsAsFactors = F)
  meta <- subset_matched_clin()
  
  exo.ra2 <- exo.ra %>%
    mutate(totcounts = sum(counts), .by = sample) %>%
    mutate(counts = (counts/totcounts) * 1e5) %>%
    select(-totcounts)
  
  exo.w <- exoToDF(data = exo.ra2)
  mics <- colnames(exo.w)[-1]
  
  modin <- meta %>%
    mutate(sample = id) %>%
    inner_join(exo.w)
  
  l2fc <- modin %>%
    select(id, EO.stat, all_of(mics)) %>%
    pivot_longer(-c("id", "EO.stat"), names_to = "microbe", values_to = "counts") %>%
    mutate(counts = counts + 1) %>%
    group_by(microbe,EO.stat) %>%
    summarise(gmean = mean(counts)) %>%
    pivot_wider(names_from = EO.stat, values_from = gmean) %>%
    mutate(l2fc = log(lo/eo, base = 2)) %>%
    # filter(is.finite(l2fc)) %>%
    select(microbe, l2fc)
  
  return(l2fc)
}