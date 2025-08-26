loadDESeqRes <- function(){
  des.res <- list(
  res.ORIEN = read.csv("../data/deseq2_all-cohort_ORIEN-data.csv"),
  res.TCGA = read.csv("../data/deseq2_all-cohort.csv"),
  res.16s = read.csv("../data/16s_modelling_eo-lo.csv")
  )
  
  return(des.res)
}
