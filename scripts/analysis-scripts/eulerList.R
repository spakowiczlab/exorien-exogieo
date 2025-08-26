eulerList <- function(des.res){
  taxalist <- list(
    ORIEN = des.res$res.ORIEN$Taxa,
    TCGA = des.res$res.TCGA$Taxa,
    "16s" = des.res$res.16s$microbe
  )
  
  return(taxalist)
}
