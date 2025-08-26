correlationPoolSource <- function(){
  micim.TCGA <- read.csv("/fs/ess/PAS1695/projects/exogieo/data/correlation-data/TCGA-results_1.csv") %>%
    rename("im.cell" = "V2",
           "TCGA.est" = "estimate",
           "TCGA.p" = "p.value") %>%
    select(microbe, im.cell, onset, TCGA.est, TCGA.p) 
  # mutate(datset = "TCGA")
  
  micim.ORIEN <- read.csv("/fs/ess/PAS1695/projects/exogieo/data/correlation-data/ORIEN-results_1.csv") %>%
    rename("im.cell" = "V2",
           "ORIEN.est" = "estimate",
           "ORIEN.p" = "p.value") %>%
    select(microbe, im.cell, onset, ORIEN.est, ORIEN.p)
  
  oldlabs <- c("B.cells.naive", "B.cells.memory", 
               "Plasma.cells", "T.cells.CD8", "T.cells.CD4.naive", 
               "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated", "T.cells.follicular.helper",
               "T.cells.gamma.delta", "T.cells.regulatory..Tregs.", "NK.cells.resting", 
               "NK.cells.activated", "Monocytes", "Macrophages.M0", 
               "Macrophages.M1", "Macrophages.M2", "Dendritic.cells.resting", 
               "Dendritic.cells.activated", "Mast.cells.resting", "Mast.cells.activated", 
               "Eosinophils", "Neutrophils", "TILs")
  
  newlabs <- c("Naive B Cells", "Memory B Cells", 
               "Plasma Cells", "CD8 T Cells", "Naive CD4 T Cells", 
               "Resting Memory CD4 T Cells", "Activated Memory CD4 T Cells", "Follicular Helper T Cells",
               "Gamma Delta T Cells", "Regulatory T Cells", "Resting NK Cells", 
               "Activated NK Cells", "Monocytes", "M0 Macrophages",
               "M1 Macrophages", "M2 Macrophages", "Resting Dendritic Cells", 
               "Activated Dendritic Cells", "Resting Mast Cells", "Activated Mast Cells", 
               "Eosonophils", "Neutrophils", "TILs")
  
  lab.df <- as.data.frame(cbind(im.cell = oldlabs, newlabs))
  
  micim.mics <- intersect(micim.ORIEN$microbe, micim.TCGA$microbe)
  
  micim.full <- inner_join(micim.ORIEN, micim.TCGA) %>%
    # filter(microbe %in% micim.mics) %>%
    left_join(lab.df) %>%
    drop_na() %>%
    mutate(sigcode = case_when(ORIEN.p < 0.05 & TCGA.p < 0.05 ~ "Both",
                               ORIEN.p < 0.05 & TCGA.p >= 0.05 ~ "ORIEN",
                               ORIEN.p >= 0.05 & TCGA.p < 0.05 ~ "TCGA",
                               ORIEN.p >= 0.05 & TCGA.p >= 0.05 ~ "Neither"),
           mlab = ifelse(sigcode == "Both", microbe, NA)) 
  
  return(micim.full)
}