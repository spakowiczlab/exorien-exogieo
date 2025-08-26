library(tidyverse)

poolCorrsHeat <- function(micim.full){

micim.full <- micim.full %>%
  filter(sigcode != "Neither")


# micim.sortmics <- micim.full %>%
#   filter(onset == "LO") %>%
#   select(microbe, im.cell, TCGA.est) %>%
#   pivot_wider(names_from = im.cell, values_from = TCGA.est, values_fill = 0) %>%
#   column_to_rownames(var = "microbe") %>%
#   dist() %>% 
#   hclust()
# micord <- micim.sortmics$labels[micim.sortmics$order]
# 
# micim.sortcells <- micim.full %>%
#   filter(onset == "LO") %>%
#   select(microbe, im.cell, TCGA.est) %>%
#   pivot_wider(names_from = microbe, values_from = TCGA.est, values_fill = 0) %>%
#   column_to_rownames(var = "im.cell") %>%
#   dist() %>% 
#   hclust()
# cellord <- micim.sortcells$labels[micim.sortcells$order]

micim.sortmics <- micim.full %>%
  select(microbe, onset, ORIEN.est, TCGA.est, im.cell) %>%
  pivot_longer(c("ORIEN.est", "TCGA.est"), names_to = "source",
               values_to = "correlation") %>%
  group_by(microbe, im.cell) %>%
  summarise(correlation = mean(correlation)) %>%
  pivot_wider(names_from = im.cell, values_from = correlation, values_fill = 0) %>%
  column_to_rownames(var = "microbe") %>%
  dist() %>%
  hclust()
micord <- micim.sortmics$labels[micim.sortmics$order]

micim.sortcells <- micim.full %>%
  select(microbe, onset, ORIEN.est, TCGA.est, newlabs) %>%
  pivot_longer(c("ORIEN.est", "TCGA.est"), names_to = "source",
               values_to = "correlation") %>%
  group_by(microbe, newlabs) %>%
  summarise(correlation = mean(correlation)) %>%
  pivot_wider(names_from = microbe, values_from = correlation, values_fill = 0) %>%
  column_to_rownames(var = "newlabs") %>%
  dist() %>%
  hclust()
cellord <- micim.sortcells$labels[micim.sortcells$order]

micim.heat <- micim.full %>%
  mutate(correlation.ORIEN = ifelse(ORIEN.p < 0.05, ORIEN.est, NA),
         correlation.TCGA = ifelse(TCGA.p < 0.05, TCGA.est, NA)) %>%
  pivot_longer(c("correlation.ORIEN", "correlation.TCGA"),
               names_to = "corr.source", values_to = "correlation") %>%
  mutate(newlabs = fct_relevel(newlabs, cellord),
         microbe = fct_relevel(microbe, micord))

return(micim.heat)

}
