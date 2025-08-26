### Holds functions for gathering plotted data for exogieo

analysis.scripts <- list.files("analysis-scripts/", full.names = T)
lapply(analysis.scripts, source)

### Data related to 16s
dist.16s <- distanceBoxplot()
stackedbar.16s <- validate16sStackedBar("p__Proteobacteria")

save(dist.16s, file = "../data/prepared-figure-data/distances_16s.rda")
save(stackedbar.16s, file = "../data/prepared-figure-data/stackedbar_16s.rda")

### Results of EO/LO microbe enrichment analysis
des.res <- loadDESeqRes()
euler.list <- eulerList(des.res)
l2fc16s <- calcL2FC()

vol.lda.dat <- prepareVolLDA(des.res, l2fc16s, 10)

save(euler.list, file = "../data/prepared-figure-data/euler_taxalist.rda")
save(vol.lda.dat, file = "../data/prepared-figure-data/vol-lda.rda")

### Figures related to correlations. The correlation compute is large, see processing scripts for that
micim.full <- correlationPoolSource()
micim.heat <- poolCorrsHeat(micim.full = micim.full)

save(micim.full, file = "../data/prepared-figure-data/scatterplot_corrs.rda")
save(micim.heat, file = "../data/prepared-figure-data/heatmap_corrs.rda")
