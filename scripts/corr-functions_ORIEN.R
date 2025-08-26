args = commandArgs(trailingOnly = T)

library(tidyverse)
library(broom)
library(rlist)

# Load data

load("/fs/ess/PAS1695/projects/exogieo/data/correlation-data/ORIEN-inputs.rda")

# Define function

mass_correlate <- function(inputdf, names1, names2, lab1, lab2){
  cor.out <- list()
  for(eo in c("EO", "LO")){
    input.tmp <- inputdf %>%
      filter(agestat == eo)
    cor.out[[eo]] <- lapply(names1, function(m) 
      lapply(names2,  function(i) try(cor.test(input.tmp[[m]],
                                               input.tmp[[i]],
                                               method = "spearman") %>%
                                        tidy() %>%
                                        mutate(names2 = i, 
                                               names1 = m, 
                                               onset = eo))))
    cor.out[[eo]] <- lapply(cor.out[[eo]], function(x) list.clean(x, is.character))
  }
  
  test <- lapply(cor.out, function(x)
    lapply(x, function(y) bind_rows(y)))
  test2 <- lapply(test, function(x) bind_rows(x))
  cor.out.df <- bind_rows(test2) %>%
    rename(!!lab1 := names1,
           !!lab2 := names2)
  
  return(cor.out.df)
}


cor.can <- mass_correlate(micimexp.all, microbes, corr.vars[[args[1]]], "microbe", "V2")

write.csv(cor.can, paste0("/fs/ess/PAS1695/projects/exogieo/data/correlation-data/ORIEN-results_",
                         args[1], ".csv"), row.names = F)
