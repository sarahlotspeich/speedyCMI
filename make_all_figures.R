setwd("~/Documents/")
fig_files = paste0("speedyCMI/figures/", list.files(path = "speedyCMI/figures/", pattern = ".R"))
sapply(X = fig_files, FUN = source)
