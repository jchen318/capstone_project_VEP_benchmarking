# install.R
pkgs <- c("yaml","data.table","readr","dplyr","tidyr","tibble",
          "purrr","stringr","ggplot2","magrittr","parallel","optparse")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")



