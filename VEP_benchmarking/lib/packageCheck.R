# ---- Safe defaults so packageCheck.R works standalone -----------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!exists("CRAN_LIBRARIES")) {
  CRAN_LIBRARIES <- c(
    "yaml","data.table","readr","dplyr","tidyr","tibble",
    "purrr","stringr","ggplot2","magrittr","optparse",
    # extras the pipeline uses / you saw missing:
    "ggtext","forcats","gridExtra","Rmisc","reshape2","logr","remotes","retry"
  )
}

if (!exists("BIOCONDUCTOR_LIBRARIES")) {
  BIOCONDUCTOR_LIBRARIES <- c("qvalue")
}

ensure_cran <- function(pkgs) {
  need <- setdiff(pkgs, rownames(installed.packages()))
  if (length(need)) install.packages(need, dependencies = TRUE)
}

ensure_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
}
# ---------------------------------------------------------------------------

CRAN_LIBRARIES = c("yaml", "data.table", "dplyr", "tidyr", "stringr",
                   "ggtext", "glue", "forcats", "gridExtra", "fs", "Rmisc",
                   "reshape2", "stringr", "logr", "BiocManager", "qvalue", "remotes", "retry")
#BIOCONDUCTOR_LIBRARIES = c("EnsDb.Hsapiens.v86")
GITHUB_LIBRARIES = c("jweile/yogiroc")

# Install needed libraries
cat("Checking missing libraries\n")
missingLibs = CRAN_LIBRARIES[!CRAN_LIBRARIES %in% installed.packages()[, "Package"]]
if (length(missingLibs)) {
  cat("Installing missing CRAN libraries:", paste0(missingLibs, collapse = ", "), "\n")
  install.packages(missingLibs, lib = Sys.getenv("R_LIBS_USER"))
}

missingLibs = BIOCONDUCTOR_LIBRARIES[!BIOCONDUCTOR_LIBRARIES %in% installed.packages()[, "Package"]]
if (length(missingLibs)) {
  cat("Installing missing Bioconductor libraries:", paste0(missingLibs, collapse = ", "), "\n")
  BiocManager::install(missingLibs, lib = Sys.getenv("R_LIBS_USER"))
}

library(stringr)
missingLibs = str_split(GITHUB_LIBRARIES, fixed("/"), simplify = T)[,2]
missingLibs = GITHUB_LIBRARIES[which(!missingLibs %in% installed.packages()[, "Package"])]
if (length(missingLibs)) {
  cat("Installing missing Github libraries:", paste0(missingLibs, collapse = ", "), "\n")
  remotes::install_github(missingLibs, lib = Sys.getenv("R_LIBS_USER"))
}