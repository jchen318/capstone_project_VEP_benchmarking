# --- Fix library paths on Windows ---
user_lib <- "C:/Users/jchen/AppData/Local/R/win-library/4.5"
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)

# Use only proper vector entries for lib paths
.libPaths(c(user_lib, "C:/Program Files/R/R-4.5.1/library"))

# Also fix env vars so child Rscript sees the right path
Sys.setenv(R_LIBS_USER = user_lib)
Sys.unsetenv("R_LIBS")  # avoid broken semicolon value
# --- end fix ---

# vep_benchmarking_smoke.r  — self-contained smoke test
# -----------------------------------------------------

# --- Pre-install GitHub dep into your user library (avoids child build) -----
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# make sure user lib is first and writable
user_lib <- "C:/Users/jchen/AppData/Local/R/win-library/4.5"
.libPaths(c(user_lib, file.path(R.home("library"))))
Sys.setenv(R_LIBS_USER = user_lib); Sys.unsetenv("R_LIBS"); Sys.unsetenv("R_LIBS_SITE")

# install yogiroc once if missing
if (!requireNamespace("yogiroc", quietly = TRUE)) {
  remotes::install_github("jweile/yogiroc", dependencies = TRUE, lib = .libPaths()[1])
}
# ---------------------------------------------------------------------------

# ---- 0) Dependencies --------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

CRAN_LIBRARIES <- c(
  "yaml","data.table","readr","dplyr","tidyr","tibble",
  "purrr","stringr","ggplot2","magrittr","optparse",
  "ggtext","forcats","gridExtra","Rmisc","reshape2","logr","remotes","retry"
)
BIOCONDUCTOR_LIBRARIES <- c("qvalue")

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

cat("Checking missing libraries\n")
ensure_cran(CRAN_LIBRARIES)
ensure_bioc(BIOCONDUCTOR_LIBRARIES)

# ---- 1) Paths ---------------------------------------------------------------
repo_dir <- "C:/Users/jchen/Desktop/Capstone_project/VEP_benchmarking"
setwd(repo_dir)

# Ensure folders exist
if (!dir.exists("logs")) dir.create("logs")
if (!dir.exists("output")) dir.create("output")

# define OUTPUT_PATH so main.R can find it
OUTPUT_PATH <- file.path(repo_dir, "output")

# Config must exist
if (!file.exists("config.local.yaml")) {
  stop("Missing config.local.yaml in: ", repo_dir)
}

# ---- 2) Run the pipeline (with child startup profile) -----------------------
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) stop("Could not find Rscript.exe at: ", rscript)

# Create a tiny startup profile for the child Rscript so it knows OUTPUT_PATH
child_profile <- file.path(repo_dir, "child_profile.R")
writeLines(c(
  'options(repos = c(CRAN = "https://cloud.r-project.org"))',
  # keep library path hygiene in the child too
  'user_lib <- "C:/Users/jchen/AppData/Local/R/win-library/4.5"',
  'if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)',
  '.libPaths(c(user_lib, file.path(R.home("library"))))',
  'Sys.setenv(R_LIBS_USER = user_lib); Sys.unsetenv("R_LIBS"); Sys.unsetenv("R_LIBS_SITE")',
  # make OUTPUT_PATH visible to main.R in the child process
  sprintf('OUTPUT_PATH <- "%s"', file.path(repo_dir, "output"))
), child_profile)

args <- c("main.R", "config.local.yaml", "logs")
cat("Running: ", rscript, " ", paste(args, collapse=" "), "\n", sep = "")

# Pass clean env + the child profile to Rscript
child_env <- c(
  paste0("R_LIBS_USER=", .libPaths()[1]),
  "R_LIBS=",
  "R_LIBS_SITE=",
  paste0("R_PROFILE_USER=", child_profile)  # <<< ensures OUTPUT_PATH is defined
)

status <- system2(rscript, args, env = child_env)
if (status != 0) stop("Pipeline run failed with exit status: ", status)
# -----------------------------------------------------------------------------

# ---- 3) Show result preview -------------------------------------------------
if (file.exists(file.path(repo_dir, "pvals.csv"))) {
  cat("\n✅ Pipeline completed. Preview of pvals.csv:\n")
  print(readr::read_csv("pvals.csv", show_col_types = FALSE, n_max = 10))
} else {
  cat("\n⚠️ Completed but pvals.csv not found. Check output/ and logs/.\n")
}




