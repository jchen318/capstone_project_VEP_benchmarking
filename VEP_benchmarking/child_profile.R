options(repos = c(CRAN = "https://cloud.r-project.org"))
user_lib <- "C:/Users/jchen/AppData/Local/R/win-library/4.5"
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, file.path(R.home("library"))))
Sys.setenv(R_LIBS_USER = user_lib); Sys.unsetenv("R_LIBS"); Sys.unsetenv("R_LIBS_SITE")
OUTPUT_PATH <- "C:/Users/jchen/Desktop/Capstone_project/VEP_benchmarking/output"
