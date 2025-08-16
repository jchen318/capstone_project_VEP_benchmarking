Warning: 2 arguments not used by format 'sample_variants_filtered.csv'

Error in fread(paste(INPUT_VARIANTs_DIR, fileName, sep = "/"), drop = 1) : File
'/sample_variants_filtered.csv' does not exist or is non-readable.
getwd()=='C:/Users/jchen/Desktop/Capstone_project/VEP_benchmarking' Calls:
setUpEnvAndRun ... apply -> FUN -> fread -> stopf -> raise_condition -> signal
Traceback:
12: signal(obj)
11: raise_condition(stop, gettextf(fmt, ..., domain = domain), c(class, 
        "simpleError", "error", "condition"))
10: stopf("File '%s' does not exist or is non-readable. getwd()=='%s'", 
        file, getwd())
9: fread(paste(INPUT_VARIANTs_DIR, fileName, sep = "/"), drop = 1)
8: FUN(newX[, i], ...)
7: apply(blocks, 1, function(row) {
       chrom = row[["chr"]]
       block = row[["block"]]
       fileName = sprintf(ALL_VARIANTs_PATH, chrom, block)
       allVariants = fread(paste(INPUT_VARIANTs_DIR, fileName, sep = "/"), 
           drop = 1)
       fileName = sprintf(UNIQUE_VARIANTs_PATH, chrom, block)
       uniqueVariants = fread(paste(INPUT_VARIANTs_DIR, fileName, 
           sep = "/"))
       mergedVariants = merge(allVariants, uniqueVariants[Gene == 
           ENSEMBL_ID], by.x = "variant_name", by.y = "#Uploaded_variation")
       return(mergedVariants)
   })
6: eval(exprs[i], envir)
5: eval(exprs[i], envir)
4: sys.source(file, envir = env, keep.source = F)
3: withCallingHandlers(expr, packageStartupMessage = function(c) tryInvokeRestart("muffleMessage"))
2: suppressPackageStartupMessages(sys.source(file, envir = env, 
       keep.source = F))
1: setUpEnvAndRun("Parsing variants", "lib/parseVariants.R")

