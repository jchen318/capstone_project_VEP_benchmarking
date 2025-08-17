Warning: cannot create dir '/CETP', reason 'Read-only file system' 

Warning: 2 arguments not used by format 'sample_variants_filtered.csv' 

Error in fread(paste(INPUT_VARIANTs_DIR, fileName, sep = "/"), drop = 1) : File '/sample_variants_filtered.csv' does not exist or is non-readable. getwd()=='/Users/andrewchen/capstone_project_VEP_benchmarking/VEP_benchmarking' Calls: setUpEnvAndRun ... apply -> FUN -> fread -> stopf -> raise_condition -> signal 
Traceback: 
 [1] "12: signal(obj)"                                                                                      
 [2] "11: raise_condition(stop, gettextf(fmt, ..., domain = domain), c(class, "                             
 [3] "        \"simpleError\", \"error\", \"condition\"))"                                                  
 [4] "10: stopf(\"File '%s' does not exist or is non-readable. getwd()=='%s'\", "                           
 [5] "        file, getwd())"                                                                               
 [6] "9: fread(paste(INPUT_VARIANTs_DIR, fileName, sep = \"/\"), drop = 1)"                                 
 [7] "8: FUN(newX[, i], ...)"                                                                               
 [8] "7: apply(blocks, 1, function(row) {"                                                                  
 [9] "       chrom = row[[\"chr\"]]"                                                                        
[10] "       block = row[[\"block\"]]"                                                                      
[11] "       fileName = sprintf(ALL_VARIANTs_PATH, chrom, block)"                                           
[12] "       allVariants = fread(paste(INPUT_VARIANTs_DIR, fileName, sep = \"/\"), "                        
[13] "           drop = 1)"                                                                                 
[14] "       fileName = sprintf(UNIQUE_VARIANTs_PATH, chrom, block)"                                        
[15] "       uniqueVariants = fread(paste(INPUT_VARIANTs_DIR, fileName, "                                   
[16] "           sep = \"/\"))"                                                                             
[17] "       mergedVariants = merge(allVariants, uniqueVariants[Gene == "                                   
[18] "           ENSEMBL_ID], by.x = \"variant_name\", by.y = \"#Uploaded_variation\")"                     
[19] "       return(mergedVariants)"                                                                        
[20] "   })"                                                                                                
[21] "6: eval(exprs[i], envir)"                                                                             
[22] "5: eval(exprs[i], envir)"                                                                             
[23] "4: sys.source(file, envir = env, keep.source = F)"                                                    
[24] "3: withCallingHandlers(expr, packageStartupMessage = function(c) tryInvokeRestart(\"muffleMessage\"))"
[25] "2: suppressPackageStartupMessages(sys.source(file, envir = env, "                                     
[26] "       keep.source = F))"                                                                             
[27] "1: setUpEnvAndRun(\"Parsing variants\", \"lib/parseVariants.R\")"                                     

