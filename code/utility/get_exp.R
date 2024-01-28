# Get expression from getGEOSuppFiles

get_exp <- function(GEO, baseDir = dir.data.raw, pattern = NULL, download = T, save = T){
  
  # get expression data with getGEOSuppFiles
  if(download){ # set download = F if the file is already downloaded
    getGEOSuppFiles(GEO, baseDir = baseDir)
  }
  
  # set dir.geo to the GEO data folder
  dir.geo <- file.path(baseDir, GEO, "/")
  
  # list all files in the folder
  exp_file <- list.files(
    dir.geo, 
    pattern = pattern, # add pattern if we want to load a specific file
    recursive = T,
    full.names = T)
  
  # read file
  exp <- data.table::fread(exp_file) %>%
    as.data.frame() %>%
    column_to_rownames("V1") %>%
    as.matrix()
  
  if(save){
    write_csv(
      exp %>% data.frame() %>% rownames_to_column("Gene"),
      file.path(dir.data.processed, paste0(GEO, "_raw_exp_data.csv"))
    )
  }
  
  return(exp)
}