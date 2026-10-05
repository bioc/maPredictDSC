getLungCancerData <- function(ano = anoLC) {
  
  ## Persistent cache directory
  cache_dir <- tools::R_user_dir("maPredictDSC", "cache")
  data_dir <- file.path(cache_dir, "lungcancer")
  
  if (!dir.exists(data_dir))
    dir.create(data_dir, recursive = TRUE)
  
  ## Get Zenodo record metadata
  api_url <- "https://zenodo.org/api/records/23108090"
  record <- jsonlite::fromJSON(api_url)
  
  ## Files required by ano
  filenames <- as.character(ano$files)
  
  ## Check that all required files are present in the Zenodo record
  missing <- setdiff(filenames, record$files$key)
  
  if (length(missing) > 0) {
    stop(
      "The following CEL files were not found in the Zenodo record: ",
      paste(missing, collapse = ", ")
    )
  }
  
  ## Match required files to their Zenodo download URLs
  ii <- match(filenames, record$files$key)
  urls <- record$files$links$self[ii]
  
  ## Download files that are not already cached
  for (i in seq_along(filenames)) {
    
    destfile <- file.path(data_dir, filenames[i])
    
    if (!file.exists(destfile)) {
      message("Downloading ", filenames[i], " ...")
      
      utils::download.file(
        urls[i],
        destfile = destfile,
        mode = "wb",
        quiet = TRUE
      )
    }
  }
  
  ## Return directory containing the CEL files
  data_dir
}