#' Download and Process GEO Data
#' @description Downloads, decompresses, and loads data from GEO FTP server
#' @param ftp_url URL of the GEO file (must be a .gz file)
#' @param rownames_col Name or index of the column to use as rownames (default: "X")
#' @param cleanup Logical, whether to remove temporary files after processing (default: TRUE)
#' @return Processed data frame
#' @import tidyverse R.utils
#' @export
#' 
download_geo_ftp <- function(ftp_url, rownames_col, cleanup = TRUE) {
  # Validate inputs
  if (!is.character(ftp_url)) {
    stop("ftp_url must be a character string")
  }
  if (!grepl("\\.gz$", ftp_url)) {
    stop("URL must point to a .gz file")
  }
  
  # Create temporary files with tryCatch for proper cleanup
  local_file <- NULL
  decompressed_file <- NULL
  
  tryCatch({
    # Create temporary files
    local_file <- tempfile(fileext = ".gz")
    decompressed_file <- tempfile(fileext = ".csv")
    
    # Download the file with progress indicator
    message("Downloading data from GEO...")
    download_status <- download.file(ftp_url, local_file, mode = "wb", quiet = FALSE)
    
    if (download_status != 0) {
      stop("Failed to download file from GEO")
    }
    
    # Decompress the file
    message("Decompressing file...")
    R.utils::gunzip(local_file, decompressed_file, remove = FALSE)
    
    # Load and process the data
    message("Loading data...")
    data <- read.csv(decompressed_file)
    
    # Set row names if specified
    if (!is.null(rownames_col)) {
      if (!(rownames_col %in% names(data)) && 
          !(is.numeric(rownames_col) && rownames_col <= ncol(data))) {
        stop("Specified rownames column not found in data")
      }
      
      data <- data %>% column_to_rownames(
        if(is.numeric(rownames_col)) names(data)[rownames_col] else rownames_col
      )
    }
    
    # Return the processed data
    message("Data processing complete")
    return(data)
    
  }, error = function(e) {
    stop("Error processing GEO data: ", e$message)
    
  }, finally = {
    # Cleanup temporary files if requested
    if (cleanup && !is.null(local_file) && file.exists(local_file)) {
      file.remove(local_file)
    }
    if (cleanup && !is.null(decompressed_file) && file.exists(decompressed_file)) {
      file.remove(decompressed_file)
    }
  })
}

# Example usage:
if (FALSE) {  # Set to TRUE to run example
  # Download GSE115978 data
  geo_data <- download_geo_data(
    ftp_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978%5Fcounts.csv.gz",
    rownames_col = "X"
  )
  
  # View the results
  head(geo_data)
}

# Function to download GEO dataset
download_geo_direct <- function(gse_id) {
  tryCatch({
    gse <- getGEO(gse_id, GSEMatrix = TRUE)
    return(gse)
  }, error = function(e) {
    stop(paste("Error downloading GEO dataset:", e$message))
  })
}

# Function to extract raw tar data files
extract_tar_files <- function(tar_file) {
  # Create a temporary directory for extraction
  temp_dir <- tempdir()
  
  tryCatch({
    # Extract files
    untar(tar_file, exdir = temp_dir)
    
    # List extracted files, excluding series files
    txt_gz_files <- list.files(
      temp_dir, 
      pattern = "\\.txt\\.gz$", 
      full.names = TRUE
    )
    txt_gz_files <- txt_gz_files[!grepl("series", basename(txt_gz_files), ignore.case = TRUE)]
    
    return(list(temp_dir = temp_dir, files = txt_gz_files))
  }, error = function(e) {
    stop(paste("Error extracting files:", e$message))
  })
}

# Function to process gsm sample files
process_sample_files <- function(files) {
  data_list <- list()
  
  for (file in files) {
    tryCatch({
      # Extract GSM code from filename
      gsm_code <- sub(".*(GSM\\d+).*", "\\1", file)
      
      # Read compressed file
      data <- read.table(gzfile(file), header = TRUE)
      
      # Rename second column with GSM code
      colnames(data)[2] <- gsm_code
      
      data_list[[gsm_code]] <- data
    }, error = function(e) {
      warning(paste("Error processing file", file, ":", e$message))
    })
  }
  
  return(data_list)
}
