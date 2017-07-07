#' Reading a number of files from a directory into a list of data frames
#'
#' \code{da_Read_Data_Files} reads files from DATA_DIR into a list of 
#' data frames. The full path of the directory should be specificied in 
#' the DATA_DIR argument. For instance, if DATA_DIR = "[YOUR_DATA]", 
#' then the function will read S files from "${HOME}/[YOUR_DATA]/" directory. Files 
#' names are given in the \code{Names} character array of length S. The data frames 
#' created are ready to feed the other package functions. 
#'
## Input:
#' @param DATA_DIR data directory to read from
#' @param Names list containing file names to read
#' @param S number of files to read
#' @param H TRUE/FALSE (with or without column names)
#' @keywords read.table, da_MacKenzie_mle.R
#' @export
#' @examples
#  (1) Importing data to check da_upgma_AIC_cedp:
#' HOME      <- Sys.getenv("HOME")
#' DATA_DIR  <- "/PROJECT_ISLAND_BIOGEOGRAPHY/DATA_ILLES_LACADIVES/"
#' DATA_DIR  <- paste(HOME, DATA_DIR , sep="");
#' Names = c("matrix_Data_156Sp_Header_TRUE.dat");
#' S <- 1
#' Data <- da_Read_Data_Files(DATA_DIR, Names, S, H=TRUE)
#
#  (2) Importing data to check da_MacKenzie_mle.R
#' Names = c("MacKenzie_Data_Matrix_AGT_156Sp_Header_TRUE.dat", 
#' "MacKenzie_Data_Matrix_KAD_156Sp_Header_TRUE.dat", 
#' "MacKenzie_Data_Matrix_KVT_156Sp_Header_TRUE.dat");
#' S <- 3 
#' Data <- da_Read_Data_Files(DATA_DIR, Names, S, H=TRUE)

da_Read_Data_Files <- function( DATA_DIR, Names, S, H=TRUE )
{
  # We can create a list of data frames from each of the files
  # to read:
  DATALIST <- list()
  for (i in 1:S)
  {
    File <- paste(DATA_DIR, Names[i], sep="")
    DATALIST[[i]] <- read.table(File, header=H)
  }
  
  # This is a data list of data frames containing all data, 
  # one member in the list per file read. 
  DATALIST;
}





