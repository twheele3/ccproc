#
#
#
#
#
#
# BlindImages
#
# An object for storing and processing CellCounter files
#
#
#

blind_files <- function( source_dir, output_dir = NA, title_header = "BLIND", overwrite_names = FALSE ){
  #'@title blind_files
  #'@description Creates blinded versions of all files in a source directory, using a title header and timestamp to create unique filenames. Additionally generates a blinding key alongside output files.
  #'@param source_dir The directory where files to be blinded are located.
  #'@param output_dir (Optional) The directory where blinded files are to be stored. Defaults to source_dir if not specified.
  #'@param title_header (Optional) A descriptive string to give context to blinded files. Defaults to "BLIND" if not specified.
  #'@param overwrite_names (Optional) Overwrites filenames in-place with blinded names. Disabled by default.
  #'@return A system message confirming
  if(is.na(output_dir) | overwrite_names){ output_dir = source_dir }
  filenames <- list.files(path = source_dir)
  filetypes <- unlist(lapply(filenames, FUN=function(x){tail(strsplit(x,"[.]")[[1]],1)}))
  blindkey <- sample(1:length(filenames),length(filenames))
  blindnum <- format(Sys.time(), "%s")
  blindnames <- paste0(title_header,"-",blindnum,"-",blindkey,".", filetypes)
  blindkeydf <- cbind("key" = blindnames, "Original File" = filenames)
  write.csv(blindkeydf, file = paste0(output_dir,"\\",title_header,"-",blindnum," Blinding Key.csv"),row.names=FALSE)
  if(overwrite_names){
    file.rename(from = paste0(source_dir,"\\",filenames), to = paste0(output_dir,"\\",blindnames))
  }else{
    file.copy(from = paste0(source_dir,"\\",filenames), to = paste0(output_dir,"\\",blindnames))
  }
  return(message(paste0("Files in directory ", source_dir, " blinded to ",title_header,"-",blindnum,"-##")))
}
