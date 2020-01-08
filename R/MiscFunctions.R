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

blind_files <- function( source_dir, output_dir = NA, key_dir = NA, title_header = "BLIND", overwrite_names = FALSE ){
  #'@title blind_files
  #'@description Creates blinded versions of all files in a source directory, using a title header and timestamp to create unique filenames. Additionally generates a blinding key alongside output files.
  #'@param source_dir The directory where files to be blinded are located.
  #'@param output_dir (Optional) The directory where blinded files are to be stored. Defaults to source_dir if not specified.
  #'@param title_header (Optional) A descriptive string to give context to blinded files. Defaults to "BLIND" if not specified.
  #'@param overwrite_names (Optional) Overwrites filenames in-place with blinded names. Disabled by default.
  #'@return A system message confirming
  if(is.na(output_dir) | overwrite_names){ output_dir = source_dir }
  if(is.na(key_dir)){ key_dir = output_dir }
  if(is.character(output_dir)){
    if(!dir.exists(output_dir)){
      tryCatch({dir.create(output_dir)},
               error = function(err){message(paste0("Could not create directory: ",output_dir,"/n",err,"/n"))})
    }
  }
  if(is.character(key_dir)){
    if(!dir.exists(key_dir)){
      tryCatch({dir.create(key_dir)},
               error = function(err){message(paste0("Could not create directory: ",key_dir,"/n",err,"/n"))})
    }
  }
  filenames <- list.files(path = source_dir)
  filetypes <- unlist(lapply(filenames, FUN=function(x){tail(strsplit(x,"[.]")[[1]],1)}))
  blindkey <- sample(1:length(filenames),length(filenames))
  blindnum <- format(Sys.time(), "%s")
  blindnames <- paste0(title_header,"-",blindnum,"-",blindkey,".", filetypes)
  blindkeydf <- cbind("key" = blindnames, "Original File" = filenames)
  write.csv(blindkeydf, file = paste0(key_dir,"\\",title_header,"-",blindnum," Blinding Key.csv"),row.names=FALSE)
  if(overwrite_names){
    file.rename(from = paste0(source_dir,"\\",filenames), to = paste0(output_dir,"\\",blindnames))
  }else{
    file.copy(from = paste0(source_dir,"\\",filenames), to = paste0(output_dir,"\\",blindnames))
  }
  return(message(paste0("Files in \"", source_dir, "\" have been blinded to ",title_header,"-",blindnum,"-## in ",c(paste0("\"",output_dir,"\"."), "place.")[1+overwrite_names],
                        "\n", "Blinding key saved to: \"",paste0(key_dir,"/",title_header,"-",blindnum," Blinding Key.csv"),"\"")))
}

add_columns_with_key = function(df, key_column, keydb, keydb_column = "key", which_rows = NA){
  # A function to add columns in bulk with a key table.
  # key_column : a string specifying a source column from df table
  # keydb : a table containing attribute info, with column titles. May be a string pointing to a csv file.
  # keydb_column : the column name within keydb that contains key information.
  # whichcells : a list of parameters to constrain attribute assignment to, as per self$which_cells function.
  # is.new : If TRUE, will add a metadata entry to parse future entries automatically.
  #TODO : integrate with add_df

  # Input error checking
  if(is.character(keydb)){
    tryCatch(
      keydb <- read.csv(keydb,header = T),
      error = function(err){message(paste0("Key database could not be parsed as data frame: /n",err,"/n"))})
  }

  if(!is.data.frame(keydb)){
    return(message("Key database could not be parsed as data frame."))
  }

  if(!any(colnames(keydb) == key_column)){
    return(message(paste0("No columns in key dataframe specified as '", key_column,"'.")))
  }

  # Experiment_regex = NA by default because may fall in key_regex
  if(is.na(whichcells)){
    row_indices <-  1:nrow(df)
  }else{
    row_indices <-  which_rows(df,which_rows)
  }

  # For future use, extract all other column names aside from key.
  valnames <- colnames(keydb)[!(colnames(keydb) %in% keydb_column)]

  # Coerce all key columns into factors.
  for(i in names(keydb)) { keydb[,i] <- factor(keydb[,i]) }

  # Add new factor levels to new columns
  for(i in valnames) {
    # Add column for new attribute types
    if(!(i %in% colnames(df))) { df[,i] <- NA }
    # Coerce target column to factor if not already
    if(class(df[,i]) != "factor") { df[,i] <- factor(df[,i]) }
    #Add new factor levels to column.
    levels(df[,i]) <- c( levels(df[,i]), levels(keydb[,i]) )
  }

  for(i in 1:nrow(keydb)){
    key <- keydb[,keydb_column][i]
    key_indices <- which(df[,key_column] %in% key)
    indices <- key_indices[which(key_indices %in% row_indices)]
    df[indices,valnames] = keydb[i,valnames]
  }


}

which_rows = function(df, listconditions){
  # Returns indices in self$cells that match conditions specified in list.
  if(!all(names(listconditions) %in% colnames(df))){
    return(message(paste0("Error in which_rows: Columns not found in table: ",paste0(names(listconditions)[!(names(listconditions) %in% colnames(df))],collapse = ", "))))
  }
  return(which(apply(sapply(names(listconditions),function(x){df[,x] %in% listconditions[[x]] }),1,prod)>0))
}


