
## CCPROC Introduction:
'ccproc' is a feature detection package that is designed for use with the Fiji Cell Counter plugin (v3.0.0+) to streamline data processing and auditing for images of structured tissue. The basic idea underlying this package is that technicians can visually inspect and score cells and markers in a structured system using an intuitive and simple plugin, to generate files which can be processed quickly into discrete data for further analysis. This can eliminate a lot of time spent switching back and forth between spreadsheets, and aid in reproducibility and auditing by allowing users to save and load marking files for inspection and correction. This processing also produces data that can be further processed using advanced data parameterization techniques, increasing the threshold of analysis rigor in classically lower-rigor systems without significantly increasing workload.

## Installation
This package can be installed most easily by using the 'devtools' package to compile the 'twheele3/ccproc' package from github.

#### Note: if you run this chunk, you may need to restart R to clear errors that may occur later.

```{r installation}
install.packages("devtools",dependencies = T)
library(devtools)
install_github("twheele3/ccproc", dependencies = T)
```

## Set-up 
After installation, the package only needs to be loaded. Working directory is also set to the Rmd file location for convenience.

```{r setup}
library(ccproc)
knitr::opts_knit$set(root.dir = getwd())
```

## Loading CellCounter files
To begin with, a new CellCounter Database object (ccdb) must be created. This is done at the same time as loading any initial files. Files are loaded here by using list.files to get all files within ./data directory ('.' denoting working directory).

```{r load_files}
ccdb <- CellCounterDB$new( list.files("./data", full.names = TRUE) ) 
```


## Removing duplicates
Duplicates can be optionally removed from files before processing. This is based on the idea that duplicates may be present due to accidental double-clicks, so cells marked within the specified pixel radius 'r' are merged.

```{r remove_duplicates}
ccdb$remove_duplicates(r=5)
```

## Processing data into features
Processing data into features means that cells are algorithmically assigned into features or "shapes" relative to origin positions specified in the CellCounter file. By default it processes all unprocessed CellCounter files into "crypt"-type features, as U-shaped cross-sections.

```{r process_features}
ccdb$processFeatures()
```

## Verify processing with debug plots
To ensure that data is being processed correctly into feature shapes, it may be helpful to create debug plots to visually inspect the processed files. This can be done by specifying which files you want to plot, based on image filenames. If any issues arise, please contact me at twheele3@uoregon.edu, sending debug plot and associated CellCounter file. 

```{r plot_debug}
ccdb$plot_debug(names(ccdb$CCFiles), savedir = "./debug/")
```

## Tagging attributes
Additional cell attributes can be automatically tagged, based on extracting them from the image filename using regex capture groups. Here I am capturing Sample names based on the pattern 'TW##'. To experiment with capturing group formatting, go to https://regexr.com/

```{r tag_attributes}
ccdb$tag_attribute( Regex = "(TW[0-9]{2})", 
                    attribute = "Sample" )
```

## Tagging additional attributes with a key.
Further attributes can be assigned to data in bulk, for the purposes of adding experimental details (such as sex, age, genotype, treatment) as well as unblinding images. This function works by choosing a source column to compare to a list of keys in a specified key_column ("key" by default) in a key database (keydb). keydb may be given as a string pointing to a CSV file, or as a dataframe. 

Additionally, a list of search criteria may be specified as 'whichcells' to constrain the application of these attributes (see formatting rules for which_cells below). This may be useful if a key is a simple alphanumeric value that is only applicable to a single experiment. 
```{r}
ccdb$add_attributes_with_key(key_attribute = "Sample",
                             keydb = "./unblindKey.csv",
                             key_column = "Sample",
                             whichcells = NA)
```

## Finding specific data with search criteria
You may be interested in pulling data with complex criteria to answer certain questions. You can get this by using the built-in 'which_cells' function and specifying a list of criteria to fill. The 'which_cells' function takes a list as input. This list should be formatted so that the names of columns are set equal to the allowed values within those columns. The return is a vector of indices of rows meeting all criteria.

ie, "Position.Crypt.Abs" = 0:10 would return the indices of all cells that are between positions 0 and 10 in crypts.

```{r calling_data}
criteria <- list("Sample" = "TW21", 
                 "Position.Crypt" = -10:10, 
                 "Reg4-dsRed" = TRUE,
                 "Lrig1" = c(TRUE,FALSE))
index <- ccdb$which_cells(criteria)
print(ccdb$cells[index,])

```

## Aggregating data by features
Data can be aggregated by specified features into count and ratio for specified columns. For instance below, counts of all markers are independently taken per-crypt per-sample. The measure.by variable allows for separate markers to be specified to aggregate them in a single command. 

If a vector of more than one marker is specified in a list entry, then this command will aggregate coexpression as well. This is compiled under MarkerCombo, which is aggregated by Marker.Expression. For instance, '1.1' would indicate that Marker 1 has expression of 1. '2.0' would indicate that Marker 2 has expression of 0. Thus in my example below, I only include MarkerCombo=="1.1" as I'm only interested in positive counts. 

```{r aggregating_data}

per.crypt.df <- ccdb$aggregate_by_features(
    group.by = c("Sample","Crypt"), 
    measure.by = as.list(ccdb$Markers),
    melt.data = TRUE)

# Only keeps rows for where Ki67 is present, as indicated by measure.by marker 1 being 1 ("1.1")
per.crypt.df <- per.crypt.df[which(per.crypt.df$MarkerCombo=="1.1"),]

print(per.crypt.df)
```

## Exporting data to CSV
Dataframes can be exported to CSV to be used in other programs like Excel and GraphPad. 
```{r}
write.csv( ccdb$cells , file ="./CCDBcells.csv", row.names = FALSE )
write.csv( per.crypt.df , file ="./CCDBpercrypt.csv", row.names = FALSE )
```

## Saving and loading database
The database object can be saved and loaded for later use. 

```{r save_and_read}
ccdb$save(file = "./ccdb.rds")
ccdb <- readRDS(file = "./ccdb.rds")
```

## Adding and removing files
Files can be removed and added later if needed. Processing will have to be done again in this case.

```{r removing_and_adding}
ccdb$removeCC("IHC009 TW30-2L Bmi1 001.lsm")
ccdb$addCC("./data/CellCounter_IHC009 TW30-2L Bmi1 001.xml")
ccdb$processFeatures("IHC009 TW30-2L Bmi1 001.lsm")
```

