library(tercen)
library(tidyverse)
library(flowCore)
library(flowClean)

data("synPerturbed")

matrix2flowset <- function(a_matrix){ 
  
  minRange <- matrixStats::colMins(a_matrix)
  maxRange <- matrixStats::colMaxs(a_matrix)
  rnge <- maxRange - minRange
  
  df_params <- data.frame(
    name = colnames(a_matrix),
    desc = colnames(a_matrix),
    range = rnge,
    minRange = minRange,
    maxRange = maxRange
  )
  
  params <- Biobase::AnnotatedDataFrame()
  Biobase::pData(params) <- df_params
  Biobase::varMetadata(params) <- data.frame(
    labelDescription = c("Name of Parameter",
                         "Description of Parameter",
                         "Range of Parameter",
                         "Minimum Parameter Value after Transformation",
                         "Maximum Parameter Value after Transformation")
  )
  flowset <- flowCore::flowFrame(a_matrix, params)
  
  return(flowset)
}

data <- exprs(synPerturbed)
fcframe <- matrix2flowset(data)
  
test <- as.matrix(flowClean::clean(fcframe, 
                 #filePrefixWithDir = "QC_output",
                 vectMarkers = c(5:17),
                 ext = "fcs",
                 binSize=0.01,
                 nCellCutoff=500,
                 announce = TRUE,
                 cutoff="median",
                 #diagnostic = TRUE,
                 fcMax=1.3,
                 returnVector = TRUE
))

# Pass cutoff is 10000
ifelse(test>= 10000, "fail" , "pass")
hist(test)


### Checking if the example output is identical to the script output. It is. ###
synPerturbed.c <- as.matrix(clean(synPerturbed, vectMarkers=c(5:17),filePrefixWithDir="sampleName", ext="fcs", returnVector = FALSE))
identical(test, synPerturbed.c)



### Trying to use the flowCore::new function to provide a better flowFrame (without error)


a_matrix <- data

minRange <- matrixStats::colMins(a_matrix)
maxRange <- matrixStats::colMaxs(a_matrix)
rnge <- maxRange - minRange

df_params <- data.frame(
  name = colnames(a_matrix),
  desc = colnames(a_matrix),
  range = rnge,
  minRange = minRange,
  maxRange = maxRange
)

params <- Biobase::AnnotatedDataFrame()
Biobase::pData(params) <- df_params
Biobase::varMetadata(params) <- data.frame(
  labelDescription = c("Name of Parameter",
                       "Description of Parameter",
                       "Range of Parameter",
                       "Minimum Parameter Value after Transformation",
                       "Maximum Parameter Value after Transformation")
)
description <- list("test")

new_flowFrame <- new("flowFrame", exprs= a_matrix, params, description  )

flowClean::clean(new_flowFrame, 
                 #filePrefixWithDir = "QC_output",
                 vectMarkers = c(5:17),
                 ext = "fcs",
                 binSize=0.01,
                 nCellCutoff=500,
                 #announce = TRUE,
                 cutoff="median",
                 #diagnostic = FALSE,
                 fcMax=1.3,
                 returnVector = FALSE
)














