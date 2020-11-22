library(tercen)
library(dplyr)
library(tibble)
library(tidyverse)
library(flowCore)
library(flowClean)
# Minimum amount of cells is 30.000
# QC cutoff value = 10.000


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

data("synPerturbed")
data <- exprs(synPerturbed)
fc_frame <- matrix2flowset(data)
colnames(data)
qc_list <- 
  flowClean::clean(
    fc_frame,
    #filePrefixWithDir = "QC_output",
    vectMarkers = c(5:16),
    #ext = "fcs",
    binSize = 0.01,
    nCellCutoff = 500,
    #announce = TRUE,
    cutoff = "median",
    #diagnostic = TRUE,
    fcMax = 1.3,
    returnVector = TRUE
  )

### vectMarkers is needed for result. 
### V-705 has the outliers.


flag <- ifelse(qc_list >= 10000, "fail", "pass")
qc_df <- data.frame(flag, .ci = (0:(length(flag)-1)))
#qc_result <- ctx$addNamespace(qc_df)
#ctx$save(qc_result)








