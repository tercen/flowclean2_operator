library(tercen)
library(dplyr)
library(tidyverse)
library(flowCore)
library(flowClean)

# Minimum amount of cells is 30.000
# QC cutoff value = 10.000
# markercolumns have to be selected

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

ctx <- tercenCtx(workflowId = "affe1383f9318d3b8141a4bfa700a359",
                 stepId = "526578e0-b152-4db0-ba75-0f32ee999a84")

time <- ctx$cselect(ctx$cnames[[1]])
data <- ctx$as.matrix() %>% t()
data <- as.matrix(cbind(data, time))

# Indicate which columns are the markers of which you want to QC. (take out FSC SSC etc.)
# FSC and SSC are in 1, 2. Time is in the last column.
markercolumns <- c(3: (ncol(data)-1))

fc_frame <- matrix2flowset(data)

qc_df <- as.data.frame(flowClean::clean(fc_frame, 
                                        #filePrefixWithDir = "QC_output",
                                        vectMarkers = markercolumns,
                                        #ext = "fcs",
                                        binSize=0.01,
                                        nCellCutoff=500,
                                        #announce = TRUE,
                                        cutoff="median",
                                        #diagnostic = TRUE,
                                        fcMax=1.3,
                                        returnVector = TRUE
                                        ))
### returnVector has to be TRUE, otherwise errormessage:
### Error in makeFCS(fF, GoodVsBad, filePrefixWithDir, numbins, nCellCutoff,  : 
### (list) object cannot be coerced to type 'double'

flag <- ifelse(qc_df >= 10000, "fail", "pass")
qc_result <- ctx$addNamespace(flag)

