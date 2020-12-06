# flowClean operator

##### Description

`flowClean` operator performs a quality control on flowcytometric data.
Number of events has to be >30.000.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. channels, markers)
`col`   | represents the observations (Use Time on top of RowID) 
`y-axis`| measurement value


Output relations|.
---|---
`flag`  |character, quality flag `pass` or `fail`


Input parameters|.
---|---
`binSize`     |  A number in [0,1]; represents the fraction of duration of collection per bin.
`nCellCutoff` |  An integer; represents the minimum number of cells a population must have tobe included in analysis.
`cutoff`      | Method for determining threshold for parameter.  Can be "median" (default) orin [0, 1], which is interpreted as a percentile. Integers > 1 will be interpreted as the fluorescence value to be used for a threshold.
`fcMax`       |  Maximum allowable increase relative to presumed ’good’ data.announceIf TRUE, will print message to screen if errors detected.
`nstable`     |  The number of stable populations required to be observed during the duration ofan experiment. Default is 5.

##### Details

`flowClean` tests  for  deviations  fromuniformity of collection. 
The collection time is discretized into periods and the quality tested. 
Returns a column with "pass" or "fail" per event. 

#### Reference
[flowClean R package]((http://bioconductor.org/packages/release/bioc/html/flowClean.html))

##### See Also

[flowAI R package]((http://bioconductor.org/packages/release/bioc/html/flowAI.html))
, 
[flowsom operator](https://github.com/tercen/flowsom_operator)
