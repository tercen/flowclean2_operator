# flowClean operator

##### Description

`flowClean` operator performs a quality control on flowcytometric data.
Number of events has to be >30.000.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. channels, markers)
`col`   | represents the observations (e.g. cells, samples, individuals) 
`y-axis`| measurement value

Output relations|.
---|---
`flag`         | character, quality flag `pass` or `fail`

Input parameters|.
---|---
`input_var`        | parameter description

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
