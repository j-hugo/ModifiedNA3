# Modified Neural Activity Cubic

## Original Publication

This repository contains four R script files that implement functionality for activity detection from calcium imaging recordings originally published in Prada, J. et al. (2018). An open source tool for automatic spatiotemporal assessment of calcium transients and local ‘signal-close-to-noise’ activity in calcium imaging data. PLoS Comput Biol, 14(3), e1006054.

## Modifications to the original code

The modifications to the original source code merely enhance usability. They *do not* influence the analysis results provided by Neural Activity Cubic. 

The following functionality was added:

* Allow CSV and AVI files as pipeline input
* Return analysis results as CSV file
* Remove dependence on Bio7 environment and Java code

## Usage

* Open `NACmodROI.R` and adjust directories and settings in the file.
* Source `NACmodROI.R`.

## Disclaimer

The provided scripts by no means constitute a functional R program/package. The scripts were merely used for a specific use case.