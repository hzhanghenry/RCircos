\name{RCircos-package}
\alias{RCircos-package}
\alias{RCircos}
\docType{package}
\title{
Circos-like Plot of Genomic Data
}
\description{
This is the R implement of basic Circos plots of genimic data.
}
\details{
\tabular{ll}{
Package: \tab RCircos\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2016-8-01\cr
License: \tab GPL (>=2) \cr
}
RCircos package provides a simple and flexible way to make Circos 2D track plots with R and could be easily integrated into other R data processing and graphic manipulation pipelines for presenting large-scale multi-sample genomic research data. It can also serve as a base tool to generate complex Circos images. Currently, RCircos supports genomic data plot for three species: human, mouse, and rat. Data plots include chromosome ideogram with cytoband, heatmap, histogram, lines, scatterplot, and tiles. Plot items for further decorations have connectors, links (lines and ribbons), and text (gene) labels.  Each plot is implemented with a specific function and input data for all functions are data frames which can be objects read from text files or generated with other R pipelines.
}
\author{
Henry H Zhang

Maintainer: Henry H Zhang <hzhang.henry@gmail.com>
}
\references{
Hongen Zhang, Paul Meltzer, and Sean Davis. RCircos: an R package for Circos 2D track plots. BMC Bioinformatics, 2013, 14:244.
}
\keyword{ package }
