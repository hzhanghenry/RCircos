\name{RCircos.Heatmap.Data}
\alias{RCircos.Heatmap.Data}
\docType{data}
\title{
Sample Data for RCircos Heatmap Plot
}
\description{
A dataset with the first three columns as chromosome names, start position, and end position followed by a column of gene names then columns of gene expression values.
}
\usage{data(RCircos.Heatmap.Data)}
\format{
  A data frame with 6660 observations on the following 10 variables.
  \describe{
    \item{\code{Chromosome}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
    \item{\code{chromStart}}{a numeric vector}
    \item{\code{chromEnd}}{a numeric vector}
    \item{\code{GeneName}}{a character vector}
    \item{\code{X786.O}}{a numeric vector}
    \item{\code{A498}}{a numeric vector}
    \item{\code{A549.ATCC}}{a numeric vector}
    \item{\code{ACHN}}{a numeric vector}
    \item{\code{BT.549}}{a numeric vector}
    \item{\code{CAKI.1}}{a numeric vector}
  }
}
\source{
Modified from NCBI GEO dataset: GSE32474.
}
\keyword{datasets}
