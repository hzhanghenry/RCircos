\name{RCircos.Rat.Expr.Data}
\alias{RCircos.Rat.Expr.Data}
\docType{data}
\title{
Sample Data of Rat Gene Expression
}
\description{
A data frame containing rat genes,  genomic position and expression values of these genes.
}
\usage{data(RCircos.Rat.Expr.Data)}
\format{
  A data frame with 11426 observations on the following 5 variables.
  \describe{
    \item{\code{Chromosome}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX}}
    \item{\code{txStart}}{a numeric vector}
    \item{\code{txEnd}}{a numeric vector}
    \item{\code{Gene}}{a character vector}
    \item{\code{Expr.Mean}}{a numeric vector}
  }
}
\source{
Modified from GEO dataset: GSE42081.
}
\keyword{datasets}
