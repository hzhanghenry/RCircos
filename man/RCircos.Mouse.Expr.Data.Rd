\name{RCircos.Mouse.Expr.Data}
\alias{RCircos.Mouse.Expr.Data}
\docType{data}
\title{
Sample Data of Mouse Gene Expression
}
\description{
A data frame containing mouse genomic position information and gene expresio
n values.
}
\usage{data(RCircos.Mouse.Expr.Data)}
\format{
  A data frame with 16499 observations on the following 5 variables.
  \describe{
    \item{\code{Chromosome}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
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
