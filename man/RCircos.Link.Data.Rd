\name{RCircos.Link.Data}
\alias{RCircos.Link.Data}
\docType{data}
\title{
Sample Data for RCircos Link Plot 
}
\description{
A data frame containing paired genomic position data for each row. The content of each row must be in the order of name, start and end position of chromosome A then name, start and end position of chromosome B. This data frame is directly used to draw link lines without extra processing.  
}
\usage{data(RCircos.Link.Data)}
\format{
  A data frame with 71 observations on the following 6 variables.
  \describe{
    \item{\code{Chromosome}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr14} \code{chr15} \code{chr17} \code{chr19} \code{chr2} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX}}
    \item{\code{chromStart}}{a numeric vector}
    \item{\code{chromEnd}}{a numeric vector}
    \item{\code{Chromosome.1}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX}}
    \item{\code{chromStart.1}}{a numeric vector}
    \item{\code{chromEnd.1}}{a numeric vector}
  }
}
\source{
Unpublished data. 
}
\keyword{datasets}

