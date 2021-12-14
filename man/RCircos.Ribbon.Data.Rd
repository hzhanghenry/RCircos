\name{RCircos.Ribbon.Data}
\alias{RCircos.Ribbon.Data}
\docType{data}
\title{
Sample Data for RCircos Ribbon Plot
}
\description{
A data frame containing paired genomic position data for each row. The contents of each row are in the order of name, start and end positions of chromosome A then name, start and end positions of chromosome B. This data frame is directly used to draw ribbons without extra processing.  
}
\usage{data(RCircos.Ribbon.Data)}
\format{
  A data frame with 4 observations on the following 6 variables.
  \describe{
    \item{\code{chromA}}{a factor with levels \code{chr1} \code{chr5} \code{chr8}}
    \item{\code{chromStartA}}{a numeric vector}
    \item{\code{chromEndA}}{a numeric vector}
    \item{\code{chromB}}{a factor with levels \code{chr10} \code{chr13} \code{chr17} \code{chr18}}
    \item{\code{chromStartB}}{a numeric vector}
    \item{\code{chromEndB}}{a numeric vector}
  }
}
\source{
Unpublished data. 
}
\keyword{datasets}
