\name{RCircos.Get.Single.Point.Positions}
\alias{RCircos.Get.Single.Point.Positions}
\title{
Calculate Plot Positions for Genomic Data
}
\description{
Calculate x and y coordinates of one genomic position for each row of a genomic data set. As a set of points for a circlar line is held in the RCircos session, we only need the index of the point for each genomic position.
}
\usage{
RCircos.Get.Single.Point.Positions(genomic.data=NULL, genomic.columns=3)
}
\arguments{
  \item{genomic.data}{
A data frame contains genomic positions (at least two or three columns for chromosome names, start and/or end positions). 
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position (chromosome name, start and/or end position).
}
}
\value{
A data frame same as input but with a new column for index of plot positions on a circular line. 
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(RCircos.Heatmap.Data)
dat <- RCircos.Get.Single.Point.Positions(RCircos.Heatmap.Data)
}
}
\keyword{ methods }

