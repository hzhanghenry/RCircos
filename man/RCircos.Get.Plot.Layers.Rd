\name{RCircos.Get.Plot.Layers}
\alias{RCircos.Get.Plot.Layers}
\title{
Get the Layer Numbers for Tile Plot
}
\description{
Check out overlaps between different genomic positions on same chromosome and set up layer numbers for each line of the genomic data.
}
\usage{
RCircos.Get.Plot.Layers(genomic.data=NULL, genomic.columns=NULL)
}
\arguments{
  \item{genomic.data}{
A data frame with genomic positions (chromosomes, start and end positions) and the positions should  be already validated and sorted by chromosome then start position.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic positions.
}
}
\value{
A non-negative integer vector with length same as the total rows of input data.
}
\author{
Henry Zhang
}
\examples{
library(RCircos)
data(RCircos.Tile.Data)
layers <- RCircos.Get.Plot.Layers(RCircos.Tile.Data, 3)
}
\keyword{ methods }
