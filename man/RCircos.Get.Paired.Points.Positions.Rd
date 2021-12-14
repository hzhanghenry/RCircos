\name{RCircos.Get.Paired.Points.Positions}
\alias{RCircos.Get.Paired.Points.Positions}
\title{
Calcualte Plot Positions for A Genomic Interval
}
\description{
Calculate x and y coordinates for a genomic interval which needs two plot positions (start and end) such as link lines/ribbons, tiles, parallel lines. As a set of points for a circlar line is held in the RCircos session, we only need the index of the point for each genomic position.
}
\usage{
RCircos.Get.Paired.Points.Positions(genomic.data=NULL, genomic.columns=3,  
    plot.type=c("link", "ribbon", "pLink", "polygon", "tile"))
}
\arguments{
  \item{genomic.data}{
A data frame contains paired genomic positions (chromosome names, start and end positions). The data does not need to be sorted.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position (chromosome name, start and/or end position).
}
  \item{plot.type}{
Chraracter vector, either "link", "ribbon", "pLink", "polygon", or "tile".
}
}
\value{
A data frame same as input but with two new columns for index of plot positions on circular line. 
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
data(RCircos.Link.Data)

RCircos.Set.Core.Components(UCSC.HG19.Human.CytoBandIdeogram,  
    chr.exclude=NULL, track.inside=10, tracks.outside=0)
linkData<-RCircos.Get.Paired.Points.Positions(RCircos.Link.Data,
    genomic.columns=3, plot.type="link")
}
}
\keyword{ methods }
