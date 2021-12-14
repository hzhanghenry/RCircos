\name{RCircos.ZoomIn.Chromosome}
\alias{RCircos.ZoomIn.Chromosome}
\title{
Zoom in One or Partial Chromosome Ideogram
}
\description{
Increase plot length for one or partial chromosome in order to plot more data points. This works on original chromosome ideogram data and should be done before set RCircos core components. Also, this function only work with one region each time. If there are more than one regions to be zoomed in, repeat this function for each region.
}
\usage{
RCircos.ZoomIn.Chromosome(ideogram=NULL, chromosome=NULL, from=NULL, 
    to=NULL, zoom.in=NULL)
}
\arguments{
  \item{ideogram}{
Data frame, chromosome ideogram data with chromosome names, start position, end position, band name and stains for each band or each chromosome. 
}
  \item{chromosome}{
Character vector, name of chromosome to be zoomed in.
}
  \item{from}{
Non-negative integer, genomic coordinates of start position of chromosome region to be zoomed in.
}
  \item{to}{
Non-negative integer, genomic coordinates of end position of chromosome region to be zoomed in.
}
  \item{zoom.in}{
Non-negative number, fold to zoom in.
}
}
\value{
Data frame, ideogram with new zoomed in region.
}
\author{
Henry Zhang
}
\examples{
\dontrun{ 
library(RCircos) 
data(UCSC.HG19.Human.CytoBandIdeogram) 
zoomedIdeo <- RCircos.ZoomIn.Chromosome( 
    ideogram=UCSC.HG19.Human.CytoBandIdeogram, 
    chromosome="chr17", from=10000, to=20000, 
    zoom.in=1000)
}
}
\keyword{ methods }

