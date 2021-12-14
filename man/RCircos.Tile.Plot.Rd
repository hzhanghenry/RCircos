\name{RCircos.Tile.Plot}
\alias{RCircos.Tile.Plot}
\title{
Plot Tiles for One Data Track
}
\description{
Draw one track of Tiles. RCircos core components and graphics device must be initialized before drawing.
}
\usage{
RCircos.Tile.Plot(tile.data=NULL, track.num=NULL, 
    side= c("in", "out"), inside.pos=NULL, outside.pos=NULL, 
    genomic.columns=3, is.sorted=TRUE)
}
\arguments{
  \item{tile.data}{
A data frame with the three columns as chromosome names, start position and end position of each data point. Other column(s) is optional.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the lines will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track .
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}
  \item{genomic.columns}{
Non-negative interger. Must be 3 for tile plot.
}
  \item{is.sorted}{
Logic, whether the data is sorted by chromosome names and start positions.
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram);

RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=0) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data(RCircos.Tile.Data);
RCircos.Tile.Plot(tile.data=RCircos.Tile.Data, 
	track.num=4, side="in")
}
}
\keyword{ methods }
