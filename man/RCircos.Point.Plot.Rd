\name{RCircos.Point.Plot}
\alias{RCircos.Point.Plot}
\title{
Point Plot for One Data Track
}
\description{
Plot points on a track without track outlines and sub-track lines. The size and location od each point can be adjusted with the data values. RCircos core components and graphic device must be initialized first.
}
\usage{
RCircos.Point.Plot(point.data=NULL, data.col=4, 
    track.num=NULL, side=c("in", "out"), 
    min.value=NULL, max.value=NULL, 
    point.type=19, with.height=TRUE, 
    with.size=FALSE, point.scale=1, 
    inside.pos=NULL, outside.pos=NULL, 
    genomic.columns=3,  is.sorted=TRUE)
}
\arguments{
  \item{point.data}{
A data frame with leading columns of chromosome names, start position and/or end position of each data point followed by data values for one or more samples
}
  \item{data.col}{
Non-negative integer, representing the ordinal number of the column in input dataset that contains the data to be plotted.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the lines will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{min.value}{
Numeric, minimum value of point data.
}
  \item{max.value}{
Numeric, maximum value of point data
}  
\item{point.type}{
Non-negative integer for pch. Default is 19.
}
  \item{with.height}{
Logic, if TURE, the location of each point will be adjusted according to the data value.
}
  \item{with.size}{
Logic, if TURE, the size of each point will be adjusted according to the data value.
}
  \item{point.scale}{
Non-negative numeric, more scale for point size, must be greater than or equal to 1.
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position in each row. Must be either 3 or 2.
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
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)

RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram,
	chr.exclude=c("chrX", "chrY"), 
	num.inside=10, num.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data(RCircos.Scatter.Data)
RCircos.Scatter.Plot(RCircos.Scatter.Data, 
	data.col=5, track.num=4, side="in", 
	with.height=TRUE, with.size=FALSE)
}
}
\keyword{ methods }

