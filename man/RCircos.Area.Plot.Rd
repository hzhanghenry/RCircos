\name{RCircos.Area.Plot}
\alias{RCircos.Area.Plot}
\title{
Paint Areas on One Data Track
}
\description{
Paint area  on one data track with different height and different locations inside of a track. Plot types includes "mountain", "curtain", and "band". Plot data should have genomic positions(chromosome names, start and end positions) as well as height values. For band type plot, two columns of height values are required.
}
\usage{
RCircos.Area.Plot(area.data=NULL, data.col=c(4,5),
    track.num=NULL, side=c("in", "out"), 
    plot.type=c("mountain", "curtain", "band"), 
    min.value=NULL, max.value=NULL, area.color="gray",
    border.col="black", inside.pos=NULL, outside.pos=NULL, 
    genomic.columns=3, is.sorted=TRUE)
}
\arguments{
  \item{area.data}{
A data frame with two or three columns for genomic positions, one or more columns for heights of each data point, and an optional column for polygon colour.
}
  \item{data.col}{
Non-negative integer, representing the ordinal number of the column in input data set that contains the data to be plotted.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the lines will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}  
  \item{plot.type}{
Character vector, either "mountain", "curtain", or "band"
}
  \item{min.value}{
Numeric, minimum value in data column of polygon data.
}
  \item{max.value}{
Numeric, maximum value in data column of polygon data.
}
  \item{area.color}{
Color names for fill of the area.
}
  \item{border.col}{
Color name for border color, default null.
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
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Polygon.Data);

RCircos.Set.Core.Components(  
    cyto.info=UCSC.HG19.Human.CytoBandIdeogram,  
    chr.exclude=NULL, tracks.inside=10, tracks.outside=5)  
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()

load("RCircos/data/RCircos.Histogram.Data.RData")

area.data <- RCircos.Histogram.Data;
adj.value <- runif(nrow(area.data), 0, 0.4)
area.data["DataT"] <- 0.5 + adj.value
area.data["DataB"] <- 0.5 - adj.value

RCircos.Area.Plot(area.data, data.col=4, plot.type="mountain", 
    inside.pos=1.2, outside.pos=1.5, is.sorted=FALSE)

RCircos.Area.Plot(area.data, data.col=4, plot.type="curtain", 
    inside.pos=0.9, outside.pos=1.1, is.sorted=FALSE)

RCircos.Area.Plot(area.data, data.col=c(5,6), plot.type="band", 
    inside.pos=0.4, outside.pos=0.7, is.sorted=FALSE)}
}
\keyword{ methods }

