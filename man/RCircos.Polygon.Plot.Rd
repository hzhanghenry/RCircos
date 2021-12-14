\name{RCircos.Polygon.Plot}
\alias{RCircos.Polygon.Plot}
\title{
Draw Polygons on One Data Track
}
\description{
Plot polygons with different height and different locations inside of a track. Polygon plot is an alternative bar plot that takes both positive and negative height values and genomic intervals with different lengths. Plot data should have genomic positions(chromosome names, start and end positions) as well as height values. Optional column for polygon colors may follow.
}
\usage{
RCircos.Polygon.Plot(polygon.data=NULL, data.col=NULL,  
    track.num=NULL, side=c("in", "out"), border.col=NULL,  
    polygon.col="pink", min.value=NULL, max.value=NULL,  
    inside.pos=NULL, outside.pos=NULL, genomic.columns=3,  
    is.sorted = TRUE)
}
\arguments{
  \item{polygon.data}{
A data frame with three columns for genomic positions, one or columns for polygon heights, and an optional column for polygon colour.
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
   \item{border.col}{
Vector of color names for border color, default null.
}
  \item{polygon.col}{
Color name for fill of polygon.
}
  \item{min.value}{
Numeric, minimum value in data column of polygon data.
}
  \item{max.value}{
Numeric, maximum value in data column of polygon data.
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
    chr.exclude=c("chrX", "chrY"), 
    tracks.inside=10, tracks.outside=5)  
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()

RCircos.Polygon.Plot(RCircos.Polygon.Data,  
    track.num=1, side="in")
}
}
\keyword{ methods }

