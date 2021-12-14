\name{RCircos.Customized.Shape.Plot}
\alias{RCircos.Customized.Shape.Plot}
\title{
Plot A Customized Shape
}
\description{
Plot one customized shape on a data track. RCircos core components and graphic device must be initialized first. The customized shape should be represented by coordinates of a polygon inside a circle with radius of 1 and default location is at 12 o'clock.  When plotting, the polygon center will be scaled and transformed for new size and location. For example, following code will plot an arrow:

polygonX <- c(0, -0.7, -0.2, -0.2, 0.2, 0.2, 0.7, 0)
polygonY <- c(-1, 0.7,  0.4,  1,   1,   0.4, 0.7, -1)
plot(c(-2, 2), c(-2, 2)) 
polygon(polygonX, polygonY, col="red")
}
\usage{
RCircos.Customized.Shape.Plot(shape.data=NULL, track.num=NULL,  
    side=c("in", "out"), location=NULL, shape.color="red",
    inside.pos=NULL, outside.pos=NULL)
}
\arguments{
  \item{shape.data}{
A two dimensional numeric matrix for coordinates of a polygon
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the lines will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{location}{
Vector with chromosome name, start position, and end position where the shap to be drawn.
}
  \item{shape.color}{
Character vector, color for the shape.
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)  

data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(UCSC.HG19.Human.CytoBandIdeogram,  
    chr.exclude=NULL, tracks.inside=10, tracks.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

shape <- cbind(c(0, -0.7, -0.2, -0.2, 0.2, 0.2, 0.7, 0),  
        c(-1, 0.7,  0.4,  1,   1,   0.4, 0.7, -1))
RCircos.Customized.Shape.Plot(shape, track.num=1, side="in")
}
}
\keyword{ methods }
