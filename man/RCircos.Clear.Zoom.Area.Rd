\name{RCircos.Clear.Zoom.Area}
\alias{RCircos.Clear.Zoom.Area}
\title{
Erase Zoom-in Plot Area
}
\description{
Erase zoom-in plot area. RCircos core components and graphic device must be already initialized.
}
\usage{
RCircos.Clear.Zoom.Area(zoom.pos=NULL, track.num=NULL,  
        inside.pos=NULL, outside.pos=NULL)
}
\arguments{
  \item{zoom.pos}{
Non-negative integer vector, representing the index of RCircos base plot positions.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track to be erased.
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
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,  
    tracks.inside=10, tracks.outside=0 )

RCircos.Set.Plot.Area()
RCircos.Clear.Zoom.Area(zoom.pos=c(1:10000), inside.pos=2, outside.pos=2.5) 
}
}
\keyword{ methods }

