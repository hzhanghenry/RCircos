\name{RCircos.Clear.Track}
\alias{RCircos.Clear.Track}
\title{
Erase One or More Data Tracks
}
\description{
Erase one or more data track or all tracks from current track to the center of plot area. RCircos core components and graphic device must be already initialized.
}
\usage{
RCircos.Clear.Track(track.num=NULL, side=NULL, to.center=FALSE, 
    inside.pos=NULL, outside.pos=NULL)
}
\arguments{
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track to be erased.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram. It could be NULL is both inside.pos and outside.pos are defined.
}
  \item{to.center}{
Logic, FALSE for erasing current track only and TRUE for erasing all tracks from current track to the center of plot plot area. 
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
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram)
data(RCircos.Link.Data)
data(RCircos.Scatter.Data)

RCircos.Set.Core.Components(UCSC.HG19.Human.CytoBandIdeogram,  
    chr.exclude=c("chrX", "chrY"),  
    tracks.inside=10, tracks.outside=0 )

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col=5,  
    track.num=1, side="in", by.fold=1)
RCircos.Link.Plot(RCircos.Link.Data, track.num=3)  

RCircos.Clear.Track(track.num=1, side="in", to.center=FALSE)
RCircos.Clear.Track(track.num=5, side="in", to.center=TRUE) 
}
}
\keyword{ methods }

