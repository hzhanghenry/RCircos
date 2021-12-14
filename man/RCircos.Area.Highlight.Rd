\name{RCircos.Area.Highlight}
\alias{RCircos.Area.Highlight}
\title{
Highlight PLot Area Cross One or More Tracks
}
\description{
Highlight a plot area with transparent color cross one or more tracks. RCircos core components and graphic device must be initialized first.
}
\usage{
RCircos.Area.Highlight(highlight.area=NULL, track.num=NULL,  
    side=c("in", "out"), hightlight.color=rgb(0.5, 0.5, 0, 0.5),  
    inside.pos=NULL, outside.pos=NULL)
}
\arguments{
  \item{highlight.area}{
Vector with chromosome name, start and end position to be highlighted.
}
  \item{track.num}{
Vector of non-negative integer, which track or track(s) to be highlighted.
}
  \item{side}{
Character vector, location relative to the ideogram, either "in" or "out".
}
  \item{hightlight.color}{
An RGB color definition, alpha value must be defined.
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) 
of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) 
of the track.
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
data(RCircos.Gene.Label.Data)

RCircos.Set.Core.Components(UCSC.HG19.Human.CytoBandIdeogram,  
    chr.exclude=NULL, tracks.inside=10, tracks.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

highlight.area <- c("chr1", 100000, 200000)
RCircos.Area.Highlight(highlight.area=c("chr1", 100000, 200000),  
    track.num=c(1:3), side="in")
RCircos.Area.Highlight(highlight.areac("chr10", 100000, 200000),  
    inside.pos=1.5, outside.pos=2)  
}
}
\keyword{ method }

