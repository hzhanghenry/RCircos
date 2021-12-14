\name{RCircos.Zoom.Area.Outline}
\alias{RCircos.Zoom.Area.Outline}
\title{
Draw Outline for Zoom Plot Area
}
\description{
Draw outline and  subtrack lines for zoomed-in area of a track in outside of chromosome ideogram. RCircos core components must be initialized first. This function is mainly for internal use.
}
\usage{
RCircos.Zoom.Area.Outline(zoom.pos=NULL, inside.pos=NULL,  
    outside.pos=NULL, num.layers=5, fill.col="white")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zoom.pos}{
Non-negative numeric vector, the index of RCircos plot position.
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}
  \item{num.layers}{
Non-negative integer, number of sub-tracks in a data track.
}
  \item{fill.col}{
Character vector for color name to fill the polygon.
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram)

RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=5) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data(RCircos.Heatmap.Data)
zoom.data <- RCircos.Get.Zoom.Data(RCircos.Heatmap.Data, 
			name.col=4, genomic.columns=3,
			target.gene="SP5", neighbor.genes=5)
zoom.range <- RCircos.Get.Zoom.Range(zoom.data, 3)
zoom.pos <- RCircos.Set.Zoom.Plot.Positions(zoom.range, 
        total.genes=11, area.length=0.25, gene.width=NULL)

RCircos.Zoom.Area.Outline(zoom.pos=NULL, inside.pos=2.5,
	outside.pos=3, num.layers=5, fill.col="white")
}
}
\keyword{ methods }
