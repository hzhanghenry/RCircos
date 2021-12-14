\name{RCircos.Track.Outline}
\alias{RCircos.Track.Outline}
\title{
Draw Outline for A Data Track
}
\description{
Draw outline with subtracks for one data track. RCircos core components and graphic device must be initialized first. This function is mainly for internal use.
}
\usage{
RCircos.Track.Outline(inside.pos=NULL, outside.pos=NULL,
    num.layers=1, chrom.list=NULL, track.colors=NULL)
}
\arguments{
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}  \item{num.layers}{
Non-negative integer, number of sub-tracks lines. 0 for no sub-track line.
}
  \item{chrom.list}{
List of chromosome names for which outline will be drawn.
}
  \item{track.colors}{
Vector of color names for outline of each chromosome.
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
	tracks.inside=10, tracks.outside=5) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

locations <- RCircos.Track.Positions(side="in", track.num=1);  
out.pos <- locations[1];
in.pos  <- locations[2];
subtracks <- 5;
RCircos.Track.Outline(inside.pos=locations[2], 
	outside.pos=locations[1], num.layers=5,
	chrom.list=NULL, track.colors=NULL) 
}
}
\keyword{ methods }

