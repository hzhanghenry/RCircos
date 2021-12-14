\name{RCircos.Parallel.Line.Plot}
\alias{RCircos.Parallel.Line.Plot}
\title{
Draw Lines between Two Genomic Positions on Same Chromosome
}
\description{
Plot link lines inside a data track between two genomic positions on same chromosome. RCircos core components and graphic device must be initialized first.
}
\usage{
RCircos.Parallel.Line.Plot(line.data=NULL, track.num=NULL, 
    side=c("in", "out"), line.width=1, inside.pos=NULL,
    outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{line.data}{
A data frame containing chromosome, start and end position for each line.
}
  \item{track.num}{
Non-negative integer, the number of track to plot
}
  \item{side}{
Character vector, relative location to ideogram, either "in" or "out".
}
  \item{line.width}{
Non-negative numeric, same as lwd in graphic package
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position in each row. Must be 3.
}
  \item{is.sorted}{
Logic, whether the data is sorted by chromosome names and start positions.
}}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Line.Data);

RCircos.Set.Core.Components(
    cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
    chr.exclude=c("chrX", "chrY"), 
    tracks.inside=10, tracks.outside=5) 
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()

RCircos.Parallel.Line.Plot(line.data, track.num=5, side="in")
RCircos.Parallel.Line.Plot(line.data, line.width=2, 
    inside.pos=2, outside.pos=2.5)
}
}
\keyword{ methods }
