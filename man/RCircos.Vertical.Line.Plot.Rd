\name{RCircos.Vertical.Line.Plot}
\alias{RCircos.Vertical.Line.Plot}
\title{
Draw Vertical Lines on One Data Track
}
\description{
Plot vertical lines on a data track for each genomic position  without track 
outlines and sub-track lines. RCircos core components and graphic devices must 
be initialized first.
}
\usage{
RCircos.Vertical.Line.Plot(line.data=NULL, track.num=NULL, 
    side=c("in", "out"), line.width=1, inside.pos=NULL,
    outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
}
\arguments{
  \item{line.data}{
A data frame with leading columns as chromosome names, start position and/or end
 position of each data point followed by data values for one or more samples.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where 
the lines will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to 
chromosome ideogram.
}
  \item{line.width}{
Non-negative integer, width of vertical lines
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of 
the track .
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of 
the track.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position in each row. 
Must be either 3 or 2.
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
data(RCircos.Line.Data);

RCircos.Set.Core.Components(
    cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
    chr.exclude=c("chrX", "chrY"), 
    tracks.inside=10, tracks.outside=5) 
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()

RCircos.Vertical.Line.Plot(RCircos.Line.Data, track.num=1, side="in")
RCircos.Vertical.Line.Plot(RCircos.Line.Data, side="in", 
    inside.pos=1.5, outside.pos=1.75)
}
}
\keyword{ methods }
