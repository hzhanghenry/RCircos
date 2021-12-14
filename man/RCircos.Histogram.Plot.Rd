\name{RCircos.Histogram.Plot}
\alias{RCircos.Histogram.Plot}
\title{
Plot Histogram for One Data Track
}
\description{
Draw one track of histogram. RCircos core components and graphic device must be initialized before drawing.
}
\usage{
RCircos.Histogram.Plot(hist.data=NULL, data.col=4,
    track.num=NULL, side=c("in", "out"), 
    min.value=NULL, max.value=NULL, inside.pos=NULL, 
    outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
}
\arguments{
  \item{hist.data}{
A data frame with the first three columns as chromosome names, start position and end position of each data point followed by data values for one or more samples.
}
  \item{data.col}{
Non-negative integer, representing the ordinal number of the column in input data that contains the data to be plotted.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the histogram will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{min.value}{
Numeric, minimum value for histogram
}
  \item{max.value}{
Numeric, maximum value for histogram
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track .
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
data(RCircos.Histogram.Data);

RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=5) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

RCircos.Histogram.Plot(RCircos.Histogram.Data, 
	data.col=4, track.num=1, side="in")
RCircos.Histogram.Plot(RCircos.Histogram.Data, 
	data.col=4, side="in", inside.pos=1.7, 
	outside.pos=1.8)
}
}
\keyword{ methods }
