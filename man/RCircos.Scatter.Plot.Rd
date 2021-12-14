\name{RCircos.Scatter.Plot}
\alias{RCircos.Scatter.Plot}
\title{
Scater Plot for One Data Track
}
\description{
Draw one track of scatter plot. RCircos core components and graphics device must be initialized before drawing.
}
\usage{
RCircos.Scatter.Plot(scatter.data=NULL, data.col=4, 
    track.num=NULL, side=c("in", "out"), by.fold=0,
    min.value=NULL, max.value=NULL, inside.pos=NULL, 
    outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
}
\arguments{
  \item{scatter.data}{
A data frame with leading columns of chromosome names, start position and/or end position of each data point followed by data values for one or more samples
}
  \item{data.col}{
Non-negative integer, representing the ordinal number of the column in input dataset that contains the data to be plotted.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the lines will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{by.fold}{
Zero or a positive float number, if it is greater than zero, any data point with a value greater than or equal to the by.fold will be plot as red color and data points with values less or equal to the negative by.fold will be plot as blue colr. Otherwise, black color is used.
}
  \item{min.value}{
Numeric, minimum value of scatter plot data.
}
  \item{max.value}{
Numeric, maximum value of scatter plot data
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
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=0) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data(RCircos.Scatter.Data)
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col=5, 
	track.num=4, side="in", by.fold=1)
}
}
\keyword{ methods }
