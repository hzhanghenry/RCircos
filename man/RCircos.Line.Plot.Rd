\name{RCircos.Line.Plot}
\alias{RCircos.Line.Plot}
\title{
Draw Lines for One Data Track
}
\description{
Draw one track of continue lines which connect each neighbor data points. RCircos core components and graphics device must bae initialized before drawing.
}
\usage{
RCircos.Line.Plot(line.data=NULL, data.col=4, track.num=NULL,  
    side=c("in", "out"), min.value=NULL, max.value=NULL,  
    inside.pos=NULL, outside.pos=NULL, genomic.columns=3,  
    is.sorted=TRUE)
}
\arguments{
  \item{line.data}{
A data frame with the first three columns for chromosome names, start and end position of each data point followed by data values for one or more samples.
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
  \item{min.value}{
Numeric, minimum value of line data.
}
  \item{max.value}{
Numeric, maximum value of line data
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
data(RCircos.Line.Data);

RCircos.Set.Core.Components(  
    cyto.info=UCSC.HG19.Human.CytoBandIdeogram,  
    chr.exclude=c("chrX", "chrY"),  
    num.inside=10, num.outside=0) 

RCircos.Set.Plot.Area() 
RCircos.Chromosome.Ideogram.Plot()  

RCircos.Line.Plot(RCircos.Line.Data, data.col=5,  
    track.num=1, side="in")
RCircos.Line.Plot(RCircos.Line.Data, data.col=5,  
    side="in", inside.pos=1.5, outside.pos=1.75)  
}
}
\keyword{ methods }
