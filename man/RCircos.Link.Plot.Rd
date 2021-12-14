\name{RCircos.Link.Plot}
\alias{RCircos.Link.Plot}
\title{
Draw Link Lines between Two or More Genomic Positions
}
\description{
Draw lines (quandratic Bezier curves) between paired genomic positions. RCircos core components and graphic device must be initialize before drawing. Link lines are always in the center of plot area.
}
\usage{
RCircos.Link.Plot(link.data=NULL, track.num=NULL,
    by.chromosome=FALSE, start.pos=NULL, genomic.columns=3,
    is.sorted=TRUE, lineWidth=rep(1, nrow(link.data)))
}
\arguments{
  \item{link.data}{
A data frame with paired genomic positions in each row. The two sets of genomic positions must have same number of columns, e.g., chromosome names, start and end positions or only chromosome names and start positions.     
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the link lines will start.
}
  \item{by.chromosome}{
Logic. If true, red color will be used for links between positions of same chromosomes and blue color for links between different chromosomes. Otherwise, user defined or random colors (rainbow) will be automatically assigned for each line.
}
  \item{start.pos}{
Non-negative numeric, scale factor relative to chromosome ideogram position. Must be smaller than 1.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position in each row. Must be either 3 or 2.
}
  \item{is.sorted}{
Logic, whether the data is sorted by chromosome names and start positions.
}  
\item{lineWidth}{
Non-negative integer vector, width for each link line.
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Link.Data);

RCircos.Set.Core.Components(
	cyto.inf=UCSC.HG19.Human.CytoBandIdeogram,  
	chr.exclude=c("chrX", "chrY"), 
	num.inside=10, num.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Link.Plot(RCircos.Link.Data, track.num=4,
	by.chromosome=TRUE)
RCircos.Link.Plot(RCircos.Link.Data, start.pos=0.75,
	by.chromosome=TRUE)
}
}
\keyword{ methods }
