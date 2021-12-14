\name{RCircos.Ribbon.Plot}
\alias{RCircos.Ribbon.Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Ribbon between Two Genomic Regions
}
\description{
Draw various ribbons instead of lines with uniform thickness. The thickness of ribbons are relative to the size of the genomic regions to be linked. 
}
\usage{
RCircos.Ribbon.Plot(ribbon.data=NULL, track.num=NULL, 
    by.chromosome=FALSE, twist=FALSE, start.pos=NULL, 
    genomic.columns=3, is.sorted=TRUE)
}
\arguments{
  \item{ribbon.data}{
A data frame with paired genomic positions in each row. The first three columns of each row are for chromosome name, start and end positions of one genomic position followed by three columns for chromosome name, start and end positions of the second genomic position. 
}
  \item{track.num}{
Non negative-integer, representing the ordinal number of the plot track where the ribbons will start.
}
  \item{by.chromosome}{
Logic, If true, red color will be used for ribbons between positions of same chromosomes and blue color for ribbons between different chromosomes. Otherwise, different colors will be automatically assigned for each ribbon.
}
  \item{twist}{
Logic, TRUE for the ribbons with twisted shape and FALSE for no twist.
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
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Link.Data);
data(RCircos.Ribbon.Data);

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
chr.exclude <- c("chrX", "chrY")
RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=0) 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

RCircos.Link.Plot(RCircos.Link.Data, track.num=4, twist=TRUE)
RCircos.Ribbon.Plot(RCircos.Ribbon.Data, 
	track.num=4, twist=TRUE, by.chromosome=FALSE)
}
}
\keyword{ methods }
