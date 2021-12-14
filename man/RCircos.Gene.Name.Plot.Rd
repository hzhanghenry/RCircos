\name{RCircos.Gene.Name.Plot}
\alias{RCircos.Gene.Name.Plot}
\title{
Draw Gene Names on a Data Track
}
\description{
Label gene names along chromosome ideogram or a data track. RCircos core components and graphic device must be initialized before drawing.
}
\usage{
RCircos.Gene.Name.Plot(gene.data=NULL, name.col=NULL,  
    track.num=NULL, side="in", inside.pos=NULL,  
    outside.pos=NULL, genomic.columns=3, is.sorted=FALSE)
}
\arguments{
  \item{gene.data}{
A data frame. The first three columns should be chromosome name, start position, end position.
}
  \item{name.col}{
Non-negative integer, representing the ordinal number of the column in input data that contains gene names.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the gene names are plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
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
data(RCircos.Gene.Label.Data)

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info,  
    chr.exclude=c("chrX", "chrY"),  
    num.inside=10, num.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

gene.data <- RCircos.Gene.Label.Data;
RCircos.Gene.Name.Plot(gene.data, name.col=4,  
    track.num=2, side="in")
RCircos.Gene.Name.Plot(gene.data, name.col=4,  
    track.num=2, side="out",  
    inside.pos=2.5, outside.pos=3)
}
}
\keyword{ methods }
