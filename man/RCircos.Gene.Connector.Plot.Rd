\name{RCircos.Gene.Connector.Plot}
\alias{RCircos.Gene.Connector.Plot}
\title{
Draw Connectors between Chromosome Ideogram and Gene Names
}
\description{
Draw a set of connectors (three lines for each connector) between chromosome ideogram or a data track and gene names. RCircos core components and graphics device must be initialized before drawing.
}
\usage{
RCircos.Gene.Connector.Plot(genomic.data=NULL, track.num=NULL, 
        side="in", inside.pos=NULL, outside.pos=NULL, 
        genomic.columns=3, is.sorted=FALSE)
}
\arguments{
  \item{genomic.data}{
A data frame with the first four columns for chromosome names, start positions, end positions, and gene names.
}
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the connectors are plotted.
}
  \item{side}{
Character vector, either "in" or "out" representing the position related to chromosome ideogram.
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
data(RCircos.Gene.Label.Data)
data(UCSC.HG19.Human.CytoBandIdeogram)

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,  
       tracks.inside=10, tracks.outside=0)

gene.data <- RCircos.Gene.Label.Data
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(gene.data, track.num=1, side="in")
RCircos.Gene.Connector.Plot(gene.data, side="out",  
    inside.pos=2.5, outside.pos=3)
}
}
\keyword{ methods }

