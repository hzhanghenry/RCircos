\name{RCircos.Pseudo.Ideogram.From.Table}
\alias{RCircos.Pseudo.Ideogram.From.Table}
\title{
Generate Pseudo Chromosome Ideogram from Plot Data  Table
}
\description{
Generate a pseudo ideogram from plot data for Circos plot. Plot data must have at least two columns serving as chromosome names and band names. An optional gene length or locations can also be provided, e.g.,

Groups  LocusName   LocusPosition
Group1  1           1.5
Group1  3           2.0
Group1  4           3.0
}
\usage{
RCircos.Pseudo.Ideogram.From.Table(plot.data=NULL, location.col=NULL, band.col=NULL)
}
\arguments{
  \item{plot.data}{
A data frame with two or three columns, the first column is for chromosome names.
}
  \item{location.col}{
Non negative integer, column number in plot.data from which end positions will be derived for each chromosome or each band. This argument must be defined and could be either length of the locus or relative locus position.
}
  \item{band.col}{
Non negative integer, column number in plot.data which serves as band names, this argument is optional.
}
}
\value{
A data frame with ideogram data.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(RCircos.Gene.Label.Data);
plot.data <- RCircos.Gene.Label.Data;
cyto.info <- RCircos.Get.Pseudo.Ideogram(plot.data=plot.data, 
	location.col=2, band.col=3)
}
}
\keyword{ methods }
