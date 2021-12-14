\name{RCircos.Pseudo.Ideogram.From.Labels}
\alias{RCircos.Pseudo.Ideogram.From.Labels}
\title{
Generate A Pseudo Chromosome Ideogram from List of Lables
}
\description{
Generate a pseudo ideogram from a character list for RCircos plot Each chromosome will have same length and no stain colors.
}
\usage{
RCircos.Pseudo.Ideogram.From.Labels(chromosomes=NULL)
}
\arguments{
  \item{chromosomes}{
A vector of labels (such as gene names) used for chromosome names and each element must be unique.
}
}
\value{
A data frame as an ideogram data table.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
geneNames <- paste0("Gene_", 1:20)
cyto.info <- RCircos.Pseudo.Ideogram.From.Labels(geneNames)
}
}
\keyword{ methods }
