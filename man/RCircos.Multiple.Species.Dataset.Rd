\name{RCircos.Multiple.Species.Dataset}
\alias{RCircos.Multiple.Species.Dataset}
\title{
Combine Multiple Data Sets form Different Genomes
}
\description{
Combine and modify the chromosome names in multiple species datasets to match the chromosomes in multiple species cytoband data.
}
\usage{
RCircos.Multiple.Species.Dataset(data.list, species)
}
\arguments{
  \item{data.list}{
List of genomic data from multiple species.
}
  \item{species}{
Character vector for prefix of chromosome names to identify different species. The order of species in species list must match the species in data.list.
}
}
\value{
A data frame same as the input data list but with modified chromosome names.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(RCircos.Mouse.Expr.Data)
data(RCircos.Rat.Expr.Data)
dataSets <- list(RCircos.Mouse.Expr.Data, RCircos.Rat.Expr.Data)
dataset <- RCircos.Get.Multiple.Species.Dataset(dataSets,  
    species.list=c("M", "L"))
}
}
\keyword{ methods }
