\name{RCircos.Get.Chromosome.Order}
\alias{RCircos.Get.Chromosome.Order}
\title{
Find the order of Chromosome Names
}
\description{
Generate an ordered chromosome names from input. For human and other mammalian animals, numeric names (integers or Roman numbers) will go first followed by chromosome X, Y, and M. If all chromosome names are character only, they will be in alphabetical order. 
}
\usage{
RCircos.Get.Chromosome.Order(chromosomes=NULL)
}
\arguments{
  \item{chromosomes}{
Character vector, names of chromosomes.
}
}
\value{
Character vector, ordered chromosome name list.
}
\author{
Henry Zhang
}
\examples{
library(RCircos)
chromosomes <- paste0("chr", c(1:10, 21,22, 11:20, "X", "Y"))  
chromosomes <- RCircos.Get.Chromosome.Order(chromosomes)  
}
\keyword{ methods }

