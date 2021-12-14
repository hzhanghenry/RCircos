\name{RCircos.Multiple.Species.Core.Components}
\alias{RCircos.Multiple.Species.Core.Components}
\title{
Setup RCircos Core Components with Multiple Species Genomes
}
\description{
Setup RCircos core components with chromosome ideograms of multiple species for visulization of similarity and difference of genomic data between different species.
}
\usage{
RCircos.Multiple.Species.Core.Components(cyto.info.list, species.list, 
    chr.exclude= NULL, tracks.inside, tracks.outside)
}
\arguments{
  \item{cyto.info.list}{
A List contains data frames for chromosome ideogram data of different species. All data frames in the list must have same columns.
}
  \item{species.list}{
Character vector represent multiple species. Each one will be used as prefix of chromosome names for relavant species. The order of species must match the order of the species in cyto.info.list.
}
  \item{chr.exclude}{
Character vector of chromosome names to be excluded from plot, e.g., chr.exclude <- c("chrX", "chrY");
}
  \item{tracks.inside}{
Integer, total number of data tracks inside of chromosome ideogram.
}
  \item{tracks.outside}{
Integer, total number of data tracks outside of chromosome ideogram.
}
}
\author{
Hongen Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram)

cyto.list <- list(UCSC.Mouse.GRCm38.CytoBandIdeogram,  
               UCSC.Baylor.3.4.Rat.cytoBandIdeogram);
species <- c("M", "R")

RCircos.Multiple.Species.Core.Components(  
       cyto.list, species, chr.exclude=NULL,  
       tracks.inside=10, tracks.outside=0)
}
}
\keyword{ methods }
