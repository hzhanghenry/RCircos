\name{RCircos.Validate.Genomic.Info}
\alias{RCircos.Validate.Genomic.Info}
\title{
Validate Genomic Information
}
\description{
Check out if the chromosome name, start and end positions in the range of ideogram data.
}
\usage{
RCircos.Validate.Genomic.Info(genomic.info=NULL)
}
\arguments{
  \item{genomic.info}{
A vector with a chromosome name, start and end position on the chromosome
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)

RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram, 
    chr.exclude=NULL, tracks.inside=10, tracks.outside=0 )
RCircos.Validate.Genomic.Info(c("chr1", 10000, 500000))
}
}
\keyword{ methods }
