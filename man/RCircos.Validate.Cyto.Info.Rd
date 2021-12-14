\name{RCircos.Validate.Cyto.Info}
\alias{RCircos.Validate.Cyto.Info}
\title{
Validate Chromosome Ideogram Information
}
\description{
Validate chromosome ideogram information for correct chromosome order including of correct order of chromosome names, chromosome start, and end positions. The prefix of "chr" in chromosome names is no longer required version 1.2.
}
\usage{
RCircos.Validate.Cyto.Info(cyto.info=NULL, chr.exclude=NULL, is.sorted=TRUE)
}
\arguments{
  \item{cyto.info}{
A data frame with chromosome ideogram data, e.g., an object returned by function of read.table() which read a file containing full information of cytoBandIdeo table from UCSC genome browser.
}
  \item{chr.exclude}{
Character vector, name(s) of chromosomes to be excluded from RCircos plot.
}
  \item{is.sorted}{
Logic, whether the cyto.info has been sorted.
}
}

\value{
A data frame with validated chromosome ideogram table.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- RCircos.Validate.Cyto.Info(
    cyto.info=UCSC.HG38.Human.CytoBandIdeogram,
    chr.exclude=NULL, is.sorted=FALSE)
}
}
\keyword{ methods }
