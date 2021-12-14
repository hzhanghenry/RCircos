\name{RCircos.Data.Point}
\alias{RCircos.Data.Point}
\title{
Convert A Genomic Position to An Index of Circos Plot Position
}
\description{
Calculate the index of Circos plot positions (a set of x- and y-coordinates for a circular line) for a genomic position (chromosome name and start  position). This function is for internal use only.
}
\usage{
RCircos.Data.Point(chromosome=NULL, start.posotion=NULL)
}
\arguments{
  \item{chromosome}{
Character vector, a chromosome name with prefix of "chr", e.g., "chr1".
}
  \item{start.posotion}{
Non-negative integer, start position of a genomic interval (e.g, a gene or a band) on the chromosome.
}
}
\value{
An integer representing the index of RCircos plot positions.
}
\author{
Henry Zhang
}
\keyword{ methods }
