\name{RCircos.Get.Gene.Label.Locations}
\alias{RCircos.Get.Gene.Label.Locations}
\title{
Calculate Plot Positions for Gene Labels
}
\description{
This function calculates new plot positions from genomic positions for gene labeling. In case there are too many genes in a genomic interval, the lables may become overlapped so that correct numbers and new plot positions are often needed to make the labels readable. If all labels cannot be plotted for the given genomic interval extra labels will be ignored. This function is for internal use only.
}
\usage{
RCircos.Get.Gene.Label.Locations(genomic.data=NULL, 
	genomic.columns=3, is.sorted = TRUE)
}
\arguments{
  \item{genomic.data}{
A data frame with leading columns for chromosome name, start position, and/or end position.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position in each row. Must be either 3 or 2.
}
  \item{is.sorted}{
Logic, whether the data is sorted by chromosome names and start positions.
}
}
\value{
All or subset of input data frame with a new column for plot positions.
}
\author{
Henry Zhang
}
\keyword{ methods }
