\name{RCircos.Get.Zoom.Data}
\alias{RCircos.Get.Zoom.Data}
\title{
Retriev Subset of Plot Data for Zoom
}
\description{
Retrieve subset of plot data based on gene name or row header.
}
\usage{
RCircos.Get.Zoom.Data(plot.data=NULL, name.col= NULL,
 	genomic.columns=3, target.gene=NULL, neighbor.genes=5)
}
\arguments{
  \item{plot.data}{
A data frame containing genomic positions, gene names, and plot values. The data should be already sorted by chromosome names then start positions.
}
  \item{name.col}{
Non-negative integer, which column is for gene names.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position in each row. Must be either 3 or 2.
}
  \item{target.gene}{
Character vector, target gene names.
}
  \item{neighbor.genes}{
Non-negative integer, number of neighbor genes on each side of taget gene to be included in subset.
}
}
\value{
Data frame, subset of input data.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(RCircos.Heatmap.Data)
zoom.data <- RCircos.Get.Zoom.Data(RCircos.Heatmap.Data, 
	name.col=4, target.gene="TP53")
}
}
\keyword{ methods }
