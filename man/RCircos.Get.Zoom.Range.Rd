\name{RCircos.Get.Zoom.Range}
\alias{RCircos.Get.Zoom.Range}
\title{
Retrieve Genomic Interval from Genomic Data for Zoom Plot
}
\description{
Extract zoom in range (chromosome name, start position, and end position) from zoom data with a gene name or row header. If  using this function to get zoom-in range, there will be no need to call RCircos.Validate.Genomic.Info().
}
\usage{
RCircos.Get.Zoom.Range(zoom.data=NULL, genomic.columns=3)
}
\arguments{
  \item{zoom.data}{
A data frame with genomic positions, gene names, and plot values for genes to be zoomed-in only.
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position in each row. Must be either 3 or 2.
}
}
\value{
A list with chromosome names, start end end positions.
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
zoom.info <- RCircos.Get.Zoom.Range(zoom.data, genomic.columns=3)
}
}
\keyword{ methods }
