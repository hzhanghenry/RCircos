\name{RCircos.Get.Start.End.Locations}
\alias{RCircos.Get.Start.End.Locations}
\title{
Calculate Start and End Positions for Plot Data
}
\description{
Generate start and end positions from single point for special plot types such as heatmap, histogram, or rectangle-like which need start and end positions along the track. This function is mainly for internal use.
}
\usage{
RCircos.Get.Start.End.Locations(plot.data, plot.width)
}

\arguments{
  \item{plot.data}{
Data frame returned from RCircos.Get.Single.Point.Positions() and the last column is the plot locations
}
  \item{plot.width}{
Non-negative integer in number of base pairs, e.g., width or heatmap or histogram of each data points
}
}
\author{
Henry Zhang
}
\keyword{ methods }
