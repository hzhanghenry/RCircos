\name{RCircos.Get.Heatmap.Data.Colors}
\alias{RCircos.Get.Heatmap.Data.Colors}
\title{
Set Colors for One Track of Heatmap
}
\description{
Assign colors to each value of a numeric vector for heatmap plot. The colour scale needs to be calculated from whole datasets when there are more than one data column so the minimum and maximum values for heatmap data must be provided.
}
\usage{
RCircos.Get.Heatmap.Data.Colors(heatmap.value=NULL, min.value=NULL, max.value=NULL)
}
\arguments{
  \item{heatmap.value}{
A numeric vector.
}
  \item{min.value}{
Numeric, the minimu value for heatmap scale.
}
  \item{max.value}{
Numeric,the maximum value for heatmap scale.
}
}
\value{
Character vector of R color names with length same as the length of data values.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(RCircos.Heatmap.Data)
hColors <- RCircos.Get.Heatmap.Data.Colors(RCircos.Heatmap.Data[,5], min.value=-3, max.value=3)}
}
\keyword{ methods }
