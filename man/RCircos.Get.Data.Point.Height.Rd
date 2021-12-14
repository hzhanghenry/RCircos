\name{RCircos.Get.Data.Point.Height}
\alias{RCircos.Get.Data.Point.Height}
\title{
Calculate the Height of Data in A Track
}
\description{
Calculate data point height inside a plot track such as scatter location, top or bottom location of a bar, layer of a tile or parallel link line.
}
\usage{
RCircos.Get.Data.Point.Height(plot.values=NULL, min.value=NULL, 
    max.value=NULL, plot.type=NULL, track.height=NULL)
}
\arguments{
  \item{plot.values}{
Numeric, the data to be plotted on a data track.
}
  \item{min.value}{
Numeric, the minimum value of data range.
}
  \item{max.value}{
Numeric, the maximum value of data range.
}
  \item{plot.type}{
Character vector, plot type, valid values are "bar", "histogram", "uniform", or "points".
}
  \item{track.height}{
Non-negative numeric, height of plot track.
}
}
\value{
Numeric vector with values between 0 ~ 1
}
\author{
Henry Zhang
}
\examples{
data.values <- runif(1000, -4, 11)
library(RCircos)
data.height <- RCircos.Get.Data.Point.Height(data.values,  
    min.value=-4, max.value=14, plot.type="points",  
    track.height=NULL)
}
\keyword{ methods }

