\name{RCircos.Get.Arrow.Shape}
\alias{RCircos.Get.Arrow.Shape}
\title{
Get Default Coordinates for An Arrow Shape
}
\description{
Get default coordinates for an arrow shape. The arrow is represented as a polygon inside of a circle with radius of 1. The cocordinates will be different based which side of chromosome ideogram the arrow will be plotted.
}
\usage{
RCircos.Get.Arrow.Shape(side="in")
}
\arguments{
  \item{side}{
Character vector, either "in" or "out".
}
}
\value{
A two dimensional numeric matrix for x and y coordinates of a polygon.
}
\author{
Henry Zhang
}
\examples{
library(RCircos)
arrow <- RCircos.Get.Arrow.Shape("in")
arrow <- RCircos.Get.Arrow.Shape("out")
}
\keyword{ methods }
