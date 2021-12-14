\name{RCircos.Get.Heatmap.Color.Scale.Location}
\alias{RCircos.Get.Heatmap.Color.Scale.Location}
\title{
Calculate X and Y Coordinates for Color Scale
}
\description{
Calculate X and Y coordinates for the topleft corner of heatmap color scale.
}
\usage{
RCircos.Get.Heatmap.Color.Scale.Location(scale.location = 1)
}
\arguments{
  \item{scale.location}{
Integer of 1 ~ 12, represents the plot location (default is 1):\cr
1:  "bottom" \cr
2:  "left" \cr
3:  "top" \cr
4.: "right" \cr
5:  "bottomleft" \cr
6:  "bottomright" \cr
7:  "leftright" \cr
8:  "leftbottom" \cr
9:  "topleft" \cr
10: "topright" \cr
11: "righttop" \cr
12: "rightbottom" \cr
}
}
\value{
Numeric verctor of length 2 for x and y coordinates
}
\author{
Henry Zhang
}
\examples{
library(RCircos)
location <- RCircos.Get.Heatmap.Color.Scale.Location(1)
}
\keyword{ methods }

