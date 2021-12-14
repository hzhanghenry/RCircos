\name{RCircos.Get.Heatmap.Color.Scale}
\alias{RCircos.Get.Heatmap.Color.Scale}
\title{
Generate Color Scales for Heatmap Plot
}
\description{
Create color scales for heatmap plot. Currently there are six color scales could be generated: \cr\cr
BlueWhiteRed:   colors from blue (lowest) to white then red (highest) \cr
GreenWhiteRed:  colors from green (lowest) to white then red (highest) \cr
GreenYellowRed: colors from green (lowest) to yellow then red (highest) \cr
GreenBlackRed:  colors from green (lowest) to black then red (highest) \cr
YellowToRed:    colors from yellow (lowest) to red (highest) \cr
BlackOnly:      black only
}
\usage{
RCircos.Get.Heatmap.Color.Scale(heatmap.color=NULL)
}
\arguments{
  \item{heatmap.color}{
Character vector, one of "BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed", "GreenBlackRed", "YellowToRed", and "BlackOnly".  
}
}
\value{
A vector or two column matrix of rgb colors.
}
\author{
Henry Zhang
}
\examples{
library(RCircos)
colorScales <- RCircos.Get.Heatmap.Color.Scale(heatmap.color="BlueWhiteRed")
}
\keyword{ methods }

