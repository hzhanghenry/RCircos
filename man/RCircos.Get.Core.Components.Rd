\name{RCircos.Get.Core.Components}
\alias{RCircos.Get.Plot.Parameters}
\alias{RCircos.Get.Plot.Positions}
\alias{RCircos.Get.Plot.Ideogram}
\title{
Methods to Retrieve RCircos Core Components Used in Current RCircos Session
}
\description{
Retrieve base plot positions, plot ideogram, and plot parameters used in current RCircos session. RCircos core components must be initialized first. 
}
\usage{
RCircos.Get.Plot.Parameters()
RCircos.Get.Plot.Positions()
RCircos.Get.Plot.Ideogram()
}
\value{
RCircos.Get.Plot.Parameters() returns a list that holding all parameters used in current RCircos session.

RCircos.Get.Plot.Positions() returns a data frame with three columns for x, and y coordinates of a circular line with radius of 1 as well as the degrees used for text rotation at each point of the positions. 

RCircos.Get.Plot.Ideogram() returns a data frame containing all information for chromosome ideogram plot. The first five column are orginal cytoband data and following are chromosome highlight colors, and plot colors, length, number of chromosome unit, and relative location on the circular layout of each band.
}
\author{
Henry Zhang
}

\keyword{ methods }

