\name{RCircos.Get.Plot.Boundary}
\alias{RCircos.Get.Plot.Boundary}
\title{
Calculate Boundary of A Data Track
}
\description{
Calculate plot track location based on track number, plot side, or user defined inside and outside position relative to plot area center. If all arguments are defined, user defined inside and outside position will be used. This function is mainly for internal use.
}
\usage{
RCircos.Get.Plot.Boundary(track.num=NULL, side=NULL, 
	inside.pos=NULL, outside.pos=NULL, erase.area=FALSE)
}
\arguments{
  \item{track.num}{
Non-negative integer, representing the ordinal number of the plot track where the lines will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram. It could be NULL if both inside.pos and outside.pos are defined.
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}
  \item{erase.area}{
Logic, if the track is going to be erased.
}
}
\value{
Numeric vector of length 2, x and y coordinates of a data track.
}
\author{
Henry Zhang
}
\keyword{ methods }
