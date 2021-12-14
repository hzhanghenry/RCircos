\name{RCircos.Validate.Track.Positions}
\alias{RCircos.Validate.Track.Positions}
\title{
Validate Track Positions
}
\description{
Validate if user defined plot locations (inside.pos and outside.pos) are inside of plot area and do not overlap with chromosome ideogram. RCircos core components must be initialized first. This function is mainly for internal use.
}
\usage{
RCircos.Validate.Track.Positions(inside.pos=0, outside.pos=0, 
	erase.area=FALSE)
}
\arguments{
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
A numeric vector with length of 2 for the real plot position of a customized track.
}
\author{
Henry Zhang
}
\keyword{ methods }
