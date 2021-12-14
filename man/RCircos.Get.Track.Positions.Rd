\name{RCircos.Get.Track.Positions}
\alias{RCircos.Get.Track.Positions}
\title{
Calculate Plot Positions for A Data Track 
}
\description{
Calculate inner and outer plot positions for a data track. RCircos core components must be initialized first. This function is mainly for internal use.
}
\usage{
RCircos.Get.Track.Positions(side=NULL, track.num=NULL)
}
\arguments{
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{track.num}{
Integer, representing the ordinal number of the plot track
}
}
\value{
Numeric vector of length 2. The outer and inner positions of a track.
}
\author{
Hongen Zhang
}
\keyword{ method }

