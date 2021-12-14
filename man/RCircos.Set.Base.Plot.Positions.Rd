\name{RCircos.Set.Base.Plot.Positions}
\alias{RCircos.Set.Base.Plot.Positions}
\title{
Setup Base Plot Positions
}
\description{
Calculate x and y coordinates for a circular line that served as base plot positions to derive positions for all other plot items. Total number of coordinates are based on total chromosome units and padding. Degrees that text should be rotated at each point are also calculated. This function is for internal use.
}
\usage{
RCircos.Set.Base.Plot.Positions(total.points=NULL)
}
\arguments{
  \item{total.points}{
Non-negative integer, total number of points to form a circular line.
}
}
\author{
Henry Zhang
}
\keyword{ methods }
