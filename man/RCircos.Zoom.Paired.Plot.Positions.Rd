\name{RCircos.Zoom.Paired.Plot.Positions}
\alias{RCircos.Zoom.Paired.Plot.Positions}
\title{
Calculate Zoomed Index of Circos Plot Positions for Genomic Intervals
}
\description{
Scale up the index of Circos plot positions of a small range for zoom-in plot with data that require two positions to plot e.g., tiles and parallele lines. This function is mainly for internal use.
}
\usage{
RCircos.Zoom.Paired.Plot.Positions(zoom.data=NULL, zoom.pos=NULL)
}
\arguments{
  \item{zoom.data}{
Data frame, data to be plotted onto zoom area.
}
  \item{zoom.pos}{
Non-negative integer, index of RCircos base plot positions.
}
}
\value{
A two columns integer matrix representing the index of RCircos plot positions.
}
\author{
Henry Zhang
}
\keyword{ methods }
