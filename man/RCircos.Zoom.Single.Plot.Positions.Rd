\name{RCircos.Zoom.Single.Plot.Positions}
\alias{RCircos.Zoom.Single.Plot.Positions}
\title{
Calculate Zoomed Index of Circos Plot Positions for Single Point Genomic Data 
}
\description{
Scale up the index of Circos plot positions of a small range for zoom-in plot with data that require only one position to plot e.g., scatter plot. This function is mainly for internal use.
}
\usage{
RCircos.Zoom.Single.Plot.Positions(zoom.data=NULL, zoom.pos=NULL)
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
An integer representing the index of RCircos plot positions.
}
\author{
Henry Zhang
}
\keyword{ methods }
