\name{RCircos.Get.Plot.Colors}
\alias{RCircos.Get.Plot.Colors}
\title{
Assign Plot Colors for Non-Link Data Plot
}
\description{
Assign plot colors for each data point of genomic data except of link plot (link lines or ribbons). Users can define their own colors by appending a column with color names to plot data and give the column name as "PlotColor". Otherwise, default colors will be used.

This function is for internal use only.
}
\usage{
RCircos.Get.Plot.Colors(plot.data, color)
}
\arguments{
  \item{plot.data}{
A data frame with the first three columns as chromosome name, start and end position followed by data values for one or more samples.
}
  \item{color}{
Character vector of color names.
}
}
\value{
A vector of color names with length same as the rows of input data.
}
\author{
Hongen Zhang
}
\keyword{ methods }

