\name{RCircos.Set.Plot.Area}
\alias{RCircos.Set.Plot.Area}
\title{
Setup Plot Area for RCircos Plot
}
\description{
Open an new window for plotting. RCircos core components must be initialized before call this function. All four margins are set to 0.25 and window size are set to plot.radius in RCircos parameters. This function is optional and could be replaced by call par(), plot.new(), and plot.window() if users want control the plot window by themselves.
}
\usage{
RCircos.Set.Plot.Area(margins = 0.25, ...)
}
\arguments{
  \item{margins}{
Non-negative numeric, value for margins around plot area.
}
  \item{\dots}{
Other arguments passed...
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram);
RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=5) 
RCircos.Set.Plot.Area() 
}
}
\keyword{ methods }

