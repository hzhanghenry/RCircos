\name{RCircos.Set.Core.Components}
\alias{RCircos.Set.Core.Components}
\title{
Initialize Core Components for RCircos Plot
}
\description{
Setup core components for RCircos Plot including plot parameters, cytoband data, and base plot positions. This function must be called before chromosome ideogram and any other data plot.
}
\usage{
RCircos.Set.Core.Components(cyto.info, chr.exclude = NULL, tracks.inside, tracks.outside)
}
\arguments{
  \item{cyto.info}{
A data frame with chromosome ideogram data and columns are ordered by chromosome names, start position, end position, band names, and stain intensity for each chromosome band.
}
  \item{chr.exclude}{
Character vector of chromosome names to be excluded from plot, e.g., chr.exclude <- c("chrX", "chrY");
}
  \item{tracks.inside}{
Non-negative integer, total number of data tracks inside of chromosome ideogram to be plotted.
}
  \item{tracks.outside}{
Non-negative integer, total number of data tracks outside of chromosome ideogram to be plotted.
}
}
\author{
Hongen Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=0) 
} 
}
\keyword{ methods }
