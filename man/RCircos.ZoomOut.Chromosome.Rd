\name{RCircos.ZoomOut.Chromosome}
\alias{RCircos.ZoomOut.Chromosome}
\title{
Zoom out Chromosome Ideogram
}
\description{
Zoom out chromosome ideogram to leave a small room between the first and the last chromosome (at 12 O'clock) so that track names could be added. This function works on the RCircos chromosome ideogram component so the core components must be set first.
}
\usage{
RCircos.ZoomOut.Chromosome(zoom.out.ratio=NULL)
}
\arguments{
  \item{zoom.out.ratio}{
A float number between 0 and 1, the ratio of target chromosome ideogram to be zoomed out.
}
}
\value{
None. New ideogram object is set to RCircos environment.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=0) RCircos.ZoomOut.Chromosome(zoom.out.ratio=0.95)
}
}
\keyword{ methods }