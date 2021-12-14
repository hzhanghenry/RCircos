\name{RCircos.List.Plot.Parameters}
\alias{RCircos.List.Plot.Parameters}
\title{
List RCircos Plot Parameters
}
\description{
Print out onto screen all plot parameters used by current RCircos session. RCircos core components must be already initialized.
}
\usage{
RCircos.List.Plot.Parameters()
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
	tracks.inside=10, tracks.outside=0) RCircos.List.Plot.Parameters()
}
}
\keyword{ methods }
