\name{RCircos.Reset.Core.Components}
\alias{RCircos.Reset.Plot.Ideogram}
\alias{RCircos.Reset.Plot.Parameters}
\alias{RCircos.Reset.Plot.Positions}
\title{
Reset Methods for Modifying RCircos Core Components
}
\description{
Reset RCircos core components: plot ideogram, plot position, and plot parameters.
}
\usage{
RCircos.Reset.Plot.Ideogram(chrom.ideo)
RCircos.Reset.Plot.Positions(plot.positions)
RCircos.Reset.Plot.Parameters(new.params)
}
\arguments{
  \item{chrom.ideo}{
Data frame, object of RCircos cytoband data returned from RCircos.Get.Plot.Ideogram(). 
}
\item{plot.positions}{
Data frame, object of RCircos plot positions returned from RCircos.Get.Plot.Positions().
}
\item{new.params}{
List, object of RCircos plot positions returned from RCircos.Get.Plot.Parameters()
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

ideo <- RCircos.Get.Plot.Ideogram()
rows <- which(ideo$BandColor=="red"
ideo$BandColor[rows] <- "green"
RCircos.Reset.Plot.Ideogram(ideo)

params <- RCircos.Get.Plot.Parameters()
params$base.per.unit <- 6000
RCircos.Reset.Plot.Paramters(params)

positions <- RCircos.Get.Plot.Positions();
positions[1, 3] <- 0;
RCircos.Reset.Plot.Positions(positions);
}
}
\keyword{ methods }
