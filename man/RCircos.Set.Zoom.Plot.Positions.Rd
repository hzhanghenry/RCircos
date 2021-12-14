\name{RCircos.Set.Zoom.Plot.Positions}
\alias{RCircos.Set.Zoom.Plot.Positions}
\title{
Set Position for Zoom-in Plot
}
\description{
Set plot position for zooming-in plot by selecting index of default plot positions. The zoomed plot will be aligned with the original plot position such as gene position or a cytoband. The zoom area was defined as fraction of circle circumference, default is 1/8.
}
\usage{
RCircos.Set.Zoom.Plot.Positions(zoom.info=NULL, total.genes=11,
    area.length=0.25, fixed.width=FALSE, gene.width=NULL)
}
\arguments{
  \item{zoom.info}{
A vector with chromosome name, start and end position to be zoomed.
}
  \item{total.genes}{
Non-negative integer from 1 ~ maximum number of genes to be plotted. The maximum number of genes is calculated based on
track height and will be adjusted inside the function.
}
  \item{area.length}{
Non-negative numeric, for better layout, it must be smaller than or equal to 1/4 (quarter of the circumstance).
}
  \item{fixed.width}{
Logic, if the plot witdth for each data point is fixed or not.
}
  \item{gene.width}{
Non-negative integer, number of units for gene width.
}
}
\value{
Numeric vector for index of RCircos.Base.Position
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG19.Human.CytoBandIdeogram);

RCircos.Set.Core.Components(
	cyto.info=UCSC.HG19.Human.CytoBandIdeogram, 
	chr.exclude=c("chrX", "chrY"), 
	tracks.inside=10, tracks.outside=0) 

zoom.info <- c("chr1", 10000, 50000)
zoom.pos <- RCircos.Set.Zoom.Plot.Positions(zoom.info)
}
}
\keyword{ methods }
