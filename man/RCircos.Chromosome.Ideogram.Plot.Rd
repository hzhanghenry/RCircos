\name{RCircos.Chromosome.Ideogram.Plot}
\alias{RCircos.Chromosome.Ideogram.Plot}
\alias{RCircos.Draw.Chromosome.Ideogram}
\alias{RCircos.Highligh.Chromosome.Ideogram}
\alias{RCircos.Ideogram.Tick.Plot}
\alias{RCircos.Label.Chromosome.Names}
\title{
Plot Chromosome Ideogram
}
\description{
Draw chromosome ideogram, chromosome highlights and chromosome names. Graphic  device must be initialized first. The original function was split to three new  ones in order to plot chromosome highlight or names in different area. RCircos core components and graphic device must be initialized first.
}
\usage{
RCircos.Chromosome.Ideogram.Plot(tick.interval=0)
RCircos.Draw.Chromosome.Ideogram(ideo.pos=NULL, ideo.width=NULL)
RCircos.Highligh.Chromosome.Ideogram(highlight.pos=NULL, highlight.width=NULL)
RCircos.Ideogram.Tick.Plot(tick.interval=50, track.for.ticks=3)
RCircos.Label.Chromosome.Names(chr.name.pos=NULL)
}
\arguments{
  \item{tick.interval}{
Non-negative integer, the distance between two chromosome ideogram ticks in million base pairs. Set to 0 (zero) will ignore the ideogram ticks.
}
\item{ideo.pos}{
Non-negative numeric, the distance of chromosome ideogram position from center of plot area.
}
\item{ideo.width}{
Non-negative numeric, the width of chromosome ideogram.
}
\item{highlight.pos}{
Non-negative numeric, the distance of chromosome highlight line position from center of plot area.
}
\item{highlight.width}{
Non-negative integer, the width of chromosome highlight lines.
}
  \item{track.for.ticks}{
Non-negative integer, total number of tracks for ideogram tick height.
}
  \item{chr.name.pos}{
Non-negative numeric, the distance of chromosome names from center of plot area.
}
}
\details{
RCircos.Chromosome.Ideogram.Plot() is the default method to plot chromosome ideogram, chromosome highlights, and chromosome names.  

RCircos.Draw.Chromosome.Ideogram() plots chromosome ideogram only.  

RCircos.Highligh.Chromosome.Ideogram() highlights chromosomes with different color lines only.  

RCircos.Ideogram.Tick.Plot() draws chromosome ideogram ticks starting from chromosome highlight positions.  

RCircos.Label.Chromosome.Names(chr.name.pos=NULL) add chromosome names for chromosome ideogram at default or customized location.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,  
    tracks.inside=10, tracks.outside=0 )
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
}
}
\keyword{ methods }

