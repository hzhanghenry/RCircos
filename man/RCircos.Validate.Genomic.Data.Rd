\name{RCircos.Validate.Genomic.Data}
\alias{RCircos.Validate.Genomic.Data}
\title{
Validate Genomic Data
}
\description{
Validate input dataset for correct chromosome names, chromosome start, and chromosome end positions. Chromosome names will be converted to character vectors if they are factor variables.
}
\usage{
RCircos.Validate.Genomic.Data(genomic.data=NULL, 
    plot.type=c("plot", "link"), genomic.columns=3)
}
\arguments{
  \item{genomic.data}{
Data frame with genomic position data. 
}
  \item{plot.type}{
Character vector, either "plot" or "link".
}
  \item{genomic.columns}{
Non-negative integer, total number of columns for genomic position (chromosome name, start and/or end position).
}
}
\value{
None. Function will exit with an error prompt.
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos)
data(RCircos.Heatmap.Data)
RCircos.Validate.Genomic.Data(RCircos.Heatmap.Data, 
    plot.type="plot", genomic.columns=3)

data(RCircos.Link.Data)
RCircos.Validate.Genomic.Data(RCircos.Link.Data, 
    plot.type="link", genomic.columns=3)
}
}
\keyword{ methods }
