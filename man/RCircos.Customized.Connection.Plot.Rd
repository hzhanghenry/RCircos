\name{RCircos.Customized.Connection.Plot}
\alias{RCircos.Customized.Connection.Plot}
\title{
Draw Connectors between Chromosome Ideogram and Gene Names
}
\description{
Draw connection lines (simple lines plot) between two sets of data points. One example of usage is to label genes at modified plot position and connect the gene label to its genomic position.os core components and graphics device must be initialized before drawing.
}
\usage{
RCircos.Customized.Connection.Plot(gene.data, label.data, 
			gene.pos=NULL, label.pos=NULL)
}
\arguments{
  \item{gene.data}{
A data frame with chromosome names and actual genomic positions for the genes to be labeled.
}
  \item{label.data}{
A data frame with chromosome name and genomic positions that will be used to label gene names.
}
  \item{gene.pos}{
Float numeric, scale factor relative to the center of plot area (0) for start points of connection lines.
}
  \item{label.pos}{
Float numeric, scale factor relative to the center of plot area (0) for end points of connection lines.
}

}
\author{
Henry Zhang
}
\keyword{ methods }

