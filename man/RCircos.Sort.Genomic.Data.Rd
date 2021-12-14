\name{RCircos.Sort.Genomic.Data}
\alias{RCircos.Sort.Genomic.Data}
\title{
Sort Genomic or Ideogram Data
}
\description{
Sort genomics/ideogram data. The order of chromosome names should be numeric names (integers or Roman numbers) first then character names. If chromosome names are all characters alphabets order will be used. This function could be used before making RCircos plot.
}
\usage{
RCircos.Sort.Genomic.Data(genomic.data=NULL, is.ideo=FALSE)
}
\arguments{
  \item{genomic.data}{
A data frame with the first two or three columns for chromosome names, start and end positions. If it is ideogram data, next two columns must be band names, and Giemsa stain status.
}
  \item{is.ideo}{
Logic, weither the genomic data is ideogram or not.
}
}
\value{
Data frame same as input data but ordered by chromosome names then start positions
}
\author{
Henry
}
\examples{
\dontrun{
library(RCircos)

data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
cyto <- UCSC.Mouse.GRCm38.CytoBandIdeogram
ideogram <- RCircos.Sort.Genomic.Data(  
    genomic.data=UCSC.Mouse.GRCm38.CytoBandIdeogram,  
    is.ideo=TURE)

data(RCircos.Line.Data)
lineData <- RCircos.Sort.Genomic.Data( 	genomic.data=RCircos.Line.Data, 
    is.ideo=FALSE)
}
}
\keyword{ methods }
