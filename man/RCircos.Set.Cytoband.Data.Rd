\name{RCircos.Set.Cytoband.Data}
\alias{RCircos.Set.Cytoband.Data}
\title{
Set Chromosome Cytoband Data for Circos Plot
}
\description{
Set chromosome cytoband data (one of three RCircos core components) for Circos plot. The cytoband data will be put into RCircos Environment. This function is for internal use and users may not directly call this function and use get and reset methods instead when need to modify ideogram data.
}
\usage{
RCircos.Set.Cytoband.Data(cyto.band.info=NULL)
}
\arguments{
  \item{cyto.band.info}{
A data frame containing the chromosome ideogram data returned from RCircos.Validate.Cyto.Info() or a sorted ideogram table with correct chromosome names,start and end positions, band names, and staing status.
}
}
\author{
Henry Zhang
}
\keyword{ methods }

