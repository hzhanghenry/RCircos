\name{RCircos.Adjust.Scatter.Values}
\alias{RCircos.Adjust.Scatter.Values}
\title{
Adjust Scatter Plot Data to Fit the Data Track
}
\description{
Adjust scatter plot values to fit the plot track. The data with lowest value will be plot on the centre of first subtrack and the one with highest value will be on the centre of top subtrack. After this adjustment, values could be added to inside position of plot track as final plot positions. This function is mainly for internal use.
}
\usage{
RCircos.Adjust.Scatter.Values(scatter.values=NULL, min.value=NULL,  
    max.value=NULL, track.height=NULL, subtrack=NULL) 
}
\arguments{
  \item{scatter.values}{
Numeric vector, the data to be plotted.
}
  \item{min.value}{
Numeric, minimum value of scatter data.
}
  \item{max.value}{
Numeric, maximum value of scatter data.
}
  \item{track.height}{
Non-negative numeric, height of plot track.
}
  \item{subtrack}{
Non-negative nummeric, the number of sub tracks.
}
}
\author{
Henry Zhang
}
\keyword{ methods }

