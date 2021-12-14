\name{RCircos.Get.Polygon.Height}
\alias{RCircos.Get.Polygon.Height}
\title{
Calculate the Height of Polygons in A Data Track
}
\description{
Calculate polygon heights (bottom and top location) inside a plot track. Polygon can be draw from inside position to outside position,  from outside position to inside position, or from the middle to both direction of the track. This is mainly for internal use.
}
\usage{
RCircos.Get.Polygon.Height(data.heights, min.value=NULL,  
    max.value=NULL, inside.pos=NULL, outside.pos=NULL)
}
\arguments{
  \item{data.heights}{
Numeric vector, the data to be plotted on a data track.
}
  \item{min.value}{
Numeric, the minimum value of data range.
}
  \item{max.value}{
Numeric, the maximum value of data range.
}
  \item{inside.pos}{
Non-negative numeric, the close position of a plot track relative to the center of plot area.
}
  \item{outside.pos}{
Non-negative numeric, the far position of a plot track relative  to the center of plot area.
}
}
\value{
Non-negative umeric matrix of column 2 representing the bottom and top positions of polygons.
}
\author{
Henry Zhang
}
\keyword{ methods }

