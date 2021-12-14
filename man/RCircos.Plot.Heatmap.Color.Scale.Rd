\name{RCircos.Plot.Heatmap.Color.Scale}
\alias{RCircos.Plot.Heatmap.Color.Scale}
\title{
Plot A Color Scale for Heatmap
}
\description{
Plot a color scale showing the color ranges of heatmap. One or more tracks should be defined for plots in outside of chromosome ideogram.
}
\usage{
RCircos.Plot.Heatmap.Color.Scale(max.value=NULL, min.value=NULL,
        color.type="BlueWhiteRed", scale.location=1, 
        scale.width=0, scale.height=0)
}
\arguments{
  \item{max.value}{
Numeric, maximum value of heatmap scale
}
  \item{min.value}{
Numeric, minimum value of heatmap scale
}
  \item{color.type}{
Character vector for color specification, either "BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed","GreenBlackRed", "YellowToRed", or "BlackOnly". Default is "BlueWhiteRed".
}
  \item{scale.location}{
Non negative integer representing color scale location, 1~4 for 
bottom, left, top, and right.
}
  \item{scale.width}{
Non-negative numeric, width of color scale. if not defined, 1/2 of x axis will be used
}
  \item{scale.height}{
Non-negative numeric, height of color scale. if not defined, 1/10 of scaleWidth axis will be used
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
plot.new()
plot.window(c(-5, 5), c(-5, 5))
RCircos.Plot.Heatmap.Color.Scale(max.value=3, min.value=-3,
    scaleLocation=1, colorType="BlueWhiteRed", 
    scaleWidth=0, scaleHeight=0)
}
}
\keyword{ methods }
