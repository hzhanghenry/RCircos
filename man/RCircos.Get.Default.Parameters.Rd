\name{RCircos.Get.Default.Parameters}
\alias{RCircos.Get.Default.Circos.Units}
\alias{RCircos.Get.Default.Base.Per.Units}
\alias{RCircos.Get.Default.Chrom.Padding}
\alias{RCircos.Get.Padding.Constant}
\alias{RCircos.Get.Supported.HeatmapColors}
\alias{RCircos.Get.Supported.Plot.Types}
\alias{RCircos.Get.Default.Char.Width}
\alias{RCircos.Get.Default.Text.Size}
\title{
Retrieve Parameter Values Stored in RCircos Environment
}
\description{
These methods are used to view the defaul parameters stored in RCircos Environments.
}

\usage{
RCircos.Get.Default.Circos.Units()
RCircos.Get.Default.Base.Per.Units()
RCircos.Get.Default.Chrom.Padding()
RCircos.Get.Padding.Constant()
RCircos.Get.Supported.HeatmapColors()
RCircos.Get.Supported.Plot.Types()
RCircos.Get.Default.Char.Width()
RCircos.Get.Default.Text.Size()
}
\value{
RCircos.Get.Default.Circos.Units() returns the default number of units (points) for a circular line. This values is derived from  human genome length (3x10^9) and one unit represents length of 3000 base pairs.

RCircos.Get.Default.Base.Per.Units() returns the default number of base pairs that a unit (point) will represent.

RCircos.Get.Default.Chrom.Padding() returns the default lenght of gaps between two neighbot chromosomes in chromosome (Circos) unit.

RCircos.Get.Padding.Constant() returns the default ratio of chromosome padding length and total default Circos (chromosome) units.

RCircos.Get.Supported.HeatmapColors() returns the types of colors for heatmap plot.

RCircos.Get.Supported.Plot.Types() returns the plot types supported by current RCircos version.

RCircos.Get.Default.Char.Width() returns the default number of units a character will span on the plot.

RCircos.Get.Default.Text.Size() returns the default size of text in RCircos.
}
\author{
Henry Zhang
}
\keyword{ methods }
