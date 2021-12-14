\name{RCircos.ZoomIn.Plot}
\alias{RCircos.Plot.Zoomed.Area}
\alias{RCircos.Plot.Zoomed.Gene.Connectors}
\alias{RCircos.Plot.Zoomed.Heatmap}
\alias{RCircos.Plot.Zoomed.Histogram}
\alias{RCircos.Plot.Zoomed.Ideogram.Ticks}
\alias{RCircos.Plot.Zoomed.Continue.Lines}
\alias{RCircos.Plot.Zoomed.Parallel.Lines}
\alias{RCircos.Plot.Zoomed.Vertical.Lines}
\alias{RCircos.Plot.Zoomed.Scatters}
\alias{RCircos.Plot.Zoomed.Tiles}
\alias{RCircos.Plot.Zoomed.Polygons}
\alias{RCircos.Label.Zoom.Region}
\alias{RCircos.Mark.Zoom.Area}
\title{
Zoomed-in Plot of Small Genomic Region
}
\description{
Zoom-in a small area on outide of chromosome ideogram to show the details for: heatmap, histogram, gene connector, gene labels, continue lines, vertical lines, parallel lines, scatters(points), areas, and tiles.
}
\details{
RCircos.Plot.Zoomed.Area() paints partial area inside the zoomed region (an irregular polygon ship). Plot types include "mountain", "curtain" (reversed mountain), and "band".

RCircos.Plot.Zoomed.Gene.Connectors() plot zoomed-in connectors for small number of genes(default 11).

RCircos.Plot.Zoomed.Heatmap() plot zoomed-in heatmap for small number of genes (default 11).

RCircos.Plot.Zoomed.Histogram() plot zoomed-in histogram for small number of genes (default 11).

RCircos.Plot.Zoomed.Ideogram.Ticks() plot chromosome ticks for a small genomic intervals.

RCircos.Plot.Zoomed.Continue.Lines() plot zoomed-in continue lines between neighbors of small number genes (default 11).

RCircos.Plot.Zoomed.Parallel.Lines() plot zoomed-in parallel lines for small number of genomic intevals (default 11).

RCircos.Plot.Zoomed.Vertical.Lines() plot zoomed-in vertical lines from bottom to top of a track for small number of genomic intevals (default 11).

RCircos.Plot.Zoomed.Scatters() plot zoomed-in scatters for small number of genes(default 11).

RCircos.Plot.Zoomed.Tiles() plot zoomed-in tiles for small number of genes(default 11).

RCircos.Plot.Zoomed.Polygons() plot zoomed-in polygons in a samll genomic region.

RCircos.Label.Zoom.Region() plot gene names for small number of genes(default 11).

RCircos.Mark.Zoom.Area() highligh the area between zoomed-in plot area and chromosome ideogram to mark the original genomic interval.
}
\usage{
RCircos.Plot.Zoomed.Area(zoom.data=NULL, 
    plot.type="mountain", data.col=NULL, track.num=NULL, 
    zoom.pos=NULL, min.value=NULL, max.value=NULL, 
    area.color="gray", border.col="black", inside.pos=NULL, 
    outside.pos=NULL, outline=TRUE)
 
RCircos.Plot.Zoomed.Gene.Connectors(zoom.data=NULL,  
    track.num=NULL, zoom.pos=NULL, line.width=1,  
    inside.pos=NULL, outside.pos=NULL)  

RCircos.Plot.Zoomed.Heatmap(zoom.data=NULL, data.col=NULL,  
    track.num=NULL, zoom.pos=NULL, min.value=NULL,  
    max.value=NULL, inside.pos=NULL, outside.pos=NULL)  

RCircos.Plot.Zoomed.Histogram(zoom.data=NULL, data.col=NULL,  
    track.num=NULL, zoom.pos=NULL, min.value=NULL,  
    max.value=NULL, inside.pos=NULL, outside.pos=NULL, 
    outline=TRUE)  

RCircos.Plot.Zoomed.Ideogram.Ticks(zoom.info=NULL,  
    track.num=NULL, zoom.pos=NULL, tick.interval=5,  
    inside.pos=NULL, outside.pos=NULL)  

RCircos.Plot.Zoomed.Continue.Lines(zoom.data=NULL,  
    data.col=NULL, track.num=NULL, zoom.pos=NULL,  
    line.width=1, min.value=NULL, max.value=NULL,  
    inside.pos=NULL, outside.pos=NULL, outline=TRUE)  

RCircos.Plot.Zoomed.Parallel.Lines(zoom.data=NULL, 
    track.num=NULL, zoom.pos=NULL, genomic.cols=3, 
    line.width=NULL,inside.pos=NULL, outside.pos=NULL, 
    outline=FALSE)

RCircos.Plot.Zoomed.Vertical.Lines(zoom.data=NULL,  
    track.num=NULL, zoom.pos=NULL, line.width=1,  
    inside.pos=NULL, outside.pos=NULL, outline=FALSE)  

RCircos.Plot.Zoomed.Scatters(zoom.data=NULL, data.col=NULL,  
    track.num=NULL, zoom.pos=NULL, min.value=NULL,  
    max.value=NULL, point.type=16, by.fold=0,  
    with.size=TRUE, with.height=FALSE, point.scale=1,  
    inside.pos=NULL, outside.pos=NULL, outline=TRUE)  

RCircos.Plot.Zoomed.Tiles(zoom.data=NULL, track.num=NULL,  
    zoom.pos=NULL, genomic.cols=3, layers=5, border.col=NULL,  
    inside.pos=NULL, outside.pos=NULL, outline=TRUE)  

RCircos.Plot.Zoomed.Polygons(zoom.data=NULL, data.col=4,   
    track.num=NULL, zoom.pos=NULL, genomic.cols=3,  
    min.value=NULL, max.value=NULL, border.col=NULL,  
    inside.pos=NULL, outside.pos=NULL, outline=TRUE) 

RCircos.Label.Zoom.Region(zoom.data=NULL, name.col=NULL,  
    track.num=NULL, zoom.pos=NULL, text.size=0.75,  
    inside.pos=NULL, outside.pos=NULL)  

RCircos.Mark.Zoom.Area(zoom.range=NULL, track.num=1,  
    zoom.pos=NULL, fill.color="yellow",  
    inside.pos=NULL, outside.pos=NULL)  
}
\arguments{
  \item{zoom.data}{
A data frame containing genomic positions, gene names, and plot values for zoom-in genes/rows.
}  \item{zoom.range}{
Character vecter, chromosome name, start and end positions of a genomic interval to be zoomed in.
}
  \item{data.col}{
Non-negative integer/vector, which column(s) is plot data values.
}
  \item{track.num}{
Non-negative integer, which track will be plotted.
}
  \item{zoom.pos}{
Non-negative numeric vector, the index of RCircos plot position.
}
  \item{min.value}{
Numeric, minimum value for heatmap scale.
}
  \item{max.value}{
Numeric, maximum value for heatmap scale.
}
  \item{inside.pos}{
Non-negative numeric, inside position (relative to the centre of plot area) of the track.
}
  \item{outside.pos}{
Non-negative numeric, outside position (relative to the centre of plot area) of the track.
}
\item{name.col}{
Non-negative integer, which column is for labels.
}
\item{genomic.cols}{
Non-negative integer, total number of columns for genomic position information (chromosome, start and/or end position).
}
\item{fill.color}{
Character vector for name of a color to fill the polygon.
}
\item{border.col}{
Character vector for name of a color to draw the polygon border.
}
\item{area.color}{
Character vector for name of a color to fill the area.
}
\item{text.size}{
Non-negative numeric, character size for text plot.
}
\item{line.width}{
Non-negative numeric, width of lines.
}
\item{point.type}{
Non-negative integer, type of character for point plot (Default 19).
}
\item{plot.type}{
Character vector, either "mountain", "curtain", or "ban", types of area plot.
}
\item{by.fold}{
Zero or a positive float number, if it is greater than zero, any data point with a value greater than or equal to the by.fold will be plot as red color and data points with values less or equal to the negative by.fold will be plot as blue colr. Otherwise, black color is used.
}
\item{with.size}{
Logic. If TRUE, the point size will match to the data value in scatter/point plot.
}
\item{with.height}{
Logic. If TURE, point location will match to the data value in scatter/point plot.
}
\item{point.scale}{
Non-negative numeric, more scale to add for scatter/point size in scatter/point plot.
}
\item{outline}{
Logic. If TURE, outline will be drawn for zoomed area.
}
\item{zoom.info}{
Character vector contains chromosome name, start and end position.
}
\item{tick.interval}{
Non-negative integer, length of genomic interval in million base pairs.
}
\item{layers}{
Non-negative integer, number of subtrack in the data track.
}
}
\author{
Henry Zhang
}
\examples{
\dontrun{
library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram)

cyto <- UCSC.HG19.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto, NULL, 10, 5)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data(RCircos.Heatmap.Data)
min.value <- min(as.matrix(RCircos.Heatmap.Data[,5:10]))
max.value <- max(as.matrix(RCircos.Heatmap.Data[,5:10]))

zoom.data <- RCircos.Get.Zoom.Data(RCircos.Heatmap.Data,
        name.col=4, genomic.columns=3,
        target.gene="SP5", neighbor.genes=5)
zoom.range <- RCircos.Get.Zoom.Range(zoom.data, 3)
zoom.pos <- RCircos.Set.Zoom.Plot.Positions(zoom.range,
        total.genes=11, area.length=0.25, gene.width=NULL)

RCircos.Plot.Zoomed.Heatmap(zoom.data, data.col=5,
        track.num=2, zoom.pos=zoom.pos,
        min.value=min.value, max.value=max.value,
        inside.pos=NULL, outside.pos=NULL)

RCircos.Mark.Zoom.Area(zoom.pos, zoom.data, color="red",
        track.num=2, inside.pos=NULL, outside.pos=NULL)

RCircos.Label.Zoom.Region(zoom.data, name.col=4,
        track.num=3, zoom.pos=zoom.pos, text.size=0.75,
        inside.pos=NULL, outside.pos=NULL)
}
}
\keyword{ methods }
