#
#   Zoom-in plot for heatmap, point, bar/histogram, vLine, cLine .......
#   Note: zoom-in plot is always done at outside of chromosome ideogram.
#
#   There will be two categories of zoom plot:
#
#   1.  Fixed data width: Examples are heatmap, histogram(bar), connector, 
#       scatters, and gene names. They do not needs to deal with genomic 
#       interval length 
#
#   2.  Segment related:  Examples have vertical lines, continue lines,
#       parallel lines, tiles, 
#
#   Functions
#
#   1.  RCircos.Get.Zoom.Data()
#   2.  RCircos.Get.Zoom.Range()
#   3.  RCircos.Set.Zoom.Plot.Positions()
#   4.  RCircos.Mark.Zoom.Area()
#   5.  RCircos.Label.Zoom.Region()
#   6.  RCircos.Plot.Zoomed.Heatmap()
#   7.  RCircos.Plot.Zoomed.Histogram()
#   8.  RCircos.Plot.Zoomed.Gene.Connectors()
#   9.  RCircos.Plot.Zoomed.Vertical.Lines()
#   10. RCircos.Plot.Zoomed.Continue.Lines()
#   11. RCircos.Plot.Zoomed.Parallel.Lines()
#   12. RCircos.Plot.Zoomed.Scatters()
#   13. RCircos.Plot.Zoomed.Tiles()
#   14. RCircos.Plot.Zoomed.Ideogram.Ticks()
#   15. RCircos.Plot.Zoomed.Polygons()
#   16. RCircos.Zoom.Area.Outline()
#   17. RCircos.Clear.Zoom.Area()
#   18. RCircos.Zoom.Single.Plot.Positions()
#   19. RCircos.Zoom.Paired.Plot.Positions()
#   20. RCircos.Plot.Zoomed.Area()
#
#   Procedures:
#
#   zoom.data <- RCircos.Get.Zoom.Data(plot.data, 
#                   name.col=4, target.gene="TP53");
#   zoom.info <- RCircos.Get.Zoom.Range(zoom.data, 
#                   genomic.cols=3)
#   zoom.pos  <- RCircos.Set.Zoom.Plot.Positions(zoom.info)
#   RCircos.Plot.Zoomed.Heatmap(zoom.data, data.col=5, 
#                   zoom.pos, track.num=2, 
#                   min.value=-3, max.value=3)
#   RCircos.Plot.Zoomed.Heatmap(zoom.data, data.col=6, 
#                   zoom.pos, track.num=3, 
#                   min.value=-3, max.value=3)
#   RCircos.Gene.Labels.On.Zoom(zoom.data, name.col=4, 
#                   zoom.pos, track.num=4)
#
#
#   Last modified on January 06, 2016
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>




#   =========================================================================
#
#   1.  RCircos.Get.Zoom.Data()
#
#   Get zoom data (subset of plot data) based on a gene name or row header
#
#   Arguments:
#
#       plot.data:          A data frame containing genomic positions, gene 
#                           names,and plot values. The data should be already  
#                           sorted by chromosome names then start positions.
#       name.col:           Non-negative integer, column for gene names
#       genomic.columns:    Non-negative integer, total number of columns for
#                           genomic position, valid values are 2 or 3.
#       target.gene:        Character vector, name of gene to be focused on
#       neighbor.genes:     Non-negative integer, how many genes will be plot
#                           on each side of target gene, default 5
#
#       Return value:       A data frame (subset of input dataset)
#       Example:            zoom.data <- RCircos.Get.Zoom.Data(plot.data,  
#                                           name.col=4, target.gene="TP53");
#

RCircos.Get.Zoom.Data <- function(plot.data=NULL, name.col=NULL, 
              genomic.columns=3, target.gene=NULL, neighbor.genes=5)
{
    if(is.null(plot.data) || is.null(name.col) || is.null(target.gene) ||
         is.null(genomic.columns) || is.null(neighbor.genes))
        stop("Missing argument(s) in RCircos.Get.Gene.Info().\n")
    if(genomic.columns < 2 || genomic.columns > 3)
        stop("Total columns for position information should be 2 or 3.\n");
    if(name.col <= genomic.columns) 
        stop("Incorrect name column defined.\n");

    gene.row <- which(as.character(plot.data[, name.col]) == target.gene);
    if(length(gene.row) > 1) gene.row <- gene.row[1];

    from <- gene.row - neighbor.genes;
    to   <- gene.row + neighbor.genes;
    target.rows <- from:to;
    zoom.data <- plot.data[target.rows, ]
    
    # in case target gene is on one end of chromosome
    target.chr <- as.character(plot.data[gene.row,  1])
    zoom.rows <- which(as.character(zoom.data[,1]) == target.chr)
    zoom.data <- zoom.data[zoom.rows,]

    return (zoom.data);
}




#   ==========================================================================
#
#   2.  RCircos.Get.Zoom.Range()
#
#   Extract chromosome name, start position, and end position from plot data
#   with a gene name or row header. If using this function to get zoom-in
#   range, there will be no need to call RCircos.Validate.Genomic.Info().
#
#   Argument:
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values
#       genomic.columns:  Non-negative integer, total number of columns for
#                       genomic position, valid values are 2 or 3.
#
#   Return value:       An vector containing chromosome name, start and end
#                       position.
#
#   Example:    zoom.data <- RCircos.Get.Zoom.Data(plot.data, name.col=4, 
#                                       target.gene="TP53");
#               zoom.info <- RCircos.Get.Zoom.Range(zoom.data, 3)
#
#

RCircos.Get.Zoom.Range <- function(zoom.data=NULL,  genomic.columns=3)
{
    if(is.null(zoom.data) || is.null(genomic.columns))
        stop("Missing argument(s) in RCircos.Get.Gene.Info().\n")
    if(genomic.columns < 2 || genomic.columns > 3)
        stop("Total columns for genomic position should be 2 or 3.\n");

    zoom.chr <- unique(as.character(zoom.data[, 1]))
    if(length(zoom.chr) > 1) 
        stop("Only data on same chromosome could be zoomed-in.\n");

    zoom.start <- min(as.numeric(zoom.data[, 2]));
    zoom.end <- max(as.numeric(zoom.data[, genomic.columns]));
    zoom.range <- c(zoom.chr, zoom.start, zoom.end);    

    RCircos.Validate.Genomic.Info(zoom.range);

    return (zoom.range);
}




#    =========================================================================
#
#   3.  RCircos.Set.Zoom.Plot.Positions()
#
#   Set plot position for zooming-in pot area by selecting index of default
#   plot positions. The zoomed plot will be aligned with the original plot 
#   position such as gene position or a cytoband. The zoom area was defined 
#   as fraction of circle circumference, default is 1/4.
#
#    Arguments: 
#
#       zoom.info:      A vector contains chromosome name, start and optional
#                       end position to be zoomed in
#       total.genes:    integer of 1 ~ maximum number of genes to be plotted.
#                       The maximum number of genes is calculated based on
#                       track height and will be adjusted inside the function
#       area.length:    Non-negative numeric, for better layout, it must be 
#                       smaller than or equal to 1/4
#       gene.width:     Non-negative integer, number of units for gene width
#
#    Return value:      Numeric vector for index of RCircos.Base.Position
#    Example:           RCircos.Set.Zoom.Plot.Positions(area.length=4)
#

RCircos.Set.Zoom.Plot.Positions <- function(zoom.info=NULL, total.genes=11, 
                    area.length=0.25, fixed.width=FALSE, gene.width=NULL)
{
    if(is.null(zoom.info)) 
        stop("Missing zoom info in RCircos.Get.Zoom.Plot.Positions().\n");
    if(total.genes < 1 ) stop("Total genes to zoom must be 1 or more.\n")
 
    RCircos.Pos  <- RCircos.Get.Plot.Positions();
    RCircos.Par  <- RCircos.Get.Plot.Parameters();

    #   Get number of points for polygon width (same as rectangle with
    #   length of track height)
    #   ===================================================================
    #
    if(fixed.width == FALSE) 
    {
        area.points <- nrow(RCircos.Pos)*area.length;
        
    } else {
        if(is.null(gene.width)) {
            gene.width <- ceiling(nrow(RCircos.Pos) * 
                (RCircos.Par$track.height/(2*pi*RCircos.Par$track.out.start)));
        } else {
            if(gene.width<0) 
                stop("Area for a zoomed gene cannot be negative.\n");
        }
        max.genes <- floor(nrow(RCircos.Pos) * area.length / gene.width);
   
        if(total.genes > max.genes)
        {
            stop(paste0("Set area.length=", area.length, 
                " can only plot ", max.genes, "genes.\n"));
        }
        area.points <- ceiling(gene.width*total.genes);
    }  

    #    Default index of RCircos.Pos
    #    ===========================================================
    #
    point.index <- 1:area.points;
    area.center <- floor(area.points/2);

    #    Rotate the index clockwise so that its mid point is the 
    #    center of segment to be zoomed
    #    =========================================================
    #
    target.index <- RCircos.Data.Point(as.character(zoom.info[1]),
        mean(as.numeric(zoom.info[2]), as.numeric(zoom.info[2])));
    distance <- abs(target.index - area.center);
    point.index <- point.index + distance;    
    
    #   In case of rotate the index counter-clockwise 
    #   ============================================================
    #
    zeros <- which(point.index <= 0)
    if(length(zeros)>0) 
        point.index[zeros] <- point.index[zeros] + nrow(RCircos.Pos);

    return (point.index);
}



#   =========================================================================
#
#   4.  RCircos.Mark.Zoom.Area()
#
#   Plot colored area between chromosome highlight and the start location of
#   zoom plot area to mark the original location of zoom-in data on chromosome.
#
#   Arguments:
#
#       zoom.pos:       Non-negative numeric vector, the index of RCircos p
#                       lot position
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows       
#       track.num:      Non-negative integer, which track will be plotted
#       color:          Character vector for a color name or a R color
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#   Return value:   None
#   Example:        RCircos.Mark.Zoom.Area(zoom.pos, 2)
#

RCircos.Mark.Zoom.Area <- function(zoom.range=NULL, track.num=1, zoom.pos=NULL, 
        fill.color="yellow", inside.pos=NULL, outside.pos=NULL)
{
    if(is.null(zoom.pos) || is.null(zoom.range)) 
        stop("Zoom in position and zoom in data must be defined.\n")

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                        inside.pos, outside.pos);
    out.pos <- boundary[1];
    in.pos  <- boundary[2];

    zoom.start <- zoom.pos[1];
    zoom.end   <- zoom.pos[length(zoom.pos)];
    
    chr.start <- RCircos.Data.Point(as.character(zoom.range[1]),
                        as.numeric(zoom.range[2]));
    chr.end   <- RCircos.Data.Point(as.character(zoom.range[1]),
                        as.numeric(zoom.range[3]));

    rgb.color <- col2rgb(fill.color, alpha=FALSE);
    mark.color <- rgb(rgb.color[1,1]/255, rgb.color[2,1]/255, 
                        rgb.color[3,1]/255, 0.1);

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    highlight.pos <- RCircos.Par$highlight.pos;

    polygonX <- c(RCircos.Pos[zoom.start:zoom.end, 1]*in.pos, 
                RCircos.Pos[chr.end:chr.start, 1]*highlight.pos);
    polygonY <- c(RCircos.Pos[zoom.start:zoom.end, 2]*in.pos, 
                RCircos.Pos[chr.end:chr.start, 2]*highlight.pos);
    polygon(polygonX, polygonY, col=mark.color);
}




#   =========================================================================
#
#   5.  RCircos.Label.Zoom.Region()
#
#   Plot gene names for zoom-in track(s).
#
#   Arguments:
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows
#       name.col:       Non-negative integer, which column is labels
#       track.num:      Non-negative integer, which track will be plotted
#       zoom.pos:       Non-negative numeric vector, the index of RCircos p
#                       lot position
#       text.size:      Size of text for labels
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#
#   Return value:       None
#   Example:            RCircos.Gene.Labels.On.Zoom(zoom.data, 5, zoom.pos, 6)
#

RCircos.Label.Zoom.Region <- function(zoom.data=NULL, name.col=NULL, 
            track.num=NULL, zoom.pos=NULL, text.size=0.75, 
            inside.pos=NULL, outside.pos=NULL)
{    
    if(is.null(zoom.data) || is.null(name.col) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Heatmap().\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                        inside.pos, outside.pos);
    out.pos <- boundary[1];
    in.pos  <- boundary[2];

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters(); 

    total.gene <- nrow(zoom.data);
    gene.width <- floor(length(zoom.pos) / total.gene);
    label.pos <- (1:nrow(zoom.data))*gene.width - floor(gene.width/2);
    label.pos <- label.pos + zoom.pos[1];

    #   Label all gene names
    #   ==================================================================

    gene.names <- as.character(zoom.data[, name.col]);
    text.start   <- in.pos;
    text.loc <- rep(4, length(label.pos));
    text.loc[which(label.pos > (nrow(RCircos.Pos)/2))] <- 2;

    for(a.gene in seq_len(total.gene))
    {
        index <- label.pos[a.gene]
        text(RCircos.Pos[index, 1]*text.start, 
                RCircos.Pos[index, 2]*text.start, offset=0,
                cex=text.size, srt=RCircos.Pos[index, 3], 
                gene.names[a.gene], pos=text.loc[a.gene]);
    }
}




#    =========================================================================
#
#   6.  RCircos.Plot.Zoomed.Heatmap()
#
#   Plot zoom-in heatmap for a small area with fixed width for each gene. 
#
#   Arguments:
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows
#       data.col:       Non-negative integer, which column of plot data
#       track.num:      Non-negative integer, which track will be plotted
#       zoom.pos:       Non-negative numeric vector, the index of RCircos  
#                       plot position
#       min.value:      Numeric, the minimum value of heatmap 
#       max.value:      Numeric, the maximum value for heatmap
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#
#   Return value:   None
#   Example:        RCircos.Plot.Zoomed.Heatmap(zoom.data, 5, zoom.pos, 4)
#

RCircos.Plot.Zoomed.Heatmap <- function (zoom.data=NULL, data.col=NULL, 
            track.num=NULL, zoom.pos=NULL, min.value=NULL, max.value=NULL,
            inside.pos=NULL, outside.pos=NULL)
{
    if(is.null(zoom.data) || is.null(data.col) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Heatmap().\n");

    if(is.null(min.value) || is.null(max.value))
        stop("Min and max heatmap value must be defiend.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                            inside.pos, outside.pos);                           
    out.pos <- boundary[1];
    in.pos  <- boundary[2];
    
    #   heatmap plot for a sample in zoom-in area
    #   ============================================================
    #
    heatmap.values <- as.numeric(zoom.data[, data.col]);
    heatmap.colors <- RCircos.Get.Heatmap.Data.Colors(heatmap.values, 
                        min.value, max.value);
    gene.width <- floor(length(zoom.pos) / nrow(zoom.data));

    RCircos.Pos <- RCircos.Get.Plot.Positions();                
    for(a.gene in seq_len(length(heatmap.values)))
    {
            first <- zoom.pos[1] + (a.gene-1)*gene.width;
            last  <- first + gene.width;

            polygon(c(RCircos.Pos[first:last, 1]*in.pos, 
                RCircos.Pos[last:first, 1]*out.pos),
                c(RCircos.Pos[first:last , 2]*in.pos, 
                RCircos.Pos[last:first, 2]*out.pos),
                col=heatmap.colors[a.gene], border=NA);
    }
}




#   =========================================================================
#
#   7.  RCircos.Plot.Zoomed.Histogram()
#
#   Plot zoom-in histogram on zoom-in area with fixed width for each gene.
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows
#       data.col:       Non-negative integer, which column is labels
#       zoom.pos:       Non-negative numeric vector, the index of RCircos
#                       plot position
#       track.num:      Non-negative integer, which track will be plotted
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#   Return value:       None
#   Example:    RCircos.Plot.Zoomed.Histogram(zoom.data, 5, zoom.pos, 1);
#               RCircos.Plot.Zoomed.Histogram(zoom.data, 5, zoom.pos, 
#                                   outside.pos=5.4, inside.pos=5);
#

RCircos.Plot.Zoomed.Histogram <- function(zoom.data=NULL, data.col=NULL, 
            track.num=NULL, zoom.pos=NULL, min.value=NULL, max.value=NULL, 
            inside.pos=NULL, outside.pos=NULL, outline=TRUE) 
{
    if(is.null(zoom.pos) || is.null(zoom.data)) 
        stop("Zoom in position and zoom in data must be defined.\n")

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", inside.pos, 
                                               outside.pos);
    out.pos <- boundary[1];
    in.pos  <- boundary[2];
    gene.width <- floor(length(zoom.pos) / nrow(zoom.data));

    RCircos.Par <- RCircos.Get.Plot.Parameters();
    hist.color  <- RCircos.Get.Plot.Colors(zoom.data, RCircos.Par$hist.color);
    
    if(outline==TRUE)
    {
        layers <- RCircos.Par$sub.tracks;
        RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, layers);
    }
    
    histValues <- as.numeric(zoom.data[, data.col]);
    if(is.null(max.value) || is.null(min.value)) {
        max.value <- max(histValues);
        min.value <- min(histValues);
    } else {
        if(min.value > max.value) stop("min.value > max.value.")
    }
    hist.height <- RCircos.Get.Data.Point.Height(histValues, min.value, 
                    max.value, plot.type="points", out.pos-in.pos);
    hist.height <- hist.height + in.pos;

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    for(a.gene in seq_len(length(hist.height)))
    {
            first <- zoom.pos[1] + (a.gene-1)*gene.width;
            last  <- first + gene.width;

            polygon(c(RCircos.Pos[first:last, 1]*in.pos, 
                RCircos.Pos[last:first, 1]*hist.height[a.gene]),
                c(RCircos.Pos[first:last, 2]*in.pos, 
                RCircos.Pos[last:first, 2]*hist.height[a.gene]),
                col=hist.color[a.gene]);
    }
}




#   ==========================================================================
#
#   8. RCircos.Plot.Zoomed.Gene.Connectors()
#
#   Plot gene connectors to mark the original location of gene names plotted
#   on the zoom-in area so that more gene labels could be done for a small
#   genomic interval. Gene names should be plotted with fixed area size. 
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and must be sorted by start positions
#       zoom.pos:       Non-negative numeric vector, the index of RCircos p
#                       lot position
#       track.num:      Non-negative integer, the number of a track where
#                       connecter will be started. Default is NULL and if 
#                       inside.pos or outside.pos is not defined, total of 3
#                       tracks will be used to draw connectors.
#       genomic.cols:   Non-negative integer, total number of columns for
#                       genomic positions. Valid values are 2 or 3.
#       line.size:      Non-negative numeric, lwd parameter for line plot.
#       inside.pos:     Non-negative numeric, the close location relative  
#                       to center of plot area
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#
#   Return value:   None
#   Example:        RCircos.Plot.Zoomed.Gene.Connectors(zoom.data, zoom.pos, 
#                           track.num=1, genomic.cols=3, line.size=1);
#                   RCircos.Plot.Zoomed.Gene.Connectors(zoom.data, zoom.pos, 
#                               genomic.cols=3, line.size=1, outside.pos=2.0, 
#                               inside.pos=1.5);
#

RCircos.Plot.Zoomed.Gene.Connectors <- function(zoom.data=NULL, track.num=NULL,  
                zoom.pos=NULL, line.width=1, inside.pos=NULL, outside.pos=NULL) 
{
    if(is.null(zoom.data) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(is.null(line.width) || line.width < 1) line.size <- 1;
   
    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", inside.pos, 
                                              outside.pos);
    out.pos <- boundary[1]; 
    in.pos  <- boundary[2];

    track.height <- out.pos - in.pos;
    if(is.null(track.num) == FALSE) {
        out.pos <- in.pos + track.height*3;
        low.end   <- in.pos + track.height;
        top.start <- low.end + track.height;
    } else {
        low.end   <- in.pos + track.height * 0.25;
        top.start <- in.pos + track.height * 0.75;
    }

    gene.width <- floor(length(zoom.pos) / nrow(zoom.data));
    label.pos  <- (1:nrow(zoom.data))*gene.width - floor(gene.width/2);
    label.pos <- label.pos + zoom.pos[1];
    gene.pos  <- RCircos.Get.Single.Point.Positions(zoom.data, 
                            genomic.columns=2);
    gene.Pos <- as.numeric(gene.pos$Location)
   
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    line.colors <- RCircos.Get.Plot.Colors(zoom.data, RCircos.Par$line.color);
    
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    for(a.gene in seq_len(length(label.pos)))
    {
        line.x <- c(RCircos.Pos[gene.Pos[a.gene], 1]*in.pos,
                    RCircos.Pos[gene.Pos[a.gene], 1]*low.end,
                    RCircos.Pos[label.pos[a.gene], 1]*top.start,
                    RCircos.Pos[label.pos[a.gene], 1]*out.pos);

        line.y <- c(RCircos.Pos[gene.Pos[a.gene], 2]*in.pos,
                    RCircos.Pos[gene.Pos[a.gene], 2]*low.end,
                    RCircos.Pos[label.pos[a.gene], 2]*top.start,
                    RCircos.Pos[label.pos[a.gene], 2]*out.pos);

        lines(line.x, line.y, lwd=line.width);
    }
}




#   ==========================================================================
#
#   9.  RCircos.Plot.Zoomed.Vertical.Lines()
#
#   Plot zoom-in vertical lines on zoom-in area.
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows
#       track.num:      Non-negative integer, which track will be plotted
#       zoom.pos:       Non-negative numeric vector, the index of RCircos
#                       plot position
#       line.width:     Non-negative integet, line width parameter
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#
#   Return value:   None
#   Example:        RCircos.Plot.Zoomed.Vertical.Lines(zoom.data, 
#                               zoom.pos, line.width=4, track.num=5)
#                    RCircos.Plot.Zoomed.Vertical.Lines(zoom.data,  
#                               zoom.pos, line.width=4, outside.pos=5.5, 
#                               inside.pos=4.5)
#

RCircos.Plot.Zoomed.Vertical.Lines <- function(zoom.data=NULL,  
                track.num=NULL, zoom.pos=NULL, line.width=1, 
                inside.pos=NULL, outside.pos=NULL, outline=FALSE) 
{
    if(is.null(zoom.data) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(line.width < 1) stop("Line width cannot be 0 or less.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", inside.pos, 
                                               outside.pos);
    out.pos <- boundary[1]; 
    in.pos  <- boundary[2];

    line.pos <- RCircos.Zoom.Single.Plot.Positions(zoom.data, zoom.pos);
    
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    line.color  <- RCircos.Get.Plot.Colors(zoom.data, RCircos.Par$line.color);
    
    if(outline==TRUE)
    {
        layers <- RCircos.Par$sub.tracks;
        RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, layers);
    }

    for(a.line in seq_len(length(line.pos)))
    {
        lines(c(RCircos.Pos[line.pos[a.line], 1]*in.pos, 
                RCircos.Pos[line.pos[a.line], 1]*out.pos),
              c(RCircos.Pos[line.pos[a.line], 2]*in.pos, 
                RCircos.Pos[line.pos[a.line], 2]*out.pos),
                col=line.color[a.line], lwd=line.width);
    }
}




#   ==========================================================================
#
#   10.  RCircos.Plot.Zoomed.Continue.Lines()
#
#   Plot continuous lines from one point to next point on zoom in area
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows
#       data.col:       Non-negative integer, which column is labels
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       line.width:     Non-negative integet, line width parameter
#       min.value:      Numeric, minimum value ifor point height
#       max.value:      Numeric, maximum value ifor point height
#       track.num:      Non-negative integer, which track will be plotted
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#   Return value:   None
#   Example:        RCircos.Plot.Zoomed.Cotinuous.Lines(zoom.data, 5, zoom.pos,
#                                           4, -5, 5, 1);
#                   RCircos.Plot.Zoomed.Cotinuous.Lines(zoom.data, data.col=5, 
#                           zoom.pos, line.width=4, min.value=-5, max.value=5, 
#                           outside.pos=5, inside.pos=4.5);
#

RCircos.Plot.Zoomed.Continue.Lines <- function(zoom.data=NULL,   
        data.col=NULL, track.num=NULL, zoom.pos=NULL, line.width=1, 
        min.value=NULL, max.value=NULL, inside.pos=NULL, 
        outside.pos=NULL, outline=TRUE) 
{
    if(is.null(zoom.data) || is.null(data.col) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(line.width < 1) stop("Line width cannot be 0 or less.\n");
    
    if(is.null(min.value) || is.null(max.value))
        stop("min.value and max.value must be defined.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                          inside.pos, outside.pos);
    out.pos <- boundary[1]; in.pos  <- boundary[2];
    track.height <- out.pos - in.pos;

    #   Points (line start and end) locations
    #   =============================================================
    #
    gene.width <- floor(length(zoom.pos) / nrow(zoom.data));
    point.pos <- (1:nrow(zoom.data))*gene.width - floor(gene.width/2);
    point.pos <- RCircos.Zoom.Single.Plot.Positions(zoom.data, zoom.pos);

    #   Point height
    #   =============================================================
    #
    point.value <- as.numeric(zoom.data[, data.col]);
    point.height <- RCircos.Get.Data.Point.Height(point.value, 
                min.value, max.value, "points", track.height);
    point.height <- point.height + in.pos;

    #   line colors
    #   ===================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    line.colors <- RCircos.Get.Plot.Colors(zoom.data, RCircos.Par$line.color); 

    #   Outline at zoom area
    #   =============================================================
    #
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    if(outline == TRUE)
    {
        layers <- RCircos.Par$sub.tracks;
        RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, layers);
    }

    for(a.point in seq_len((length(point.pos)-1)))
    {
        point.one <- point.pos[a.point];
        point.two <- point.pos[a.point+1];
        
        lines(c(RCircos.Pos[point.one, 1]*point.height[a.point], 
                RCircos.Pos[point.two, 1]*point.height[a.point+1]),
              c(RCircos.Pos[point.one, 2]*point.height[a.point], 
                RCircos.Pos[point.two, 2]*point.height[a.point+1]),
               col=line.colors[a.point]);
    }
}




#   ==========================================================================
#
#   11. RCircos.Plot.Zoomed.Parallel.Line()
#
#   Plot parallel links between two genomic positions on zoom in area.
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and layer values for zoom-in genes/rows
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       genomic.cols:   Non-negative integer, total number of columns for
#                       genomic positions. Valid values are 2 or 3.
#       track.num:      Non-negative integer, which track will be plotted
#       line.width:     Non-negative integet, line width parameter
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#   Return value:   None
#   Example:    RCircos.Plot.Zoomed.Scatters(zoom.data, 5, zoom.pos, 
#                       min.value=-10, max.value=10, track.num=1)
#               RCircos.Plot.Zoomed.Scatters(zoom.data, 5, zoom.pos,
#                       min.value=-10, max.value=10, 
#                       outside.pos=2.5, inside.pos=2);
#
RCircos.Plot.Zoomed.Parallel.Lines <- function(zoom.data=NULL, track.num=NULL, 
            zoom.pos=NULL, genomic.cols=3, line.width=NULL, 
            inside.pos=NULL, outside.pos=NULL, outline=FALSE) 
{
    if(is.null(zoom.data) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(line.width < 1) stop("Line width cannot be 0 or less.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                                inside.pos, outside.pos);
    out.pos <- boundary[1]; 
    in.pos  <- boundary[2];
    track.height <- out.pos - in.pos;

    #   Data points to link each other
    #   ================================================
    #
    line.data <- RCircos.Get.Paired.Points.Positions(zoom.data, 
                       genomic.cols, "tile");
    layers <- RCircos.Get.Plot.Layers(line.data, genomic.cols);
     if(max(layers)==1) {
        line.height <- in.pos + track.height/2;
    } else {
        line.height <- layers*(track.height/max(layers)) + in.pos;
    }
     
    line.location <- RCircos.Zoom.Paired.Plot.Positions(zoom.data, zoom.pos);
    line.start <- line.location[,1];
    line.end   <- line.location[,2];
    
    #   Get link line colors for each pair of locations
    #   ================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    line.colors <- RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color);

    #   Draw link lines for each pair of genomic positions
    #   ==========================================
    #
    if(outline == TRUE)
    { RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, max(layers));}
    for(a.line in seq_len(nrow(line.data)))
    {
        a.start <- line.start[a.line];
        a.end   <- line.end[a.line];
        
        lines(RCircos.Pos[a.start:a.end, 1]*line.height[a.line], 
              RCircos.Pos[a.start:a.end, 2]*line.height[a.line], 
               col=line.colors[a.line], lwd=line.width);
    }
}




#   ==========================================================================
#
#   12. RCircos.Plot.Zoomed.Scatters()
#
#   Plot scatters on zoom in area.
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows
#       data.col:       Non-negative integer, which column is labels
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       line.width:     Non-negative integet, line width parameter
#       min.value:      Numeric, minimum value ifor point height
#       max.value:      Numeric, maximum value ifor point height
#       track.num:      Non-negative integer, which track will be plotted
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#   Return value:   None
#   Example:    RCircos.Plot.Zoomed.Scatters(zoom.data, 5, zoom.pos, 
#                       min.value=-10, max.value=10, track.num=1)
#               RCircos.Plot.Zoomed.Scatters(zoom.data, 5, zoom.pos,
#                       min.value=-10, max.value=10, 
#                       outside.pos=2.5, inside.pos=2);
#

RCircos.Plot.Zoomed.Scatters <- function(zoom.data=NULL, data.col=NULL,  
            track.num=NULL, zoom.pos=NULL, min.value=NULL, max.value=NULL,  
            point.type=16, by.fold=0, with.size=TRUE, with.height=FALSE, 
            point.scale=1, inside.pos=NULL, outside.pos=NULL, outline=TRUE) 
{
   if(is.null(zoom.data) || is.null(data.col) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(is.null(point.scale) || point.scale < 1) point.scale <- 1;
    if(is.null(min.value) || is.null(max.value))
        stop("min.value and max.value must be defined.\n")

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                        inside.pos, outside.pos);
    out.pos <- boundary[1]; 
    in.pos  <- boundary[2];

    #   plot height inside of plot area
    #   =============================================================
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    scatter.values <- as.numeric(zoom.data[, data.col]);   
    if(with.height == TRUE) {
        plot.values <- RCircos.Adjust.Scatter.Values(scatter.values, 
                    max.value, min.value, out.pos-in.pos, 
                    RCircos.Par$sub.track );
    } else { 
        plot.values <- rep((out.pos-in.pos)/2, nrow(zoom.data)); 
    }
    point.height <- in.pos + plot.values;

    #   plot size
    #   ===========================================================
    if(with.size == TRUE)
    {
        point.size <- scatter.values/(max.value-min.value)*10;
    } else {point.size <- rep(point.scale, nrow(zoom.data)); }


    #   plot colors
    #   ===========================================================
    if(by.fold>0) {
        point.colors <- rep("black", nrow(zoom.data));
        red.rows <- which(scatter.values>by.fold);
        if(length(red.rows)>0) point.colors[red.rows] <- "red"; 
        
        blue.rows <- which(scatter.values <= -by.fold);
        if(length(blue.rows)>0) point.colors[blue.rows] <- "blue"; 
    } else {    
        point.colors <- RCircos.Get.Plot.Colors(zoom.data, 
                            RCircos.Par$scatter.color);
    }
                            
    #   Points locations
    #   =============================================================
    #
    point.pos <- RCircos.Zoom.Single.Plot.Positions(zoom.data, zoom.pos)
   
    #   Outline scatter plot at zoom area
    #   =============================================================
    #
    if(outline==TRUE) 
    {
        layers <- RCircos.Par$sub.tracks;
        RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, layers);
    }
    
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    for(a.point in seq_len(length(point.pos)))
    {
        points(RCircos.Pos[point.pos[a.point], 1]*point.height[a.point],
                RCircos.Pos[point.pos[a.point], 2]*point.height[a.point],
                pch=point.type, cex=point.size[a.point], 
                col=point.colors[a.point]);
    }
}





#   ==========================================================================
#
#   13. RCircos.Plot.Zoomed.Tiles()
#
#   Plot zoomed in tiles on zoom in area. Most of this function is same as
#   parallel link plot except replacing lines() by polygons
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and layer values for zoom-in genes/rows
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       genomic.cols:   Non-negative integer, total number of columns for
#                       genomic positions. Valid values are 2 or 3.
#       track.num:      Non-negative integer, which track will be plotted
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#   Return value:   None
#   Example:    RCircos.Plot.Zoomed.Scatters(zoom.data, 5, zoom.pos, 
#                       min.value=-10, max.value=10, track.num=1)
#               RCircos.Plot.Zoomed.Scatters(zoom.data, 5, zoom.pos,
#                       min.value=-10, max.value=10, 
#                       outside.pos=2.5, inside.pos=2);
#
RCircos.Plot.Zoomed.Tiles <- function(zoom.data=NULL, track.num=NULL, 
                zoom.pos=NULL, genomic.cols=3, layers=5, border.col=NULL, 
                inside.pos=NULL, outside.pos=NULL, outline=TRUE)
{
    if(is.null(zoom.data) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(genomic.cols != 3) 
        stop("Tile data must have three columns of genomic position.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                                inside.pos, outside.pos);
    out.pos <- boundary[1]; in.pos  <- boundary[2];
    track.height <- out.pos - in.pos;

    #   Get color, start and end position for each tile 
    #   ================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    tile.colors <- RCircos.Get.Plot.Colors(zoom.data, RCircos.Par$tile.color); 
    tile.location <- RCircos.Zoom.Paired.Plot.Positions(zoom.data, zoom.pos);

    #   Top and bottome locations of each tile
    #   =======================================================
    #
    layer.height <- track.height/layers;
    zoomed.layers <- RCircos.Get.Plot.Layers(zoom.data, genomic.columns=3);
    polygon.bot <- (zoomed.layers-1)*layer.height + in.pos;
    polygon.top <- zoomed.layers*layer.height + in.pos;
    
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    if(outline == TRUE)
    { RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, layers); }

    for(a.tile in seq_len(nrow(zoom.data)))
    {
        first <- tile.location[a.tile, 1];
        last  <- tile.location[a.tile, 2];
        
        polygon.x <- c(RCircos.Pos[first:last, 1]*polygon.bot[a.tile], 
                       RCircos.Pos[last:first, 1]*polygon.top[a.tile]);
        polygon.y <- c(RCircos.Pos[first:last, 2]*polygon.bot[a.tile], 
                       RCircos.Pos[last:first, 2]*polygon.top[a.tile]);
        polygon(polygon.x, polygon.y, col=tile.colors[a.tile]); 
    }

}



#   ==========================================================================
#
#   14. RCircos.Plot.Zoomed.Ideogram.Ticks()
#
#       Draw chromosome ticks and lables along outside of zoom area. Short 
#       ticks will use one track height and labels will start from one track
#       outside of zoom area unless the inside.pos and outside.pos are used
#       instead of track number.
#       
#       Arguments:
#
#       zoom.info:      A vector contains chromosome name, start and optional
#                       end position to be zoomed in.
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       tick.interval:  Non-negative numeric, length in million base pairs 
#                       between two chromosome ticks.
#       track.num:      Non-negative integer, which track will be plotted
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#       Returned value: None
#       Example:        zoom.pos <- c(1:125000)
#                       RCircos.Plot.Zoomed.Ideogram.Ticks(zoom.pos, )
#
#   

RCircos.Plot.Zoomed.Ideogram.Ticks <- function(zoom.info=NULL, track.num=NULL,
    zoom.pos=NULL, tick.interval=5, inside.pos=NULL, outside.pos=NULL)
{
    if(is.null(zoom.info) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(tick.interval <= 0) stop("Tick interval must be > 0")

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                                inside.pos, outside.pos);
    out.pos <- boundary[1]; 
    in.pos  <- boundary[2];

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Tick locations from interval data width 
    #   ========================================================
    interval <- tick.interval*1000000;
    chr.start <- as.numeric(zoom.info[2]);
    chr.end   <- as.numeric(zoom.info[3]);

    #   Move start position to next whole interval
    to.pass <- interval - (chr.start %% interval)
    chr.start <- chr.start + to.pass;
    first.name <- chr.start/interval;

    tick.starts <- seq(chr.start, chr.end, by=interval)
    total.ticks <- length(tick.starts);

    zoom.data <- data.frame(rep(zoom.info[1], total.ticks), tick.starts)
    tick.pos  <- RCircos.Zoom.Single.Plot.Positions(zoom.data, zoom.pos)
    mid.tick.pos <- round((tick.pos[2]-tick.pos[1])/2, digits=0);
    mid.tick.len <- in.pos + (out.pos-in.pos)/2;

    #   line where all ticks start
    line.x <- RCircos.Pos[zoom.pos, 1]*in.pos;
    line.y <- RCircos.Pos[zoom.pos, 2]*in.pos;
    lines(line.x, line.y, col="black");

    text.loc <- rep(4, total.ticks);
    text.loc[which(tick.pos > (nrow(RCircos.Pos)/2))] <- 2;

    for(a.tick in seq_len(total.ticks))
    {
        #   long tick
        lines(c(RCircos.Pos[tick.pos[a.tick],1]*in.pos, 
                RCircos.Pos[tick.pos[a.tick],1]*out.pos),
              c(RCircos.Pos[tick.pos[a.tick],2]*in.pos, 
                RCircos.Pos[tick.pos[a.tick],2]*out.pos),
                col="black");

        #   tick label
        tick.label <- paste0((first.name+a.tick-1)*tick.interval, "MB");
        text(RCircos.Pos[tick.pos[a.tick],1]*out.pos, 
             RCircos.Pos[tick.pos[a.tick],2]*out.pos, 
             tick.label, srt=RCircos.Pos[tick.pos[a.tick], 3],
             offset=0.25, pos=text.loc[a.tick], cex=RCircos.Par$text.size);

        #   short tick
        if(a.tick == total.ticks) next;
        short.tick.pos <- tick.pos[a.tick] + mid.tick.pos;
        lines(c(RCircos.Pos[short.tick.pos,1]*in.pos, 
                RCircos.Pos[short.tick.pos,1]*mid.tick.len),
              c(RCircos.Pos[short.tick.pos,2]*in.pos, 
                RCircos.Pos[short.tick.pos,2]*mid.tick.len),
                col="black");
        
    }
}



#
#   ==========================================================================
#
#   15. RCircos.Plot.Zoomed.Polygons()
#
#       Plot zoomed polygons inside of zoom area. Polygons can be drawn as:
#       1). from top to bottom, 2). from bottom to top, or 3). from middle 
#       to both top and bottom. The difference between zoomed histogram plot 
#       and the second polygon plot is the polygon width. 
#       is 
#       
#       Arguments:
#
#       zoom.data:      A data frame containing genomic positions and
#                       plot values for polygons
#       data.col:       Non-negative integer, which column is labels
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       tick.interval:  Non-negative numeric, length in million base pairs 
#                       between two chromosome ticks.
#       track.num:      Non-negative integer, which track will be plotted
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#       Returned value: None
#       Example:        zoom.pos <- c(1:125000)
#                       RCircos.Plot.Zoomed.Ideogram.Ticks(zoom.pos, )
#
RCircos.Plot.Zoomed.Polygons <- function(zoom.data=NULL, data.col=4, 
                track.num=NULL, zoom.pos=NULL, genomic.cols=3, 
                min.value=NULL, max.value=NULL, border.col=NULL, 
                inside.pos=NULL, outside.pos=NULL, outline=TRUE)
{
    if(is.null(zoom.data) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(genomic.cols != 3) 
        stop("Polygon data must have three columns of genomic position.\n");
    if(data.col < 4) 
        stop("Data for polygon height must be after genomic positions.\n");
        
    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                                inside.pos, outside.pos);
    out.pos <- boundary[1]; in.pos  <- boundary[2];
    track.height <- out.pos - in.pos;

    #   Get color, start and end position for each tile 
    #   ================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    polygon.colors <- RCircos.Get.Plot.Colors(zoom.data, RCircos.Par$tile.color); 
    polygon.location <- RCircos.Zoom.Paired.Plot.Positions(zoom.data, zoom.pos);

    #   Top and bottom locations of each tile
    #   =======================================================
    #
    data.heights <- as.numeric(zoom.data[,data.col])
    if(is.null(min.value) || is.null(max.value)) {
        min.value <- min(data.heights);
        max.value <- max(data.heights);
    }
    
    polygon.loc <- RCircos.Get.Polygon.Height(zoom.data[, data.col],
        min.value, max.value, inside.pos=in.pos, outside.pos=out.pos);
    bottom <- polygon.loc[,1];
    top    <- polygon.loc[,2];

    if(is.null(border.col)) 
        border.col <- rep(NA, nrow(zoom.data))
    else border.col <- rep(border.col, nrow(zoom.data))
    
    RCircos.Pos <- RCircos.Get.Plot.Positions();

    if(outline == TRUE)
    { 
        if(min.value >= 0 || max.value <= 0) {
            layers <- RCircos.Par$sub.tracks;
        } else { layers <- 10; }   
        RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, num.layers=layers); 
    }

    for(a.row in seq_len(nrow(zoom.data)))
    {
        a.start <- polygon.location[a.row, 1];
        a.end   <- polygon.location[a.row, 2];
        
        polygon.x <- c(RCircos.Pos[a.start:a.end, 1]*top[a.row], 
                        RCircos.Pos[a.end:a.start, 1]*bottom[a.row]);
        polygon.y <- c(RCircos.Pos[a.start:a.end, 2]*top[a.row], 
                        RCircos.Pos[a.end:a.start, 2]*bottom[a.row]);

        polygon(polygon.x, polygon.y, col=polygon.colors[a.row], 
                        border=border.col);
    }

}




#   ==========================================================================
#
#   16. RCircos.Zoom.Area.Outline()
#
#   Draw outline for zoom area with lines of sub-tracks. This function is
#   mainly for internal use.
#
#   Arguments:
#   
#       outside.pos:    Non-negative numeric, outside position of plot area
#       inside.pos:     Non-negative numeric, inside position of plot area
#       num.layers:     Non-negative integer, total number of sub tracks
#       zoom.positions: Non-negative numeric vector, the index of RCircos 
#                       plot position
#       fill.col:       Character vector of color names
#
#   Returned value:     None
#   Example:            RCircos.Zoom.Area.Outline(zoom.pos=NULL, 
#                           inside.pos=2.5, outside.pos=3, 
#                           num.layers=5, fill.col="yellow")                
#
RCircos.Zoom.Area.Outline <- function(zoom.pos=NULL, inside.pos=NULL, 
        outside.pos=NULL, num.layers=5, fill.col="white")
{
    if(is.null(outside.pos) || is.null(inside.pos) || is.null(zoom.pos))
        stop("Missing arguments for RCircos.Zoom.Area.Outline() ")

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    
    start.pos <- zoom.pos[1];
    end.pos <- zoom.pos[length(zoom.pos)];
    
    polygonX <- c(RCircos.Pos[start.pos:end.pos,1]*outside.pos,
                  RCircos.Pos[end.pos:start.pos,1]*inside.pos);
    polygonY <- c(RCircos.Pos[start.pos:end.pos,2]*outside.pos,
                  RCircos.Pos[end.pos:start.pos,2]*inside.pos);
    polygon(polygonX, polygonY, col=fill.col);
    
    if(num.layers > 1) {
        if(fill.col=="white") { 
            grid.col <- RCircos.Par$grid.line.color;
        } else { grid.col <- fill.col; }
    
        sub.height <- (outside.pos - inside.pos)/ num.layers;
        for(a.line in seq_len(num.layers)[-num.layers])
        {
            line.height <- inside.pos + sub.height*a.line;
            lines(RCircos.Pos[start.pos:end.pos, 1] * line.height, 
                  RCircos.Pos[start.pos:end.pos, 2] * line.height,
                  col=grid.col);  
        }
    }
}



#   =========================================================================
#
#   17. RCircos.Clear.Zoom.Area
#
#   Erase everything plotted on zoomed area for a data track. Note: marking
#   area cannot be erased since it spans from chromosome highlight to zoom 
#   area and erase this area may remove part or all of a chromosome name.
#
#   Arguments:
# 
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       inside.pos:     Non-negative numeric, inside position of plot area
#       outside.pos:    Non-negative numeric, outside position of plot area
#
#   Returned value:     None
#   Example:            zoom.position <- c(10000:20000)
#                       RCircos.Clear.Zoom.Area(zoom.pos=zoom.position, 
#                                   inside.pos=2.5, outside.pos=3)
#

RCircos.Clear.Zoom.Area <- function(zoom.pos=NULL, track.num=NULL, 
                            inside.pos=NULL, outside.pos=NULL)
{
    if(is.null(zoom.pos))
        stop("Zoom in position must be defined.\n")

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                        inside.pos, outside.pos);
    out.pos <- boundary[1];
    in.pos  <- boundary[2];

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    start.pos <- zoom.pos[1]-10;
    end.pos <- zoom.pos[length(zoom.pos)]+10;

    polygonX <- c(RCircos.Pos[start.pos:end.pos ,1]*out.pos,
                    RCircos.Pos[end.pos:start.pos ,1]*in.pos);
    polygonY <- c(RCircos.Pos[start.pos:end.pos ,2]*out.pos,
                    RCircos.Pos[end.pos:start.pos ,2]*in.pos);
    polygon(polygonX, polygonY, col="white", border="white");
}




#   =========================================================================
#
#   18. RCircos.Zoom.Single.Plot.Positions()
#
#   Scale the original plot locations to plot the data on zoom area. This is
#   used to zoom the plot position of data that require only single genomic 
#   positions such as the heatmap, histogram, points and is not suitable for 
#   data with two genomic positions such as tiles, polygons, parallele lines. 
#
#   Arguments:
#
#       zoom.data:  data frame, data to be plotted onto zoom area
#       zoom.pos:   Non-negative numeric vector, index of RCircos base plot
#                   positions
#
#   Returned value: Non-negative numeric vector, index of RCircos base plot
#                   positions
#
#   Example:    data(RCircos.Scatter.Data.RData)
#               data.rows <- which(RCircos.Scatter.Data$Chromosome=="chr11")
#               zoom.data <- RCircos.Scatter.Data[data.rows[11:21], ];
#               zoom.range <- RCircos.Get.Zoom.Range(zoom.data, 3)
#               zoom.pos <- RCircos.Set.Zoom.Plot.Positions(zoom.range, 
#                               total.genes=11)
#               plot.location <- RCircos.Zoom.Original.Plot.Positions(zoom.data,
#                                        zoom.pos)
#

RCircos.Zoom.Single.Plot.Positions <- function(zoom.data=NULL, zoom.pos=NULL)
{
    if(is.null(zoom.data) || is.null(zoom.pos))
        stop("Missing argument in RCircos.Zoom.Original.Plot.Positions()")

    #   Orignal plot locations
    original.loc <- RCircos.Get.Single.Point.Positions(zoom.data, 
                            genomic.columns=2);
    original.loc <- as.numeric(original.loc$Location);
    original.len <- max(original.loc) - min(original.loc);
    
    zoom.factor <- length(zoom.pos)/original.len;
    if(zoom.factor < 1) 
        stop("Length of orignal plot location is greater than zoom position");
   
    zoomed.loc <- original.loc - min(original.loc);
    zoomed.loc <- zoomed.loc * zoom.factor + zoom.pos[1];
    zoomed.loc <- round(zoomed.loc, digits=0);

    return (zoomed.loc);
}




#   =========================================================================
#
#   19. RCircos.Zoom.Paired.Plot.Positions()
#
#   Scale the original paired locations to plot the data on zoom area. This
#   function is used to zoom the plot position of data that require two 
#   genomic positions such as tiles, polygons, and parallele lines. 
#
#   Arguments:
#
#       zoom.data:  data frame, data to be plotted onto zoom area. Must
#                   have three columns for genomic positions
#       zoom.pos:   Non-negative numeric vector, index of RCircos base plot
#                   positions
#
#   Returned value: Non-negative numeric matrix with two columns represents
#                   index of RCircos base plot positions
#
#   Example:    data(RCircos.Tile.Data.RData)
#           data.rows <- which(RCircos.Scatter.Data$Chromosome=="chr11")
#           zoom.data <- RCircos.Scatter.Data[data.rows[11:21], ];
#           zoom.range <- RCircos.Get.Zoom.Range(zoom.data, 3)
#           zoom.pos <- RCircos.Set.Zoom.Plot.Positions(zoom.range, 
#                               total.genes=11)
#           plot.location <- RCircos.Zoom.Paired.Plot.Positions(zoom.data,
#                                        zoom.pos)
#

RCircos.Zoom.Paired.Plot.Positions <- function(zoom.data=NULL, zoom.pos=NULL)
{
    if(is.null(zoom.data) || is.null(zoom.pos))
        stop("Missing argument in RCircos.Zoom.Original.Plot.Positions()")

    RCircos.Validate.Genomic.Data(genomic.data=zoom.data, 
            plot.type="plot", genomic.columns=3)

    #   Convert genomic positions from three columns to two columns 
    #   Note: genomic positions of different data points may overlap
    #
    genome.loc <- c(zoom.data[,2], zoom.data[,3]);
    loc.order <- order(genome.loc);
    chromosomes <- rep(zoom.data[1,1], length(genome.loc));

    plot.data <- data.frame(chromosomes, genome.loc);
    
    #   Get zoomed position for two columns genomic position
    #
    plot.loc <- RCircos.Get.Single.Point.Positions(plot.data, 
                            genomic.columns=2);
    original.loc <- as.numeric(plot.loc$Location);
    original.len <- max(original.loc) - min(original.loc);
    
    zoom.factor <- length(zoom.pos)/original.len;
    if(zoom.factor < 1) 
        stop("Length of orignal plot location is greater than zoom position");
   
    zoom.loc <- original.loc - min(original.loc);
    zoom.loc <- floor(zoom.loc * zoom.factor) + zoom.pos[1];
    zoom.loc <- zoom.loc[order(loc.order)];
    zoom.loc <- round(zoom.loc, digits=0);
    
    #   Covert the vector to two column matrix
    start.loc <- c(1:nrow(zoom.data))
    end.loc   <- start.loc + nrow(zoom.data)
    zoomed.location <- cbind(zoom.loc[start.loc], zoom.loc[end.loc])
    return (zoomed.location);
}




#   =========================================================================
#
#   20. RCircos.Plot.Zoomed.Area()
#   Plot continuous lines from one point to next point on zoom in area
#
#       zoom.data:      A data frame containing genomic positions, gene 
#                       names, and plot values for zoom-in genes/rows
#       data.col:       Non-negative integer, which column is labels
#       zoom.pos:       Non-negative numeric vector, the index of RCircos 
#                       plot position
#       line.width:     Non-negative integet, line width parameter
#       min.value:      Numeric, minimum value ifor point height
#       max.value:      Numeric, maximum value ifor point height
#       track.num:      Non-negative integer, which track will be plotted
#       outside.pos:    Non-negative numeric, the far location relative
#                       to center of plot area
#       inside.pos:     Non-negative numeric, the cloase location relative  
#                       to center of plot area
#
#   Return value:   None
#   Example:        RCircos.Plot.Zoomed.Cotinuous.Lines(zoom.data, 5, zoom.pos,
#                                           4, -5, 5, 1);
#                   RCircos.Plot.Zoomed.Cotinuous.Lines(zoom.data, data.col=5, 
#                           zoom.pos, line.width=4, min.value=-5, max.value=5, 
#                           outside.pos=5, inside.pos=4.5);
#

RCircos.Plot.Zoomed.Area <- function(zoom.data=NULL, plot.type="mountain",
        data.col=NULL, track.num=NULL, zoom.pos=NULL, min.value=NULL, 
        max.value=NULL, area.color="gray", border.col="black",
        inside.pos=NULL, outside.pos=NULL, outline=TRUE) 
{
    if(is.null(zoom.data) || is.null(data.col) || is.null(zoom.pos)) 
        stop("Missing argument in RCircos.Plot.Zoomed.Lines().\n");
    if(! plot.type %in% c("mountain","curtain", "band")) 
        stop("Plot type must be either 'mountain', 'curtain', or 'band'.\n");
    
    if(is.null(min.value) || is.null(max.value))
        stop("min.value and max.value must be defined.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side="out", 
                          inside.pos, outside.pos);
    out.pos <- boundary[1]; in.pos  <- boundary[2];
    track.height <- out.pos - in.pos;

    #   Location of each data point
    #   =============================================================
    #
    point.pos <- RCircos.Zoom.Single.Plot.Positions(zoom.data, zoom.pos)
    area.color <- RCircos.Get.Plot.Colors(zoom.data, area.color);

    #   area height
    #   =============================================================
    #
    if(plot.type != "band")
    {
        area.values <- as.numeric(zoom.data[, data.col]);
        point.height <- RCircos.Get.Data.Point.Height(area.values, min.value, 
            max.value, plot.type="points", track.height=track.height);
        if(plot.type == "mountain")
        {
            area.bot <- rep(in.pos, length(point.height));
            area.top <- in.pos + point.height;
        } else {
            area.top <- rep(out.pos, length(point.height));
            area.bot <- out.pos - point.height;
        }
    } else {
        area.values <- as.matrix(zoom.data[, data.col]);
        differ <- sum(area.values[,1] > area.values[, 2])
        if(differ == 0) {
            bottom <- 1; top <- 2; 
        } else if(differ == nrow(area.values)) {
            top <- 1; bottom <- 2;
        } else { stop("Band top cross area bottom.\n")}

        band.top <- as.numeric(area.values[, top]);
        area.top <- RCircos.Get.Data.Point.Height(band.top, min.value, 
            max.value, plot.type="points", track.height=track.height);
        area.top <- area.top + in.pos;
        
        band.bot <- as.numeric(area.values[, bottom]);
        area.bot <- RCircos.Get.Data.Point.Height(band.bot, min.value, 
            max.value, plot.type="points", track.height=track.height);
        area.bot <- area.bot + in.pos;
    }

    #   Outline at zoom area
    #   =============================================================
    #
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    if(outline == TRUE)
    {
        RCircos.Par <- RCircos.Get.Plot.Parameters();
        layers <- RCircos.Par$sub.tracks;
        RCircos.Zoom.Area.Outline(zoom.pos, in.pos, out.pos, layers);
    }

    polygon.x <- c(RCircos.Pos[point.pos, 1]*area.top, 
                    RCircos.Pos[rev(point.pos), 1]*rev(area.bot));
    polygon.y <- c(RCircos.Pos[point.pos, 2]*area.top, 
                    RCircos.Pos[rev(point.pos), 2]*rev(area.bot));

    polygon(polygon.x, polygon.y, col=area.color, border=border.col);
}


#   End of RCircosZoomPlot.R
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
