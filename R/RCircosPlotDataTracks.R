#    
#    Functions of RCircos data plot on different tracks
#
#		 1.		RCircos.Gene.Connector.Plot()
#		 2.		RCircos.Gene.Name.Plot()
#		 3.		RCircos.Heatmap.Plot()
#		 4.		RCircos.Histogram.Plot()
#		 5.		RCircos.Line.Plot()
#		 6.		RCircos.Scatter.Plot()
#		 7.		RCircos.Tile.Plot()
#		 8.		RCircos.Link.Plot()
#		 9.		RCircos.Ribbon.Plot()
#		10.		RCircos.Clear.Track()
#		11.		RCircos.Track.Outline()
#		12.		RCircos.Vertical.Line.Plot()
#		13.		RCircos.Point.Plot()
#		14.		RCircos.Parallele.Link.Plot()
#		15.		RCircos.Polygon.Plot()
#		16.		RCircos.Area.Highlight()
#		17.		RCircos.Customized.Shape.Plot()
#		18.		RCircos.Get.Start.End.Locations()
#		19.		RCircos.Adjust.Scatter.Values()
#		20.		RCircos.Area.plot()
#		21.		RCircos.Customized.Connection.Plot()
#
#   New arguments have been added for version 1.2 to allow advanced users
#
#   1.  to define track locations directly
#   2.  accept different number (2 or 3) of columns for genomic position
#   3.  pre-sorting is no required for genomic data
#
#   Old arguments will be first and in original order so new arguments
#   could be ignored for use of old version
#
#    Last modified on August 21, 2017
#    ________________________________________________________________________
#    <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>





#     ==========================================================================
#
#   1.  RCircos.Gene.Connector.Plot()
#
#   Draw connectors between chromosome ideogram and gene labels. The plot data 
#   (connectData) has two paired point locations on two neighbour tracks 
#   (chromosomes ideogram and gene label track). The first column is for outer
#   points, and the second column is inner points. 
#
#   The points are sorted by relative positions on their chromosomes and are 
#   held in the last two columns of connectData.
#
#   Arguments:
#
#       genomic.data:   Data frame with the first four columns for chromosome
#                       name, start and end position, and names of genes.
#       track.num:      Non-negative integer, number of the track from 
#                       chromosome ideogram.
#       side:           Character vector, must be either "in" or "out".
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be either 3 or 2. 
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#   Return value:    None
#
#   Example:    RCircos.Gene.Connector.Plot(genomic.data, 1, "in")
#               RCircos.Gene.Connector.Plot(genomic.data, 0.9, 0.85)
#

RCircos.Gene.Connector.Plot <- function(genomic.data=NULL, 
            track.num=NULL, side="in", inside.pos=NULL, 
            outside.pos=NULL, genomic.columns=3, is.sorted=FALSE)
{
    if(is.null(genomic.data)) 
        stop("Genomic data missing for RCircos.Gene.Connector.Plot().\n");
    
    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, erase.area=FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #    Construct Connector data from gene name data
    #    =======================================================
    #
    geneData    <- RCircos.Get.Single.Point.Positions(genomic.data, 
                             genomic.columns);
    labelData   <- RCircos.Get.Gene.Label.Locations(geneData, 
                            genomic.columns, is.sorted);
    connectData <- data.frame(labelData$Location, labelData$LabelPosition);

    #    Switch the columns of genomic location and label 
    #    location for inside or outside
    #    =================================================
    #
    if(outerPos < RCircos.Par$chr.ideo.pos) {
        genomicCol <- ncol(connectData) - 1;
        labelCol <- ncol(connectData);
    } else {
        genomicCol <- ncol(connectData);
        labelCol <- ncol(connectData) - 1;
    }

    #    Heights for the two vertical lines of connectors and
    #    the horizontal line range
    #    ====================================================
    #
    vHeight <- round((outerPos-innerPos)/10, digits=4);
    hRange <- outerPos - innerPos - 2*vHeight;

    topLoc <- outerPos - vHeight;
    botLoc <- innerPos + vHeight;

    #    Connector colors
    #    ===============================================
    #
    lineColors <- RCircos.Get.Plot.Colors(labelData, RCircos.Par$text.color);

    #    Plot Connectors
    #    ===============================================
    #
    chroms <- unique(connectData[,1]);
    for(aChr in seq_along(chroms))
    {
        chrRows <- which(connectData[,1]==chroms[aChr]);
        total <- length(chrRows);

        for(aPoint in seq_len(total))
        {
            p1 <- connectData[chrRows[aPoint], genomicCol];
            p2 <- connectData[chrRows[aPoint], labelCol];

            #    draw top vertical line
            #    ======================================
            #
            lines(c(RCircos.Pos[p1, 1]*outerPos, RCircos.Pos[p1,  1]*topLoc),
                    c(RCircos.Pos[p1,2]*outerPos, RCircos.Pos[p1, 2]*topLoc),
                    col=lineColors[chrRows[aPoint]]);

            #    draw bottom vertical line
            #    ======================================
            #
            lines(c(RCircos.Pos[p2, 1]*botLoc, RCircos.Pos[p2, 1]*innerPos),
                    c(RCircos.Pos[p2,2]*botLoc, RCircos.Pos[p2, 2]*innerPos),
                    col=lineColors[chrRows[aPoint]]);

            #    draw horizontal line
            #    ======================================
            #
            lines(c(RCircos.Pos[p1,  1]*topLoc, RCircos.Pos[p2, 1]*botLoc),
                    c(RCircos.Pos[p1, 2]*topLoc, RCircos.Pos[p2, 2]*botLoc),
                    col=lineColors[chrRows[aPoint]]);
        }
    }
}




#     ==========================================================================
#
#   2.  RCircos.Gene.Name.Plot()
#
#    Label genes beside of track. This is only suitable for small number of 
#    labels. When cex=0.4, each character of the label will occupy about 5000 
#    units. This is the best visualization for a 8x8 inch image.
#
#    Arguments:
#
#       gene.data:      A data frame returned from the function call to
#                       RCircos.Get.Gene.Label.Locations(genomic.data). 
#                       The genomic.data has three leading columns for  
#                       genomic positions followed by gene names.
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be either 3 or 2. 
#       name.col:       Non-negative integer, number of column in gene.data 
#                       for data to be plotted (gene names)
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram.
#       side:           Character vector, must be either "in" or "out"
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#    Return value:    None
#
#    Example:   RCircos.Gene.Label(gene.data, name.col=4, track.num=3, "in")
#               RCircos.Gene.Label(gene.data, inside.pos=0.7, outside.pos=0.8) 
#
#

RCircos.Gene.Name.Plot <- function(gene.data=NULL, name.col=NULL, 
            track.num=NULL, side="in", inside.pos=NULL, outside.pos=NULL,
            genomic.columns=3, is.sorted=FALSE)
{
    if(is.null(gene.data)) 
        stop("Genomic data missing in RCircos.Gene.Name.Plot().\n");
    if(is.null(genomic.columns)) 
        stop("Missing number of columns for genomic position.\n");
     if( is.null(name.col) ||name.col <= genomic.columns)
        stop("Data column must be ", genomic.columns+1, " or bigger.\n"); 

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    textColors <- RCircos.Get.Plot.Colors(gene.data, RCircos.Par$text.color);

    #    Convert raw data to plot data. The raw data will be validated
    #    first during the conversion
    #    =============================================================
    #
    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);   
    gene.data <- RCircos.Get.Single.Point.Positions(gene.data, 
                             genomic.columns);
    gene.data <- RCircos.Get.Gene.Label.Locations(gene.data,  genomic.columns,
                            is.sorted);
    #    Label positions
    #    =============================================================
    #
    rightSide <- nrow(RCircos.Pos)/2;
    thePoints <- as.numeric(gene.data[, ncol(gene.data)]);

    if(side=="in") {
        labelPos <- boundary[1]; 
        textSide <- rep(4, nrow(gene.data));
        textSide[thePoints <= rightSide] <- 2;    
    } else {
        labelPos  <- boundary[2];
        textSide <- rep(2, nrow(gene.data));
        textSide[thePoints <= rightSide] <- 4; 
    }

    #    Plot labels
    #    =============================================================
    #
    for(aText in seq_len(nrow(gene.data)))
    {
        geneName <- as.character(gene.data[aText, name.col]);
        rotation <- RCircos.Pos$degree[thePoints[aText]];

        text(RCircos.Pos[thePoints[aText],1]*labelPos,
             RCircos.Pos[thePoints[aText],2]*labelPos,
            label=geneName, pos=textSide[aText], 
            cex=RCircos.Par$text.size, srt=rotation, 
            offset=0, col=textColors[aText]);
    }
}




#   =========================================================================
#
#   3.  RCircos.Heatmap.Plot()
#
#   Draw one track of heatmap with blue and red colors. The first four columns 
#   of headmap data must be chromosome, chromStart, chromEnd, and gene names.
# 
#   Arguments:
#
#       heatmap.data:   A data frame with returned from function call to
#                       RCircos.Get.Plot.Data(genomic.data, plot.type). 
#                       Heatmap data must have leading columns for gene 
#                       position and gene names.
#       data.col:       Non-negative integer, number of column in heatmap.data 
#                       for data to be plotted
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram.
#       side:           Character vector, must be either "in" or "out"
#       min.value:      Numeric, minimum value for heatmap scale
#       max.value:      Numeric, maximum value for heatmap scale
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be either 3 or 2.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#    Return value:    None
#
#    Example:   RCircos.Heatmap.Plot(heatmap.data, 3, 3, "in")
#               RCircos.Heatmap.Plot(heatmap.data, data.col=3, 
#                               inside.pos=0.8, outside.pos=0.9)
#

RCircos.Heatmap.Plot <- function(heatmap.data=NULL, data.col=NULL, 
    track.num=NULL, side=c("in", "out"), min.value=NULL, max.value=NULL, 
    inside.pos=NULL, outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(heatmap.data)) 
        stop("Genomic data missing in RCircos.Heatmap.Plot().\n");
    if(is.null(genomic.columns)) 
        stop("Missing number of columns for genomic position.\n");
     if( is.null(data.col) || data.col <= genomic.columns)  
        stop("Data column must be ", genomic.columns+1, " or bigger.\n"); 

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos <- boundary[2];

    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Colors for different data values
    #   ===========================================================
    #
    colorMap <- RCircos.Get.Heatmap.Color.Scale(RCircos.Par$heatmap.color);

    if(is.null(min.value) || is.null(max.value)) 
    {
        columns <- (genomic.columns+2):ncol(heatmap.data);
        min.value <- min(as.matrix(heatmap.data[, columns]));
        max.value <- max(as.matrix(heatmap.data[, columns]));
    }
    colorLevel  <- seq(min.value, max.value, length=length(colorMap));

    #   Each heatmap cell starts from data start location and span for
    #   heatmap.width. Make sure each one will be in the range of thire 
    #   chromosome
    #   ===============================================================
    #
    heatmap.data <- RCircos.Get.Single.Point.Positions(heatmap.data, 
                                genomic.columns);
    plotLocations <- RCircos.Get.Start.End.Locations(heatmap.data, 
                                RCircos.Par$heatmap.width);

    #   outline of chromosomes. No lines inside.
    #   ===============================================================
    #
    chromosomes <- unique(as.character(RCircos.Cyto$Chromosome));
    outlineColors <- rep("white", length(chromosomes));
    RCircos.Track.Outline(outerPos, innerPos, num.layers=1, 
                chrom.list=chromosomes, track.colors=outlineColors);

    #   Plot heatmap for each gene.
    #   ===============================================================
    #
    heatmapValues <- as.numeric(heatmap.data[, data.col]);
    for(aPoint in 1:length(heatmapValues))
    {
        theLevel <- which(colorLevel >= heatmapValues[aPoint]);
        cellColor <- colorMap[min(theLevel)];
        
        theStart <- plotLocations[aPoint, 1];
        theEnd   <- plotLocations[aPoint, 2];

        polygonX <- c(RCircos.Pos[theStart:theEnd,1]*outerPos, 
                RCircos.Pos[theEnd:theStart,1]*innerPos);
        polygonY <- c(RCircos.Pos[theStart:theEnd,2]*outerPos, 
                RCircos.Pos[theEnd:theStart,2]*innerPos);
        polygon(polygonX, polygonY, col=cellColor, border=NA);
    }
}




#     =========================================================================
#
#   4.  RCircos.Histogram.Plot()
#
#    Arguments:
#
#       hist.data:      A data frame with returned from Circos.Get.Plot.Data(). 
#                       Histogram data has three leading columns for genomic 
#                       positions.
#       data.col:       Non-negative integer, number of column in heatmap.data 
#                       for data to be plotted
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram.
#       side:           Character vector, must be either "in" or "out"
#       min.value:      Numeric, minimum value for histogram height
#       max.value:      Numeric, maximum value for histogram height
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be either 3 or 2.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#    Return value:    None
#
#    Example:   RCircos.Histogram.Plot(hist.data, 4, 2, "in")
#               RCircos.Histogram.Plot(hist.data, data.col=4, 
#                               inside.pos=0.8, outside.pos=0.9)
#

RCircos.Histogram.Plot <- function(hist.data=NULL, data.col=4, 
    track.num=NULL, side=c("in", "out"), min.value=NULL, 
    max.value=NULL, inside.pos=NULL, outside.pos=NULL,
    genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(hist.data)) 
        stop("Genomic data missing in RCircos.Histogram.Plot().\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];

    if(is.null(genomic.columns) || genomic.columns<2 || genomic.columns>3) 
        stop("Incorrect number of columns for genomic position.\n");
     if( is.null(data.col) || data.col <= genomic.columns)  
        stop("Hist data column must be greater than", genomic.columns, ".\n"); 

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();

    #    Convert raw data to plot data. The raw data will be validated
    #    first during the convertion
    #    ============================================================
    #
    hist.data <- RCircos.Get.Single.Point.Positions(hist.data,
                            genomic.columns);
    locations <- RCircos.Get.Start.End.Locations(hist.data, 
                            RCircos.Par$hist.width)

    #    histgram colors and height
    #    =========================================================
    #
    histColors <- RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color); 

    histValues <- as.numeric(hist.data[, data.col]);
    if(is.null(max.value) || is.null(min.value)) {
        max.value <- max(histValues);
        min.value <- min(histValues);
    } else {
        if(min.value > max.value) stop("min.value > max.value.")
    }
    histHeight <- RCircos.Get.Data.Point.Height(histValues, min.value, 
            max.value, plot.type="points", outerPos-innerPos);
  
    #    Draw histogram
    #    =============================================================
    #
    RCircos.Track.Outline(outerPos, innerPos, RCircos.Par$sub.tracks);

    for(aPoint in seq_len(nrow(hist.data)))
    {
        height <- innerPos + histHeight[aPoint];
        theStart <- locations[aPoint, 1];
        theEnd <- locations[aPoint, 2];

        #    Plot rectangle with specific height for each data point
        #    =========================================================
        #
        polygonX <- c(RCircos.Pos[theStart:theEnd,1]*height, 
                        RCircos.Pos[theEnd:theStart,1]*innerPos);
        polygonY <- c(RCircos.Pos[theStart:theEnd,2]*height, 
                        RCircos.Pos[theEnd:theStart,2]*innerPos);
        polygon(polygonX, polygonY, col=histColors[aPoint], border=NA);
    }
}




#     =========================================================================
#
#   5.  RCircos.Line.Plot()
#
#       This will draw lines between two neibors for all points
# 
#    Arguments:
#
#       line.data:      A data frame returned from function call of 
#                       RCircos.Get.Plot.Data(genomic.data). line.data has  
#                       three leading columns for genomic positions.
#       data.col:       Non-negative integer, number of column in heatmap.data 
#                       for the data to be plotted
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram
#       side:           Character vector, must be either "in" or "out"
#       min.value:      Numeric, minimum value of line data
#       max.value:      Numeric, maximum value of line data
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be either 3 or 2.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#    Return value:    None
#
#    Example:    RCircos.Line.Plot(line.data, 4, 3, "in")
#
#
#    Last revised on June 16, 2015
#
#

RCircos.Line.Plot <- function(line.data=NULL, data.col=4, track.num=NULL, 
        side=c("in", "out"), min.value=NULL, max.value=NULL, 
        inside.pos=NULL, outside.pos=NULL, genomic.columns=3,is.sorted=TRUE)
{
    if(is.null(line.data)) 
        stop("Genomic data missing in RCircos.Line.Plot().\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];

    if(is.null(genomic.columns)) 
        stop("Missing number of columns for genomic position.\n");
     if( is.null(data.col) || data.col <= genomic.columns)  
        stop("Line data column must be ", genomic.columns+1, " or bigger.\n"); 

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #    Convert raw data to plot data. The raw data will 
    #    be validated first during the convertion
    #    ===================================================
    #
    line.data <- RCircos.Get.Single.Point.Positions(line.data,
                        genomic.columns);

    pointValues <- as.numeric(line.data[, data.col]);
    if(is.null(min.value) || is.null(max.value)){
        min.value <- min(pointValues);
        max.value <- max(pointValues);
    } else {
        if(min.value > max.value) stop("min.value > max.value.")
    }
    pointHeight <- RCircos.Get.Data.Point.Height(pointValues, min.value, 
            max.value, plot.type="points", outerPos-innerPos);
    pointHeight <- pointHeight + innerPos;

    #     Line colors
    #    ============================================================
    #     
    line.colors <- RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color); 

    #    Start plotting. Line plot is connecting two neighbor points
    #    so no exception catch needed.
    #    ===========================================================
    #
    RCircos.Track.Outline(outerPos, innerPos, RCircos.Par$sub.tracks)

    for(aPoint in seq_len((nrow(line.data)-1)))
    {
        #   chromosome changed....
        if(line.data[aPoint, 1]!= line.data[aPoint+1, 1]) { next;}

        point.one <- line.data[aPoint, ncol(line.data)];
        point.two <- line.data[aPoint+1, ncol(line.data)];

        #    Draw lines
        #    =====================================================
        #
        lines(c(RCircos.Pos[point.one , 1]*pointHeight[aPoint],
                RCircos.Pos[point.two , 1]*pointHeight[aPoint+1]),
              c(RCircos.Pos[point.one , 2]*pointHeight[aPoint],
                RCircos.Pos[point.two , 2]*pointHeight[aPoint+1]),
                col=line.colors[aPoint]);
    }
}




#   =========================================================================
#
#   6.  RCircos.Scatter.Plot()
#
#   Draw one track of scatterplot. 
# 
#   Arguments:
#
#       scatter.data:   A data frame returned from function call of
#                       RCircos.Get.Plot.Data(genomic.data, plot.type). 
#                       The scatter.data has three leading columns for 
#                       genomic positions.
#       data.col:       Non-negative integer, number of column in heatmap.data 
#                       for the data to be plotted
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram
#       side:           Character vector, must be either "in" or "out"
#       by.fold:
#       min.value:      Numeric, minimum value of scatter plot data
#       max.value:      Numeric, maximum value of scatter plot data
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be either 3 or 2.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#   Return value:   None
#
#   Example:    RCircos.Scatter.Plot(scatter.data, 5, 3, "in", 1)
#               RCircos.Scatter.Plot(scatter.data, data.col=5, by.fold=0, 
#                               inside.pos=1.5, outside.pos=1.6)
#

RCircos.Scatter.Plot <- function(scatter.data=NULL, data.col=4, 
        track.num=NULL, side=c("in", "out"), by.fold=0, 
        min.value=NULL, max.value=NULL, inside.pos=NULL, 
        outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(scatter.data)) 
        stop("Genomic data missing in RCircos.Scatter.Plot().\n");

    if(is.null(genomic.columns)) 
        stop("Missing number of columns for genomic position.\n");
    if( is.null(data.col) || data.col <= genomic.columns)  
        stop("Data column must be ", genomic.columns+1, " or bigger.\n"); 
    if(by.fold<0) stop("Value of by.fold cannot be negative.\n")
    
    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];
    
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #    Convert raw data to plot data then adjust the data value to 
    #    avoid overflow. The raw data will be validated first 
    #    =============================================================
    scatter.values <- as.numeric(scatter.data[,data.col]);
    scatter.data <- RCircos.Get.Single.Point.Positions(scatter.data,
                      genomic.columns);

    #   scatter colors 
    #   =====================================================
    if(by.fold>0) {
        scatter.colors <- rep("black", nrow(scatter.data));
        red.rows <- which(scatter.values>by.fold);
        if(length(red.rows)>0) scatter.colors[red.rows] <- "red"; 
        
        blue.rows <- which(scatter.values <= -by.fold);
        if(length(blue.rows)>0) scatter.colors[blue.rows] <- "blue"; 
    } else {    
        scatter.colors <- RCircos.Get.Plot.Colors(scatter.data, 
                            RCircos.Par$scatter.color);
    }
    
    if(is.null(max.value) || is.null(min.value))
    { 
        max.value <- max(scatter.values); 
        min.value <- min(scatter.values) 
    }
    plot.values <- RCircos.Adjust.Scatter.Values(scatter.values, 
        min.value=min.value, max.value=max.value, 
        track.height=outerPos-innerPos, subtrack=RCircos.Par$sub.track)
    plot.height <- innerPos + plot.values;

    #    Start plotting
    #    ============================================================
    RCircos.Track.Outline(outerPos, innerPos, RCircos.Par$sub.tracks)

    for(a.point in seq_len(nrow(scatter.data)))
    {
        the.point <- scatter.data[a.point, ncol(scatter.data)];
        points(RCircos.Pos[the.point ,1]*plot.height[a.point], 
                RCircos.Pos[the.point ,2]*plot.height[a.point],
                col=scatter.colors[a.point], pch=RCircos.Par$point.type, 
                cex=RCircos.Par$point.size);
    }
}




#   =========================================================================
#
#   7.  RCircos.Tile.Plot()
#
#   Draws one track of tiles. Note: Tile plot needs genomic position only
#   and data column is not requied.
#
#   Arguments:
#
#       tile.data:      A data frame returned from the function call to 
#                       RCircos.Get.Plot.Data(genomic.data, plot.type). The
#                       tile.data has three columns for genomic positions only
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram.
#       side:           Character vector, must be either "in" or "out"
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be 3.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#   Return value:  None
#
#   Example:    RCircos.Tile.Plot(tile.data, track.num=3, side="in")
#               RCircos.Tile.Plot(tile.data, inside.pos=1.2, outside.pos=1.3)
#

RCircos.Tile.Plot <- function(tile.data=NULL, track.num=NULL, 
                side=c("in", "out"), inside.pos=NULL, outside.pos=NULL,
                genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(tile.data)) 
        stop("Genomic data missing in RCircos.Tile.Plot().\n");

    if(is.null(genomic.columns) || genomic.columns < 3) 
        stop(paste("Genomic position must include chromosome name,",
            "start and end position.\n")); 

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Convert raw data to plot data. The raw data will be validated
    #   first during the conversion
    #   =============================================================
    #
    tile.data <- RCircos.Get.Paired.Points.Positions(tile.data,
                genomic.columns, plot.type="tile");

    #   Assign a layer number to each data point and find the maximum
    #   layer number
    #   ============================================================
    #
    tile.layers <- RCircos.Get.Plot.Layers(tile.data, genomic.columns);
    layer.height <- RCircos.Par$track.height/RCircos.Par$max.layers;
    num.layers <- max(tile.layers);

    if(num.layers>RCircos.Par$max.layers) 
    { 
        if(side=="in")
        {  innerPos <- outerPos - layer.height*num.layers;        
        } else { outerPos <- innerPos + layer.height*num.layers; }

        message(paste("Tiles plot may use more than one track.",
            "Please select correct area for next track if necessary.\n"));
    }

    if(num.layers<RCircos.Par$max.layers) 
    { layer.height <- RCircos.Par$track.height/num.layers; }

    #   Tile colors
    #   =====================
    #
    tile.colors <- RCircos.Get.Plot.Colors(tile.data, RCircos.Par$tile.color); 

    #   Start plotting
    #   =============================================================
    #
    RCircos.Track.Outline(outerPos, innerPos, num.layers);

    the.loc <- ncol(tile.data);
    for(a.row in seq_len(nrow(tile.data)))
    {
        tile.start <- tile.data$LinkStart[a.row]
        tile.end   <- tile.data$LinkEnd[a.row];

        layer.bot <- innerPos + layer.height*(tile.layers[a.row]-1);
        layer.top <- layer.bot + layer.height*0.8;

        polygon.x<- c(RCircos.Pos[tile.start:tile.end,1]*layer.top, 
                RCircos.Pos[tile.end:tile.start,1]*layer.bot);
        polygon.y<- c(RCircos.Pos[tile.start:tile.end,2]*layer.top, 
                RCircos.Pos[tile.end:tile.start,2]*layer.bot);
        polygon(polygon.x, polygon.y, col=tile.colors[a.row]);
    }
}




#   =========================================================================
#
#   8.  RCircos.Link.Plot()
#
#   Draw link lines between two chromosomes. Link line is always inside 
#   of ideogram and data has only paired genomic positions
#  
#   Arguments:
#
#       link.data:      Data frame with paired genomic positions in each row
#       track.num:      Non-negative number, number of the track from 
#                       chromosome ideogram.
#       by.chromosome:  Logical, if true, intra chromosome will be red, 
#                       otherwise random color will be used
#       start.pos:      Non-negative number, scale factor based on chromosome
#                       ideogram position. Must be smaller than 1.
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be 3.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#       lineWidth:      Non-negative numeric vector, width for each link line.
#
#   Example:    RCircos.Link.Plot(link.data, 9, TRUE)
#               RCircos.Link.Plot(link.data, by.chromosome=FALSE, start.pos=0.6)
#

RCircos.Link.Plot <- function(link.data=NULL, track.num=NULL, by.chromosome=FALSE, 
                start.pos=NULL, genomic.columns=3, is.sorted=TRUE, 
                lineWidth=rep(1, nrow(link.data)))
{
    if(is.null(link.data)) stop("Link data missing in RCircos.Link.Plot().\n");
    if(by.chromosome!=TRUE && by.chromosome!=FALSE)
        stop("Error: by.chromosome must be either TRUE or FALSE.\n");
		
    if(length(which(lineWidth < 0)) > 1 ) 
		stop("Line width must be positive.");
	if(length(lineWidth) > 1 && length(lineWidth) != nrow(link.data))
		stop("Length of lineWidth must match rows of link.data");

    RCircos.Par <- RCircos.Get.Plot.Parameters();
    
    if(is.null(start.pos)) {
        locations <- RCircos.Get.Plot.Boundary(track.num, side="in", 
                        inside.pos=NULL, outside.pos=NULL, FALSE);
        line.start <- locations[1];
    } else {
        if(start.pos>=1) stop("Link line must be inside chromosome ideogram");
        line.start <- RCircos.Par$chr.ideo.pos * start.pos;
    }

     if(is.null(genomic.columns) || genomic.columns < 3) 
        stop("Incorrect number of columns for genomic position.\n");
 
    #   Start and  end point for each link line
    #   ================================================
    #
    link.data <- RCircos.Get.Paired.Points.Positions(link.data, 
                    genomic.columns, plot.type="link");

    #   Get link line colors for each pair of locations
    #   ================================================
    # 
    link.colors <- RCircos.Get.Link.Colors(link.data, genomic.columns, 
                        by.chromosome);

    #   Draw link lines for each pair of locations
    #   ==========================================
    #
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    base.positions <- RCircos.Pos[, 1:2]*line.start;
    for(a.link in seq_len(nrow(link.data)))
    {  
        point.one <- as.numeric(link.data$LinkStart[a.link]);
        point.two <- as.numeric(link.data$LinkEnd[a.link]);
        if(point.one > point.two)
        { 
            point.one <- link.data$LinkEnd[a.link];
            point.two <- link.data$LinkStart[a.link];
        }

        P0 <- as.numeric(base.positions[point.one, ]);
        P2 <- as.numeric(base.positions[point.two, ]);
        links <- RCircos.Link.Line(P0, P2); 
        lines(links$pos.x, links$pos.y, type="l", 
            lwd=lineWidth[a.link], col=link.colors[a.link] ); 
    }
}




#   =========================================================================
# 
#   9.  RCircos.Ribbon.Plot()
#
#   Ribbon link plot. Ribbons are wide link line between two chromosomes
#
#   Arguments:
#
#       link.data:      Data frame with paired genomic positions in each row
#       track.num:      Non-negative number,  number of the track from  
#                       chromosome ideogram.
#       by.chromosome:  Logic, if true, intra chromosome will be red, 
#                       otherwise random color will be used
#       twist:          Logic, if the ribbon is twisted
#       start.pos:      Non-negative number, scale factor based on chromosome
#                       ideogram position. Must be smaller than 1.
#
#   Return value:    None
#
#   Example:    RCircos.Ribbon.Plot(ribbon.data, 10, FALSE, FALSE)
#               RCircos.Ribbon.Plot(ribbon.data, start.pos=0.75)
#
#    Last revised on June 11, 2014
#
#

RCircos.Ribbon.Plot <- function(ribbon.data=NULL, track.num=NULL, 
            by.chromosome=FALSE, twist=FALSE, start.pos=NULL, 
            genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(ribbon.data)) 
        stop("Ribbon data missing in RCircos.Ribbon.Plot().\n");
  
    if(by.chromosome!=TRUE && by.chromosome!=FALSE)
        stop("Error: by.chromosome must be either TRUE or FALSE.\n");   
    if(twist!=TRUE &&twist!=FALSE) 
        stop("Error: twist must be either TRUE or FALSE.\n"); 
        
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    
    if(is.null(start.pos)) {
        locations <- RCircos.Get.Plot.Boundary(track.num, side="in", 
                        inside.pos=NULL, outside.pos=NULL, FALSE);
        ribbon.start <- locations[1];
    } else {
        if(start.pos>=1) stop("Link line must be inside of ideogram");
        ribbon.start <- RCircos.Par$chr.ideo.pos * start.pos;
    }

    #   Check chromosome names, start, and end positions
    #   =================================================
    #
    RCircos.Validate.Genomic.Data(ribbon.data, plot.type="link", 
                        genomic.columns=genomic.columns);

    #   Coordinates of the four conner of each ribbon (polygon)
    #   =======================================================
    #
    data.points <- matrix(rep(0, nrow(ribbon.data)*4), ncol=4);
    for(a.link in seq_len(nrow(ribbon.data)))
    {
        data.points[a.link, 1] <- RCircos.Data.Point(
            ribbon.data[a.link, 1], ribbon.data[a.link, 2]);
        data.points[a.link, 2] <- RCircos.Data.Point(
            ribbon.data[a.link, 1], ribbon.data[a.link, 3]);
        data.points[a.link, 3]<- RCircos.Data.Point(
            ribbon.data[a.link, 4], ribbon.data[a.link, 5]);
        data.points[a.link, 4]<- RCircos.Data.Point(
            ribbon.data[a.link, 4], ribbon.data[a.link, 6]);

        if(data.points[a.link, 1]==0 || data.points[a.link, 2]==0 ||
            data.points[a.link, 3]==0 || data.points[a.link, 4]==0)
            stop("Error in chromosome locations ..."); 
    }

    #   Ribbon colors
    #   ====================================================
    #
    ribbon.colors <- RCircos.Get.Link.Colors(ribbon.data, genomic.columns,
                            by.chromosome);

    #   Draw each ribbon (polygon)
    #   ============================
    #
    RCircos.Pos <- RCircos.Get.Plot.Positions();;
    base.positions <- RCircos.Pos*ribbon.start;

    for(a.ribbon in seq_len(nrow(ribbon.data)))
    {
        start.one <- data.points[a.ribbon, 1];
        end.one <- data.points[a.ribbon, 2];

        if(twist==FALSE) {
            start.two <- data.points[a.ribbon, 3];
            end.two <- data.points[a.ribbon, 4];
        } else {
            start.two <- data.points[a.ribbon, 4];
            end.two <- data.points[a.ribbon, 3];
        }

        P0 <- as.numeric(base.positions[end.one,]);
        P2 <- as.numeric(base.positions[start.two,]);
        line.one <- RCircos.Link.Line(P0, P2);

        P0 <- as.numeric(base.positions[end.two,]);
        P2 <- as.numeric(base.positions[start.one,]);
        line.two <- RCircos.Link.Line(P0, P2);

        polygon.x<- c(base.positions[start.one:end.one,1], line.one$pos.x,
                    base.positions[start.two:end.two,1], line.two$pos.x );
        polygon.y<- c(base.positions[start.one:end.one,2], line.one$pos.y,
                    base.positions[start.two:end.two,2], line.two$pos.y );
        polygon(polygon.x, polygon.y, border=NA, col=ribbon.colors[a.ribbon]);
    }
}




#   =========================================================================
#
#   10. RCircos.Clear.Track()
#
#   Erase one track or center area
#
#   Arguments:
#
#       track.num:      Non-negative number, the numbers of track to be erased 
#                       (e.g., 2, 3, 4, even 3.5)
#       side:           Character vector, relative position to chromosome 
#                       ideogram, must be either "in" or "out".
#       to.center:      Logical, TRUE for erase inner area including current 
#                       track and FALSE for clear current track only
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#
#   Returned value: None
#
#   Example:    RCircos.Clear.Track(track.num=2, side="in", to.center=FALSE);
#

RCircos.Clear.Track <- function(track.num=NULL, side=NULL, 
                to.center=FALSE, inside.pos=NULL, outside.pos=NULL)
{
    if(to.center!=FALSE && to.center!=TRUE)
    {  stop("to.center must be either TRUE or FALSE.\n"); }

    #   Adjust the far position relative to chromosome ideogram
    #   =======================================================
    #
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    boundary <- RCircos.Get.Plot.Boundary(track.num=track.num, 
        side=side, inside.pos=inside.pos, outside.pos=outside.pos, 
        erase.area=TRUE);
    if(boundary[1] <= RCircos.Par$chr.ideo.pos) {
        outerPos <- boundary[1];
        innerPos  <- boundary[2] - RCircos.Par$track.padding;
    } else {
        outerPos <- boundary[1] + RCircos.Par$track.padding;
        innerPos  <- boundary[2];
    }

    #   Clear all inner area including current track 
    #   ===================================================
    #
    if(to.center == TRUE)
    {
        polygon.x <- RCircos.Pos[,1]*outerPos;
        polygon.y <- RCircos.Pos[,2]*outerPos;

    #   Clear current track area only
    #   ==============================
    #
    } else {
        start <- 1;  end <- nrow(RCircos.Pos);
        polygon.x <- c(RCircos.Pos[start:end,1]*outerPos, 
                       RCircos.Pos[end:start,1]*innerPos);
        polygon.y <- c(RCircos.Pos[start:end,2]*outerPos, 
                       RCircos.Pos[end:start,2]*innerPos);
    }

    polygon(polygon.x, polygon.y, col="white", border="white");
}




#   =========================================================================
#
#    11.    RCircos.Track.Outline()
#
#   Draw outline of one plot track. 
# 
#    Arguments:
#
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       num.layers:     Non-negative integer, number of sub-tracks to plot, 
#                       0 for no sub-track line
#       chrom.list:     Character vector of chromosome name list or NULL
#       track.colors:   Color names for each chromosome in chrom.list or all.
#
#    Return value:    None
#
#    Example:    RCircos.Track.Outline(0.75, 0.65, 5, NULL)
#

RCircos.Track.Outline <- function(inside.pos=NULL, outside.pos=NULL, 
                    num.layers=1, chrom.list=NULL, track.colors=NULL)
{
    if(is.null(outside.pos) || is.null(inside.pos)) 
        stop("Missing outside.pos/inside.pos in RCircos.Track.Outline().\n")

    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Sub-track height. Note: Some times there may have 
    #   more or less subtracks than default maximum layers
    #   ===================================================
    #
    subtrack.height <- (outside.pos-inside.pos)/num.layers;
    chromosomes <- unique(as.character(RCircos.Cyto$Chromosome));

    #   In case one one or more but not all chromosome outlines
    #   need to be drawn
    #   ========================================================
    #
    if(!is.null(chrom.list)) { 
        if(sum(chrom.list %in% chromosomes)!= length(chrom.list))
        {  
            stop(paste("One or more chromosome is not",
                "in chromosome ideogram data.\n"));  
        }
        chromosomes <- chrom.list;
    } 

    if(is.null(track.colors)) {
        track.colors <- rep(RCircos.Par$track.background, length(chromosomes))
    } else {
        if(length(track.colors) != length(chromosomes))
            track.colors <- rep(track.colors, length(chromosomes));
    }

    for(aChr in seq_len(length(chromosomes)))
    {
        chr.rows <- which(RCircos.Cyto$Chromosome == chromosomes[aChr]);
        the.chr  <- RCircos.Cyto[chr.rows,];
 
        plot.start <- min(RCircos.Cyto$StartPoint[chr.rows]);
        plot.end   <- max(RCircos.Cyto$EndPoint[chr.rows]);

        polygon.x<- c(RCircos.Pos[plot.start:plot.end, 1]*outside.pos, 
                        RCircos.Pos[plot.end:plot.start, 1]*inside.pos);
        polygon.y<- c(RCircos.Pos[plot.start:plot.end,2]*outside.pos, 
                        RCircos.Pos[plot.end:plot.start,2]*inside.pos);
        polygon(polygon.x, polygon.y, col=track.colors[aChr]);

        if(num.layers>1) {
            for(a.line in seq_len(num.layers-1))
            {
                height <- outside.pos - a.line * subtrack.height;
                lines(RCircos.Pos[plot.start:plot.end,1] * height, 
                        RCircos.Pos[plot.start:plot.end,2] * height,
                        col=RCircos.Par$grid.line.color);
            }
        }
    }
}




#   =========================================================================
#
#   12. RCircos.Vertical.Line.Plot()
#
#   Plot vertical lines on a track without track outlines and sub-track lines
#
#   Arguments:
#
#       line.data:      A data frame with four columns for start and/or end 
#                       position, height in the data track, and plot colors
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram.
#       side:           Character vector, must be either "in" or "out"
#       line.width:     Non-negative integer, width of lines
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be 3.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#   Return value:   None
#   Example:        data(RCircos.Line.Data);
#                   line.data <- RCircos.Get.Track.Plot.Data(line.data,
#                               plot.type="vLine");
#                   RCircos.Vertical.Line.Plot(line.data, 1, "in")
#                   RCircos.Vertical.Line.Plot(line.data, 1, inside.pos=2, 
#                               outside.pos=2.2)
#

RCircos.Vertical.Line.Plot <- function(line.data=NULL, track.num=NULL, 
        side=c("in", "out"), line.width=1, inside.pos=NULL, outside.pos=NULL,
        genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(line.data)) 
        stop("Genomic data missing in RCircos.Vertical.Line.Plot.\n");
    if(is.null(genomic.columns) || genomic.columns < 2) 
        stop("Missing number of columns for genomic position.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];

    if(length(line.width) == 1)
        line.width <- rep(line.width, nrow(line.data));

    line.data <- RCircos.Get.Single.Point.Positions(line.data, 
                        genomic.columns);
    point.index <- as.numeric(line.data$Location);
    
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    line.colors <- RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color);

    for(a.point in seq_len((nrow(line.data))))
    {
        lines(  c(RCircos.Pos[point.index[a.point] ,1]*innerPos,
                    RCircos.Pos[point.index[a.point], 1]*outerPos),
                c(  RCircos.Pos[point.index[a.point], 2]*innerPos,
                    RCircos.Pos[point.index[a.point], 2]*outerPos),
                col=line.colors[a.point], lwd=line.width[a.point]);
    }
}




#   =========================================================================
#
#   13. RCircos.Point.Plot()
#
#   Plot points on a track without track outlines and sub-track lines. This 
#   function only use one column of data and does not consider the data 
#   range of whole datasets even there are more than one data columns.
#
#   Argument:
#
#       point.data:     A data frame containing genomic position and
#                       plot data values such as NGS read depth 
#       data.col:       Non-negative integer, the column number of data
#                       to be plotted 
#       track.num:      Non-negative integer, the number of track to plot
#       side:           Character vector, either "in" or "out" for inside
#                       outside of chromosome ideogram

#       point.type:     Non-negative integer for pch, default is 16
#       with.height:    Logical, if the point heights vary based on data values
#       with.size:      Logical, if the point sizes vary based on data values
#       point.scale:    Non-negative numeric, more scale for point size,
#                       must be greater than or equal to 1.
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area

#
#       Return value:   None
#       Example:        RCircos.Point.Plot(point.data=readDepth, data.col=4, 
#                               inside.pos=1.5, outside.pos=2)
#
#       Last modified on: July 10, 2015
#

RCircos.Point.Plot <- function(point.data=NULL, data.col=4, track.num=NULL, 
        side=c("in", "out"), min.value=NULL, max.value=NULL, 
        point.type=19, with.height=TRUE, with.size=FALSE, point.scale=1,
        inside.pos=NULL, outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(point.data) || !is.data.frame(point.data))
        stop("Missing or incorrect genomic data in RCircos.Point.Plot().\n");
    if( genomic.columns < 2 || genomic.columns > 3) 
        stop("Incorrect number of columns for genomic position.\n");
    if(data.col <= genomic.columns || !is.numeric(data.col))  
        stop("Plot data column must be ", genomic.columns+1, " or greater.\n");
    if(point.scale < 1) stop("Point scale must be greater than 1.\n")

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];
    trackHeight <- outerPos - innerPos;

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   Point size, color, and height from inner position of the track
    #   ==============================================================
    #
    point.values <- as.numeric(point.data[, data.col]);
    if(is.null(min.value) || is.null(max.value))
    {   min.value <- min(point.values);
        max.value <- max(point.values);
    }
    if(with.size == TRUE) {
        point.size <- point.values/(max.value-min.value)*point.scale;
    } else (point.size <- rep(point.scale, nrow(point.data)))

    if(with.height == TRUE) {
        point.height <- RCircos.Get.Data.Point.Height(point.values, min.value, 
            max.value, plot.type="points", track.height=trackHeight);
        point.height <- point.height + innerPos;
    } else { point.height <- rep(mean(boundary), nrow(point.data)); }
    
    point.color <- RCircos.Get.Plot.Colors(point.data, RCircos.Par$line.color);

    #   Plot each point
    #   =================================================================
    #
    RCircos.Track.Outline(outside.pos, inside.pos, num.layers=5)
    point.location <- RCircos.Get.Single.Point.Positions(point.data, 
                        genomic.columns);
    for(a.point in seq_len(length(point.values)))
    {
        index <- point.location$Location[a.point];
        points(RCircos.Pos[index, 1]*point.height[a.point],
                RCircos.Pos[index, 2]*point.height[a.point],
                pch=point.type, col=point.color[a.point], 
                cex=point.size[a.point]);
    }
}




#   =========================================================================
#
#   14. RCircos.Parallele.Link.Plot()
#
#   Plot horizontal line inside a track to show data overlaps. Horizontal 
#   lines are link lines between two genomic positions on same chromosome 
#
#   Arguments:
#
#       line.data:      A data frame containing chromosome, start and end 
#                       position for each line.
#       track.num:      Non-negative integer, the number of track to plot
#       side:           Character vector, either "in" or "out"
#       line.width:     lwd in graphic package
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns:   Non-negative integer, total number of colums for
#                           genomic position
#       is.sorted:      Logical, if the data have been sorted
#
#   Return value:   None
#
#   Example:        RCircos.Parallel.Line.Plot(line.data, 5, "in")
#                   RCircos.Parallel.Line.Plot(line.data, inside.pos=2,
#                        outside.pos=2.5)
#

RCircos.Parallel.Line.Plot <- function(line.data=NULL, track.num=NULL, 
         side=c("in", "out"), line.width=1, inside.pos=NULL, outside.pos=NULL, 
         genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(line.data)) stop("Line data missing in RCircos.Link.Plot.\n");
    if(is.null(genomic.columns) || genomic.columns < 3) 
        stop("Total number of columns for genomic position must be 3.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];
    track.height <- outerPos - innerPos;

    line.data <- RCircos.Get.Paired.Points.Positions(line.data, 
                    genomic.columns, plot.type="tile");

    layers <- RCircos.Get.Plot.Layers(line.data, genomic.columns);
    line.height <- layers*(track.height/max(layers)) + innerPos;

    #   Get link line colors for each pair of locations
    #   ================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    line.colors <- RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color);

    #   Draw link lines for each pair of genomic positions
    #   ==========================================
    #
    RCircos.Track.Outline(outside.pos, inside.pos, max(layers));
    for(a.line in seq_len(nrow(line.data)))
    {
        lines(c(RCircos.Pos[line.data$LinkStart[a.line],1]*line.height[a.line], 
            RCircos.Pos[line.data$LinkEnd[a.line],1]*line.height[a.line]), 
            c(RCircos.Pos[line.data$LinkStart[a.line],2]*line.height[a.line], 
               RCircos.Pos[line.data$LinkEnd[a.line],2]*line.height[a.line]),      
            type="l", col=line.colors[a.line], lwd=line.width); 
    }
}




#   =========================================================================
#
#   15. RCircos.polygons.Plot()
#
#   Plot polygons with different height and different locations inside of a 
#   track. Polygon plot is an alternative bar plot that takes both positive 
#   and negative height values and genomic intervals with different lengths   
#   Plot data should have genomic positions(chromosome names, start and end
#   positions) as well as height values. Optional column for polygon colors 
#   may follow.
#
#   Arguments:
#
#       polygon.data:   A data frame with three columns for genomic positions,
#                       one or more column for polygon heights. An optional 
#                       column for polygon colour may be in the last column.
#       data.col:       Non-negative integer, column number in polygon data for  
#                       polygon height. 
#       track.num:      Non-negative integer, which track to plot
#       side:           Character vector, must be either "in" or "out"
#       border.col:     Vector of color names for border color, default null
#       polygon.col:    Color name for fill of polygon
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area

#       genomic.columns:   total number of columns for genomic position, 
#                           must be 2 or 3.
#       is.sorted:      Logic, is the data is sorted or not
#
#   Return value:       None
#
#   Example:    RCircos.Polygon.Plot(polygon.data, 2, "in")
#               RCircos.Polygon.Plot(polygon.data, inside.pos=5, 
#                                       outside.pos=5.2)
#

RCircos.Polygon.Plot <- function(polygon.data=NULL, data.col=NULL,
        track.num=NULL, side=c("in", "out"), border.col=NULL, 
        polygon.col="pink", min.value=NULL, max.value=NULL,
        inside.pos=NULL, outside.pos=NULL, genomic.columns=3,  
        is.sorted=TRUE)
{
    if(is.null(polygon.data) || ncol(polygon.data) < 4) 
        stop("Polygon data with at least 4 columns must be defined.\n");
    if(is.null(genomic.columns) || genomic.columns != 3) 
        stop("Chromosome names, start and end positions are needed.\n");

    #   Get track position from track number or from user defined
    #   =============================================================
    #
    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];
    track.height <- outerPos - innerPos;

    data.heights <- as.numeric(polygon.data[,data.col])
    if(is.null(min.value) || is.null(max.value)) {
        min.value <- min(data.heights);
        max.value <- max(data.heights);
    }
    polygon.heights <- RCircos.Get.Polygon.Height(data.heights,
        min.value, max.value, inside.pos=innerPos, outside.pos=outerPos);
    bottom <- polygon.heights[,1];
    top    <- polygon.heights[,2];

    polygon.colors <- RCircos.Get.Plot.Colors(polygon.data, polygon.col)
    if(is.null(border.col)) 
        border.col <- rep(NA, nrow(polygon.data))
    else border.col <- rep(border.col, nrow(polygon.data))

    polygon.data <- RCircos.Get.Paired.Points.Positions(polygon.data, 
            genomic.columns, plot.type="polygon");
    the.start <- polygon.data$LinkStart;
    the.end <- polygon.data$LinkEnd;

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    if(min.value >= 0 || max.value <= 0) {
        layers <- RCircos.Par$sub.tracks;
    } else { layers <- 10; }
    RCircos.Track.Outline(outerPos, innerPos, layers);
    
    for(a.row in seq_len(nrow(polygon.data)))
    {
        a.start <- the.start[a.row];
        a.end   <-  the.end[a.row];
        polygon.x <- c(RCircos.Pos[a.start:a.end, 1]*top[a.row], 
                        RCircos.Pos[a.end:a.start, 1]*bottom[a.row]);
        polygon.y <- c(RCircos.Pos[a.start:a.end, 2]*top[a.row], 
                        RCircos.Pos[a.end:a.start, 2]*bottom[a.row]);

        polygon(polygon.x, polygon.y, col=polygon.colors[a.row], 
                        border=border.col);
    }
}




#   =========================================================================
#
#   16. RCircos.Area.Highlight()
#
#   Highlight a plot area with transparent color
#
#   Argument:   
#
#       highlight.area: Vector with chromosome name, start position, and
#                       end position to be highlighted
#       track.nums:     Vector of non-negative integer, which track(s) to
#                       be highlighted.
#       side:           Character vector, either "in" or "out"
#       color:          An RGB color definition, alpha value must be defined
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#
#   Return value:   None
#
#   Example:    highlight.area <- c("chr1", 100000, 200000);
#               RCircos.Area.Highlight(highlight.area, c(1:3), "in")
#               RCircos.Area.Highlight(highlight.area, inside.pos=3, 
#                                       outside.pos=3.5)
#

RCircos.Area.Highlight <- function(highlight.area=NULL, track.num=NULL, 
                side=c("in", "out"), hightlight.color=rgb(0.5, 0.5, 0, 0.5), 
                inside.pos=NULL, outside.pos=NULL)
{
    if(is.null(highlight.area))
        stop("Missing highlight area definition.\n");

    RCircos.Validate.Genomic.Info(highlight.area);
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    
    if(length(track.num) == 1) {
        boundary <- RCircos.Get.Plot.Boundary(track.num, side, 
                    inside.pos, outside.pos, FALSE);
    } else {
        track.one <- RCircos.Get.Plot.Boundary(track.num[1], side, 
                                inside.pos, outside.pos, FALSE);
        track.two <- RCircos.Get.Plot.Boundary(max(track.num), side, 
                                inside.pos, outside.pos, FALSE);
        if(side == "in") {
            boundary <- c(track.one[1], track.two[2]);
        } else {boundary <- c(track.one[2], track.two[1]);}
    }
    outerPos <- boundary[1];
    innerPos  <- boundary[2];

    start.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                as.numeric(highlight.area[2]));
    end.pos <- RCircos.Data.Point(as.character(highlight.area[1]), 
                                as.numeric(highlight.area[3]));

    polygon.x <- c(RCircos.Pos[start.pos:end.pos, 1]*outerPos, 
                        RCircos.Pos[end.pos:start.pos, 1]*innerPos);
    polygon.y <- c(RCircos.Pos[start.pos:end.pos, 2]*outerPos, 
                        RCircos.Pos[end.pos:start.pos, 2]*innerPos);

    polygon(polygon.x, polygon.y, col=hightlight.color);    
}




#   =========================================================================
#
#   17. RCircos.Customized.Shape.Plot()
#
#   Plot one customized shape on a track. The customized shape 
#   should be represented by coordinates of a polygon inside a circle 
#   with radius of 1 and default location of 12 o'clock. When plotting, 
#   the polygon center will be scaled and transformed for new size and 
#   location. For example, following code will plot an arrow:
#
#       plot(c(-2, 2), c(-2, 2))
#       polygonX <- c(0, -0.7, -0.2, -0.2, 0.2, 0.2, 0.7, 0)
#       polygonY <- c(-1, 0.7,  0.4,  1,   1,   0.4, 0.7, -1)   
#
#       nodeRadius <- 2;
#       polygon(polygonX*nodeRadius, polygonY*nodeRadius, col="red")
#       polygon(polygonX, polygonY, col="white")
#
#       nodeRadius <- 0.5
#       polygon(polygonX*nodeRadius, polygonY*nodeRadius, col="blue")
#
#   Arguments:
#
#       shape.data:     A two dimensional numeric matrix for coordinates 
#                       of a polygon
#       track.num:     Non-negative numeric vector for tracks to plot
#       side:           Character vector either "in" or "out"
#       location:       Vector with chromosome name, start position, and
#                       end position where the shap to be drawn
#       shape.color:    Character vector, color for the shape
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#
#   Return valuse:  None
#   Example:    RCircos.Customized.Shape.Plot()
#

RCircos.Customized.Shape.Plot <- function(shape.data=NULL,  
        track.num=NULL, side=c("in", "out"), location=NULL, 
        shape.color="red", inside.pos=NULL, outside.pos=NULL)
{
    if(is.null(shape.data) || is.null(location))
        stop("Missing argument in RCircos.Customized.Shape.Plot().\n");
    if(!is.numeric(shape.data) || !is.matrix(shape.data))
        stop("Shape data should be in a two dimension numeric matrix.\n");

    RCircos.Validate.Genomic.Info(location);
    
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    
    if(length(track.num) == 1) {
        boundary <- RCircos.Get.Plot.Boundary(track.num, side, 
                    inside.pos, outside.pos, TRUE);
    } else {
        track.one <- RCircos.Get.Plot.Boundary(track.num[1], side, 
                                inside.pos, outside.pos, TRUE);
        track.two <- RCircos.Get.Plot.Boundary(max(track.num), side, 
                                inside.pos, outside.pos, TRUE);
        if(side == "in") {
            boundary <- c(track.one[1], track.two[2]);
        } else {boundary <- c(track.one[2], track.two[1]);}
    }
    outerPos <- boundary[1];
    innerPos  <- boundary[2];
    center <- mean(c(outerPos, innerPos));
    scale.factor <- (outerPos - innerPos);

    #   Rotate shape data first
    #   ===========================================
    #
    point.index <- RCircos.Data.Point(as.character(location[1]), 
                        as.numeric(location[2]));
    angle <- 2*pi / nrow(RCircos.Pos) * point.index * -1;

    shapeX <- as.numeric(shape.data[,1]);
    shapeY <- as.numeric(shape.data[,2]);

    newX <- shapeX*cos(angle) -  shapeY*sin(angle);
    newY <- shapeX*sin(angle) +  shapeY*cos(angle);

    #   Scale and transform shape data 
    #   ===========================================
    #
    scaledX <- newX * scale.factor;
    scaledY <- newY * scale.factor;

    polygonX <- scaledX + RCircos.Pos[point.index, 1] * center;
    polygonY <- scaledY + RCircos.Pos[point.index, 2] * center;

    #   Plot the shape 
    #   ===========================================
    #
    polygon(polygonX, polygonY, col=shape.color);
}




#   =========================================================================
#
#   18. RCircos.Get.Start.End.Locations()
#
#   Generate plot start and end positions for special plot types such as
#   heatmap, histogram, or rectangle-like which needs start and end positions
#   along the track
#
#   Arguments:
#
#       plot.data:  Data frame returned from the function of 
#                   RCircos.Get.Single.Point.Positions() and the last
#                   column is the plot locations
#       plot.width: Non-negative integer in number of base pairs, e.g.,
#                   width or heatmap or histogram of each data points
#
#   Returned values:    Numeric matrix of two columns for start and end
#                       positions
#
#   Example:    RCircos.Par <- RCircos.Get.Plot.Parameters()
#               plot.locations <- RCircos.Get.Start.End.Locations(plot.data,
#                                   RCircos.Par$hist.width)
#

RCircos.Get.Start.End.Locations <- function(plot.data, plot.width)
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();

    dataChroms <- as.character(plot.data[,1]);
    chromosomes <- unique(dataChroms);
    cyto.chroms <- as.character(RCircos.Cyto$Chromosome);

    point.loc <- as.numeric(plot.data$Location)
    locations <- cbind(point.loc-plot.width, point.loc+plot.width)

    for(aChr in seq_len(length(chromosomes)))
    {
        cyto.rows <- which(cyto.chroms==chromosomes[aChr]);
        chr.start <- min(RCircos.Cyto$StartPoint[cyto.rows]);
        chr.end   <- max(RCircos.Cyto$EndPoint[cyto.rows]);

        data.rows <- which(dataChroms==chromosomes[aChr]);

        start.outliers <- which(locations[data.rows, 1] < chr.start)
        if(length(start.outliers)>0) 
            locations[data.rows[start.outliers], 1] <- chr.start;

        end.outliers <- which(locations[data.rows, 2] > chr.end)
        if(length(end.outliers)>0) 
            locations[data.rows[end.outliers], 2] <- chr.end;
    }

    return (locations);
}




#   =========================================================================
#
#   19. RCircos.Adjust.Scatter.Values()
#
#   Adjust scatter plot values to fit the plot track. The scatter with
#   lowest value will be plot on the centre of first subtrack and the one
#   with highest value will be on the centre of last subtrack. After this
#   adjustment, values could be added to inside position of plot track as
#   final plot positions.
#
#   Argument:
#
#       scatter.values: Numeric vector, the data to be plotted
#       max.value:      Numeric, the highest value to be plotted
#       min.value:      Numeric, the lowest value to be plotted
#       track.height:   Non-negative numeric, the height of plot area
#       subtrack:       Non-negative nummeric, the number of sub tracks
#
#   Returned value:     Adjusted scatter values
#
#   Example:    data(RCircos.Scatter.Data)
#               scatter.data <- as.numeric(RCircos.Scatter.Data[, 5])
#               plot.data <- RCircos.Adjust.Scatter.Values(scatter.data,
#                               track.heigh=0.2)
#
RCircos.Adjust.Scatter.Values <- function(scatter.values=NULL, min.value=NULL, 
             max.value=NULL, track.height=NULL, subtrack=NULL)
{
    if(is.null(scatter.values) || is.null(track.height) )
        stop("Missing argument for RCircos.Adjust.Scatter.Values().\n");
    if(is.null(max.value) || is.null(min.value))
    {
        min.value=min(scatter.values);
        max.value=max(scatter.values);
    }
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    if(is.null(subtrack) || subtrack < RCircos.Par$sub.tracks)
        subtrack <- RCircos.Par$sub.tracks;

    subtrack.height <- track.height/subtrack;
    plot.range <- track.height - subtrack.height;
    value.range <- max.value - min.value; 

    scatter.values <- (scatter.values - min.value)/value.range;
    scatter.values <- scatter.values*plot.range + subtrack.height/2;
    
    return (scatter.values);
}




#   =========================================================================
#
#   20.     RCircos.Area.plot()
#
#   Highligh part area of a data track. There will be three area types
#   to be plotted:
#
#   1.  Bottom is along inside of data track and top is a continue lines
#       connecting each data point with different height (mountain-type)
#
#   2.  Reversed mountain-type plot (Bottom and top in 1 is switched)
#
#   3.  Both bottom and top of the area are continue lines connecting each
#       each data point with different height (confidence region)
#
#		Arguments:
#
#       area.data:      A data frame with leading columns for genomic positions
#                       followed by data column(s). Data must be sorted first  
#                       by chromosome names then start positions and data values
#                       are better in range of 0 and 1.
#       data.col:       Non-negative integer, number of column in heatmap.data 
#                       for the data to be plotted
#       track.num:      Non-negative integer, number of the track from  
#                       chromosome ideogram
#       side:           Character vector, must be either "in" or "out"
#       plot.type:      Character vector, either "mountain", "curtain", or 
#                       "region"
#       min.value:      Numeric, minimum value of line data
#       max.value:      Numeric, maximum value of line data
#       area.color:     Character vector of R color name
#       border.col:     Character vector of R color name
#       inside.pos:     Non-negative number, inside position of the track 
#                       relative to the centre of plot area
#       outside.pos:    Non-negative number, outside position of the track  
#                       relative to the centre of plot area
#       genomic.columns: Non-negative integer, total number of columns for 
#                       genomic position in each row. Must be either 3 or 2.   
#       is.sorted:      Logic, whether the data is sorted by chromosome names
#                       and start position
#
#    Return value:    None
#
#    Example:    RCircos.Line.Plot(line.data, 4, 3, "in")
#
#   Last updated on:
#
#

RCircos.Area.Plot <- function(area.data=NULL, data.col=c(4,5), track.num=NULL, 
        side=c("in", "out"), plot.type=c("mountain", "curtain", "band"),
        min.value=NULL, max.value=NULL, area.color="gray", border.col="black",
        inside.pos=NULL, outside.pos=NULL, genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(area.data) || !is.data.frame(area.data))
        stop("Missing or incorrect genomic data in RCircos.Point.Plot().\n");
    if( genomic.columns < 2 || genomic.columns > 3) 
        stop("Incorrect number of columns for genomic position.\n");
    if(data.col <= genomic.columns || !is.numeric(data.col))  
        stop("Plot data column must be ", genomic.columns+1, " or greater.\n");
    if(plot.type == "band" && length(data.col) != 2)
            stop("Two columns of data required to plot band.\n");
    if(is.sorted == FALSE) area.data <- RCircos.Sort.Genomic.Data(area.data);
    
    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];
    trackHeight <- outerPos - innerPos;

    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    #   height of area top and area bottom inside a data track
    #   ==============================================================
    #
    area.values <- as.matrix(area.data[, data.col]);
    if(is.null(min.value) || is.null(max.value))
    {   min.value <- min(area.values);
        max.value <- max(area.values);
    }
    if(plot.type != "band")
    {
        area.values <- as.numeric(area.data[, data.col]);
        point.height <- RCircos.Get.Data.Point.Height(area.values, min.value, 
            max.value, plot.type="points", track.height=trackHeight);
        if(plot.type == "mountain")
        {
            area.bot <- rep(innerPos, length(point.height));
            area.top <- innerPos + point.height;
        } else {
            area.top <- rep(outerPos, length(point.height));
            area.bot <- outerPos - point.height;
        }
    } else {
        differ <- sum(area.values[,1] > area.values[, 2])
        if(differ == 0) {
            bottom <- 1; top <- 2; 
        } else if(differ == nrow(area.values)) {
            top <- 1; bottom <- 2;
        } else { stop("Band top cross area bottom.\n")}

        band.top <- as.numeric(area.values[, top]);
        area.top <- RCircos.Get.Data.Point.Height(band.top, min.value, 
            max.value, plot.type="points", track.height=trackHeight);
        area.top <- area.top + innerPos;
        
        band.bot <- as.numeric(area.values[, bottom]);
        area.bot <- RCircos.Get.Data.Point.Height(band.bot, min.value, 
            max.value, plot.type="points", track.height=trackHeight);
        area.bot <- area.bot + innerPos;
    }

    #   Location of each data point
    #   ==========================================================
    point.location <- RCircos.Get.Single.Point.Positions(area.data, 
                        genomic.columns);
    area.color <- RCircos.Get.Plot.Colors(area.data, area.color);
    
    #   Plot area by chromosome
    #   =================================================================
    #
    RCircos.Track.Outline(outside.pos, inside.pos, num.layers=5)
    chromosomes <- unique(area.data[,1])
    for(a.chr in seq_along(chromosomes))
    {
        data.index <- which(area.data[,1] == chromosomes[a.chr]);
        point.index <- as.numeric(point.location$Location[data.index]);

        polygon.x <- c(RCircos.Pos[point.index, 1]*area.top[data.index], 
            RCircos.Pos[rev(point.index), 1]*area.bot[rev(data.index)]);
        polygon.y <- c(RCircos.Pos[point.index, 2]*area.top[data.index], 
            RCircos.Pos[rev(point.index), 2]*area.bot[rev(data.index)]);

        polygon(polygon.x, polygon.y, col=area.color[data.index[1]], 
                    border=border.col);
    }
}




#   =========================================================================
#
#'		21.		RCircos.Customized.Connection.Plot()
#'
#'		Plot connection lines (simple lines plot) between two set 
#'		of data points. One example of usage is to label genes at
#'		modified plot position and connect the gene label to its
#'		genomic position.
#'
#'		Prerequisite:		
#'		RCircos core components and graphic device must be 
#'		initialized.
#'
#"
#'		Arguments:
#'
#'		gene.data:
#'		Data frame with chromosome name and actual genomic positions 
#'		for the genes to be labeled.
#'
#'		label.data: 
#'		Data frame with chromosome name and genomic positions that
#'		was used to label gene names
#
#'		gene.pos, label.pos:
#'		float numeric, scale factors relative to the center of plot   
#'		area (0). These two locations represent the start and end 
#'		location of lines.
#'
#'
#'		Returned value:	None.
#'

RCircos.Customized.Connection.Plot <- function(gene.data, label.data, 
			gene.pos=NULL, label.pos=NULL)
{
	gene_loc  <- RCircos.Get.Single.Point.Positions(gene.data, 3);
	label_loc <- RCircos.Get.Single.Point.Positions(label.data,3);
	
	Positions <- RCircos.Get.Plot.Positions(); 
	start_pos <- Positions[gene_loc$Location, 1:2]*gene.pos;
	  end_pos <- Positions[label_loc$Location, 1:2]*label.pos;

	for(a_loc in 1:nrow(gene.data)) {
		lines(c(start_pos[a_loc, 1], end_pos[a_loc,1]),
				c(start_pos[a_loc, 2], end_pos[a_loc,2]));
	}
}



#   End of RCircosPlotDataTracks.R
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
