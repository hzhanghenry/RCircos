#
#   Functions to initialize and reset RCircos plot positions and calculate
#   relevant positions for plot data
#
#   1.  RCircos.Set.Base.Plot.Positions()
#   2.  RCircos.Data.Point()
#   3.  RCircos.Get.Track.Positions()
#   4.  RCircos.Validate.Track.Positions()
#   5.  RCircos.Get.Gene.Label.Locations()
#   6.  RCircos.Link.Line()
#   7.  RCircos.Reset.Plot.Positions()
#   8.  RCircos.Get.Plot.Boundary()
#   9.  RCircos.Get.Gene.Name.Plot.Parameters()
#   10. RCircos.Get.Polygon.Height
#
#   Last debug done on September 14, 2016
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>




#   ========================================================================
# 
#   1.  RCircos.Set.Base.Plot.Positions()
#
#   Calculate the x and y coordinates for a circular line. These values will 
#   be used as base values to calculate the coordinates of plot tracks, data 
#   points, and chromosome band. Rotation degrees are also attached for text
#   labelling at each point.
#
#   Argument:
#
#       total.points:   Non-negative integer, how many points will be used
#                       to draw the circle. Default is null, means calculating
#                       it from chromosome ideogram data. This argument is
#                       reserved for advanced user only
# 
#   Return value:   None
#

RCircos.Set.Base.Plot.Positions<-function(total.points=NULL)
{
	#   If total points are not provided, calculating it from  
	#   chromosome ideogram data and add one padding length to 
	#	the last chromosome. 
	#   =========================================================
	#
    if(is.null(total.points)) {
        RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
        RCircos.Par <- RCircos.Get.Plot.Parameters();    
        total.points <- RCircos.Cyto$EndPoint[nrow(RCircos.Cyto)] +  
                        RCircos.Par$chrom.paddings
    }

    #   x and y coordinates for a circular line with radius of 1,
    #   circumference of 2*PI, and interval of 2PI/total.points
    #   ========================================================
    #
    interval <- 2*pi/total.points;
    baseVal <- seq(0, 2*pi, interval);

    coordinateX <- sin(baseVal);
    coordinateY <- cos(baseVal);

    #   Degrees for text rotating at each point
    #   ========================================
    #
    degree <- rep(0, length(baseVal));
    mid <- round((length(baseVal)-1)/2, digits=0) + 1;

    numOfPoints <- length(baseVal);
    degree[1:mid] <- 90 - (baseVal[1:mid]*180/pi);
    degree[(mid+1):numOfPoints] <- 270 - (baseVal[(mid+1):numOfPoints]*180/pi);

    #   Put the plot position data in RCircos environment
    #   =================================================
    #
    plotPostions <- data.frame(coor.x=coordinateX, coor.y=coordinateY, degree);

    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Base.Position"]] <- plotPostions;
}




#   ==========================================================================
# 
#   2.  RCircos.Data.Point()
#
#   Calculate plot coordinates for a chromosome position associated to a data 
#   point (e.g., the chromosome name and start position of a gene).
#
#   Arguments:
#
#       chromosome: Character vector, a chromosome name, e.g., "chr1"
#       start:      Start position of the gene on chromosomes
#
#   Return value:   An integer representing the index of x and y coordinates 
#                   of RCircos base plot positions
#   Example:        data.point <-  RCircos.Data.Point("chr1", 100000)
#

RCircos.Data.Point<-function(chromosome=NULL, start.posotion=NULL)
{
    if(is.null(chromosome) || is.null(start.posotion))
        stop("Missing argument for RCircos.Data.Point().\n")

    thePoint <- 0;
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Par  <- RCircos.Get.Plot.Parameters();

    #   Which band the start position is in. The start and end position
    #   is  are not overlap.
    #   ============================================================
    #
    chromRows <- grep(paste("^", chromosome, "$", sep=""), 
                        RCircos.Cyto$Chromosome);
    theRow <- which(RCircos.Cyto$ChromStart[chromRows] <= start.posotion & 
                    RCircos.Cyto$ChromEnd[chromRows] >= start.posotion)[1];
    
    #   How far from the chromosome start the point is (by units)
    #   ==============================================================
    #
    theBases <- start.posotion - RCircos.Cyto$ChromStart[chromRows[theRow]];
    theUnits  <- theBases/RCircos.Par$base.per.unit;
    
    thePoint <- RCircos.Cyto$StartPoint[chromRows[theRow]] + theUnits;

    #   Return the index. Arguments are validated outside of this 
    #   function so that there is no need to catch exception.
    #   =========================================================
    #
    return (round(thePoint, digits=0));
}




#   ========================================================================
#
#   3.  RCircos.Get.Track.Positions()
#
#   Calculate inner and outer positions for a data track
#
#   Arguments:
#
#       side:       Character vector, must be either "in" or "out".
#       track.num:  Non-negative integer, number of the track from chromosome 
#                   ideogram.
#
#   Return value:   The outer and inner positions of a track
#   Example:        track.loc <- RCircos.Track.Positions("in", 2)
#

RCircos.Get.Track.Positions<-function(side=NULL, track.num=NULL)
{
    #   argument validation
    #   ===========================================================
    #
    if(is.null(side) || is.null(track.num))
        stop("Missing function argument in RCircos.Track.Positions().\n")
    if(track.num < 1) stop("track number cannot be smaller than 1.\n"); 

    side <- tolower(side);
    if(side != "in" && side != "out") stop("side must be in or out.\n");

    RCircos.Par <- RCircos.Get.Plot.Parameters();
    oneTrack <- RCircos.Par$track.height + RCircos.Par$track.padding;

    #   Positions based on side
    #   ==============================================================
    #
    side <- tolower(side);
    if(side == "in") 
    {
        outPos <- RCircos.Par$track.in.start - (track.num-1)*oneTrack;
        inPos  <- outPos - RCircos.Par$track.height;

        if(outPos <= 0 || inPos <= 0) 
        { stop("Incorrect track location.\n"); }

    } else if(side == "out") {
        inPos  <- RCircos.Par$track.out.start +(track.num-1)*oneTrack;
        outPos <- inPos + RCircos.Par$track.height;
    
        RCircos.Par <- RCircos.Get.Plot.Parameters();
        if(outPos > RCircos.Par$plot.radius) 
            stop("Track position is out of plot area.\n");
    } else {  
        stop("Incorrect track location. It must be \"in\" or \"out\".\n"); 
    }

    #   The position needs to be held for the RCircos session
    #   =====================================================
    #
    return (locations=c(out.pos=outPos, in.pos=inPos));
}




#   =========================================================================
#
#   4.  RCircos.Validate.Track.Positions()
#
#   Calculate customized plot locations
#
#   Arguments:
#   
#       inside.pos:     Non-negative number, scale factor based on chromosome
#                       ideogram position for the close position relative to 
#                       the centre of plot area
#       outside.pos:    Non-negative number, scale factor based on chromosome
#                       ideogram position for the far position relative to 
#                       the centre of plot area
#
#   Return value:   A numeric vector with length of 2 for the real plot 
#                   position of a customized track
#
#   Example:    locations <- RCircos.Get.Customized.Positions(0.95, 0.90)
#

RCircos.Validate.Track.Positions <- function(inside.pos=0, outside.pos=0, 
                erase.area=FALSE)
{
    #   Plot positions must be positive numbers, inside.pos should
    #   be always smaller than outside.pos
    #   ==========================================================
    #
    if(inside.pos <= 0 || outside.pos <= 0)
        stop("Position <= 0 in RCircos.Validate.Track.Positions().\n");
    if(inside.pos > outside.pos) 
        stop("Outside position must be greater than inside position.");

    #   Inside and outside position cannot overlap with chromosome 
    #   ideogram and chromosome names unless to erase the area
    #   ==========================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    ideoPos <- RCircos.Par$chr.ideo.pos;
    namePos <- RCircos.Par$chr.name.pos;

    if(erase.area != TRUE) {
        if(inside.pos > RCircos.Par$track.in.start && 
            inside.pos < RCircos.Par$track.out.start)
            stop("Plot position overlap with chromosome ideogram.\n")
        if(outside.pos > RCircos.Par$track.in.start && 
            outside.pos < RCircos.Par$track.out.start)
            stop("Plot position overlap with chromosome ideogram.\n");
        if(inside.pos <= RCircos.Par$chr.ideo.pos && 
            outside.pos >= RCircos.Par$chr.ideo.pos)
            stop("Plot position overlap with chromosome ideogram.\n");
    }

    #   make adjustments to match the first track position inside 
    #   or outside chromosome ideogram if necessary
    #   =========================================================
    #
    customizedHeight <- outside.pos - inside.pos;

    if(inside.pos > RCircos.Par$chr.name.pos 
        && inside.pos < RCircos.Par$track.out.start ) {
            inPos <- RCircos.Par$track.out.start;
            outPos <- inPos + customizedHeight;
    } else if (outside.pos < RCircos.Par$chr.ideo.pos 
        && outside.pos > RCircos.Par$track.in.start) {
            outPos <- RCircos.Par$track.in.start
            inPos <- outPos - customizedHeight ;
    } else { 
        inPos <- inside.pos;  
        outPos <- outside.pos; 
    }

    return (locations=c(out.pos=outPos, in.pos=inPos));
}




#   =========================================================================
#
#   5.  RCircos.Get.Gene.Label.Locations()
#
#   In case of too many gene along one chromosome or a genomic region, the 
#   gene names may become very crowded so that we need reset plot positions 
#   to make the image more readable. This function is for plot gene names
#   and connectors next to chromosome ideogram track. For more gene names,
#   use zoomed plot instead on outside of chromosome ideogram.
#
#   Arguments:
#
#   genomic.data:   A data frame with the four columns for chromosome 
#                   name, start position, end position, and gene name
#   genomic.columns: Non-negative integer, total number of columns for
#                   genomic positions, must be 2 or 3.
#   plot.pos:       Non-negative numeric, inside position of plot track
#   is.sorted:      Logic, if the data is sorted by chromosomes
#
#   Return value:   All or subset of input data frame with a new column 
#                   for plot positions.
#
#   Example:    labelData <- RCircos.Get.Gene.Label.Locations(genomic.data)
#

RCircos.Get.Gene.Label.Locations <- function(genomic.data=NULL,
                genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(genomic.data)) 
        stop("Missing genomic.data in RCircos.Get.Gene.Label.Locations().\n");
    if(genomic.columns <2 || genomic.columns > 3)
        stop("Incorrect number for columns of genomic positions.");

    if(is.sorted == FALSE)
        genomic.data <- RCircos.Sort.Genomic.Data(genomic.data, is.ideo=FALSE);
        
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Pos <- RCircos.Get.Plot.Positions();

 
    #   Attach a new column to plot data for label locations. 
    #   ======================================================
    #
    genomic.data <- RCircos.Get.Single.Point.Positions(genomic.data, 
                            genomic.columns);
    genomic.data["LabelPosition"] <- genomic.data$Location;
    dataChr <- as.character(genomic.data$Chromosome);

    geneNameParameters <- RCircos.Get.Gene.Name.Plot.Parameters();
    cytoChr <- as.character(geneNameParameters[, 1]);

    
    #   Reset label locations
    #   ======================================================
    #
    toMuchChrom <- FALSE;   
    labelData <- NULL;  
    labelWidth <- as.numeric(geneNameParameters[1, 5]);
    for(aChr in seq_len(length(cytoChr)))
    {
        index <- which(dataChr == cytoChr[aChr]);
        if(length(index) == 0) { next; }

        #   If there too many gene labels, remove extra 
        #   genes for best visualization
        #   ===================================================
        #
        if(length(index) > geneNameParameters[aChr, 2])
        {
            toMuchChrom <- TRUE;
            theChr <- genomic.data[index[1:geneNameParameters[aChr, 2]],];
            theChr <- theChr[order(theChr$Location),];

            for(aGene in seq_len(nrow(theChr)))
            {
                newLoc <- geneNameParameters[aChr, 3] + (aGene-1) * labelWidth;  
                theChr$LabelPosition[aGene] <- newLoc;
            }
        } else {
            #   modify label locations if necessary
            #   ===========================================
            #
            theChr <- genomic.data[index,];
            theChr <- theChr[order(theChr$Location),];
            allGene <- nrow(theChr);
            lastPos <- theChr$LabelPosition[1];
            
            endLoc <- as.numeric(geneNameParameters$endLoc);
            for(aGene in seq_len(allGene))
            {
                currLoc <-  theChr$LabelPosition[aGene];
                lenNeeded <- currLoc + (allGene - aGene)*labelWidth;
                if(lenNeeded > endLoc[aChr]) 
                {  currLoc <- currLoc - (lenNeeded - endLoc[aChr]); }

                if(aGene > 1 && (currLoc-lastPos) < labelWidth)
                { currLoc  <- lastPos + labelWidth;  }

                theChr$LabelPosition[aGene] <- currLoc;
                lastPos <- currLoc;
            }
        }
        labelData <- rbind(labelData, theChr); 
    }
    colnames(labelData) <- colnames(genomic.data);

    if(toMuchChrom == TRUE) {
        message("Not all labels will be plotted.");
        message("Type RCircos.Get.Gene.Name.Plot.Parameters()");
        message("to see the number of labels for each chromosome.");
    }
    
    #   The position needs to be held for the RCircos session
    #   =====================================================
    #
    return (labelData);
}




#   =========================================================================
#
#   6.  RCircos.Link.Line()
#
#   Calculate x and y coordinates for a quandratic Bezier curve between two 
#   chromosome locations with the equation:  
#
#       B(t) = (1-t) ((1-t)P0 + tP1) + t((1-t)P1 + tP2)
#
#   where P0 is the start point, P2 is the end point, and P1 is the control 
#   point. Since we set P1 to (0,0), the equation become: 
#
#       B(t) =(1-t)^2*P0 + t^2*P2
#
#   Arguments:
#
#       line.start: The point where Bezier line starts
#       line.end:   The point where Bezier line ends
#
#   Return value:   A list containing x and y coordinates for a quandratic 
#                   Bezier curve
#
#   Example:    bzLine <- RCircos.Link.Line(c(0, -1), c(1, 1))

RCircos.Link.Line <- function(line.start=NULL, line.end=NULL)
{
    if(is.null(line.start) || is.null(line.end))
        stop("Missing function argument in RCircos.Link.Line().\n");

    #   Set up the points for Bezure curve
    #   ==================================
    #
    P0 <- line.start;
    P2 <- line.end;

    #   Calculate total number of points for the Bezier curve
    #   =====================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    numOfBCPoint <- RCircos.Par$Bezier.point;
    t <- seq(0, 1, 1/numOfBCPoint);

    #   Calculate the point values for Bezuer curve
    #   ===========================================
    #
    linkX <- (1-t)^2*P0[1] + t^2*P2[1];
    linkY <- (1-t)^2*P0[2] + t^2*P2[2];

    #   Return the coordinates
    #   ===========================================
    #
    return (list(pos.x=linkX, pos.y=linkY));
}




#   =========================================================================
#
#   7.  RCircos.Reset.Plot.Positions()
#
#   Arguments:
#
#       plot.positions: A Data frame of three columns for: x and y coordinates 
#                       for points which form a circular line and degrees at 
#                       each point the text should be rotated.
#
#   Returned values: None
#
#   Example:    Reserved for advanced usage.
#
#

RCircos.Reset.Plot.Positions<-function(plot.positions=NULL)
{
    if(is.null(plot.positions)) 
        stop("Missing function argument in RCircos.Reset.Plot.Positions().\n")

    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Base.Position"]] <- plot.positions;
}




#   =========================================================================
#
#   8.  RCircos.Get.Plot.Boundary()
#
#   Calculate plot track location based on track number, plot side, or user
#   defined in boundary and out boundary relative to plot area center
#
#   Argument:
#
#       track.num:      Non-negative integer, which track to plot
#       side:           Character vector, either "in" or "out"
#       inside.pos:     Non-negative numeric, a point that is close to 0
#       outside.pos:    Non-negative numeric, a point that is far from 0
#
#   Return value:       A non-negative numeric vector with length of 2 for
#                       outside and inside boundary of plot area relative
#                       to the center of whole plot area
#
#   Example:    locations <- RCircos.Get.Plot.Boundary(2, "in");
#               
#               RCircos.Par <- RCircos.Get.Plot.Parameters();
#               start.pos <- RCircos.Par$track.out.start;
#               locations <- RCircos.Get.Plot.Boundary(inside.pos=start.pos, 
#                               outside.pos=start.pos + 0.5);
#

RCircos.Get.Plot.Boundary <- function(track.num=NULL, side=NULL,
                inside.pos=NULL, outside.pos=NULL, erase.area=FALSE)
{
    if(is.null(inside.pos) == FALSE && is.null(outside.pos) == FALSE) {
        locations <- RCircos.Validate.Track.Positions(inside.pos, 
                        outside.pos, erase.area);
    } else {
        if(is.null(side) || is.null(track.num))
            stop("track.num and side must be defined.\n")
        locations <- RCircos.Get.Track.Positions(side, track.num);
    }
    
    return (locations);
}





#   =========================================================================
#
#   9.  RCircos.Get.Gene.Name.Plot.Parameters() 
#
#       Calculate maximum number of gene names to be plotted for each
#       chromosome based on plot RCircos plot parameters
#
#   Argument:   None
#   Returned values:    A matrix with total number of gene names, start
#                       and end locations for each chromosome
#   
#   Example:    genePlotPar <- RCircos.Get.Gene.Name.Plot.Parameters()
#
#

RCircos.Get.Gene.Name.Plot.Parameters <- function()
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    defaultUnits <- RCircos.Get.Default.Base.Per.Units();   # 30000 units
    defaultWidth <- RCircos.Get.Default.Char.Width();   # 500 units
    defaultSize  <- RCircos.Get.Default.Text.Size();    # cex = 0.4

    sizeFactor  <- RCircos.Par$text.size/defaultSize;
    unitFactor  <- RCircos.Par$base.per.unit/defaultUnits;
    widthFactor <- RCircos.Par$char.width/defaultWidth;
    
    labelWidth <- defaultWidth*sizeFactor*unitFactor*widthFactor;
    
    #   Get maximum number of labels for each chromosome. 
    #   ======================================================
    #
    cyto.chroms <- as.character(RCircos.Cyto$Chromosome);
    chromosomes <- unique(cyto.chroms);   
    bandLength <- as.numeric(RCircos.Cyto$ChromEnd) - 
                    as.numeric(RCircos.Cyto$ChromStart);
    bandUnits  <- round(bandLength/RCircos.Par$base.per.unit, digits=0);

    maxLabels <- rep(0, length(chromosomes));
    startLoc  <- rep(0, length(chromosomes));
    endLoc    <- rep(0, length(chromosomes));
    for(aChr in seq_along(chromosomes))
    {
        index <- which(cyto.chroms == chromosomes[aChr]);
        totalUnits <- sum(bandUnits[index]);
        maxLabels[aChr] <- floor(totalUnits/labelWidth);
        startLoc[aChr] <- min(RCircos.Cyto$StartPoint[index]);
        endLoc[aChr]   <- max(RCircos.Cyto$EndPoint[index]);
    }
    
    labelWidth <- rep(labelWidth, length(chromosomes))
    GeneNamePars <- data.frame(chromosomes, maxLabels, 
                        startLoc, endLoc, labelWidth);

    return (GeneNamePars);
}




#   =========================================================================
#
#   10.  RCircos.Get.Polygon.Height() 
#
#       Calculate height of polygons for different bottom and top positions
#
#   Argument:
#
#       data.heights:   Numeric vector for polygon heights
#       min.value:      Minimum value of polygon heights
#       max.value:      Maximum value of polygon heights
#       inside.pos:     Non-negative numeric, location of the track close
#                       to the center of plot area
#       outside.pos:    Non-negative numeric, location of the track far 
#                       from the center of plot area
#
#   Returned values:    A matrix with two columns for bottom and top height
#   
#   Example:    heights <- RCircos.Get.Polygon.Height(data.heights, -1, 1, 
#                               1.0, 1.5)
#
#
RCircos.Get.Polygon.Height <- function(data.heights, min.value=NULL, 
        max.value=NULL, inside.pos=NULL, outside.pos=NULL)
{
    #
    if(is.null(data.heights) || is.null(inside.pos) || is.null(outside.pos))
        stop("Missing function arguments in RCircos.Get.Polygon.Height()!")
    
    if(is.null(inside.pos) || is.null(outside.pos))
    {
        min.value <- min(data.heights);
        max.value <- max(data.heights)
    }
    track.height <- outside.pos - inside.pos;
    
    #   Polygon heights are all positive, draw polygon 
    #   from bottom of the data track
    if( min.value >= 0)
    {
        polygon.heights <- RCircos.Get.Data.Point.Height(data.heights, 
                min.value, max.value, plot.type="points", track.height);
        bottom <- rep(inside.pos, length(polygon.heights));
        top <- inside.pos + polygon.heights;

    #   Polygon heights are all negative, draw polygon 
    #   from top of the data track
    } else if( max.value <= 0) {
        polygon.heights <- RCircos.Get.Data.Point.Height(data.heights*-1, 
                min.value=0, max.value, plot.type="points", track.height);
        top <- rep(outside.pos, length(polygon.heights));
        bottom <- top - polygon.heights;

    #   Polygon heights have both negative and positive values, 
    #   draw polygon from the middle of the data track 
    } else {
        min.value <- 0; max.values <- max(abs(data.heights));

        positive.heit <- data.heights;
        negatives <- which(positive.heit < 0)
        positive.heit[negatives] <- 0;

        negative.heit <- data.heights;
        positives <- which(negative.heit > 0)
        negative.heit[positives] <- 0;

        positive.heights <- RCircos.Get.Data.Point.Height(positive.heit, 
                min.value, max.value, plot.type="points", track.height/2);
        negative.heights <- RCircos.Get.Data.Point.Height(negative.heit*-1,
                min.value, max.value, plot.type="points", track.height/2);
        negative.heights <- negative.heights*-1;

        polygon.heights <- positive.heights + negative.heights;
        bottom <- rep(inside.pos + track.height/2, length(polygon.heights));
        top <- bottom + polygon.heights;
    }

    return(cbind(bottom, top));
}


#   End of RCircosPlotPositions.R
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>



