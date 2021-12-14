#
#   Functions to initialize, reset, and plot RCircos chromosome ideogram
#
#   1.  RCircos.Validate.Cyto.Info()
#   2.  RCircos.Set.Cytoband.Data()
#   3.  RCircos.Chromosome.Ideogram.Plot()
#   4.  RCircos.Draw.Chromosome.Ideogram()
#   5.  RCircos.Highligh.Chromosome.Ideogram()
#   6.  RCircos.Label.Chromosome.Names()
#   7.  RCircos.Ideogram.Tick.Plot()
#   8.  RCircos.Reset.Plot.Ideogram()
#   9.  RCircos.ZoomIn.Chromosome()
#   10. RCircos.ZoomOut.Chromosome()
#   11. RCircos.Get.Chromosome.Order()
#   12. RCircos.Pseudo.Ideogram.From.Labels()
#   13. RCircos.Pseudo.Ideogram.From.Table()
#   14. RCircos.Validate.Genomic.Info()
#
#   Last debugged on December 7, 2016
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>




#   =========================================================================
#
#   1.  RCircos.Validate.Cyto.Info()
#
#       Validate chromosome ideogram information for correct chromosome order 
#       including of correct order of chromosome names, chromosome start, and 
#       end positions. The prefix of "chr" in chromosome names is no longer 
#       requiring since this version.
#
#   Argument:   
#
#       cyto.info   A data frame which containing the chromosome ideogram 
#                   data, e.g., an object returned by calls to the function
#                   of read.table() which read a file containing full 
#                   information of cytoBandIdeo table from UCSC genome 
#                   browser. This is not A RCircos component.
#
#       chr.exclude: Character vector, name(s) of chromosomes to be excluded
#                   from RCircos plot.
#       is.sorted:  Logic, whether the ideogram has been sorted.
#
#   Return value:   data frame with validated chromosome ideogram data  
#   Example:        ideo <- RCircos.Validate.Cyto.Info(cyto.info, "chrY")
#

RCircos.Validate.Cyto.Info <- function(cyto.info=NULL, chr.exclude=NULL,
                is.sorted=TRUE)
{
    if(is.null(cyto.info) || is.null(is.sorted)) 
            stop("Missing argument cyto.info.\n")
    if(is.data.frame(cyto.info) == FALSE)
        stop("Ideogram data must be in data frame.\n");

    #   Set standard headers for cyto.info
    #   =============================================================
    #   
    if(ncol(cyto.info) < 5) { 
        stop( paste("Cytoband data must have columns for Chromosome,",
                "ChromStart, ChromEnd, Band, Stain\n",
                "Current cyto.info columns:", ncol(cyto.info)));  
    }
    colnames(cyto.info)[1:5] <- c("Chromosome", "ChromStart", 
                "ChromEnd", "Band", "Stain");

    #   Remove rows for chromosomes to be excluded if any
    #   =================================================
    #
    if(length(chr.exclude) > 0) 
    {
        toExclude <- which(chr.exclude %in% cyto.info$Chromosome)
        if(length(toExclude) != length(chr.exclude) )
            stop("Not all chr.exclude found in ideogram data.\n");

        ignore <- which(cyto.info$Chromosome %in% chr.exclude);
        cyto.info <- cyto.info[-ignore, ]; 
    }

    #   Any ChromEnd must be greater than its ChromStart
    #   ================================================
    #
    endGTstart <- cyto.info$ChromEnd>cyto.info$ChromStart;
    if(sum(endGTstart) < length(endGTstart)) 
    {
        print(cyto.info[which(endGTstart == FALSE), ]);
        stop("ChromEnd position less than ChromStart position.\n"); 
    }

    #   Sort the ideogram table if need.
    #   ===============================
    #
    if(is.sorted != TRUE) 
        cyto.info <- RCircos.Sort.Genomic.Data(cyto.info, is.ideo=TRUE);

    #   Check start and end positions for each chromosome band. The first
    #   start position for each chromosome must be 0 or 1 and all other
    #   start positions must be greater than its previous end position.
    #   ==================================================================
    #
    chromosomes <- unique(cyto.info$Chromosome);
    for(aChr in seq_len(length(chromosomes)))
    {
        theRow <- which(cyto.info$Chromosome==chromosomes[aChr]);
        theInfo <- cyto.info[theRow,];
        theStart <- as.numeric(theInfo$ChromStart);
        theEnd   <- as.numeric(theInfo$ChromEnd);

        if(theStart[1] > 1) stop("Cytoband start should be 0.");

        #   If a chromosome has more than one bands, start positions 
        #   of each band must be greater than the end position of its
        #   previous band. 
        #   =========================================================
        #
        if(length(theRow) > 1) 
        { 
            for(aRow in seq_len(nrow(theInfo))[-1]) {
                if(theStart[aRow] < theEnd[(aRow-1)]) 
                { 
                    stop(paste("Cytoband start position cannot be ",
                            "less than previous end position.")); 
                }
            }
        }
    }

    #   Return the validated cyto.info 
    #   ===============================
    #
    return (cyto.info);
}




#   ==========================================================================
#
#   2.  RCircos.Set.Cytoband.Data()
#
#       Set chromosome cytoband data for ideogram plot. The cytoband data will
#       be put into RCircos Environment.
#
#   Arguments:
#
#       cyto.band.info: A data frame which containing the chromosome ideogram 
#                       data returned from RCircos.Validate.Cyto.Info() or
#                       a sorted ideogram table with correct chromosome names,
#                       start and end positions, band names, and staing status.
#
#   Returned values: None
#   Example:    data(UCSC.HG19.Human.CytoBandIdeogram);
#               RCircos.Set.Cytoband.data(UCSC.HG19.Human.CytoBandIdeogram);

RCircos.Set.Cytoband.Data <- function(cyto.band.info=NULL)
{
    if(is.null(cyto.band.info)) stop("Missing ideogram data.\n");
    if(is.data.frame(cyto.band.info) == FALSE)
        stop("Ideogram data must be in data frame or matrix.\n");
        
    #   Reset colors for chromosome bands. Use yellow color for unknown 
    #   ===============================================================
    #
    stain2color <- as.character(cyto.band.info$Stain);
    bandColor <- rep(colors()[652], length(stain2color));

    stains <- c("gneg", "acen", "stalk", "gvar", "gpos", "gpos100", 
        "gpos75", "gpos66", "gpos50", "gpos33", "gpos25");
    colorIndex <- c(1, 552, 615, 418, 24, 24, 193, 203, 213, 223, 233);

    for(aStain in seq_len(length(stains)))
    {
        bands <- which(stain2color == stains[aStain]);
        if(length(bands) > 0) 
        { bandColor[bands] <- colors()[colorIndex[aStain]]; }
    }
    cyto.band.info["BandColor"] <- bandColor;

    #   Assign colors to chromosome highlight. There are total 50
    #   colors and the last 26 colors are reserved for future.
    #   =========================================================
    #
    chromColor <- c(552, 574, 645, 498, 450, 81, 26, 584, 524, 472,
            32, 57, 615, 635, 547, 254, 100, 72, 630, 589,
            8, 95, 568, 52);
    chrom2color <- as.character(cyto.band.info$Chromosome);
    chromosomes <- unique(chrom2color);

    #   In case of multiple ideogram plot, recycle the colors
    #   ======================================================
    #
    numOfChrom <- length(chromosomes);
    numOfColor <- length(chromColor);
    if(numOfChrom > numOfColor)
    {
        recycleTime <- floor(numOfChrom/numOfColor) + 1;
        moreColors  <- rep(chromColor, recycleTime);
        chromColor <- moreColors[1:numOfChrom]
    }

    for(aChr in seq_len(length(chromosomes)))
    {
        rows <- which(chrom2color==chromosomes[aChr]);
        if(length(rows) > 0) 
        { chrom2color[rows] <- colors()[chromColor[aChr]]; }
    }
    cyto.band.info["ChrColor"] <- chrom2color;

    #   Total base pairs, relative length and location of each band
    #   are replaced with start point and end point in version 1.2
    #   ===========================================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    chromosomeStart <- as.numeric(cyto.band.info$ChromStart);
    chromosomeEnd   <- as.numeric(cyto.band.info$ChromEnd);

    #   chromosome 1
    #   ==========================================================
    #
    chromosomes <- unique(cyto.band.info$Chromosome);
    chrRows <- which(cyto.band.info$Chromosome == chromosomes[1])
    
    startIndex <- chromosomeStart[chrRows]/RCircos.Par$base.per.unit;
    endIndex   <- chromosomeEnd[chrRows]/RCircos.Par$base.per.unit;

    startIndex[1] <- 1;
    lastEnd  <- endIndex[length(chrRows)] + RCircos.Par$chrom.paddings;

    #   Other chromosomes
    #   ==========================================================
    #
    for(aChr in seq_len(length(chromosomes))[-1])
    {
        chrRows <- which(cyto.band.info$Chromosome == chromosomes[aChr])

        theStart <- chromosomeStart[chrRows]/RCircos.Par$base.per.unit;
        theEnd   <- chromosomeEnd[chrRows]/RCircos.Par$base.per.unit;

        theStart <- theStart +  lastEnd;
        theEnd  <- theEnd + lastEnd;

        startIndex <- c(startIndex, theStart);
        endIndex  <- c(endIndex, theEnd);

        lastEnd  <- theEnd[length(chrRows)] + RCircos.Par$chrom.paddings;
    }
    cyto.band.info["StartPoint"] <- round(startIndex, digits=0);
    cyto.band.info["EndPoint"] <- round(endIndex, digits=0);
    
    #   Put the cyto.band.info data in RCircos environment
    #   ==================================================
    #
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Cytoband"]] <- cyto.band.info;
}




#   =========================================================================
#
#   3.  RCircos.Chromosome.Ideogram.Plot()
#
#   Draw chromosome ideogram. Graphic device must be initialized first. The 
#   original function was split to three new ones in order to plot chromosome
#   highlight or names in different area.
#
#   Argument:
#
#       tick.interval:  Non-negative integer, length between two ticks
#                       in million base pairs. Set to 0 to ignore the 
#                       ticks.
#
#   Return value:   None
#
#   Example:    RCircos.Chromosome.Ideogram.Plot();
#

RCircos.Chromosome.Ideogram.Plot<-function(tick.interval=0)
{
    RCircos.Draw.Chromosome.Ideogram();
    RCircos.Highligh.Chromosome.Ideogram();

    if(tick.interval>0) {
        RCircos.Ideogram.Tick.Plot(tick.interval);
    }

    RCircos.Label.Chromosome.Names();
}




#   =========================================================================
#
#   4.  RCircos.Draw.Chromosome.Ideogram()
#
#   Draw chromosome ideogram only (outlines and dark bands). Graphic device 
#   must be initialized first
#
#   Argument:       None
#   Return value:   None
#
#   Example:        RCircos.Drawn.Chromosome.Ideogram()
#

RCircos.Draw.Chromosome.Ideogram <- function (ideo.pos=NULL, ideo.width=NULL)
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos  <- RCircos.Get.Plot.Positions()
    RCircos.Par  <- RCircos.Get.Plot.Parameters()
    
    if(is.null(ideo.pos)) ideo.pos <- RCircos.Par$chr.ideo.pos;
    if(is.null(ideo.width)) ideo.width <- RCircos.Par$chrom.width;
    
    #   Plot outlines for each chromosome
    #   =================================
    #
    outerPos <- ideo.pos + ideo.width;
    innerPos <- ideo.pos;

    chromosomes <- unique(RCircos.Cyto$Chromosome);
    RCircos.Track.Outline(outerPos, innerPos, num.layers=1, chromosomes,
            track.colors=rep("white", length(chromosomes)));

    #   Add chromosome bands (Giemsa stain positive only)
    #   ================================================
    #
    whiteBands <- which(RCircos.Cyto$BandColor == "white");
	darkBands <- RCircos.Cyto;	
	if(length(whiteBands)>0) darkBands <- darkBands[-whiteBands, ];

    for(aBand in seq_len(nrow(darkBands)))
    {
        aColor <- darkBands$BandColor[aBand];
        aStart <- darkBands$StartPoint[aBand];
        aEnd   <- darkBands$EndPoint[aBand];

        posX <- c(RCircos.Pos[aStart:aEnd,1]*outerPos, 
                RCircos.Pos[aEnd:aStart,1]*innerPos);
        posY <- c(RCircos.Pos[aStart:aEnd,2]*outerPos, 
                RCircos.Pos[aEnd:aStart,2]*innerPos);
        polygon(posX, posY, col=aColor, border=NA);
    }
}




#   ==========================================================================
#
#   5.  RCircos.Highligh.Chromosome.Ideogram()
#
#   Highlight chromosomes (draw a color line outside of each chromosome).
#   Graphic device must be initialized first
#
#   Argument:       None
#   Return value:   None
#
#   Example:        RCircos.Highligh.Chromosome.Ideogram()
#

RCircos.Highligh.Chromosome.Ideogram <- function(highlight.pos=NULL,
                highlight.width=NULL)
{
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos  <- RCircos.Get.Plot.Positions()
    RCircos.Par  <- RCircos.Get.Plot.Parameters()

    if(is.null(highlight.pos)) 
        highlight.pos <- RCircos.Par$highlight.pos;
    if(is.null(highlight.width))
        highlight.width <- RCircos.Par$highlight.width;

    chroms <- unique(RCircos.Cyto$Chromosome);
    for(aChr in seq_len(length(chroms)))
    {
        chrRows   <- which(RCircos.Cyto$Chromosome==chroms[aChr]);
        lineStart <- min(RCircos.Cyto$StartPoint[chrRows]);
        lineEnd   <- max(RCircos.Cyto$EndPoint[chrRows]);

        chrColor <- RCircos.Cyto$ChrColor[chrRows[1]];
        lines(RCircos.Pos[lineStart:lineEnd, 1:2]*highlight.pos, 
                col=chrColor, lwd=highlight.width);
    }
}




#   ==========================================================================
#
#   6.  RCircos.Label.Chromosome.Names()
#
#   Label chromosome names.Graphic device must be initialized first.
#
#   Argument:
#       chr.name.pos:   Non-negative numeric, plot position of chromosome
#                       names relative to the center of plot area
#
#   Return value:   None
#
#   Example:        RCircos.Label.Chromosome.Names()
#

RCircos.Label.Chromosome.Names <- function (chr.name.pos=NULL)
{
    RCircos.Par  <- RCircos.Get.Plot.Parameters()
    if(is.null(chr.name.pos)) {
        chr.name.pos <- RCircos.Par$chr.name.pos;
    } else {
        if(chr.name.pos < RCircos.Par$track.in.start)
            stop("Chromosome name positions overlap with inside track.\n");
        if(chr.name.pos > RCircos.Par$track.out.start)
            message("May not plot data tracks outside of chromosomes.\n")
    }
    
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos  <- RCircos.Get.Plot.Positions()
    chroms <- unique(RCircos.Cyto$Chromosome);
    
    rightSide <- nrow(RCircos.Pos)/2;
    for(aChr in seq_len(length(chroms)))
    {
        chrRows <- which(RCircos.Cyto$Chromosome == chroms[aChr]);
        chrStart <- min(RCircos.Cyto$StartPoint[chrRows]);
        chrEnd   <- max(RCircos.Cyto$EndPoint[chrRows]);
        mid <- round((chrEnd - chrStart + 1) / 2, digits=0) + chrStart;

        chrName <- sub(pattern="chr", replacement="", chroms[aChr]);
        text(RCircos.Pos[mid, 1]*RCircos.Par$chr.name.pos,
             RCircos.Pos[mid, 2]*RCircos.Par$chr.name.pos,
             label=chrName, pos=ifelse(mid <= rightSide, 4, 2),
             srt=RCircos.Pos$degree[mid]);
    }
}




#   =========================================================================
#
#   7. RCircos.Ideogram.Tick.Plot()
#
#   Draw ticks along chromosome highlight lines and reset plot parameters
#   with new chromosome name position and outside track start position.
#
#    Argument:    
#       tick.interval: Non-negative integer, distance between two ticks in 
#                      millions base pairs
#
#    Return value:    None
#
#    Example:    RCircos.Tick.Plot(50);
#

RCircos.Ideogram.Tick.Plot <- function(tick.interval=50, track.for.ticks=3)
{
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();

    #   Check if there is enough space for ticks and labels. From outside 
    #   of chromosome ideogram, ticks start at highlight position and take
    #   one track height, tick label takes three tracks, and chromosome 
    #   names use two tracks. There will be total of 6 tracks needed.
    #   ===================================================================
    #
    track.height <- RCircos.Par$track.height;
    tick.height  <- track.height*track.for.ticks;
    ticks.span   <- RCircos.Par$chr.ideo.pos + tick.height*2;
    
    if(RCircos.Par$plot.radius < ticks.span ) 
    { 
        stop(paste0("There is no enough room to draw ticks.\n",
            "Please reset plot radius and redraw chromosome ideogram.\n") ); 
    }

    #   Draw ticks and labels. Positions are calculated based on
    #   chromosome highlight positions
    #   ========================================================
    #
    start.pos <- RCircos.Par$highlight.pos;
    innerPos <- RCircos.Pos[, 1:2]*start.pos;
    outerPos <- RCircos.Pos[, 1:2]*(start.pos + track.height);
    mid.pos  <- RCircos.Pos[, 1:2]*(start.pos + track.height/2);

    the.interval <- tick.interval*1000000;
    short.tick <- round(the.interval/RCircos.Par$base.per.unit, digits=0);
    long.tick  <- short.tick*2;

    lab.pos  <- RCircos.Pos[, 1:2]*(start.pos + tick.height);
    chroms <- unique(RCircos.Cyto$Chromosome);
    for(aChr in seq_len(length(chroms)))
    {
        the.chr  <- RCircos.Cyto[RCircos.Cyto[,1]==chroms[aChr],];
        chr.start <- the.chr$StartPoint[1];
        chr.end   <- the.chr$EndPoint[nrow(the.chr)];

        total.ticks <- ceiling((chr.end-chr.start)/long.tick);
        for(a.tick in seq_len(total.ticks))
        {
            tick.pos <- chr.start + (a.tick-1)*long.tick;
            if(tick.pos < chr.end)
            {
                lines(c(innerPos[tick.pos,1], outerPos[tick.pos,1]), 
                    c(innerPos[tick.pos,2], outerPos[tick.pos,2]),
                    col=the.chr$ChrColor[1]);

                lab.text <- paste0((a.tick-1)*tick.interval*2, "MB");
                text(lab.pos[tick.pos,1] , lab.pos[tick.pos,2], 
                    lab.text, cex=RCircos.Par$text.size,
                    srt=RCircos.Pos$degree[tick.pos]);
            }

            tick.pos <- tick.pos + short.tick;
            if(tick.pos < chr.end)
            {
                lines(c(innerPos[tick.pos,1], mid.pos[tick.pos,1]), 
                    c(innerPos[tick.pos,2], mid.pos[tick.pos,2]), 
                    col=the.chr$ChrColor[1]);
            }
        }
    }
    
    #   Reset plot parameters with new chromosome name position and 
    #   outside track start position. As chr.name.pos is a read-only
    #   parameter, direct work with RCircosEnvironment is needed.
    #   =======================================================

    old.name.pos <- RCircos.Par$chr.name.pos;
    old.out.pos  <- RCircos.Par$track.out.start
    old.distance <- old.out.pos - old.name.pos;
    
    RCircos.Par$chr.name.pos <- ticks.span;
    RCircos.Par$track.out.start <- RCircos.Par$chr.name.pos + old.distance;

    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv());
    RCircosEnvironment[["RCircos.PlotPar"]] <- NULL;
    RCircosEnvironment[["RCircos.PlotPar"]] <- RCircos.Par;
}




#   ==========================================================================
# 
#   8.  RCircos.Reset.Plot.Ideogram()
#
#   Reset chromosome ideogram plot data
#
#   Arguments:
#
#       chromIdeo:  Data frame, object of RCircos cytoband data returned 
#                   from RCircos.Get.Plot.Ideogram(). 
#
#   Returned values: None
#
#   Example:    chromIdeo <- RCircos.Get.Plot.Ideogram();
#               chromIdeo$Location <- round(chromIdeo$Location*0.95);
#               RCircos.Reset.Plot.Ideogram(chromIdeo);
#

RCircos.Reset.Plot.Ideogram <- function(chrom.ideo)
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.Cytoband"]] <- chrom.ideo;
}




#   ==========================================================================
#
#   9.  RCircos.ZoomIn.Chromosome()
#
#   Zoom in one or partial chromosome ideogram before set up core components.
#   When zoom in more than one regions, call this function more times. This
#   function works only with original cyto band table not the plot ideogram
#   object in RCircos environment.
#
#   Arguments:
#
#       ideogram:   Data frame, chromosome ideogram data 
#       chromosome: Character vector, name of chromosome to be zoomed
#       from:       Non-negative integer, genomic coordinates of start  
#                   position of chromosome region to be zoomed
#       to:         Non-negative integer, genomic coordinates, end position
#                    of chromosome region to be zoomed
#       zoomIn:     Non-negative number, fold to zoom in 
#
#   Returned values: Data frame of ideogram with one chromosome region zoomed
#
#   Example:    data(UCSC.HG19.Human.CytoBandIdeogram)
#               oldIdeo <- UCSC.HG19.Human.CytoBandIdeogram
#               zoomedIdeo <- RCircos.ZoomIn.Chromosome(oldIdeo, "chr17",
#                                       10000, 20000, 1000)
#

RCircos.ZoomIn.Chromosome <- function(ideogram=NULL, chromosome=NULL, 
                    from=NULL, to=NULL, zoom.in=NULL)
{
    if(is.null(ideogram) || is.null(chromosome) || 
            is.null(from) || is.null(to) || is.null(zoom.in) )
        stop("Missing function argument(s).");
    if(is.data.frame(ideogram) == FALSE)
        stop("Ideogram data must be in data frame.\n");
    if(length(chromosome) > 1)
        stop("Only one chromosome can be zoomed each time.\n");
    if(from < 0 || to < 0) stop("Numeric arguments cannot be negative.");
    if(from > to) stop("From value is greater the to value.");
    if(zoom.in <= 1) stop("Zoom in fold should be greater than 1.");

    ideo <- RCircos.Validate.Cyto.Info(ideogram);
    chromNames <- unique(as.character(ideo$Chromosome));
    if(!chromosome %in% chromNames)
        stop(message(chromosome, " not found in ideogram."));

    #   which band will be zoomed in
    #   ============================
    #
    chromRows <- which(as.character(ideo$Chromosome)== chromosome);
    startRow <- max(which(ideo$ChromStart[chromRows] <= from));
    endRow   <- min(which(ideo$ChromEnd[chromRows] >= to));
    zoomRows <- chromRows[startRow:endRow];
    newWidth <- (ideo$ChromEnd[zoomRows]-ideo$ChromStart[zoomRows])*zoom.in;
    
    #   scale up the length of the bands to be zoomed. The first start row
    #   will have no change. The band after zoomed rows, if any, are also
    #   need be modified.
    #   ===================================================================
    #
    for(aRow in seq_along(zoomRows))
    {
        ideo$ChromEnd[zoomRows[aRow]] <- newWidth[aRow] +
                        ideo$ChromStart[zoomRows[aRow]];
         ideo$ChromStart[zoomRows[aRow+1]] <- ideo$ChromEnd[zoomRows[aRow]];
    }

    if(endRow < length(chromRows))
    {
        others <- chromRows[-c(1:endRow)];
        differ <- ideo$ChromEnd[chromRows[endRow]] - ideo$ChromStart[others[1]]
        ideo$ChromStart[others] <- ideo$ChromStart[others] + differ;
        ideo$ChromEnd[others]   <- ideo$ChromEnd[others] + differ;
    }

    #   Return zoomed ideogram data
    #   ===========================
    #
    return (ideo);
}




#   ==========================================================================
# 
#   10.  RCircos.ZoomOut.Chromosome()
#
#   Zoom out chromosome ideogram to leave small room between the first and the
#   last chromosome (at 12 O'clock) so that track names could be added.
#
#   Arguments:
# 
#       zoom.out.ratio: A float number between 0 and 1, percentage of target
#                       chromosome ideogram
#
#   Returned value: None. New ideogram object is set to RCircos environment.
#
#   Example:    RCircos.ZoomOut.Chromosome(0.95)
#

RCircos.ZoomOut.Chromosome <- function(zoom.out.ratio=NULL)
{
    if(is.null(zoom.out.ratio)) stop("Zoom out ratio not defined.\n")
    if(zoom.out.ratio >= 1 || zoom.out.ratio < 0) 
        stop("Zoom out ratio should be greatet than 0 and less than 1.");

    ideogram <- RCircos.Get.Plot.Ideogram();
    if(is.null(ideogram)) stop("RCircos ideogram has not been set.\n");

    RCircos.Pos  <- RCircos.Get.Plot.Positions()
    totalUnits <- dim(RCircos.Pos)[1];

    emptySpace <- 1 - zoom.out.ratio;
    spaceTotal <- round(totalUnits*emptySpace, digits=0);
    spaceHalf  <- round(spaceTotal/2, digits=0);

    # Reset chromosome ideogram object
    # ================================
    #
    ideogram$StartPoint <- round(ideogram$StartPoint*zoom.out.ratio, digits=0);
    ideogram$EndPoint   <- round(ideogram$EndPoint*zoom.out.ratio, digits=0);
    ideogram$StartPoint <- ideogram$StartPoint + spaceHalf;
    ideogram$EndPoint <- ideogram$EndPoint + spaceHalf;

    RCircos.Reset.Plot.Ideogram(ideogram);
}




#   =========================================================================
#
#   11. RCircos.Get.Chromosome.Order()
#
#   Generate an ordered chromosome names from input. For human and other
#   mammalian animals, numeric names (integers or Roman numbers) will go 
#   first followed by chromosome X, Y, and M. Character names will be in
#   alphabetical order.  
#
#   Argument:       chromosomes: Character vector, names of chromosomes
#   Return value:   ordered chromosome names
#
#   Example:    chrNames <- paste0("chr", c(14:22, 1:14, "Y", "X"))
#               chrNames <- RCircos.Get.Chromosome.Name.Order(chrNames)
#

RCircos.Get.Chromosome.Order <- function(chromosomes=NULL)
{
    if(is.null(chromosomes)) stop("Missing function argument.\n")
    if(length(chromosomes) != length(unique(chromosomes))) 
        stop("Chromosome names is not an unique set.\n");

    #   Remove prefix "chr" from chromosome names, if any
    # 
    chromNames <- chromosomes;
    chromWithPrefix <- length(grep("^chr", chromosomes));
    if(chromWithPrefix > 0 )
    { 
        if(chromWithPrefix != length(chromosomes)) {
            stop("Not all chromosome name have prefix 'chr'");
        } else { chromNames <- gsub("^chr", "", chromNames); }
    }

    genericChroms <- c(1:100, "X", "Y", "M");
    RomanChroms <- c(as.character(as.roman(1:100)), "M");
    
    totalGeneric <- which(chromNames %in% genericChroms);    
    totalRomans  <- which(chromNames %in% RomanChroms);

    #   Mammalian genomes have both integer and character
    #   chromosome names with or without prefix of "chr"
    #   =================================================
    #
    if(length(totalGeneric) == length(chromosomes))
    {
        orderedChroms <- chromNames[grep("[0-9]", chromNames)];
        orderedChroms <- orderedChroms[order(as.numeric(orderedChroms))];
        
        if("X" %in% chromNames) orderedChroms <- c(orderedChroms, "X");
        if("Y" %in% chromNames) orderedChroms <- c(orderedChroms, "Y");
        if("M" %in% chromNames) orderedChroms <- c(orderedChroms, "M"); 

    } else if(length(totalRomans) == length(chromosomes)) {
    
        chrM <- which(chromNames == "M");
        if(length(chrM) > 0) {
            orderedChroms <- chromNames[-chrM];
        } else { orderedChroms <- chromNames; }
        
        romanIntegers <- as.integer(as.roman(orderedChroms));
        orderedChroms <- orderedChroms[order(romanIntegers)];
        
        if(length(chrM) > 0) orderedChroms <- c(orderedChroms, "M");

    } else {
        message("There is unsupported chromosome names in ideogram\n",
                "and chromosomes are sorted in alphabetical order.")
        orderedChroms <- chromNames[order(chromNames)]; 
    }

    if(chromWithPrefix > 0 ) orderedChroms <- paste0("chr", orderedChroms);

    return (orderedChroms);
}




#   =========================================================================
#
#   12. RCircos.Pseudo.Ideogram.From.Labels()
#
#   Generate a pseudo ideogram from a character list for Circos plot without 
#   cyto band information. Each chromosome will have same length.
#
#   Argument:
#
#       chromsomes: a vector for chromosome names and each element must
#               be unique
#
#   Return value:   Data frame as an ideogram data table
#
#   Example:    geneNames <- paste0("Gene_", 1:20)
#               cyto.info <- RCircos.Pseudo.Ideogram.From.Labels(geneNames)
#

RCircos.Pseudo.Ideogram.From.Labels <- function(chromosomes=NULL)
{
    if(is.null(chromosomes)) stop("Missing function argument.\n");
    if(is.vector(chromosomes) == FALSE) 
        stop("Names for chromosomes must to be a character vector.\n");
    if(length(chromosomes) < 2) 
        stop("Names for chromosomes have length less than 2.\n");

    totalRows <- length(chromosomes)
    if(length(unique(chromosomes)) != totalRows)
        stop("Names for chromosomes must be unique.\n")

    defaultUnits <- RCircos.Get.Default.Circos.Units();
    defaultBases <- RCircos.Get.Default.Base.Per.Units();
    segmentLen <- floor(defaultUnits/totalRows)*defaultBases;

    pseudoIdeogram <- data.frame(Chromosome=chromosomes, 
                        chromStart=rep(0, totalRows),
                        chromEnd=rep(segmentLen, totalRows), 
                        Bands=chromosomes,
                        Stain=rep("gneg", totalRows));

    return (pseudoIdeogram);
}




#   =========================================================================
#
#   13. RCircos.Pseudo.Ideogram.From.Table()
#
#   Generate a pseudo ideogram from plot data for Circos plot. Plot data must
#   have at least two columns serving as chromosome names and band names. An
#   optional gene length or locations can also be provided, e.g.,
#
#   Groups  LocusName   LocusPosition
#   Group1  1           1.5
#   Group1  3           2.0
#   Group1  4           3.0
#
#   Argument:
#
#       plot.data:      A data frame with two or three columns, the first  
#                       column serves as chromosome names
#       location.col:   Column number in plot.data from which end positions 
#                       will be derived for each chromosome or each band.  
#                       This argument must be defined and could be either 
#                       length of the locus or relative locus position.
#       locus.col:      Column number in plot.data which serves as band names, 
#                       this argument is optional
#
#   Return value:   Data frame as an ideogram data table
#
#   Example:    data(RCircos.Gene.Label.Data);
#               plot.data <- RCircos.Gene.Label.Data;
#               cyto.info <- RCircos.Get.Pseudo.Ideogram(plot.data, 2, 3)
#

RCircos.Pseudo.Ideogram.From.Table <- function(plot.data=NULL, 
                location.col=NULL, band.col=NULL)
{
    if(is.null(plot.data)) stop("Missing data argument.\n");
    if(is.null(location.col)) 
        stop("Missing column number for locus location.\n");
    if(is.numeric(location.col) == FALSE)
        stop("The column for locus position must be numeric.\n");

    totalRows <- nrow(plot.data);
    if(totalRows <2) stop("Row number in plot data less than 2.")

    #   If band name is not defined, chromosome names must be unique
    #   and use plot.data[, 1] as band names
    #   ============================================================
    #
    if(is.null(band.col)) {
        if(length(unique(plot.data[, 1]) != totalRows))
            stop("The first column must be unique if band.col is null.\n");
        band.col <- 1;
    }

    #   Scale pseudo genome length to default ideogram (hg19)
    #   ====================================================
    #
    defaultBasePerUnit <- RCircos.Get.Default.Base.Per.Units();
    defaultUnits <- RCircos.Get.Default.Circos.Units();
    defaultTotalBases <- defaultUnits*defaultBasePerUnit;
    
    ideoTotalBases <- sum(plot.data[, location.col]);
    fold <- defaultTotalBases/ideoTotalBases;
    plot.data[, location.col] <- floor(plot.data[, location.col]*fold);

    pseudoIdeogram <- data.frame(
        Chromosome = as.character(plot.data[, 1]),
        ChromStart = rep(0, totalRows), 
        ChromEnd   = as.numeric(plot.data[, location.col]),
        Band = as.character(plot.data[, band.col]),
        Stain = rep("gneg", totalRows)  );

    #   if there are more than one bands on a chromosome, start
    #   and end positions of bands should be continuous.
    #   =================================================================
    #
    chromosomes <- as.character(pseudoIdeogram$Chromosome);
    chromNames <- unique(chromosomes);
    if(length(chromNames) < length(chromosomes) )
    {
        for(aChr in seq_len(length(chromNames)))
        {
            chrRows <- which(chromosomes==chromNames[aChr]);
            theIdeo <- pseudoIdeogram[chrRows, ];
            if(length(chrRows) > 1) 
            {
                band.len <- theIdeo$ChromEnd - theIdeo$ChromStart;
                for(bChr in seq_len(length(chrRows))[-1]){
                    theIdeo$ChromStart[bChr] <- theIdeo$ChromEnd[bChr-1] + 1;
                    theIdeo$ChromEnd[bChr] <- theIdeo$ChromStart[bChr] +
                        band.len[bChr];
                }
            }

            if(aChr==1) {  tmpIdeo <- theIdeo;
            } else { tmpIdeo <- rbind(tmpIdeo, theIdeo); }
        }
        pseudoIdeogram <- tmpIdeo;
    }

    return (pseudoIdeogram);
}




#    =========================================================================
#
#   14. RCircos.Validate.Genomic.Info()
#
#    Check out if the chromosome name, start and end positions match
#    the ideogram data
#
#    Argument:
#
#       genomic.info:   A vector with a chromosome name, start and end 
#                       position on the chromosome
#
#   Return value:   None.
#
#   Example:    RCircos.Validate.Genomic.Info(c("chr1", 10000, 500000))
#

RCircos.Validate.Genomic.Info <- function(genomic.info=NULL)
{
    if(is.null(genomic.info)) 
        stop("Missing plot.data in RCircos.Get.Area.Info().\n");
    if(length(genomic.info) < 2) stop("Missing items in area.info.\n")

    #   Check out if the chromosome is in ideogram data
    #   ===============================================
    #
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    ideoChrom <- as.character(RCircos.Cyto$Chromosome);
    chromosome <- as.character(genomic.info[1]);
    if(!chromosome %in% ideoChrom)
        stop("The chromosome name is not in ideogram data.\n")

    #   Check out if the start and end position is in the
    #   chromosome boundary
    #   ===============================================
    #
    chrRows <- which(ideoChrom == chromosome)
    chrEnd   <- max(as.numeric(RCircos.Cyto$ChromEnd[chrRows]));

    if(is.na(suppressWarnings(as.numeric(genomic.info[2]))))
        stop("Genomic start position must be numeric.\n")
    startPos <- as.numeric(genomic.info[2]);
    if(startPos < 0 || startPos > chrEnd)
       stop("Start position is out of chromosome boundary.\n");

    if(length(genomic.info) >= 3)  {
        if(is.na(suppressWarnings(as.numeric(genomic.info[3]))))
            stop("Genomic start position must be numeric.\n")
        endPos <- as.numeric(genomic.info[3]);
        if(endPos <= 0 || endPos >  chrEnd ) 
            stop("End position is out of chromosome boundary.\n");
        if(startPos > endPos) stop("start position > end position.\n");  
    }
}

#   End of RCircosChromosomeIdeogram.RCircos
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>