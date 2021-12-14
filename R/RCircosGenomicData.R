#
#   Functions for plot data preparation
#
#   1.  RCircos.Validate.Genomic.Data()
#   2.  RCircos.Get.Single.Point.Positions()
#   3.  RCircos.Get.Paired.Points.Positions()
#   4.  RCircos.Multiple.Species.Dataset()
#   5.  RCircos.Sort.Genomic.Data()
#   6.  RCircos.Get.Data.Point.Height()
#   7.  RCircos.Get.Plot.Layers()
#
#   Last debug done on September 14, 2016
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>


#   ========================================================================
#
#   1.  RCircos.Validate.Genomic.Data()
#
#   Validate input dataset for correct chromosome names, chromosome start, 
#   and chromosome end positions. Chromosome names will be converted to 
#   character vectors if they are factor variables.
#
#   Arguments:
#
#       genomic.data:   Data frame with genomic position data. Sorting is 
#                       not required.
#       genomic.columns: Non-negative integer, total number of columns for
#                       genomic position (chromosome name, start and/or end
#                       position)
#       plot.type:      Character vector, either "plot" or "link"
#
#   Return value:       None. If there is error, exit.
#
#   Example:    RCircos.Validate.Genomic.Data(genomic.data, "plot")
#

RCircos.Validate.Genomic.Data <- function(genomic.data=NULL,
            plot.type=c("plot", "link"), genomic.columns=3) 
{
    if(is.null(genomic.data)) stop("Missing genomic data.\n");
    if(genomic.columns < 2 || genomic.columns > 3)
        stop("Incorrect number of genomic position columns defined.\n");

    #   Plot data has only one chromosome column and link data has two
    #   ===============================================================
    #
    plot.type <- tolower(plot.type);
    if(plot.type=="plot") { chromCol <- 1;
    } else if(plot.type=="link") { 
            chromCol <- c(1, genomic.columns + 1);
    } else { stop("Plot type must be \"plot\" or \"link\"") }
    
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
    cytoChroms <- unique(as.character(RCircos.Cyto$Chromosome));
    for (aCol in seq_len(length(chromCol)))
    {
        #   Make sure chromosomes in input data are all included 
        #   in chromosome ideogram data
        #   ====================================================
        #
        theCol <- chromCol[aCol];
        dataChroms <- unique(as.character(genomic.data[, theCol]));

        if(sum(dataChroms %in% cytoChroms) < length(dataChroms)) 
            stop("Some chromosomes in plot data is not in ideogram.");

        #   Make sure chromosome start and end positions in genomic
        #   data are not negative.
        #   ==============================================
        #
        if(min(genomic.data[, theCol+1]) < 0) 
        { stop("One or more chromStart position less than 0."); }

        #   if there are three columns for genomic positions, end
        #   position must be greater or equal to start position
        #   ======================================================
        if(genomic.columns == 3)
        { 
            if(min(genomic.data[, theCol+2]) < 0)
                stop("One or more chromEnd position less than 0."); 

            #   Make sure in genomic data all chromosome start positions 
            #   are smaller than their paired chromosome end positions
            #   ========================================================
            #
            startPos <- as.numeric(genomic.data[, theCol+1]);
            endPos <- as.numeric(genomic.data[, theCol+2]);
            if(length(which(endPos<startPos))>0) 
                stop("\nOne or more chromStart greater than chromEnd.\n");
        }

        #   Make sure chromosome start and end locations in genomic
        #   data are not out of chromosome length
        #   ========================================================
        # 
        startCol <- theCol + 1;
        if(genomic.columns == 3) {
            endCol <- startCol + 1;
        } else { endCol <- startCol; }
        
        for(aChr in seq_len(length(dataChroms)))
        {
            theChr <- dataChroms[aChr];
            inData <- genomic.data[genomic.data[,theCol] == theChr,];

            #   Be careful for the cases of grep("1", "1")  and ("1", "12")
            #   ===========================================================
            #
            cytoData <- RCircos.Cyto[grep(paste(theChr, "$", sep=""), 
                                RCircos.Cyto$Chromosome),];

            if(max(inData[, startCol]) > max(cytoData[,3]) ||
               max(inData[, endCol]) > max(cytoData[,3]))
            {
                stop("One or more genomic position in plot data is\n",
                    "outside of chromosome length for ", theChr, ".\n"); 
            }
        }
    }
}




#   ========================================================================
#
#   2.  RCircos.Get.Single.Point.Positions()
#
#   Convert input genomic data to plot data by adding x and y coordinates 
#   for each row of a data set. A set of points for a circle is held in the 
#   RCircos session. We only need the index of the point for each chromosome 
#   position
#
#   Auguments:
#
#       genomic.data:   A data frame contains genomic positions (at least
#                       two or three columns for chromosome names, start 
#                       and/or end positions). It does not need to be sorted.
#       genomic.columns:    Non-negative integer, total number of columns for
#                       genomic position (chromosome name, start and/or end
#                       position).
#
#   Returned value:     A data frame same as input but with a new column
#                       for index of plot positions on circular line. 
#
#   Example:    data(RCircos.Heatmap.Data)
#               plot.data<-RCircos.Get.Track.Plot.Position(RCircos.Heatmap.Data)
#
#   Last revised on July 6, 2015
#

RCircos.Get.Single.Point.Positions <- function(genomic.data=NULL, 
                genomic.columns=3)
{
    if(is.null(genomic.data))
        stop("Missing genomic.data in RCircos.Get.Track.Plot.Position().\n");

    #   Check chromosome names, chromStart, and chromEnd positions, 
    #   if failed, function will exit here. No more processing.
    #   ==========================================================
    #
    RCircos.Validate.Genomic.Data(genomic.data, "plot", genomic.columns);

    #   Calculate the point index for each chromosome location. If both
    #   start and end positions are defined, the location will be in the
    #   mid of start and end positions
    #   _______________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dataPoints <- rep(0, nrow(genomic.data));
    
    if(genomic.columns == 3) {
         location <- round((genomic.data[, 2] + genomic.data[, 3])/2, 
                        digits=0);
    } else if(genomic.columns == 2) {
        location <- as.numeric(genomic.data[, 2]);
    } else { stop("Incorrect genomic column numbers.\n"); }

    for(aRow in seq_len(nrow(genomic.data)))
    {
        chromosome <- as.character(genomic.data[aRow, 1]);
        dataPoints[aRow] <- RCircos.Data.Point(chromosome, location[aRow]);
    }
    genomic.data["Location"] <- dataPoints;

    #   Sort the data by plot position
    #   ===============================================
    #   
    genomic.data <- genomic.data[order(genomic.data$Location),];

    #   The data needs to be held for the RCircos session
    #    =================================================
    #
    return (genomic.data);
}




#    ========================================================================
#
#   3.  RCircos.Get.Paired.Points.Positions()
#
#    Convert genomic link data to plot data by adding two columns for start
#    and end plot positions as start and end of link line/ribbon/hLine. A
#    set of points for a circle is held in the RCircos session. We only need
#    the index of the point for each chromosome position.
#
#    Auguments:
#
#       genomic.data:       A data frame contains paired genomic positions (
#                           chromosome names, start and end positions). The  
#                           data does not need to be sorted.
#       genomic.columns:    Non-negative integer, total number of columns for
#                           genomic position (chromosome name, start and/or
#                           end position).
#       plot.type:          Chraracter vector, either "link", "ribbon", "pLink", 
#                           "polygon", or "tile",
#
#    Returned value:    A data frame same as input but with two new columns
#                       for index of plot positions on circular line. 
#
#    Example:   data(RCircos.Heatmap.Data)
#               linkData<-RCircos.Get.Paired.Points.Positions(RCircos.Link.Data,
#                       genomic.columns=3, plot.type="link")
#
#    Last revised on June 29, 2015
#

RCircos.Get.Paired.Points.Positions <- function(genomic.data=NULL, 
            genomic.columns=3, 
            plot.type=c("link", "ribbon", "pLink", "polygon", "tile"))
{
    if(is.null(genomic.data) || is.null(plot.type))
        stop("Missing genomic.data or plot.type in ",
            "RCircos.Get.Paired.Points.Positions().\n");

    #    Calculate the point index for each chromosome location
    #    _______________________________________________________________
    #    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    linkStart <- rep(0, nrow(genomic.data));
    linkEnd <- rep(0, nrow(genomic.data));

    startChrCol <- 1;
    startPosCol <- 2;
    
    #   Plot types are "pLink", "polygon", "tile", which need two points
    #   on same chromosome. The first point uses start position and the 
    #   second point uses the end position.
    #
    supportedTypes <- RCircos.Get.Supported.Plot.Types();
    if(plot.type %in% supportedTypes[3:5]) 
    {
        if(genomic.columns != 3) 
            stop("Incorrect column number for genomic position.\n");
        RCircos.Validate.Genomic.Data(genomic.data, "plot", genomic.columns);
        endChrCol <- 1; endPosCol <- 3;

    #   plot types are "link" and "ribbon" which need two points but not 
    #   necessarily be on same chromosome. Position columns could be 2
    #   and 4 or 2 and 5.
    #
    } else if(plot.type %in% supportedTypes[1:2]){ 
        if(genomic.columns == 2) {
            endChrCol  <- 3; endPosCol  <- 4; 
        } else if(genomic.columns == 3){ 
            endChrCol  <- 4; endPosCol  <- 5; 
        } else { stop("Incorrect number for genomic position.\n"); }
        
        RCircos.Validate.Genomic.Data(genomic.data, "link", genomic.columns);
    } else { stop("Incorrect plot type for paired positions."); }

    for(aRow in seq_len(nrow(genomic.data)))
    {
        chromosome <- as.character(genomic.data[aRow,  startChrCol]);
        location <- as.numeric(genomic.data[aRow, startPosCol]);
        linkStart[aRow] <- RCircos.Data.Point(chromosome, location);

        chromosome <- as.character(genomic.data[aRow, endChrCol]);
        location <- as.numeric(genomic.data[aRow, endPosCol]);
        linkEnd[aRow] <- RCircos.Data.Point(chromosome, location);    
    }
    genomic.data["LinkStart"] <- linkStart;
    genomic.data["LinkEnd"] <- linkEnd;

    #    Sort the data by plot position
    #    ===============================================
    #    
    genomic.data <- genomic.data[order(genomic.data$LinkStart),];

    #    The data needs to be held for the RCircos session
    #    =================================================
    #    
    return (genomic.data);
}




#   ========================================================================
#
#   4.  RCircos.Multiple.Species.Dataset()
#
#   Combine and modify the chromosome names in multiple species datasets to 
#   match the chromosomes in multiple species cytoband data
#
#   Arguments:
#
#       data.list:  List of genomic data from multiple species
#       species:    Character vector for prefix of chromosome names to 
#                   identify different species
#
#       Note:   The order of each dataset in data list and species in 
#               species list must match.
#
#   Return value:   A data frame contain all data in the input data list with 
#                   chromosome names modified
#
#   Example:    dataSets <- list(mouse.data, rat.data)
#               species.list  <- c("m", "r")
#               dataset <- RCircos.Get.Multiple.Species.Dataset(dataSets, 
#                               species.list)
#

RCircos.Multiple.Species.Dataset <- function(data.list, species)
{
    #   Number of datasets and species must be same
    #   ===========================================
    #
    if(length(data.list) != length(species)) 
    { stop("Error! Number of datasets and species must be same") }

    #   Modify chromosome names in each dataset then combine them as one
    #   ================================================================
    #
    numOfData <- length(data.list);
    numOfColumns <- ncol(data.frame(data.list[1]));

    for(aData in seq_len(numOfData))
    {
        dataset <- data.frame(data.list[aData]);
        prefix <- species[aData];

        #   Number of columns of each dataset must be same
        #   ==============================================
        #
        if(ncol(dataset)!= numOfColumns) 
        { stop("Error! Datasets have different columns.") }

        dataset[,1] <- paste(prefix, dataset[,1], sep="")
        if(aData==1) { 
            newData <- dataset;
        } else { newData <- rbind(newData, dataset) }
    }

    return (newData);
}




#   =========================================================================
#
#   5.  RCircos.Sort.Genomic.Data()
#
#   Sort genomics/ideogram data. The order of chromosome names should be
#   numeric names (integers or Roman numbers) first then character names.
#   If chromosome names are all characters alphabets order will be used. This
#   function could be used before RCircos plot.
#
#   Argument:
#
#       genomic.data:   A data frame with the first three columns for 
#                       chromosome names, start and end positions. If it is 
#                       ideogram data, next two columns must be band names, 
#                       and Giemsa stain status
#       is.ideo:        Logic, if the data is ideogram or not
#
#   Return value:   The input data ordered by chromosome names then start
#                   positions
#
#   Example:    ideogram <- RCircos.Sort.Genomic.Data(cyto.info, TURE);
#               line.data <- RCircos.Sort.Genomic.Data(line.data, FALSE);
#
#    Last debugged on June 27, 2015
#

RCircos.Sort.Genomic.Data <- function(genomic.data=NULL, is.ideo=FALSE)
{

    if(is.null(genomic.data)) 
        stop("Missing argument in RCircos.Sort.Genomic.Data().\n");
    if(!is.data.frame(genomic.data))
        stop("Input data must be in data frame.\n");

    #   Check out the ideogram data since this function could
    #   be called before RCircos core component initialization
    # 
    if(is.ideo) {
        if( ncol(genomic.data) <5 ) {
            stop( paste("Cytoband data must have columns for Chromosome,",
                    "ChromStart, ChromEnd, Band, Stain\n",
                    "Current cyto.info columns:", ncol(genomic.data)));
        } else {
                colnames(genomic.data)[1:5] <- c("Chromosome", "ChromStart", 
                "ChromEnd", "Band", "Stain");
        }
    } else {
        if( ncol(genomic.data) < 2 )
            stop("Genomic data must have two or more columns.\n");
        colnames(genomic.data)[1:2] <- c("Chromosome", "ChromStart");
    }

    #   Get chromosome order
    #   ========================================================
    #
    chromosomes <- unique(genomic.data$Chromosome);
    if(length(chromosomes)>=2)
        chromosomes <- RCircos.Get.Chromosome.Order(chromosomes);

    #   Reconstruct ideogram data and sort by chromStart for 
    #   each chromosome
    #   =====================================================
    #
    rows <- which(genomic.data$Chromosome %in% chromosomes[1]);
    sortedData <- genomic.data[rows, ];
    sortedData <- sortedData[order(sortedData$ChromStart), ];
    
    for(aChr in seq_len(length(chromosomes))[-1])
    {
        rows <- which(genomic.data$Chromosome %in% chromosomes[aChr]);
        theData <- genomic.data[rows, ];
        theData <- theData[order(theData$ChromStart), ];
    
        sortedData <- rbind(sortedData, theData);
    }

    return (sortedData);
}




#    =========================================================================
#
#   6.  RCircos.Get.Data.Point.Height()
#
#   Calculate data point height inside a plot track such as scatter location,
#   top or bottom location of a bar, layer of a tile or parallel link line.
#   Note: if user does not provide min.value and max.value, the smallest value
#   in data will be plot on the bottom border of data track and the highest one 
#   will be on the top border of data track.
#
#    Argument:
#
#       plot.values:    Numeric, the data to be plotted on a data track
#       min.value:       Numeric, the minimum value of data range
#       max.value:       Numeric, the maximum value of data range
#       plot.type:      Character vector, plot type, valid values are
#                       "bar", "histogram", "uniform", or "points"
#       height.range:   Non-negative numeric, height of plot track
#
#   Return value:       Numeric vector with values between 0 ~ 1
#   Example:    data.values <- runif(1000, -4, 11)
#               data.height <- RCircos.Get.Data.Point.Height(data.values, 
#                           min.value=-4, max.value=14, plot.type="point", 
#                           height.range=NULL)
#

RCircos.Get.Data.Point.Height <- function(plot.values=NULL, min.value=NULL, 
            max.value=NULL, plot.type=NULL, track.height=NULL)
{
    if(is.null(plot.values) || is.null(plot.type))
        stop("Missing argument in RCircos.Get.Data.Point.Height().\n");
    if(!is.vector(plot.values) || !is.numeric(plot.values))
        stop("Plot values must be in a numeric vector.\n");
    
    supportedTypes <- RCircos.Get.Supported.Plot.Types();
    if(!plot.type %in% supportedTypes[c(3, 8, 9, 10, 12)]) 
        stop("The plot type has no height property.\n");

    #   If height.range is not defined from customized track 
    #   height, use track height for height range base.
    #   ===================================================
    #
    if(is.null(track.height)) {
        RCircos.Par <- RCircos.Get.Plot.Parameters();
        track.height  <- RCircos.Par$track.height;
    } 

    #   when plot values is between 0 and 1. No converting needed.
    #   =========================================================
    if(min(plot.values) >= 0 && max(plot.values) <= 1) {
        dataHeight <- plot.values*track.height;

    #   Converting plot values to ratios to track.height
    #   ===========================================================
    } else {
        if(is.null(min.value) || is.null(max.value)) {
            min.value <- min(plot.values);
            max.value <- max(plot.values);
        } else { 
            if(is.numeric(min.value) == FALSE || 
                is.numeric(max.value) == FALSE || 
                length(min.value) > 1 || length(max.value) > 1)
                stop("Values for data height must be one numeric value.\n");

            outliers <- which(plot.values > max.value);
            if(length(outliers) > 0 ) plot.values[outliers] <- max.value;

            outliers <- which(plot.values < min.value);
            if(length(outliers) > 0 ) plot.values[outliers] <- min.value;
        } 

        dataScale <- max.value - min.value;
        dataHeight <- (plot.values-min.value)/dataScale*track.height;
    }

    return (dataHeight);
}




#   ===========================================================================
#
#   7.  RCircos.Get.Plot.Layers()
#
#   Check out overlaps between different genomic positions on same chromosome
#   and get layer numbers for each line
#
#   Argument:   
#       genomic.data:   A data frame with genomic positions (chromosomes,  
#                       start and end positions) and the positions should  
#                       be already validated and sorted by chromosome then 
#                       start position.
#       genomic.columns: Non-negative integer, total number of columns for
#                       genomic positions.
#
#   Return value:   A non-negative integer vector with length same as the
#                   total rows of input data
#
#   Example:    data(RCircos.Tile.Data)
#               layers <- RCircos.Get.Plot.Layers(RCircos.Tile.Data, 3);
#

RCircos.Get.Plot.Layers <- function(genomic.data=NULL, genomic.columns=NULL) 
{
    if(is.null(genomic.data)) 
        stop("Missing argument in RCircos.Check.Position.Overlaps().\n");
    if(is.null(genomic.columns) || genomic.columns != 3)
        stop("Start and end position are needed for layer assignment.\n");

    theLayer <- 1;
    theChr   <- as.character(genomic.data[1, 1]);
    theStart <- as.numeric(genomic.data[1, 2]);  
    theEnd   <- as.numeric(genomic.data[1, 3]);

    segLayers <- rep(1, nrow(genomic.data));
    for(aRow in seq_len(nrow(genomic.data))[-1])
    {
        #   Meet a new region without overlap with previous or
        #   a different chromosome, reset relevant variables
        #   ==================================================
        #
        if (genomic.data[aRow, 2] >= theEnd ) {
            theLayer <- 1; 
            theStart <- genomic.data[aRow, 2];
            theEnd   <- genomic.data[aRow, 3];
        } else if (genomic.data[aRow, 1] != theChr) {
            theLayer <- 1;  
            theChr   <- genomic.data[aRow, 1];
            theStart <- genomic.data[aRow, 2];
            theEnd   <- genomic.data[aRow, 3];
        } else {  
            theLayer <- theLayer + 1; 
            if(genomic.data[aRow, 3] > theEnd) 
            { theEnd <- genomic.data[aRow, 3]; }
        }
        segLayers[aRow] <- theLayer;
    }

    return(segLayers);
}


#   End of RCircosGenomicData.R
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>