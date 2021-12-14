#
#   Source code for RCircos package
#
#   Date created:   January 2, 2013
#   Version:        RCircos v.1.0.0
#
#   Hongen Zhang, Ph.D. (hzhang@mail.nih.gov)
#
#   Genetics Branch
#   Center for Cancer Research 
#   US National Cancer Institute
#   US National Institutes of Health
#   Bethesda, Maryland 20892
# 
#   Last debug done on September 14, 2016
#   Version:    RCircos v.1.2
#   =================================================
#
#   Functions in this file:
#
#   1.  RCircos.Workflow()
#   2.  RCircos.Set.Core.Components()
#   3.  RCircos.Get.Plot.Parameters()
#   4.  RCircos.Get.Plot.Ideogram()
#   5.  RCircos.Get.Plot.Positions()
#   6.  RCircos.Get.Default.Circos.Units()
#   7.  RCircos.Get.Default.Base.Per.Units()
#   8.  RCircos.Get.Padding.Constant()
#   9.  RCircos.Get.Supported.HeatmapColors()
#   10. RCircos.Get.Supported.Plot.Types()
#   11. RCircos.Get.Default.Char.Width()
#   12. RCircos.Get.Default.Text.Size()
#   13. RCircos.Set.Plot.Area()
#   14. RCircos.Multiple.Species.Core.Components()
#   15. RCircos.Get.Plot.Colors()
#   16. RCircos.Get.Link.Colors()
#   17. RCircos.Get.Arrow.Shape()
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>




#   ========================================================================
#
#   Working environment for RCircos to hold RCircos core components and other
#   default global variables
#

RCircos.Env <- new.env();
RCircos.defaultCircosUnits <- 103190;
RCircos.defaultBasePerUnits <- 30000;
RCircos.defaultChromPadding <- 300;
RCircos.paddingConst     <- 300/103190;
RCircos.defaultCharWidth <- 500;
RCircos.defaultTextSize  <- 0.4;

RCircos.heatmapColors <- c("BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed",
                        "GreenBlackRed", "YellowToRed", "BlackOnly");

RCircos.plotTypes <- c("link", "ribbon", "pLink", "polygon", "tile", "ideogram",  
                        "heatmap", "bar", "histogram", "cLine", "vLine",
                        "points",  "connector", "ticks",  "text", "area");


#   ========================================================================
#
#   1.  RCircos.Workflow()
#
#   Print out work flow on screen as a quick guide.
#
#   Argument:       None
#   Return value:   None
#
#   Example:    library(RCircos);
#               RCircos.Workflow();
#

RCircos.Workflow <- function()
{
    helpText <- c("\n1. Load RCircos library:\n\n",
        "   library(RCircos)\n\n",
        "2. Load chromosome cytoband data:\n\n",
        "   data(UCSC.HG19.Human.CytoBandIdeogram);\n",
        "   cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;\n\n",
        "   # Other chromosome ideogram data installed:\n",
        "   # UCSC.Mouse.GRCm38.CytoBandIdeogram\n",
        "   # UCSC.Baylor.3.4.Rat.cytoBandIdeogram\n\n",
        "3. Setup RCircos core components:\n\n",
        "   RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,
                    tracks.inside=10, tracks.outside=0);\n\n",
        "4. Load input data:\n\n",
        "   heatmap.data <- read.table(\"/path/Heatmap.data.txt\", 
                        sep=\"\\t\", quote=\"\", head=T);\n",
        "   hist.data <- read.table(\"/path/histgram.data.txt\", 
                        sep=\"\\t\", quote=\"\", head=T);\n",
        "   link.data <- read.table(\"/path/link.data.txt\", 
                        sep=\"\\t\", quote=\"\", head=T);\n\n",
        "5. Modify plot parameters if necessary:\n\n",
        "   rcircos.params <- RCircos.Get.Plot.Parameters()\n",
        "   rcircos.params$radiu.len <- 1.5;\n",
        "   RCircos.Reset.Plot.Parameters(rcircos.params);\n\n",
        "6. Open graphic device:\n\n",
        "   RCircos.Set.Plot.Area();\n\n",
        "   or submit your own code. For example: \n\n",
        "   par(mai=c(0.25, 0.25, 0.25, 0.25));\n",
        "   plot.new();\n",
        "   plot.window(c(-2.5,2.5), c(-2.5, 2.5));\n\n",
        "7. Call plot function to plot each data track:\n\n",
        "   RCircos.Chromosome.Ideogram.Plot();\n",
        "   RCircos.Heatmap.Plot(heatmap.data, data.col=5, 
                        track.num=1, side=\"in\");\n",
        "   RCircos.Histogram.Plot(hist.data, data.col=4, 
                        track.num=4, side=\"in\");\n",
        "   RCircos.Link.Plot(link.data, track.num=5, 
                        by.chromosome=FALSE);\n\n",
        "8. Close the graphic device if you was plotting to file:\n\n",
        "   dev.off();\n\n"               
        );

    message(helpText);
}




#   ========================================================================
#
#   2.  RCircos.Set.Core.Components()
#
#   Initialize RCircos core components including of:
#
#   1.  Data frame for chromosome ideogram data
#   2.  x and y coordinates for a circular line
#   3.  Parameters to control the RCircos plot
#  
#   Argument:
#
#       cyto.info:      A data frame contains chromosome ideogram data
#       chr.exclude:    Character vector, chromosome name(s) to be excluded
#       tracks.inside:  Non-negative integer, number of data tracks inside of 
#                       chromosome ideogram
#       tracks.outside: Non-negative integer, number of data tracks outside of 
#                       chromosome ideogram
#
#   Return value:   None. 
#
#   Example:    library(RCircos);
#               RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL, 
#                       tracks.inside=10, tracks.outside=0);
#

RCircos.Set.Core.Components<-function(cyto.info=NULL, chr.exclude=NULL, 
            tracks.inside=10, tracks.outside=0) {

    if(tracks.inside<0 || tracks.outside<0) 
    { stop("Track number cannot be smaller than 0.\n") }

    #   Step 1. validate cyto.info for correct chromosome start and
    #   end positions of each chromosome band. The data will not be 
    #   hold for the RCircos environment
    #   ===========================================================
    #
    cytoBandData <- RCircos.Validate.Cyto.Info(cyto.info, chr.exclude)

    #   Step 2. Initialize RCircos core components
    #   ===========================================
    #
    RCircos.Initialize.Plot.Parameters(tracks.inside,tracks.outside);
    RCircos.Set.Cytoband.Data(cytoBandData);
    RCircos.Set.Base.Plot.Positions();

    #   User friendly notice.
    #   ===============================================
    #
    message("\nRCircos.Core.Components initialized.\n",
            "Type ?RCircos.Reset.Plot.Parameters to see",
            " how to modify the core components.\n\n");
}




#   ========================================================================
#
#   3 ~ 12. Methods to retrieve RCircos core components stored in RCircos 
#           environment. No calculations in each method.
#
#   Argument:       None
#   Return value:   Parameters, cytoband, or default plot positions stored in  
#                   RCircos environment.
#
#   Example:        plot.param <- RCircos.Get.Plot.Parameters();
#                   plot.cyto  <- RCircos.Get.Plot.Ideogram();
#                   plot.pos   <- RCircos.Get.Plot.Positions();
#

RCircos.Get.Plot.Parameters<-function()
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    return (RCircosEnvironment[["RCircos.PlotPar"]]);
}

RCircos.Get.Plot.Ideogram<-function()
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    return (RCircosEnvironment[["RCircos.Cytoband"]]);
}

RCircos.Get.Plot.Positions<-function()
{
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    return (RCircosEnvironment[["RCircos.Base.Position"]]);
}

RCircos.Get.Default.Circos.Units <- function()
{
    return (RCircos.defaultCircosUnits);
}

RCircos.Get.Default.Base.Per.Units <- function()
{
    return (RCircos.defaultBasePerUnits);
}

RCircos.Get.Default.Chrom.Padding <- function()
{
    return(RCircos.defaultChromPadding);
}

RCircos.Get.Padding.Constant <- function()
{
    return (RCircos.paddingConst);
}

RCircos.Get.Supported.HeatmapColors <- function()
{
    return (RCircos.heatmapColors);
}

RCircos.Get.Supported.Plot.Types <- function()
{
    return (RCircos.plotTypes);
}

RCircos.Get.Default.Char.Width <- function()
{
    return (RCircos.defaultCharWidth);
}

RCircos.Get.Default.Text.Size <- function()
{
    return (RCircos.defaultTextSize);
}



#   =========================================================================
#
#   13. RCircos.Set.Plot.Area()
#
#   Open a new graphic device with the build-in plot.radius. If user want an 
#   image file, correct image type must be created and closed from command
#   line or other script. This function call can also be replaced with script
#   from command line if user known how much of area are needed.
#
#   Argument:       None
#   Return value:   None
#
#   Example:    RCircos.Set.Plot.Area();
#

RCircos.Set.Plot.Area <- function(margins=0.25, ...)
{
    if(!is.numeric(margins) || margins < 0)
        stop("Incorrect margins value.\n");
        
    RCircos.Par <- RCircos.Get.Plot.Parameters();

    par(mai=c(margins, margins, margins, margins));
    plot.new();
    plot.window(c(-1*RCircos.Par$plot.radius, RCircos.Par$plot.radius), 
                c(-1*RCircos.Par$plot.radius, RCircos.Par$plot.radius));
}




#   ========================================================================
# 
#   14. RCircos.Multiple.Species.Core.Components()
#
#   Generate customized RCircos core components for multiple species plot 
#   with all or part of chromosomes
#
#   Arguments:
#
#       cyto.info.list: List of multiple chromosome ideogram data 
#       species:        Character vector for prefix of chromosome names to 
#                       identify different species
#       chr.exclude:    Chromosomes which should be excluded from dataset
#       tracks.inside:  Non-negative integer, number of data tracks inside 
#                       of chromosome ideogram
#       tracks.outside: Non-negative integer, number of data tracks outside 
#                       of chromosome ideogram
#
#   Return value:   None
#
#   Example:    ideos <- list(mouse.cyto, rat.cyto);
#               species   <- c("m", "r");
#               RCircos.Multipl.Species.Core.Components(cyto.info.list=ideos, 
#                            species.list=species, chr.exclude=NULL, 
#                            tracks.inside=10, tracks.outside=0);
#

RCircos.Multiple.Species.Core.Components <- function(cyto.info.list=NULL, 
    species.list=NULL, chr.exclude=NULL, tracks.inside=10, tracks.outside=0)
{
    #   cyto.info.list and species must have same length
    #   ================================================
    #
    if(is.null(cyto.info.list) || is.null(species.list))
        stop("Missing cyto info or species list.\n");
    if(length(cyto.info.list) != length(species.list)) 
        stop("Number of cyto band data and species must be same");

    if(!is.numeric(tracks.inside) || !is.numeric(tracks.outside))
        stop("tracks.inside and tracks.outside must be numeric.\n");
    if(tracks.inside<0 || tracks.outside<0) 
        stop("Track.inside and track.outside cannot be negative.\n");

    #   Validte each chromosome ideogram data then combine them as one
    #   ==============================================================
    #
    numOfSpecies <- length(cyto.info.list);
    for(aCyto in seq_len(numOfSpecies))
    {
        cytoInfo <- data.frame(cyto.info.list[aCyto]);
        prefix <- species.list[aCyto];

        cytoInfo <- RCircos.Validate.Cyto.Info(cytoInfo, NULL);
        cytoInfo$Chromosome <- paste(prefix, cytoInfo$Chromosome, sep="");

        if(aCyto==1) { 
            newCytoInfo <- cytoInfo;
        } else {
            newCytoInfo <- rbind(newCytoInfo, cytoInfo);
        }
    }

    #   Initialize RCircos core components
    #   ==================================
    #
    RCircos.Initialize.Plot.Parameters(tracks.inside,tracks.outside);
    RCircos.Set.Cytoband.Data(newCytoInfo);
    RCircos.Set.Base.Plot.Positions();
}




#   ========================================================================
# 
#   15. RCircos.Get.Plot.Colors()
#
#   Retrieve plot color for each data point which is stored in plot data
#
#   Arguments:
#
#       plot.data:  Data frame of genomic data
#       color:      character vector, default color names
#
#   Returned value: Vector of color names with length same as number of rows 
#                   of plot data
#   Example:    plotColors <- RCircos.Get.Plot.Colors(plotData, "red")
#
#   Last modified on June 30, 2015
#

RCircos.Get.Plot.Colors <- function(plot.data=NULL, color="black")
{
    if(is.null(plot.data))
        stop("Missing argument in RCircos.Get.Plot.Colors().\n");
    
    colorCol <- grep("PlotColor", colnames(plot.data));

    if(length(colorCol)==0) {
        plotColors <- rep(color, nrow(plot.data));
    } else if(length(colorCol)==1) {
        plotColors <- as.character(plot.data[, colorCol]);
    } else {
        stop("Incorrect plot colors defined in dataset.\n"); 
    }
        
    return (plotColors)
}




#   ========================================================================
#
#   16. RCircos.Get.Link.Colors()
#
#   Get plot color for each link line or ribbon
#
#   Arguments:
#
#       plot.data:      Data frame of paired genomic position data
#       chrom.columns:  Numeric vector of 2, column numbers for chromosome
#                       names in link data.
#       by.chromosome:  Logic, if TRUE, read color is used for links on same 
#                       chromosome and blue color is used for link between 
#                       different chromosomes. If FALSE, user defined colors 
#                       or rainbow colors will be used.
#
#   Returned value:     vector of color names with length same as number of 
#                       rows in input data
#
#   Example:        linkColor <- RCircos.Get.Link.Colors(link.data, TRUE)
#

RCircos.Get.Link.Colors <- function(link.data, genomic.columns=3,
                        by.chromosome=TRUE)
{
    if(is.null(link.data)) stop("Missing link data.\n");
    if(genomic.columns < 2 || genomic.columns > 3)
        stop("Genomic position columns must be 2 or 3.\n")

    redColor  <- rgb(1, 0, 0, alpha=0.5);
    blueColor <- rgb(0, 0, 1, alpha=0.5);

    #   Default colors
    #   =======================================
    #
    linkColors <- rep(blueColor, nrow(link.data));

    #   If by.chromosome is set to true, red color will be used 
    #   for links in same chromosome and blue color for links 
    #   between different chromosomes
    #   =======================================================
    #
    if(by.chromosome==TRUE) 
    {
        if(genomic.columns == 2) {
            chromColumns <- c(1, 3)
        } else { chromColumns <- c(1, 4); }

        startChroms <- as.character(link.data[, chromColumns[1]]);
        endChroms   <- as.character(link.data[, chromColumns[2]]);

        for(aRow in seq_len(nrow(link.data)))
        {
            if(startChroms[aRow] == endChroms[aRow]) 
                linkColors[aRow] <- redColor;
        }
        #   If the plot color is provided in dataset, use it to
        #   replace the default one (rainbow)
        #   ====================================================
        #
        } else {
            colorCol <- grep("PlotColor", colnames(link.data));

            if(length(colorCol==1)) 
            {
                theColor <- as.character(link.data[, colorCol]);
                for(aRow in seq_len(length(theColor)))
                {
                    rgbVal <- as.vector(col2rgb(theColor[aRow]))/255;
                    linkColors[aRow] <- rgb(red=rgbVal[1], 
                        green=rgbVal[2], blue=rgbVal[3], alpha=0.5);
                }
            } else {
                for(aRow in seq_len(nrow(link.data)))
                {
                    rgbVal <- as.vector(col2rgb(aRow+1))/255;
                    linkColors[aRow] <- rgb(red=rgbVal[1], 
                        green=rgbVal[2], blue=rgbVal[3], alpha=0.5);
                }
            }
        }

        return (linkColors);
}




#   ========================================================================
#
#   17. RCircos.Get.Arrow.Shape()
#
#   Get default coordinates for an arrow shape. The arrow is represented
#   as a polygon inside of a circle with radius of 1.
#
#   Argument:       None
#   Return value:   A two dimensional numeric matrix for x and y coordinates
#                   of polygon.
#
#   example:    arrow <- RCircos.Get.Arrow.Shape()
#
#   Last updated on July 15, 2015
#

RCircos.Get.Arrow.Shape <- function(side="in")
{
    coor.x <- c( 0, -0.7, -0.2,  -0.2,  0.2, 0.2,  0.7,  0);
    coor.y <- c(-1,  0.7,  0.4,   1,    1,    0.4,  0.7, -1);
  
    if(side == "in") coor.y <- coor.y*-1;
    arrow <- cbind(coor.x, coor.y);
    
    return (arrow);
}

#   End of RCircosMain.R
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>