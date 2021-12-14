#   Functions for heatmap plot
#
#   1.  RCircos.Get.Heatmap.Color.Scale()
#   2.  RCircos.Plot.Heatmap.Color.Scale()
#   3.  RCircos.Get.Heatmap.Data.Colors()
#   4.  RCircos.Get.Heatmap.Color.Scale.Location()
#
#   Last debugged on September 14, 2016
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>


#   ========================================================================
# 
#   1.  RCircos.Get.Heatmap.Color.Scale()
#
#   Create color map for heatmap plot
#
#   Arguments:  heatmap.color, character vector, one of 
#
#       BlueWhiteRed:   colors from blue to white then red
#       GreenWhiteRed:  colors from green to white then red
#       GreenYellowRed: colors from green to yellow then red
#       GreenBlackRed:  colors from green to black then red
#       YellowToRed:    colors from yellow to red
#       BlackOnly:      default
#
#   Example:    RCircos.Get.Heatmap.Color.Scale("BlueWhiteRed")
#

RCircos.Get.Heatmap.Color.Scale <- function(heatmap.color=NULL)
{
    if(is.null(heatmap.color))
        stop("Missing argument in RCircos.Get.Heatmap.Color.Scale().\n");

    allOnes  <- seq(1, 1, length=256);
    allZeros <- seq(0, 0, length=256);
    one2zeor <- seq(1, 0, length=256);
    zero2one <- seq(0, 1, length=256);

    #   Blue, White, and Red
    #   ==============================================
    #
    if(heatmap.color=="BlueWhiteRed") {

        RedRamp  <- rgb(allOnes, one2zeor, one2zeor);
        BlueRamp <- rgb(zero2one, zero2one, allOnes);
        ColorRamp <- cbind(BlueRamp, RedRamp);

    #   Green, White, and Red
    #   ==============================================
    #
    } else if (heatmap.color=="GreenWhiteRed") {

        RedRamp   <- rgb(allOnes, one2zeor, one2zeor);
        GreenRamp <- rgb(zero2one, allOnes, zero2one);
        ColorRamp <- cbind(GreenRamp, RedRamp);

    #   Green, Yellow, and Red
    #   ==============================================
    #
    } else if (heatmap.color=="GreenYellowRed"){

        RedRamp <- rgb(allOnes, one2zeor, allZeros);
        GreenRamp <- rgb(zero2one, allOnes, allZeros);
        ColorRamp <- cbind(GreenRamp, RedRamp);

    #   Green, Black, and Red
    #   ==============================================
    #
    } else if (heatmap.color=="GreenBlackRed"){

        RedRamp <- rgb( zero2one, allZeros, allZeros);
        GreenRamp <- rgb(allZeros, one2zeor, allZeros);
        ColorRamp <- cbind(GreenRamp, RedRamp);

    #   Yellow to Red
    #   ==============================================
    #
    } else if (heatmap.color=="YellowToRed") {
    
        ColorRamp <- rgb(allOnes, one2zeor, allZeros);

    #   black only
    #   ==============================================
    #
    } else {
        ColorRamp <- rgb(one2zeor, one2zeor, one2zeor);
    }

    return (ColorRamp);
}




#   ========================================================================
#
#   2.  RCircos.Plot.Heatmap.Color.Scale()
#
#   Plot a color scale for heatmap plot
#
#   Arguments:
#
#       maxValue:       Numeric, maximum value in heatmap scale
#       minValue:       Numeric, minimum value in heatmap scale
#       colorType:      Character vector for color specification, either 
#                       "BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed",
#                       "GreenBlackRed", "YellowToRed", or "BlackOnly".
#                       default is "BlueWhiteRed"
#       scaleLocation:  Integer representing color scale location, 1~4 for 
#                       bottom, left, top, right, bottomleft, bottomright,
#                       leftright, leftbottom, topleft, topright, righttop,
#                       and rightbottom.
#       scaleWidth:     Non-negative numeric, width of color scale. if not 
#                       defined, 1/2 of x axis will be used
#       scaleHeight:    Non-negative numeric, height of color scale. if not 
#                       defined, 1/10 of scaleWidth axis will be used
#
#   Returned value:     None
#
#   Example:    RCircos.Plot.Heatmap.Color.Scale(max.value=3, min.value=-3,
#                       scaleLocation=1, colorType="BlueWhiteRed", 
#                       scaleWidth=0, scaleHeight=0)
#

RCircos.Plot.Heatmap.Color.Scale <- function(max.value=NULL, min.value=NULL,
    color.type="BlueWhiteRed", scale.location=1, scale.width=0, scale.height=0) 
{
    #   Argument checking
    #   =====================================================
    #
    if(is.null(max.value) || is.null(min.value)) 
        stop("Max or Min value for color scale is missing.\n");
    if(!is.numeric(max.value) || !is.numeric(min.value) || 
        !is.numeric(scale.location) || !is.numeric(scale.width) || 
        !is.numeric(scale.height))
        stop("All arguments except of color.type must be numeric.\n")
    if(scale.location > 4 || scale.location < 1)
        stop("scale.location must be 1 ~ 4.\n");
    if(scale.width<0 || scale.height<0)
        stop("Incorrect scale.width or scale.height.");

    if(scale.width == 0 || scale.height == 0)
    {
        RCircos.Par <- RCircos.Get.Plot.Parameters();
        scale.width <- round(RCircos.Par$plot.radius/2);
        scale.height <- signif(scale.width/10, digits=2);
    }
    scaleCorr <- RCircos.Get.Heatmap.Color.Scale.Location(scale.location)
    coorX <- scaleCorr[1]; coorY <- scaleCorr[2];

    #   Standard color map
    #   ==================================================
    #
    colorRamp <- RCircos.Get.Heatmap.Color.Scale(color.type)
    totalRect <- nrow(colorRamp)*ncol(colorRamp)

    #   Plot horizontal color scale at bottom or top
    #   ===========================================
    #
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    if (scale.location %in% c(1, 3)) 
    {
        rectWidth <- scale.width/totalRect;
        rectHeight <- scale.height;

        if(scale.location == 1) {
            yTop <- -1*RCircos.Par$plot.radius; 
        } else { yTop <- RCircos.Par$plot.radius + rectHeight; }
        yBottom <- yTop - scale.height;
     
        for(aRect in seq_len(totalRect)) {
            xLeft <-  coorX + (aRect-1)*rectWidth;
            xRight <- xLeft + rectWidth;

            rect(xLeft, yBottom, xRight, yTop, col=colorRamp[aRect],  
                        border = NA);
        }

        text(coorX, yTop-(scale.height/2), min.value, pos=2);
        text(coorX+scale.width, yTop-(scale.height/2), max.value, pos=4);

    #   Plot vertical color scale at left or right
    #   ===========================================
    #
    } else {
        rectWidth  <- scale.height;
        rectHeight <- scale.width/totalRect;

        if(scale.location == 2) {
            xLeft <-  -1*RCircos.Par$plot.radius - rectWidth;
        } else {xLeft <-  RCircos.Par$plot.radius;}
        xRight <- xLeft + rectWidth;


        for(aRect in rev(seq_len(totalRect))) 
        {
            yTop <- coorY - (aRect-1)*rectHeight;
            yBottom <- yTop - rectHeight;
            rect(xLeft, yBottom, xRight, yTop, col=colorRamp[aRect],  
                border = NA);
        }

        text(xLeft+rectWidth/2, coorY, min.value, pos=3);
        text(xLeft+rectWidth/2, coorY-scale.width, max.value, pos=1);
    }
}




#   =========================================================================
#
#   3.  RCircos.Get.Heatmap.Data.Colors()
#
#   Calculate heatmap colors for one column of heatmap data. As the colour
#   map needs to be calculated from whole datasets when there are more than
#   one data column.
#
#   Argument:
#
#       heatmap.value:  A nemeric vector.
#       min.value:      Numeric, the minimum value for heatmap scale.
#       max.value:      Numeric, the maximum value for heatmap scale.
#
#   Return value:       Character vector of R color names with length same
#                       as the length of data values
#
#   Example:    data(RCircos.Heatmap.Data)
#               heatmap.data <- RCircos.Heatmap.Data;
#               hColors <- RCircos.Get.Heatmap.Data.Colors(heatmap.data, -3, 3);
#
#   Last modified on June 30, 2015
#

RCircos.Get.Heatmap.Data.Colors <- function(heatmap.value=NULL, 
                                    min.value=NULL, max.value=NULL) 
{
    if(is.null(heatmap.value) || is.null(min.value) || is.null(max.value))
        stop("Missing argument in RCircos.Get.Heatmap.Data.Colors().\n");

    if(!is.numeric(heatmap.value) || !is.numeric(min.value) || 
        !is.numeric(max.value))
        stop("Arguments must be numeric.\n");

    RCircos.Par <- RCircos.Get.Plot.Parameters();
    colorMap <- RCircos.Get.Heatmap.Color.Scale(RCircos.Par$heatmap.color);
    colorLevel <- seq(min.value, max.value, length=length(colorMap));

    cellColor <- rep(colorMap[1], length(heatmap.value));
    for(aRow in seq_len(length(heatmap.value)))
    {
        theLevel <- which(colorLevel>=heatmap.value[aRow]);
        cellColor[aRow] <- colorMap[min(theLevel)];
    }
    
    return (cellColor);
}




#   ===================================================================
#
#   4.  RCircos.Get.Heatmap.Color.Scale.Location()
#
#   Get heatmap color scale coordinates.
#
#   Argument:   scaleLocation, integer of 1:12, represents the plot 
#               location, default 1.
#               1:  "bottom"
#               2:  "left"
#               3:  "top"
#               4.: "right"
#               5:  "bottomleft"
#               6:  "bottomright"
#               7:  "leftright"
#               8:  "leftbottom"
#               9:  "topleft"
#               10: "topright"
#               11: "righttop"
#               12: "rightbottom"
#
#   Returned value: numeric verctor of length 2 for x and y coordinates
#
#   Example:    scale.location <- RCircos.Get.Heatmap.Color.Scale.Location(1)
#

RCircos.Get.Heatmap.Color.Scale.Location <- function(scale.location=1)
{
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    scaleWidth  <- RCircos.Par$plot.radius;
    scaleHeight <- scaleWidth/10;
    
    scaleX <- c(-0.5*scaleWidth, -1*scaleWidth-scaleHeight, 
                -0.5*scaleWidth, scaleWidth,
                0, -1*scaleWidth, 
                -1*scaleWidth-scaleHeight, -1*scaleWidth-scaleHeight,
                -1*scaleWidth, 0, 
                 scaleWidth, scaleWidth);
    scaleY <- c(-1*scaleWidth, 0.5*scaleWidth,
                scaleWidth + scaleHeight, 0.5*scaleWidth,
                -1*scaleWidth, -1*scaleWidth,
                0, scaleWidth, 
                scaleWidth + scaleHeight, scaleWidth + scaleHeight,
                scaleWidth, 0);
                
    coorX <- scaleX[scale.location];
    coorY <- scaleY[scale.location];
    
    return (c(coorX, coorY));
}



#   End of RCircosHeatmap.R
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>
