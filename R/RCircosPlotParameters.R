#
#   Functions to initialize and reset RCircos plot parameters
#
#   1.  RCircos.Initialize.Plot.Parameters()
#   2.  RCircos.Validate.Plot.Parameters()
#   3.  RCircos.List.Plot.Parameters()
#   4.  RCircos.Reset.Plot.Parameters()
#
#   Last debug done on September 14, 2016
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>




#   ========================================================================
#
#   1.  RCircos.Initialize.Plot.Parameters()
#
#   Arguments:
#
#       tracks.inside:  non-negative  integer, number of plot tracks
#                       inside of chromosome ideogram
#       tracks.outside: non-negative  integer, number of plot tracks
#                       outside of chromosome ideogram
#   Return value: None
#
#   Example:    tracks.inside  <- 10 
#               tracks.outside <- 5
#               RCircos.Initialize.Parameters(tracks.inside, tracks.outside)
#
#   Last revised on June 11, 2015
#   

RCircos.Initialize.Plot.Parameters <- function(tracks.inside=NULL, 
                        tracks.outside=NULL)
{
    if(is.null(tracks.inside) || is.null(tracks.outside<0) )
    { stop("tracks.inside and tracks.outside must be defined.\n");}

    if(tracks.inside < 0 ||tracks.outside < 0) 
    { stop("tracks.inside and tracks.outside cannot be negative.\n");}

    radius.default        <- 1.0;
    tracks.inside.default <- 4;
    track.thick.default   <- 0.12;
    plot.radius.default   <- 1.8;

    #   Total of four data tracks are allowed inside of chromosome
    #   ideogram when radius.len is 1.0. Remains are for link lines. 
    #   If there are more than four tracks, increase the radius 
    #   length for extra plot area.
    #   ===========================================================
    #
    moreInside <- tracks.inside - tracks.inside.default;
    if(moreInside > 0)
    { 
        moreArea <- (moreInside + 1)* track.thick.default;
        radius.default <- radius.default + moreArea;
        plot.radius.default <- plot.radius.default + moreArea; 
    }

    #   If there will be data tracks outside of chromosome 
    #   ideogram, increase plot radius to get more room
    #   ==================================================
    #
    if(tracks.outside>0)
    { 
        moreArea <- (tracks.outside + 1)* track.thick.default;
        plot.radius.default <- plot.radius.default + moreArea; 
    }

    #   Set default plot parameters to a list. Color parameters 
    #   go first then character and numeric parameters
    #   =======================================================
    #
    plot.param <- list(
    
        #   part 1: these two parameters will affect plot 
        #   locations of cyto band and genomic data points
        #
        base.per.unit=RCircos.Get.Default.Base.Per.Units(),
        chrom.paddings=RCircos.Get.Default.Chrom.Padding(),

        #   part 2: These parameters are layout related
        #
        radius.len=radius.default,
        tracks.inside=tracks.inside,
        tracks.outside=tracks.outside,

        chr.ideo.pos=radius.default + 0.1,
        highlight.pos=radius.default + 0.25,
        chr.name.pos=radius.default  + 0.3,

        track.out.start=radius.default + 0.65, 
        track.in.start=radius.default + 0.05, 
        
        plot.radius=plot.radius.default,

        chrom.width=0.1,
        track.padding=0.02, 
        track.height=0.1, 
 
        #   part 3: These parameters can be reset freely 
        #   with no need to reset other parameters
        #
        hist.width=100,
        heatmap.width=100,
        text.size=0.4,
        char.width=500,
        highlight.width=round(radius.default, digits=0),
        point.size=1,

        Bezier.point=1000,
        max.layers=5,
        sub.tracks=5,

        text.color="black",
        hist.color="red",
        line.color="black",
        scatter.color="black",
        tile.color="black",
        track.background="wheat",
        grid.line.color="gray",

        heatmap.color="BlueWhiteRed",
        point.type="." 
    );

    #   Put the plot parameter in RCircos environment
    #   ==============================================
    #
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir=globalenv());
    RCircosEnvironment[["RCircos.PlotPar"]] <- plot.param;
}




#   ========================================================================
#
#   2.  RCircos.Validate.Plot.Parameters()
#
#   All numeric parameters must be positive and colors must be supported 
#   by R and RCircos. No ajustment for numeric values here.
#
#   Argument:       new.params, a list holding all plot parameters
#   Return values:  None
#
#   Example:    new.params <- RCircos.Get.Plot.Parameters()
#               new.params$base.per.unit <- 1500
#               RCircos.Validate.Plot.Parameters(new.params) 
#
RCircos.Validate.Plot.Parameters <- function(parameters=NULL)
{
    if(is.null(parameters)) stop("Missing function argument.\n");

    #   No negative numeric parameters allowed
    #   ==============================================
    #
    for(aParam in seq_along(parameters))
    {
        if(class(parameters[[aParam]]) == "numeric" &&
            parameters[[aParam]] < 0) {
            stop("Plot parameters cannot have negative values.");
        }
    }

    #   Check out heatmap colors
    #   ==============================================
    #
    heatmapColors <- RCircos.Get.Supported.HeatmapColors();
    if((parameters$heatmap.color %in% heatmapColors) == FALSE)
    {   stop("\nUnsupported heatmap color defined.\nType 
            RCircos.Get.Supported.HeatmapColors() for more info.");
    }

    #   Check out other plot colors
    #   ==============================================
    #
    colorNames <- colors();
    if( !parameters$text.color %in% colorNames ||
        !parameters$hist.color  %in% colorNames ||
        !parameters$line.color  %in% colorNames ||
        !parameters$scatter.color  %in% colorNames ||
        !parameters$tile.color  %in% colorNames  )
    { stop("Unsupported R plot color defined.") }

    #   track background and grid line can be NULL or NA
    #   ================================================
    background <- parameters$track.background;
    if(!is.null(background) && !is.na(background)){
        if(!background  %in% colorNames)
            stop("Unsupported color for track background.")
    }

    grid.color <- parameters$grid.line.color;
    if(!is.null(grid.color) && !is.na(grid.color)){
        if(!grid.color %in% colorNames)
            stop("Unsupported color for grid line.")
    }
}




#   ========================================================================
#
#   3.  RCircos.Reset.Plot.Parameters()
#
#   Argument:       new.params, a list holding all plot parameters
#   Return values:  None
#
#   Example:    new.params <- RCircos.Get.Plot.Parameters()
#               new.params$base.per.unit <- 1500
#               RCircos.Reset.Plot.Parameters(new.params) 
#

RCircos.Reset.Plot.Parameters <- function (new.params=NULL) 
{
    if(is.null(new.params)) stop("Missing function argument.\n");
    old.params <- RCircos.Get.Plot.Parameters();

    #   1.  If parameters related to total number of data tracks 
    #       need reset, use reset core components instead.
    #   ==========================================================
    if( new.params$radius.len != old.params$radius.len ||
        new.params$plot.radius != old.params$plot.radius ||
        new.params$chr.ideo.pos != old.params$chr.ideo.pos ||
        new.params$tracks.inside != old.params$tracks.inside ||
        new.params$tracks.outside != old.params$tracks.outside )
    { stop("Please use RCircos.Set.Core.Components() instead.\n") }

    #   2.  If parameters related to chromosome ideogram plot 
    #       need reset, use customized plot methods instead.
    #   ==========================================================
    if( new.params$chr.ideo.pos != old.params$chr.ideo.pos ||
        new.params$highlight.pos != old.params$highlight.pos ||
        new.params$chr.name.pos  != old.params$chr.name.pos )
    { stop("Please use customized ideogram plot methods instead.\n")}

    #   3.  Validate the parameter values in case of multiple
    #       parameters were reset in new.params such as nemeric
    #       values and color values, and ideogram layout values.
    #   ==========================================================
    RCircos.Validate.Plot.Parameters(new.params);

    #   4.  Parameters related to ideogram width change.  
    #       Note: chr.ideo.pos is a read-only parameter
    #   ========================================================
    if( new.params$chrom.width != old.params$chrom.width )
    {
        differ <- new.params$chrom.width - old.params$chrom.width;
        new.params$highlight.pos <- old.params$highlight.pos + differ;

        new.name.pos <- old.params$chr.name.pos + differ;
        if(new.params$chr.name.pos < new.name.pos)
            new.params$chr.name.pos  <- new.name.pos;
    }
    
    #   5.  In case user midified track.in.start and track.out.start
    #   ===========================================================
    if(new.params$track.in.start >= new.params$chr.ideo.pos)
        new.params$track.in.start <- new.params$chr.ideo.pos - 0.05;

    new.name.end <- new.params$chr.name.pos + 0.3;
    if(new.params$track.out.start < new.name.end)
        new.params$track.out.start <- new.name.end;

    #   6.  Parameters related to data track layout. If reset, total
    #       tracks will be different. Just validate and give a prompt
    #   ============================================================
    if(new.params$track.padding !=  old.params$track.padding ||
            new.params$track.height != old.params$track.height )
    { 
        message(paste0("Track height and/or track padding have been ",
            "reset\n. Actual total data track plotted may differ.\n"));
    }

    #   7.  Adjust chromosome padding parameter with default constant
    #       if base.per.unit was rest but chrom.padding was unchanged
    #       or the new chrom.paddings is too big
    #   =============================================================
    if(old.params$base.per.unit != new.params$base.per.unit &&
        old.params$chrom.paddings == new.params$chrom.paddings )
    {
        RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
        band.len <- RCircos.Cyto$ChromEnd - RCircos.Cyto$ChromStart;
        genome.len <- sum(as.numeric(band.len));
        padding.const <- RCircos.Get.Padding.Constant();

        total.units <- genome.len/new.params$base.per.unit;
        new.padding <- round(padding.const*total.units, digits=0);

        if(new.padding != new.params$base.per.unit) {
            message(paste("\nNote: chrom.padding", 
                new.params$chrom.paddings,
                " was reset to", new.padding, "\n"));
            new.params$chrom.paddings <- new.padding;
        }
    }

    #   Save new parameters to RCircos Environment
    #   =====================================================
    RCircosEnvironment <- NULL;
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv());
    RCircosEnvironment[["RCircos.PlotPar"]] <- NULL;
    RCircosEnvironment[["RCircos.PlotPar"]] <- new.params;

    #   Ideogram/band positions are binded to base.per.unit and
    #   chromosome padding so have to be reset if base.per.unit 
    #   and/or chrom.paddings are reset.
    #   ======================================================
    if(old.params$base.per.unit != new.params$base.per.unit ||
        old.params$chrom.paddings != new.params$chrom.paddings) 
    {
        RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
        RCircos.Cyto <- RCircos.Cyto[,1:5];

        RCircosEnvironment[["RCircos.Cytoband"]] <- NULL
        RCircos.Set.Cytoband.Data(RCircos.Cyto);

        RCircosEnvironment[["RCircos.Base.Position"]] <- NULL
        RCircos.Set.Base.Plot.Positions();
    }
}






#   ========================================================================
#
#   4.  RCircos.List.Plot.Parameters()
#
#   Print out all parameters. This could be ran any time to checking out the 
#   current values of RCircos plot parameters
#
#   Argument:       None
#   Return value:   None
#
#   Example:        RCircos.List.Parameters();
#

RCircos.List.Plot.Parameters <- function()
{
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    message("Parameters for current RCircos session.\n\n",

    "Parameters in inch:\n",
    "==============================\n",
    "radius.len:\t\t",     RCircos.Par$radius.len, "\n",
    "chr.ideo.pos:\t\t",  RCircos.Par$chr.ideo.pos, "\n",
    "highlight.pos:\t\t",  RCircos.Par$highlight.pos, "\n",
    "chr.name.pos:\t\t",   RCircos.Par$chr.name.pos, "\n",
    "plot.radius:\t\t",    RCircos.Par$plot.radius, "\n",
    "track.in.start:\t\t", RCircos.Par$track.in.start, "\n",
    "track.out.start:\t",  RCircos.Par$track.out.start, "\n",
    "chrom.width:\t\t",    RCircos.Par$chrom.width, "\n",
    "track.padding:\t\t",  RCircos.Par$track.padding, "\n",
    "track.height:\t\t",   RCircos.Par$track.height, "\n\n",

    "Parameters in chromosome unit:\n",
    "==============================\n",
    "base.per.unit:\t\t",  RCircos.Par$base.per.unit, "\n",
    "chrom.paddings:\t\t", RCircos.Par$chrom.paddings, "\n",
    "heatmap.width:\t\t",  RCircos.Par$heatmap.width, "\n",
    "hist.width:\t\t",     RCircos.Par$hist.width, "\n",
    "gene name char. width:\t", RCircos.Par$char.width, "\n\n",
   
    "General R graphic parameters:\n",
    "==============================\n",
    "text.size:\t\t",      RCircos.Par$text.size, "\n",
    "highlight.width:\t",  RCircos.Par$highlight.width, "\n",
    "point.type:\t\t",     RCircos.Par$point.type, "\n",
    "point.size:\t\t",     RCircos.Par$point.size, "\n",

    "text.color:\t\t",    RCircos.Par$text.color, "\n",
    "heatmap.color:\t\t", RCircos.Par$heatmap.color, "\n",
    "hist.color:\t\t",    RCircos.Par$hist.color, "\n",
    "line.color:\t\t",    RCircos.Par$line.color, "\n",
    "scatter.color:\t\t", RCircos.Par$scatter.color, "\n",
    "tile.color:\t\t",    RCircos.Par$tile.color, "\n",

    "track.background:\t",  RCircos.Par$track.background, "\n",
    "grid.line.color:\t",   RCircos.Par$grid.line.color, "\n",

    "Bezier.point:\t\t",  RCircos.Par$Bezier.point, "\n",
    "max.layers:\t\t",    RCircos.Par$max.layers, "\n",
    "sub.tracks:\t\t",    RCircos.Par$sub.tracks, "\n\n",
    
    "Data track numbers:\n",
    "==============================\n",
    "tracks.inside:\t\t", RCircos.Par$tracks.inside, "\n",
    "tracks.outside:\t\t", RCircos.Par$tracks.outside, "\n\n"
    );


    #   User friendly notice
    #   ===========================================================
    #
    message("Following are procedures to change RCircos plot parameters:\n",
        "params <- RCircos.Get.Plot.Parameters();\n",
        "params$radius.len <- 2.0;\n",
        "params$base.per.unit <- 5000;\n",
        "RCircos.Reset.Plot.Parameters(params)\n\n",
        "Chromosome ideogram data were automatically modified.\n\n");
}


#   End of RCircosPlotParameters.R
#   ________________________________________________________________________
#   <RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos><RCircos>