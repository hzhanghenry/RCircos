    #
    #   This demo draw chromosome ideogram with padding between 
    #   chromosomes, highlights, chromosome names, and tile plot. 
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Tile.Plot.Demo");
    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos library
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);


    #   Load human cytoband data and gene expression data
    #   _________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(RCircos.Tile.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


    #   Open the graphic device (here a pdf file)
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    pdf(file="RCircos.Tile.Plot.Demo.pdf", height=8, width=8);
    RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   _________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Tile Plot Demo");


    #   Tile plot. Note: tile plot data have chromosome 
    #   locations only and each data file is for one track
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    tile.data <- RCircos.Tile.Data;
    tile.colors <- rainbow(nrow(tile.data));
    tile.data["PlotColor"] <- tile.colors;

    track.num <- 9;
    RCircos.Tile.Plot(tile.data, track.num, "in");


    #   Close the graphic device and clear memory
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("RCircos Tile Plot Demo Done!");

