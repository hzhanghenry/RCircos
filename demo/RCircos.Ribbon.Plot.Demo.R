    #
    #   This demo draw ribbon links between chromosomes
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Tile.Plot.Demo");
    #
    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos library
    #   _________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);


    #   Load human cytoband data and link data
    #   _________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(RCircos.Ribbon.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


    #   Open the graphic device (here a pdf file)
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    pdf(file="RCircos.Ribbon.Plot.Demo.pdf", height=8, width=8);
    RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Ribbon Plot Demo");


    #   Link lines. Link data has only paired chromosome 
    #   locations in each row and link lines are always 
    #   drawn inside of chromosome ideogram.
    #   __________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    ribbon.data <- RCircos.Ribbon.Data;
    ribbon.data["PlotColor"] <- c("red", "blue", "green", "cyan");

    track.num <- 11;
    RCircos.Ribbon.Plot(ribbon.data, track.num, FALSE);


    #   Close the graphic device and clear memory
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("RCircos Ribbon Plot Demo Done!");





