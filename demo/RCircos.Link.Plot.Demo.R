    #
    #   This demo draw chromosome ideogram with padding between chromosomes, 
    #   highlights, chromosome names, and link lines. 
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Link.Plot.Demo");
    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO> 



    #   Load RCircos library
    #   _______________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);


    #   Load human cytoband data and link data
    #   _______________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(RCircos.Link.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


    #   Open the graphic device (here a pdf file)
    #   _____________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    out.file <- "RCircos.Link.Plot.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   _____________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Draw chromosome ideogram ...\n");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Link Plot Demo");


    #   Link lines. Link data has paired chromosome 
    #   locations only in each row and link lines 
    #   are always drawn inside chromosome ideogram.
    #   ____________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add link track ...\n");

    link.data <- RCircos.Link.Data;
    link.colors <- rep("blue", nrow(link.data));
    rows <- seq(1, nrow(link.data), by=5);
    link.colors[rows] <- "red";
    link.data["PlotColor"] <- link.colors;

    track.num <- 2;
    RCircos.Link.Plot(link.data, track.num, FALSE);


    #   Close the graphic device and clear memory
    #   ___________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("RCircos Link Plot Demo Done!");


