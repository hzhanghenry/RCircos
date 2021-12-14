    #
    #   This demo draw chromosome ideogram with padding between
    #   chromosomes, highlights, chromosome names, and line plot. 
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Line.Plot.Demo");
    #   ______________________________________________________________________
    #   <RCircos Demo><RCircos Demo><RCircos Demo><RCircos Demo><RCircos Demo>



    #   Load RCircos library
    #   ___________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);


    #   Load human cytoband data and line plot data
    #   ____________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(RCircos.Line.Data);
    RCircos.Line.Data$chromosome <- paste0("chr", 
        RCircos.Line.Data$chromosome);

    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   _____________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


    #   Open the graphic device (here a pdf file)
    #   _____________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    out.file <- "RCircos.Line.Plot.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   _____________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Draw chromosome ideogram ...\n");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Line Plot Demo");


    #   Plot lines
    #   _____________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    line.data <- RCircos.Line.Data;
    line.colors <- rep("blue", nrow(line.data))
    rows <- which(abs(line.data$seg.mean)>=1.5);
    line.colors[rows] <- "red";
    line.data["PlotColor"] <- line.colors;

    data.col <- 5; 
    track.num <- 1;
    direction <- "in";
    RCircos.Line.Plot(line.data, data.col, track.num, "in");


    #   Close the graphic device and clear memory
    #   ___________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("RCircos Line Plot Demo Done!");




