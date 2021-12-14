    #
    #   This demo draw chromosome ideogram with padding between  
    #   chromosomes, highlights, chromosome names, and scatters. 
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Scatter.Plot.Demo");
    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos library
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);


    #   Load human cytoband data and scatterplot data
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(RCircos.Scatter.Data);
	if(length(grep("chr", RCircos.Scatter.Data$chromosome)) == 0)
		RCircos.Scatter.Data$chromosome <- paste0("chr", 
        RCircos.Scatter.Data$chromosome);
    
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


    #   Open the graphic device (here a pdf file)
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    pdf(file="RCircos.Scatter.Plot.Demo.pdf", height=8, width=8);
    RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Scatter Plot Demo");


    #   Scatterplot 
    #   _________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    scatter.data <- RCircos.Scatter.Data;
    scatter.colors <- rep("cyan", nrow(scatter.data));
    scatter.colors[which(scatter.data$seg.mean>=2)] <- "red";
    scatter.colors[which(scatter.data$seg.mean<=-2)] <- "blue";
    scatter.data["PlotColor"] <- scatter.colors;

    data.col <- 5;
    track.num <- 6; 
    RCircos.Scatter.Plot(scatter.data, data.col, track.num, "in", 1);



    #   Close the graphic device
    #   ___________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("RCircos Scatter Plot Demo Done!");

