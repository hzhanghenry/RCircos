    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>
    #
    #   This demo draw chromosome ideogram with padding between chromosomes, 
    #   highlights, chromosome names, and histogram. 
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Polygon.Demo");
    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos package and defined parameters
    #   ________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);


    #   Load human cytoband data and gene expression data
    #   ________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(RCircos.Polygon.Data);
    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info, NULL, 10, 10);


    #   Open the graphic device (here a pdf file)
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    pdf(file="RCircos.Polygon.Demo.pdf", height=8, width=8);
    RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Draw chromosome ideogram ...\n");

    RCircos.Chromosome.Ideogram.Plot();
    title("RCircos Polygon Plot Demo");


    #   Plot histogram Inside of chromosome ideogram
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    polygon.data <- RCircos.Polygon.Data;
    RCircos.Polygon.Plot(polygon.data, track.num=NULL, data.col=4,
        genomic.columns=3, side="in", inside.pos=1.6, outside.pos=1.8, 
        border.col="black", polygon.col="red", is.sorted=TRUE)

    polygon.data$Data <- polygon.data$Data*-1
    RCircos.Polygon.Plot(polygon.data, track.num=NULL, data.col=4,
        genomic.columns=3, side="in", inside.pos=1.3, outside.pos=1.5, 
        border.col="black", polygon.col="green", is.sorted=TRUE)

    rows <- seq(1, 77, by=2)

    plot.colors <- rep("red", nrow(polygon.data))
    plot.colors[rows] <- "green"
    polygon.data["PlotColor"] <- plot.colors
    polygon.data$Data[rows] <- polygon.data$Data[rows]*-1

    RCircos.Polygon.Plot(polygon.data, track.num=NULL, data.col=4,
        genomic.columns=3, side="in", inside.pos=0.8, outside.pos=1.2, 
        border.col="black", polygon.col=NULL, is.sorted=TRUE)

    title("RCircos Polygon Plot Demo")



    #   Close the graphic device and clear memory
    #   _________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("RCircos Polygon Plot Demo Done!");

