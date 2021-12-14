    #
    #   This demo draw human chromosome ideogram and data tracks for:
    #
    #   1.  Connectors
    #   2.  Gene lables
    #   3.  Heatmap
    #   4.  Scatterplot 
    #   5.  Line plot
    #   6.  Histogram
    #   7.  Tile plot
    #   8.  Link lines
    #
    # ______________________________________________________________________
    # <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos package
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);



    #   Load human cytoband data 
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(UCSC.HG19.Human.CytoBandIdeogram);
    hg19.cyto <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info=hg19.cyto, chr.exclude=NULL, 
                tracks.inside=10, tracks.outside=0);



    #   Open the graphic device (here a pdf file)
    #
    #   png(file="RCircos.Demo.Human.png", height=8, width=8, unit="in", 
    #       type="cairo", res=300);
    #
    #   tiff(file="RCircos.Demo.Human.tif", height=8, width=8, unit="in", 
    #           type="cairo", res=300);
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Open graphic device and start plot ...\n");
    pdf(file="RCircos.Demo.Human.pdf", height=8, width=8);

    RCircos.Set.Plot.Area();
    title("RCircos 2D Track Plot with Human Genome");



    #   Draw chromosome ideogram
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Draw chromosome ideogram ...\n");
    RCircos.Chromosome.Ideogram.Plot();



    #   Connectors in first track and gene names in the second track. 
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add Gene and connector tracks ...\n");
    data(RCircos.Gene.Label.Data);

    RCircos.Gene.Connector.Plot(genomic.data=RCircos.Gene.Label.Data, 
                track.num=1, side="in");
    RCircos.Gene.Name.Plot(gene.data=RCircos.Gene.Label.Data, name.col=4, 
                track.num=2, side="in");


    #   Heatmap plot.  Since some gene names plotted above are longer 
    #   than one track height, we skip two tracks 
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add heatmap track ...\n");

    data(RCircos.Heatmap.Data);
    RCircos.Heatmap.Plot(heatmap.data=RCircos.Heatmap.Data, data.col=6, 
                track.num=5, side="in");



    #   Scatterplot. 
    #   Note:   There is no prefix of "chr" in chromosome names 
    #           of RCircos.Scatter.Data 
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add scatterplot track ...\n");
    data(RCircos.Scatter.Data);
    RCircos.Scatter.Plot(scatter.data=RCircos.Scatter.Data, data.col=5, 
                track.num=6, side="in", by.fold=1);



    #   Line plot. 
    #   Note:   There is no prefix of "chr" in chromosome names 
    #           of RCircos.Line.Data
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add line plot track ...\n");
    data(RCircos.Line.Data);
    RCircos.Line.Data[,1] <- paste0("chr", RCircos.Line.Data[,1])
    RCircos.Line.Plot(line.data=RCircos.Line.Data, data.col=5, 
                track.num=7, side="in");



    #   Histogram plot
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add histogram track ...\n");

    data(RCircos.Histogram.Data);
    RCircos.Histogram.Plot(hist.data=RCircos.Histogram.Data, data.col=4, 
                track.num=8, side="in");



    #   Tile plot. Note: tile plot data have chromosome locations and each
    #   data file is for one track
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add tile track ...\n");

    data(RCircos.Tile.Data);
    RCircos.Tile.Plot(tile.data=RCircos.Tile.Data, track.num=9, side="in");



    #   Link lines. Link data has only paired chromosome locations in
    #   each row and link lines are always drawn inside of chromosome 
    #   ideogram.
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Add link track ...\n");

    data(RCircos.Link.Data);
    RCircos.Link.Plot(link.data=RCircos.Link.Data, track.num=11, 
                by.chromosome=FALSE);



    #   Add ribbon link to the center of plot area (link lines).
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   message("Add ribbons to link track ...\n");

    data(RCircos.Ribbon.Data);
    RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, track.num=11, 
                by.chromosome=FALSE, twist=FALSE);





    #   Close the graphic device and clear memory
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("R Circos Demo Done ...\n\n");














