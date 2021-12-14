    #
    #   This demo draw chromosome ideogram with padding between  
    #   chromosomes, highlights, chromosome names, and three 
    #   empty tracks inside and outside of chromosome ideogram.
    #
    #   Usage:
    #
    #   library(RCircos);
    #   demo("RCircos.Layout.Demo");
    #   ______________________________________________________________________
    #   <RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO>



    #   Load RCircos library
    #   _______________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    library(RCircos);


    #   Load human cytoband data 
    #   _______________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data(UCSC.HG19.Human.CytoBandIdeogram);
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


    #   Setup RCircos core components:
    #   _______________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Set.Core.Components(cyto.info, NULL, 5, 5);



    #   Open the graphic device (here a pdf file)
    #
    #   png(file="RCircos.Layout.Demo.png", height=8, 
    #      width=8, unit="in", type="cairo", res=300);
    #   ________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    out.file <- "RCircos.Layout.Demo.pdf";
    pdf(file=out.file, height=8, width=8);

    RCircos.Set.Plot.Area();


    #   Draw chromosome ideogram
    #   _________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    message("Draw chromosome ideogram ...\n");

    RCircos.Chromosome.Ideogram.Plot();
    title("R Circos Layout Demo");



    #   Marking plot areas both inside and outside of
    #   chromosome ideogram ( 3 for each)
    #   ********************************************
    total.track <- 3;
    subtrack <- 5;


    #   Outside of chromosome ideogram
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    RCircos.Par <- RCircos.Get.Plot.Parameters();
    one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;

    for(a.track in 1:total.track)
    {   
        in.pos    <- RCircos.Par$track.out.start  + (a.track-1)*one.track;
        out.pos <- in.pos + RCircos.Par$track.height;
        RCircos.Track.Outline(out.pos, in.pos, subtrack);
    }

    #   Inside of chromosome ideogram
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    for(a.track in 1:total.track)
    {
        out.pos <- RCircos.Par$track.in.start  - (a.track-1)*one.track;
        in.pos    <- out.pos - RCircos.Par$track.height;
        RCircos.Track.Outline(out.pos, in.pos, subtrack);
    }


    #   Close the graphic device and clear memory
    #   _________________________________________________________________
    #   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    dev.off();
    message("RCircos Layout Demo Done!");




