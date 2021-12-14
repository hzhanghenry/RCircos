\name{UCSC.HG19.Human.CytoBandIdeogram}
\alias{UCSC.HG19.Human.CytoBandIdeogram}
\docType{data}
\title{
Cytoband Information for Human Chromosome Ideogram 
}
\description{
A data frame containing chromosome name, start and end position, band name, and intensity of Giemsa stains for each cytoband. 
}
\usage{data(UCSC.HG19.Human.CytoBandIdeogram)}
\format{
  A data frame with 862 observations on the following 5 variables.
  \describe{
    \item{\code{Chromosome}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
    \item{\code{ChromStart}}{a numeric vector}
    \item{\code{ChromEnd}}{a numeric vector}
    \item{\code{Band}}{a factor with levels \code{p11} \code{p11.1} \code{p11.11} \code{p11.12} \code{p11.2} \code{p11.21} \code{p11.22} \code{p11.23} \code{p11.3} \code{p11.31} \code{p11.32} \code{p11.4} \code{p12} \code{p12.1} \code{p12.2} \code{p12.3} \code{p12.31} \code{p12.32} \code{p12.33} \code{p13} \code{p13.1} \code{p13.11} \code{p13.12} \code{p13.13} \code{p13.2} \code{p13.3} \code{p13.31} \code{p13.32} \code{p13.33} \code{p14} \code{p14.1} \code{p14.2} \code{p14.3} \code{p15} \code{p15.1} \code{p15.2} \code{p15.3} \code{p15.31} \code{p15.32} \code{p15.33} \code{p15.4} \code{p15.5} \code{p16.1} \code{p16.2} \code{p16.3} \code{p21} \code{p21.1} \code{p21.2} \code{p21.3} \code{p21.31} \code{p21.32} \code{p21.33} \code{p22} \code{p22.1} \code{p22.11} \code{p22.12} \code{p22.13} \code{p22.2} \code{p22.3} \code{p22.31} \code{p22.32} \code{p22.33} \code{p23} \code{p23.1} \code{p23.2} \code{p23.3} \code{p24.1} \code{p24.2} \code{p24.3} \code{p25.1} \code{p25.2} \code{p25.3} \code{p26.1} \code{p26.2} \code{p26.3} \code{p31.1} \code{p31.2} \code{p31.3} \code{p32.1} \code{p32.2} \code{p32.3} \code{p33} \code{p34.1} \code{p34.2} \code{p34.3} \code{p35.1} \code{p35.2} \code{p35.3} \code{p36.11} \code{p36.12} \code{p36.13} \code{p36.21} \code{p36.22} \code{p36.23} \code{p36.31} \code{p36.32} \code{p36.33} \code{q11} \code{q11.1} \code{q11.2} \code{q11.21} \code{q11.22} \code{q11.221} \code{q11.222} \code{q11.223} \code{q11.23} \code{q12} \code{q12.1} \code{q12.11} \code{q12.12} \code{q12.13} \code{q12.2} \code{q12.3} \code{q13} \code{q13.1} \code{q13.11} \code{q13.12} \code{q13.13} \code{q13.2} \code{q13.3} \code{q13.31} \code{q13.32} \code{q13.33} \code{q13.4} \code{q13.41} \code{q13.42} \code{q13.43} \code{q13.5} \code{q14} \code{q14.1} \code{q14.11} \code{q14.12} \code{q14.13} \code{q14.2} \code{q14.3} \code{q15} \code{q15.1} \code{q15.2} \code{q15.3} \code{q16.1} \code{q16.2} \code{q16.3} \code{q21} \code{q21.1} \code{q21.11} \code{q21.12} \code{q21.13} \code{q21.2} \code{q21.21} \code{q21.22} \code{q21.23} \code{q21.3} \code{q21.31} \code{q21.32} \code{q21.33} \code{q22} \code{q22.1} \code{q22.11} \code{q22.12} \code{q22.13} \code{q22.2} \code{q22.3} \code{q22.31} \code{q22.32} \code{q22.33} \code{q23} \code{q23.1} \code{q23.2} \code{q23.3} \code{q23.31} \code{q23.32} \code{q23.33} \code{q24} \code{q24.1} \code{q24.11} \code{q24.12} \code{q24.13} \code{q24.2} \code{q24.21} \code{q24.22} \code{q24.23} \code{q24.3} \code{q24.31} \code{q24.32} \code{q24.33} \code{q25} \code{q25.1} \code{q25.2} \code{q25.3} \code{q25.31} \code{q25.32} \code{q25.33} \code{q26} \code{q26.1} \code{q26.11} \code{q26.12} \code{q26.13} \code{q26.2} \code{q26.3} \code{q26.31} \code{q26.32} \code{q26.33} \code{q27} \code{q27.1} \code{q27.2} \code{q27.3} \code{q28} \code{q28.1} \code{q28.2} \code{q28.3} \code{q29} \code{q31.1} \code{q31.2} \code{q31.21} \code{q31.22} \code{q31.23} \code{q31.3} \code{q31.31} \code{q31.32} \code{q31.33} \code{q32} \code{q32.1} \code{q32.11} \code{q32.12} \code{q32.13} \code{q32.2} \code{q32.3} \code{q32.31} \code{q32.32} \code{q32.33} \code{q33} \code{q33.1} \code{q33.2} \code{q33.3} \code{q34} \code{q34.1} \code{q34.11} \code{q34.12} \code{q34.13} \code{q34.2} \code{q34.3} \code{q35} \code{q35.1} \code{q35.2} \code{q35.3} \code{q36.1} \code{q36.2} \code{q36.3} \code{q37.1} \code{q37.2} \code{q37.3} \code{q41} \code{q42.11} \code{q42.12} \code{q42.13} \code{q42.2} \code{q42.3} \code{q43} \code{q44}}
    \item{\code{Stain}}{a factor with levels \code{acen} \code{gneg} \code{gpos100} \code{gpos25} \code{gpos50} \code{gpos75} \code{gvar} \code{stalk}}
  }
}
\source{
http://genome.ucsc.edu 
}
\keyword{datasets}

