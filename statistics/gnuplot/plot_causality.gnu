#! /usr/bin/gnuplot/
set size square
set term pngcairo size 480, 480
set output "TE.png"

set colorsequence classic

data = "../np_data/info/TE.d"

set palette defined (0 "white", 1 "slategray", 2 "navy")

set cbtics 0.1

set xtics ("{/=12 {/Times-New-Roman:Italic u'}}" 0, \
           "{/=12 {/Times-New-Roman:Italic v'}}" 1, \
           "{/=12 {/Times-New-Roman:Italic w'}}" 2)
set ytics ("{/=12 {/Times-New-Roman:Italic u'}}" 0, \
           "{/=12 {/Times-New-Roman:Italic v'}}" 1, \
           "{/=12 {/Times-New-Roman:Italic w'}}" 2)

set cblabel "{/=18 {/Times-New-Roman:Italic Transfer entropy}}"
plot data matrix with image

