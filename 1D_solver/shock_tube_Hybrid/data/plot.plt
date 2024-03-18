set term pngcairo size 1280, 480
set output "output.png"

set format xy "%4.2f"
set key off

set xlabel "{/=12 {/Times-New-Roman:Italic x} (m)}"

set xrange [0:1]

set tics font "Times-New-Roman,10"
set xtics nomirror
set ytics nomirror

datafile = "Q00010.d"

set multiplot layout 1, 4

set ylabel "{/=12 {/Symbol:Italic r} (kg / m^3)}"
set size square
set mxtics 10
set mytics 10
# stats datafile using 2 nooutput
# set yrange [0 : STATS_max * 1.1]
plot datafile using 1:2 with points pt 2 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic u} (m/s)}"
set size square
set mxtics 10
set mytics 10
# stats datafile using 3 nooutput
# set yrange [0 : STATS_max * 1.1]
plot datafile using 1:3 with points pt 2 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic P} (Pa)}"
set size square
set mxtics 10
set mytics 10
# stats datafile using 4 nooutput
# set yrange [0 : STATS_max * 1.1]
plot datafile using 1:4 with points pt 2 lc rgb "blue"

set ylabel "{/=12 {/Symbol:Italic f}}"
set size square
set mxtics 10
set mytics 10
# stats datafile using 5 nooutput
# set yrange [0 : STATS_max * 1.1]
plot datafile using 1:5 with points pt 2 lc rgb "blue"

unset multiplot

set term x11
replot

