set term pngcairo size 1280, 960
set output "output.png"

set format xy "%4.1f"
set key off

set xlabel "{/=12 {/Times-New-Roman:Italic x}}"

set xrange [0:1]

set tics font "Times-New-Roman,10"
set xtics nomirror
set ytics nomirror

datafile = "Q00100.d"

set multiplot layout 2, 3

set ylabel "{/=12 {/Symbol:Italic r}}"
set size square
set mxtics 5
set mytics 5
plot datafile using 1:2 with points pt 1 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic u}}"
set size square
set mxtics 5
set mytics 5
plot datafile using 1:3 with points pt 1 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic P}}"
set size square
set mxtics 5
set mytics 5
plot datafile using 1:4 with points pt 1 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic T}}"
set size square
set mxtics 5
set mytics 5
plot datafile using 1:5 with points pt 1 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic M}}"
set size square
set mxtics 5
set mytics 5
plot datafile using 1:6 with points pt 1 lc rgb "blue"

set ylabel "{/=12 {/Symbol:Italic m}}"
set size square
set mxtics 5
set mytics 5
plot datafile using 1:7 with points pt 1 lc rgb "blue"

unset multiplot

set term x11
replot

