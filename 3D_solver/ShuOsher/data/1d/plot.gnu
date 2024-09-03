#! /usr/bin/gnuplot/
set term pngcairo size 480, 480
set output "ShuOsher.png"

set format xy "%4.1f"
set key off

set xlabel "{/=12 {/Times-New-Roman:Italic x}}"

set xrange [0:1]

set tics font "Times-New-Roman,10"
set xtics nomirror
set ytics nomirror

data0 = "Q00000.d"
data1 = "Q00100.d"

set ylabel "{/=12 {/Symbol:Italic r}}"
set yrange [0.0:5.0]
set size square
set mxtics 5
set mytics 5
plot data0 using 1:2 with lines lw 2 lc rgb "black", \
     data1 using 1:2 with lines lw 2 lc rgb "red"

