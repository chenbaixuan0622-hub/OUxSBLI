#! /usr/bin/gnuplot/
set key off
set size square

set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 5

data = "wall_prop.d"

delta = 2.e-3

set output "wall_prop.png"
set term pngcairo size 960, 480

set multiplot layout 1, 2

set xlabel "{/=18 {/Times-New-Roman:Italic x /{/Symbol:Italic d}}}"
set format x "%.1f"
set xtics 5
set mxtics 2
set xrange [5:20]

set format y "%.1f"
set ytics  1
set mytics 2
set yrange [-1:3]
plot data using ($1/delta):($2*1e3) with lines lw 2 lc rgb "blue"

set format y "%.1f"
set ytics 0.5
set mytics 2
set yrange [0.8:2]
plot data using ($1/delta):3 with lines lw 2 lc rgb "blue"

unset multiplot

