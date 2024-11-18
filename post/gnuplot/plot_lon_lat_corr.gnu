#! /usr/bin/gnuplot/
set format x "%.0f"

set key off
set xlabel "{/=18 {/Times-New-Roman:Italic r}}"
set xrange [0:300]

set tics font "Times-New-Roman, 16"
set xtics nomirror
set mxtics 2

corr = "vel_corr.d"
TE   = "TE_lon_lat.d"

set output "lon_lat_corr.png"
set term pngcairo size 480, 960

set multiplot layout 2, 1

set ylabel "{/=18 {/Times-New-Roman:Italic f}}"
set format y "%.2f"
set yrange [0.6:1]
set ytics nomirror
set mytics 1
plot corr using 1:2 with lines lw 2 lc rgb "blue" 

set ylabel "{/=18 {/Times-New-Roman:Italic g}}"
set format y "%.1f"
set yrange [-0.2:1]
set ytics nomirror
set mytics 1
plot corr using 1:3 with lines lw 2 lc rgb "blue" 

unset multiplot

set output "TE.png"
set term pngcairo size 480, 960

set multiplot layout 2, 1

set ylabel "{/=18 {/Times-New-Roman:Italic TE}}"
set format y "%.2f"
set yrange [0:0.3]
set ytics nomirror
set mytics 1
plot TE using 1:2 with lines lw 2 lc rgb "red", \
     TE using 1:3 with lines lw 2 lc rgb "blue"

set ylabel "{/=18 {/Times-New-Roman:Italic TE}}"
set format y "%.1f"
set yrange [0:0.7]
set ytics nomirror
set mytics 1
plot TE using 1:4 with lines lw 2 lc rgb "red", \
     TE using 1:5 with lines lw 2 lc rgb "blue"

unset multiplot


