#! /usr/bin/gnuplot/
set tics font "Times-New-Roman, 16"
set key off
set xtics nomirror
set ytics nomirror
set mxtics 5
set mytics 5
set logscale xy
set size square

PSD_t_1 = "PSD_time_1delta.d"
PSD_t_2 = "PSD_time_05delta.d"

set output "PSD_time.png"
set term pngcairo size 960, 480

set format y "10^{%L}"
set format x "10^{%L}"

set xlabel "{/=18 {/Times-New-Roman:Italic k} {/Times-New-Roman: [Hz]}}"
set ylabel "{/=18 {/Times-New-Roman:Italic PSD} {/Times-New-Roman: [m^2/s]}}"

set xrange [1e3:1e6]

plot PSD_t_1 using 1:2 with lines dt 1 lw 2 lc rgb "black", \
     PSD_t_1 using 1:3 with lines dt 1 lw 2 lc rgb "blue", \
     PSD_t_1 using 1:4 with lines dt 1 lw 2 lc rgb "red", \
     PSD_t_1 using 1:5 with lines dt 1 lw 2 lc rgb "green", \
     PSD_t_2 using 1:2 with lines dt 2 lw 2 lc rgb "black", \
     PSD_t_2 using 1:3 with lines dt 2 lw 2 lc rgb "blue", \
     PSD_t_2 using 1:4 with lines dt 2 lw 2 lc rgb "red", \
     PSD_t_2 using 1:5 with lines dt 2 lw 2 lc rgb "green"

