#! /usr/bin/gnuplot/
set key off

set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 5

data = "turb_stat.d"

ruu = "ruu.csv"
rvv = "rvv.csv"
ruv = "ruv.csv"

urms  = "urms.csv"
vrms  = "vrms.csv"
uvrms = "uvrms.csv"

set output "turb_stat.png"
set term pngcairo size 2500, 960

set multiplot layout 2, 3

set size square
set ylabel "{/=18{/Times-New-Roman:Italic Reynolds stress}}"
set format y "%.0f"
set ytics 1
set mytics 2
set logscale x
set format x "10^{%L}"
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set xrange [1:1e4]
set yrange [-2:3]
plot ruu  using 1:2  with points pt 4 lc rgb "black", \
     rvv  using 1:2  with points pt 4 lc rgb "black", \
     ruv  using 1:2  with points pt 4 lc rgb "black", \
     data using 2:8  with lines dt 1 lw 2 lc rgb "blue", \
     data using 2:9  with lines dt 2 lw 2 lc rgb "blue", \
     data using 2:10 with lines dt 3 lw 2 lc rgb "blue"

set size square
set ylabel "{/=18{/Times-New-Roman:Italic Reynolds stress}}"
set format y "%.1f"
set ytics 0.2
set mytics 2
set logscale x
set format x "10^{%L}"
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set xrange [1:1e4]
set yrange [0:2]
plot data using 2:(-$10) with lines lw 2 lc rgb "blue", \
     data using 2:11 with lines lw 2 lc rgb "red", \
     data using 2:(-$10+$11) with lines lw 2 lc rgb "black"

set size square
set ylabel "{/=18{/Symbol:Italic r_{/Times-New-Roman:Italic RMS}}}"
set format y "%.0f"
set ytics 5
set mytics 2
set logscale x
set format x "10^{%L}"
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set xrange [1:1e4]
set yrange [0:15]
plot data using 2:12 with lines lw 2 lc rgb "blue"

set size square
set ylabel "{/=18{/Times-New-Roman:Italic P_{RMS}}}"
set ytics 1
set mytics 2
set logscale x
set format x "10^{%L}"
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set xrange [1:1e4]
set yrange [0:5]
plot data using 2:13 with lines lw 2 lc rgb "blue"

set size square
set ylabel "{/=18{/Times-New-Roman:Italic K^+}}"
set ytics 1
set mytics 2
set logscale x
set format x "10^{%L}"
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set xrange [1:1e4]
set yrange [0:6]
plot data using 2:7 with lines lw 2 lc rgb "blue"

set size ratio 0.5
set ylabel "{/=18{/Times-New-Roman:Italic Turbulence intensity}}"
set ytics 1
set mytics 2
unset logscale x
set format x "%.1f"
set xlabel "{/=18 {/Times-New-Roman:Italic y /{/Symbol:Italic d}}}"
set xrange [0:1.2]
set yrange [-2:3]
plot urms  using 1:2 with points pt 4 lc rgb "black", \
     vrms  using 1:2 with points pt 4 lc rgb "black", \
     uvrms using 1:2 with points pt 4 lc rgb "black", \
     data  using 1:4 with lines dt 1 lw 2 lc rgb "blue", \
     data  using 1:5 with lines dt 2 lw 2 lc rgb "blue", \
     data  using 1:6 with lines dt 3 lw 2 lc rgb "blue"

unset multiplot

