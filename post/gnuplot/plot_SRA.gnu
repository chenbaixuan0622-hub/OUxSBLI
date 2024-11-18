#! /usr/bin/gnuplot/
set key off

set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 5

data = "SRA_upstream.d"

set output "SRA_upstream.png"
set term pngcairo size 960, 480

set multiplot layout 1, 2

set size square
set format y "%.1f"
set ytics 0.2
set mytics 2
set format x "%.1f"
set xlabel "{/=18 {/Times-New-Roman:Italic y /{/Symbol:Italic d}}}"
set xtics 0.2
set mxtics 2
set xrange [0:1]
set yrange [-0.5:1]
set key right bottom
set key spacing 1.5
plot data  using 1:(-$2) with lines dt 1 lw 2 lc rgb "blue"  title "{/=16{/Times-New-Roman:Italic  -R_{uv}}}", \
     data  using 1:3      with lines dt 2 lw 2 lc rgb "blue"  title "{/=16{/Times-New-Roman:Italic  R_{vT}}}"

set size square
set format y "%.1f"
set ytics 0.2
set mytics 2
set format x "%.1f"
set xlabel "{/=18 {/Times-New-Roman:Italic y /{/Symbol:Italic d}}}"
set xtics 0.2
set mxtics 2
set xrange [0:1]
set yrange [-0.1:1]
set key right bottom
set key spacing 1.5
plot data using 1:(-$4) with lines dt 1 lw 2 lc rgb "blue" title "{/=16{/Times-New-Roman:Italic -R_{uT}}}"

unset multiplot

