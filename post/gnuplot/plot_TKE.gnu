#! /usr/bin/gnuplot/
set key off

set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 5

UP1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/TKE_Budget_upstream.d"
RR1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/TKE_Budget_relaxation_region.d"
UP2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/TKE_Budget_upstream.d"
RR2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/TKE_Budget_relaxation_region.d"

set output "TKE_Budget.png"
set term pngcairo size 960, 480

set multiplot layout 1, 2

set size square
set ylabel "{/=18{/Times-New-Roman:Italic Budget terms}}"
set format y "%.1f"
set yrange [-0.3:0.3]
set ytics 0.1
set mytics 2
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set format x "%.0f"
set key right top
set key spacing 1.2

set xtics 10
set mxtics 2
set xrange [1:50]
plot UP1 using 1:2     with lines dt 1 lw 2 lc rgb "black" title "{/=16{/Times-New-Roman:Italic P}}", \
     UP1 using 1:3     with lines dt 1 lw 2 lc rgb "blue" title "{/=16{/Times-New-Roman:Italic T}}", \
     UP1 using 1:5     with lines dt 1 lw 2 lc rgb "dark-green" title "{/=16{/Times-New-Roman:Italic D}}", \
     UP1 using 1:(-$6) with lines dt 1 lw 2 lc rgb "0x6400B4" title "{/=16{/Symbol:Italic -e}}", \
     UP2 using 1:2     with lines dt 2 lw 2 lc rgb "black" notitle, \
     UP2 using 1:3     with lines dt 2 lw 2 lc rgb "blue" notitle, \
     UP2 using 1:5     with lines dt 2 lw 2 lc rgb "dark-green" notitle, \
     UP2 using 1:(-$6) with lines dt 2 lw 2 lc rgb "0x6400B4" notitle

set xtics 100
set mxtics 2
set xrange [1:400]
plot RR1 using 1:2     with lines dt 1 lw 2 lc rgb "black" title "{/=16{/Times-New-Roman:Italic P}}", \
     RR1 using 1:3     with lines dt 1 lw 2 lc rgb "blue" title "{/=16{/Times-New-Roman:Italic T}}", \
     RR1 using 1:5     with lines dt 1 lw 2 lc rgb "dark-green" title "{/=16{/Times-New-Roman:Italic D}}", \
     RR1 using 1:(-$6) with lines dt 1 lw 2 lc rgb "0x6400B4" title "{/=16{/Symbol:Italic -e}}", \
     RR2 using 1:2     with lines dt 2 lw 2 lc rgb "black" notitle, \
     RR2 using 1:3     with lines dt 2 lw 2 lc rgb "blue" notitle, \
     RR2 using 1:5     with lines dt 2 lw 2 lc rgb "dark-green" notitle, \
     RR2 using 1:(-$6) with lines dt 2 lw 2 lc rgb "0x6400B4" notitle

unset multiplot

