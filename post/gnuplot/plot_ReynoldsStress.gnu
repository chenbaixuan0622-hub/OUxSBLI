#! /usr/bin/gnuplot/
set key off

set output "ReynoldsStress.png"
set term pngcairo size 960, 960
set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 5
set size square
set format y "%.0f"
set ytics 1
set mytics 2
set logscale x
set format x "10^{%L}"
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set xrange [1:2e3]

up1 = "ReynoldsStress_upstream.d"
dw1 = "ReynoldsStress_downstream.d"
rr1 = "ReynoldsStress_relaxation_region.d"
up2 = "ReynoldsStress_upstream.d"
dw2 = "ReynoldsStress_downstream.d"
rr2 = "ReynoldsStress_relaxation_region.d"

set multiplot layout 2, 2

set ylabel "{/=18{/Times-New-Roman:Italic uu}}"
plot up1 using 1:2 with lines dt 1 lw 3 lc rgb "black", \
     dw1 using 1:2 with lines dt 1 lw 3 lc rgb "blue", \
     rr1 using 1:2 with lines dt 1 lw 3 lc rgb "dark-green", \
     up2 using 1:2 with lines dt 2 lw 3 lc rgb "black", \
     dw2 using 1:2 with lines dt 2 lw 3 lc rgb "blue", \
     rr2 using 1:2 with lines dt 2 lw 3 lc rgb "dark-green", \

set ylabel "{/=18{/Times-New-Roman:Italic vv}}"
plot up1 using 1:3 with lines dt 1 lw 3 lc rgb "black", \
     dw1 using 1:3 with lines dt 1 lw 3 lc rgb "blue", \
     rr1 using 1:3 with lines dt 1 lw 3 lc rgb "dark-green", \
     up2 using 1:3 with lines dt 2 lw 3 lc rgb "black", \
     dw2 using 1:3 with lines dt 2 lw 3 lc rgb "blue", \
     rr2 using 1:3 with lines dt 2 lw 3 lc rgb "dark-green", \

set ylabel "{/=18{/Times-New-Roman:Italic ww}}"
plot up1 using 1:4 with lines dt 1 lw 3 lc rgb "black", \
     dw1 using 1:4 with lines dt 1 lw 3 lc rgb "blue", \
     rr1 using 1:4 with lines dt 1 lw 3 lc rgb "dark-green", \
     up2 using 1:4 with lines dt 2 lw 3 lc rgb "black", \
     dw2 using 1:4 with lines dt 2 lw 3 lc rgb "blue", \
     rr2 using 1:4 with lines dt 2 lw 3 lc rgb "dark-green", \

set ylabel "{/=18{/Times-New-Roman:Italic uv}}"
plot up1 using 1:5 with lines dt 1 lw 3 lc rgb "black", \
     dw1 using 1:5 with lines dt 1 lw 3 lc rgb "blue", \
     rr1 using 1:5 with lines dt 1 lw 3 lc rgb "dark-green", \
     up2 using 1:5 with lines dt 2 lw 3 lc rgb "black", \
     dw2 using 1:5 with lines dt 2 lw 3 lc rgb "blue", \
     rr2 using 1:5 with lines dt 2 lw 3 lc rgb "dark-green", \

unset multiplot

