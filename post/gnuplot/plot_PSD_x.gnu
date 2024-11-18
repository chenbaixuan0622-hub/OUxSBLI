#! /usr/bin/gnuplot/
set tics font "Times-New-Roman, 16"
set key off
set xtics nomirror
set ytics nomirror
set mxtics 5
set mytics 5
set logscale xy
set size square

PSD_x_1_up = "PSD_x_1delta_upstream.d"
PSD_x_1_dw = "PSD_x_1delta_downstream.d"
PSD_x_1_sb = "PSD_x_1delta_separation_bubble.d"
PSD_x_1_rr = "PSD_x_1delta_relaxation_region.d"
PSD_x_2_up = "PSD_x_05delta_upstream.d"
PSD_x_2_dw = "PSD_x_05delta_downstream.d"
PSD_x_2_sb = "PSD_x_05delta_separation_bubble.d"
PSD_x_2_rr = "PSD_x_05delta_relaxation_region.d"
PSD_z_1_up = "PSD_z_1delta_upstream.d"
PSD_z_1_dw = "PSD_z_1delta_downstream.d"
PSD_z_1_sb = "PSD_z_1delta_separation_bubble.d"
PSD_z_1_rr = "PSD_z_1delta_relaxation_region.d"
PSD_z_2_up = "PSD_z_05delta_upstream.d"
PSD_z_2_dw = "PSD_z_05delta_downstream.d"
PSD_z_2_sb = "PSD_z_05delta_separation_bubble.d"
PSD_z_2_rr = "PSD_z_05delta_relaxation_region.d"

set output "PSD_xz.png"
set term pngcairo size 1440, 480

set format y "10^{%L}"
set format x "10^{%L}"

set multiplot layout 1, 4

set xlabel "{/=18 {/Times-New-Roman:Italic k} {/Times-New-Roman: [1/m]}}"
set ylabel "{/=18 {/Times-New-Roman:Italic PSD} {/Times-New-Roman: [m^3/s^2]}}"
set xrange [2e2:4e4]
plot PSD_x_1_up using 1:2 with lines dt 1 lw 2 lc rgb "black", \
     PSD_x_1_dw using 1:2 with lines dt 1 lw 2 lc rgb "blue", \
     PSD_x_1_rr using 1:2 with lines dt 1 lw 2 lc rgb "dark-green", \
     PSD_x_2_up using 1:2 with lines dt 2 lw 2 lc rgb "black", \
     PSD_x_2_dw using 1:2 with lines dt 2 lw 2 lc rgb "blue", \
     PSD_x_2_rr using 1:2 with lines dt 2 lw 2 lc rgb "dark-green"

plot PSD_x_1_up using 1:3 with lines dt 1 lw 2 lc rgb "black", \
     PSD_x_1_dw using 1:3 with lines dt 1 lw 2 lc rgb "blue", \
     PSD_x_1_rr using 1:3 with lines dt 1 lw 2 lc rgb "dark-green", \
     PSD_x_2_up using 1:3 with lines dt 2 lw 2 lc rgb "black", \
     PSD_x_2_dw using 1:3 with lines dt 2 lw 2 lc rgb "blue", \
     PSD_x_2_rr using 1:3 with lines dt 2 lw 2 lc rgb "dark-green"

set xrange [4e2:1e5]
plot PSD_z_1_up using 1:2 with lines dt 1 lw 2 lc rgb "black", \
     PSD_z_1_dw using 1:2 with lines dt 1 lw 2 lc rgb "blue", \
     PSD_z_1_rr using 1:2 with lines dt 1 lw 2 lc rgb "dark-green", \
     PSD_z_2_up using 1:2 with lines dt 2 lw 2 lc rgb "black", \
     PSD_z_2_dw using 1:2 with lines dt 2 lw 2 lc rgb "blue", \
     PSD_z_2_rr using 1:2 with lines dt 2 lw 2 lc rgb "dark-green"

plot PSD_z_1_up using 1:3 with lines dt 1 lw 2 lc rgb "black", \
     PSD_z_1_dw using 1:3 with lines dt 1 lw 2 lc rgb "blue", \
     PSD_z_1_rr using 1:3 with lines dt 1 lw 2 lc rgb "dark-green", \
     PSD_z_2_up using 1:3 with lines dt 2 lw 2 lc rgb "black", \
     PSD_z_2_dw using 1:3 with lines dt 2 lw 2 lc rgb "blue", \
     PSD_z_2_rr using 1:3 with lines dt 2 lw 2 lc rgb "dark-green"

unset multiplot

