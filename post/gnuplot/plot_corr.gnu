#! /usr/bin/gnuplot/
set tics font "Times-New-Roman, 16"
set key off
set xtics nomirror
set ytics nomirror
set mxtics 5

#cxup1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_x_upstream.d"
#cxdw1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_x_downstream.d"
#cxsb1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_x_separation_bubble.d"
#cxlr1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_x_relaxation_region.d"
#czup1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_z_upstream.d"
#czdw1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_z_downstream.d"
#czsb1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_z_separation_bubble.d"
#czlr1 = "../SBLI_1delta/stat03ms_05ms_SLAU/post/corr_z_relaxation_region.d"
#aup1  = "../SBLI_1delta/stat03ms_05ms_SLAU/post/auto_corr_upstream.d"
#adw1  = "../SBLI_1delta/stat03ms_05ms_SLAU/post/auto_corr_downstream.d"
#asb1  = "../SBLI_1delta/stat03ms_05ms_SLAU/post/auto_corr_separation_bubble.d"
#alr1  = "../SBLI_1delta/stat03ms_05ms_SLAU/post/auto_corr_relaxation_region.d"

cxup1 = "../corr/corr_x_upstream.d"
cxdw1 = "../corr/corr_x_downstream.d"
cxsb1 = "../corr/corr_x_separation_bubble.d"
cxlr1 = "../corr/corr_x_relaxation_region.d"
czup1 = "../corr/corr_z_upstream.d"
czdw1 = "../corr/corr_z_downstream.d"
czsb1 = "../corr/corr_z_separation_bubble.d"
czlr1 = "../corr/corr_z_relaxation_region.d"
aup1  = "../corr/auto_corr_upstream.d"
adw1  = "../corr/auto_corr_downstream.d"
asb1  = "../corr/auto_corr_separation_bubble.d"
alr1  = "../corr/auto_corr_relaxation_region.d"

cxup2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_x_upstream.d"
cxdw2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_x_downstream.d"
cxsb2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_x_separation_bubble.d"
cxlr2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_x_relaxation_region.d"
czup2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_z_upstream.d"
czdw2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_z_downstream.d"
czsb2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_z_separation_bubble.d"
czlr2 = "../SBLI_05delta/stat03ms_09ms_SLAU/post/corr_z_relaxation_region.d"
aup2  = "../SBLI_05delta/stat03ms_09ms_SLAU/post/auto_corr_upstream.d"
adw2  = "../SBLI_05delta/stat03ms_09ms_SLAU/post/auto_corr_downstream.d"
asb2  = "../SBLI_05delta/stat03ms_09ms_SLAU/post/auto_corr_separation_bubble.d"
alr2  = "../SBLI_05delta/stat03ms_09ms_SLAU/post/auto_corr_relaxation_region.d"

delta  = 2.e-3

up1    = 6
dw1    = 15
sb1    = 16
lr1    = 21

up2    = 6
dw2    = 14
sb2    = 16
lr2    = 21

up4    = 6
dw4    = 14
sb4    = 22
lr4    = 30

t1     = 1e3
t2     = 3e3
t3     = 4e3

set output "corr.png"
set term pngcairo size 2400, 480

set size square
set ylabel "{/=18{/Times-New-Roman:Italic Correlation}}"
set format y "%.1f"
set ytics 0.2
set mytics 2
set format x "%.1f"
set xtics 0.5
set mxtics 2
set yrange [-0.3:1.1]

set multiplot layout 1, 5

set xlabel "{/=18 {/Times-New-Roman:Italic x /{/Symbol:Italic d}}}"
set xrange [0:2]
plot cxup1 using ($1/delta-up1):2 with lines dt 1 lw 2 lc rgb "black", \
     cxdw1 using ($1/delta-dw1):2 with lines dt 1 lw 2 lc rgb "blue", \
     cxlr1 using ($1/delta-lr1):2 with lines dt 1 lw 2 lc rgb "dark-green", \
     cxup2 using ($1/delta-up2):2 with lines dt 2 lw 2 lc rgb "black", \
     cxdw2 using ($1/delta-dw2):2 with lines dt 2 lw 2 lc rgb "blue", \
     cxlr2 using ($1/delta-lr2):2 with lines dt 2 lw 2 lc rgb "dark-green"

plot cxup1 using ($1/delta-up1):3 with lines dt 1 lw 2 lc rgb "black", \
     cxdw1 using ($1/delta-dw1):3 with lines dt 1 lw 2 lc rgb "blue", \
     cxlr1 using ($1/delta-lr1):3 with lines dt 1 lw 2 lc rgb "dark-green", \
     cxup2 using ($1/delta-up2):3 with lines dt 2 lw 2 lc rgb "black", \
     cxdw2 using ($1/delta-dw2):3 with lines dt 2 lw 2 lc rgb "blue", \
     cxlr2 using ($1/delta-lr2):3 with lines dt 2 lw 2 lc rgb "dark-green"

set xlabel "{/=18 {/Times-New-Roman:Italic z /{/Symbol:Italic d}}}"
set format x "%.2f"
set xrange [0:0.25]
set xtics 0.05
plot czup1 using ($1/delta):2 with lines dt 1 lw 2 lc rgb "black", \
     czdw1 using ($1/delta):2 with lines dt 1 lw 2 lc rgb "blue", \
     czlr1 using ($1/delta):2 with lines dt 1 lw 2 lc rgb "dark-green", \
     czup2 using ($1/delta):2 with lines dt 2 lw 2 lc rgb "black", \
     czdw2 using ($1/delta):2 with lines dt 2 lw 2 lc rgb "blue", \
     czlr2 using ($1/delta):2 with lines dt 2 lw 2 lc rgb "dark-green"

plot czup1 using ($1/delta):3 with lines dt 1 lw 2 lc rgb "black", \
     czdw1 using ($1/delta):3 with lines dt 1 lw 2 lc rgb "blue", \
     czlr1 using ($1/delta):3 with lines dt 1 lw 2 lc rgb "dark-green", \
     czup2 using ($1/delta):3 with lines dt 2 lw 2 lc rgb "black", \
     czdw2 using ($1/delta):3 with lines dt 2 lw 2 lc rgb "blue", \
     czlr2 using ($1/delta):3 with lines dt 2 lw 2 lc rgb "dark-green"

set ylabel "{/=18{/Times-New-Roman:Italic ACF}}"
set xlabel "{/=18 {/Times-New-Roman:Italic t}}"
set format x "%.3f"
set xrange [0:0.01]
set xtics 0.0025
plot aup1 using ($1*t1):2 with lines dt 1 lw 2 lc rgb "black", \
     adw1 using ($1*t1):2 with lines dt 1 lw 2 lc rgb "blue", \
     alr1 using ($1*t1):2 with lines dt 1 lw 2 lc rgb "dark-green", \
     aup2 using ($1*t2):2 with lines dt 2 lw 2 lc rgb "black", \
     adw2 using ($1*t2):2 with lines dt 2 lw 2 lc rgb "blue", \
     alr2 using ($1*t2):2 with lines dt 2 lw 2 lc rgb "dark-green"

unset multiplot

