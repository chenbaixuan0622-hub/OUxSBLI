set size square
set output "loglayer.png"
set term pngcairo size 480, 480
set logscale x
set format x "10^{%L}"
set format y "%.0f"
set key off
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set ylabel "{/=18{/Times-New-Roman:Italic u^+}}"
set xrange [1:1e4]
set yrange [0:30]

set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 5
set mytics 2

data = "../np_data/Qvd/u_y_plus.d"

plot data using 1:2 with points pt 2 lc rgb "blue", \
     x with lines dashtype 2 lw 1 lc rgb "black", \
     log(x) / 0.38 + 4.1 with lines dashtype 2 lw 1 lc rgb "black"

