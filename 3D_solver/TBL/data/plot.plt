set term pngcairo size 640, 640
set output "loglayer.png"
set size square
set logscale x
set format x "10^{%L}"
set format y "%.0f"
set key off  # 修正済み
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"  # 修正済み
set ylabel "{/=18{/Times-New-Roman:Italic u^+}}"  # 修正済み
set xrange [1:1e4]
set yrange [0:30]

set tics font "Times-New-Roman, 16"  # 修正済み
set xtics nomirror
set ytics nomirror
set mxtics 2
set mytics 2

data = "yplus.d"

plot data using 1:2 with points pt 2 lc rgb "red"

