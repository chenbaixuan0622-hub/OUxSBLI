# plot_normal_distribution.gnu
set format x "%2.1f"
set format y "%.1f"
set key off

set xlabel "{/=18 {/Times-New-Roman:Italic x}}"
set ylabel "{/=18 {/Times-New-Roman:Italic f(x)}}"

set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 2
set mytics 2

pdf1 = "pdf1.d"
pdf2 = "pdf2.d"

plot pdf1 using 1:2 with lines lc rgb "blue"
replot pdf2 using 1:2 with lines lc rgb "red"

set term pngcairo size 800, 600
set output 'pdf.png'
replot

