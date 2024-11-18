set size square
set term pngcairo size 480, 480
set output "vrms.png"
set format x "%.0f"
set format y "%.0f"
set key off
set xlabel "{/=18 {/Times-New-Roman:Italic y^+}}"
set ylabel "{/=18 {/Times-New-Roman:Italic u^+_{rms}, v^+_{rms}}}"
#set xrange [1:150]
#set yrange [0:10]

set tics font "Times-New-Roman, 16"
set xtics nomirror
set ytics nomirror
set mxtics 5
set mytics 2

data = "../np_data/Qrms/velocity_rms.d"

plot data using 1:2 with points pt 2 lc rgb "blue", \
     data using 1:3 with points pt 2 lc rgb "red"

