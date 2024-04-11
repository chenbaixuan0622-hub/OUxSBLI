set term pngcairo size 1280, 480
set output "output.png"

set format xy "%4.2f"
set key off

set xlabel "{/=12 {/Times-New-Roman:Italic tM_0}}"

set xrange [0:140]

set tics font "Times-New-Roman,10"
set xtics nomirror
set ytics nomirror

entropy_3 = "entropy_3.d"
entropy_4 = "entropy_4.d"

ke_3 = "kinetic_energy_3.d"
ke_4 = "kinetic_energy_4.d"

set multiplot layout 1, 2

set ylabel "{/=12 {/Times-New-Roman:Italic ({/Symbol:Italic r}s - {/Symbol:Italic r_0}s_0)/{/Symbol:Italic}s_0}}"
set size square
set mxtics 10
set mytics 10
plot entropy_3 using 1:2 with points pt 2 lc rgb "blue", entropy_4 using 1:2 with points pt 2 lc rgb "red"

set ylabel "{/=12 {/Times-New-Roman:Italic {/Symbol:Italic r}k / {/Symbol:Italic r_0}k_0}}"
set size square
set mxtics 10
set mytics 10
plot ke_3 using 1:2 with points pt 2 lc rgb "blue", ke_4 using 1:2 with points pt 2 lc rgb "red"

unset multiplot

set term x11
replot

