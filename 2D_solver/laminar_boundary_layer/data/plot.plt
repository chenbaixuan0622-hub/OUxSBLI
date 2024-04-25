set format xy "%4.2f"
set key off

set size ratio 1.5

set xlabel "{/=12 {/Times-New-Roman:Italic u / U_0}}"

set xrange [0:1]
set yrange [0:7]

set tics font "Times-New-Roman,10"
set xtics nomirror
set ytics nomirror

numerical = "boundary_layer.d"
theoretical = "blasius.d"

# boundary layer thickness
delta = 0.00175

set ylabel "{/=12 {/Symbol:Italic h} = {/Times-New-Roman:Italic y} / {/Symbol:Italic d}}"
set mxtics 2
set mytics 2
plot numerical using 2:($1/delta) with points pt 2 lc rgb "blue"
replot theoretical using 2:1 with lines linecolor rgb "red"

set term pngcairo size 480, 640
set output "output.png"
replot

