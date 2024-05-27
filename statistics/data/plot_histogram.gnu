# plot_normal_distribution.gnu
set terminal pngcairo size 800, 600
set output 'normal_distribution.png'

set xlabel "Value"
set ylabel "Density"

# Get the number of data
stats 'output.d' nooutput

# Bin width and range for the histgram
binwidth = 0.1
bin(x,width) = width*floor(x/width) + binwidth/2

# Plot the histgram
plot 'output.d' using (bin($1,binwidth)):(1.0/(STATS_records*binwidth)) smooth freq with boxes lc rgb "blue"
