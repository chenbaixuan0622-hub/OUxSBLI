set format x "%4.0f"
set key off

set xlabel "{/=18 {/Times-New-Roman:Italic t / t_c}}"
set xrange [0:20]
set yrange [0:]

set tics font "Times-New-Roman,16"
set xtics nomirror
set ytics nomirror

# parameters
T = 530 * 5 / 9
S = 111
mu = 1.716e-5 * (273.2 + S) / (T + S) * (T / 273.2)**1.5
Re = 1600
L = 1.524e-3
M = 0.1
gamma = 1.4
R = 287.03
V = M * sqrt(gamma * R * T)
RHO = mu * Re / (V * L)
nu = mu / RHO

# data
e = "enstrophy.d"
KE = "kinetic_energy.d"

# convective time scale
tc = L / V

Lx = 2 * 3.14 * L
CFL = 0.03
nx = 128
dx = Lx / nx
# dimensional time
dt0 = CFL * dx / V
# non-dimensional time
dt1 = V * dt0 / L
nt = 20 / (100 * dt1)
dt = nt * dt0  

# set init value
y = 0

set ylabel "{/=18 {/Times-New-Roman:Italic Dissipation Rate}}"
set size square
set mxtics 5
set mytics 5
Lt = 22714.3668

set term pngcairo size 480, 480
set output "ke.png"
set ylabel "{/=18 {/Times-New-Roman:Italic Kinetic Energy}}"
set format y "%4.2f"
plot KE using ($1*Lt):($2/(RHO*V**2)) with points pt 2 lc rgb "blue"

set term pngcairo size 480, 480
set output "enstrophy.png"
set ylabel "{/=18 {/Times-New-Roman:Italic Enstrophy}}"
set format y "%4.0f"
plot e using ($1*Lt):($2*(tc**2)/RHO) with points pt 2 lc rgb "red"

set term pngcairo size 480, 480
set output "dissipation_rate.png"
set ylabel "{/=18 {/Times-New-Roman:Italic Dissipation Rate}}"
set format y "%4.2f"
plot KE using ($1*Lt):((dy=($2-y)/dt, y=$2), -tc*dy/dt/(RHO*V**2)) every ::1 smooth unique with points pt 2 lc rgb "blue"

