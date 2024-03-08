<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
</head>

<body>
<h1>
Compressible solver
</h1>
<p>
2D and 3D solvers are accelerated by GPU
</p>

<h2>
KEEP scheme
</h2>
<p>
Taylor-Green vortex was calculated. 
</p>
<img src="./visuals/Vorticity.gif" alt="Taylor-Green vortex">

<h2>
SLAU scheme
</h2>
<p>
1 dimensional shock tube was calculated.
</p>
<img src="./visuals/shock_tube_SLAU.png" alt="shock-tube">

<h2>
Dependency
</h2>
<ul>
<li>nvfortran</li>
<li>ParaView</li>
<li>gnuplot</li>
</ul>

<h2>
Usage
</h2>
<ul>
<li> Edit input.d file grid, physical properties, and simulation time etc. </li>
<li> Edit mod_globals.f90 file to choose accuracy. 2nd-order accuracy and 4th-order accuracy are available. </li> 
<li> Edit mod_globals.f90 file to choose Euler solver or Navier-Stokes solver. Moreover, turbulent model can be used. </li>
<li> Edit calc_time_dev.f90 to choose appropriate gridDim and blockDim. The number of threads should be a multiple of 32 to make the most of GPU. </li>
<li><pre>$ make</pre></li>
</ul>

<h2>
Reference
</h2>
<ul>
<li><a href="https://www.sciencedirect.com/science/article/abs/pii/S0021999118305916">
Yuichi Kuya, Kosuke Totani, Soshi Kawai, Kinetic energy and entropy preserving schemes for compressible flows by split convective forms, Journal of Computational Physics, 2018</a>
</li>
<li><a href="https://www.sciencedirect.com/science/article/abs/pii/S0021999121003776">
Yuichi Kuya, Soshi Kawai, High-order accurate kinetic-enrgy and entropy preserving (KEEP) schemes on curvilinear grids, Journal of Computational Physics, 2021</a>
</li>
<li><a href="https://www.sciencedirect.com/science/article/pii/S187775031630299X?via%3Dihub">
Christian T. Jacobs, Satya P. Jammy, Neil D. Sandham, OpenSBLI: A framework for the automated derivation and parallel execution of finite difference solvers on a range of computer architectures, 2017</a>
</li>
<li><a href="https://arc.aiaa.org/doi/abs/10.2514/6.2009-3797">
Yves Allaneau, Antony Jameson, Direct Numerical Simulations of a Two-Dimensional Viscous Flow in a Shocktube Using Kinetic Energy Preserving Scheme, 2012</a>
</li>
<li><a href="https://docs.nvidia.com/hpc-sdk/pgi-compilers/2017/pgi17cudaforug.pdf">
CUDA FORTRAN PROGRAMMING GUIDE AND REFERENCE, 2017</a>
</li>
</ul>
</body>
</html>
