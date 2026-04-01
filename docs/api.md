# Project API Index
> **AI Agent Instruction:** This is an auto-generated API reference for the CUDA Fortran project. Use this to understand the project structure, modules, subroutines, and their arguments without reading the full implementation details.

## File: `calc_time_dev.f90`

### `module calc_time_dev`
- **Description:** Implements 3rd-order and 4th-order Runge-Kutta time stepping with MPI/GPU support

### `module procedure`

### `subroutine RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)`
- **Description:** GPU computation: Each rank manages one GPU asynchronously; MPI sync only for I/O
- **Arguments & Variables:**
  - `id_RungeKutta` (`integer(2), intent(in)`) : time integration method ID
  - `id_rescale` (`integer(2), intent(in)`) : rescaling method ID
  - `myrank` (`integer, intent(in)`) : MPI rank
  - `mygpu` (`integer, intent(in)`) : GPU index for this rank
  - `nx` (`integer, intent(in)`) : x grid dimension
  - `ny` (`integer, intent(in)`) : y grid dimension
  - `nz` (`integer, intent(in)`) : z grid dimension
  - `x(nx)` (`real(8), intent(in)`) : x coordinates
  - `dx_cpu(nx-1)` (`real(8), intent(in)`) : inverse x spacing (host)
  - `y(ny)` (`real(8), intent(in)`) : y coordinates
  - `dy_cpu(ny-1)` (`real(8), intent(in)`) : inverse y spacing (host)
  - `z(nz)` (`real(8), intent(in)`) : z coordinates
  - `dz_cpu(nz-1)` (`real(8), intent(in)`) : inverse z spacing (host)
  - `Jacobian_cpu(nx,ny)` (`real(8), intent(in)`) : Jacobian determinant (host)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`) : conservative variables on host
  - `ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)` (`real(8), allocatable, device`)
  - `T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)` (`real(8), allocatable, device`)
  - `dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)` (`real(8), allocatable, device`)
  - `ke0 = 1.d0, entropy0 = 1.d0` (`real(4)`)

### `subroutine RungeKutta_3rd_rescale(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)`
- **Description:** Similar TVD RK3 stages as above, plus calls to step_rescale() for non-conservative correction
- **Arguments & Variables:**
  - `id_RungeKutta` (`integer(2), intent(in)`)
  - `id_rescale` (`integer(4), intent(in)`)
  - `myrank, mygpu, nx, ny, nz` (`integer, intent(in)`)
  - `x(nx), dx_cpu(nx-1)` (`real(8), intent(in)`)
  - `y(ny), dy_cpu(ny-1)` (`real(8), intent(in)`)
  - `z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)` (`real(8), intent(in)`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`)
  - `step = 1, flag_re = 0, flag_req` (`integer`)
  - `Qre(:), Qm(:)` (`real(8), allocatable, device`)
  - `Qm_cpu(:)` (`real(8), allocatable, pinned`)
  - `ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)` (`real(8), allocatable, device`)
  - `T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)` (`real(8), allocatable, device`)
  - `dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)` (`real(8), allocatable, device`)
  - `ke0 = 1.d0, entropy0 = 1.d0` (`real(4)`)

### `subroutine RungeKutta_4th(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)`
- **Description:** More accurate than RK3 but lacks TVD property; requires smaller CFL (~0.8 vs 1.0)
- **Arguments & Variables:**
  - `id_RungeKutta` (`integer(4), intent(in)`)
  - `id_rescale` (`integer(2), intent(in)`)
  - `myrank, mygpu, nx, ny, nz` (`integer, intent(in)`)
  - `x(nx), dx_cpu(nx-1)` (`real(8), intent(in)`)
  - `y(ny), dy_cpu(ny-1)` (`real(8), intent(in)`)
  - `z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)` (`real(8), intent(in)`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`)
  - `ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)` (`real(8), allocatable, device`)
  - `T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)` (`real(8), allocatable, device`)
  - `dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)` (`real(8), allocatable, device`)
  - `ke0 = 1.d0, entropy0 = 1.d0` (`real(4)`)

### `subroutine RungeKutta_4th_rescale(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)`
- **Description:** Step-rescale calls handle inter-rank communication for marking critical planes
- **Arguments & Variables:**
  - `id_RungeKutta` (`integer(4), intent(in)`)
  - `id_rescale` (`integer(4), intent(in)`)
  - `myrank, mygpu, nx, ny, nz` (`integer, intent(in)`)
  - `x(nx), dx_cpu(nx-1)` (`real(8), intent(in)`)
  - `y(ny), dy_cpu(ny-1)` (`real(8), intent(in)`)
  - `z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)` (`real(8), intent(in)`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`)
  - `step = 1, flag_re = 0, flag_req` (`integer`)
  - `Qre(:), Qm(:)` (`real(8), allocatable, device`)
  - `Qm_cpu(:)` (`real(8), allocatable, pinned`)
  - `ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)` (`real(8), allocatable, device`)
  - `T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)` (`real(8), allocatable, device`)
  - `dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)` (`real(8), allocatable, device`)
  - `ke0 = 1.d0, entropy0 = 1.d0` (`real(4)`)

## File: `calc_visc4.f90`

### `module calc_visc4`
- **Description:** Generally more accurate but requires larger stencils than 2nd-order

### `subroutine calc_Ev4(nx, ny, nz, dx, dy, dz, Q, T, mu, E)`
- **Description:** High-order accurate computation of viscous stress and heat flux
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(inout), device`) : viscous flux components in x direction
  - `u(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `v(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `w(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `uy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `vy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `uz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)` (`real(8), shared`)
  - `wz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)` (`real(8), shared`)
  - `txx, txy, txz, utxx, vtxy, wtxz, kTx` (`real(8)`)
  - `mu3(3)` (`real(8), device`)
  - `kTx3(3)` (`real(8), device`)
  - `my(2)` (`real(8), device`)
  - `mz(2)` (`real(8), device`)

### `subroutine calc_Ev_LES4(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E)`
- **Description:** CUDA Fortran kernel for 4th-order viscous flux with LES SGS model in x direction
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `mut(nx,ny,nz)` (`real(8), intent(in), device`) : turbulent eddy viscosity (LES model)
  - `qc2(nx,ny,nz)` (`real(8), intent(in), device`) : quadratic constitutive relation correction
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(inout), device`) : viscous + SGS flux in x direction
  - `u(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `v(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `w(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `uy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `vy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)` (`real(8), shared`)
  - `uz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)` (`real(8), shared`)
  - `wz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)` (`real(8), shared`)
  - `txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs` (`real(8)`)
  - `mu3(3), mut3(3)` (`real(8), device`)
  - `kTx3(3)` (`real(8), device`)
  - `H(4)` (`real(8)`)
  - `my, mysgs, mz, mzsgs` (`real(8), dimension(2), device`)
  - `my(2)` (`real(8), device`)
  - `mz(2)` (`real(8), device`)
  - `H(2)` (`real(8), device`)

### `subroutine calc_Fv4(nx, ny, nz, dy, dx, dz, Q, T, mu, F)`
- **Description:** CUDA Fortran kernel for 4th-order viscous flux in y direction
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(inout), device`) : viscous flux components in y direction
  - `u(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `v(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `w(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `ux(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `vx(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `vz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `wz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `tyx, tyy, tyz, utyx, vtyy, wtyz, kTy` (`real(8)`)
  - `mu3(3)` (`real(8), device`)
  - `kTy3(3)` (`real(8), device`)
  - `mx(2)` (`real(8), device`)
  - `mz(2)` (`real(8), device`)

### `subroutine calc_Fv_LES4(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F)`
- **Description:** CUDA Fortran kernel for 4th-order viscous flux with LES SGS model in y direction
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `mut(nx,ny,nz)` (`real(8), intent(in), device`) : turbulent eddy viscosity (LES model)
  - `qc2(nx,ny,nz)` (`real(8), intent(in), device`) : quadratic constitutive relation correction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(inout), device`) : viscous + SGS flux in y direction
  - `u(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `v(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `w(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `ux(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `vx(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `vz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `wz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)` (`real(8), shared`)
  - `tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs` (`real(8)`)
  - `mu3(3), mut3(3)` (`real(8), device`)
  - `kTy3(3)` (`real(8), device`)
  - `H(4)` (`real(8), device`)
  - `u2, v2, w2, mz, mzsgs, mx, mxsgs` (`real(8), dimension(2), device`)
  - `mx(2)` (`real(8), device`)
  - `mz(2)` (`real(8), device`)
  - `H(2)` (`real(8), device`)

### `subroutine calc_Gv4(nx, ny, nz, dx, dy, dz, Q, T, mu, G)`
- **Description:** CUDA Fortran kernel for 4th-order viscous flux in z direction
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(inout), device`) : viscous flux components in z direction
  - `u(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `v(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `w(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `wx(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `wy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `ux(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `vy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `tzx, tzy, tzz, utzx, vtzy, wtzz, kTz` (`real(8)`)
  - `mu3(3)` (`real(8), device`)
  - `kTz3(3)` (`real(8), device`)
  - `mx(2)` (`real(8), device`)
  - `my(2)` (`real(8), device`)

### `subroutine calc_Gv_LES4(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)` (`real(8), intent(in), device`)
  - `mut(nx,ny,nz), qc2(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(inout), device`)
  - `u(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `v(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `w(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `wx(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `wy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `ux(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `vy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)` (`real(8), shared`)
  - `tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs` (`real(8)`)
  - `mu3(3), mut3(3)` (`real(8), device`)
  - `kTz3(3)` (`real(8), device`)
  - `H(4)` (`real(8), device`)
  - `mx, mxsgs, my, mysgs` (`real(8), dimension(2), device`)
  - `mx(2)` (`real(8), device`)
  - `my(2)` (`real(8), device`)
  - `H(2)` (`real(8), device`)

## File: `set_bc_common.f90`

### `module set_bc_common`

### `module procedure`

### `subroutine set_bc_cyclic2_init(id_accuracy, nx, ny, nz, Q)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`)

### `subroutine set_bc_cyclic4_init(id_accuracy, nx, ny, nz, Q)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`)

### `subroutine set_bc_cyclic6_init(id_accuracy, nx, ny, nz, Q)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`)

### `subroutine set_bc_cyclic2(id_accuracy, nx, ny, nz, Q)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout), device`)

### `subroutine set_bc_cyclic4(id_accuracy, nx, ny, nz, Q)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout), device`)

### `subroutine set_bc_cyclic6(id_accuracy, nx, ny, nz, Q)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout), device`)

### `subroutine set_bc_cyclic_z(nx, ny, nz, QJ)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `QJ(5,nx,ny,nz)` (`real(8), intent(inout), device`)

### `subroutine set_bc_mut_common(nx, ny, nz, mut, qc2)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `mut(nx,ny,nz), qc2(nx,ny,nz)` (`real(8), intent(inout), device`)

## File: `calc_hybrid.f90`

### `module calc_hybrid`
- **Description:** Computes Ducros sensor for automatic scheme switching between KEEP and SLAU

## File: `calc_les.f90`

### `module calc_les`
- **Description:** Implements DNS/RANS and dynamic Smagorinsky LES turbulence models

## File: `calc_rescale.f90`

### `module calc_rescale`

### `subroutine calc_mean(step, ireq, flag_re, nx, ny, nz, Jacobian, QJ, Qm)`
- **Arguments & Variables:**
  - `step, ireq` (`integer, intent(inout)`)
  - `flag_re, nx, ny, nz` (`integer, intent(in)`)
  - `Jacobian(nx,ny), QJ(5,nx,ny,nz)` (`real(8), intent(in), device`)
  - `Qm(ny*5)` (`real(8), intent(inout), device`)

### `subroutine copy(nx, ny, nz, QJ, Qre)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in)`)
  - `QJ(5,nx,ny,nz)` (`real(8), intent(in), device`)
  - `Qre(ny*(nz-6)*5)` (`real(8), intent(out), device`)

### `subroutine step_rescale(num, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ, Qm, Qre)`
- **Arguments & Variables:**
  - `num, myrank, nx, ny, nz` (`integer, intent(in)`)
  - `step, flag_re, flag_req, ireq, ireq2(2)` (`integer, intent(inout)`)
  - `Jacobian(nx,ny), QJ(5,nx,ny,nz)` (`real(8), intent(in), device`)
  - `Qm(ny*5), Qre(ny*(nz-6)*5)` (`real(8), intent(inout), device`)

### `subroutine wait_rescale(myrank, ireq, ireq2, istat, istat2)`
- **Arguments & Variables:**
  - `myrank` (`integer, intent(in)`)
  - `ireq, ireq2(2), istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)` (`integer, intent(inout)`)

### `subroutine rescale_recv_send(num, flag_re, nx, ny, nz, step, y, Jacobian, Qm_cpu)`
- **Arguments & Variables:**
  - `num` (`integer, intent(in)`)
  - `flag_re` (`integer, intent(inout)`)
  - `nx, ny, nz, step` (`integer, intent(in)`)
  - `y(ny), Jacobian(nx,ny)` (`real(8), intent(in)`)
  - `Qm_cpu(ny*5)` (`real(8), intent(inout)`)
  - `Qre_cpu(ny*(nz-6)*5), bltre` (`real(8)`)
  - `Qre(ny*(nz-6)*5), Qm(ny*5)` (`real(8), device`)

### `subroutine write_Qm(ny, step, y, Qm)`
- **Arguments & Variables:**
  - `ny, step` (`integer, intent(in)`)
  - `y(ny), Qm(ny*5)` (`real(8), intent(in)`)

### `subroutine set_rescale(flag_re, step, nx, ny, nz, y, Jacobian, Qm, bltre, Qre)`
- **Arguments & Variables:**
  - `flag_re` (`integer, intent(inout)`)
  - `y(ny), Jacobian(nx,ny), Qm(ny*5)` (`real(8), intent(in)`)
  - `bltre` (`real(8), intent(out)`)
  - `jj_y, jj_e` (`integer, dimension(ny)`)
  - `Um, Vm, Wm, rhom, Tm, pm` (`real(8), dimension(ny)`)
  - `ufre, vfre, wfre, Tfre, pfre` (`real(8), dimension(ny,nz)`)
  - `ypre, ypin, etre, etin` (`real(8), dimension(ny)`)
  - `ufin, vfin, wfin, Tfin, pfin` (`real(8), dimension(ny,nz)`)
  - `ufout, vfout, wfout, Tfout, pfout` (`real(8), dimension(ny,nz)`)
  - `Umin, Vmin, pmin, Tmin, Umout, Vmout, pmout, Tmout` (`real(8), dimension(ny)`)
  - `weight` (`real(8), dimension(ny)`)
  - `u_tmp, v_tmp, p_tmp, T_tmp, weight_tmp, Jacobian_tmp, over_rhore, utin_nu, utre_nu` (`real(8)`)
  - `over_blt = 1.d0 / blt` (`real(8)`)
  - `blt_min = 0.5d0 * blt, blt_max = 2.d0 * blt` (`real(8)`)

## File: `calc_visc2.f90`

### `module calc_visc2`
- **Description:** Computes viscous stresses and heat flux for compressible Navier-Stokes equations

### `subroutine calc_Ev2(nx, ny, nz, dx, dy, dz, Q, T, mu, E)`
- **Description:** Accounts for molecular viscosity and thermal conductivity
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables: rho, rho*u, rho*v, rho*w, E
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(inout), device`) : viscous flux components in x direction
  - `u(threadsEv%x+1,0:threadsEv%y+1,0:threadsEv%z+1)` (`real(8), shared`)
  - `v(threadsEv%x+1,0:threadsEv%y+1,threadsEv%z)` (`real(8), shared`)
  - `w(threadsEv%x+1,0:threadsEv%z+1,threadsEv%y)` (`real(8), shared`)
  - `txx, txy, txz, utxx, vtxy, wtxz, kTx` (`real(8)`)

### `subroutine calc_Ev_LES2(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E)`
- **Description:** Computes molecular + subgrid-scale viscous stresses and heat flux
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `mut(nx,ny,nz)` (`real(8), intent(in), device`) : turbulent eddy viscosity (LES model)
  - `qc2(nx,ny,nz)` (`real(8), intent(in), device`) : quadratic constitutive relation correction
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(inout), device`) : viscous + SGS flux in x direction
  - `txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs` (`real(8)`)
  - `my, mysgs, mz, mzsgs` (`real(8), dimension(2), device`)
  - `my(2)` (`real(8), device`)
  - `mz(2)` (`real(8), device`)
  - `H(2)` (`real(8), device`)

### `subroutine calc_Fv2(nx, ny, nz, dy, dx, dz, Q, T, mu, F)`
- **Description:** Computes stress tensor components and heat flux at cell faces
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(inout), device`) : viscous flux components in y direction
  - `u(threadsFv%y+1,0:threadsFv%x+1,threadsFv%z)` (`real(8), shared`)
  - `v(threadsFv%y+1,0:threadsFv%x+1,0:threadsFv%z+1)` (`real(8), shared`)
  - `w(threadsFv%y+1,0:threadsFv%z+1,threadsFv%x)` (`real(8), shared`)
  - `tyx, tyy, tyz, utyx, vtyy, wtyz, kTy` (`real(8)`)

### `subroutine calc_Fv_LES2(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F)`
- **Description:** Computes molecular + subgrid-scale viscous stresses and heat flux
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `mut(nx,ny,nz)` (`real(8), intent(in), device`) : turbulent eddy viscosity (LES model)
  - `qc2(nx,ny,nz)` (`real(8), intent(in), device`) : quadratic constitutive relation correction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(inout), device`) : viscous + SGS flux in y direction
  - `tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs` (`real(8)`)
  - `u2, v2, w2, mz, mzsgs, mx, mxsgs` (`real(8), dimension(2), device`)
  - `mx(2)` (`real(8), device`)
  - `mz(2)` (`real(8), device`)
  - `H(2)` (`real(8), device`)

### `subroutine calc_Gv2(nx, ny, nz, dx, dy, dz, Q, T, mu, G)`
- **Description:** Computes stress tensor components and heat flux at cell faces
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing in x (1/dx)
  - `dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing in y (1/dy)
  - `dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing in z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature at grid points
  - `mu(nx,ny,nz)` (`real(8), intent(in), device`) : molecular viscosity coefficient
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(inout), device`) : viscous flux components in z direction
  - `u(threadsGv%z+1,0:threadsGv%x+1,threadsGv%y)` (`real(8), shared`)
  - `v(threadsGv%z+1,0:threadsGv%y+1,threadsGv%x)` (`real(8), shared`)
  - `w(threadsGv%z+1,0:threadsGv%x+1,0:threadsGv%y+1)` (`real(8), shared`)
  - `tzx, tzy, tzz, utzx, vtzy, wtzz, kTz` (`real(8)`)

### `subroutine calc_Gv_LES2(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)` (`real(8), intent(in), device`)
  - `mut(nx,ny,nz), qc2(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(inout), device`)
  - `tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs` (`real(8)`)
  - `mx, mxsgs, my, mysgs` (`real(8), dimension(2), device`)
  - `mx(2)` (`real(8), device`)
  - `my(2)` (`real(8), device`)
  - `H(2)` (`real(8), device`)

## File: `calc_flux_base.f90`

### `module calc_flux_base`
- **Description:** Groups all GPU kernel calls for computing E, F, G flux components

### `module procedure`

### `module procedure`

### `subroutine calc_conv_keep(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)`
- **Description:** Provides 4th-5th order accuracy by minimizing dispersive errors in smooth regions
- **Arguments & Variables:**
  - `id_scheme` (`integer(2), intent(in), value`) : ID for scheme: int 2 means KEEP
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `inv_dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing x (1/dx)
  - `inv_dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing y (1/dy)
  - `inv_dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature field
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : convective flux in x direction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : convective flux in y direction
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : convective flux in z direction

### `subroutine calc_conv_slau(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)`
- **Description:** Dissipation modulated by Ducros sensor: f_d controls blend ratio
- **Arguments & Variables:**
  - `id_scheme` (`real(2), intent(in), value`) : ID for scheme: real 2 means SLAU
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `inv_dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing x (1/dx)
  - `inv_dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing y (1/dy)
  - `inv_dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature field
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : convective flux in x direction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : convective flux in y direction
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : convective flux in z direction
  - `sensor(nx,ny,nz)` (`real(8), device`)

### `subroutine calc_conv_roe(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)`
- **Description:** Entropy fix via sensor prevents expansion shocks at sonic points
- **Arguments & Variables:**
  - `id_scheme` (`real(4), intent(in), value`) : ID for scheme: real 4 means Roe
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `inv_dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing x (1/dx)
  - `inv_dy(ny-1)` (`real(8), intent(in), device`) : inverse grid spacing y (1/dy)
  - `inv_dz(nz-1)` (`real(8), intent(in), device`) : inverse grid spacing z (1/dz)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : conservative variables
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : temperature field
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : convective flux in x direction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : convective flux in y direction
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : convective flux in z direction
  - `sensor(nx,ny,nz)` (`real(8), device`)

### `subroutine calc_conv_hybrid(id_scheme, nx, ny, nz, inv_dx, inv_dy, inv_dz, Q, T, E, F, G)`
- **Description:** Ducros shock sensor: f_d = (∇·u)²/[(∇·u)² + (∇×u)² + ε]
- **Arguments & Variables:**
  - `id_scheme` (`real(8), intent(in), value`) : ID for scheme: real 8 means Hybrid
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `inv_dx(nx-1)` (`real(8), intent(in), device`) : inverse grid spacing x (1/dx)
  - `sensor(nx,ny,nz)` (`real(8), device`)

### `subroutine calc_EFG_Euler(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)`

### `subroutine calc_EFG_visc(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)`
- **Description:** and heat flux via Fourier's law: q = -k*dT/dx where k depends on Prandtl number

### `subroutine calc_EFG_LES(id_visc, nx, ny, nz, inv_dx, inv_dy, inv_dz, Jacobian, QJ, Q, T, mu, mut, qc2, E, F, G)`
- **Description:** Filters out subgrid scales: nu_t = (C_s * Delta)^2 * |S_ij| where Delta is grid filter width

## File: `calc_slau_kernel.f90`

### `module calc_slau_kernel`

### `module procedure`

### `module procedure`

### `module procedure`

### `module procedure`

### `subroutine calc_slau_x6(id_accuracy, nx, ny, nz, Q, sensor, E)`
- **Description:** CUDA Fortran kernel for 6 points SLAU scheme in x direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(8), intent(in), value`) : ID for accuracy, 8 means 6 ppoints
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : Flux in x direction
  - `sx  = threadsE%x + 5` (`integer, parameter`) : tile size in x direction to calc high-order interpolation
  - `sxr = threadsE%x` (`integer, parameter`) : tile size in x direction to store the results
  - `sy  = threadsE%y` (`integer, parameter`) : tile size in y direction
  - `sz  = threadsE%z` (`integer, parameter`) : tile size in z direction
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:sx*sy*sz-2), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension(sxr*sy*sz), shared`)

### `subroutine calc_slau_y6(id_accuracy, nx, ny, nz, Q, sensor, F)`
- **Description:** CUDA Fortran kernel for 6 points SLAU scheme in y direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(8), intent(in), value`) : ID for accuracy, 8 means 6 ppoints
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : Flux in y direction
  - `sx  = threadsF%x` (`integer, parameter`) : tile size in x direction
  - `sy  = threadsF%y + 5` (`integer, parameter`) : tile size in y direction to calc high-order interpolation
  - `syr = threadsF%y` (`integer, parameter`) : tile size in y direction to store the results
  - `sz  = threadsF%z` (`integer, parameter`) : tile size in z direction
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:sx*sy*sz-2), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension(sx*syr*sz), shared`)

### `subroutine calc_slau_z6(id_accuracy, nx, ny, nz, Q, sensor, G)`
- **Description:** CUDA Fortran kernel for 6 points SLAU scheme in z direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(8), intent(in), value`) : ID for accuracy, 8 means 6 ppoints
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : Flux in z direction
  - `sx  = threadsG%x` (`integer, parameter`) : tile size in x direction
  - `sy  = threadsG%y` (`integer, parameter`) : tile size in y direction
  - `sz  = threadsG%z + 5` (`integer, parameter`) : tile size in z direction to calc high-order interpolation
  - `szr = threadsG%z` (`integer, parameter`) : tile size in z direction to store the results
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:sx*sy*sz-2), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension(sx*sy*szr), shared`)

### `subroutine calc_slau_x4(id_accuracy, nx, ny, nz, Q, sensor, E)`
- **Description:** CUDA Fortran kernel for 4 points SLAU scheme in x direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(4), intent(in), value`) : ID for accuracy, 4 means 4 ppoints
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : Flux in x direction
  - `sx  = threadsE%x + 3` (`integer, parameter`) : tile size in x direction to calc high-order interpolation
  - `sxr = threadsE%x` (`integer, parameter`) : tile size in x direction to store the results
  - `sy  = threadsE%y` (`integer, parameter`) : tile size in y direction
  - `sz  = threadsE%z` (`integer, parameter`) : tile size in z direction
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:sx*sy*sz-1), shared`) : smem to calc high-order interpolation
  - `rhor, ur, vr, wr, pr` (`real(8), dimension(sxr*sy*sz), shared`) : smem to store the results

### `subroutine calc_slau_y4(id_accuracy, nx, ny, nz, Q, sensor, F)`
- **Description:** CUDA Fortran kernel for 4 points SLAU scheme in y direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(4), intent(in), value`) : ID for accuracy, 4 means 4 ppoints
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : Flux in y direction
  - `sx  = threadsF%x` (`integer, parameter`) : tile size in x direction
  - `sy  = threadsF%y + 3` (`integer, parameter`) : tile size in y direction to calc high-order interpolation
  - `syr = threadsF%y` (`integer, parameter`) : tile size in y direction to store the results
  - `sz  = threadsF%z` (`integer, parameter`) : tile size in z direction
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:sx*sy*sz-1), shared`) : smem to calc high-order interpolation
  - `rhor, ur, vr, wr, pr` (`real(8), dimension(sx*syr*sz), shared`) : smem to store the results

### `subroutine calc_slau_z4(id_accuracy, nx, ny, nz, Q, sensor, G)`
- **Description:** CUDA Fortran kernel for 4 points SLAU scheme in z direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(4), intent(in), value`) : ID for accuracy, 4 means 4 ppoints
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : Flux in z direction
  - `sx  = threadsG%x` (`integer, parameter`) : tile size in x direction
  - `sy  = threadsG%y` (`integer, parameter`) : tile size in y direction
  - `sz  = threadsG%z + 3` (`integer, parameter`) : tile size in z direction to calc high-order interpolation
  - `szr = threadsG%z` (`integer, parameter`) : tile size in z direction to store the results
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:sx*sy*sz-1), shared`) : smem to calc high-order interpolation
  - `rhor, ur, vr, wr, pr` (`real(8), dimension(sx*sy*szr), shared`) : smem to store the results

### `subroutine calc_slau_x2(id_accuracy, nx, ny, nz, Q, sensor, E)`
- **Description:** CUDA Fortran kernel for 2 points SLAU scheme in x direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(2), intent(in), value`) : ID for accuracy, 2 means 2nd-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : Flux in x direction
  - `sx = threadsE%x+1` (`integer, parameter`) : tile size in x direction
  - `sy = threadsE%y` (`integer, parameter`) : tile size in y direction
  - `sz = threadsE%z` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p` (`real(8), dimension(sx*sy*sz), shared`)

### `subroutine calc_slau_y2(id_accuracy, nx, ny, nz, Q, sensor, F)`
- **Description:** CUDA Fortran kernel for 2 points SLAU scheme in y direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(2), intent(in), value`) : ID for accuracy, 2 means 2nd-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : Flux in y direction
  - `sx = threadsF%x` (`integer, parameter`) : tile size in x direction
  - `sy = threadsF%y+1` (`integer, parameter`) : tile size in y direction
  - `sz = threadsF%z` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p` (`real(8), dimension(sx*sy*sz), shared`)

### `subroutine calc_slau_z2(id_accuracy, nx, ny, nz, Q, sensor, G)`
- **Description:** CUDA Fortran kernel for 2 points SLAU scheme in z direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(2), intent(in), value`) : ID for accuracy, 2 means 2nd-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `sensor(nx,ny,nz)` (`real(8), intent(in), device`) : shock sensor
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : Flux in z direction
  - `sx = threadsG%x` (`integer, parameter`) : tile size in x direction
  - `sy = threadsG%y` (`integer, parameter`) : tile size in y direction
  - `sz = threadsG%z+1` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p` (`real(8), dimension(sx*sy*sz), shared`)

## File: `calc_hybrid_kernel.f90`

### `module calc_hybrid_kernel`

### `module procedure`

### `module procedure`

### `module procedure`

### `module procedure`

### `subroutine calc_hybrid_x6(id_accuracy, nx, ny, nz, Q, T, sensor, E)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:threadsE%x+3,threadsE%y,threadsE%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsE%x, threadsE%y,threadsE%z), shared`)

### `subroutine calc_hybrid_y6(id_accuracy, nx, ny, nz, Q, T, sensor, F)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:threadsF%y+3,threadsF%x,threadsF%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsF%y, threadsF%x,threadsF%z), shared`)

### `subroutine calc_hybrid_z6(id_accuracy, nx, ny, nz, Q, T, sensor, G)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:threadsG%z+3,threadsG%y,threadsG%x), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsG%z, threadsG%y,threadsG%x), shared`)

### `subroutine calc_hybrid_x4(id_accuracy, nx, ny, nz, Q, T, sensor, E)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsE%x, threadsE%y,threadsE%z), shared`)

### `subroutine calc_hybrid_y4(id_accuracy, nx, ny, nz, Q, T, sensor, F)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsF%y+2,threadsF%x,threadsF%z), shared`)

### `subroutine calc_hybrid_z4(id_accuracy, nx, ny, nz, Q, T, sensor, G)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsG%z+2,threadsG%y,threadsG%x), shared`)

### `subroutine calc_hybrid_x2(id_accuracy, nx, ny, nz, Q, T, sensor, E)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`)
  - `rho, u, v, w, p` (`real(8), dimension(threadsE%x+1,threadsE%y,threadsE%z), shared`)

### `subroutine calc_hybrid_y2(id_accuracy, nx, ny, nz, Q, T, sensor, F)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`)
  - `rho, u, v, w, p` (`real(8), dimension(threadsF%y+1,threadsF%x,threadsF%z), shared`)

### `subroutine calc_hybrid_z2(id_accuracy, nx, ny, nz, Q, T, sensor, G)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`)
  - `rho, u, v, w, p` (`real(8), dimension(threadsG%z+1,threadsG%y,threadsG%x), shared`)

## File: `preprocess.f90`

### `module preprocess`
- **Description:** Handles device setup, memory allocation, and initial data transfers

### `subroutine check_gpu(mygpu)`
- **Description:** Prints device name and capability information
- **Arguments & Variables:**
  - `mygpu` (`integer, intent(in)`) : GPU device ID to check

### `subroutine allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)`
- **Description:** Size and allocation depends on viscosity model selection
- **Arguments & Variables:**
  - `myrank` (`integer, intent(in)`) : MPI rank
  - `nx` (`integer, intent(in)`) : x grid dimension
  - `ny` (`integer, intent(in)`) : y grid dimension
  - `nz` (`integer, intent(in)`) : z grid dimension
  - `dx(:)` (`real(8), intent(out), allocatable, device`) : inverse grid spacing x
  - `dy(:)` (`real(8), intent(out), allocatable, device`) : inverse grid spacing y
  - `dz(:)` (`real(8), intent(out), allocatable, device`) : inverse grid spacing z
  - `xix(:)` (`real(8), intent(out), allocatable, device`) : coordinate transform metric in x
  - `etay(:)` (`real(8), intent(out), allocatable, device`) : coordinate transform metric in y
  - `zetaz(:)` (`real(8), intent(out), allocatable, device`) : coordinate transform metric in z
  - `Jacobian(:,:)` (`real(8), intent(out), allocatable, device`) : Jacobian determinant for coordinate transform
  - `ruvwp(:,:,:,:)` (`real(8), intent(out), allocatable, device`) : work array for momentum/velocities
  - `T(:,:,:)` (`real(8), intent(out), allocatable, device`) : temperature field
  - `mu(:,:,:)` (`real(8), intent(out), allocatable, device`) : molecular viscosity
  - `mut(:,:,:)` (`real(8), intent(out), allocatable, device`) : turbulent viscosity (LES)
  - `qc2(:,:,:)` (`real(8), intent(out), allocatable, device`) : quadratic constitutive terms
  - `E(:,:,:,:)` (`real(8), intent(out), allocatable, device`) : flux in x direction
  - `F(:,:,:,:)` (`real(8), intent(out), allocatable, device`) : flux in y direction
  - `G(:,:,:,:)` (`real(8), intent(out), allocatable, device`) : flux in z direction

### `subroutine pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap,                        dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)`
- **Description:** Divides computational domain across MPI ranks
- **Arguments & Variables:**
  - `nx` (`integer, intent(in)`) : x grid dimension
  - `ny` (`integer, intent(in)`) : y grid dimension
  - `nz` (`integer, intent(in)`) : z grid dimension
  - `myrank` (`integer, intent(in)`) : MPI rank of this process
  - `nranks` (`integer, intent(in)`) : total number of MPI ranks
  - `x(nx)` (`real(8), intent(in)`) : x coordinate array (host)
  - `dx_cpu(nx-1)` (`real(8), intent(in)`) : inverse x spacing (host)
  - `y(ny)` (`real(8), intent(in)`) : y coordinate array (host)
  - `dy_cpu(ny-1)` (`real(8), intent(in)`) : inverse y spacing (host)
  - `z(nz)` (`real(8), intent(in)`) : z coordinate array (host)
  - `dz_cpu(nz-1)` (`real(8), intent(in)`) : inverse z spacing (host)
  - `Jacobian_cpu(nx,ny)` (`real(8), intent(in)`) : Jacobian determinant (host)
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout)`) : conservative variables on host
  - `overlap` (`integer, intent(out)`) : ghost cell width for MPI halo exchange
  - `dx(nx-1)` (`real(8), intent(out), device`) : inverse x spacing (device)
  - `dy(ny-1)` (`real(8), intent(out), device`) : inverse y spacing (device)
  - `dz(nz-1)` (`real(8), intent(out), device`) : inverse z spacing (device)
  - `xix(nx-1)` (`real(8), intent(out), device`) : x coordinate metric (device)
  - `etay(ny-1)` (`real(8), intent(out), device`) : y coordinate metric (device)
  - `zetaz(nz-1)` (`real(8), intent(out), device`) : z coordinate metric (device)
  - `Jacobian(nx,ny)` (`real(8), intent(out), device`) : Jacobian determinant (device)
  - `QJ(5,nx,ny,nz)` (`real(8), intent(out), device`) : Q divided by Jacobian (device)
  - `ke0` (`real(4), intent(inout)`) : reference kinetic energy
  - `entropy0` (`real(4), intent(inout)`) : reference entropy

### `subroutine pre_rescale(myrank, flag_re, flag_req, ny, nz, Qre, Qm, Qm_cpu)`
- **Arguments & Variables:**
  - `myrank, ny, nz` (`integer, intent(in)`)
  - `flag_re, flag_req` (`integer, intent(inout)`)
  - `Qre(:), Qm(:)` (`real(8), intent(inout), allocatable, device`)
  - `Qm_cpu(:)` (`real(8), intent(inout), allocatable`)

## File: `calc_roe_kernel.f90`

### `module calc_roe_kernel`

### `module procedure`

### `module procedure`

### `module procedure`

### `subroutine calc_roe_x6(id_accuracy, nx, ny, nz, Q, sensor, E)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:threadsE%x+3,threadsE%y,threadsE%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsE%x, threadsE%y,threadsE%z), shared`)

### `subroutine calc_roe_y6(id_accuracy, nx, ny, nz, Q, sensor, F)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:threadsF%y+3,threadsF%x,threadsF%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsF%y, threadsF%x,threadsF%z), shared`)

### `subroutine calc_roe_z6(id_accuracy, nx, ny, nz, Q, sensor, G)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=8), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(-1:threadsG%z+3,threadsG%y,threadsG%x), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsG%z, threadsG%y,threadsG%x), shared`)

### `subroutine calc_roe_x4(id_accuracy, nx, ny, nz, Q, sensor, E)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsE%x, threadsE%y,threadsE%z), shared`)

### `subroutine calc_roe_y4(id_accuracy, nx, ny, nz, Q, sensor, F)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsF%y+2,threadsF%x,threadsF%z), shared`)

### `subroutine calc_roe_z4(id_accuracy, nx, ny, nz, Q, sensor, G)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=4), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`)
  - `rho,  u,  v,  w,  p` (`real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared`)
  - `rhor, ur, vr, wr, pr` (`real(8), dimension( threadsG%z+2,threadsG%y,threadsG%x), shared`)

### `subroutine calc_roe_x2(id_accuracy, nx, ny, nz, Q, sensor, E)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`)
  - `rho, u, v, w, p` (`real(8), dimension(threadsE%x+1,threadsE%y,threadsE%z), shared`)

### `subroutine calc_roe_y2(id_accuracy, nx, ny, nz, Q, sensor, F)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`)
  - `rho, u, v, w, p` (`real(8), dimension(threadsF%y+1,threadsF%x,threadsF%z), shared`)

### `subroutine calc_roe_z2(id_accuracy, nx, ny, nz, Q, sensor, G)`
- **Arguments & Variables:**
  - `id_accuracy` (`integer(kind=2), intent(in), value`)
  - `nx, ny, nz` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz), sensor(nx,ny,nz)` (`real(8), intent(in), device`)
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`)
  - `rho, u, v, w, p` (`real(8), dimension(threadsG%z+1,threadsG%y,threadsG%x), shared`)

## File: `set_init_common.f90`

### `module set_init_common`

### `subroutine calc_Gaussian_filter_x(nx, ny, nz, n, x, phi)`
- **Arguments & Variables:**
  - `nx, ny, nz, n` (`integer, intent(in)`)
  - `x(nx)` (`real(8), intent(in)`)
  - `phi(nx,ny,nz)` (`real(8), intent(inout)`)
  - `i, j, k, ix` (`integer`)
  - `phi_tmp(:,:,:), wg(:)` (`real(8), allocatable`)

### `subroutine calc_Gaussian_filter_y(nx, ny, nz, n, y, phi)`
- **Arguments & Variables:**
  - `nx, ny, nz, n` (`integer, intent(in)`)
  - `y(ny)` (`real(8), intent(in)`)
  - `phi(nx,ny,nz)` (`real(8), intent(inout)`)
  - `i, j, k, iy` (`integer`)
  - `phi_tmp(:,:,:), wg(:)` (`real(8), allocatable`)

### `subroutine calc_Gaussian_filter_z(nx, ny, nz, n, z, phi)`
- **Arguments & Variables:**
  - `nx, ny, nz, n` (`integer, intent(in)`)
  - `z(nz)` (`real(8), intent(in)`)
  - `phi(nx,ny,nz)` (`real(8), intent(inout)`)
  - `i, j, k, iz` (`integer`)
  - `phi_tmp(:,:,:), wg(:)` (`real(8), allocatable`)

### `subroutine set_bc_cyclic_x_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in)`)
  - `ustd(nx,ny,nz), vstd(nx,ny,nz), wstd(nx,ny,nz), Tstd(nx,ny,nz)` (`real(8), intent(inout)`)

### `subroutine set_bc_cyclic_z_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in)`)
  - `ustd(nx,ny,nz), vstd(nx,ny,nz), wstd(nx,ny,nz), Tstd(nx,ny,nz)` (`real(8), intent(inout)`)

### `subroutine calc_rms(nx, ny, nz, phi, rms)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in)`)
  - `phi(nx,ny,nz)` (`real(8), intent(in)`)
  - `rms` (`real(8), intent(out)`)
  - `phi2(:,:,:)` (`real(8), allocatable`)

### `subroutine set_init_tbl(nx, ny, nz, x, y, z, rand, blt0, blt, u0, p0, T0, M0, Q)`
- **Arguments & Variables:**
  - `nx, ny, nz` (`integer, intent(in)`)
  - `x(nx), y(ny), z(nz)` (`real(8), intent(in)`)
  - `rand, blt0, blt, u0, p0, T0, M0` (`real(8), intent(in)`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(out)`)
  - `seed(:)` (`integer(4), allocatable`)
  - `p_wall` (`real(8)`)
  - `fd, pi = acos(-1.d0)` (`real(8)`)
  - `rho(:), u(:), v(:), T(:), randum(:,:,:,:), ustd(:,:,:), vstd(:,:,:), wstd(:,:,:), Tstd(:,:,:)` (`real(8), allocatable`)

## File: `calc_para.f90`

### `module calc_para`
- **Description:** Handles halo exchange, data flattening, and ghost cell synchronization

### `module procedure`

### `subroutine flatten(nx, ny, nz, overlap, Q, Q1d_left, Q1d_right)`
- **Description:** Used for MPI halo exchange preparation
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : x dimension
  - `ny` (`integer, intent(in), value`) : y dimension
  - `nz` (`integer, intent(in), value`) : z dimension
  - `overlap` (`integer, intent(in), value`) : ghost cell width
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : 3D conservative variables
  - `Q1d_left(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(out), device`) : left boundary 1D array
  - `Q1d_right(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(out), device`) : right boundary 1D array

### `subroutine flatten_left(nx, ny, nz, overlap, Q, Q1d_left)`
- **Description:** Flatten 3D Q data into 1D array for left boundary only
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : x dimension
  - `ny` (`integer, intent(in), value`) : y dimension
  - `nz` (`integer, intent(in), value`) : z dimension
  - `overlap` (`integer, intent(in), value`) : ghost cell width
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : 3D conservative variables
  - `Q1d_left(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(out), device`) : left boundary 1D array

### `subroutine flatten_right(nx, ny, nz, overlap, Q, Q1d_right)`
- **Description:** Flatten 3D Q data into 1D array for right boundary only
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : x dimension
  - `ny` (`integer, intent(in), value`) : y dimension
  - `nz` (`integer, intent(in), value`) : z dimension
  - `overlap` (`integer, intent(in), value`) : ghost cell width
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : 3D conservative variables
  - `Q1d_right(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(out), device`) : right boundary 1D array

### `subroutine flatten_rescale(nx, ny, nz, nre, overlap, Q, Q1d_right)`
- **Arguments & Variables:**
  - `nx, ny, nz, nre, overlap` (`integer, intent(in), value`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`)
  - `Q1d_right(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(out), device`)

### `subroutine reconstruct(nx, ny, nz, overlap, Q1d_left, Q1d_right, Q)`
- **Arguments & Variables:**
  - `nx, ny, nz, overlap` (`integer, intent(in), value`)
  - `Q1d_left(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(in), device`)
  - `Q1d_right(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(in), device`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(out), device`)

### `subroutine reconstruct_left(nx, ny, nz, overlap, Q1d_left, Q)`
- **Arguments & Variables:**
  - `nx, ny, nz, overlap` (`integer, intent(in), value`)
  - `Q1d_left(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(in), device`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(out), device`)

### `subroutine reconstruct_right(nx, ny, nz, overlap, Q1d_right, Q)`
- **Arguments & Variables:**
  - `nx, ny, nz, overlap` (`integer, intent(in), value`)
  - `Q1d_right(overlap*(ny-2)*(nz-6)*5)` (`real(8), intent(in), device`)
  - `Q(5,nx,ny,nz)` (`real(8), intent(out), device`)

### `subroutine reconstruct_sbli_inlet(nx, ny1, ny2, nz, overlap, Q1d, Q)`
- **Arguments & Variables:**
  - `nx, ny1, ny2, nz, overlap` (`integer, intent(in), value`)
  - `Q1d(overlap*(ny1-2)*(nz-6)*5)` (`real(8), intent(in), device`)
  - `Q(5,nx,ny2,nz)` (`real(8), intent(inout), device`)

### `subroutine exchange_cyclic(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)`
- **Arguments & Variables:**
  - `id_rescale` (`integer(kind=2), intent(in), value`)
  - `myrank, nranks, overlap, nx, ny, nz` (`integer, intent(in), value`)
  - `Qs_left,   Qs_right,   Qr_left,   Qr_right` (`real(8), dimension(overlap*(ny-2)*(nz-6)*5)`)
  - `Qs1d_left, Qs1d_right, Qr1d_left, Qr1d_right` (`real(8), dimension(overlap*(ny-2)*(nz-6)*5), device`)

### `subroutine exchange_rescale(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)`
- **Arguments & Variables:**
  - `id_rescale` (`integer(kind=4), intent(in), value`)
  - `myrank, nranks, overlap, nx, ny, nz` (`integer, intent(in), value`)
  - `Qs_left,   Qs_right,   Qr_left,   Qr_right` (`real(8), dimension(overlap*(ny-2)*(nz-6)*5)`)
  - `Qs1d_left, Qs1d_right, Qr1d_left, Qr1d_right` (`real(8), dimension(overlap*(ny-2)*(nz-6)*5), device`)

## File: `set_bc_tbl_sbli.f90`

### `module set_bc_tbl_sbli`

### `subroutine set_bc_Neumann_tbl_top_down(nx, ny, nz, offset, istart, iend, Jacobian, QJ)`
- **Arguments & Variables:**
  - `nx, ny, nz, offset, istart, iend` (`integer, intent(in), value`)
  - `Jacobian(nx,ny)` (`real(8), intent(in), device`)
  - `QJ(5,nx,ny,nz)` (`real(8), intent(inout), device`)
  - `Jacobian_tmp, over_QJ1, p_wall, rhob, ub, vb, wb, pb` (`real(8)`)

### `subroutine set_bc_Riemann_tbl_top_down(nx, ny, nz, offset, istart, iend, Jacobian, QJ)`
- **Arguments & Variables:**
  - `nx, ny, nz, offset, istart, iend` (`integer, intent(in), value`)
  - `Jacobian(nx,ny)` (`real(8), intent(in), device`)
  - `QJ(5,nx,ny,nz)` (`real(8), intent(inout), device`)
  - `T       = Taw - rf * u0**2 / (2.d0 * Cp)` (`real(8), parameter`)
  - `rho0    = p0 / (R * T)` (`real(8), parameter`)
  - `c0      = sqrt(gamma * p0 / rho0)` (`real(8), parameter`)
  - `over_c0 = 1.d0 / c0` (`real(8), parameter`)
  - `Jacobian_tmp, p_wall, rhoin, pin, cin, vin, Rp, Rm, rhob, vb, cb, pb, v0 = 0.d0` (`real(8)`)

## File: `calc_steps.f90`

### `module calc_steps`

### `subroutine calc_step1(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2)`
- **Description:** CUDA Fortran kernel for 1st step of 3-3 TVD Runge-Kutta
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `coef` (`real(8), intent(in), value`) : coefficient for Runge-Kutta
  - `dx(nx-1)` (`real(8), intent(in), device`) : grid size in x direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : grid size in y direction
  - `dz(nz-1)` (`real(8), intent(in), device`) : grid size in z direction
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(in), device`) : Flux in x direction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(in), device`) : Flux in y direction
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(in), device`) : Flux in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : present Q(rho, rhou, rhov, rhow, E) / Jacobian
  - `Q2(5,nx,ny,nz)` (`real(8), intent(out), device`) : next    Q(rho, rhou, rhov, rhow, E) / Jacobian

### `subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)`
- **Description:** CUDA Fortran kernel for 1st~3rd step of 4-4 Runge-Kutta
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `coef1` (`real(8), intent(in), value`) : coefficient for Runge-Kutta
  - `coef2` (`real(8), intent(in), value`) : coefficient for Runge-Kutta
  - `dx(nx-1)` (`real(8), intent(in), device`) : grid size in x direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : grid size in y direction
  - `dz(nz-1)` (`real(8), intent(in), device`) : grid size in z direction
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(in), device`) : Flux in x direction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(in), device`) : Flux in y direction
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(in), device`) : Flux in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : present Q(rho, rhou, rhov, rhow, E) / Jacobian
  - `Q2(5,nx,ny,nz)` (`real(8), intent(out), device`) : next    Q(rho, rhou, rhov, rhow, E) / Jacobian
  - `Rs(5,nx-2,ny-2,nz-2)` (`real(8), intent(inout), device`) : accumulation for 4-4 Runge-Kutta

### `subroutine calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Qin, Qout)`
- **Description:** TVD RK3 Stage 2 & 3: Q^(n+1) = (α*Q^n + β*Q^(*) - γ*R)/(α+β)
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `coef1` (`real(8), intent(in), value`) : α coefficient (weight of original Q^n)
  - `coef2` (`real(8), intent(in), value`) : β coefficient (weight of Q^(*))
  - `coef3` (`real(8), intent(in), value`) : γ coefficient (weight of flux residual)
  - `coef4` (`real(8), intent(in), value`) : 1/(α+β) normalization factor
  - `dx(nx-1)` (`real(8), intent(in), device`) : grid size in x direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : grid size in y direction
  - `dz(nz-1)` (`real(8), intent(in), device`) : grid size in z direction
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(in), device`) : Flux in x direction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(in), device`) : Flux in y direction
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(in), device`) : Flux in z direction
  - `Qin(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q^n (original from previous step)
  - `Qout(5,nx,ny,nz)` (`real(8), intent(inout), device`) : Q^(*) on input, Q^(n+1) on output

### `subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)`
- **Description:** Final RK4 Stage: Q^n+1 = Q^n - (1/6)·∑(R_ᵢ) where R_ᵢ indexed over 4 stages
- **Arguments & Variables:**
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `dx(nx-1)` (`real(8), intent(in), device`) : grid size in x direction
  - `dy(ny-1)` (`real(8), intent(in), device`) : grid size in y direction
  - `dz(nz-1)` (`real(8), intent(in), device`) : grid size in z direction
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(in), device`) : Flux in x direction
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(in), device`) : Flux in y direction
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(in), device`) : Flux in z direction
  - `Rs(5,nx-2,ny-2,nz-2)` (`real(8), intent(inout), device`) : accumulated residuals from stages 1-3
  - `Q(5,nx,ny,nz)` (`real(8), intent(inout), device`) : Q^n on input, Q^n+1 on output

## File: `calc_keep_kernel.f90`

### `module calc_keep_kernel`

### `module procedure`

### `module procedure`

### `module procedure`

### `subroutine calc_keep_x6(id_accuracy, nx, ny, nz, Q, T, E)`
- **Description:** CUDA Fortran kernel for 6th-order KEEP scheme in x direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(8), intent(in), value`) : ID for accuracy, 8 means 6th-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : Flux in x direction
  - `sx = threadsE%x + 5` (`integer, parameter`) : tile size in x direction
  - `sy = threadsE%y` (`integer, parameter`) : tile size in y direction
  - `sz = threadsE%z` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(-1:sx*sy*sz-2), shared`)

### `subroutine calc_keep_y6(id_accuracy, nx, ny, nz, Q, T, F)`
- **Description:** CUDA Fortran kernel for 6th-order KEEP scheme in y direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(8), intent(in), value`) : ID for accuracy, 8 means 6th-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : Flux in y direction
  - `sx = threadsF%x` (`integer, parameter`) : tile size in x direction
  - `sy = threadsF%y + 5` (`integer, parameter`) : tile size in y direction
  - `sz = threadsF%z` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(-1:sx*sy*sz-2), shared`)

### `subroutine calc_keep_z6(id_accuracy, nx, ny, nz, Q, T, G)`
- **Description:** CUDA Fortran kernel for 6th-order KEEP scheme in z direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(8), intent(in), value`) : ID for accuracy, 8 means 6th-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : Flux in z direction
  - `sx = threadsG%x` (`integer, parameter`) : tile size in x direction
  - `sy = threadsG%y` (`integer, parameter`) : tile size in y direction
  - `sz = threadsG%z + 5` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(-1:sx*sy*sz-2), shared`)

### `subroutine calc_keep_x4(id_accuracy, nx, ny, nz, Q, T, E)`
- **Description:** CUDA Fortran kernel for 4th-order KEEP scheme in x direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(4), intent(in), value`) : ID for accuracy, 4 means 4th-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : Flux in x direction
  - `sx = threadsE%x + 3` (`integer, parameter`) : tile size in x direction
  - `sy = threadsE%y` (`integer, parameter`) : tile size in y direction
  - `sz = threadsE%z` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(0:sx*sy*sz-1), shared`)

### `subroutine calc_keep_y4(id_accuracy, nx, ny, nz, Q, T, F)`
- **Description:** CUDA Fortran kernel for 4th-order KEEP scheme in y direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(4), intent(in), value`) : ID for accuracy, 4 means 4th-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : Flux in y direction
  - `sx = threadsF%x` (`integer, parameter`) : tile size in x direction
  - `sy = threadsF%y + 3` (`integer, parameter`) : tile size in y direction
  - `sz = threadsF%z` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(0:sx*sy*sz-1), shared`)

### `subroutine calc_keep_z4(id_accuracy, nx, ny, nz, Q, T, G)`
- **Description:** CUDA Fortran kernel for 4th-order KEEP scheme in z direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(4), intent(in), value`) : ID for accuracy, 4 means 4th-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : Flux in z direction
  - `sx = threadsG%x` (`integer, parameter`) : tile size in x direction
  - `sy = threadsG%y` (`integer, parameter`) : tile size in y direction
  - `sz = threadsG%z + 3` (`integer, parameter`) : tile size in z direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(0:sx*sy*sz-1), shared`)

### `subroutine calc_keep_x2(id_accuracy, nx, ny, nz, Q, T, E)`
- **Description:** CUDA Fortran kernel for 2nd-order KEEP scheme in x direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(2), intent(in), value`) : ID for accuracy, 2 means 2nd-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `E(5,nx-1,ny-2,nz-2)` (`real(8), intent(out), device`) : Flux in x direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(2), device`)

### `subroutine calc_keep_y2(id_accuracy, nx, ny, nz, Q, T, F)`
- **Description:** CUDA Fortran kernel for 2nd-order KEEP scheme in y direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(2), intent(in), value`) : ID for accuracy, 2 means 2nd-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `F(5,nx-2,ny-1,nz-2)` (`real(8), intent(out), device`) : Flux in y direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(2), device`)

### `subroutine calc_keep_z2(id_accuracy, nx, ny, nz, Q, T, G)`
- **Description:** CUDA Fortran kernel for 2nd-order KEEP scheme in z direction
- **Arguments & Variables:**
  - `id_accuracy` (`integer(2), intent(in), value`) : ID for accuracy, 2 means 2nd-order
  - `nx` (`integer, intent(in), value`) : number of grid points in x direction
  - `ny` (`integer, intent(in), value`) : number of grid points in y direction
  - `nz` (`integer, intent(in), value`) : number of grid points in z direction
  - `Q(5,nx,ny,nz)` (`real(8), intent(in), device`) : Q(rho, u, v, w, p)
  - `T(nx,ny,nz)` (`real(8), intent(in), device`) : Temperature
  - `G(5,nx-2,ny-2,nz-1)` (`real(8), intent(out), device`) : Flux in z direction
  - `rho, u, v, w, p, tmp` (`real(8), dimension(2), device`)

