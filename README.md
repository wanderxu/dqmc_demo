# DQMC demonstration code for spin-fermion model

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

```
FORTRAN Compiler such as Intel or Gfortran. Lapack and Blas libraries needed.
```

### Installation

Firstly git clone the project with


> git clone git@github.com:wanderxu/dqmc_demo.git

or

> git clone https://github.com/wanderxu/dqmc_demo.git

To compile it, you should choose the right make.sys depends on your compiler. Let's take Gfortran as an example.

Compile lib


> cd lib/
>
> cp make.sys.gfortran make.sys
>
>  make

Compile src

> cd src/
>
> cp make.sys.gfortran make.sys
>
>  make

### Input file format

We use namelist in the input file. The file src/ftdqmc.in is an example.

## Code Introduction

### ftdqmc variables

- **l**: System size in terms of primitive cell
- **lq**: l^2, Number of primitive cells
- **norb**: Number of orbital
- **ndim**: Total number of sites namely lq*norb
- **nfam**: Checkerboard decomposition families of hopping
- **mu**: Chemical potential
- **beta**: Inversed temperature
- **ltrot**: Number of imaginary time slices
- **dtau**:  beta/ltrot
- **xmag**: z-flux
- **flux_x, flux_y**: Twisted boundary condition
- **nwrap**: Frequency of numerical stablization
- **nsweep**: Number of sweeps of one bin
- **nbin**: Number of bins
- **nwarnup**: Number of warm up sweeps
- **nsw_stglobal**: When set to  1, do the global update
- **rt**: nearest neighbor hopping
- **rt2**: next nearest neighbor hopping
- **rt3**: third nearest neighbor hopping
- **js**: exchange coupling of Ising spins
- **hx**: transverse field for Ising spins
- **rhub**: spin-fermion coupling strength

### ftdqmc flags

- **ltau**: Whether to do time-displaced measurements
- **lstglobal**: Whether to do the global update such as the Wolff update or self-learning
- **llocal**: Whether to do the local update
- **lwarnup**: Whether to do the warm up

### ftdqmc set lattice

- **sli.f90**: set lists for lattice
- **sltpf.f90**: set up tables for checkboard decomposition
- **generate_neighbor.f90**: generate neighbor lists

### ftdqmc set hopping matrices

- **sthop.f90**: set hopping matrices
- **thop_mag**: z_flux and twisted boundary conditions

### ftdqmc auxiliary fields related

- **salph.f90**: set auxiliary fields related variables
- **inconfc.f90**: initial auxiliary fields
- **outconfc.f90**: output auxiliary fields

### ftdqmc matrix operation subroutines

Subroutines and modules associated with matrix operations in DQMC algorithm.

- **data_tmp** :  Temporary matrices used in matrix operations
- **mmuur** : Right multiply by exp(V) **V**: Interacting matrix
- **mmuurH**: Right multiply by hermitian of exp(V)
- **muurm1**: Right division by exp(V)
- **mmuul**:  Left multiply by exp(V)
- **mmuulm1**:  Left division by exp(V)
- **mmthr**: Right multiply by exp(-dtau*T) **T**: hopping matrix
- **mmthrH**: Right multiply by hermitian of exp(-datu*T)
- **mmthrm1**: Right division by exp(-dtau*T)
- **mmthl**: Left  multiply by exp(-dtau*T)
- **mmthlm1**: Right division by exp(-dtau*T)

### ftdqmc_core.f90
DQMC sweep

- **ftdqmc_sweep_start_0b**: Sweep from 0 -> beta to prepare UDV.
- **ftdqmc_sweep_start_b0**: Sweep from beta -> 0 to prepare UDV.
- **ftdqmc_sweep_b0**: Sweep from beta -> 0 and doing measurements.
- **ftdqmc_sweep_0b**: Sweep from 0 -> beta and doing measurements.
- **Bmat_tau_R**:  Compute B(tau1, tau2) tau1 > tau2
- **Bmat_tau_RH**: Compute B(tau1, tau2) * tau1 > tau2
- **Bmat_tau_L**: Compute * B(tau1, tau2) tau1 > tau2
- **Bmatinv_tau_L**: Compute * B(tau1, tau2)^-1 tau1 > tau2

Numerical stablization using ASvQRD method

- **ftdqmc_stablize_0b_svd** : Numerical stablization during Monte Carlo local update sweep from 0 -> beta
- **ftdqmc_stablize_b0_svd** : Numerical stablization during Monte Carlo local update sweep from beta -> 0

Subroutines associated with Green's function (equal time and time displaced).

- **green_equaltime**: Calculate G(tau, tau).
- **green_equaltime00**: Calculate G(0, 0).
- **green_equaltimebb**: Calculate G(beta, beta).
- **green_tau**: Calculate g00, gt0, g0t, gtt.

### ftdqmc measurements

- **obser.f90**: perform measurements
- **preq.f90**: mpi reduce equal-time measurements and output
- **prtau.f90**: mpi reduce time-displaced measurements and output

### ftdqmc_main.f90
Main program.



## Running the tests

There are some scripts in example/ that can start a test locally directly or submit jobs on server with a little change according to the server operating system.



## Authors

* **Xiao Yan Xu**   [wanderxu@gmail.com](mailto:wanderxu@gmail.com)
* **Zi Hong Liu**   [zihongliu@iphy.ac.cn](mailto:zihongliu@iphy.ac.cn)
* **Chuang Chen**   [chenchuang@iphy.ac.cn](mailto:chenchuang@iphy.ac.cn)


## Acknowledgments

* Zi Yang Meng
* Gaopei Pan
