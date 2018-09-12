# DQMC demonstration code for spin-fermion model

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

```
FORTRAN Compiler such as Intel or Gfortran. Lapack and Blas libraries needed.
```

### Installation

Firstly git clone the project with


> git clone git@gitlab.com:xyxu/slmc.git

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

## Code Introduction

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


### ftdqmc_asvqrd.f90
Numerical stablization using ASvQRD method.

- **ftdqmc_asvqrd_alloc** : Allocate memory for matrices used in numerical stablization process.
- **ftdqmc_asvqrd_free** : Free up space
- **ftdqmc_asvqrd_initst** : Initialize Ust, Dst, Vst and logdetQst
- **ftdqmc_asvqrd_initR** : Initialize UR, DRvec, VR and logdetQR
- **ftdqmc_asvqrd_initL** : Initialize UL, DLvec, VL and logdetQL
- **ftdqmc_stablize_0b_svd** : Numerical stablization during Monte Carlo local update sweep from 0 -> beta
- **ftdqmc_stablize_b0_svd** : Numerical stablization during Monte Carlo local update sweep from beta -> 0

### ftdqmc_gfun.f90
Subroutines associated with Green's function (equal time and time displaced).

- **ftdqmc_gfun_alloc** : Allocate memory.
- **ftdqmc_gfun_free** : Free up space.
- **ftdqmc_gfun_initgt** : Initialize g00, gt0 and g0t.
- **green_equaltime**: Calculate G(tau, tau).
- **green_equaltime00**: Calculate G(0, 0).
- **green_equaltimebb**: Calculate G(beta, beta).
- **green_tau**: Calculate g00, gt0, g0t, gtt.
- **push_stage** : Save UDV and Green's function.
- **pop_stage**: Restore UDV and Green's function.

### ftdqmc_hamilt.f90
Set up global parameters for the model and runtime flags.

- **ftdqmc_hamilt_init**: Read parameters from file ftdqmc.in and broadcast them to other MPI processes.
- **ftdqmc_hamilt_free**: Free up space

### ftdqmc_latt.f90
Set up the lattice.

- **ftdqmc_latt_alloc**: Allocate memory.
- **ftdqmc_latt_free**: Free up space.
- **ftdqmc_latt_sli**: Set up tables such list, invlist, listk and nearest-neighbor tables as well as calculate the distance.
- **ftdqmc_latt_sltpf**: Set up tables for checkboard decomposition.
- **ftdqmc_latt_sthop**: Set up hopping matrix.

### ftdqmc_auxfield.f90
Class for auxiliary field or bosonic field. Currently one abstract superclass **ftdqmc_auxfield** and one subclass **ftdqmc_auxfield_f1** (phonon field)

- **ftdqmc_auxfield_alloc**: Allocate memory.
- **ftdqmc_auxifield_free**: Free up space
- **ftdqmc_auxfield_inconfc**: Read previous configurations from file **confin** or initialize the field.
- **ftdqmc_auxfield_outconfc**: Output configurations to file **confout**.
- **ftdqmc_auxfield_outconfc_bin**: Output configurations and logweight after every sweep for learning the effective model.
- **ftdqmc_auxfield_calc_ediff**: Method that calculate the logweight difference of the bosonic action when doing the local update. This method must be overrided for every subclass.

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


### ftdqmc_phy0.f90
Equal time measurements

- **ftdqmc_phy0_meas**: Scalar quanties and correlation function measurement.
- **ftdqmc_phy0_getavg**: Get average of scalar of one bin and output it.
- **ftdqmc_phy0_corFT**: Fourier transformation of equal time corrlation function.
- **ftdqmc_phy0_ft**: Detailed Fourier transformation.

### ftdqmc_tdm.f90
Time-displaced correlations

- **ftdqmc_tdm_meas**: Time-displaced correlation function measurements.
- **ftdmqc_tdm_corFT**: Fourier transformation of time-displaced correlation function.
- **ftdqmc_tdm_ft**: Detailed Fourier transformation.

### ftdqmc_main.f90
Main program.


### Flags and Variables
#### Flags

- **ltau**: Whether to do time-displaced measurements
- **lstglobal**: Whether to do the global update such as the block update in the Holstein model
- **llocal**: Whether to do the local update
- **lwarnup**: Whether to do the warm up
- **lwrapph**: g_hol > 0 then update the phonon field.

#### Variables

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
- **nsw_stglobal**: When set to  1, do the block update
- **g_hol**: Electron-phonon coupling
- **Omg**: Phonon frequency
- **Mass**: Mass of phonon
- **dx**: Step length of local update of phonon field
- **order**: 0 start from disorder phase, 1 start from randomly disorder checkboard order configurations
- **gcount**: Number of block updates in one sweep
- **phonon_f**: exp(-dtau*g_hol*ran_x)
- **ran_x**: phonon field)



## Running the tests

There are some scripts in run/ that can start a test or submit jobs on server.



## Authors

* **Xiao Yan Xu**   [wanderxu@gmail.com](mailto:wanderxu@gmail.com)
* **Zi Hong Liu**   [zihongliu@iphy.ac.cn](mailto:zihongliu@iphy.ac.cn)
* **Chuang Chen**   [chenchuang@iphy.ac.cn](mailto:chenchuang@iphy.ac.cn)


## Acknowledgments

* Zi Yang Meng
* Gaopei Pan
