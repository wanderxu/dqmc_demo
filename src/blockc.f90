module blockc
  integer, parameter :: dp = 8
  integer, parameter :: fout = 50
  integer, save :: irank, isize, ierr
  real(dp), parameter :: zero = 0.d0
  complex(dp), parameter :: cone = dcmplx( 1.d0, 0.d0 )
  complex(dp), parameter :: czi = dcmplx( 0.d0, 1.d0 )
  complex(dp), parameter :: czero = dcmplx( 0.d0, 0.d0 )
  complex(dp), parameter :: chalf = dcmplx( 0.5d0, 0.d0 )
  complex(dp), parameter :: cquarter = dcmplx( 0.25d0, 0.d0 )
  real(dp), save :: pi = dacos(-1.d0)
 
  ! lattice
  integer, save :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)
  real(dp), save :: b1_p(2), b2_p(2)
  integer, save :: ndim
  integer, save :: l
  integer, save :: lq
  integer, save :: nfam
  integer, save :: lfam
  integer, allocatable, dimension(:,:), save :: list
  integer, allocatable, dimension(:,:), save :: invlist
  integer, allocatable, dimension(:,:), save :: nnlist
  integer, allocatable, dimension(:,:), save :: list_plaq
  integer, allocatable, dimension(:,:), save :: ltpf
  integer, allocatable, dimension(:,:), save :: lthf
  integer, allocatable, dimension(:,:), save :: lthf2 ! next nearest
  integer, allocatable, dimension(:,:), save :: lthf3 ! 3rd nearest

  integer, allocatable, dimension(:), save :: orblist

  integer, allocatable, dimension(:,:), save :: latt_imj
  integer, allocatable, dimension(:,:), save :: listk

  integer, allocatable, dimension(:), save :: imjdeg
  integer, allocatable, dimension(:), save :: distance_index, equ_distance, irre_distance_deg
  real, allocatable, dimension(:), save :: distance_len, irre_distance_len
  integer, save :: num_equ_distance

  ! model
  integer, save :: ne
  real(dp), parameter :: rt = 1.d0
  real(dp), parameter :: rt2 = -0.32d0
  real(dp), parameter :: rt3 = 0.128d0
  real(dp), save :: beta
  real(dp), save :: mu
  real(dp), save :: muA ! chemical potentail for A sublattice
  real(dp), save :: muB ! chemical potential for B sublattice
  real(dp), save :: rhub
  real(dp), save :: rj
  real(dp), save :: js
  real(dp), save :: hx
  real(dp), save :: fill
  real(dp), save :: xmag
  real(dp), save :: dimer
  real(dp), save :: flux_x
  real(dp), save :: flux_y
  complex(dp), save :: alpha_u

  ! tmp 
  real(dp), save :: tanhdth, cothdth, gamma_s



  ! dqmc relative
  logical, save :: lwrapu, lwrapj, llocal, lstglobal
  integer(dp), save :: iseed
  real(dp), save :: dtau
  integer, save :: ltrot
  complex(dp), save :: phase
  real(dp), save :: weight_track
  real(dp), save :: logweightf_old, logweightf_new, logweights_old, logweights_new
  complex(dp), save :: logweightf_up, logweightf_dn
  complex(dp), allocatable, dimension(:,:), save :: Imat
  complex(dp), allocatable, dimension(:,:), save :: grup, grdn, grupc, grdnc
  real(dp), allocatable, dimension(:), save :: Ivec

  ! cal. control
  integer, save :: nbin
  integer, save :: nsweep
  integer, save :: nskip
  integer, save :: nst
  integer, save :: nwrap
  logical, save :: lwarnup
  integer, save :: nwarnup
  integer, save :: n_outconf_pace 
  integer, allocatable, dimension(:) :: iwrap_nt
  integer, allocatable, dimension(:,:) :: wrap_step

  ! for dynamical
  logical, save :: lsstau
  logical, save :: ltau
  integer, save :: nuse
  real(dp), allocatable, dimension(:), save :: wn
  complex(dp), allocatable, dimension(:,:), save :: zexpiwt
  complex(dp), allocatable, dimension(:,:), save :: zexpiqr

  ! for obs
  integer, save :: nobs
  integer, save :: obs_segment_len 
  integer, save :: nmeas_bin
  complex(dp), dimension(10), save :: main_obs, mpi_main_obs


  ! for DQMC
  integer, dimension(-1:1,1), save :: nflipl
  complex(dp), dimension(-1:1), save :: xsigp2, xsigm2
  complex(dp), dimension(-1:1), save :: xsigma_u_up, xsigma_u_dn
  complex(dp), dimension(-1:1,1), save :: delta_u_up, delta_u_dn 
  complex(dp), dimension(-1:1,1), save :: dellp2, dellm2, dgaml, deta_u

  complex(dp), dimension(2,2), save :: ur_k, urt_k, ur_j, urt_j

#ifdef BREAKUP_T
  complex(dp), allocatable, dimension(:,:,:), save :: urt, urtm1
  complex(dp), allocatable, dimension(:,:,:), save :: urtc, urtcm1 ! next nearest
  complex(dp), allocatable, dimension(:,:,:), save :: urtd, urtdm1 ! 3rd nearest
#ifdef SPINDOWN
  complex(dp), allocatable, dimension(:,:,:), save :: urt_dn, urtm1_dn
  complex(dp), allocatable, dimension(:,:,:), save :: urtc_dn, urtcm1_dn
  complex(dp), allocatable, dimension(:,:,:), save :: urtd_dn, urtdm1_dn
#endif
#else
  complex(dp), allocatable, dimension(:,:), save :: urt, urtm1
#ifdef SPINDOWN
  complex(dp), allocatable, dimension(:,:), save :: urt_dn, urtm1_dn
#endif
#endif
  complex(dp), allocatable, dimension(:,:), save :: hopping_tmp
#ifdef SPINDOWN
  complex(dp), allocatable, dimension(:,:), save :: hopping_tmp_dn
#endif
  complex(dp), allocatable, dimension(:), save :: hop_plusx, hop_minusx

  integer, allocatable, dimension(:,:), save :: nsigl_u
  integer, allocatable, dimension(:,:,:), save :: nsigl_k, nsigl_j
  real(dp) :: wsxsz(128), wjs(32)
  real(dp), save :: max_wrap_error, max_wrap_error_tmp 
  real(dp), save :: xmax_dyn, xmax_dyn_tmp

  ! for stglobal update
  integer, save :: icount_nsw_stglobal, nsw_stglobal
  integer, save :: nstcluster
  integer, save :: ntentacle
  integer, save :: ntentacle_old
  integer, save :: num_st_nn
  integer, allocatable, dimension(:,:), save :: tentacle, tentacle_old
  integer, allocatable, dimension(:,:), save :: stcluster
  integer, dimension(:,:,:,:), allocatable, save :: stbonds_neib
  real(dp), save :: ratio_nn_st(6)

  ! delay update
  integer, save :: nublock
  
  contains

  subroutine make_tables
#ifdef MPI
    use mpi
#endif
    implicit none

    !include 'mpif.h'

    integer :: i, nwrap_mid
    logical :: exists

    namelist /model_para/ l, beta, dtau, mu, muA, muB, rhub, rj, js, hx, xmag, flux_x, flux_y
    namelist /ctrl_para/ nwrap, nsweep, nbin, llocal, nsw_stglobal, lsstau, ltau, nuse, nublock

    ! default parameters
    l    = 2
    beta = 20
    dtau = 0.05d0
    mu   = 0.d0 ! default is half filling
    muA   = 0.d0
    muB   = 0.d0
    rhub = 1.0d0
    rj   = 0.d0
    js   = 1.d0
    hx   = 5.d0
    nwrap = 10
    n_outconf_pace = 1

    ltau = .false.
    lsstau = .false.
    nuse = 0

    xmag = 0.d0
    dimer = 0.d0
    flux_x = 0.d0
    flux_y = 0.d0

    nsweep = 20
    nskip = 1
    nbin = 10
    obs_segment_len = 10

    lwrapu = .false.
    lwrapj = .false.

    num_st_nn = 6
    nsw_stglobal = -1
    icount_nsw_stglobal = 0
    lstglobal = .false.
    llocal = .true.

    nublock = 16

    ! read parameters
    if ( irank.eq.0 ) then
        exists = .false.
        inquire (file = 'ftdqmc.in', exist = exists)
        if ( exists .eqv. .true. ) then
            open(unit=100, file='ftdqmc.in',status='unknown')
            read(100, model_para)
            read(100, ctrl_para)
            close(100)
        end if
    end if

#ifdef MPI
    call mpi_bcast( l,            1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( beta,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( dtau,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( mu,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( muA,          1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( muB,          1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rhub,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( rj,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( js,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( hx,           1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( xmag,         1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( flux_x,       1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( flux_y,       1, mpi_real8,    0, mpi_comm_world, ierr )
    call mpi_bcast( nwrap,        1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nsweep,       1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nbin,         1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( llocal,       1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( nsw_stglobal, 1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( lsstau,       1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( ltau,         1, mpi_logical,  0, mpi_comm_world, ierr )
    call mpi_bcast( nuse,         1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_bcast( nublock,      1, mpi_integer,  0, mpi_comm_world, ierr )
    call mpi_barrier(mpi_comm_world,ierr)
#endif

    ! check parameters
    if(mod(l,4) .ne. 0 ) then
        write(fout, '(a)') 'error in make_tables : requre system size l to be fourfold, namely mod(l,4)=0'
        stop
    end if

    ! tune parameters
    if( rhub .gt. 0.d0 ) lwrapu = .true.
    if( rj .gt. 0.d0 ) lwrapj = .true.
    if( nsw_stglobal .gt. 0 .and. lwrapu ) lstglobal = .true.
    if( .not. llocal .and. lstglobal ) nsw_stglobal = 1
    lq = l*l
    lfam = max(lq/2,1)
    nfam = 1
    ndim = lq   ! the dimension of matrix inside determinant
    ltrot = nint( beta / dtau )

    ! tune para for delay update
    if( lq/5 .lt. 16) then
        nublock = 4
    else if( lq/5 .lt. 32 ) then
        nublock = 8
    else if( lq/5 .lt. 64 ) then
        nublock = 16
    else if( lq/5 .lt. 256 ) then
        nublock = 32
    else ! equal to or greater than 256
        nublock = 64
    end if


    allocate( iwrap_nt(0:ltrot) )
    iwrap_nt(0:ltrot) = 0
    ! set nst, and wrap_step
    if( ltrot .lt. nwrap ) then
        if(irank.eq.0) write(fout,'(a,i3,a)')  " WARNNING, ltrot is less than nwrap ", ltrot, ', do not need stablization '
        nst = 0
    else
        if( mod(ltrot,nwrap) .eq. 0 ) then
            nst = ltrot/nwrap
            allocate( wrap_step(2,nst) )
            do i = 1, nst
                iwrap_nt(i*nwrap) = i
                wrap_step(1,i) = 1+(i-1)*nwrap
                wrap_step(2,i) = i*nwrap
            end do
        else
            nst = ltrot/nwrap + 1
            allocate( wrap_step(2,nst) )
            do i = 1, nst-1
                iwrap_nt(i*nwrap) = i
                wrap_step(1,i) = 1+(i-1)*nwrap
                wrap_step(2,i) = i*nwrap
            end do
            i = nst
            iwrap_nt(ltrot) = i
            wrap_step(1,i) = (i-1)*nwrap+1
            wrap_step(2,i) = ltrot
        end if
    end if

    nmeas_bin = 2*(2*obs_segment_len+1)*nsweep*isize
    weight_track = 0.d0

    a1_p(1) = 1 ; a1_p(2) =  0
    a2_p(1) = 0 ; a2_p(2) =  1
    L1_p = l*a1_p
    L2_p = l*a2_p

    b1_p(1) = 2.d0*pi/dble(l) ; b1_p(2) = 0.d0
    b2_p(1) = 0.d0            ; b2_p(2) = 2.d0*pi/dble(l)

    ! allocate tables
    allocate( list(lq,2) )
    allocate( invlist(l,l) )
    allocate( nnlist(lq, 0:8) )
    allocate( list_plaq(lq, 1:5) )
    allocate( ltpf(max(lq/2,1), 4) )
    allocate( lthf(max(lq/4,1), 2) )
    allocate( lthf2(max(lq/4,1), 2) )
    allocate( lthf3(max(lq/16,1), 8) )

    allocate( orblist(ndim) )

    allocate( latt_imj(lq,lq) )
    allocate( listk(lq,2) )

    allocate( imjdeg(lq) )
    allocate( distance_index(lq), equ_distance(lq) )
    allocate( distance_len(lq) )
    allocate( irre_distance_deg(lq) )
    allocate( irre_distance_len(lq) )

    allocate( wn(-nuse:nuse) )
    do i = -nuse, nuse
        wn(i) = 2.d0*dble(i)*pi/beta
    end do
    allocate( zexpiwt(0:ltrot-1,0:nuse) )
    allocate( zexpiqr(lq,lq) )

#ifdef BREAKUP_T
    allocate( urt(max(lq/2,1),4,4), urtm1(max(lq/2,1),4,4) )
    allocate( urtc(max(lq,1),2,2), urtcm1(max(lq,1),2,2) )
    allocate( urtd(max(lq/2,1),4,4), urtdm1(max(lq/2,1),4,4) )
#ifdef SPINDOWN
    allocate( urt_dn(max(lq/2,1),4,4), urtm1_dn(max(lq/2,1),4,4) )
    allocate( urtc_dn(max(lq,1),2,2), urtcm1_dn(max(lq,1),2,2) )
    allocate( urtd_dn(max(lq/2,1),4,4), urtdm1_dn(max(lq/2,1),4,4) )
#endif
#else
    allocate( urt(ndim,ndim), urtm1(ndim,ndim) )
#ifdef SPINDOWN
    allocate( urt_dn(ndim,ndim), urtm1_dn(ndim,ndim) )
#endif
#endif
    !allocate( hopping_tmp(4,max(lq/2,1)) )
    allocate( hopping_tmp(6,max(lq,1)) )
#ifdef SPINDOWN
    !allocate( hopping_tmp_dn(4,max(lq/2,1)) )
    allocate( hopping_tmp_dn(6,max(lq,1)) )
#endif
    allocate( hop_plusx(lq) )
    allocate( hop_minusx(lq) )

    allocate(grup(ndim,ndim), grdn(ndim,ndim), grupc(ndim,ndim), grdnc(ndim,ndim))

    allocate( Imat(ndim,ndim) )
    call s_identity_z(ndim,Imat)
    allocate( Ivec(ndim) )
    do i = 1, ndim
        Ivec(i) = 1.d0
    end do

    if( lstglobal ) then
        allocate( tentacle(2,lq*ltrot), tentacle_old(2,lq*ltrot) )
        allocate( stcluster(lq,ltrot) )
        allocate( stbonds_neib(2,6,lq,ltrot) )
    end if
  
  end subroutine make_tables

  subroutine deallocate_tables
    if( lstglobal ) then
        deallocate( stbonds_neib )
        deallocate( stcluster )
        deallocate( tentacle_old, tentacle )
    end if
    deallocate( Ivec, Imat )
    deallocate( grdnc, grupc, grdn, grup )
    deallocate( hop_minusx )
    deallocate( hop_plusx )
#ifdef SPINDOWN
    deallocate( hopping_tmp_dn )
#endif
    deallocate( hopping_tmp )

#ifdef SPINDOWN
    deallocate( urtm1_dn, urt_dn )
#endif
    deallocate( urtm1, urt )
    deallocate( zexpiqr )
    deallocate( zexpiwt )
    deallocate( wn )

    deallocate( irre_distance_len )
    deallocate( irre_distance_deg )
    deallocate( distance_len )
    deallocate( equ_distance, distance_index )
    deallocate( imjdeg )

#ifdef BREAKUP_T
    deallocate(lthf, lthf2, lthf3)
    deallocate( urtcm1, urtc )
    deallocate( urtdm1, urtd )
#ifdef SPINDOWN                      
    deallocate( urtcm1_dn, urtc_dn )
    deallocate( urtdm1_dn, urtd_dn )
#endif
#endif


    deallocate( listk, latt_imj )
    deallocate( orblist, ltpf, list_plaq, nnlist, invlist, list )
    if(allocated(wrap_step) ) deallocate( wrap_step )
    deallocate( iwrap_nt )
  end subroutine deallocate_tables

end module blockc
