module blockc
  integer, parameter :: dp = 8
  integer, parameter :: fout = 50
  integer, save :: irank, isize, ierr
  real(dp), parameter :: zero = 0.d0
  complex(dp), parameter :: cone = dcmplx( 1.d0, 0.d0 )
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
  real(dp), save :: beta
  real(dp), save :: mu
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
  real(dp), save :: tanhdth, cothdth



  ! dqmc relative
  logical, save :: lupdateu, lupdatej, lstglobal
  integer(dp), save :: iseed
  real(dp), save :: dtau
  integer, save :: ltrot
  complex(dp), save :: phase
  real(dp), save :: weight_track
  complex(dp), allocatable, dimension(:,:), save :: Imat
  complex(dp), allocatable, dimension(:,:), save :: grup, grdn, grupc, grdnc
  real(dp), allocatable, dimension(:), save :: Ivec

  ! cal. control
  integer, save :: nbin
  integer, save :: nsweep
  integer, save :: nst
  integer, save :: nwrap
  integer, save :: lwarnup
  integer, save :: nwarnup
  integer, save :: n_outconf_pace 
  integer, allocatable, dimension(:) :: iwrap_nt
  integer, allocatable, dimension(:,:) :: wrap_step

  ! for dynamical
  logical, save :: lsstau
  logical, save :: lsstau0r
  logical, save :: ltau
  logical, save :: ltauall
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

#IFDEF BREAKUP_T
  complex(dp), allocatable, dimension(:,:,:), save :: urt, urtm1
#IFDEF SPINDOWN
  complex(dp), allocatable, dimension(:,:,:), save :: urt_dn, urtm1_dn
#ENDIF
#ELSE
  complex(dp), allocatable, dimension(:,:), save :: urt, urtm1
#IFDEF SPINDOWN
  complex(dp), allocatable, dimension(:,:), save :: urt_dn, urtm1_dn
#ENDIF
#ENDIF
  complex(dp), allocatable, dimension(:,:), save :: hopping_tmp
#IFDEF SPINDOWN
  complex(dp), allocatable, dimension(:,:), save :: hopping_tmp_dn
#ENDIF

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
  
  contains

  subroutine make_tables
    use parser, only : p_create, p_parse, p_get, p_get_vec, p_destroy
    use mmpi, only : mp_bcast, mp_barrier
    implicit none

    include 'mpif.h'

    integer :: i, nwrap_mid
    logical :: exists

    ! default parameters
    l    = 2
    beta = 20
    dtau = 0.05d0
    mu   = 0.d0 ! default is half filling
    rhub = 1.0d0
    rj   = 0.d0
    js   = 1.d0
    hx   = 5.d0
    nwrap = 10
    n_outconf_pace = 1

    lsstau = .false.
    lsstau0r = .false.
    ltau = .false.
    ltauall = .false.
    nuse = 0

    xmag = 0.d0
    dimer = 0.d0
    flux_x = 0.d0
    flux_y = 0.d0

    nsweep = 20
    nbin = 10
    obs_segment_len = 10

    lupdateu = .false.
    lupdatej = .false.

    num_st_nn = 6
    nsw_stglobal = -1
    icount_nsw_stglobal = 0
    lstglobal = .false.
    
    ! read parameters
    if ( irank.eq.0 ) then
        exists = .false.
        inquire (file = 'ftdqmc.in', exist = exists)
        if ( exists .eqv. .true. ) then
            call p_create()
            call p_parse('ftdqmc.in')
            call p_get( 'L'        , l       )            ! 1
            call p_get( 'beta'     , beta    )            ! 2
            call p_get( 'dtau'     , dtau    )            ! 3
            call p_get( 'mu'       , mu      )            ! 3.5
            call p_get( 'rhub'     , rhub    )            ! 4
            call p_get( 'rj'       , rj      )            ! 5
            call p_get( 'js'       , js      )            ! 6
            call p_get( 'hx'       , hx      )            ! 7
            call p_get( 'xmag'     , xmag    )            ! 7
            call p_get( 'flux_x'   , flux_x  )            ! 7
            call p_get( 'flux_y'   , flux_y  )            ! 7
            call p_get( 'nwrap'    , nwrap   )            ! 8
            call p_get( 'nsweep'   , nsweep  )            ! 9
            call p_get( 'nbin'     , nbin    )            ! 10
            call p_get( 'nsw_stglobal', nsw_stglobal )    ! 11
            call p_get( 'lsstau'     , lsstau    )
            call p_get( 'lsstau0r'     , lsstau0r    )
            call p_get( 'ltau'     , ltau    )
            call p_get( 'ltauall'  , ltauall )
            call p_get( 'nuse'     , nuse    )
            call p_destroy()
        end if
    end if

    call mp_bcast( l,    0 )                ! 1
    call mp_bcast( beta, 0 )                ! 2
    call mp_bcast( dtau, 0 )                ! 3
    call mp_bcast( mu,   0 )                ! 3.5
    call mp_bcast( rhub, 0 )                ! 4
    call mp_bcast( rj,   0 )                ! 5
    call mp_bcast( js,   0 )                ! 6
    call mp_bcast( hx,   0 )                ! 7
    call mp_bcast( xmag, 0 )                ! 7
    call mp_bcast( flux_x, 0 )                ! 7
    call mp_bcast( flux_y, 0 )                ! 7
    call mp_bcast( nwrap, 0 )               ! 8
    call mp_bcast( nsweep, 0 )              ! 9
    call mp_bcast( nbin, 0 )                ! 10
    call mp_bcast( nsw_stglobal, 0 )    ! 11
    call mp_bcast( lsstau, 0 )
    call mp_bcast( lsstau0r, 0 )
    call mp_bcast( ltau, 0 )
    call mp_bcast( ltauall, 0 )
    call mp_bcast( nuse, 0 )
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! tune parameters
    if( rhub .gt. 0.d0 ) lupdateu = .true.
    if( rj .gt. 0.d0 ) lupdatej = .true.
    if( nsw_stglobal .gt. 0 .and. lupdateu ) lstglobal = .true.
    lq = l*l
    lfam = max(lq/2,1)
    nfam = 1
    ndim = lq   ! the dimension of matrix inside determinant
    ltrot = nint( beta / dtau )

    allocate( iwrap_nt(ltrot) )
    iwrap_nt(:) = 0
    ! set nst, and wrap_step
    if( mod(ltrot,2) .eq. 0 ) then
        if( mod(ltrot/2,nwrap) .eq. 0 ) then
            nst = (ltrot/2/nwrap)*2
            allocate( wrap_step(2,nst) )
            do i = 1, nst
                iwrap_nt(i*nwrap) = i
                wrap_step(1,i) = 1+(i-1)*nwrap
                wrap_step(2,i) = i*nwrap
            end do
        else
            nst = (ltrot/2/nwrap+1)*2
            nwrap_mid = ltrot/2-(nst/2-1)*nwrap
            allocate( wrap_step(2,nst) )
            do i = 1, nst/2-1
                iwrap_nt(i*nwrap) = i
                wrap_step(1,i) = 1+(i-1)*nwrap
                wrap_step(2,i) = i*nwrap
            end do

            i = nst/2
            iwrap_nt(ltrot/2) = i
            wrap_step(1,i) = (i-1)*nwrap+1
            wrap_step(2,i) = ltrot/2

            i = nst/2+1
            iwrap_nt(ltrot/2+nwrap_mid) = i
            wrap_step(1,i) = ltrot/2+1
            wrap_step(2,i) = ltrot/2+nwrap_mid

            do i = 1, nst/2-1
                iwrap_nt(ltrot/2+nwrap_mid+i*nwrap) = nst/2+1+i
                wrap_step(1,i+nst/2+1) = ltrot/2+nwrap_mid+(i-1)*nwrap+1
                wrap_step(2,i+nst/2+1) = ltrot/2+nwrap_mid+i*nwrap
            end do
        end if
    else
        write(*,*)  " ltrot should be even number ", ltrot
        stop
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

#IFDEF BREAKUP_T
    allocate( urt(max(lq/2,1),4,4), urtm1(max(lq/2,1),4,4) )
#IFDEF SPINDOWN
    allocate( urt_dn(max(lq/2,1),4,4), urtm1_dn(max(lq/2,1),4,4) )
#ENDIF
#ELSE
    allocate( urt(ndim,ndim), urtm1(ndim,ndim) )
#IFDEF SPINDOWN
    allocate( urt_dn(ndim,ndim), urtm1_dn(ndim,ndim) )
#ENDIF
#ENDIF
    allocate( hopping_tmp(4,max(lq/2,1)) )
#IFDEF SPINDOWN
    allocate( hopping_tmp_dn(4,max(lq/2,1)) )
#ENDIF

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
#IFDEF SPINDOWN
    deallocate( urtm1_dn, urt_dn )
#ENDIF
    deallocate( urtm1, urt )
    deallocate( zexpiqr )
    deallocate( zexpiwt )
    deallocate( wn )

    deallocate( irre_distance_len )
    deallocate( irre_distance_deg )
    deallocate( distance_len )
    deallocate( equ_distance, distance_index )
    deallocate( imjdeg )

    deallocate( listk, latt_imj )
    deallocate( orblist, lthf, ltpf, list_plaq, nnlist, invlist, list )
    deallocate( wrap_step )
    deallocate( iwrap_nt )
  end subroutine deallocate_tables

end module blockc
