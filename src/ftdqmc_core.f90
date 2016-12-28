module ftdqmc_core
  use spring
  use blockc
  use data_tmp
  use obser
  implicit none

  complex(dp), allocatable, dimension(:,:,:), save :: Ust_up, Vst_up, Ust_up_tmp, Vst_up_tmp
  complex(dp), allocatable, dimension(:,:,:), save :: Ust_dn, Vst_dn, Ust_dn_tmp, Vst_dn_tmp
  real(dp), allocatable, dimension(:,:), save :: Dst_up, Dst_up_tmp
  real(dp), allocatable, dimension(:,:), save :: Dst_dn, Dst_dn_tmp
  complex(dp), dimension(:,:), allocatable, save :: UR_up, VR_up, VL_up, UL_up, Bdtau1_up, grup_tmp
  complex(dp), dimension(:,:), allocatable, save :: UR_dn, VR_dn, VL_dn, UL_dn, Bdtau1_dn, grdn_tmp
  real(dp), dimension(:), allocatable, save :: DRvec_up, DLvec_up
  real(dp), dimension(:), allocatable, save :: DRvec_dn, DLvec_dn
  complex(dp), dimension(:,:), allocatable, save :: Bt2t1_up, gt0up, g0tup, g00up
  complex(dp), dimension(:,:), allocatable, save :: Bt2t1_dn, gt0dn, g0tdn, g00dn


  contains

    subroutine allocate_core
      implicit none
      if(nst.gt.0) then
          allocate( Ust_up(ndim,ndim,0:nst) )     ! 1
          allocate( Dst_up(ndim,0:nst) )          ! 2
          allocate( Vst_up(ndim,ndim,0:nst) )     ! 3
          allocate( UR_up(ndim,ndim) )             ! 4
          allocate( DRvec_up(ndim) )               ! 5
          allocate( VR_up(ndim,ndim) )             ! 6
          allocate( VL_up(ndim,ndim) )             ! 7
          allocate( DLvec_up(ndim) )               ! 8
          allocate( UL_up(ndim,ndim) )             ! 9
      end if
      if( nst.gt.0 .or. llocal ) then
          allocate( Bdtau1_up(ndim,ndim) )       ! 16
      end if
      if(ltau) then
          allocate( Bt2t1_up(ndim,ndim) )       ! 16
          allocate( gt0up(ndim,ndim) )
          allocate( g0tup(ndim,ndim) )
          allocate( g00up(ndim,ndim) )
      end if

      if(nst.gt.0) then
          allocate( Ust_up_tmp(ndim,ndim,0:nst) ) ! 17
          allocate( Dst_up_tmp(ndim,0:nst) )      ! 18
          allocate( Vst_up_tmp(ndim,ndim,0:nst) ) ! 19
      end if
      if( llocal ) then
          allocate( grup_tmp(ndim,ndim) )         ! 20
      end if

      if( nst.gt.0 .or. llocal ) then
          allocate( Bdtau1_dn(ndim,ndim) )       ! 16
      end if
#IFDEF SPINDOWN
      if(nst.gt.0) then
          allocate( Ust_dn(ndim,ndim,0:nst) )     ! 1
          allocate( Dst_dn(ndim,0:nst) )          ! 2
          allocate( Vst_dn(ndim,ndim,0:nst) )     ! 3
          allocate( UR_dn(ndim,ndim) )             ! 4
          allocate( DRvec_dn(ndim) )               ! 5
          allocate( VR_dn(ndim,ndim) )             ! 6
          allocate( VL_dn(ndim,ndim) )             ! 7
          allocate( DLvec_dn(ndim) )               ! 8
          allocate( UL_dn(ndim,ndim) )             ! 9
      end if
      if(ltau) then
          allocate( Bt2t1_dn(ndim,ndim) )       ! 16
          allocate( gt0dn(ndim,ndim) )
          allocate( g0tdn(ndim,ndim) )
          allocate( g00dn(ndim,ndim) )
      end if

      if(nst.gt.0) then
          allocate( Ust_dn_tmp(ndim,ndim,0:nst) ) ! 17
          allocate( Dst_dn_tmp(ndim,0:nst) )      ! 18
          allocate( Vst_dn_tmp(ndim,ndim,0:nst) ) ! 19
      end if
      if( llocal ) then 
          allocate( grdn_tmp(ndim,ndim) )         ! 20
      end if
#ENDIF

    end subroutine allocate_core

    subroutine deallocate_core
      implicit none
      if(allocated(Bdtau1_dn)) deallocate( Bdtau1_dn )      ! 16
#IFDEF SPINDOWN
      if(allocated(grdn_tmp)) deallocate( grdn_tmp )         ! 20
      if(nst.gt.0) then
          deallocate( Vst_dn_tmp )       ! 19
          deallocate( Dst_dn_tmp )       ! 18
          deallocate( Ust_dn_tmp )       ! 17
      end if
      if(ltau) then
          deallocate( g00dn )
          deallocate( g0tdn )
          deallocate( gt0dn )
          deallocate( Bt2t1_dn )      ! 16
      end if
      if(nst.gt.0) then
          deallocate( UL_dn )             ! 9
          deallocate( DLvec_dn )          ! 8
          deallocate( VL_dn )             ! 7
          deallocate( VR_dn )             ! 6
          deallocate( DRvec_dn )          ! 5
          deallocate( UR_dn )             ! 4
          deallocate( Vst_dn )         ! 3
          deallocate( Dst_dn )         ! 2
          deallocate( Ust_dn )         ! 1
      end if
#ENDIF
      if(allocated(grup_tmp)) deallocate( grup_tmp )         ! 20
      if(nst.gt.0) then
          deallocate( Vst_up_tmp )       ! 19
          deallocate( Dst_up_tmp )       ! 18
          deallocate( Ust_up_tmp )       ! 17
      end if
      if(ltau) then
          deallocate( g00up )
          deallocate( g0tup )
          deallocate( gt0up )
          deallocate( Bt2t1_up )      ! 16
      end if
      if(allocated(Bdtau1_up)) deallocate( Bdtau1_up )      ! 16
      if(nst.gt.0) then
          deallocate( UL_up )             ! 9
          deallocate( DLvec_up )          ! 8
          deallocate( VL_up )             ! 7
          deallocate( VR_up )             ! 6
          deallocate( DRvec_up )          ! 5
          deallocate( UR_up )             ! 4
          deallocate( Vst_up )         ! 3
          deallocate( Dst_up )         ! 2
          deallocate( Ust_up )         ! 1
      end if
    end subroutine deallocate_core
  
    subroutine ftdqmc_stablize_0b_svd(n)
      ! B( n*tau1, 0 ) = B( n*tau1, (n-1)*tau1 ) * B( (n-1)*tau1, 0 )
      implicit none
      integer, intent(in) :: n

      ! local
      integer :: i
      integer, allocatable, dimension(:) :: jpvt
      complex(dp), allocatable, dimension(:,:) :: Umat2, Vmat1, Vmat2
      real(dp), allocatable, dimension(:) :: Dvec1, Dvec2

      allocate( jpvt(ndim) )
      allocate( Umat2(ndim,ndim), Vmat1(ndim,ndim), Vmat2(ndim,ndim) )
      allocate( Dvec1(ndim), Dvec2(ndim) )
#IFDEF TEST
      write(fout,'(a)') ' '
      write(fout,'(a)') ' ------------------- '
      write(fout,'(a,i4)') ' In stablize_0b_svd, n = ', n
#ENDIF
#IFDEF TEST_LEVEL3
      Bdtau1_up(:,:) = Imat(:,:)
      Bdtau1_dn(:,:) = Imat(:,:)
      call Bmat_tau_R( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn )
      write(fout, '(a)') ' '
      write(fout, '(a)') 'in stablize_0b_svd, Bmat_tau_R give Bdtau1_up(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') Bdtau1_up(i,:)
      end do
      write(fout, '(a)') 'in stablize_0b_svd, Bmat_tau_R give Bdtau1_dn(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') Bdtau1_dn(i,:)
      end do
#ENDIF

      Bdtau1_up(:,:) = Ust_up(:,:,n-1)
      Bdtau1_dn(:,:) = Ust_dn(:,:,n-1)
      ! Bdtau1_up = Bup(tau+dtau,tau)*U_up
      ! Bdtau1_dn = Bdn(tau+dtau,tau)*U_dn
      call Bmat_tau_R( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn )

      if( n.eq.1 ) then
#IFDEF TEST
          write(fout, '(a)') 'B*U*D*V ='
          call s_v_d_u(ndim, Bdtau1_up, Ivec, Imat, Btmp)
          do i = 1, 1
              write(fout,'(16(2e12.4))') Btmp(i,1:ndim)
          end do
#ENDIF
          call s_zgeQRPT(ndim, ndim, Bdtau1_up, Umat2, Vtmp, jpvt)
#IFDEF TEST_LEVEL3
          write(fout, '(a)') 'B*U*D = Q * R * jpvt, with '
          write(fout, '(a)') 'Q ='
          do i = 1, ndim
              write(fout,'(16(2e12.4))') Umat2(i,:)
          end do
          write(fout, '(a)') 'R ='
          do i = 1, ndim
              write(fout,'(16(2e12.4))') Vtmp(i,:)
          end do
#ENDIF
#IFDEF TEST
          write(fout,'(a,16(i6))') 'jpvt = ', jpvt(:)
#ENDIF
          do i = 1, ndim
              Dvec2(i) = Vtmp(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vtmp, Vmat2 )
          call s_zmp(ndim,ndim,Vmat2,jpvt)
      else
          Dvec1(:)   = Dst_up(:,n-1)
          Vmat1(:,:) = Vst_up(:,:,n-1)
#IFDEF TEST
          write(fout, '(a)') 'B*U*D*V ='
          call s_v_d_u(ndim, Bdtau1_up, Dvec1, Vmat1, Btmp)
          do i = 1, 1
              write(fout,'(16(2e12.4))') Btmp(i,1:ndim)
          end do
#ENDIF
          call s_z_x_diag_d(ndim,Bdtau1_up,Dvec1,Btmp) ! Btmp = Bdtau1_up * Dmat1
          call s_zmcpt(ndim,ndim,Btmp,jpvt)
#IFDEF TEST
          write(fout,'(a)') 'after s_zmcpt, jpvt = '
          write(fout,*) jpvt(:)
#ENDIF
          call s_zgeQR(ndim,ndim,Btmp,Umat2,Vmat2)
#IFDEF TEST_LEVEL3
          write(fout, '(a)') 'B*U*D = Q * R, with '
          write(fout, '(a)') 'Q ='
          do i = 1, ndim
              write(fout,'(16(2e12.4))') Umat2(i,:)
          end do
          write(fout, '(a)') 'R ='
          do i = 1, ndim
              write(fout,'(16(2e12.4))') Vmat2(i,:)
          end do
#ENDIF
          do i = 1, ndim
              Dvec2(i) = Vmat2(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vmat2, Vtmp )
#IFDEF TEST_LEVEL3
          write(fout, '(a)') 'Vmat1 ='
          do i = 1, ndim
              write(fout,'(16(2e12.4))') Vmat1(i,:)
          end do
#ENDIF
          call s_zpm(ndim,ndim,jpvt,Vmat1)
#IFDEF TEST_LEVEL3
          write(fout, '(a)') 'after s_zpm, Vmat1 ='
          do i = 1, ndim
              write(fout,'(16(2e12.4))') Vmat1(i,:)
          end do
#ENDIF
          call zgemm('n','n',ndim,ndim,ndim,cone,Vtmp,ndim,Vmat1,ndim,czero,Vmat2,ndim)
      end if
#IFDEF TEST
      write(fout, '(a)') 'B*U*D*V = Up*Dp*Vp, with '
      write(fout, '(a)') 'Up*Dp*Vp ='
      call s_v_d_u(ndim, Umat2, Dvec2, Vmat2, Btmp)
      do i = 1, 1
          write(fout,'(16(2e12.4))') Btmp(i,1:ndim)
      end do
#ENDIF
#IFDEF TEST_LEVEL3
      write(fout, '(a)') 'Up ='
      do i = 1, 4
          write(fout,'(4(2e12.4))') Umat2(i,1:4)
      end do
      write(fout, '(a,4e12.4)') 'Dp =', Dvec2(1:4)
      write(fout, '(a)') 'Vp ='
      do i = 1, 4
          write(fout,'(4(2e12.4))') Vmat2(i,1:4)
      end do
#ENDIF
      Ust_up(:,:,n) = Umat2(:,:)
      Dst_up(:,n)   = Dvec2(:)
      Vst_up(:,:,n) = Vmat2(:,:)
#IFDEF SPINDOWN
      if( n.eq.1 ) then
          call s_zgeQRPT(ndim, ndim, Bdtau1_dn, Umat2, Vtmp, jpvt)
          do i = 1, ndim
              Dvec2(i) = Vtmp(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vtmp, Vmat2 )
          call s_zmp(ndim,ndim,Vmat2,jpvt)
      else
          Dvec1(:)   = Dst_dn(:,n-1)
          Vmat1(:,:) = Vst_dn(:,:,n-1)
          call s_z_x_diag_d(ndim,Bdtau1_dn,Dvec1,Btmp) ! Btmp = Bdtau1_dn * Dmat1
          call s_zmcpt(ndim,ndim,Btmp,jpvt)
          call s_zgeQR(ndim,ndim,Btmp,Umat2,Vmat2)
          do i = 1, ndim
              Dvec2(i) = Vmat2(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vmat2, Vtmp )
          call s_zpm(ndim,ndim,jpvt,Vmat1)
          call zgemm('n','n',ndim,ndim,ndim,cone,Vtmp,ndim,Vmat1,ndim,czero,Vmat2,ndim)
      end if
      Ust_dn(:,:,n) = Umat2(:,:)
      Dst_dn(:,n)   = Dvec2(:)
      Vst_dn(:,:,n) = Vmat2(:,:)
#ENDIF
#IFDEF TEST
      write(fout,'(a)') ' '
      write(fout,'(a)') ' '
#ENDIF
      deallocate( Dvec2, Dvec1 )
      deallocate( Vmat2, Vmat1, Umat2 )
      deallocate( jpvt )

    end subroutine ftdqmc_stablize_0b_svd
  
    subroutine ftdqmc_stablize_b0_svd(n)
      ! B( beta, (n-1)*tau1 ) = B( beta, n*tau1 ) * B( n*tau1, (n-1)*tau1 )
      implicit none
      integer, intent(in) :: n

      ! local
      integer :: i
      integer, allocatable, dimension(:) :: jpvt
      complex(dp), allocatable, dimension(:,:) :: Umat2, Vmat1, Vmat2
      real(dp), allocatable, dimension(:) :: Dvec1, Dvec2

      allocate( jpvt(ndim) )
      allocate( Umat2(ndim,ndim), Vmat1(ndim,ndim), Vmat2(ndim,ndim) )
      allocate( Dvec1(ndim), Dvec2(ndim) )
#IFDEF TEST
      write(fout,'(a)') ' '
      write(fout,'(a)') ' ------------------- '
      write(fout,'(a,i4)') ' In stablize_b0_svd, n = ', n
#ENDIF
#IFDEF TEST_LEVEL3
      Bdtau1_up(:,:) = Imat(:,:)
      Bdtau1_dn(:,:) = Imat(:,:)
      call Bmat_tau_L( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn )
      write(fout, '(a)') 'in stablize_b0_svd, Bmat_tau_L give Bdtau1_up(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') Bdtau1_up(i,:)
      end do
      write(fout, '(a)') 'in stablize_b0_svd, Bmat_tau_L give Bdtau1_dn(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') Bdtau1_dn(i,:)
      end do

      Bdtau1_up(:,:) = Imat(:,:)
      Bdtau1_dn(:,:) = Imat(:,:)
      call Bmat_tau_RH( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn )
      write(fout, '(a)') 'in stablize_b0_svd, Bmat_tau_RH give Bdtau1_up(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') Bdtau1_up(i,:)
      end do
      write(fout, '(a)') 'in stablize_b0_svd, Bmat_tau_RH give Bdtau1_dn(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') Bdtau1_dn(i,:)
      end do
#ENDIF

      Bdtau1_up(:,:) = Ust_up(:,:,n)
      Bdtau1_dn(:,:) = Ust_dn(:,:,n)
      ! Bdtau1_up = U_up*Bup(tau+dtau,tau)
      ! Bdtau1_dn = U_dn*Bdn(tau+dtau,tau)
      ! call Bmat_tau_L( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn )
      call Bmat_tau_RH( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn ) ! we need to get ( U*B(n*tau1,(n-1)*tau1) )^H = B^H * U^H
                                                                              ! for special model here, (B1B2)^H = B2*B1
      if( n.eq.nst ) then
#IFDEF TEST
          write(fout, '(a)') 'B*U*D*V ='
          call s_v_d_u(ndim, Bdtau1_up, Ivec, Imat, Btmp)
          do i = 1, 4
              write(fout,'(4(2e12.4))') Btmp(i,1:4)
          end do
#ENDIF
          call s_zgeQRPT(ndim, ndim, Bdtau1_up, Umat2, Vtmp, jpvt) ! we have Vtmp^+ and Umat^+, actually
          do i = 1, ndim
              Dvec2(i) = Vtmp(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vtmp, Vmat2 )
          call s_zmp(ndim,ndim,Vmat2,jpvt)  ! we have Vmat2^+, actually
      else
          Dvec1(:)   = Dst_up(:,n)
          Vmat1(:,:) = Vst_up(:,:,n)
#IFDEF TEST
          write(fout, '(a)') 'B*U*D*V ='
          call s_v_d_u(ndim, Bdtau1_up, Dvec1, Vmat1, Btmp)
          do i = 1, 4
              write(fout,'(4(2e12.4))') Btmp(i,1:4)
          end do
#ENDIF
          call s_z_x_diag_d(ndim,Bdtau1_up,Dvec1,Btmp) ! Btmp = Bdtau1_up * Dmat1
          call s_zmcpt(ndim,ndim,Btmp,jpvt)
#IFDEF TEST
          write(fout,'(a)') 'after s_zmcpt, jpvt = '
          write(fout,*) jpvt(:)
#ENDIF
          call s_zgeQR(ndim,ndim,Btmp,Umat2,Vmat2)
          do i = 1, ndim
              Dvec2(i) = Vmat2(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vmat2, Vtmp )


#IFDEF TEST
          write(fout, '(a)') 'Vmat1 ='
          do i = 1, 16
              write(fout,'(16(2e12.4))') Vmat1(i,1:16)
          end do
#ENDIF
          call s_zpm(ndim,ndim,jpvt,Vmat1)
#IFDEF TEST
          write(fout, '(a)') 'after s_zpm, Vmat1 ='
          do i = 1, 16
              write(fout,'(16(2e12.4))') Vmat1(i,1:16)
          end do
          write(fout,'(a,16(i6))') ' jpvt = ', jpvt(:)
#ENDIF

          call zgemm('n','n',ndim,ndim,ndim,cone,Vtmp,ndim,Vmat1,ndim,czero,Vmat2,ndim)
      end if
#IFDEF TEST
      write(fout, '(a)') 'B*U*D*V = Up*Dp*Vp, with '
      write(fout, '(a)') 'Up*Dp*Vp ='
      call s_v_d_u(ndim, Umat2, Dvec2, Vmat2, Btmp)
      do i = 1, 4
          write(fout,'(4(2e12.4))') Btmp(i,1:4)
      end do
#ENDIF
      Ust_up(:,:,n-1) = Umat2(:,:)  ! note we have Umat2^+ and Vmat2^+
      Dst_up(:,n-1)   = Dvec2(:)
      Vst_up(:,:,n-1) = Vmat2(:,:)
#IFDEF SPINDOWN
      if( n.eq.nst ) then
          call s_zgeQRPT(ndim, ndim, Bdtau1_dn, Umat2, Vtmp, jpvt) ! we have Vtmp^+ and Umat^+, actually
          do i = 1, ndim
              Dvec2(i) = Vtmp(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vtmp, Vmat2 )
          call s_zmp(ndim,ndim,Vmat2,jpvt)  ! we have Vmat2^+, actually
      else
          Dvec1(:)   = Dst_dn(:,n)
          Vmat1(:,:) = Vst_dn(:,:,n)
          call s_z_x_diag_d(ndim,Bdtau1_dn,Dvec1,Btmp) ! Btmp = Bdtau1_dn * Dmat1
          call s_zmcpt(ndim,ndim,Btmp,jpvt)
          call s_zgeQR(ndim,ndim,Btmp,Umat2,Vmat2)
          do i = 1, ndim
              Dvec2(i) = Vmat2(i,i)
          end do
          call s_invdiag_d_x_zr( ndim, Dvec2, Vmat2, Vtmp )
          call s_zpm(ndim,ndim,jpvt,Vmat1)
          call zgemm('n','n',ndim,ndim,ndim,cone,Vtmp,ndim,Vmat1,ndim,czero,Vmat2,ndim)
      end if
      Ust_dn(:,:,n-1) = Umat2(:,:)  ! note we have Umat2^+ and Vmat2^+
      Dst_dn(:,n-1)   = Dvec2(:)
      Vst_dn(:,:,n-1) = Vmat2(:,:)
#ENDIF
#IFDEF TEST
      write(fout,'(a)') ' '
      write(fout,'(a)') ' '
#ENDIF
      deallocate( Dvec2, Dvec1 )
      deallocate( Vmat2, Vmat1, Umat2 )
      deallocate( jpvt )
    end subroutine ftdqmc_stablize_b0_svd
  
    subroutine ftdqmc_sweep_start_0b
      implicit none
      integer :: n, i, info
      real(dp) :: tmp

      IF ( nst .gt. 0 ) THEN

      ! at tau = 0
      grup(:,:) = Imat(:,:)
      Ust_up(:,:,0) = Imat(:,:)
      Dst_up(:,0)   = Ivec(:)
      Vst_up(:,:,0) = Imat(:,:)

#IFDEF SPINDOWN
      grdn(:,:) = Imat(:,:)
      Ust_dn(:,:,0) = Imat(:,:)
      Dst_dn(:,0)   = Ivec(:)
      Vst_dn(:,:,0) = Imat(:,:)
#ENDIF

      do n = 1, nst
          ! at tau = n * tau1
          call ftdqmc_stablize_0b_svd(n)
#IFDEF TEST
          write(fout, '(a,i4,a)') ' in ftdqmc_sweep_start_0b, Dst_up(:,', n, ' ) = '
          write(fout,'(4(e16.8))') Dst_up(:,n)
#ENDIF
      end do
  
      ! at tau = beta
      UR_up(:,:) = Ust_up(:,:,nst)
      DRvec_up(:)= Dst_up(:,nst)
      VR_up(:,:) = Vst_up(:,:,nst)
      !call green_equaltime( nst, ndim, UR_up, DRvec_up, VR_up, Imat, Ivec, Imat, grup, info )
      call green_equaltimebb( nst, ndim, UR_up, DRvec_up, VR_up, grup, info )

#IFDEF SPINDOWN
      UR_dn(:,:) = Ust_dn(:,:,nst)
      DRvec_dn(:)= Dst_dn(:,nst)
      VR_dn(:,:) = Vst_dn(:,:,nst)
      !call green_equaltime( nst, ndim, UR_dn, DRvec_dn, VR_dn, Imat, Ivec, Imat, grdn, info )
      call green_equaltimebb( nst, ndim, UR_dn, DRvec_dn, VR_dn, grdn, info )
#ENDIF

      if( info .eq. -1 ) then
          write(fout,'(a)') ' WRONG in sweep_start, exit '
          stop
      end if

#IFDEF TEST_LEVEL3
      write(fout, '(a)') ' After sweep_start_0b, grup(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') grup(i,:)
      end do
#ENDIF


#IFDEF TEST

      tmp = 0.d0
      do i = 1, ndim
          tmp = tmp + real( cone - grup(i,i) )
      end do
      write(fout, '(a)') ' After sweep_start_0b '
      write(fout,'(a,2e12.4)') ' grup(1,1) = ', grup(1,1)
      write(fout,'(a,e12.4)') ' ne_up = ', tmp

#IFDEF SPINDOWN
      tmp = 0.d0
      do i = 1, ndim
          tmp = tmp + real( cone - grdn(i,i) )
      end do
      write(fout,'(a,2e12.4)') ' grdn(1,1) = ', grdn(1,1)
      write(fout,'(a,e12.4)') ' ne_dn = ', tmp
#ENDIF

#ENDIF

     ELSE

     if ( llocal ) then
         Bdtau1_up(:,:) = Imat(:,:)
         Bdtau1_dn(:,:) = Imat(:,:)
         call Bmat_tau_R( ltrot, 1, Bdtau1_up, Bdtau1_dn )
         do  i = 1, ndim
             Bdtau1_up(i,i) = Bdtau1_up(i,i) + cone
         end do
         call s_inv_z( ndim, Bdtau1_up )
         grup(:,:) = Bdtau1_up
#IFDEF SPINDOWN
         do  i = 1, ndim
             Bdtau1_dn(i,i) = Bdtau1_dn(i,i) + cone
         end do
         call s_inv_z( ndim, Bdtau1_dn )
         grdn(:,:) = Bdtau1_dn
#ENDIF
     else
         grup(:,:) = Imat(:,:)
         grdn(:,:) = Imat(:,:)
         call Bmat_tau_R( ltrot, 1, grup, grdn)
         do  i = 1, ndim
             grup(i,i) = grup(i,i) + cone
             grdn(i,i) = grdn(i,i) + cone
         end do
     end if

     END IF
  
    end subroutine ftdqmc_sweep_start_0b

    subroutine ftdqmc_sweep_start_b0
      implicit none
      integer :: n, i, info
      real(dp) :: tmp

      IF ( nst .gt. 0 ) THEN

      ! at tau = beta
      Vst_up(:,:,nst) = Imat(:,:)
      Dst_up(:,nst)   = Ivec(:)
      Ust_up(:,:,nst) = Imat(:,:)
#IFDEF SPINDOWN
      Vst_dn(:,:,nst) = Imat(:,:)
      Dst_dn(:,nst)   = Ivec(:)
      Ust_dn(:,:,nst) = Imat(:,:)
#ENDIF

      do n = nst, 1, -1
          ! at tau = (n-1) * tau1
          ! calculate B(n*tau1,(n-1)*tau1), and set Vst(:,:,n-1), Dst(:,:,n-1), Ust(:,:,n-1)
          call ftdqmc_stablize_b0_svd(n)
#IFDEF TEST
          write(fout, '(a,i4,a)') ' Dst_up(:,', n, ' ) = '
          write(fout,'(4(e16.8))') Dst_up(:,n)
#ENDIF
      end do

      ! at tau = 0
      UL_up(:,:) = Ust_up(:,:,0)
      DLvec_up(:)= Dst_up(:,0)
      VL_up(:,:) = Vst_up(:,:,0)
      !call green_equaltime( nst, ndim, Imat, Ivec, Imat, VL_up, DLvec_up, UL_up, grup, info )
      call green_equaltime00( nst, ndim, VL_up, DLvec_up, UL_up, grup, info )

#IFDEF SPINDOWN
      UL_dn(:,:) = Ust_dn(:,:,0)
      DLvec_dn(:)= Dst_dn(:,0)
      VL_dn(:,:) = Vst_dn(:,:,0)
      !call green_equaltime( nst, ndim, Imat, Ivec, Imat, VL_dn, DLvec_dn, UL_dn, grdn, info )
      call green_equaltime00( nst, ndim, VL_dn, DLvec_dn, UL_dn, grdn, info )
#ENDIF

      if( info .eq. -1 ) then
          write(fout,'(a)') ' WRONG in sweep_start, exit '
          stop
      end if

#IFDEF TEST_LEVEL3
      write(fout, '(a)') ' After sweep_start, grup(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e12.4))') grup(i,:)
      end do
#ENDIF


#IFDEF TEST

      tmp = 0.d0
      do i = 1, ndim
          tmp = tmp + real( cone - grup(i,i) )
      end do
      write(fout,'(a,2e12.4)') ' grup(1,1) = ', grup(1,1)
      write(fout,'(a,e12.4)') ' ne_up = ', tmp

#IFDEF SPINDOWN
      tmp = 0.d0
      do i = 1, ndim
          tmp = tmp + real( cone - grdn(i,i) )
      end do
      write(fout,'(a,2e12.4)') ' grdn(1,1) = ', grdn(1,1)
      write(fout,'(a,e12.4)') ' ne_dn = ', tmp
#ENDIF

#ENDIF

     ELSE

     if ( llocal ) then
         Bdtau1_up(:,:) = Imat(:,:)
         Bdtau1_dn(:,:) = Imat(:,:)
         call Bmat_tau_R( ltrot, 1, Bdtau1_up, Bdtau1_dn )
         do  i = 1, ndim
             Bdtau1_up(i,i) = Bdtau1_up(i,i) + cone
         end do
         call s_inv_z( ndim, Bdtau1_up )
         grup(:,:) = Bdtau1_up
#IFDEF SPINDOWN
         do  i = 1, ndim
             Bdtau1_dn(i,i) = Bdtau1_dn(i,i) + cone
         end do
         call s_inv_z( ndim, Bdtau1_dn )
         grdn(:,:) = Bdtau1_dn
#ENDIF
     else
         grup(:,:) = Imat(:,:)
         grdn(:,:) = Imat(:,:)
         call Bmat_tau_R( ltrot, 1, grup, grdn)
         do  i = 1, ndim
             grup(i,i) = grup(i,i) + cone
             grdn(i,i) = grdn(i,i) + cone
         end do
     end if
     END IF
  
    end subroutine ftdqmc_sweep_start_b0
  
    subroutine ftdqmc_sweep_b0(lupdate, lmeasure_equaltime)
      implicit none
      logical, intent(in) :: lupdate, lmeasure_equaltime
      ! local variables
      integer :: nt, n, nf, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
      logical :: lterminate 
      real(dp) :: tmp, ratiof, ratiofi

      ! at tau = beta
      Vst_up(:,:,nst) = Imat(:,:)
      Dst_up(:,nst)   = Ivec(:)
      Ust_up(:,:,nst) = Imat(:,:)
#IFDEF SPINDOWN
      Vst_dn(:,:,nst) = Imat(:,:)
      Dst_dn(:,nst)   = Ivec(:)
      Ust_dn(:,:,nst) = Imat(:,:)
#ENDIF
  
      nt_ob = ceiling( spring_sfmt_stream() * ltrot )

      do nt = ltrot, 1, -1
#IFDEF TEST
          write(fout,*)
          write(fout,'(a)') " ----------------"
          write(fout, '(a,i8)') ' |=> nt = ',  nt
          write(fout,'(a)') " ----------------"
          write(fout,*)
#ENDIF
          ! obser
          if( lmeasure_equaltime .and. ( abs(nt-nt_ob) .le. obs_segment_len .or. abs(nt-nt_ob) .ge. (ltrot-obs_segment_len) ) ) then
             call obser_equaltime(nt)
          end if
  
          !! update
          ! updateu
          if( lwrapu ) then
              nflag = 3 ! onsite
              if(lupdate) call upgradeu( nt, grup, grdn )
              call mmuul  ( grup, grdn, nf, nt, nflag )
              call mmuurm1( grup, grdn, nf, nt, nflag )
          end if
  
          ! updatej
          if( lwrapj ) then
              do nf = nfam,1,-1
                  nflag = 2
                  call mmuul  ( grup, grdn, nf, nt, nflag )
                  call mmuurm1( grup, grdn, nf, nt, nflag )
                  if(lupdate) call upgradej(nt,nf,grup,grdn)
                  nflag = 1
                  call mmuul  ( grup, grdn, nf, nt, nflag )
                  call mmuurm1( grup, grdn, nf, nt, nflag )
              enddo
          end if
  
          ! wrap H0
          call mmthl  (grup, grdn)
          call mmthrm1(grup, grdn)

          if ( (iwrap_nt(nt-1) .gt. 0 .and. nt.ne.ltrot) .or. ( nt.eq.1 .and. nst .gt. 0 ) ) then
              n = iwrap_nt(nt-1)
              ! at tau = n * tau1
              UR_up(:,:) = Ust_up(:,:,n)
              DRvec_up(:)= Dst_up(:,n)
              VR_up(:,:) = Vst_up(:,:,n)

#IFDEF SPINDOWN
              UR_dn(:,:) = Ust_dn(:,:,n)
              DRvec_dn(:)= Dst_dn(:,n)
              VR_dn(:,:) = Vst_dn(:,:,n)
#ENDIF

              call ftdqmc_stablize_b0_svd(n+1)

              ! for spin up
              UL_up(:,:)  = Ust_up(:,:,n)
              DLvec_up(:) = Dst_up(:,n)  
              VL_up(:,:)  = Vst_up(:,:,n)
              !if( .not. ltau ) then
                  call green_equaltime( n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, grtmp, info )
              !else
              !    call green_tau(n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, g00up, gt0up, g0tup, grtmp, info )
              !end if
              call s_compare_max_z( ndim, grtmp, grup, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
#IFDEF TEST_LEVEL3
              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' progating grup(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grup(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch grup(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grtmp(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_up(:,', n, ' ) before wrap = '
              write(fout,'(4(e16.8))') DRvec_up(:)

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_up(:,', n, ' ) after wrap = '
              write(fout,'(4(e16.8))') Dst_up(:,n)
              tmp = 0.d0
              do i = 1, ndim
                  tmp = tmp + real( cone - grup(i,i) )
              end do
              write(fout,'(a,e12.4)') ' progating  ne_up = ', tmp
              tmp = 0.d0
              do i = 1, ndim
                  tmp = tmp + real( cone - grup(i,i) )
              end do
              write(fout,'(a,e12.4)') ' scratch  ne_up = ', tmp
#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' grup, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
              ! whether use the scrath grup
              if( info .eq. 0 ) grup(:,:) = grtmp(:,:)

#IFDEF SPINDOWN
              ! for spin down
              UL_dn(:,:)  = Ust_dn(:,:,n)
              DLvec_dn(:) = Dst_dn(:,n)  
              VL_dn(:,:)  = Vst_dn(:,:,n)
              call green_equaltime( n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, grtmp, info )
              call s_compare_max_z( ndim, grtmp, grdn, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
#IFDEF TEST_LEVEL3
              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' progating grdn(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grdn(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch grdn(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grtmp(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_dn(:,', n, ' ) before wrap = '
              write(fout,'(4(e16.8))') DRvec_dn(:)

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_dn(:,', n, ' ) after wrap = '
              write(fout,'(4(e16.8))') Dst_dn(:,n)

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_dn(:,', n, ' ) after wrap = '
              write(fout,'(4(e16.8))') Dst_dn(:,n)
              tmp = 0.d0
              do i = 1, ndim
                  tmp = tmp + real( cone - grdn(i,i) )
              end do
              write(fout,'(a,e12.4)') ' progating  ne_dn = ', tmp
              tmp = 0.d0
              do i = 1, ndim
                  tmp = tmp + real( cone - grdn(i,i) )
              end do
              write(fout,'(a,e12.4)') ' scratch  ne_dn = ', tmp
#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' grdn, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
              ! whether use the scrath grdn
              if( info .eq. 0 ) grdn(:,:) = grtmp(:,:)
#ENDIF
          end if
  
  
      end do
    end subroutine ftdqmc_sweep_b0
  
    subroutine ftdqmc_sweep_0b(lupdate, lmeasure_equaltime, lmeasure_dyn )
      implicit none
      logical, intent(in) :: lupdate, lmeasure_equaltime, lmeasure_dyn
      ! local variables
      integer :: nt, n, nf, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
      logical :: lterminate 
      real(dp) :: tmp, ratiof, ratiofi
      ! at tau = 0
      Ust_up(:,:,0) = Imat(:,:)
      Dst_up(:,0)   = Ivec(:)
      Vst_up(:,:,0) = Imat(:,:)
#IFDEF SPINDOWN
      Ust_dn(:,:,0) = Imat(:,:)
      Dst_dn(:,0)   = Ivec(:)
      Vst_dn(:,:,0) = Imat(:,:)
#ENDIF


!!!!#include "stglobal.f90"

      if( ltau ) then
          g00up = grup
          gt0up = grup
          g0tup = grup-Imat
#IFDEF SPINDOWN
          g00dn = grdn
          gt0dn = grdn
          g0tdn = grdn-Imat
#ENDIF
      end if
  
      nt_ob = ceiling( spring_sfmt_stream() * ltrot )
      do nt = 1, ltrot, 1

#IFDEF TEST
          write(fout,*)
          write(fout,'(a)') " ----------------"
          write(fout, '(a,i8)') ' |=> nt = ',  nt
          write(fout,'(a)') " ----------------"
          write(fout,*)
#ENDIF
  
          ! wrap H0
          call mmthr  (grup, grdn)
          call mmthlm1(grup, grdn)

          ! update
          ! updatej
          if( lwrapj ) then
              do nf = 1, nfam
                  nflag = 2
                  call mmuur   (grup, grdn, nf, nt, nflag)
                  call mmuulm1 (grup, grdn, nf, nt, nflag)
                  if( lupdate ) call upgradej(nt,nf,grup,grdn)
                  nflag = 1
                  call mmuur  (grup, grdn, nf, nt, nflag)
                  call mmuulm1(grup, grdn, nf, nt, nflag)
              end do
          end if
          ! updateu
          if( lwrapu ) then
              nflag = 3 ! onsite
              call mmuur  ( grup, grdn, nf, nt, nflag )
              call mmuulm1( grup, grdn, nf, nt, nflag )
              if( lupdate ) call upgradeu( nt, grup, grdn )
          end if
  
          ! obser
          if( lmeasure_equaltime .and. ( abs(nt-nt_ob) .le. obs_segment_len .or. abs(nt-nt_ob) .ge. (ltrot-obs_segment_len) ) ) then
             call obser_equaltime(nt)
          end if
  
          if ( iwrap_nt(nt) .gt. 0 ) then
              n = iwrap_nt(nt)
              ! at tau = n * tau1
              VL_up(:,:) = Vst_up(:,:,n)
              DLvec_up(:)= Dst_up(:,n)
              UL_up(:,:) = Ust_up(:,:,n)
#IFDEF SPINDOWN
              VL_dn(:,:) = Vst_dn(:,:,n)
              DLvec_dn(:)= Dst_dn(:,n)
              UL_dn(:,:) = Ust_dn(:,:,n)
#ENDIF
              call ftdqmc_stablize_0b_svd(n)

              ! for spin up
              UR_up(:,:)  = Ust_up(:,:,n)
              DRvec_up(:) = Dst_up(:,n)
              VR_up(:,:)  = Vst_up(:,:,n)
              if( .not. ltau .or. .not. lmeasure_dyn ) then
                  call green_equaltime( n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, grtmp, info )
              else
              ! only when we need measure dynamical quantities, we will call green_tau
#IFDEF DYNERROR
                  ! B(nt1,nt2) with nt1 >= nt2
                  nt1 = nt
                  nt2 = nt
                  ! G(t',0) = B(t',t) * G(t,0)
                  call Bmat_tau_R( nt1, nt2, gt0up, gt0dn)

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call Bmatinv_tau_L( nt1, nt2, g0tup, g0tdn)
#ENDIF

                  !call green_tau(n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, g00up, gt0up,   g0tup,   grtmp, info )
                  call  green_tau(n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, g00up, gt0tmp,  g0ttmp,  grtmp, info )
#IFDEF TEST_LEVEL3
                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' progating gt0up(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') gt0up(i,:)
                  end do

                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch gt0tmp(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') gt0tmp(i,:)
                  end do

                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' progating g0tup(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') g0tup(i,:)
                  end do

                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch g0ttmp(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') g0ttmp(i,:)
                  end do
#ENDIF

#IFDEF DYNERROR
                  call s_compare_max_z( ndim, gt0up, gt0tmp, xmax_dyn_tmp )
                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
#ENDIF
                  gt0up = gt0tmp
#IFDEF TEST
                  write(fout, '(a,e16.8)') 'gt0up, xmax_dyn_tmp = ',  xmax_dyn_tmp
#ENDIF
#IFDEF DYNERROR
                  call s_compare_max_z( ndim, g0tup, g0ttmp, xmax_dyn_tmp )
                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
#ENDIF
                  g0tup = g0ttmp
#IFDEF TEST
                  write(fout, '(a,e16.8)') 'g0tup, xmax_dyn_tmp = ',  xmax_dyn_tmp
#ENDIF
              end if
              call s_compare_max_z( ndim, grtmp, grup, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
#IFDEF TEST_LEVEL3
              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' progating grup(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grup(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch grup(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grtmp(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_up(:,', n, ' ) before wrap = '
              write(fout,'(4(e16.8))') DLvec_up(:)

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_up(:,', n, ' ) after wrap = '
              write(fout,'(4(e16.8))') Dst_up(:,n)
#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' grup, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
              if( info .eq. 0 ) grup(:,:) = grtmp(:,:)

#IFDEF SPINDOWN
              ! for spin down
              UR_dn(:,:)  = Ust_dn(:,:,n)
              DRvec_dn(:) = Dst_dn(:,n)
              VR_dn(:,:)  = Vst_dn(:,:,n)
              if( .not. ltau .or. .not. lmeasure_dyn ) then
                  call green_equaltime( n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, grtmp, info )
              else
                  !call green_tau(n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, g00dn, gt0dn,  g0tdn,  grtmp, info )
                  call  green_tau(n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, g00dn, gt0tmp, g0ttmp, grtmp, info )
#IFDEF TEST_LEVEL3
                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' progating gt0dn(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') gt0dn(i,:)
                  end do

                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch gt0tmp(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') gt0tmp(i,:)
                  end do

                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' progating g0tdn(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') g0tdn(i,:)
                  end do

                  write(fout,*)
                  write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch g0ttmp(:,:) = '
                  do i = 1, ndim
                      write(fout,'(4(2e12.4))') g0ttmp(i,:)
                  end do
#ENDIF

#IFDEF DYNERROR
                  call s_compare_max_z( ndim, gt0dn, gt0tmp, xmax_dyn_tmp )
                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
#ENDIF
                  gt0dn = gt0tmp
#IFDEF TEST
                  write(fout, '(a,e16.8)') 'gt0dn, xmax_dyn_tmp = ',  xmax_dyn_tmp
#ENDIF
#IFDEF DYNERROR
                  call s_compare_max_z( ndim, g0tdn, g0ttmp, xmax_dyn_tmp )
                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
#ENDIF
                  g0tdn = g0ttmp
#IFDEF TEST
                  write(fout, '(a,e16.8)') 'g0tdn, xmax_dyn_tmp = ',  xmax_dyn_tmp
#ENDIF
              end if
              call s_compare_max_z( ndim, grtmp, grdn, max_wrap_error_tmp )
              if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
#IFDEF TEST_LEVEL3
              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' progating grdn(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grdn(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i5,a)') 'nt = ', nt, ' scratch grdn(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grtmp(i,:)
              end do

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_dn(:,', n, ' ) before wrap = '
              write(fout,'(4(e16.8))') DLvec_dn(:)

              write(fout,*)
              write(fout, '(a,i4,a)') ' Dst_dn(:,', n, ' ) after wrap = '
              write(fout,'(4(e16.8))') Dst_dn(:,n)
#ENDIF

#IFDEF TEST
              write(fout,*)
              write(fout, '(a,e16.8)') ' grdn, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
              if( info .eq. 0 ) grdn(:,:) = grtmp(:,:)
#ENDIF

          end if

#include "dyn.f90"
  
      end do
    end subroutine ftdqmc_sweep_0b
  
    subroutine green_equaltime( n, ndm, ure, dre, vre, vle, dle, ule, gtt, infoe )
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(in) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      complex(dp), dimension(ndm,ndm), intent(out) :: gtt
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      complex(dp), allocatable, dimension(:,:) :: ulinv_tmp, vlhtr_tmp, urinv_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin

      allocate( ulinv_tmp(ndm,ndm), vlhtr_tmp(ndm,ndm), urinv_tmp(ndm,ndm) )
      allocate( drmax(ndm), drmin(ndm), dlmax(ndm), dlmin(ndm) )

      ! breakup dre = drmax * drmin
      !         dle = dlmax * dlmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dre,drmax,drmin)
      call s_dvec_min_max(ndm,dle,dlmax,dlmin)

#IFDEF TEST
      write(fout,*)
      write(fout,'(a,i6)') ' in green_equaltime, n = ', n
      write(fout,'(a)') '     dre(i)        drmax(i)        drmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dre(i), drmax(i), drmin(i)
      end do
      write(fout,*)
      write(fout,'(a)') '     dle(i)        dlmax(i)        dlmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dle(i), dlmax(i), dlmin(i)
      end do
#ENDIF

!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ulinv_tmp(i,j) = dconjg(ule(j,i))
              vlhtr_tmp(i,j) = dconjg(vle(j,i))
              urinv_tmp(i,j) = dconjg(ure(j,i))
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL

      ! uutmp = ule*ure
      call zgemm('n','n',ndm,ndm,ndm,cone,ulinv_tmp,ndm,ure,ndm,czero,uutmp,ndm)  ! uutmp = ule*ure
      ! vvtmp = vre*vle
      call zgemm('n','n',ndm,ndm,ndm,cone,vre,ndm,vlhtr_tmp,ndm,czero,vvtmp,ndm)  ! vvtmp = vre*vle

      !! >> g(t,t)
      ! drmax^-1 * ( ule * ure )^-1 dlmax^-1
      ! drmin * ( vre * vle ) * dlmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ! note uutmp^-1 = uutmp ^ +
              dvvdtmp(i,j) = dconjg(uutmp(j,i)) / ( drmax(i)*dlmax(j) ) + vvtmp(i,j) * drmin(i) * dlmin(j)
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL

      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, ule, dlmax, dvvdtmp, Btmp )
      call s_v_invd_u( ndm, Btmp, drmax, urinv_tmp, gtt )

      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( urinv_tmp )
      deallocate( vlhtr_tmp )
      deallocate( ulinv_tmp )
    end subroutine green_equaltime

    subroutine green_equaltime00( n, ndm, vle, dle, ule, gtt, infoe )
      ! calcultate G(0,0), can save 3 matrix products when compare with green_equaltime
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(in) :: vle, ule
      real(dp), dimension(ndm), intent(in) :: dle
      complex(dp), dimension(ndm,ndm), intent(out) :: gtt
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      real(dp), allocatable, dimension(:) :: dlmax, dlmin

      allocate( dlmax(ndm), dlmin(ndm) )

      ! breakup dle = dlmax * dlmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dle,dlmax,dlmin)

#IFDEF TEST
      write(fout,'(a,i6)') ' in green_equaltime00, n = ', n
      write(fout,'(a)') '     dle(i)        dlmax(i)        dlmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dle(i), dlmax(i), dlmin(i)
      end do
#ENDIF

      !! >> g(0,0)
      ! ule^-1 * dlmax^-1
      ! vle * dlmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ! note uutmp^-1 = uutmp ^ +
              dvvdtmp(i,j) = ule(i,j) / dlmax(j) + dconjg( vle(j,i) ) * dlmin(j)
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, ule, dlmax, dvvdtmp, gtt )

      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
    end subroutine green_equaltime00

    subroutine green_equaltimebb( n, ndm, ure, dre, vre, gtt, infoe )
      ! calcultate G(beta,beta), can save 3 matrix products when compare with green_equaltime
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(in) :: ure, vre
      real(dp), dimension(ndm), intent(in) :: dre
      complex(dp), dimension(ndm,ndm), intent(out) :: gtt
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      complex(dp), allocatable, dimension(:,:) :: urinv_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin

      allocate( urinv_tmp(ndm,ndm) )
      allocate( drmax(ndm), drmin(ndm) )

      ! breakup dre = drmax * drmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dre,drmax,drmin)

#IFDEF TEST
      write(fout,*)
      write(fout,'(a,i6)') ' in green_equaltimebb, n = ', n
      write(fout,'(a)') '     dre(i)        drmax(i)        drmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dre(i), drmax(i), drmin(i)
      end do
#ENDIF

      !! >> g(beta,beta)
      ! drmax^-1 * ure^-1
      ! drmin * vre
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ! note uutmp^-1 = uutmp ^ +
              dvvdtmp(i,j) = dconjg(ure(j,i)) / drmax(i) + vre(i,j) * drmin(i)
              urinv_tmp(i,j) = dconjg(ure(j,i))
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, dvvdtmp, drmax, urinv_tmp, gtt )

      infoe = 0

      deallocate( drmin )
      deallocate( drmax )
      deallocate( urinv_tmp )
    end subroutine green_equaltimebb

    subroutine green_tau(n, ndm, ure, dre, vre, vle, dle, ule, g00, gt0, g0t, gtt, infoe )
      implicit none
      integer, intent(in) :: n, ndm
      complex(dp), dimension(ndm,ndm), intent(inout) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      complex(dp), dimension(ndm,ndm), intent(out) :: g00, gt0, g0t, gtt
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      complex(dp), allocatable, dimension(:,:) :: ulinv_tmp, urinv_tmp, vlhtr_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin

      allocate( ulinv_tmp(ndm,ndm), urinv_tmp(ndm,ndm), vlhtr_tmp(ndm,ndm) )
      allocate( drmax(ndm), drmin(ndm), dlmax(ndm), dlmin(ndm) )

      ! breakup dre = drmax * drmin
      !         dle = dlmax * dlmin
      ! with <1 value set to 1 in dmax,  >1 value to 1 in dmin
      call s_dvec_min_max(ndm,dre,drmax,drmin)
      call s_dvec_min_max(ndm,dle,dlmax,dlmin)

#IFDEF TEST
      write(fout,*)
      write(fout,'(a)') '     dre(i)        drmax(i)        drmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dre(i), drmax(i), drmin(i)
      end do
      write(fout,*)
      write(fout,'(a)') '     dle(i)        dlmax(i)        dlmin(i) '
      do i = 1, ndm
          write(fout, '(3e16.8)') dle(i), dlmax(i), dlmin(i)
      end do
#ENDIF
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ulinv_tmp(i,j) = dconjg(ule(j,i))
              vlhtr_tmp(i,j) = dconjg(vle(j,i))
              urinv_tmp(i,j) = dconjg(ure(j,i))
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL

      ! uutmp = ule*ure
      call zgemm('n','n',ndm,ndm,ndm,cone,ulinv_tmp,ndm,ure,ndm,czero,uutmp,ndm)  ! uutmp = ule*ure
      ! vvtmp = vre*vle
      call zgemm('n','n',ndm,ndm,ndm,cone,vre,ndm,vlhtr_tmp,ndm,czero,vvtmp,ndm)  ! vvtmp = vre*vle

      !! >> g(t,t)
      ! drmax^-1 * ( ule * ure )^-1 dlmax^-1
      ! drmin * ( vre * vle ) * dlmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              ! note uutmp^-1 = uutmp ^ +
              dvvdtmp(i,j) = dconjg(uutmp(j,i)) / ( drmax(i)*dlmax(j) ) + vvtmp(i,j) * drmin(i) * dlmin(j)
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, ule, dlmax, dvvdtmp, Btmp )
      call s_v_invd_u( ndm, Btmp, drmax, urinv_tmp, gtt )

      !! >> g(t,0)
      call s_v_d_u( ndm, Btmp, drmin, vre, gt0 )


      call s_inv_z(ndm,vre)       ! vre has been rewritten here.
      call s_inv_z(ndm,vlhtr_tmp) ! vlhtr_tmp has been rewritten here.
      ! vvtmp = (vre*vle)^-1 = vle^-1 * vre^-1
      call zgemm('n','n',ndm,ndm,ndm,cone,vlhtr_tmp,ndm,vre,ndm,czero,vvtmp,ndm)  ! vvtmp = vle^-1 * vre^-1
      !! >> g(0,0)
      ! dlmax^-1 * ( vre * vle )^-1 drmax^-1
      ! dlmin * ( ule * ure ) * drmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              dvvdtmp(i,j) = vvtmp(i,j) / ( dlmax(i)*drmax(j) ) + uutmp(i,j) * dlmin(i) * drmin(j)
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, vre, drmax, dvvdtmp, Btmp )
      call s_v_invd_u( ndm, Btmp, dlmax, vlhtr_tmp, g00 )

      !! >> g(0,t)
      dlmin(:)=-dlmin(:)  ! note minus sign here
      call s_v_d_u( ndm, Btmp, dlmin, ulinv_tmp, g0t )

      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( vlhtr_tmp )
      deallocate( urinv_tmp )
      deallocate( ulinv_tmp )

    end subroutine green_tau


    !!!subroutine green_tau(nt, ndm, ure, dre, vre, vle, dle, ule, g00, gt0, g0t, gtt, infoe )
    !!!  implicit none
    !!!  integer, intent(in) :: nt, ndm
    !!!  complex(dp), dimension(ndm,ndm), intent(inout) :: ure, vre, vle, ule
    !!!  real(dp), dimension(ndm), intent(in) :: dre, dle
    !!!  complex(dp), dimension(ndm,ndm), intent(out) :: g00, gt0, g0t, gtt
    !!!  integer, intent(out) :: infoe

    !!!  ! local
    !!!  integer :: i, j
    !!!  complex(dp), allocatable, dimension(:,:) :: vrulmat, v2mat, u2mat, vlurmat, udvmat, vuvmat, uvumat
    !!!  real(dp), allocatable, dimension(:) :: d2vec

    !!!  allocate( vrulmat( 2*ndm, 2*ndm ) )   ! 1
    !!!  allocate(   v2mat( 2*ndm, 2*ndm ) )   ! 2
    !!!  allocate(   u2mat( 2*ndm, 2*ndm ) )   ! 3
    !!!  allocate( vlurmat( 2*ndm, 2*ndm ) )   ! 4
    !!!  allocate(  udvmat( 2*ndm, 2*ndm ) )   ! 5
    !!!  allocate(  vuvmat( 2*ndm, 2*ndm ) )   ! 6
    !!!  allocate(  uvumat( 2*ndm, 2*ndm ) )   ! 7
    !!!  allocate(   d2vec( 2*ndm ) )          ! 8

    !!!  call zgemm('n','n',ndm,ndm,ndm,cone,vre,ndm,vle,ndm,czero,vvtmp,ndm)  ! vvtmp = vre*vle
    !!!  call s_inv_z(ndm,vvtmp)
    !!!  call zgemm('n','n',ndm,ndm,ndm,cone,ule,ndm,ure,ndm,czero,Atmp,ndm)  ! Atmp = ule*ure
    !!!  call s_inv_z(ndm,Atmp)
    !!!  udvmat = czero
    !!!  do j = 1, ndm
    !!!      do i = 1, ndm
    !!!          udvmat(i,     j     ) =  vvtmp(i,j)
    !!!          udvmat(i+ndm, j+ndm ) =  Atmp(i,j)
    !!!          if(i.eq.j) then
    !!!              udvmat(i+ndm, j     ) = dcmplx( -dre(i), 0.d0 )
    !!!              udvmat(i    , j+ndm ) = dcmplx(  dle(i), 0.d0 )
    !!!          end if
    !!!      end do
    !!!  end do
    !!!  call s_svd_zg(2*ndm, 2*ndm, 2*ndm, udvmat, u2mat, d2vec, v2mat)

    !!!  ! check d2vec
    !!!  infoe = 0
    !!!  do i = 1, 2*ndm
    !!!      if( d2vec(i) .eq. 0.d0 ) then
    !!!          infoe = -1
    !!!          write(fout,'(a,i5,a)') 'WARNING!!! d2vec(',i, ' ) = 0 in green_tau !!! '
    !!!      end if
    !!!  end do

    !!!  call s_inv_z(2*ndm, v2mat)
    !!!  call s_inv_z(2*ndm, u2mat)

    !!!  !! attention here, we are now changing vre, ule, vle, ure
    !!!  call s_inv_z(ndm,vre)
    !!!  call s_inv_z(ndm,ule)
    !!!  call s_inv_z(ndm,vle)
    !!!  call s_inv_z(ndm,ure)

    !!!  vrulmat = czero
    !!!  vrulmat(1:ndm,1:ndm) = vre(1:ndm,1:ndm)
    !!!  vrulmat(ndm+1:2*ndm, ndm+1:2*ndm) = ule(1:ndm,1:ndm)

    !!!  vlurmat = czero
    !!!  vlurmat(1:ndm,1:ndm) = vle(1:ndm,1:ndm)
    !!!  vlurmat(ndm+1:2*ndm, ndm+1:2*ndm) = ure(1:ndm,1:ndm)

    !!!  call zgemm('n','n',2*ndm,2*ndm,2*ndm,cone,vrulmat,2*ndm,v2mat,  2*ndm,czero,vuvmat,2*ndm)
    !!!  call zgemm('n','n',2*ndm,2*ndm,2*ndm,cone,u2mat,  2*ndm,vlurmat,2*ndm,czero,uvumat,2*ndm)

    !!!  call s_v_invd_u( 2*ndm, vuvmat, d2vec, uvumat, udvmat )

    !!!  do j = 1, ndm
    !!!      do i = 1, ndm
    !!!          g00(i,j) = udvmat(i,j)
    !!!          gt0(i,j) = udvmat(i+ndm,j)
    !!!          g0t(i,j) = udvmat(i,j+ndm)
    !!!          gtt(i,j) = udvmat(i+ndm,j+ndm)
    !!!      end do
    !!!  end do

    !!!  deallocate(   d2vec )   ! 8
    !!!  deallocate(  uvumat )   ! 7
    !!!  deallocate(  vuvmat )   ! 6
    !!!  deallocate(  udvmat )   ! 5
    !!!  deallocate( vlurmat )   ! 4
    !!!  deallocate(   u2mat )   ! 3
    !!!  deallocate(   v2mat )   ! 2
    !!!  deallocate( vrulmat )   ! 1
    !!!end subroutine green_tau

    subroutine Bmat_tau_R( nt1, nt2, bmat_up, bmat_dn )
      ! B(tau1,tau2) * 
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_dn

      ! local
      integer :: nt, nf, nflag
      complex(dp) :: phaseu

      phaseu = cone
      do nt = nt2, nt1
          call mmthr(bmat_up,bmat_dn)
          if( lwrapj ) then
            do nf = 1, nfam
                nflag = 2
                call mmuur(bmat_up, bmat_dn, nf, nt, nflag)
                nflag = 1
                call mmuur(bmat_up, bmat_dn, nf, nt, nflag)
            end do
          end if
          if( lwrapu ) then
            nflag = 3 ! onsite
            call mmuur(bmat_up, bmat_dn, nf, nt, nflag )
            !call get_phase_u(nt, phaseui)
            !phaseu = phaseu * phaseui
          end if
      end do
      !phaseu = dcmplx( (0.5d0*exp(-dtau*rhub*0.5d0) ) ** (lq*(abs(nt1-nt2)+1)), 0.d0 )
      !write(fout,'(a,2e16.8)') ' phaseu = ', phaseu
      !bmat(:,:) = bmat(:,:) * phaseu
    end subroutine Bmat_tau_R

    subroutine Bmat_tau_RH( nt1, nt2, bmat_up, bmat_dn )
      ! B(tau1,tau2) * 
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_dn

      ! local
      integer :: nt, nf, nflag
      complex(dp) :: phaseu

      phaseu = cone
      do nt = nt1, nt2, -1
          if( lwrapu ) then
            nflag = 3 ! onsite
            call mmuurH(bmat_up, bmat_dn, nf, nt, nflag )
          end if
          !!!!if( lwrapj ) then
          !!!!  do nf = nfam, 1, -1
          !!!!      nflag = 2
          !!!!      call mmuur(bmat_up, bmat_dn, nf, nt, nflag)
          !!!!      nflag = 1
          !!!!      call mmuur(bmat_up, bmat_dn, nf, nt, nflag)
          !!!!  end do
          !!!!end if
          call mmthrH(bmat_up,bmat_dn)
      end do
    end subroutine Bmat_tau_RH

    subroutine Bmat_tau_L( nt1, nt2, bmat_up, bmat_dn )
      ! * B(tau1,tau2)
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_dn

      ! local
      integer :: nt, nf, nflag
      complex(dp) :: phaseu

      phaseu = cone
      do nt = nt1, nt2, -1
          if( lwrapu ) then
            nflag = 3 ! onsite
            call mmuul(bmat_up, bmat_dn, nf, nt, nflag )
          end if
          if( lwrapj ) then
            do nf = nfam, 1, -1
                nflag = 2
                call mmuul(bmat_up, bmat_dn, nf, nt, nflag)
                nflag = 1
                call mmuul(bmat_up, bmat_dn, nf, nt, nflag)
            end do
          end if
          call mmthl(bmat_up,bmat_dn)
      end do
    end subroutine Bmat_tau_L

    subroutine Bmatinv_tau_L( nt1, nt2, bmat_up, bmat_dn )
      ! *B(tau1,tau2)^-1
      ! make sure nt1 > nt2
      implicit none
      integer, intent(in) :: nt1, nt2
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up
      complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_dn

      ! local
      integer :: nt, nf, nflag
      complex(dp) :: phaseu

      phaseu = cone
      do nt = nt2, nt1
          call mmthlm1(bmat_up,bmat_dn)
          if( lwrapj ) then
            do nf = 1, nfam
                nflag = 2
                call mmuulm1(bmat_up, bmat_dn, nf, nt, nflag)
                nflag = 1
                call mmuulm1(bmat_up, bmat_dn, nf, nt, nflag)
            end do
          end if
          if( lwrapu ) then
            nflag = 3 ! onsite
            call mmuulm1(bmat_up, bmat_dn, nf, nt, nflag )
          end if
      end do
    end subroutine Bmatinv_tau_L

#IFDEF CUMC
#include "stglobal_sl.f90"
#ELSE
#include "stglobal_wolff.f90"
#ENDIF

    subroutine push_stage
      implicit none
      Ust_up_tmp(:,:,:) = Ust_up(:,:,:)
      Dst_up_tmp(:,:)   = Dst_up(:,:)
      Vst_up_tmp(:,:,:) = Vst_up(:,:,:)
      grup_tmp(:,:)     = grup(:,:)
#IFDEF SPINDOWN
      Ust_dn_tmp(:,:,:) = Ust_dn(:,:,:)
      Dst_dn_tmp(:,:)   = Dst_dn(:,:)
      Vst_dn_tmp(:,:,:) = Vst_dn(:,:,:)
      grdn_tmp(:,:)     = grdn(:,:)
#ENDIF
    end subroutine

    subroutine pop_stage
      implicit none
      Ust_up(:,:,:) =  Ust_up_tmp(:,:,:)
      Dst_up(:,:)   =  Dst_up_tmp(:,:)
      Vst_up(:,:,:) =  Vst_up_tmp(:,:,:)
      grup(:,:)     =  grup_tmp(:,:)
#IFDEF SPINDOWN
      Ust_dn(:,:,:) =  Ust_dn_tmp(:,:,:)
      Dst_dn(:,:)   =  Dst_dn_tmp(:,:)
      Vst_dn(:,:,:) =  Vst_dn_tmp(:,:,:)
      grdn(:,:)     =  grdn_tmp(:,:)
#ENDIF
    end subroutine pop_stage

end module ftdqmc_core
