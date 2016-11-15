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
      allocate( Ust_up(ndim,ndim,0:nst) )     ! 1
      allocate( Dst_up(ndim,0:nst) )          ! 2
      allocate( Vst_up(ndim,ndim,0:nst) )     ! 3
      allocate( UR_up(ndim,ndim) )             ! 4
      allocate( DRvec_up(ndim) )               ! 5
      allocate( VR_up(ndim,ndim) )             ! 6
      allocate( VL_up(ndim,ndim) )             ! 7
      allocate( DLvec_up(ndim) )               ! 8
      allocate( UL_up(ndim,ndim) )             ! 9
      allocate( Bdtau1_up(ndim,ndim) )       ! 16
      if(ltau) then
          allocate( Bt2t1_up(ndim,ndim) )       ! 16
          allocate( gt0up(ndim,ndim) )
          allocate( g0tup(ndim,ndim) )
          allocate( g00up(ndim,ndim) )
      end if

      allocate( Ust_up_tmp(ndim,ndim,0:nst) ) ! 17
      allocate( Dst_up_tmp(ndim,0:nst) )      ! 18
      allocate( Vst_up_tmp(ndim,ndim,0:nst) ) ! 19
      allocate( grup_tmp(ndim,ndim) )         ! 20

      allocate( Bdtau1_dn(ndim,ndim) )       ! 16
#IFDEF SPINDOWN
      allocate( Ust_dn(ndim,ndim,0:nst) )     ! 1
      allocate( Dst_dn(ndim,0:nst) )          ! 2
      allocate( Vst_dn(ndim,ndim,0:nst) )     ! 3
      allocate( UR_dn(ndim,ndim) )             ! 4
      allocate( DRvec_dn(ndim) )               ! 5
      allocate( VR_dn(ndim,ndim) )             ! 6
      allocate( VL_dn(ndim,ndim) )             ! 7
      allocate( DLvec_dn(ndim) )               ! 8
      allocate( UL_dn(ndim,ndim) )             ! 9
      if(ltau) then
          allocate( Bt2t1_dn(ndim,ndim) )       ! 16
          allocate( gt0dn(ndim,ndim) )
          allocate( g0tdn(ndim,ndim) )
          allocate( g00dn(ndim,ndim) )
      end if

      allocate( Ust_dn_tmp(ndim,ndim,0:nst) ) ! 17
      allocate( Dst_dn_tmp(ndim,0:nst) )      ! 18
      allocate( Vst_dn_tmp(ndim,ndim,0:nst) ) ! 19
      allocate( grdn_tmp(ndim,ndim) )         ! 20
#ENDIF

    end subroutine allocate_core

    subroutine deallocate_core
      implicit none
      deallocate( Bdtau1_dn )      ! 16
#IFDEF SPINDOWN
      deallocate( grdn_tmp )         ! 20
      deallocate( Vst_dn_tmp )       ! 19
      deallocate( Dst_dn_tmp )       ! 18
      deallocate( Ust_dn_tmp )       ! 17
      if(ltau) then
          deallocate( g00dn )
          deallocate( g0tdn )
          deallocate( gt0dn )
          deallocate( Bt2t1_dn )      ! 16
      end if
      deallocate( UL_dn )             ! 9
      deallocate( DLvec_dn )          ! 8
      deallocate( VL_dn )             ! 7
      deallocate( VR_dn )             ! 6
      deallocate( DRvec_dn )          ! 5
      deallocate( UR_dn )             ! 4
      deallocate( Vst_dn )         ! 3
      deallocate( Dst_dn )         ! 2
      deallocate( Ust_dn )         ! 1
#ENDIF
      deallocate( grup_tmp )         ! 20
      deallocate( Vst_up_tmp )       ! 19
      deallocate( Dst_up_tmp )       ! 18
      deallocate( Ust_up_tmp )       ! 17
      if(ltau) then
          deallocate( g00up )
          deallocate( g0tup )
          deallocate( gt0up )
          deallocate( Bt2t1_up )      ! 16
      end if
      deallocate( Bdtau1_up )      ! 16
      deallocate( UL_up )             ! 9
      deallocate( DLvec_up )          ! 8
      deallocate( VL_up )             ! 7
      deallocate( VR_up )             ! 6
      deallocate( DRvec_up )          ! 5
      deallocate( UR_up )             ! 4
      deallocate( Vst_up )         ! 3
      deallocate( Dst_up )         ! 2
      deallocate( Ust_up )         ! 1
    end subroutine deallocate_core
  
    subroutine ftdqmc_stablize_0b_svd(n)
      ! B( n*tau1, 0 ) = B( n*tau1, (n-1)*tau1 ) * B( (n-1)*tau1, 0 )
      implicit none
      integer, intent(in) :: n

      ! local
      integer :: i
      complex(dp), allocatable, dimension(:,:) :: Umat1, Umat2, Vmat1, Vmat2
      real(dp), allocatable, dimension(:) :: Dvec1, Dvec2

      allocate( Umat1(ndim,ndim), Umat2(ndim,ndim), Vmat1(ndim,ndim), Vmat2(ndim,ndim) )
      allocate( Dvec1(ndim), Dvec2(ndim) )

      call Bmat_tau( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn )
      
#IFDEF TEST_LEVEL3
      write(fout, '(a)') ' Bdtau1_up(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2f8.3))') Bdtau1_up(i,:)
      end do
#ENDIF

      Umat1(:,:) = Ust_up(:,:,n-1)
      Dvec1(:)   = Dst_up(:,n-1)
      Vmat1(:,:) = Vst_up(:,:,n-1)

      ! Btmp = ( Bdtau1_up * Umat1 ) * Dmat1
      call zgemm('n','n',ndim,ndim,ndim,cone,Bdtau1_up,ndim,Umat1,ndim,czero,Atmp,ndim)  ! Atmp = Bdtau1_up * Umat1
      call s_z_x_diag_d(ndim,Atmp,Dvec1,Btmp) ! Btmp = Atmp * Dmat1

#IFDEF TEST_LEVEL3
      write(fout, '(a)') ' Btmp(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2e16.8))') Btmp(i,:)
      end do
#ENDIF
      call s_svd_zg(ndim, ndim, ndim, Btmp, Umat2, Dvec2, Vtmp)
!!!#IFDEF TEST
!!!      ! test SVD
!!!      call s_z_x_diag_d(ndim,Umat2,Dvec2,Atmp)  ! Atmp = Umat2 * Dmat2
!!!      call zgemm('n','n',ndim,ndim,ndim,cone,Atmp,ndim,Vtmp,ndim,czero,Btmp,ndim)    ! Btmp = Atmp * Vtmp
!!!      write(fout, '(a)') ' Btmp(:,:) = '
!!!      do i = 1, ndim
!!!          write(fout,'(4(2f8.3))') Btmp(i,:)
!!!      end do
!!!      stop
!!!#ENDIF
      call zgemm('n','n',ndim,ndim,ndim,cone,Vtmp,ndim,Vmat1,ndim,czero,Vmat2,ndim)  ! Vmat2 = Vtmp * Vmat1
      Ust_up(:,:,n) = Umat2(:,:)
      Dst_up(:,n)   = Dvec2(:)
      Vst_up(:,:,n) = Vmat2(:,:)

#IFDEF SPINDOWN
      Umat1(:,:) = Ust_dn(:,:,n-1)
      Dvec1(:)   = Dst_dn(:,n-1)
      Vmat1(:,:) = Vst_dn(:,:,n-1)
      call zgemm('n','n',ndim,ndim,ndim,cone,Bdtau1_dn,ndim,Umat1,ndim,czero,Atmp,ndim)  ! Atmp = Bdtau1_dn * Umat1
      call s_z_x_diag_d(ndim,Atmp,Dvec1,Btmp) ! Btmp = Atmp * Dmat1
      call s_svd_zg(ndim, ndim, ndim, Btmp, Umat2, Dvec2, Vtmp)
      call zgemm('n','n',ndim,ndim,ndim,cone,Vtmp,ndim,Vmat1,ndim,czero,Vmat2,ndim)  ! Vmat2 = Vtmp * Vmat1
      Ust_dn(:,:,n) = Umat2(:,:)
      Dst_dn(:,n)   = Dvec2(:)
      Vst_dn(:,:,n) = Vmat2(:,:)
#ENDIF

      deallocate( Dvec2, Dvec1 )
      deallocate( Vmat2, Vmat1, Umat2, Umat1 )

    end subroutine ftdqmc_stablize_0b_svd
  
    subroutine ftdqmc_stablize_b0_svd(n)
      ! B( beta, (n-1)*tau1 ) = B( beta, n*tau1 ) * B( n*tau1, (n-1)*tau1 )
      implicit none
      integer, intent(in) :: n

      ! local
      integer :: i
      complex(dp), allocatable, dimension(:,:) :: Umat1, Umat2, Vmat1, Vmat2
      real(dp), allocatable, dimension(:) :: Dvec1, Dvec2

      allocate( Umat1(ndim,ndim), Umat2(ndim,ndim), Vmat1(ndim,ndim), Vmat2(ndim,ndim) )
      allocate( Dvec1(ndim), Dvec2(ndim) )

      call Bmat_tau( wrap_step(2,n), wrap_step(1,n), Bdtau1_up, Bdtau1_dn )

#IFDEF TEST_LEVEL3
      write(fout, '(a)') ' Bdtau1_up(:,:) = '
      do i = 1, ndim
          write(fout,'(4(2f8.3))') Bdtau1_up(i,:)
      end do
#ENDIF

      Vmat1(:,:) = Vst_up(:,:,n)
      Dvec1(:)   = Dst_up(:,n)
      Umat1(:,:) = Ust_up(:,:,n)

      ! Btmp = Dmat1 * Umat1 * Bdtau1_up
      call zgemm('n','n',ndim,ndim,ndim,cone,Umat1,ndim,Bdtau1_up,ndim,czero,Atmp,ndim)  ! Atmp = Umat1 * Bdtau1_up
      call s_diag_d_x_z(ndim,Dvec1,Atmp,Btmp) ! Btmp = Dmat1 * Atmp
      call s_svd_zg(ndim, ndim, ndim, Btmp, Vtmp, Dvec2, Umat2)  ! Btmp = Vtmp * Dmat2 * Umat2
      call zgemm('n','n',ndim,ndim,ndim,cone,Vmat1,ndim,Vtmp,ndim,czero,Vmat2,ndim)  ! Vmat2 = Vmat1 * Vtmp 
      Vst_up(:,:,n-1) = Vmat2(:,:)
      Dst_up(:,n-1) = Dvec2(:)
      Ust_up(:,:,n-1) = Umat2(:,:)

#IFDEF SPINDOWN
      Vmat1(:,:) = Vst_dn(:,:,n)
      Dvec1(:)   = Dst_dn(:,n)
      Umat1(:,:) = Ust_dn(:,:,n)
      call zgemm('n','n',ndim,ndim,ndim,cone,Umat1,ndim,Bdtau1_dn,ndim,czero,Atmp,ndim)  ! Atmp = Umat1 * Bdtau1_dn
      call s_diag_d_x_z(ndim,Dvec1,Atmp,Btmp) ! Btmp = Dmat1 * Atmp
      call s_svd_zg(ndim, ndim, ndim, Btmp, Vtmp, Dvec2, Umat2)  ! Btmp = Vtmp * Dmat2 * Umat2
      call zgemm('n','n',ndim,ndim,ndim,cone,Vmat1,ndim,Vtmp,ndim,czero,Vmat2,ndim)  ! Vmat2 = Vmat1 * Vtmp 
      Vst_dn(:,:,n-1) = Vmat2(:,:)
      Dst_dn(:,n-1) = Dvec2(:)
      Ust_dn(:,:,n-1) = Umat2(:,:)
#ENDIF

      deallocate( Dvec2, Dvec1 )
      deallocate( Vmat2, Vmat1, Umat2, Umat1 )

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
          write(fout, '(a,i4,a)') ' Dst_up(:,', n, ' ) = '
          write(fout,'(4(e16.8))') Dst_up(:,n)
#ENDIF
      end do
  
      ! at tau = beta
      UR_up(:,:) = Ust_up(:,:,nst)
      DRvec_up(:)= Dst_up(:,nst)
      VR_up(:,:) = Vst_up(:,:,nst)
      call green_equaltime( nst, ndim, UR_up, DRvec_up, VR_up, Imat, Ivec, Imat, grup, info )

#IFDEF SPINDOWN
      UR_dn(:,:) = Ust_dn(:,:,nst)
      DRvec_dn(:)= Dst_dn(:,nst)
      VR_dn(:,:) = Vst_dn(:,:,nst)
      call green_equaltime( nst, ndim, UR_dn, DRvec_dn, VR_dn, Imat, Ivec, Imat, grdn, info )
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

     call Bmat_tau( ltrot, 1, Bdtau1_up, Bdtau1_dn )
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
      call green_equaltime( nst, ndim, Imat, Ivec, Imat, VL_up, DLvec_up, UL_up, grup, info )

#IFDEF SPINDOWN
      UL_dn(:,:) = Ust_dn(:,:,0)
      DLvec_dn(:)= Dst_dn(:,0)
      VL_dn(:,:) = Vst_dn(:,:,0)
      call green_equaltime( nst, ndim, Imat, Ivec, Imat, VL_dn, DLvec_dn, UL_dn, grdn, info )
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

     call Bmat_tau( ltrot, 1, Bdtau1_up, Bdtau1_dn )
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
     END IF
  
    end subroutine ftdqmc_sweep_start_b0
  
    subroutine ftdqmc_sweep_b0(lupdate, lmeasure)
      implicit none
      logical, intent(in) :: lupdate, lmeasure
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
          if( lmeasure .and. ( abs(nt-nt_ob) .le. obs_segment_len .or. abs(nt-nt_ob) .ge. (ltrot-obs_segment_len) ) ) then
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
  
    subroutine ftdqmc_sweep_0b(lupdate, lmeasure)
      implicit none
      logical, intent(in) :: lupdate, lmeasure
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
                  if( lupdate .and. (.not.ltau .or. .not.lmeasure)) then
                      call upgradej(nt,nf,grup,grdn)
                  end if
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
              if( lupdate .and. (.not.ltau .or. .not.lmeasure)) then
                  call upgradeu( nt, grup, grdn )
              end if
          end if
  
          ! obser
          if( lmeasure .and. ( abs(nt-nt_ob) .le. obs_segment_len .or. abs(nt-nt_ob) .ge. (ltrot-obs_segment_len) ) ) then
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
              if(.not. ltau .or. .not. lmeasure ) then
              !!!if(.not. ltau .or. nt .gt. ltrot/2) then
                  call green_equaltime( n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, grtmp, info )
              else
#IFDEF DYNERROR
                  ! B(nt1,nt2) with nt1 >= nt2
                  nt1 = nt
                  nt2 = nt
                  call Bmat_tau( nt1, nt2, Bdtau1_up, Bdtau1_dn )

                  ! G(t',0) = B(t',t) * G(t,0)
                  Btmp = gt0up
                  call zgemm('n','n',ndim,ndim,ndim,cone,Bdtau1_up,ndim,Btmp,ndim,czero,gt0up,ndim)

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call s_inv_z(ndim,Bdtau1_up)
                  Btmp = g0tup
                  call zgemm('n','n',ndim,ndim,ndim,cone,Btmp,ndim,Bdtau1_up,ndim,czero,g0tup,ndim)
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
              if( .not. ltau .or. .not. lmeasure ) then
              !!!if(.not. ltau .or. nt .gt. ltrot/2) then
                  call green_equaltime( n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, grtmp, info )
              else
#IFDEF DYNERROR
                  ! B(nt1,nt2) with nt1 >= nt2
                  !!!nt1 = nt
                  !!!nt2 = nt
                  !!!call Bmat_tau( nt1, nt2, Bdtau1_up, Bdtau1_dn )

                  ! G(t',0) = B(t',t) * G(t,0)
                  Btmp = gt0dn
                  call zgemm('n','n',ndim,ndim,ndim,cone,Bdtau1_dn,ndim,Btmp,ndim,czero,gt0dn,ndim)

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call s_inv_z(ndim,Bdtau1_dn)
                  Btmp = g0tdn
                  call zgemm('n','n',ndim,ndim,ndim,cone,Btmp,ndim,Bdtau1_dn,ndim,czero,g0tdn,ndim)
#ENDIF

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
  
    subroutine green_equaltime( nt, ndm, ure, dre, vre, vle, dle, ule, gtt, infoe )
      implicit none
      integer, intent(in) :: nt, ndm
      complex(dp), dimension(ndm,ndm), intent(in) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      complex(dp), dimension(ndm,ndm), intent(out) :: gtt
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      complex(dp), allocatable, dimension(:,:) :: ulinv_tmp, urinv_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin

      allocate( ulinv_tmp(ndm,ndm), urinv_tmp(ndm,ndm) )
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

      ! uutmp = ule*ure
      call zgemm('n','n',ndm,ndm,ndm,cone,ule,ndm,ure,ndm,czero,uutmp,ndm)  ! uutmp = ule*ure
      ! vvtmp = vre*vle
      call zgemm('n','n',ndm,ndm,ndm,cone,vre,ndm,vle,ndm,czero,vvtmp,ndm)  ! vvtmp = vre*vle

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
              ulinv_tmp(i,j) = dconjg(ule(j,i))
              urinv_tmp(i,j) = dconjg(ure(j,i))
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, ulinv_tmp, dlmax, dvvdtmp, Btmp )
      call s_v_invd_u( ndm, Btmp, drmax, urinv_tmp, gtt )

      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( urinv_tmp )
      deallocate( ulinv_tmp )

!!!      ! local
!!!      integer :: i
!!!      complex(dp), allocatable, dimension(:,:) :: Umat1, Vmat1
!!!      real(dp), allocatable, dimension(:) :: Dvec1
!!!
!!!      allocate( Umat1(ndm,ndm), Vmat1(ndm,ndm) )
!!!      allocate( Dvec1(ndm) )
!!!
!!!!!!#IFDEF TEST
!!!!!!      write(fout,*)
!!!!!!      write(fout, '(a)') ' before zgemm, ure(:,:) = '
!!!!!!      do i = 1, ndm
!!!!!!          write(fout,'(4(2e16.8))') ure(i,:)
!!!!!!      end do
!!!!!!#ENDIF
!!!      call zgemm('n','n',ndm,ndm,ndm,cone,ule,ndm,ure,ndm,czero,Atmp,ndm)  ! Atmp = ule*ure
!!!!!!#IFDEF TEST
!!!!!!      write(fout,*)
!!!!!!      write(fout, '(a)') ' after zgemm, ure(:,:) = '
!!!!!!      do i = 1, ndm
!!!!!!          write(fout,'(4(2e16.8))') ure(i,:)
!!!!!!      end do
!!!!!!#ENDIF
!!!      call s_inv_z(ndm,Atmp)
!!!      call zgemm('n','n',ndm,ndm,ndm,cone,vre,ndm,vle,ndm,czero,vvtmp,ndm)  ! vvtmp = vre*vle
!!!      call s_diag_dvd(ndim, dre, vvtmp, dle, dvvdtmp )
!!!      !!!if( nt .lt. nst/2 ) then
!!!      !!!    call s_diag_d_x_z(ndm,dre,vvtmp,dvvtmp)  ! dvvtmp = dre * vvtmp
!!!      !!!    call s_z_x_diag_d(ndm,dvvtmp,dle,dvvdtmp)  ! dvvdtmp = dvvtmp * dle
!!!      !!!else
!!!      !!!    call s_z_x_diag_d(ndm,vvtmp,dle,dvvtmp)  ! dvvtmp = vvtmp * dle
!!!      !!!    call s_diag_d_x_z(ndm,dre,dvvtmp,dvvdtmp)  ! dvvdtmp = dre * dvvtmp
!!!      !!!end if
!!!
!!!      Btmp(:,:) = Atmp(:,:) + dvvdtmp(:,:)
!!!
!!!#IFDEF TEST_LEVEL3
!!!      write(fout,*)
!!!      write(fout, '(a)') ' before svd, Btmp(:,:) = '
!!!      do i = 1, ndm
!!!          write(fout,'(4(2e16.8))') Btmp(i,:)
!!!      end do
!!!#ENDIF
!!!
!!!      call s_svd_zg(ndm, ndm, ndm, Btmp, Umat1, Dvec1, Vmat1)
!!!
!!!#IFDEF TEST
!!!      write(fout,*)
!!!      write(fout,'(a)') ' in green_equaltime, after svd, Dvec1 = '
!!!      write(fout,'(4(e16.8))') Dvec1(:)
!!!      write(fout,*)
!!!#ENDIF
!!!      ! check Dvec1
!!!      infoe = 0
!!!      do i = 1, ndim
!!!          if( Dvec1(i) .eq. 0.d0 ) then
!!!              infoe = -1
!!!              write(fout,'(a,i5,a)') 'WARNING!!! Dvec1(',i, ' ) = 0 in green_equaltime !!! '
!!!              return
!!!          end if
!!!      end do
!!!
!!!#IFDEF TEST_LEVEL3
!!!      ! test SVD
!!!      call s_z_x_diag_d(ndim,Umat1,Dvec1,Atmp)  ! Atmp = Umat1 * Dmat1
!!!      call zgemm('n','n',ndim,ndim,ndim,cone,Atmp,ndim,Vmat1,ndim,czero,Btmp,ndim)    ! Btmp = Atmp * Vmat1
!!!      write(fout, '(a)') ' after svd, Btmp(:,:) = '
!!!      do i = 1, ndim
!!!          write(fout,'(4(2e16.8))') Btmp(i,:)
!!!      end do
!!!#ENDIF
!!!
!!!      call zgemm('n','n',ndm,ndm,ndm,cone,Vmat1,ndm,ule,ndm,czero,Atmp,ndm)  ! Atmp = Vmat1 * ule
!!!      call zgemm('n','n',ndm,ndm,ndm,cone,ure,ndm,Umat1,ndm,czero,Btmp,ndm)  ! Btmp = ure * Umat1
!!!      call s_inv_z(ndm, Atmp)
!!!      call s_inv_z(ndm, Btmp)
!!!
!!!      call s_v_invd_u(ndim, Atmp, Dvec1, Btmp, gtt)
!!!
!!!      !!!do i = 1, ndm
!!!      !!!    Dvec1(i) = 1.d0 / Dvec1(i)
!!!      !!!end do
!!!      !!!call s_z_x_diag_d(ndm,Atmp,Dvec1,Vtmp) ! Vtmp = Atmp * Dmat2
!!!      !!!call zgemm('n','n',ndm,ndm,ndm,cone,Vtmp,ndm,Btmp,ndm,czero,gtt,ndm)  ! gtt = Vtmp * Btmp
!!!
!!!!!!#IFDEF TEST
!!!!!!      ! test 
!!!!!!      call s_z_x_diag_d(ndm,ure,dre,Atmp)                                       ! Atmp = ure*dre
!!!!!!      call zgemm('n','n',ndm,ndm,ndm,cone,Atmp,ndm,vre,ndm,czero,Btmp,ndm) ! Btmp = Atmp*vre
!!!!!!      call zgemm('n','n',ndm,ndm,ndm,cone,Btmp,ndm,vle,ndm,czero,Atmp,ndm) ! Atmp = Btmp*vle
!!!!!!      call s_z_x_diag_d(ndm,Atmp,dle,Btmp)                                      ! Btmp = Atmp*dle
!!!!!!      call zgemm('n','n',ndm,ndm,ndm,cone,Btmp,ndm,ule,ndm,czero,Atmp,ndm) ! Atmp = Btmp*ule
!!!!!!      gtt(:,:) = Imat(:,:) + Atmp(:,:)
!!!!!!      call s_inv_z(ndm,gtt)
!!!!!!#ENDIF
!!!      deallocate( Dvec1 )
!!!      deallocate( Vmat1, Umat1 )
    end subroutine green_equaltime

    subroutine green_tau(nt, ndm, ure, dre, vre, vle, dle, ule, g00, gt0, g0t, gtt, infoe )
      implicit none
      integer, intent(in) :: nt, ndm
      complex(dp), dimension(ndm,ndm), intent(inout) :: ure, vre, vle, ule
      real(dp), dimension(ndm), intent(in) :: dre, dle
      complex(dp), dimension(ndm,ndm), intent(out) :: g00, gt0, g0t, gtt
      integer, intent(out) :: infoe

      ! local
      integer :: i, j
      complex(dp), allocatable, dimension(:,:) :: ulinv_tmp, urinv_tmp, vrinv_tmp, vlinv_tmp
      real(dp), allocatable, dimension(:) :: drmax, drmin, dlmax, dlmin

      allocate( ulinv_tmp(ndm,ndm), urinv_tmp(ndm,ndm) )
      allocate( vrinv_tmp(ndm,ndm), vlinv_tmp(ndm,ndm) )
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
      ! uutmp = ule*ure
      call zgemm('n','n',ndm,ndm,ndm,cone,ule,ndm,ure,ndm,czero,uutmp,ndm)  ! uutmp = ule*ure
      ! vvtmp = vre*vle
      call zgemm('n','n',ndm,ndm,ndm,cone,vre,ndm,vle,ndm,czero,vvtmp,ndm)  ! vvtmp = vre*vle

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
              ulinv_tmp(i,j) = dconjg(ule(j,i))
              urinv_tmp(i,j) = dconjg(ure(j,i))
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, ulinv_tmp, dlmax, dvvdtmp, Btmp )
      call s_v_invd_u( ndm, Btmp, drmax, urinv_tmp, gtt )

      !! >> g(t,0)
      call s_v_d_u( ndm, Btmp, drmin, vre, gt0 )

      !!!!! >> g(t,0)
      !!!!  g(t,0) = g(t,t)B(t,0) = g(t,t) * ure * dre * vre
      !!!call zgemm('n','n',ndm,ndm,ndm,cone,gtt,ndm,ure,ndm,czero,Atmp,ndm)  ! Atmp = gtt*ure
      !!!call s_v_d_u( ndm, Atmp, dre, vre, gt0 )


      !! >> g(0,0)
      ! dlmax^-1 * ( vre * vle )^-1 drmax^-1
      ! dlmin * ( ule * ure ) * drmin
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
      do j = 1, ndm
          do i = 1, ndm
              dvvdtmp(i,j) = dconjg(vvtmp(j,i)) / ( dlmax(i)*drmax(j) ) + uutmp(i,j) * dlmin(i) * drmin(j)
              vrinv_tmp(i,j) = dconjg(vre(j,i))
              vlinv_tmp(i,j) = dconjg(vle(j,i))
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      call s_inv_z(ndm,dvvdtmp)
      call s_v_invd_u( ndm, vrinv_tmp, drmax, dvvdtmp, Btmp )
      call s_v_invd_u( ndm, Btmp, dlmax, vlinv_tmp, g00 )

      !! >> g(0,t)
      dlmin(:)=-dlmin(:)
      call s_v_d_u( ndm, Btmp, dlmin, ule, g0t )

      !!!!! >> g(0,t)
      !!!!  g(0,t) = -B(t,0)^-1 * ( 1 - g(t,t) ) = - vre^-1 * dre^-1 * ure^-1 * ( 1- g(t,t) )
      !!!Atmp(:,:) = gtt(:,:) - Imat(:,:)  ! note minus sign here
      !!!call zgemm('n','n',ndm,ndm,ndm,cone,uinv_tmp,ndm,Atmp,ndm,czero,Btmp,ndm)  ! Btmp = ure^-1 * ( g(t,t) - 1 )
      !!!call s_v_invd_u( ndm, vre, dre, Btmp, g0t )

      infoe = 0

      deallocate( dlmin )
      deallocate( dlmax )
      deallocate( drmin )
      deallocate( drmax )
      deallocate( vlinv_tmp )
      deallocate( vrinv_tmp )
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

    subroutine Bmat_tau( nt1, nt2, bmat_up, bmat_dn )
      ! B(tau1,tau2)
      implicit none
      integer, intent(in) :: nt1, nt2
      complex(dp), dimension(ndim,ndim), intent(out) :: bmat_up
      complex(dp), dimension(ndim,ndim), intent(out) :: bmat_dn

      ! local
      integer :: nt, nf, nflag
      complex(dp) :: phaseu

      phaseu = cone
      bmat_up(:,:) = Imat(:,:)
      bmat_dn(:,:) = Imat(:,:)
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
    end subroutine Bmat_tau

#IFDEF CUMC
#include "stglobal_sl.f90"
#ELSE
#include "stglobal_wolff.f90"
#ENDIF

end module ftdqmc_core
