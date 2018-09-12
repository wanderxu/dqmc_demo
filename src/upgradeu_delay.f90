subroutine upgradeu(ntau, green_up, green_dn)

#ifdef _OPENMP
  USE OMP_LIB
#endif
  use spring
  use blockc
#ifdef CUMC
  use mod_cumulate, only: heff, nei_cord, nei_Jeff, num_nei
#endif

  implicit none

  !arguments
  integer,intent(in) :: ntau
  complex(dp), intent(inout), dimension(ndim,ndim) :: green_up, green_dn

  !local
  complex(dp) ::  ratioup, ratiodn, ratiotot, del44_up, del44_dn
  integer :: i4, nl, nl1, nl2, nrflip, nfb, id, ntm1, nta1
  real(dp) :: accm, ratio_re, ratio_re_abs, random
#ifdef CUMC
  integer :: inn, j, ntj
#endif
  complex(dp), allocatable, dimension(:) :: diagg_up
  complex(dp), allocatable, dimension(:) :: diagg_dn
  complex(dp), allocatable, dimension(:,:) :: avec_up, bvec_up
  complex(dp), allocatable, dimension(:,:) :: avec_dn, bvec_dn
  complex(dp) :: alpha_up
  complex(dp) :: alpha_dn
  integer :: i, ik, m

  allocate( diagg_up(ndim) )
#ifdef SPINDOWN
  allocate( diagg_dn(ndim) )
#endif
  allocate( avec_up(ndim,nublock) )
  allocate( bvec_up(ndim,nublock) )
#ifdef SPINDOWN
  allocate( avec_dn(ndim,nublock) )
  allocate( bvec_dn(ndim,nublock) )
#endif

  accm  = 0.d0
  ik = 0
  ! initial diag G
  do i = 1, ndim
      diagg_up(i) = green_up(i,i)
#ifdef SPINDOWN
      diagg_dn(i) = green_dn(i,i)
#endif
  end do
  ! intial avec, bvec
  avec_up = czero
  bvec_up = czero
#ifdef SPINDOWN
  avec_dn = czero
  bvec_dn = czero
#endif
  do i4 = 1,  lq
      ! delay update: after nublock steps of local update, perform a whole update of Green function
      ! calculate weight ratio, fermion part
      nrflip = 1
      del44_up   =  delta_u_up( nsigl_u(i4,ntau), nrflip )
      ratioup = dcmplx(1.d0,0.d0) + del44_up * ( cone - diagg_up(i4) )
#ifdef TEST
      write(fout,'(a,2e16.8)') 'in upgradeu, ratioup = ', ratioup
#endif
#ifdef SPINDOWN
      del44_dn   =  delta_u_dn( nsigl_u(i4,ntau), nrflip )
      ratiodn = dcmplx(1.d0,0.d0) + del44_dn * ( cone - diagg_dn(i4) )
#ifdef TEST
      write(fout,'(a,2e16.8)') 'in upgradeu, ratiodn = ', ratiodn
#endif
#endif
      ! calculate weight ratio, boson part
      id = 0
      ntm1 = ntau - 1
      if ( ntm1 .lt. 1 ) ntm1 = ntm1 + ltrot
      nta1 = ntau + 1
      if ( nta1 .gt. ltrot ) nta1 = nta1 - ltrot
      if( nsigl_u( i4, ntau ) .eq. 1 ) id = ibset( id, 6 )
      if( nsigl_u( i4, ntm1 ) .eq. 1 ) id = ibset( id, 5 )
      if( nsigl_u( i4, nta1 ) .eq. 1 ) id = ibset( id, 4 )
      do nfb = 1, 4
          if( nsigl_u( nnlist(i4,nfb), ntau ) .eq. 1 ) id = ibset( id, 4-nfb )
      end do
      id = id + 1
      ! total ratio
#ifdef SPINDOWN
      ratiotot = (ratioup*ratiodn)*dconjg(ratioup*ratiodn) * wsxsz(id) !* deta_u( nsigl_u(i4,ntau), nrflip )
#else
      ratiotot = ratioup*dconjg(ratioup) * wsxsz(id) ! * deta_u( nsigl_u(i4,ntau), nrflip )
#endif
      ! set alpha, will be used during update avec, bvec
      alpha_up = del44_up/ratioup
#ifdef SPINDOWN
      alpha_dn = del44_dn/ratiodn
#endif
      ! real part of ratio
      ratio_re = dble( ratiotot ) ! * dgaml(nsigl_u(i4,ntau),nrflip)
#ifdef TEST
      write(fout,'(a,2e16.8)') 'in upgradeu, ratio_re = ', ratio_re
#endif
      ratio_re_abs = ratio_re
      ! absolute ratio
      if (ratio_re .lt. 0.d0 )  ratio_re_abs = - ratio_re 
      random = spring_sfmt_stream()
      !! perform update
      if ( ratio_re_abs .gt. random ) then
      ! update accepted
          accm  = accm + 1.d0
          weight_track = weight_track + log( ratio_re_abs )
          logweightf_old = logweightf_old + log( (ratioup*ratiodn)*dconjg(ratioup*ratiodn) )
          logweights_old = logweights_old + log( wsxsz(id) )

          ik = ik + 1
          ! store avec(:,ik) and bvec(:,ik)
          avec_up(:,ik) = green_up(:,i4)
          bvec_up(:,ik) = green_up(i4,:)
#ifdef SPINDOWN
          avec_dn(:,ik) = green_dn(:,i4)
          bvec_dn(:,ik) = green_dn(i4,:)
#endif
          do m = 1, ik-1
              avec_up(:,ik) = avec_up(:,ik) + bvec_up(i4,m)*avec_up(:,m)
              bvec_up(:,ik) = bvec_up(:,ik) + avec_up(i4,m)*bvec_up(:,m)
#ifdef SPINDOWN
              avec_dn(:,ik) = avec_dn(:,ik) + bvec_dn(i4,m)*avec_dn(:,m)
              bvec_dn(:,ik) = bvec_dn(:,ik) + avec_dn(i4,m)*bvec_dn(:,m)
#endif
          end do
          avec_up(:,ik) =avec_up(:,ik)*alpha_up
          bvec_up(i4,ik)=bvec_up(i4,ik) - cone
#ifdef SPINDOWN
          avec_dn(:,ik) =avec_dn(:,ik)*alpha_dn
          bvec_dn(i4,ik)=bvec_dn(i4,ik) - cone
#endif
          ! update diag G
          do i = 1, ndim
              diagg_up(i) = diagg_up(i) + avec_up(i,ik)*bvec_up(i,ik)
#ifdef SPINDOWN
              diagg_dn(i) = diagg_dn(i) + avec_dn(i,ik)*bvec_dn(i,ik)
#endif
          end do
#ifdef CUMC
          ! update heff
          do inn = 1, num_nei
              j = nei_cord(1,inn,i4,ntau)
              ntj = nei_cord(2,inn,i4,ntau)
              heff(j,ntj) = heff(j,ntj) - 2.d0*nei_Jeff(inn,i4,ntau)*nsigl_u(i4,ntau)
          end do
#endif
          ! flip filed
          nsigl_u(i4,ntau) =  nflipl(nsigl_u(i4,ntau), nrflip)
      end if

      if( (ik.eq.nublock) .or. (i4.eq.lq) ) then
          ik = 0
          ! delay update: update the whole Green function
          call  zgemm('N', 'T', ndim, ndim, nublock, cone, avec_up, ndim, bvec_up, ndim, cone, green_up, ndim)
#ifdef SPINDOWN
          call  zgemm('N', 'T', ndim, ndim, nublock, cone, avec_dn, ndim, bvec_dn, ndim, cone, green_dn, ndim)
#endif
          if( i4.lt.lq) then
              ! initial diag G
              do i = 1, ndim
              diagg_up(i) = green_up(i,i)
#ifdef SPINDOWN
                  diagg_dn(i) = green_dn(i,i)
#endif
              end do
              ! intial avec, bvec
              avec_up = czero
              bvec_up = czero
#ifdef SPINDOWN
              avec_dn = czero
              bvec_dn = czero
#endif
          end if
      end if
  end do
  main_obs(1) = main_obs(1) + dcmplx( accm, dble(lq) )

#ifdef SPINDOWN
  deallocate( bvec_dn )
  deallocate( avec_dn )
#endif
  deallocate( bvec_up )
  deallocate( avec_up )
#ifdef SPINDOWN
  deallocate( diagg_dn )
#endif
  deallocate( diagg_up )
end subroutine upgradeu
