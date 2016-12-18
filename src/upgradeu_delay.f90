subroutine upgradeu(ntau, green_up, green_dn)

#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use spring
  use blockc
#IFDEF CUMC
  use mod_cumulate, only: heff, nei_cord, nei_Jeff, num_nei
#ENDIF

  implicit none

  !arguments
  integer,intent(in) :: ntau
  complex(dp), intent(inout), dimension(ndim,ndim) :: green_up, green_dn

  !local
  complex(dp) ::  ratioup, ratiodn, ratiotot, del44_up, del44_dn
  integer :: i4, nl, nl1, nl2, nrflip, nfb, id, ntm1, nta1
  real(dp) :: accm, ratio_re, ratio_re_abs, random
#IFDEF CUMC
  integer :: inn, j, ntj
#ENDIF
  complex(dp), allocatable, dimension(:) :: diagg_up
  complex(dp), allocatable, dimension(:) :: diagg_dn
  complex(dp), allocatable, dimension(:,:) :: avec_up, bvec_up
  complex(dp), allocatable, dimension(:,:) :: avec_dn, bvec_dn
  complex(dp) :: alpha_up
  complex(dp) :: alpha_dn
  integer :: i, iblock, ik

  allocate( diagg_up(ndim) )
#IFDEF SPINDOWN
  allocate( diagg_dn(ndim) )
#ENDIF
  allocate( avec_up(ndim,ublock_Nsites) )
  allocate( bvec_up(ndim,ublock_Nsites) )
#IFDEF SPINDOWN
  allocate( avec_dn(ndim,ublock_Nsites) )
  allocate( bvec_dn(ndim,ublock_Nsites) )
#ENDIF

  accm  = 0.d0
  do iblock = 1, num_ublock
      ! initial diag G
      do i = 1, ndim
          diagg_up(i) = green_up(i,i)
#IFDEF SPINDOWN
          diagg_dn(i) = green_dn(i,i)
#ENDIF
      end do
      ! perform ublock_Nsites steps of local update for one delay update
      do ik = 1, ublock_Nsites
          i4 = (iblock-1)*ublock_Nsites + ik
          ! calculate weight ratio, fermion part
          nrflip = 1
          del44_up   =  delta_u_up( nsigl_u(i4,ntau), nrflip )
          ratioup = dcmplx(1.d0,0.d0) + del44_up * ( cone - diagg_up(i4) )
#IFDEF SPINDOWN
          del44_dn   =  delta_u_dn( nsigl_u(i4,ntau), nrflip )
          ratiodn = dcmplx(1.d0,0.d0) + del44_dn * ( cone - diagg_dn(i4) )
#ENDIF
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
#IFDEF SPINDOWN
          ratiotot = (ratioup*ratiodn)*dconjg(ratioup*ratiodn) * wsxsz(id) !* deta_u( nsigl_u(i4,ntau), nrflip )
#ELSE
          ratiotot = ratioup*dconjg(ratioup) * wsxsz(id) ! * deta_u( nsigl_u(i4,ntau), nrflip )
#ENDIF
          ! set alpha, will be used during update avec, bvec
          alpha_up = del44_up/ratioup
#IFDEF SPINDOWN
          alpha_dn = del44_dn/ratiodn
#ENDIF
          ! real part of ratio
          ratio_re = dble( ratiotot ) ! * dgaml(nsigl_u(i4,ntau),nrflip)
          ratio_re_abs = ratio_re
          ! absolute ratio
          if (ratio_re .lt. 0.d0 )  ratio_re_abs = - ratio_re 
          random = spring_sfmt_stream()
          !! perform update
          if ( ratio_re_abs .gt. random ) then
          ! update accepted
              accm  = accm + 1.d0
              weight_track = weight_track + log( ratio_re_abs )
              ! store avec(:,ik) and bvec(:,ik)
              avec_up(:,ik) = greep_up(:,i4)
              bvec_up(:,ik) = green_up(i4,:)
#IFDEF SPINDOWN
              avec_dn(:,ik) = greep_dn(:,i4)
              bvec_dn(:,ik) = green_dn(i4,:)
#ENDIF
              do m = 1, ik-1
                  avec_up(:,ik) = avec_up(:,ik) + bvec_up(i4,m)*avec_up(:,m)
                  bvec_up(:,ik) = bvec_up(:,ik) + avec_up(i4,m)*bvec_up(:,m)
#IFDEF SPINDOWN
                  avec_dn(:,ik) = avec_dn(:,ik) + bvec_dn(i4,m)*avec_dn(:,m)
                  bvec_dn(:,ik) = bvec_dn(:,ik) + avec_dn(i4,m)*bvec_dn(:,m)
#ENDIF
              end do
              avec_up(:,ik) =avec_up(:,ik)*alpha_up
              bvec_up(i4,ik)=bvec_up(i4,ik) - cone
#IFDEF SPINDOWN
              avec_dn(:,ik) =avec_dn(:,ik)*alpha_up
              bvec_dn(i4,ik)=bvec_dn(i4,ik) - cone
#ENDIF
              ! update diag G
              do i = 1, ndim
                  diagg_up(i) = diagg_up(i) + avec_up(i)*bvec_up(i)
#IFDEF SPINDOWN
                  diagg_dn(i) = diagg_dn(i) + avec_dn(i)*bvec_dn(i)
#ENDIF
              end do
#IFDEF CUMC
              ! update heff
              do inn = 1, num_nei
                  j = nei_cord(1,inn,i4,ntau)
                  ntj = nei_cord(2,inn,i4,ntau)
                  heff(j,ntj) = heff(j,ntj) - 2.d0*nei_Jeff(inn,i4,ntau)*nsigl_u(i4,ntau)
              end do
#ENDIF
              ! flip filed
              nsigl_u(i4,ntau) =  nflipl(nsigl_u(i4,ntau), nrflip)
          else
          ! update rejected
              ! set avec(:,ik) and bvec(:,ik) to zero
              avec_up(:,ik) = czero
              bvec_up(:,ik) = czero
#IFDEF SPINDOWN
              avec_dn(:,ik) = czero
              bvec_dn(:,ik) = czero
#ENDIF
              ! do not need update diag G
          end if
      end do
      ! delay update: update the whole Green function
      call  zgemm('N', 'T', ndim, ndim, ublock_Nsites, cone, avec_up, ndim, bvec_up, ndim, cone, green_up)
#IFDEF SPINDOWN
      call  zgemm('N', 'T', ndim, ndim, ublock_Nsites, cone, avec_dn, ndim, bvec_dn, ndim, cone, green_dn)
#ENDIF
  end do
  main_obs(1) = main_obs(1) + dcmplx( accm, dble(lq) )

#IFDEF SPINDOWN
  deallocate( bvec_dn )
  deallocate( avec_dn )
#ENDIF
  deallocate( bvec_up )
  deallocate( avec_up )
#IFDEF SPINDOWN
  deallocate( diagg_dn )
#ENDIF
  deallocate( diagg_up )
end subroutine upgradeu
