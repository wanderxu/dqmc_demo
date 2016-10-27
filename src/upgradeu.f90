subroutine upgradeu(ntau, green_up, green_dn)

#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use spring
  use blockc
  use data_tmp

  implicit none

  !arguments
  integer,intent(in) :: ntau
  complex(dp), intent(inout), dimension(ndim,ndim) :: green_up, green_dn

  !local
  complex(dp) ::  ratioup, ratiodn, ratiotot, del44_up, del44_dn
  integer :: i4, nl, nl1, nl2, nrflip, nfb, id, ntm1, nta1
  real(dp) :: accm, ratio_re, ratio_re_abs, random

  accm  = 0.d0
  do i4 = 1,lq
     nrflip = 1

     del44_up   =  delta_u_up( nsigl_u(i4,ntau), nrflip )
     ratioup = dcmplx(1.d0,0.d0) + del44_up * ( cone - green_up(i4,i4) )

#IFDEF TEST
     write(fout,'(a,2e16.8)') 'in upgradeu, ratioup = ', ratioup
#ENDIF

#IFDEF SPINDOWN
     del44_dn   =  delta_u_dn( nsigl_u(i4,ntau), nrflip )
     ratiodn = dcmplx(1.d0,0.d0) + del44_dn * ( cone - green_dn(i4,i4) )
#IFDEF TEST
     write(fout,'(a,2e16.8)') 'in upgradeu, ratiodn = ', ratiodn
#ENDIF
#ENDIF

     ! get Ising part ratio
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

#IFDEF SPINDOWN
     ratiotot = (ratioup*ratiodn)*dconjg(ratioup*ratiodn) * wsxsz(id) !* deta_u( nsigl_u(i4,ntau), nrflip )
#ELSE
     ratiotot = ratioup*dconjg(ratioup) * wsxsz(id) ! * deta_u( nsigl_u(i4,ntau), nrflip )
#ENDIF

     !ratio_re = dgaml(nsigl_u(i4,ntau),nrflip ) * dble( ratiotot * phaseu )/    dble( phaseu )
     ratio_re = dble( ratiotot ) ! * dgaml(nsigl_u(i4,ntau),nrflip)

#IFDEF TEST
     write(fout,'(a,2e16.8)') 'in upgradeu, ratio_re = ', ratio_re
#ENDIF

     ratio_re_abs = ratio_re
     if (ratio_re .lt. 0.d0 )  ratio_re_abs = - ratio_re 

     random = spring_sfmt_stream()
     if ( ratio_re_abs .gt. random ) then
        accm  = accm + 1.d0
        weight_track = weight_track + dlog( ratio_re_abs )
        ! update greep_up
        do nl = 1, ndim
            u1(nl) = green_up(nl,i4)/ratioup
            v1(nl) = ( -Imat(nl,i4) + green_up(i4,nl) ) * del44_up
        end do
!$OMP PARALLEL &
!$OMP PRIVATE ( nl2, nl1 )
!$OMP DO
        do nl2 = 1,ndim
        do nl1 = 1,ndim
           green_up(nl1,nl2) = green_up(nl1,nl2) + u1(nl1)*v1(nl2)
        enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL

#IFDEF SPINDOWN
        ! update greep_dn
        do nl = 1, ndim
            u1(nl) = green_dn(nl,i4)/ratiodn
            v1(nl) = ( -Imat(nl,i4) + green_dn(i4,nl) ) * del44_dn
        end do
!$OMP PARALLEL &
!$OMP PRIVATE ( nl2, nl1 )
!$OMP DO
        do nl2 = 1,ndim
        do nl1 = 1,ndim
           green_dn(nl1,nl2) = green_dn(nl1,nl2) + u1(nl1)*v1(nl2)
        enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
#ENDIF

        ! flip
        nsigl_u(i4,ntau) =  nflipl(nsigl_u(i4,ntau), nrflip)
        
     endif
  enddo
  main_obs(1) = main_obs(1) + dcmplx( accm, dble(lq) )
end subroutine upgradeu
