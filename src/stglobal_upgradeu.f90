subroutine stglobal_upgradeu(ntau,ni,green_up,green_dn, ratiofi)

#ifdef _OPENMP
  USE OMP_LIB
#endif
  use spring
  use blockc
  use data_tmp

  implicit none

  !arguments
  integer,intent(in) :: ntau, ni
  complex(dp), intent(inout), dimension(ndim,ndim) :: green_up, green_dn
  real(dp), intent(out) :: ratiofi

  !local
  complex(dp) ::  ratioup, ratiodn, ratiotot, del44_up, del44_dn
  integer :: i4, nl, nl1, nl2, nrflip
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight

  accm  = 0.d0
  !do i4 = 1,lq
     i4 = ni
     nrflip = 1

     del44_up   =  delta_u_up( nsigl_u(i4,ntau), nrflip )
     ratioup = dcmplx(1.d0,0.d0) + del44_up * ( cone - green_up(i4,i4) )

#ifdef SPINDOWN
     del44_dn   =  delta_u_dn( nsigl_u(i4,ntau), nrflip )
     ratiodn = dcmplx(1.d0,0.d0) + del44_dn * ( cone - green_dn(i4,i4) )
#endif

#ifdef SPINDOWN
     ratiotot = (ratioup*ratiodn)*dconjg(ratioup*ratiodn) !* deta_u( nsigl_u(i4,ntau), nrflip )
#else
     ratiotot = ratioup*dconjg(ratioup) * deta_u( nsigl_u(i4,ntau), nrflip )
#endif

     !ratio_re = dgaml(nsigl_u(i4,ntau),nrflip ) * dble( ratiotot * phaseu )/    dble( phaseu )
     ratio_re = dble( ratiotot ) ! * dgaml(nsigl_u(i4,ntau),nrflip)
     ratiofi = ratio_re

#ifdef TEST
     write(fout,'(a,2e16.8)') 'in upgradeu, ratio_re = ', ratio_re
#endif

     ratio_re_abs = ratio_re
     if (ratio_re .lt. 0.d0 )  ratio_re_abs = - ratio_re 

     !!!random = spring_sfmt_stream()
     !!!if ( ratio_re_abs .gt. random ) then
     ! accept it with ratio = 1 temporarily

        ! upgrade the inverse.

        accm  = accm + 1.d0
            ! upgrade phaseu.
        weight = dsqrt(dble(ratiotot*dconjg(ratiotot)))
        !phaseu =  phaseu*ratiotot/dcmplx(weight,0.d0)

         
        ! update greep_up
        do nl = 1, ndim
            u1(nl) = green_up(nl,i4)/ratioup
            v1(nl) = green_up(i4,nl)
        end do
        v1(i4) = v1(i4) - cone  ! note the sign
        do nl = 1, ndim
            v1(nl) = del44_up * v1(nl)
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

        ! update greep_dn
        do nl = 1, ndim
            u1(nl) = green_dn(nl,i4)/ratiodn
            v1(nl) = green_dn(i4,nl)
        end do
        v1(i4) = v1(i4) - cone  ! note the sign
        do nl = 1, ndim
            v1(nl) = del44_dn * v1(nl)
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

        !	      flip:
        nsigl_u(i4,ntau) =  nflipl(nsigl_u(i4,ntau), nrflip)
        
     !!!endif
     ! write(50,*)
  !enddo
  !obs(27) = obs(27) + dcmplx(accm/dble(lq),0.d0)
  !obs(28) = obs(28) + dcmplx(1.d0,0.d0)
  !write(6,*) 'upgradec: acc: ', accm/dble(lq)
end subroutine stglobal_upgradeu
