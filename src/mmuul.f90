subroutine mmuul(a_up, a_dn, nf, ntau, nflag)

  !in a  out a* exp(d(nf)) * ut(nf)  if nflag = 1
  !in a  out a* u(nf)                if nflag = 2

#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  use data_tmp
  implicit none

  !arguments:
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_up
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_dn
  integer, intent(in) :: nf,ntau,nflag

  !	local
  integer :: nl, i, j, nf1, nn, i1, i2
  complex (dp) :: ut(2,2), u(2,2)

  if (nflag.eq.3) then
!$OMP PARALLEL &
!$OMP PRIVATE ( j, nl )
!$OMP DO
     do j  = 1,ndim
        do nl = 1, ndim
           a_up(nl,j) = a_up(nl,j) * xsigma_u_up(nsigl_u(j,ntau))
#IFDEF SPINDOWN
           a_dn(nl,j) = a_dn(nl,j) * xsigma_u_dn(nsigl_u(j,ntau))
#ENDIF
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
     return
  endif
  
  if (nf.gt.4) then
     !  current.
     nf1 = nf -4
     do i = 1,2
        do j = 1,2
           u(i,j)  = ur_j (i,j)
           ut(i,j) = urt_j(i,j)
        enddo
     enddo
  else
     !kinetic.
     nf1 = nf
     do i = 1,2
        do j = 1,2
           u(i,j)  = ur_k (i,j)
           ut(i,j) = urt_k(i,j)
        enddo
     enddo
  endif
  if (nf1.eq.1) nn = 1
  if (nf1.eq.2) nn = 1
  if (nf1.eq.3) nn = 2
  if (nf1.eq.4) nn = 2
  
  
  if ( nflag.eq.2 ) then
     do i = lq,1,-1
        i1 = i
        i2 = i+lq
        do j = 1,ndim
           v1(j)   =  a_up(j,i1) * u(1,1) + a_up(j,i2) * u(2,1) 
           v2(j)   =  a_up(j,i1) * u(1,2) + a_up(j,i2) * u(2,2) 
        enddo
        do j = 1,ndim
           a_up(j,i1) = v1(j)
           a_up(j,i2) = v2(j)
        enddo
     enddo
  endif
  
  
  if ( nflag.eq.1 ) then
     do i = lq,1,-1
        i1 = i
        i2 = i+lq
        if (nf.gt.4) then
           ! current.
           do j = 1,ndim
              a_up(j,i1) = xsigp2(nsigl_j(i1,nn,ntau)) * a_up(j,i1)
              a_up(j,i2) = xsigm2(nsigl_j(i1,nn,ntau)) * a_up(j,i2)
           enddo
        else
           ! kenitic
           do j = 1,ndim
              a_up(j,i1) = xsigp2(nsigl_k(i1,nn,ntau))*a_up(j,i1)
              a_up(j,i2) = xsigm2(nsigl_k(i1,nn,ntau))*a_up(j,i2)
           enddo
        endif
        do j = 1,ndim
           v1(j) = a_up(j,i1) * ut(1,1) +  a_up(j,i2) * ut(2,1)
           v2(j) = a_up(j,i1) * ut(1,2) +  a_up(j,i2) * ut(2,2)
        enddo
        do j = 1,ndim
           a_up(j,i1) = v1(j)
           a_up(j,i2) = v2(j)
        enddo
     enddo
  endif

end subroutine mmuul
