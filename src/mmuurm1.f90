subroutine mmuurm1(a_up, a_dn, nf, ntau, nflag)
  !	in a out u(nf) * a                if nflag = 1
  !	in a out exp(d(nf)) * ut(nf) * a  if nflag = 2

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

  !local
  integer :: nl, i, j, nf1, nn, i1, i2
  complex (dp) :: ut(2,2), u(2,2)

  if (nflag.eq.3) then
!$OMP PARALLEL &
!$OMP PRIVATE ( i, nl )
!$OMP DO
     do nl= 1, ndim
        do i = 1, ndim
           a_up(i,nl) =  a_up(i,nl) / xsigma_u_up( nsigl_u(i,ntau))
#IFDEF SPINDOWN
           a_dn(i,nl) =  a_dn(i,nl) / xsigma_u_dn( nsigl_u(i,ntau))
#ENDIF
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
     return
  endif

  if (nf.gt.4) then ! current.
     nf1 = nf - 4
     do i = 1,2
        do j = 1,2
           u (i,j) = ur_j (i,j)
           ut(i,j) = urt_j(i,j)
        enddo
     enddo
  else             ! kinetic.
     nf1 = nf
     do i = 1,2
        do j = 1,2
           u (i,j) = ur_k (i,j)
           ut(i,j) = urt_k(i,j)
        enddo
     enddo
  endif
  if (nf1.eq.1) nn = 1
  if (nf1.eq.2) nn = 1
  if (nf1.eq.3) nn = 2
  if (nf1.eq.4) nn = 2
  
  if (nflag.eq.2) then
     do i = lq, 1, -1
        i1 = i
        i2 = i+lq
        do j = 1,ndim
           v1(j)   =  ut(1,1) * a_up(i1,j) + ut(1,2) * a_up(i2,j)
           v2(j)   =  ut(2,1) * a_up(i1,j) + ut(2,2) * a_up(i2,j) 
        enddo
        do j = 1,ndim
           a_up(i1,j) = v1(j)
           a_up(i2,j) = v2(j)
        enddo
        
     enddo
  endif
  
  if (nflag.eq.1) then
     do i = lq, 1, -1
        i1 = i
        i2 = i+lq
        if (nf.gt.4) then ! current.
           do j = 1,ndim
              a_up(i1,j) =  a_up(i1,j) /  xsigp2(nsigl_j(i1,nn,ntau))
              a_up(i2,j) =  a_up(i2,j) /  xsigm2(nsigl_j(i1,nn,ntau))
           enddo
        else              ! kinetic
           do j = 1,ndim
              a_up(i1,j) =  a_up(i1,j) / xsigp2(nsigl_k(i1,nn,ntau))
              a_up(i2,j) =  a_up(i2,j) / xsigm2(nsigl_k(i1,nn,ntau))
           enddo
        endif
        do j = 1,ndim
           v1(j)   =  u(1,1) * a_up(i1,j) +  u(1,2) * a_up(i2,j) 
           v2(j)   =  u(2,1) * a_up(i1,j) +  u(2,2) * a_up(i2,j) 
        enddo
        do j = 1,ndim
           a_up(i1,j) = v1(j)
           a_up(i2,j) = v2(j)
        enddo
     enddo
  endif
end subroutine mmuurm1
