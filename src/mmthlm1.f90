subroutine mmthlm1(a_up, a_dn)

  !	in  a 
  !	out a* (exp(-dtau*t2) * exp(-dtau*t1))^{-1}
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  use data_tmp
  implicit none

  !arguments:
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_up
  complex(dp), dimension(ndim,ndim), intent(inout) :: a_dn

#IFDEF BREAKUP_T
  ! local
  integer :: nf, i, j, i1, i2, i3, i4, ist

  if ( rt.gt.zero) then
      do nf = 1,2
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, i3, i4, j, v1, v2, v3, v4 )
!$OMP DO
         do i = 1,lq/4
            ist =  i + (nf - 1)*lq/4
            i1 = lthf(i,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            do j = 1,lq
               v1(j) = a_up(j,i1)*urtm1(ist,1,1) + a_up(j,i2)*urtm1(ist,2,1) + &
                    &  a_up(j,i3)*urtm1(ist,3,1) + a_up(j,i4)*urtm1(ist,4,1)
               v2(j) = a_up(j,i1)*urtm1(ist,1,2) + a_up(j,i2)*urtm1(ist,2,2) + &
                    &  a_up(j,i3)*urtm1(ist,3,2) + a_up(j,i4)*urtm1(ist,4,2)
               v3(j) = a_up(j,i1)*urtm1(ist,1,3) + a_up(j,i2)*urtm1(ist,2,3) + &
                    &  a_up(j,i3)*urtm1(ist,3,3) + a_up(j,i4)*urtm1(ist,4,3)
               v4(j) = a_up(j,i1)*urtm1(ist,1,4) + a_up(j,i2)*urtm1(ist,2,4) + &
                    &  a_up(j,i3)*urtm1(ist,3,4) + a_up(j,i4)*urtm1(ist,4,4)
            enddo
            do j = 1,ndim
               a_up(j,i1) = v1(j)
               a_up(j,i2) = v2(j)
               a_up(j,i3) = v3(j)
               a_up(j,i4) = v4(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo
  endif

#IFDEF SPINDOWN
  if ( rt.gt.zero) then
      do nf = 1,2
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, i3, i4, j, v1, v2, v3, v4 )
!$OMP DO
         do i = 1,lq/4
            ist =  i + (nf - 1)*lq/4
            i1 = lthf(i,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            do j = 1,lq
               v1(j) = a_dn(j,i1)*urtm1_dn(ist,1,1) + a_dn(j,i2)*urtm1_dn(ist,2,1) + &
                    &  a_dn(j,i3)*urtm1_dn(ist,3,1) + a_dn(j,i4)*urtm1_dn(ist,4,1)
               v2(j) = a_dn(j,i1)*urtm1_dn(ist,1,2) + a_dn(j,i2)*urtm1_dn(ist,2,2) + &
                    &  a_dn(j,i3)*urtm1_dn(ist,3,2) + a_dn(j,i4)*urtm1_dn(ist,4,2)
               v3(j) = a_dn(j,i1)*urtm1_dn(ist,1,3) + a_dn(j,i2)*urtm1_dn(ist,2,3) + &
                    &  a_dn(j,i3)*urtm1_dn(ist,3,3) + a_dn(j,i4)*urtm1_dn(ist,4,3)
               v4(j) = a_dn(j,i1)*urtm1_dn(ist,1,4) + a_dn(j,i2)*urtm1_dn(ist,2,4) + &
                    &  a_dn(j,i3)*urtm1_dn(ist,3,4) + a_dn(j,i4)*urtm1_dn(ist,4,4)
            enddo
            do j = 1,ndim
               a_dn(j,i1) = v1(j)
               a_dn(j,i2) = v2(j)
               a_dn(j,i3) = v3(j)
               a_dn(j,i4) = v4(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo
  endif
#ENDIF

#ELSE
  call zgemm('n','n',ndim,ndim,ndim,cone,a_up,ndim,urtm1,ndim,czero,Atmp,ndim)
  a_up(:,:) = Atmp(:,:)
#IFDEF SPINDOWN
  call zgemm('n','n',ndim,ndim,ndim,cone,a_dn,ndim,urtm1_dn,ndim,czero,Atmp,ndim)
  a_dn(:,:) = Atmp(:,:)
#ENDIF

#ENDIF
end subroutine mmthlm1
