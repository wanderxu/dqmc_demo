subroutine mmthrm1(a_up,a_dn)

  !     in  a 
  !     out exp( dtau*t3) * exp( dtau*t4) * exp( dtau*t5) * exp( dtau*t6) *
  !	    exp( dtau*t1) * exp( dtau*t2) * a
  !     t1 and t2 is 4x4 matrix, t3, t4, t5 and t6 is 2x2 matrix

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
  integer :: nf_tmp, nf, i, j, i1, i2, i3, i4, ist

  if (rt.gt.zero) then
      do nf = 2,1,-1
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, i3, i4, j, v1, v2, v3, v4 )
!$OMP DO
         do i = 1,lq/4
            ist = i + (nf - 1)*lq/4 
            i1 = lthf(i,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            
            do j = 1,lq
               v1(j) = urtm1(ist,1,1)*a_up(i1,j) + urtm1(ist,1,2)*a_up(i2,j) + &
                    &  urtm1(ist,1,3)*a_up(i3,j) + urtm1(ist,1,4)*a_up(i4,j)
               v2(j) = urtm1(ist,2,1)*a_up(i1,j) + urtm1(ist,2,2)*a_up(i2,j) + &
                    &  urtm1(ist,2,3)*a_up(i3,j) + urtm1(ist,2,4)*a_up(i4,j)
               v3(j) = urtm1(ist,3,1)*a_up(i1,j) + urtm1(ist,3,2)*a_up(i2,j) + &
                    &  urtm1(ist,3,3)*a_up(i3,j) + urtm1(ist,3,4)*a_up(i4,j)
               v4(j) = urtm1(ist,4,1)*a_up(i1,j) + urtm1(ist,4,2)*a_up(i2,j) + &
                    &  urtm1(ist,4,3)*a_up(i3,j) + urtm1(ist,4,4)*a_up(i4,j)
            enddo
            do j = 1, lq
               a_up(i1,j) = v1(j)
               a_up(i2,j) = v2(j)
               a_up(i3,j) = v3(j)
               a_up(i4,j) = v4(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo

      do nf = 8,1,-1
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, i3, i4, j, v1, v2, v3, v4 )
!$OMP DO
         do i = 1,lq/16
            ist = i + (nf - 1)*lq/16
            i1 = lthf3(i,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            i2 = nnlist(i2,1)
            i3 = nnlist(i3,5)
            i4 = nnlist(i4,2)
            
            do j = 1,lq
               v1(j) = urtdm1(ist,1,1)*a_up(i1,j) + urtdm1(ist,1,2)*a_up(i2,j) + &
                    &  urtdm1(ist,1,3)*a_up(i3,j) + urtdm1(ist,1,4)*a_up(i4,j)
               v2(j) = urtdm1(ist,2,1)*a_up(i1,j) + urtdm1(ist,2,2)*a_up(i2,j) + &
                    &  urtdm1(ist,2,3)*a_up(i3,j) + urtdm1(ist,2,4)*a_up(i4,j)
               v3(j) = urtdm1(ist,3,1)*a_up(i1,j) + urtdm1(ist,3,2)*a_up(i2,j) + &
                    &  urtdm1(ist,3,3)*a_up(i3,j) + urtdm1(ist,3,4)*a_up(i4,j)
               v4(j) = urtdm1(ist,4,1)*a_up(i1,j) + urtdm1(ist,4,2)*a_up(i2,j) + &
                    &  urtdm1(ist,4,3)*a_up(i3,j) + urtdm1(ist,4,4)*a_up(i4,j)
            enddo
            do j = 1, lq
               a_up(i1,j) = v1(j)
               a_up(i2,j) = v2(j)
               a_up(i3,j) = v3(j)
               a_up(i4,j) = v4(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo

      do nf = 4,1,-1
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, nf_tmp, i1, i2, i3, i4, j, v1, v2 )
!$OMP DO
         do i = 1, lq/4
            nf_tmp = (nf-1)/2 + 1
            ist = i + (nf - 1)*lq/4
            i1 = lthf2(i,nf_tmp)
            i2 = nnlist(i1,5)
            i3 = nnlist(i1,1)
            i4 = nnlist(i3,6)

            if ( mod(nf,2) .ne. 0 ) then
                 i1 = i3
                 i2 = i4
            end if
            
            do j = 1,lq
               v1(j) = urtcm1(ist,1,1) * a_up(i1,j) + urtcm1(ist,1,2) * a_up(i2,j)
               v2(j) = urtcm1(ist,2,1) * a_up(i1,j) + urtcm1(ist,2,2) * a_up(i2,j)
            enddo
            do j = 1, lq
               a_up(i1,j) = v1(j)
               a_up(i2,j) = v2(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo
  endif

#IFDEF SPINDOWN
  if (rt.gt.zero) then
      do nf = 2,1,-1
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, i3, i4, j, v1, v2, v3, v4 )
!$OMP DO
         do i = 1,lq/4
            ist = i + (nf - 1)*lq/4 
            i1 = lthf(i,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            
            do j = 1,lq
               v1(j) = urtm1_dn(ist,1,1)*a_dn(i1,j) + urtm1_dn(ist,1,2)*a_dn(i2,j) + &
                    &  urtm1_dn(ist,1,3)*a_dn(i3,j) + urtm1_dn(ist,1,4)*a_dn(i4,j)
               v2(j) = urtm1_dn(ist,2,1)*a_dn(i1,j) + urtm1_dn(ist,2,2)*a_dn(i2,j) + &
                    &  urtm1_dn(ist,2,3)*a_dn(i3,j) + urtm1_dn(ist,2,4)*a_dn(i4,j)
               v3(j) = urtm1_dn(ist,3,1)*a_dn(i1,j) + urtm1_dn(ist,3,2)*a_dn(i2,j) + &
                    &  urtm1_dn(ist,3,3)*a_dn(i3,j) + urtm1_dn(ist,3,4)*a_dn(i4,j)
               v4(j) = urtm1_dn(ist,4,1)*a_dn(i1,j) + urtm1_dn(ist,4,2)*a_dn(i2,j) + &
                    &  urtm1_dn(ist,4,3)*a_dn(i3,j) + urtm1_dn(ist,4,4)*a_dn(i4,j)
            enddo
            do j = 1, lq
               a_dn(i1,j) = v1(j)
               a_dn(i2,j) = v2(j)
               a_dn(i3,j) = v3(j)
               a_dn(i4,j) = v4(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo

      do nf = 8,1,-1
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, i1, i2, i3, i4, j, v1, v2, v3, v4 )
!$OMP DO
         do i = 1,lq/16
            ist = i + (nf - 1)*lq/16
            i1 = lthf3(i,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            i2 = nnlist(i2,1)
            i3 = nnlist(i3,5)
            i4 = nnlist(i4,2)
            
            do j = 1,lq
               v1(j) = urtdm1_dn(ist,1,1)*a_dn(i1,j) + urtdm1_dn(ist,1,2)*a_dn(i2,j) + &
                    &  urtdm1_dn(ist,1,3)*a_dn(i3,j) + urtdm1_dn(ist,1,4)*a_dn(i4,j)
               v2(j) = urtdm1_dn(ist,2,1)*a_dn(i1,j) + urtdm1_dn(ist,2,2)*a_dn(i2,j) + &
                    &  urtdm1_dn(ist,2,3)*a_dn(i3,j) + urtdm1_dn(ist,2,4)*a_dn(i4,j)
               v3(j) = urtdm1_dn(ist,3,1)*a_dn(i1,j) + urtdm1_dn(ist,3,2)*a_dn(i2,j) + &
                    &  urtdm1_dn(ist,3,3)*a_dn(i3,j) + urtdm1_dn(ist,3,4)*a_dn(i4,j)
               v4(j) = urtdm1_dn(ist,4,1)*a_dn(i1,j) + urtdm1_dn(ist,4,2)*a_dn(i2,j) + &
                    &  urtdm1_dn(ist,4,3)*a_dn(i3,j) + urtdm1_dn(ist,4,4)*a_dn(i4,j)
            enddo
            do j = 1, lq
               a_dn(i1,j) = v1(j)
               a_dn(i2,j) = v2(j)
               a_dn(i3,j) = v3(j)
               a_dn(i4,j) = v4(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo

      do nf = 4,1,-1
!$OMP PARALLEL &
!$OMP PRIVATE ( i, ist, nf_tmp, i1, i2, i3, i4, j, v1, v2 )
!$OMP DO
         do i = 1, lq/4
            nf_tmp = (nf-1)/2 + 1
            ist = i + (nf - 1)*lq/4
            i1 = lthf2(i,nf_tmp)
            i2 = nnlist(i1,5)
            i3 = nnlist(i1,1)
            i4 = nnlist(i3,6)

            if ( mod(nf,2) .ne. 0 ) then
                 i1 = i3
                 i2 = i4
            end if
            
            do j = 1,lq
               v1(j) = urtcm1_dn(ist,1,1) * a_dn(i1,j) + urtcm1_dn(ist,1,2) * a_dn(i2,j)
               v2(j) = urtcm1_dn(ist,2,1) * a_dn(i1,j) + urtcm1_dn(ist,2,2) * a_dn(i2,j)
            enddo
            do j = 1, lq
               a_dn(i1,j) = v1(j)
               a_dn(i2,j) = v2(j)
            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo
  endif
#ENDIF

#ELSE
  call zgemm('n','n',ndim,ndim,ndim,cone,urtm1,ndim,a_up,ndim,czero,Atmp,ndim)
  a_up(:,:) = Atmp(:,:)

#IFDEF SPINDOWN
  call zgemm('n','n',ndim,ndim,ndim,cone,urtm1_dn,ndim,a_dn,ndim,czero,Atmp,ndim)
  a_dn(:,:) = Atmp(:,:)
#ENDIF

#ENDIF

end subroutine mmthrm1
