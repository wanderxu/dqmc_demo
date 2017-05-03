subroutine prtau
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  use obser
  implicit none

  complex(dp) :: znorm, chiszsz_qwn_tmp, chijxjx_qwn_tmp, chi0jxjx_qwn_tmp, zexpiwtqr
  complex(dp), dimension(:), allocatable :: collect1
  complex(dp), dimension(:,:), allocatable :: collect2
  character(40) :: filek
  real(dp) :: qvec(2)
  integer :: n, iq, itau, imj, i, j

  interface
     subroutine fourier_trans_tau(gr,filek)
       complex(kind=8), dimension(:,:) :: gr
       character (40) :: filek
     end subroutine fourier_trans_tau
  end interface

  include 'mpif.h'
  znorm = cone / dcmplx( dble(nsweep), 0.d0 )
  gtau_up = znorm * gtau_up
#IFDEF SPINDOWN
  gtau_dn = znorm * gtau_dn
#ENDIF
  chiszsz = znorm * chiszsz
  chijxjx = znorm * chijxjx
  chijxjxaa = znorm * chijxjxaa
  chijxjxab = znorm * chijxjxab
  chijxjxba = znorm * chijxjxba
  chijxjxbb = znorm * chijxjxbb
  if( ltau ) then
      allocate(collect2(lq,ltrot))
      n = lq*ltrot
      collect2 = czero
      call mpi_reduce(gtau_up,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      gtau_up = collect2/dcmplx( dble(isize), 0.d0 )
#IFDEF SPINDOWN
      call mpi_reduce(gtau_dn,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      gtau_dn = collect2/dcmplx( dble(isize), 0.d0 )
#ENDIF

      call mpi_reduce(chiszsz,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chiszsz = collect2/dcmplx( dble(isize), 0.d0 )

      call mpi_reduce(chijxjx,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chijxjx = collect2/dcmplx( dble(isize), 0.d0 )

      call mpi_reduce(chijxjxaa,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chijxjxaa = collect2/dcmplx( dble(isize), 0.d0 )

      call mpi_reduce(chijxjxab,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chijxjxab = collect2/dcmplx( dble(isize), 0.d0 )

      call mpi_reduce(chijxjxba,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chijxjxba = collect2/dcmplx( dble(isize), 0.d0 )

      call mpi_reduce(chijxjxbb,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chijxjxbb = collect2/dcmplx( dble(isize), 0.d0 )

      deallocate(collect2)
  end if

  if (irank.eq.0) then

     if(ltau) then

         filek = "gtau_up.bin"
         call fourier_trans_tau(gtau_up,filek)
#IFDEF SPINDOWN
         filek = "gtau_dn.bin"
         call fourier_trans_tau(gtau_dn,filek)
#ENDIF
          open (unit=188,file='szsztaur_corrlt.bin',status='unknown', action="write", position="append")
          do itau = 1, ltrot
              do imj = 1, lq
                  write(188,'(e22.12)') dble(chiszsz(imj,itau))
              end do
          end do
          close(188)

          open (unit=222,file='chijxjxtaur.bin',status='unknown', action="write", position="append")
          open (unit=223,file='chijxjxaataur.bin',status='unknown', action="write", position="append")
          open (unit=224,file='chijxjxabtaur.bin',status='unknown', action="write", position="append")
          open (unit=225,file='chijxjxbataur.bin',status='unknown', action="write", position="append")
          open (unit=226,file='chijxjxbbtaur.bin',status='unknown', action="write", position="append")
          do itau = 1, ltrot
              do imj = 1, lq
                  write(222,'(e22.12)') dble(chijxjx(imj,itau))
                  write(223,'(e22.12)') dble(chijxjxaa(imj,itau))
                  write(224,'(e22.12)') dble(chijxjxab(imj,itau))
                  write(225,'(e22.12)') dble(chijxjxba(imj,itau))
                  write(226,'(e22.12)') dble(chijxjxbb(imj,itau))
              end do
          end do
          close(222)
          close(223)
          close(224)
          close(225)
          close(226)
     end if

  endif
end subroutine prtau

subroutine fourier_trans_tau(gr,filek)
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  implicit none
  complex(dp), dimension(:,:) :: gr
  integer :: imj, no, no1, no2, nt, nk
  character (40) :: filek
  real(dp) :: xk_p(2), aimj_p(2)
  complex(dp), allocatable , dimension(:,:) :: gk

  allocate (gk(lq,ltrot))

  gk = dcmplx(0.d0,0.d0)
  do imj = 1,lq
     aimj_p = dble(list(imj,1)*a1_p) + dble(list(imj,2)*a2_p)
!$OMP PARALLEL &
!$OMP PRIVATE ( nt, nk )
!$OMP DO
     do nt = 1,ltrot
        do nk = 1,lq
           gk(nk,nt) = gk(nk,nt) +  gr(imj,nt)/zexpiqr(imj,nk)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
  enddo
  gk = gk/dcmplx(dble(lq),0.d0)

  open (unit=20,file=filek,status='unknown', action="write", position="append")
  do nk = 1,lq
     xk_p = dble(listk(nk,1))*b1_p + dble(listk(nk,2))*b2_p
     write(20,*) xk_p(1), xk_p(2)
     do nt = 1,ltrot
         write(20,*) gk(nk,nt)
     enddo
  enddo
  close(20)
  deallocate (gk)
end subroutine fourier_trans_tau
