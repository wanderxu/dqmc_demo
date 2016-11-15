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
  !integer status(mpi_status_size)
  !write(6,*) ' in prtau : ', nobst
  !!!znorm = dcmplx(1.d0,0.d0) / dcmplx( dble(nobst), 0.d0 )
  znorm = cone / dcmplx( dble(nsweep), 0.d0 )
  gtau_up = znorm * gtau_up
#IFDEF SPINDOWN
  gtau_dn = znorm * gtau_dn
#ENDIF
  chiszsz = znorm * chiszsz
  chijxjx = znorm * chijxjx
!!!  jttbin  = znorm * jttbin
!!!  j00bin  = znorm * j00bin
  if( ltauall ) then
      allocate(collect2(lq,ltrot))
      n = lq*ltrot
      collect2 = czero
      call mpi_reduce(gtau_up,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      !!!green_tau = collect2/dcmplx(dble(isize),0.d0)
      gtau_up = collect2/dcmplx( dble(isize), 0.d0 )
#IFDEF SPINDOWN
      call mpi_reduce(gtau_dn,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      gtau_dn = collect2/dcmplx( dble(isize), 0.d0 )
#ENDIF

      call mpi_reduce(chiszsz,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chiszsz = collect2/dcmplx( dble(isize), 0.d0 )

      call mpi_reduce(chijxjx,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      chijxjx = collect2/dcmplx( dble(isize), 0.d0 )

!!!      call mpi_reduce(jttbin,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
!!!      jttbin = collect2/dcmplx( dble(isize), 0.d0 )
!!!
!!!      call mpi_reduce(j00bin,collect2,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
!!!      j00bin = collect2/dcmplx( dble(isize), 0.d0 )

      deallocate(collect2)
  else 
      allocate( collect1(lq) )
      n = lq
      collect1 = czero
      call mpi_reduce(gtau_up(:,ltrot/2),collect1,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      !!!green_tau = collect1/dcmplx(dble(isize),0.d0)
      gtau_up(:,ltrot/2) = collect1(:)/dcmplx( dble(isize), 0.d0 )
#IFDEF SPINDOWN
      call mpi_reduce(gtau_dn(:,ltrot/2),collect1,n,mpi_complex16,mpi_sum,0,mpi_comm_world,ierr)
      gtau_dn(:,ltrot/2) = collect1(:)/dcmplx( dble(isize), 0.d0 )
      deallocate(collect1)
#ENDIF
  end if

  if (irank.eq.0) then

     if(ltauall) then

         filek = "gtau_up.bin"
         call fourier_trans_tau(gtau_up,filek)
#IFDEF SPINDOWN
         filek = "gtau_dn.bin"
         call fourier_trans_tau(gtau_dn,filek)
#ENDIF

          open(unit=177,file='chiszsz.bin',status='unknown', action="write", position="append")
          do iq = 1, lq
              qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
              do n = 0, nuse
                  chiszsz_qwn_tmp = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( itau, imj, zexpiwtqr )
!$OMP DO REDUCTION ( + : chiszsz_qwn_tmp )
                  do itau = 1, ltrot
                      do imj = 1, lq
                          zexpiwtqr = zexpiwt(itau-1,n) / zexpiqr(imj,iq)
                          chiszsz_qwn_tmp = chiszsz_qwn_tmp + chiszsz(imj,itau) * zexpiwtqr
                      end do
                  end do
!$OMP END DO
!$OMP END PARALLEL
                  !chiszsz_qwn(n,iq) = chiszsz_qwn_tmp
                  write(177, '(2f16.8,f16.8,2e16.8)') qvec(:), wn(n), chiszsz_qwn_tmp*dcmplx(dtau/dble(lq),0.d0)
              end do
          end do
          close(177)

          open (unit=188,file='szsztaur_corrlt.bin',status='unknown', action="write", position="append")
          do itau = 1, ltrot
              do imj = 1, lq
                  write(188,'(e22.12)') dble(chiszsz(imj,itau))
              end do
          end do
          close(188)

          open(unit=199,file='chijxjx.bin',status='unknown', action="write", position="append")
          do iq = 1, lq
              qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
              !do n = 0, nuse ! only need zero frequency needed
                  chijxjx_qwn_tmp = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( itau, imj, zexpiwtqr )
!$OMP DO REDUCTION ( + : chijxjx_qwn_tmp )
                  do itau = 1, ltrot
                      do imj = 1, lq
                          !zexpiwtqr = zexpiwt(itau-1,n) / zexpiqr(imj,iq)
                          zexpiwtqr = cone / zexpiqr(imj,iq)
                          chijxjx_qwn_tmp = chijxjx_qwn_tmp + chijxjx(imj,itau) * zexpiwtqr
                      end do
                  end do
!$OMP END DO
!$OMP END PARALLEL
                  !chijxjx_qwn(n,iq) = chijxjx_qwn_tmp
                  write(199, '(2f16.8,f16.8,2e16.8)') qvec(:), wn(0), chijxjx_qwn_tmp*dcmplx(dtau/dble(lq),0.d0)
              !end do
          end do
          close(199)

          open (unit=222,file='chijxjxtaur.bin',status='unknown', action="write", position="append")
          do itau = 1, ltrot
              do imj = 1, lq
                  write(222,'(e22.12)') dble(chijxjx(imj,itau))
              end do
          end do
          close(222)

!!!          open(unit=299,file='chi0jxjx.bin',status='unknown', action="write", position="append")
!!!          do iq = 1, lq
!!!              qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
!!!              !do n = 0, nuse ! only need zero frequency needed
!!!                  chi0jxjx_qwn_tmp = czero
!!!!$OMP PARALLEL &
!!!!$OMP PRIVATE ( itau, j, i, imj, zexpiwtqr )
!!!!$OMP DO REDUCTION ( + : chi0jxjx_qwn_tmp )
!!!                  do itau = 1, ltrot
!!!                      do j = 1, lq
!!!                        do i = 1, lq
!!!                          imj = latt_imj(i,j)
!!!                          !zexpiwtqr = zexpiwt(itau-1,n) / zexpiqr(imj,iq)
!!!                          zexpiwtqr = cone / zexpiqr(imj,iq)
!!!                          chi0jxjx_qwn_tmp = chi0jxjx_qwn_tmp - jttbin(i,itau)*j00bin(j,itau) * zexpiwtqr  ! (-1) come from (it)^2
!!!                        end do
!!!                      end do
!!!                  end do
!!!!$OMP END DO
!!!!$OMP END PARALLEL
!!!                  !chijxjx_qwn(n,iq) = chijxjx_qwn_tmp
!!!                  write(299, '(2f16.8,f16.8,2e16.8)') qvec(:), wn(0), chi0jxjx_qwn_tmp*dcmplx(dtau/dble(lq*lq),0.d0)
!!!              !end do
!!!          end do
!!!          close(299)

     else

         filek = "gtau_up_halfbeta.bin"
         call fourier_trans_tau(gtau_up,filek)
#IFDEF SPINDOWN
         filek = "gtau_dn_halfbeta.bin"
         call fourier_trans_tau(gtau_dn,filek)
#ENDIF

     end if

  endif
     !!!filek = "diractau_tot"
     !!!call dirac_point(spinup_tau , spindo_tau ,filek)
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

  if(ltauall) then

  gk = dcmplx(0.d0,0.d0)
  do imj = 1,lq
     aimj_p = dble(list(imj,1)*a1_p) + dble(list(imj,2)*a2_p)
!$OMP PARALLEL &
!$OMP PRIVATE ( nt, nk )
!$OMP DO
     do nt = 1,ltrot
        do nk = 1,lq
           !xk_p = dble(listk(nk,1))*b1_p + dble(listk(nk,2))*b2_p
           !gk(nk,nt) = gk(nk,nt) + exp( dcmplx( 0.d0, aimj_p(1)*xk_p(1) + aimj_p(2)*xk_p(2)) ) * gr(imj,nt)
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
     !! convert ffa convention to correct convention
     write(20,*) xk_p(1), xk_p(2)
     do nt = 1,ltrot
         write(20,*) gk(nk,nt)
     enddo
  enddo
  close(20)

  else
      gk = dcmplx(0.d0,0.d0)
      do imj = 1,lq
         aimj_p = dble(list(imj,1)*a1_p) + dble(list(imj,2)*a2_p)
         !do nt = 1,ltrot
         nt = ltrot/2
            do nk = 1,lq
               !xk_p = dble(listk(nk,1))*b1_p + dble(listk(nk,2))*b2_p
               !gk(nk,nt) = gk(nk,nt) + exp( dcmplx( 0.d0, aimj_p(1)*xk_p(1) + aimj_p(2)*xk_p(2)) ) * gr(imj,nt)
               gk(nk,nt) = gk(nk,nt) + gr(imj,nt)/zexpiqr(imj,nk) 
            enddo
         !enddo
      enddo
      gk(:,ltrot/2) = gk(:,ltrot/2)/dcmplx(dble(lq),0.d0)

      open (unit=20,file=filek,status='unknown', action="write", position="append")
      do nk = 1,lq
         xk_p = dble(listk(nk,1))*b1_p + dble(listk(nk,2))*b2_p
         !! convert ffa convention to correct convention
         write(20,*) xk_p(1), xk_p(2)
         !do nt = 1,ltrot
         nt = ltrot/2
             write(20,*) gk(nk,nt)
         !enddo
      enddo
      close(20)
  end if

  deallocate (gk)
end subroutine fourier_trans_tau
