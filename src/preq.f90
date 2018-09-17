subroutine preq
#ifdef _OPENMP
  USE OMP_LIB
#endif
#ifdef MPI
  use mpi
#endif
  use blockc
  use obser
  implicit none

  complex(dp) :: mpi_obs_bin(10), mpi_pair_bin(19), mpi_high_pair_bin(4)
  integer, allocatable, dimension(:,:) :: mpi_i2
  real(dp), allocatable, dimension(:,:) :: mpi_r2

  real(dp) :: qvec(2), rij(2)
  integer :: imj, iq, itau, n
  complex(dp) :: sq_ising_qwn_tmp, zexpiwtqr

  if(lsstau) allocate(mpi_i2(lq,ltrot) )
  if(lsstau) allocate(mpi_r2(lq,ltrot))

#ifdef MPI
  call mpi_reduce( obs_bin, mpi_obs_bin, 10, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  obs_bin(:) = mpi_obs_bin(:)
  call mpi_reduce( pair_bin, mpi_pair_bin, 19, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  pair_bin(:) = mpi_pair_bin(:)
  call mpi_reduce( high_pair_bin, mpi_high_pair_bin, 4, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  high_pair_bin(:) = mpi_high_pair_bin(:)
  if(lsstau) then
  call mpi_reduce( isingzztau_corrlt, mpi_i2, lq*ltrot, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr )
  isingzztau_corrlt(:,:) = mpi_i2(:,:)
  end if
#endif

  if( irank .eq. 0 ) then
      ! calculate obs_bin
      obs_bin(:) = obs_bin(:) / dcmplx( isize * nobs )
      obs_bin(8) = obs_bin(9) / (obs_bin(8)**2)
      open (unit=90,file='ener1.bin',status='unknown', action="write", position="append")
      !write(90, '(8(2e16.8,4x))') obs_bin(1)/dcmplx(dble(ndim),0.d0), obs_bin(2), obs_bin(8), obs_bin(3)/dcmplx(dble(ndim),0.d0), obs_bin(4), obs_bin(5), obs_bin(6), obs_bin(7)
      write(90, '(9(e16.8,2x))') dble(obs_bin(1))/dble(ndim), dble(obs_bin(2)), dble(obs_bin(8)), obs_bin(3)/dcmplx(dble(ndim),0.d0), dble(obs_bin(4:7))
      close(90)
      open (unit=91,file='pair.bin',status='unknown', action="write", position="append")
      pair_bin(:) = pair_bin(:) / dcmplx( isize * nobs )
#ifdef TEST
      write(91, '(19(2e16.8,4x))') pair_bin(1:19)/dcmplx( dble(ndim), 0.d0 )
#else
      write(91, '(19(e14.6,2x))') dble(pair_bin(1:19))/dble(ndim)
#endif
      close(91)

      high_pair_bin(:) = high_pair_bin(:) / dcmplx( isize * nobs )
      open (unit=92,file='highpair.bin',status='unknown', action="write", position="append")
      write(92, '(4(e14.6,2x))') dble(high_pair_bin(1:4))/dble(ndim)
      close(92)

      ! calculate Sising(q,iwn)
      if( lsstau ) then
          mpi_r2(:,:) = dble(isingzztau_corrlt(:,:)) / dble( isize * nobs )
          open (unit=85,file='isingzztau_corrlt.bin',status='unknown', action="write", position="append")
          do iq = 1, lq
              qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
              do n = 0, nuse
                  sq_ising_qwn_tmp = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( itau, imj, zexpiwtqr )
!$OMP DO REDUCTION ( + : sq_ising_qwn_tmp )
                  do itau = 1, ltrot
                      do imj = 1, lq
                          zexpiwtqr = zexpiwt(itau-1,n) / zexpiqr(imj,iq)
                          sq_ising_qwn_tmp = sq_ising_qwn_tmp + dcmplx( mpi_r2(imj,itau), 0.d0 ) * zexpiwtqr
                      end do
                  end do
!$OMP END DO
!$OMP END PARALLEL
                  !sq_ising_qwn(n,iq) = sq_ising_qwn_tmp
                  write(85, '(2f16.8,f16.8,2e16.8)') qvec(:), wn(n), sq_ising_qwn_tmp*dcmplx(dtau/dble(lq),0.d0)
              end do
          end do
          close(85)

          open (unit=89,file='sstaur_corrlt.bin',status='unknown', action="write", position="append")
          do itau = 1, ltrot
              do imj = 1, lq
                  write(89,'(e22.12)') mpi_r2(imj,itau)
              end do
          end do
          close(89)
      end if
  end if
#ifdef MPI
  call mpi_barrier( mpi_comm_world, ierr )
#endif
  if(lsstau) deallocate(mpi_r2)
  if(lsstau) deallocate(mpi_i2)
end subroutine preq
