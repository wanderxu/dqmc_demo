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

#ifdef MPI
  call mpi_reduce( obs_bin, mpi_obs_bin, 10, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  obs_bin(:) = mpi_obs_bin(:)
  call mpi_reduce( pair_bin, mpi_pair_bin, 19, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  pair_bin(:) = mpi_pair_bin(:)
  call mpi_reduce( high_pair_bin, mpi_high_pair_bin, 4, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  high_pair_bin(:) = mpi_high_pair_bin(:)
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
#IFDEF TEST
      write(91, '(19(2e16.8,4x))') pair_bin(1:19)/dcmplx( dble(ndim), 0.d0 )
#ELSE
      write(91, '(19(e14.6,2x))') dble(pair_bin(1:19))/dble(ndim)
#ENDIF
      close(91)

      high_pair_bin(:) = high_pair_bin(:) / dcmplx( isize * nobs )
      open (unit=92,file='highpair.bin',status='unknown', action="write", position="append")
      write(92, '(4(e14.6,2x))') dble(high_pair_bin(1:4))/dble(ndim)
      close(92)
  end if
  call mpi_barrier( mpi_comm_world, ierr )
end subroutine preq
