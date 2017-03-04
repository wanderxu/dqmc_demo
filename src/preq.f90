subroutine preq
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  use obser
  implicit none

  include 'mpif.h'

  complex(dp) :: mpi_obs_bin(10)

  call mpi_reduce( obs_bin, mpi_obs_bin, 10, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  obs_bin(:) = mpi_obs_bin(:)

  if( irank .eq. 0 ) then
      ! calculate obs_bin
      obs_bin(:) = obs_bin(:) / dcmplx( isize * nobs )
      open (unit=90,file='ener1.bin',status='unknown', action="write", position="append")
#IFDEF TEST
      write(90, '(5(2e16.8,4x))') obs_bin(1)/dcmplx(dble(ndim),0.d0), obs_bin(2:5)
#ELSE
      write(90, '(5(e16.8,4x))') dble(obs_bin(1))/dble(ndim), dble(obs_bin(2:5))
#ENDIF
      close(90)
  end if
  call mpi_barrier( mpi_comm_world, ierr )
end subroutine preq
