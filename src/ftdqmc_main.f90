program ftdqmc_main
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  use data_tmp
#IFDEF CUMC
  use mod_cumulate
#ENDIF
  use ftdqmc_core
  implicit none

  include 'mpif.h'

  ! local
  integer :: nbc, nsw
  character (len = 20) :: date_time_string
  real(dp) :: start_time, end_time, time1, time2
#IFDEF CAL_AUTO
  integer :: i, n, totsz
  integer, allocatable, dimension(:) :: totsz_bin
#ENDIF

  call MPI_INIT(ierr)                             
  call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr) 

#IFDEF _OPENMP
  start_time = omp_get_wtime()
#ELSE
  call cpu_time(start_time)
#ENDIF

  if( irank.eq.0 ) then
      open( unit=fout, file='ftdqmc.out', status='unknown' )
  end if

#IFDEF TEST
  write(fout,'(a,e32.15)')  ' zero  = ', zero
  write(fout,'(a,e32.15)')  ' pi    = ', pi
  write(fout,'(a,2e32.15)') ' czero = ', czero
  write(fout,'(a,2e32.15)') ' cone  = ', cone
#ENDIF

  main_obs(:) = czero

  call ftdqmc_initial

  call make_tables
  call sli
  call sltpf

  call ftdqmc_initial_print

  ! prepare for the DQMC
  call salph
  call inconfc
  call sthop
  


  call allocate_data_tmp
  call allocate_core
  call allocate_obs

#IFDEF CUMC
  call set_neighbor
  call initial_heff
#ENDIF

#IFDEF CAL_AUTO
  if( llocal .and. .not. lstglobal ) then
      allocate( totsz_bin(2*nsweep) )
  else
      allocate( totsz_bin(nsweep) )
  end if
#ENDIF

  max_wrap_error = 0.d0
  if(ltau) xmax_dyn = 0.d0

  call ftdqmc_sweep_start_0b
  call ftdqmc_calculate_weight( logweightf_old, logweights_old )

  if( irank .eq. 0 ) then
      write(fout,'(a)') ' ftdqmc_sweep_start done '
  end if

  ! warnup
  if( lwarnup ) then
      ! set nwarnup
      !nwarnup = nint( beta ) + 3
      nwarnup = 0
      if(rhub.le.0.d0) nwarnup = 0
#IFDEF TEST
      nwarnup = 0
#ENDIF
      if( irank.eq.0 ) then
          write(fout,'(a,i8)') ' nwarnup = ', nwarnup
      end if
      do nsw = 1, nwarnup
          if(llocal) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.false.)
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.false.)
          end if
          call ftdqmc_stglobal(lmeas=.false.)
      end do
      if(irank.eq.0) write(fout, '(a,e16.8)') 'after wanrup, max_wrap_error = ', max_wrap_error
      if(irank.eq.0 .and. ltau) write(fout,'(a,e16.8)')'after wanrup  xmax_dyn = ', xmax_dyn
      xmax_dyn = 0.d0 ! in warnup, xmax_dyn is not right, reset it here
  end if

#IFDEF _OPENMP
  time1 = omp_get_wtime()
#ELSE
  call cpu_time(time1)
#ENDIF
  do nbc =  1, nbin

      call obser_init
#IFDEF CAL_AUTO
#include 'sweep_auto.f90'
#ELSE
#include 'sweep.f90'
#ENDIF

      !!! --- Timming and outconfc
      if( nbc .eq. 1 )  then
#IFDEF _OPENMP
          time2 = omp_get_wtime()
#ELSE
          call cpu_time(time2)
#ENDIF
          if(irank.eq.0) then
              n_outconf_pace = nint( dble( 3600 * 12 ) / ( time2-time1 ) )
              if( n_outconf_pace .lt. 1 ) n_outconf_pace = 1
              write(fout,'(a,e16.8,a)') ' time for 1 bin: ', time2-time1, ' s'
              write(fout,'(a,i12)') ' n_out_conf_pace = ', n_outconf_pace
          end if
          call mpi_bcast( n_outconf_pace, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr )
      end if

      if( n_outconf_pace .lt. nbin/3 ) then
          if( mod(nbc,n_outconf_pace) .eq. 0 ) then
              call outconfc
          end if
      else if( mod( nbc, max(nbin/3,1) ) .eq. 0 ) then
          call outconfc
      end if

      if( irank.eq.0 .and. mod(nbc,max(nbin/10,1) ).eq.0 ) then
          write( fout, '(i5,a,i5,a)' ) nbc, '  /', nbin, '   finished '
      end if
      !!! --- END Timming and outconfc
  end do

  if(irank.eq.0) write(fout, '(a,e16.8)') ' max_wrap_error = ', max_wrap_error
  if(irank.eq.0 .and. ltau) write(fout,'(a,e16.8)')' >>> xmax_dyn = ', xmax_dyn

  call outconfc

  call mpi_reduce(main_obs, mpi_main_obs, size(main_obs), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  if(irank.eq.0) then
      if(lwrapu)  write(fout,'(a,e16.8)') ' >>> accep_u  = ', dble(main_obs(1))/aimag(main_obs(1))
      if(lwrapj)  write(fout,'(a,e16.8)') ' >>> accep_j  = ', dble(main_obs(2))/aimag(main_obs(2))
      if(lstglobal) then
          write(fout,'(a,e16.8)') ' >>> accep_st = ', dble(main_obs(3))/aimag(main_obs(3))
          write(fout,'(a,e16.8)') ' >>> cluster_size = ', dble(main_obs(4))/aimag(main_obs(4))*dble(ltrot*lq)
      end if
  end if

#IFDEF CAL_AUTO
  deallocate( totsz_bin )
#ENDIF

#IFDEF CUMC
  call deallocate_cumulate
#ENDIF

  call deallocate_core
  call deallocate_data_tmp

  call deallocate_tables

  if( irank.eq.0 ) then
#IFDEF _OPENMP
      end_time = omp_get_wtime()
#ELSE
      call cpu_time(end_time)
#ENDIF
      call s_time_builder(date_time_string)
      write(fout,*)
      write(fout,'(a,f10.2,a)') ' >>> Total time spent:', end_time-start_time, 's'
      write(fout,'(a)') ' >>> Happy ending at '//date_time_string
      write(fout,*)
      write(fout,'(a)') ' The simulation done !!! '
      write(fout,*)
      write(fout,'(a)') '        o         o    '
      write(fout,'(a)') '       o o       o o   '
      write(fout,'(a)') '       o o       o o   '
      write(fout,'(a)') '        o         o    '
      write(fout,'(a)') '       o o       o o   '
      write(fout,'(a)') '       o o       o o   '
      write(fout,'(a)') '        o         o    '
  end if

  close(fout)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

end program ftdqmc_main
