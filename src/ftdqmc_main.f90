program ftdqmc_main
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  use data_tmp
  use ftdqmc_core
  implicit none

  include 'mpif.h'

  ! local
  integer :: nbc, nsw
  character (len = 20) :: date_time_string
  real(dp) :: start_time, end_time, time1, time2

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

  max_wrap_error = 0.d0
  if(ltau) xmax_dyn = 0.d0

  call ftdqmc_sweep_start_0b

  if( irank .eq. 0 ) then
      write(fout,'(a)') ' ftdqmc_sweep_start done '
  end if

  ! warnup
  if( lwarnup ) then
      ! set nwarnup
      !nwarnup = nint( beta ) + 3
      nwarnup = nint( beta )*lq + 1
      if(rhub.le.0.d0) nwarnup = 0
      if( irank.eq.0 ) then
          write(fout,'(a,i8)') ' nwarnup = ', nwarnup
      end if
      do nsw = 1, nwarnup
          call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.false.)
          call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.false.)
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

      do nsw = 1, nsweep

          call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.true.)
#IFDEF GEN_CONFC_LEARNING
      ! output configuration for learning
      call outconfc_bin(weight_track)
      call preq
      call obser_init
#ENDIF
          call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.true.)
#IFDEF GEN_CONFC_LEARNING
       ! output configuration for learning
       call outconfc_bin(weight_track)
       call preq
       call obser_init
#ENDIF
#IFDEF TEST
          if( irank .eq. 0 ) then
              write(fout,'(a,i4,i4,a)') ' ftdqmc_sweep ', nbc, nsw,  '  done'
          end if
#ENDIF

      end do

#IFNDEF GEN_CONFC_LEARNING
      call preq  ! reduce
#ENDIF
      if(ltau) call prtau

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

  end do

  if(irank.eq.0) write(fout, '(a,e16.8)') ' max_wrap_error = ', max_wrap_error
  if(irank.eq.0 .and. ltau) write(fout,'(a,e16.8)')' >>> xmax_dyn = ', xmax_dyn

  call outconfc

  call mpi_reduce(main_obs, mpi_main_obs, size(main_obs), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  if(irank.eq.0) then
      if(lupdateu)  write(fout,'(a,e16.8)') ' >>> accep_u  = ', dble(main_obs(1))/aimag(main_obs(1))
      if(lupdatej)  write(fout,'(a,e16.8)') ' >>> accep_j  = ', dble(main_obs(2))/aimag(main_obs(2))
      if(lstglobal) write(fout,'(a,e16.8)') ' >>> accep_st = ', dble(main_obs(3))/aimag(main_obs(3))
  end if


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