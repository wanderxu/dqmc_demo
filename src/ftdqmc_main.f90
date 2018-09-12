program ftdqmc_main
#ifdef _OPENMP
  USE OMP_LIB
#endif
#ifdef MPI
  use mpi
#endif
  use blockc
  use data_tmp
#ifdef CUMC
  use mod_cumulate
#endif
  use ftdqmc_core
  implicit none

  ! local
  integer :: nbc, nsw
  character (len = 24) :: date_time_string
  real(dp) :: start_time, end_time, time1, time2
#ifdef CAL_AUTO
  integer :: i, j, imj, n, totsz, nti, ntj
  integer, allocatable, dimension(:) :: totsz_bin
  integer, allocatable, dimension(:,:) :: jjcorr_Rtau, mpi_jjcorr_Rtau
  real(dp), allocatable, dimension(:,:) :: jjcorr_Rtau_real
#endif

#ifdef MPI
  call MPI_INIT(ierr)                             
  call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr) 
#else
  irank = 0
  isize = 1
#endif

#ifdef _OPENMP
  start_time = omp_get_wtime()
#else
  call cpu_time(start_time)
#endif

  if( irank.eq.0 ) then
      open( unit=fout, file='ftdqmc.out', status='unknown' )
  end if

#ifdef TEST
  write(fout,'(a,e32.15)')  ' zero  = ', zero
  write(fout,'(a,e32.15)')  ' pi    = ', pi
  write(fout,'(a,2e32.15)') ' czero = ', czero
  write(fout,'(a,2e32.15)') ' cone  = ', cone
#endif

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
#ifndef CAL_AUTO
  call allocate_obs
#endif

#ifdef CUMC
  call set_neighbor
  call initial_heff
#endif

#ifdef CAL_AUTO
#ifdef GEN_CONFC_LEARNING
  if( llocal .and. .not. lstglobal ) then
      allocate( totsz_bin(2*nsweep) )
  else
      allocate( totsz_bin(nsweep) )
  end if
#else
  allocate(jjcorr_Rtau(lq,ltrot/2+1))
  allocate(jjcorr_Rtau_real(lq,ltrot/2+1))
  allocate(mpi_jjcorr_Rtau(lq,ltrot/2+1))
#endif
#endif

  max_wrap_error = 0.d0
  if(ltau) xmax_dyn = 0.d0

  call ftdqmc_sweep_start_0b
  logweightf_old = dble( logweightf_up + logweightf_dn )*2.d0
  call ftdqmc_calculate_weights( logweights_old )
  weight_track = logweightf_old + logweights_old

  if( irank .eq. 0 ) then
      write(fout,'(a)') ' ftdqmc_sweep_start done '
  end if

  ! warnup
  if( lwarnup ) then
      ! set nwarnup
      !nwarnup = nint( beta ) + 3
      nwarnup = 300
      if(rhub.le.0.d0) nwarnup = 0
#ifdef TEST
      nwarnup = 0
#endif
      if( irank.eq.0 ) then
          write(fout,'(a,i8)') ' nwarnup = ', nwarnup
      end if
      do nsw = 1, nwarnup
          if(llocal) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure_equaltime=.false.)
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure_equaltime=.false., lmeasure_dyn=.false.)
          end if
          call ftdqmc_stglobal(lmeas=.false.)
      end do
      if(irank.eq.0) write(fout, '(a,e16.8)') 'after wanrup, max_wrap_error = ', max_wrap_error
      if(irank.eq.0 .and. ltau) write(fout,'(a,e16.8)')'after wanrup  xmax_dyn = ', xmax_dyn
      xmax_dyn = 0.d0 ! in warnup, xmax_dyn is not right, reset it here
  end if

#ifdef _OPENMP
  time1 = omp_get_wtime()
#else
  call cpu_time(time1)
#endif
  do nbc =  1, nbin

#ifdef CAL_AUTO
#include "sweep_auto.f90"
#else
#include "sweep.f90"
#endif

      !!! --- Timming and outconfc
      if( nbc .eq. 1 )  then
#ifdef _OPENMP
          time2 = omp_get_wtime()
#else
          call cpu_time(time2)
#endif
          if(irank.eq.0) then
              n_outconf_pace = nint( dble( 3600 * 12 ) / ( time2-time1 ) )
              if( n_outconf_pace .lt. 1 ) n_outconf_pace = 1
              write(fout,'(a,e16.8,a)') ' time for 1 bin: ', time2-time1, ' s'
              write(fout,'(a,i12)') ' n_out_conf_pace = ', n_outconf_pace
          end if
#ifdef MPI
          call mpi_bcast( n_outconf_pace, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr )
#endif
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

#ifdef MPI
  call mpi_reduce(main_obs, mpi_main_obs, size(main_obs), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
#endif
  if(irank.eq.0) then
      if(lwrapu)  write(fout,'(a,e16.8)') ' >>> accep_u  = ', dble(main_obs(1))/aimag(main_obs(1))
      if(lwrapj)  write(fout,'(a,e16.8)') ' >>> accep_j  = ', dble(main_obs(2))/aimag(main_obs(2))
      if(lstglobal) then
          write(fout,'(a,e16.8)') ' >>> accep_st = ', dble(main_obs(3))/aimag(main_obs(3))
          write(fout,'(a,e16.8)') ' >>> cluster_size = ', dble(main_obs(4))/aimag(main_obs(4))*dble(ltrot*lq)
      end if
  end if

#ifdef CAL_AUTO
#ifdef GEN_CONFC_LEARNING
  deallocate( totsz_bin )
#else
  deallocate(mpi_jjcorr_Rtau)
  deallocate(jjcorr_Rtau_real)
  deallocate(jjcorr_Rtau)
#endif
#endif

#ifdef CUMC
  call deallocate_cumulate
#endif

#ifndef CAL_AUTO
  call deallocate_obs
#endif

  call deallocate_core
  call deallocate_data_tmp

  call deallocate_tables

  if( irank.eq.0 ) then
#ifdef _OPENMP
      end_time = omp_get_wtime()
#else
      call cpu_time(end_time)
#endif
      call fdate( date_time_string )
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

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
#endif

end program ftdqmc_main
