subroutine ftdqmc_initial
  use spring
  use blockc

  integer :: system_time
  integer :: stream_seed
  character (len = 20) :: date_time_string

  !================================================
  !%% inital the pseudo random number generator   $
  !------------------------------------------------
  call system_clock(system_time)
  stream_seed = abs( system_time - ( irank * 1981 + 2008 ) * 951049 )
#IFDEF TEST
  stream_seed = abs( 0 - ( irank * 1981 + 2008 ) * 951049 )
  write(fout, '(a,i20)') ' stream_seed = ', stream_seed
#ENDIF
  call spring_sfmt_init(stream_seed)

  call s_time_builder(date_time_string)

  ! print head
  if(irank.eq.0) then

  write(fout,'(a)') ' ===================================================================================='
  write(fout,*)
  write(fout,'(a)') '        The finite temperature determinant quantum monte carlo (DQMC) package '
  write(fout,*)
  write(fout,'(a)') '            FFFF   TTTTT   DDD      QQQ     M   M      CCCC                    '
  write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
  write(fout,'(a)') '            FFFF     T     D   D   Q   Q   M M M M    C                        '
  write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
  write(fout,'(a)') '            F        T     DDD      QQQ   M   M   M    CCCC                    '
  write(fout,'(a)') '                                      \                                       '
  write(fout,*)
  write(fout,*)
  write(fout,'(a)') ' written by Xiao Yan Xu ( wanderxu@gmail.com )                                '
  write(fout,*)
  write(fout,'(a)') ' history: '
  write(fout,*)
  write(fout,'(a)') '     22/01/2016,  version 1.0  '
  write(fout,*)
  write(fout,'(a)') ' ------------------------------------------------------------------------------------'
  write(fout,*)
  write(fout,'(a)') ' >>> The simulation start running at '//date_time_string
  if( isize .gt. 1 ) then
      write(fout,'(a,i6,a)') ' >>> Parallelism running with', isize, '  processes'
  else
      write(fout,'(a)') ' >>> Serial running '
  end if

  end if
end subroutine ftdqmc_initial

subroutine ftdqmc_initial_print
  use blockc
  implicit none

  integer :: i, j

  IF(irank.eq.0) THEN

  write(fout,*)
  write(fout,'(a)')' --------------------- '
  write(fout,'(a)')' The input parameters  '
  write(fout,'(a)')' --------------------- '
  write(fout,*)
  write(fout,'(a,f6.2)')    ' t      = ', rt
  write(fout,'(a,f6.2)')    ' t_2    = ', rt2
  write(fout,'(a,f6.2)')    ' t_3    = ', rt3
  write(fout,'(a,f6.2)')    ' U      = ', rhub
  write(fout,'(a,f6.2)')    ' js     = ', js
  write(fout,'(a,f6.2)')    ' rj     = ', rj
  write(fout,'(a,f6.2)')    ' hx     = ', hx
  write(fout,'(a,f6.2)')    ' B      = ', xmag
  write(fout,'(a,f6.2)')    ' dimer  = ', dimer
  write(fout,'(a,f8.5)')    ' flux_x = ', flux_x
  write(fout,'(a,f8.5)')    ' flux_y = ', flux_y
  write(fout,'(a,i4)')      ' L      = ', l
  write(fout,'(a,i4)')      ' LQ     = ', lq
  write(fout,'(a,i4)')      ' NE     = ', ne
  write(fout,'(a,f6.2)')    ' beta   = ', beta
  write(fout,'(a,f7.3)')    ' dtau   = ', dtau
  write(fout,'(a,f7.3)')    ' mu     = ', mu
  write(fout,'(a,i6)')      ' ltrot  = ', ltrot
  write(fout,'(a,i6)')      ' nwrap  = ', nwrap
  write(fout,'(a,i6)')      ' nsweep = ', nsweep
  write(fout,'(a,i6)')      ' nbin   = ', nbin
  write(fout,'(a,i6)')      ' nst    = ', nst
  write(fout,'(a,i6)')      ' nsw_stglobal = ', nsw_stglobal
  write(fout,'(a,i6)')      ' nublock = ', nublock
  write(fout,*)  ' lwrapu = ', lwrapu
  write(fout,*)  ' lwrapj = ', lwrapj
  write(fout,*)  ' llocal = ', llocal
  write(fout,*)  ' lstglobal = ', lstglobal
  write(fout,*)  ' lsstau = ', lsstau
  write(fout,*)  ' ltau = ', ltau
  write(fout,*)  ' ltauall = ', ltauall
  write(fout,*) ' number of meas for 1 bin = ', nmeas_bin
  write(fout,*) ' total number of measurements = ', nbin*nmeas_bin

  write(fout,*)
  write(fout,'(a)')' --------------------- '
  write(fout,'(a)')' wrapping coordinates  '
  write(fout,'(a)')' --------------------- '
  write(fout,'(a)')'       wrap_step(1,i)   wrap_step(2,i)   iwrap_nt(nt) '
  do i = 1, nst
      write( fout, '(3i16)') wrap_step(1,i), wrap_step(2,i), iwrap_nt( wrap_step(2,i) )
  end do

  write(fout,*)
  write(fout,'(a)')' --------------------- '
  write(fout,'(a)')' The lattice sites list '
  write(fout,'(a)')' --------------------- '
  write(fout,'(a)') '   i     list(i,:) '
  do i = 1, lq
      write(fout,'(i6,2i4)') i,  list(i,:)
  end do

  write(fout, *)
  write(fout,'(a)') '-------------------------'
  write(fout,'(a)') ' imj distance info       '
  write(fout,'(a)') '-------------------------'
  write(fout, '(a)') '  site      distance    deg '
  do i = 1, lq
      j = distance_index(i)
      write(fout, '(i5,f16.8,i5)')j, distance_len(i), imjdeg(j)
  end do
  write(fout, *)

  write(fout, *)
  write(fout,'(a)') '----------------------------'
  write(fout,'(a)') ' irreducible distance info       '
  write(fout,'(a)') '----------------------------'
  write(fout, '(a)') '  index      distance    deg '
  do i = 1, num_equ_distance
      write(fout, '(i5,f16.8,i5)')i, irre_distance_len(i), irre_distance_deg(i)
  end do
  write(fout, *)


  END IF

end subroutine ftdqmc_initial_print
