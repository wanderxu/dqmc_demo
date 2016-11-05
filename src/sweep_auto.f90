      totsz_bin(:) = 0
      do nsw = 1, nsweep
          if(lstglobal .and. llocal ) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.false.)
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.false.)
              call ftdqmc_stglobal(lmeas=.false.)
              totsz=0
!$OMP PARALLEL &
!$OMP PRIVATE ( n, i )
!$OMP DO REDUCTION ( + : totsz )
              do n = 1, ltrot
                  do i = 1, ndim
                      totsz = totsz + nsigl_u(i,n)
                  end do
              end do
!$OMP END DO
!$OMP END PARALLEL
              totsz_bin(nsw) = totsz
          else if( lstglobal ) then
              call ftdqmc_stglobal(lmeas=.false.)
              totsz=0
!$OMP PARALLEL &
!$OMP PRIVATE ( n, i )
!$OMP DO REDUCTION ( + : totsz )
              do n = 1, ltrot
                  do i = 1, ndim
                      totsz = totsz + nsigl_u(i,n)
                  end do
              end do
!$OMP END DO
!$OMP END PARALLEL
              totsz_bin(nsw) = totsz
          else if ( llocal ) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.false.)
              totsz=0
!$OMP PARALLEL &
!$OMP PRIVATE ( n, i )
!$OMP DO REDUCTION ( + : totsz )
              do n = 1, ltrot
                  do i = 1, ndim
                      totsz = totsz + nsigl_u(i,n)
                  end do
              end do
!$OMP END DO
!$OMP END PARALLEL
              totsz_bin(nsw*2-1) = totsz
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.false.)
              totsz=0
!$OMP PARALLEL &
!$OMP PRIVATE ( n, i )
!$OMP DO REDUCTION ( + : totsz )
              do n = 1, ltrot
                  do i = 1, ndim
                      totsz = totsz + nsigl_u(i,n)
                  end do
              end do
!$OMP END DO
!$OMP END PARALLEL
              totsz_bin(nsw*2) = totsz
          else
              stop ' lstglobal and llocal should  not both false '
          end if
#IFDEF TEST
          if( irank .eq. 0 ) then
              write(fout,'(a,i4,i4,a)') ' ftdqmc_sweep ', nbc, nsw,  '  done'
          end if
#ENDIF
      end do
      if( llocal .and. .not. lstglobal ) then
          open (unit=9091,file='totsz.bin',status='unknown', action="write", position="append")
          do i = 1, 2*nsweep
            write(9091, '(e16.8)') dble(abs(totsz_bin(i)))/dble(ltrot*lq)
          end do
          close(9091)
      else
          open (unit=9091,file='totsz.bin',status='unknown', action="write", position="append")
          do i = 1, nsweep
            write(9091, '(e16.8)') dble(abs(totsz_bin(i)))/dble(ltrot*lq)
          end do
          close(9091)
      end if
