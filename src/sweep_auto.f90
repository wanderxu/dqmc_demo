#IFDEF GEN_CONFC_LEARNING
      totsz_bin(:) = 0
#ELSE
      jjcorr_R(:) = 0
#ENDIF
      do nsw = 1, nsweep
          if(lstglobal .and. llocal ) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.false.)
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.false.)
              call ftdqmc_stglobal(lmeas=.false.)
#IFDEF GEN_CONFC_LEARNING
              call outconfc_bin(weight_track)
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
#ENDIF
          else if( lstglobal ) then
              call ftdqmc_stglobal(lmeas=.false.)
#IFDEF GEN_CONFC_LEARNING
              call outconfc_bin(weight_track)
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
#ENDIF
          else if ( llocal ) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.false.)
#IFDEF GEN_CONFC_LEARNING
              call outconfc_bin(weight_track)
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
#ENDIF
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.false.)
#IFDEF GEN_CONFC_LEARNING
              call outconfc_bin(weight_track)
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
#ENDIF
          else
              stop ' lstglobal and llocal should  not both false '
          end if
#IFDEF TEST
          if( irank .eq. 0 ) then
              write(fout,'(a,i4,i4,a)') ' ftdqmc_sweep ', nbc, nsw,  '  done'
          end if
#ENDIF
#IFNDEF GEN_CONFC_LEARNING
          !! first average over time
          nsiglR(:) = 0
          do n = 1, ltrot
              do i = 1, lq
                  nsiglR(i) = nsiglR(i) + nsigl_u(i,n)
              end do
          end do
          !! calculate spin-spin interaction
          do j = 1, lq
              do i = 1, lq
                  imj = latt_imj(i,j)
                  jjcorr_R(imj) = jjcorr_R(imj) + dble(nsiglR(i))*dble(nsiglR(j))/dble(ltrot*ltrot)
              end do
          end do
#ENDIF
      end do
#IFDEF GEN_CONFC_LEARNING
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
#ELSE
      call mpi_reduce( jjcorr_R, mpi_jjcorr_R, lq, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr )
      if( irank .eq. 0 ) then
          jjcorr_R(:) = mpi_jjcorr_R(:) / dble( isize*nsweep )
          open (unit=9095,file='jjcorrR.bin',status='unknown', action="write", position="append")
          do i = 1, lq
            write(9095, '(e16.8)') jjcorr_R(i)
          end do
      end if
#ENDIF
