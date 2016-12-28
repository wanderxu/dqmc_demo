#IFDEF GEN_CONFC_LEARNING
      totsz_bin(:) = 0
#ELSE
      jjcorr_Rtau(:,:) = 0
#ENDIF
      do nsw = 1, nsweep
          if(lstglobal .and. llocal ) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure_equaltime=.false.)
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure_equaltime=.false.,lmeasure_dyn=.false.)
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
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure_equaltime=.false.)
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
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure_equaltime=.false.,lmeasure_dyn=.false.)
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
          !! calculate spin-spin interaction
          do ntj = 1, ltrot
            do nti = 1, ltrot
              n = mod(nti-ntj + ltrot, ltrot) + 1
              if( n .le. (ltrot/2+1) ) then
                do j = 1, lq
                  do i = 1, lq
                    imj = latt_imj(i,j)
                    jjcorr_Rtau(imj,n) = jjcorr_Rtau(imj,n) + nsigl_u(i,nti)*nsigl_u(j,ntj)
                  end do
                end do
              end if
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
      call mpi_reduce( jjcorr_Rtau, mpi_jjcorr_Rtau, lq*(ltrot/2+1), mpi_integer, mpi_sum, 0, mpi_comm_world, ierr )
      if( irank .eq. 0 ) then
          jjcorr_Rtau_real(:,:) = dble( mpi_jjcorr_Rtau(:,:) ) / dble( isize*nsweep )
          open (unit=9095,file='jjcorrRtau.bin',status='unknown', action="write", position="append")
          do n = 1, ltrot/2+1
              do i = 1, lq
                write(9095, '(e16.8)') jjcorr_Rtau_real(i,n)
              end do
          end do
      end if
#ENDIF
