      call obser_init
      do nsw = 1, nsweep
          if(lstglobal .and. llocal ) then
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure=.false.)
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure=.false.)
              call ftdqmc_stglobal(lmeas=.true.)
#IFDEF GEN_CONFC_LEARNING
              ! output configuration for learning
              call outconfc_bin(weight_track)
              call preq
              call obser_init
#ENDIF
          else if( lstglobal ) then
              call ftdqmc_stglobal(lmeas=.true.)
#IFDEF GEN_CONFC_LEARNING
              ! output configuration for learning
              call outconfc_bin(weight_track)
              call preq
              call obser_init
#ENDIF
          else if ( llocal ) then
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
          else
              stop ' lstglobal and llocal should  not both false '
          end if
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
