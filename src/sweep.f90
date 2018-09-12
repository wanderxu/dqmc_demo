      call obser_init
      do nsw = 1, nsweep
          if(lstglobal .and. llocal ) then
              !! perform local and global update, only measure after global update
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure_equaltime=.false.)
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure_equaltime=.false.,lmeasure_dyn=.false.)
              call ftdqmc_stglobal(lmeas=.true.)
#ifdef GEN_CONFC_LEARNING
              ! output configuration for learning
              call outconfc_bin(weight_track)
              call preq
              call obser_init
#endif
          else if( lstglobal ) then
              call ftdqmc_stglobal(lmeas=.true.)
#ifdef GEN_CONFC_LEARNING
              ! output configuration for learning
              call outconfc_bin(weight_track)
              call preq
              call obser_init
#endif
          else if ( llocal ) then
              !! only perform local update, measure equaltime quantities during sweeps, measure dyn quantities when turnning off updates
              call ftdqmc_sweep_b0(lupdate=.true., lmeasure_equaltime=.true.)
#ifdef GEN_CONFC_LEARNING
              ! output configuration for learning
              call outconfc_bin(weight_track)
              call preq
              call obser_init
#endif
              if(ltau) then
                  call push_stage
                  call ftdqmc_sweep_0b(lupdate=.false., lmeasure_equaltime=.false., lmeasure_dyn=ltau)
                  call pop_stage
              end if
              
              call ftdqmc_sweep_0b(lupdate=.true., lmeasure_equaltime=.true., lmeasure_dyn=.false.)
#ifdef GEN_CONFC_LEARNING
              ! output configuration for learning
              call outconfc_bin(weight_track)
              call preq
              call obser_init
#endif
          else
              stop ' lstglobal and llocal should  not both false '
          end if
#ifdef TEST
          if( irank .eq. 0 ) then
              write(fout,'(a,i4,i4,a)') ' ftdqmc_sweep ', nbc, nsw,  '  done'
          end if
#endif
      end do
#ifndef GEN_CONFC_LEARNING
      call preq  ! reduce
#endif
      if(ltau) call prtau
