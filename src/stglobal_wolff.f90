    subroutine ftdqmc_stglobal( lmeas )
      implicit none
      logical, intent(in) :: lmeas
      ! local variables
      integer :: nt, n, nf, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
      logical :: lterminate 
      real(dp) :: ratiof, logratiof
      complex(dp) :: logwf_up_tmp, logwf_dn_tmp

      ! perform global update
      if( lstglobal ) then
         icount_nsw_stglobal = icount_nsw_stglobal  + 1
         if( icount_nsw_stglobal .eq. nsw_stglobal ) then
             icount_nsw_stglobal = 0

#IFDEF TEST
             write(fout,*)
             write(fout, '(a)') ' >>>>>>>>>> '
             write(fout, '(a)') ' in space time global update '
             write(fout, '(a)') ' >>>>>>>>>> '
             write(fout,*)
#ENDIF

             !! build the space-time cluster to be performed global update on
             !! use the Wolff algorithm

             ! first store UDV matrix and Green function
             call push_stage

             ! start build cluster
             tentacle_old(:,:) = 0
             tentacle(:,:) = 0
             stcluster(:,:) = 0

             ! choose an inital site randomly
             ntentacle = 1
             ilq = ceiling( spring_sfmt_stream()*dble(lq)    )
             it  = ceiling( spring_sfmt_stream()*dble(ltrot) )
             tentacle(1,1) = ilq
             tentacle(2,1) = it
             stcluster( ilq, it ) = 1
             nstcluster = 1
             lterminate = .false.

             do while ( .not. lterminate )
                 ntentacle_old = ntentacle
                 tentacle_old(:,:) = tentacle(:,:)
                 ntentacle = 0
                 tentacle(:,:) = 0
                 do i = 1, ntentacle_old
                     ilq = tentacle_old(1,i)
                     it  = tentacle_old(2,i)
                     do inn_st = 1, num_st_nn
                         nn_ilq =  stbonds_neib(1, inn_st, ilq, it )
                         nn_it  =  stbonds_neib(2, inn_st, ilq, it )
                         ! here ratio_nn_st(1:4) = 1 - exp(-2*dtau*|js|);     ratio_nn_st(5:6) = 1 - tanh(dtau*h)
                         if (inn_st .lt. 5) then 
                         ! antiferromagnetic interaction in space direction if js < 0
                         ! ferromagnetic interaction in space direction if js > 0
                             if ( ( stcluster( nn_ilq, nn_it ) .eq. 0 ) .and. &
                                  ( nsigl_u( nn_ilq, nn_it) .eq. int(sign(1.d0,js))*nsigl_u(ilq, it) ) .and. &
                                  ( ratio_nn_st(inn_st) .gt. spring_sfmt_stream() ) ) then
                                 stcluster( nn_ilq, nn_it ) = 1
                                 nstcluster = nstcluster + 1
                                 ntentacle = ntentacle + 1
                                 tentacle(1,ntentacle) = nn_ilq
                                 tentacle(2,ntentacle) = nn_it
                             end if
                         else
                         ! ferromagnetic interaction in imaginary time direction
                             if ( ( stcluster( nn_ilq, nn_it ) .eq. 0 ) .and. &
                                  ( nsigl_u( nn_ilq, nn_it) .eq. nsigl_u(ilq, it) ) .and. &
                                  ( ratio_nn_st(inn_st) .gt. spring_sfmt_stream() ) ) then
                                 stcluster( nn_ilq, nn_it ) = 1
                                 nstcluster = nstcluster + 1
                                 ntentacle = ntentacle + 1
                                 tentacle(1,ntentacle) = nn_ilq
                                 tentacle(2,ntentacle) = nn_it
                             end if
                         end if
                     end do
                 end do
                 if( ntentacle.eq.0 ) lterminate = .true.
             end do

             ! flip first
             do nt = 1, ltrot
             do i = 1, lq
                 nf = 1
                 if( stcluster(i,nt) .eq. 1 ) nsigl_u(i,nt) = nflipl( nsigl_u(i,nt), 1 )
             end do
             end do

             logweightf_old = dble( logweightf_up + logweightf_dn )*2.d0
             ! Store weight
             logwf_up_tmp = logweightf_up
             logwf_dn_tmp = logweightf_dn

             call ftdqmc_sweep_start_b0   ! update B(beta,0)
             logweightf_new = dble( logweightf_up + logweightf_dn )*2.d0
             call ftdqmc_calculate_weights( logweights_new )

             !!============================================================================================
             !! calculate accept ratio
             logratiof = logweightf_new - logweightf_old
             if( logratiof .gt. 0 ) then
                 ratiof = 1.001d0
             else
                 ratiof = exp(logratiof)
             end if

             !!============================================================================================
             ! update
             if( ratiof .gt. spring_sfmt_stream() ) then
                 ! global update is accepted
                 weight_track = logweightf_new + logweights_new
                 main_obs(3) = main_obs(3) + dcmplx(1.d0,1.d0)
                 main_obs(4) = main_obs(4) + dcmplx(dble(nstcluster),dble(ltrot*lq))
#IFDEF TEST
                 write(fout,'(a,e24.16,a,i8)') ' global update accepted, logratiof = ', logratiof, '  nstcluster = ',  nstcluster
                 write(fout,'(a,e24.16)') ' logweights_old = ', logweights_old
                 write(fout,'(a,e24.16)') ' logweights_new = ', logweights_new
                 write(fout,'(a,e24.16)') ' logweightf_old = ', logweightf_old
                 write(fout,'(a,e24.16)') ' logweightf_new = ', logweightf_new
                 write(fout,'(a,e24.16)') ' weight_track = ', weight_track
#ENDIF
                 ! update logweight
                 logweightf_old = logweightf_new
                 logweights_old = logweights_new
                 ! perform measurement
                 if( lmeas .or. llocal ) then ! when we have local update, we also need perform an sweep back to beta
                     call ftdqmc_sweep_0b(lupdate=.false., lmeasure_equaltime=lmeas, lmeasure_dyn=ltau)
                 end if
             else
                 ! global update is rejected
                 weight_track = logweightf_old + logweights_old
                 main_obs(3) = main_obs(3) + dcmplx(0.d0,1.d0)
#IFDEF TEST
                 write(fout,'(a,e24.16,a,i8)') ' global update rejected, logratiof = ', logratiof, '  nstcluster = ',  nstcluster
                 write(fout,'(a,e24.16)') ' logweights_old = ', logweights_old
                 write(fout,'(a,e24.16)') ' logweights_new = ', logweights_new
                 write(fout,'(a,e24.16)') ' logweightf_old = ', logweightf_old
                 write(fout,'(a,e24.16)') ' logweightf_new = ', logweightf_new
                 write(fout,'(a,e24.16)') ' weight_track = ', weight_track
#ENDIF
                 ! global update is rejected, you need flip back the spin
                 do nt = 1, ltrot
                 do i = 1, lq
                     nf = 1
                     if( stcluster(i,nt) .eq. 1 ) nsigl_u(i,nt) = nflipl( nsigl_u(i,nt), 1 )
                 end do
                 end do

                 ! Reset weight
                 logweightf_up = logwf_up_tmp
                 logweightf_dn = logwf_dn_tmp

                 ! perform measurement  ! note no matter whether the update is aceepted, you should do measrement
                 if( lmeas ) then
                     call ftdqmc_sweep_start_b0 ! first, you should recover old B matrix at tau=0
                     call ftdqmc_sweep_0b(lupdate=.false., lmeasure_equaltime=lmeas, lmeasure_dyn=ltau)
                 else
                     ! if no meas, just recover old UDV matrix and Green functions
                     call pop_stage
                 end if
             end if

         end if
      end if
    end subroutine ftdqmc_stglobal

    subroutine ftdqmc_calculate_weights( logweights )
      implicit none
      real(dp), intent(out) :: logweights
      integer :: nt, nf, i, ijs, i_1, i0, i1, i2, i3, i4
      logweights = 0.d0
      !!============================================================================================
      !!! calculate boson part ratio
      ijs = 0
!$OMP PARALLEL &
!$OMP PRIVATE ( nt, nf, i_1, i0, i1, i2, i3, i4 )
!$OMP DO REDUCTION ( + : ijs )
      do nt = 1, ltrot
          do nf = 1,2
              do i_1 = 1,lq/4
                  i0 = lthf(i_1,nf)
                  i1 = nnlist(i0,1)
                  i2 = nnlist(i0,2)
                  i3 = nnlist(i0,3)
                  i4 = nnlist(i0,4)
                  ijs = ijs + nsigl_u(i0,nt)*nsigl_u(i1,nt) + &
                              nsigl_u(i0,nt)*nsigl_u(i2,nt) + &
                              nsigl_u(i0,nt)*nsigl_u(i3,nt) + &
                              nsigl_u(i0,nt)*nsigl_u(i4,nt)
              end do
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      logweights = dtau*js*dble(ijs)

      ijs = 0
!$OMP PARALLEL &
!$OMP PRIVATE ( nt, i )
!$OMP DO REDUCTION ( + : ijs )
      do nt = 1, ltrot
          do i = 1, lq
              ijs = ijs + nsigl_u(i,nt)*nsigl_u(i,mod(nt,ltrot)+1)
          end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      logweights = logweights + gamma_s*dble(ijs)
    end subroutine ftdqmc_calculate_weights
