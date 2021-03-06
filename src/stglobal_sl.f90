    subroutine ftdqmc_stglobal( lmeas )
      use mod_cumulate
      implicit none
      logical, intent(in) :: lmeas
      ! local variables
      integer :: nt, n, nf, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
      logical :: lterminate 
      real(dp) :: ratiof, logratiof
      integer ::  icum, inn, ntj
      integer, allocatable, dimension(:,:) :: nsigl_u_old
      real(dp), allocatable, dimension(:,:) :: heff_old
      real(dp) :: ediff, local_ratio, Heff_diff
      integer :: ltrot_lq, iltlq

      ! perform global update
      if( lstglobal ) then
         icount_nsw_stglobal = icount_nsw_stglobal  + 1
         if( icount_nsw_stglobal .eq. nsw_stglobal ) then
             icount_nsw_stglobal = 0

#ifdef TEST
             write(fout,*)
             write(fout, '(a)') ' >>>>>>>>>> '
             write(fout, '(a)') ' in space time global update '
             write(fout, '(a)') ' >>>>>>>>>> '
             write(fout,*)
#endif

             !!============================================================================================
             !! build the space-time cluster to be performed global update on
             !! use the cummulate update algorithm
             
             ! store old field
             allocate( nsigl_u_old(lq,ltrot) )
             allocate( heff_old(lq,ltrot) )
             nsigl_u_old(:,:) = nsigl_u(:,:)
             heff_old(:,:) = heff(:,:)
             ! also store UDV matrix and Green function
             call push_stage

             ! cumulate update
             ltrot_lq = ltrot*lq
             Heff_diff = 0.d0
             do icum = 1, ncumulate*ltrot_lq
                 !!do nt = 1, ltrot
                 !!    do i = 1, lq
                        iltlq = ceiling( spring_sfmt_stream() * ltrot*lq )
                        nt =  (iltlq-1)/lq + 1
                        i = mod(iltlq-1, lq) + 1
                        ediff=-2.d0*nsigl_u(i,nt)*heff(i,nt) 
                        if( ediff .gt. 0 ) then
                            local_ratio = 1.001d0
                        else
                            local_ratio = exp(ediff)
                        end if
                        if( local_ratio .gt. spring_sfmt_stream() ) then
                            ! update heff
                            do inn = 1, num_nei
                                j = nei_cord(1,inn,i,nt)
                                ntj = nei_cord(2,inn,i,nt)
                                heff(j,ntj) = heff(j,ntj) - 2.d0*nei_Jeff(inn,i,nt)*nsigl_u(i,nt)
                            end do
                            ! flip field
                            nsigl_u(i,nt) = -nsigl_u(i,nt)
                            ! add ediff to Heff_diff
                            Heff_diff = Heff_diff + ediff
                        end if
                 !!    end do
                 !!end do
             end do

             nstcluster = 0
             do nt = 1, ltrot
                 do i = 1, lq
                     stcluster(i,nt) = abs(nsigl_u(i,nt) - nsigl_u_old(i,nt))/2
                     nstcluster = nstcluster +  stcluster(i,nt)
                 end do
             end do

             logweightf_old = dble( logweightf_up + logweightf_dn )*2.d0
             call ftdqmc_sweep_start_b0   ! update B(beta,0)
             logweightf_new = dble( logweightf_up + logweightf_dn )*2.d0
             call ftdqmc_calculate_weights( logweights_new )

             !!============================================================================================
             !! calculate accept ratio
             ! log(accep_ratio) = (E'-E) - (Eff'-Eff)
             logratiof = logweightf_new + logweights_new - logweightf_old - logweights_old - Heff_diff
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
#ifdef TEST
                 write(fout,'(a,e24.16,a,i8)') ' global update accepted, logratiof = ', logratiof, '  nstcluster = ',  nstcluster
                 write(fout,'(a,e24.16)') ' logweights_old = ', logweights_old
                 write(fout,'(a,e24.16)') ' logweights_new = ', logweights_new
                 write(fout,'(a,e24.16)') ' logweightf_old = ', logweightf_old
                 write(fout,'(a,e24.16)') ' logweightf_new = ', logweightf_new
                 write(fout,'(a,e24.16)') ' weight_track = ', weight_track
                 write(fout,'(a,e24.16)') ' Heff_diff = ', Heff_diff
#endif
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
#ifdef TEST
                 write(fout,'(a,e24.16,a,i8)') ' global update rejected, logratiof = ', logratiof, '  nstcluster = ',  nstcluster
                 write(fout,'(a,e24.16)') ' logweights_old = ', logweights_old
                 write(fout,'(a,e24.16)') ' logweights_new = ', logweights_new
                 write(fout,'(a,e24.16)') ' logweightf_old = ', logweightf_old
                 write(fout,'(a,e24.16)') ' logweightf_new = ', logweightf_new
                 write(fout,'(a,e24.16)') ' weight_track = ', weight_track
                 write(fout,'(a,e24.16)') ' Heff_diff = ', Heff_diff
#endif
                 ! global update is rejected, you need flip back the spin
                 nsigl_u(:,:) = nsigl_u_old(:,:)
                 heff(:,:) = heff_old(:,:)
                 !!!do nt = 1, ltrot
                 !!!do i = 1, lq
                 !!!    nf = 1
                 !!!    if( stcluster(i,nt) .eq. 1 ) nsigl_u(i,nt) = nflipl( nsigl_u(i,nt), 1 )
                 !!!end do
                 !!!end do

                 ! perform measurement  ! note no matter whether the update is aceepted, you should do measrement
                 if( lmeas ) then
                     call ftdqmc_sweep_start_b0 ! first, you should recover old B matrix at tau=0
                     call ftdqmc_sweep_0b(lupdate=.false., lmeasure_equaltime=lmeas, lmeasure_dyn=ltau)
                 else
                     ! if no meas, just recover old UDV matrix and Green functions
                     if(nst.gt.0 .or. llocal) then
                         call pop_stage
                     end if
                 end if
             end if

         end if
      end if

      if(allocated(heff_old)) deallocate( heff_old )
      if(allocated(nsigl_u_old)) deallocate( nsigl_u_old )
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
