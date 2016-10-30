    subroutine ftdqmc_stglobal( lmeas )
      implicit none
      logical, intent(in) :: lmeas
      ! local variables
      integer :: nt, n, nf, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
      logical :: lterminate 
      real(dp) :: logweightf_old, logweightf_new, logweights_old, logweights_new, ratiof, logratiof
      complex(dp) :: logweightf_up, logweightf_dn
      integer :: ijs, i_1, i0, i1, i2, i3, i4

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
             ! WARNNING, make sure you are at tau=beta

             ! calculate det(1+B(beta,0))
             ! det( I + UDV ) = det( I + DVU )
             ! at tau = beta
             UR_up(:,:) = Ust_up(:,:,nst)
             DRvec_up(:)= Dst_up(:,nst)
             VR_up(:,:) = Vst_up(:,:,nst)
             call zgemm('n','n',ndim,ndim,ndim,cone,VR_up,ndim,UR_up,ndim,czero,Atmp,ndim)  ! Atmp = V*U
             call s_diag_d_x_z(ndim,DRvec_up,Atmp,Btmp) ! Btmp = D * Atmp = DVU
             do  i = 1, ndim
                 Btmp(i,i) = Btmp(i,i) + cone
             end do
             call s_logdet_z(ndim, Btmp, logweightf_up)
#IFDEF SPINDOWN
             ! at tau = beta
             UR_dn(:,:) = Ust_dn(:,:,nst)
             DRvec_dn(:)= Dst_dn(:,nst)
             VR_dn(:,:) = Vst_dn(:,:,nst)
             call zgemm('n','n',ndim,ndim,ndim,cone,VR_dn,ndim,UR_dn,ndim,czero,Atmp,ndim)  ! Atmp = V*U
             call s_diag_d_x_z(ndim,DRvec_dn,Atmp,Btmp) ! Btmp = D * Atmp = DVU
             do  i = 1, ndim
                 Btmp(i,i) = Btmp(i,i) + cone
             end do
             call s_logdet_z(ndim, Btmp, logweightf_dn)
#ENDIF
             logweightf_old = dble( logweightf_up + logweightf_dn )*2.d0

#IFDEF GEN_CONFC_LEARNING
             ! calculate boson part ratio
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
             logweights_old = dtau*js*dble(ijs)

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
             logweights_old = logweights_old + gamma_s*dble(ijs)
#ENDIF

             !! build the space-time cluster to be performed global update on
             !! use the Wolff algorithm

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
                                  ( nsigl_u( nn_ilq, nn_it) .eq. int(sign(js,1.d0))*nsigl_u(ilq, it) ) .and. &
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

#IFDEF GEN_CONFC_LEARNING
             ! calculate boson part ratio
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
#IFDEF TEST
             write(fout,'(a,i20)') 'new ijs for spatial = ', ijs
#ENDIF 
             logweights_new = dtau*js*dble(ijs)

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
#IFDEF TEST
             write(fout,'(a,i20)') 'new ijs for temporal = ', ijs
#ENDIF 
             logweights_new = logweights_new + gamma_s*dble(ijs)
#ENDIF

             ! calculate the fermion part ratio
             ! ratio = weight_new / weight_old

             call ftdqmc_sweep_start_b0   ! update B(beta,0)
             ! calculate  new_det(1+B(beta,0))
             ! at tau = 0
             UR_up(:,:) = Ust_up(:,:,0)
             DRvec_up(:)= Dst_up(:,0)
             VR_up(:,:) = Vst_up(:,:,0)
             call zgemm('n','n',ndim,ndim,ndim,cone,VR_up,ndim,UR_up,ndim,czero,Atmp,ndim)  ! Atmp = V*U
             call s_diag_d_x_z(ndim,DRvec_up,Atmp,Btmp) ! Btmp = D * Atmp = DVU
             do  i = 1, ndim
                 Btmp(i,i) = Btmp(i,i) + cone
             end do
             call s_logdet_z(ndim, Btmp, logweightf_up)
#IFDEF SPINDOWN
             ! at tau = 0
             UR_dn(:,:) = Ust_dn(:,:,0)
             DRvec_dn(:)= Dst_dn(:,0)
             VR_dn(:,:) = Vst_dn(:,:,0)
             call zgemm('n','n',ndim,ndim,ndim,cone,VR_dn,ndim,UR_dn,ndim,czero,Atmp,ndim)  ! Atmp = V*U
             call s_diag_d_x_z(ndim,DRvec_dn,Atmp,Btmp) ! Btmp = D * Atmp = DVU
             do  i = 1, ndim
                 Btmp(i,i) = Btmp(i,i) + cone
             end do
             call s_logdet_z(ndim, Btmp, logweightf_dn)
#ENDIF
             logweightf_new = dble( logweightf_up + logweightf_dn )*2.d0

             logratiof = logweightf_new - logweightf_old
             if( logratiof .gt. 0 ) then
                 ratiof = 1.001d0
             else
                 ratiof = exp(logratiof)
             end if

             if( ratiof .gt. spring_sfmt_stream() ) then
                 ! global update is accepted
#IFDEF GEN_CONFC_LEARNING
                 weight_track = logweightf_new + logweights_new
#ENDIF
                 main_obs(3) = main_obs(3) + dcmplx(1.d0,1.d0)
                 main_obs(4) = main_obs(4) + dcmplx(dble(nstcluster),dble(ltrot*lq))
#IFDEF TEST
                 write(fout,'(a,e16.8,a,i8)') ' global update accepted, logratiof = ', logratiof, '  nstcluster = ',  nstcluster
                 write(fout,'(a,e16.8)') ' logweights_old = ', logweights_old
                 write(fout,'(a,e16.8)') ' logweights_new = ', logweights_new
                 write(fout,'(a,e16.8)') ' weight_track = ', weight_track
#ENDIF
                 ! perform measurement
                 call ftdqmc_sweep_0b(lupdate=.false., lmeasure=lmeas)
             else
                 ! global update is rejected
#IFDEF GEN_CONFC_LEARNING
                 weight_track = logweightf_old + logweights_old
#ENDIF
                 main_obs(3) = main_obs(3) + dcmplx(0.d0,1.d0)
#IFDEF TEST
                 write(fout,'(a,e16.8,a,i8)') ' global update rejected, logratiof = ', logratiof, '  nstcluster = ',  nstcluster
                 write(fout,'(a,e16.8)') ' logweights_old = ', logweights_old
                 write(fout,'(a,e16.8)') ' logweights_new = ', logweights_new
                 write(fout,'(a,e16.8)') ' weight_track = ', weight_track
#ENDIF
                 ! global update is rejected, you need flip back the spin
                 do nt = 1, ltrot
                 do i = 1, lq
                     nf = 1
                     if( stcluster(i,nt) .eq. 1 ) nsigl_u(i,nt) = nflipl( nsigl_u(i,nt), 1 )
                 end do
                 end do

                 ! perform measurement  ! note no matter whether the update is aceepted, you should do measrement
                 call ftdqmc_sweep_start_b0 ! recover old B matrix
                 call ftdqmc_sweep_0b(lupdate=.false., lmeasure=lmeas)
             end if

         end if
      end if
    end subroutine ftdqmc_stglobal
