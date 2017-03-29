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

             ! store the current state
             Ust_up_tmp(:,:,:) = Ust_up(:,:,:)
             Dst_up_tmp(:,:)   = Dst_up(:,:)
             Vst_up_tmp(:,:,:) = Vst_up(:,:,:)
             grup_tmp(:,:) = grup(:,:)

#IFDEF SPINDOWN
             Ust_dn_tmp(:,:,:) = Ust_dn(:,:,:)
             Dst_dn_tmp(:,:)   = Dst_dn(:,:)
             Vst_dn_tmp(:,:,:) = Vst_dn(:,:,:)
             grdn_tmp(:,:) = grdn(:,:)
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

             ! calculate the fermion part ratio
             ratiof = 1.d0

             do nt = 1, ltrot, 1
  
                 ! wrap H0
                 call mmthr  (grup, grdn)
                 call mmthlm1(grup, grdn)

                 ! update
                 ! updateu
                 if( lupdateu ) then
                     nflag = 3 ! onsite
                     call mmuur  ( grup, grdn, nf, nt, nflag )
                     call mmuulm1( grup, grdn, nf, nt, nflag )
                     do i = 1, lq
                         if( stcluster(i,nt) .eq. 1 ) then
                             !call stglobal_upgradej(nt,nf,i,grup,grdn,ratiofi)
                             call stglobal_upgradeu(nt,i,grup,grdn,ratiofi)
                             ratiof = ratiof * ratiofi
                         end if
                     end do
                 end if
  
                 if ( iwrap_nt(nt) .gt. 0 ) then
                     n = iwrap_nt(nt)
                     ! at tau = n * tau1
                     VL_up(:,:) = Vst_up(:,:,n)
                     DLvec_up(:)= Dst_up(:,n)
                     UL_up(:,:) = Ust_up(:,:,n)
#IFDEF SPINDOWN
                     VL_dn(:,:) = Vst_dn(:,:,n)
                     DLvec_dn(:)= Dst_dn(:,n)
                     UL_dn(:,:) = Ust_dn(:,:,n)
#ENDIF
                     call ftdqmc_stablize_0b_svd(n)

                     ! for spin up
                     UR_up(:,:)  = Ust_up(:,:,n)
                     DRvec_up(:) = Dst_up(:,n)
                     VR_up(:,:)  = Vst_up(:,:,n)
                     call green_equaltime( n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, grtmp, info )
                     call s_compare_max_z( ndim, grtmp, grup, max_wrap_error_tmp )
                     if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
                     if( info .eq. 0 ) grup(:,:) = grtmp(:,:)
#IFDEF TEST
                     write(fout,*)
                     write(fout, '(a,e16.8)') ' grup, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF

#IFDEF SPINDOWN
                     ! for spin down
                     UR_dn(:,:)  = Ust_dn(:,:,n)
                     DRvec_dn(:) = Dst_dn(:,n)
                     VR_dn(:,:)  = Vst_dn(:,:,n)
                     call green_equaltime( n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, grtmp, info )
                     call s_compare_max_z( ndim, grtmp, grdn, max_wrap_error_tmp )
                     if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
                     if( info .eq. 0 ) grdn(:,:) = grtmp(:,:)
#IFDEF TEST
                     write(fout,*)
                     write(fout, '(a,e16.8)') ' grdn, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
#ENDIF
                 end if
  
             end do

             if( ratiof .gt. spring_sfmt_stream() ) then
                 ! global update is accepted, perfrom an sweep from beta to 0
                 main_obs(3) = main_obs(3) + dcmplx(1.d0,1.d0)
#IFDEF TEST
                 write(fout,'(a,e16.8,a,i8)') ' global update accepted, ratiof = ', ratiof, '  nstcluster = ',  nstcluster
#ENDIF

                 Vst_up(:,:,nst) = Imat(:,:)
                 Dst_up(:,nst)   = Ivec(:)
                 Ust_up(:,:,nst) = Imat(:,:)
#IFDEF SPINDOWN
                 Vst_dn(:,:,nst) = Imat(:,:)
                 Dst_dn(:,nst)   = Ivec(:)
                 Ust_dn(:,:,nst) = Imat(:,:)
#ENDIF
  
                 nt_ob = ceiling( spring_sfmt_stream() * ltrot )

                 do nt = ltrot, 1, -1
                     if ( iwrap_nt(nt) .gt. 0 .and. nt.ne.ltrot ) then
                         n = iwrap_nt(nt)
                         ! at tau = n * tau1
                         UR_up(:,:) = Ust_up(:,:,n)
                         DRvec_up(:)= Dst_up(:,n)
                         VR_up(:,:) = Vst_up(:,:,n)

#IFDEF SPINDOWN
                         UR_dn(:,:) = Ust_dn(:,:,n)
                         DRvec_dn(:)= Dst_dn(:,n)
                         VR_dn(:,:) = Vst_dn(:,:,n)
#ENDIF

                         call ftdqmc_stablize_b0_svd(n+1)

                         ! for spin up
                         UL_up(:,:)  = Ust_up(:,:,n)
                         DLvec_up(:) = Dst_up(:,n)  
                         VL_up(:,:)  = Vst_up(:,:,n)
                         call green_equaltime( n, ndim, UR_up, DRvec_up, VR_up, VL_up, DLvec_up, UL_up, grtmp, info )
                         call s_compare_max_z( ndim, grtmp, grup, max_wrap_error_tmp )
                         if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp

#IFDEF TEST
                         write(fout,*)
                         write(fout, '(a,e16.8)') ' grup, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
                         ! whether use the scrath grup
                         if( info .eq. 0 ) grup(:,:) = grtmp(:,:)

#IFDEF SPINDOWN
                         ! for spin down
                         UL_dn(:,:)  = Ust_dn(:,:,n)
                         DLvec_dn(:) = Dst_dn(:,n)  
                         VL_dn(:,:)  = Vst_dn(:,:,n)
                         call green_equaltime( n, ndim, UR_dn, DRvec_dn, VR_dn, VL_dn, DLvec_dn, UL_dn, grtmp, info )
                         call s_compare_max_z( ndim, grtmp, grdn, max_wrap_error_tmp )
                         if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp

#IFDEF TEST
                         write(fout,*)
                         write(fout, '(a,e16.8)') ' grdn, max_wrap_error_tmp = ',  max_wrap_error_tmp
#ENDIF
                         ! whether use the scrath grdn
                         if( info .eq. 0 ) grdn(:,:) = grtmp(:,:)
#ENDIF
                     end if
  
                     !! update
                     ! updateu
                     if( lupdateu ) then
                         nflag = 3 ! onsite
                         call upgradeu( nt, grup, grdn )
                         call mmuul  ( grup, grdn, nf, nt, nflag )
                         call mmuurm1( grup, grdn, nf, nt, nflag )
                     end if
  
                     ! updatej
                     if( lupdatej ) then
                         do nf = nfam,1,-1
                             nflag = 2
                             call mmuul  ( grup, grdn, nf, nt, nflag )
                             call mmuurm1( grup, grdn, nf, nt, nflag )
                             call upgradej(nt,nf,grup,grdn)
                             nflag = 1
                             call mmuul  ( grup, grdn, nf, nt, nflag )
                             call mmuurm1( grup, grdn, nf, nt, nflag )
                         enddo
                     end if
  
                     ! wrap H0
                     call mmthl  (grup, grdn)
                     call mmthrm1(grup, grdn)
  
                 end do
  
                 ! at tau = 0
                 n = 0
                 Ust_up(:,:,n) = Imat(:,:)
                 Dst_up(:,n)   = Ivec(:)
                 Vst_up(:,:,n) = Imat(:,:)
#IFDEF SPINDOWN
                 Ust_dn(:,:,n) = Imat(:,:)
                 Dst_dn(:,n)   = Ivec(:)
                 Vst_dn(:,:,n) = Imat(:,:)
#ENDIF

             else
                 ! global update is rejected, recover the stored state
                 main_obs(3) = main_obs(3) + dcmplx(0.d0,1.d0)
#IFDEF TEST
                 write(fout,'(a,e16.8,a,i8)') ' global update rejected, ratiof = ', ratiof, '  nstcluster = ',  nstcluster
#ENDIF
                 Ust_up(:,:,:) =  Ust_up_tmp(:,:,:)
                 Dst_up(:,:)   =  Dst_up_tmp(:,:)
                 Vst_up(:,:,:) =  Vst_up_tmp(:,:,:)
                 grup(:,:)     =  grup_tmp(:,:)
#IFDEF SPINDOWN
                 Ust_dn(:,:,:) =  Ust_dn_tmp(:,:,:)
                 Dst_dn(:,:)   =  Dst_dn_tmp(:,:)
                 Vst_dn(:,:,:) =  Vst_dn_tmp(:,:,:)
                 grdn(:,:)     =  grdn_tmp(:,:)
#ENDIF
                 ! also you need flip back the spin
                 do nt = 1, ltrot
                 do i = 1, lq
                     nf = 1
                     if( stcluster(i,nt) .eq. 1 ) nsigl_u(i,nt) = nflipl( nsigl_u(i,nt), 1 )
                 end do
                 end do
             end if

         end if
      end if
