             ! order_branch
             order_branch = 0
             do i = 1, ndim
                 order_branch = order_branch + nsigl_u(i,nt)
             end do

             if(order_branch .lt. 0 ) then

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
                         !call stglobal_upgradej(nt,nf,i,grup,grdn,ratiofi)
                         call stglobal_upgradeu(nt,i,grup,grdn,ratiofi)
                         !ratiof = ratiof * ratiofi
                     end do
                 end if
  
                 if ( mod(nt, nwrap) .eq. 0 ) then
                     n = nt/nwrap
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

                 Vst_up(:,:,nst) = Imat(:,:)
                 Dst_up(:,nst)   = Ivec(:)
                 Ust_up(:,:,nst) = Imat(:,:)
#IFDEF SPINDOWN
                 Vst_dn(:,:,nst) = Imat(:,:)
                 Dst_dn(:,nst)   = Ivec(:)
                 Ust_dn(:,:,nst) = Imat(:,:)
#ENDIF
  
                 !nt_ob = ceiling( spring_sfmt_stream() * ltrot )

                 do nt = ltrot, 1, -1
                     if ( mod(nt, nwrap) .eq. 0 .and. (nt/nwrap) .lt. nst ) then
                         n = nt/nwrap
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
                         !call upgradeu( nt, grup, grdn )  ! no update
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
         end if  ! if(order_branch .lt. 0 ) then
