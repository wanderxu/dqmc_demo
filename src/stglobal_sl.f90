    subroutine ftdqmc_stglobal( lmeas )
      use mod_cumulate
      implicit none
      logical, intent(in) :: lmeas
      ! local variables
      integer :: nt, n, nf, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
      logical :: lterminate 
      real(dp) :: logweightf_old, logweightf_new, logweights_old, logweights_new, ratiof, logratiof
      complex(dp) :: logweightf_up, logweightf_dn
      integer :: ijs, i_1, i0, i1, i2, i3, i4, icum, inn, ntj
      integer, allocatable, dimension(:,:) :: nsigl_u_old
      real(dp), allocatable, dimension(:,:) :: heff_old
      real(dp) :: ediff, local_ratio, Heff_diff

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

#IFDEF TEST
             call Bmat_tau( ltrot, 1, Bdtau1_up, Bdtau1_dn )
             Btmp(:,:) = Bdtau1_up(:,:)
             do  i = 1, ndim
                 Btmp(i,i) = Bdtau1_up(i,i) + cone
             end do
             call s_logdet_z(ndim, Btmp, logweightf_up)

             Btmp(:,:) = Bdtau1_dn(:,:)
             do  i = 1, ndim
                 Btmp(i,i) = Bdtau1_dn(i,i) + cone
             end do
             call s_logdet_z(ndim, Btmp, logweightf_dn)
             write(fout,'(a,2e24.12)') ' without stablize, logweightf_up = ', logweightf_up
             write(fout,'(a,2e24.12)') ' without stablize, logweightf_dn = ', logweightf_dn
#ENDIF

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
#IFDEF TEST
             write(fout,'(a,2e24.12)') ' with stablize, logweightf_up = ', logweightf_up
             write(fout,'(a,2e24.12)') ' with stablize, logweightf_dn = ', logweightf_dn
#ENDIF

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

             !! build the space-time cluster to be performed global update on
             !! use the cummulate update algorithm
             
             ! store old field
             allocate( nsigl_u_old(lq,ltrot) )
             allocate( heff_old(lq,ltrot) )
             nsigl_u_old(:,:) = nsigl_u(:,:)
             heff_old(:,:) = heff(:,:)

             ! cumulate update
             Heff_diff = 0.d0
             do icum = 1, ncumulate
                 do nt = 1, ltrot
                     do i = 1, lq
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
                     end do
                 end do
             end do

             nstcluster = 0
             do nt = 1, ltrot
                 do i = 1, lq
                     stcluster(i,nt) = abs(nsigl_u(i,nt) - nsigl_u_old(i,nt))/2
                     nstcluster = nstcluster +  stcluster(i,nt)
                 end do
             end do

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

             !! calculate accept ratio
             ! log(accep_ratio) = (E'-E) - (Eff'-Eff)
             logratiof = logweightf_new + logweights_new - logweightf_old - logweights_old - Heff_diff
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
                 write(fout,'(a,e16.8)') ' logweightf_old = ', logweightf_old
                 write(fout,'(a,e16.8)') ' logweightf_new = ', logweightf_new
                 write(fout,'(a,e16.8)') ' Heff_diff = ', Heff_diff
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
                 write(fout,'(a,e16.8)') ' logweightf_old = ', logweightf_old
                 write(fout,'(a,e16.8)') ' logweightf_new = ', logweightf_new
                 write(fout,'(a,e16.8)') ' Heff_diff = ', Heff_diff
                 write(fout,'(a,e16.8)') ' weight_track = ', weight_track
#ENDIF
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
                 call ftdqmc_sweep_start_b0 ! recover old B matrix
                 call ftdqmc_sweep_0b(lupdate=.false., lmeasure=lmeas)
             end if

         end if
      end if

      deallocate( heff_old )
      deallocate( nsigl_u_old )
    end subroutine ftdqmc_stglobal
