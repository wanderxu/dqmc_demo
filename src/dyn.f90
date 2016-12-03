          ! dyn
          ! g00up g00dn
          ! gt0up gt0dn
          if( ( ltau .and. ltauall ) .or. ( ltau .and. nt .eq. ltrot/2 ) ) then
              IF ( lmeasure ) THEN
              !IF ( nt .le. ltrot/2 .or.  nt .eq. ltrot ) THEN
              !!!!!!!!n = nt/nwrap
!!!              if( nt .eq. 1 ) then
!!!                  Bt2t1_up = Imat
!!!#IFDEF SPINDOWN
!!!                  Bt2t1_dn = Imat
!!!#ENDIF
!!!              end if
              if( iwrap_nt(nt) .gt. 0 ) then

!!!                  Bt2t1_up = Imat
!!!                  UR_up(:,:) = Ust_up(:,:,n)
!!!                  DRvec_up(:) = Dst_up(:,n)
!!!                  VR_up(:,:) = Vst_up(:,:,n)
!!!
!!!                  ! gt0up = gttup * UDV
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,grup,ndim,UR_up,ndim,czero,Vtmp,ndim)  ! Vtmp = gttup * U
!!!                  call s_v_d_u( ndim, Vtmp, DRvec_up, VR_up, gt0up )
!!!
!!!                  ! g0tup = - V^-1 * D^-1 * U^-1 * (1-gttup)
!!!                  grtmp(:,:) = grup(:,:) - Imat(:,:)  ! - ( 1 - G(t,t) ), note here already include the sign
!!!                  call s_inv_z(ndim,VR_up)
!!!                  call s_inv_z(ndim,UR_up)
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,UR_up,ndim,grtmp,ndim,czero,Vtmp,ndim)  ! Vtmp = U^-1 * (1-gttup)
!!!                  call s_v_invd_u( ndim, VR_up, DRvec_up, Vtmp, g0tup )
!!!#IFDEF SPINDOWN
!!!                  Bt2t1_dn = Imat
!!!                  UR_dn(:,:) = Ust_dn(:,:,n)
!!!                  DRvec_dn(:) = Dst_dn(:,n)
!!!                  VR_dn(:,:) = Vst_dn(:,:,n)
!!!
!!!                  ! gt0dn = gttdn * UDV
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,grdn,ndim,UR_dn,ndim,czero,Vtmp,ndim)  ! Vtmp = gttdn * U
!!!                  call s_v_d_u( ndim, Vtmp, DRvec_dn, VR_dn, gt0dn )
!!!#ENDIF
              else
                  ! B(nt1,nt2) with nt1 >= nt2
                  nt1 = nt
                  nt2 = nt
                  ! G(t',0) = B(t',t) * G(t,0)
                  call Bmat_tau_R( nt1, nt2, gt0up, gt0dn)

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call Bmatinv_tau_L( nt1, nt2, g0tup, g0tdn)

!!!                  ! Bt2t1_up = Bdtau1_up * Bt2t1_up
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Bdtau1_up,ndim,Bt2t1_up,ndim,czero,Btmp,ndim)  ! Btmp =  Bdtau1_up * Bt2t1_up
!!!                  Bt2t1_up = Btmp
!!!
!!!                  UR_up(:,:) = Ust_up(:,:,n)
!!!                  DRvec_up(:) = Dst_up(:,n)
!!!                  VR_up(:,:) = Vst_up(:,:,n)
!!!
!!!                  ! gt0up = gttup * Bt2t1_up * UDV
!!!
!!!                  ! -> Bt2t1_up * U
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Bt2t1_up,ndim,UR_up,ndim,czero,Atmp,ndim)  ! Atmp =  Bt2t1_up * UR_up
!!!
!!!                  ! -> Bt2t1_up * UD = U1 * D1 * V1
!!!                  call s_z_x_diag_d(ndim,Atmp,DRvec_up,Btmp) ! Btmp = Atmp * D = Bt2t1_up * UD
!!!                  call s_svd_zg(ndim, ndim, ndim, Btmp, Umat1, Dvec1, Vmat1)
!!!
!!!                  ! -> then gt0up = gttup * U1 * D1 * V1 * V
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,grup,ndim,Umat1,ndim,czero,Vtmp,ndim)  ! Vtmp = gttup * Umat1
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Vmat1,ndim,VR_up,ndim,czero,Btmp,ndim)  ! Btmp = V1 * V
!!!                  call s_v_d_u( ndim, Vtmp, Dvec1, Btmp, gt0up )
!!!#IFDEF DYNERROR
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,grup,ndim,Atmp,ndim,czero,Vtmp,ndim)  ! Vtmp = gttup * Atmp
!!!                  call s_v_d_u( ndim, Vtmp, DRvec_up, VR_up, Btmp )  ! Btmp = gt0up_tmp
!!!                  call s_compare_max_z( ndim, gt0up, Btmp, xmax_dyn_tmp )
!!!                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
!!!#IFDEF TEST
!!!                  write(fout,*)
!!!                  write(fout, '(a,e16.8)') ' gt0up, gt0up_tmp, xmax_dyn_tmp = ',  xmax_dyn_tmp
!!!#ENDIF
!!!#ENDIF
!!!                  ! g0tup = - ( Bt2t1_up * UDV )^-1 * ( 1 - gttup )
!!!                  grtmp(:,:) = grup(:,:) - Imat(:,:)  ! - ( 1 - G(t,t) ), note here already include the sign
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Vmat1,ndim,VR_up,ndim,czero,Vtmp,ndim)  ! Vtmp =  V1 * V
!!!                  call s_inv_z(ndim,Vtmp)
!!!                  call s_inv_z(ndim,Umat1)
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Umat1,ndim,grtmp,ndim,czero,Btmp,ndim)  ! Btmp = U1^-1 * (1-gttup)
!!!                  call s_v_invd_u( ndim, Vtmp, Dvec1, Btmp, g0tup )
!!!
!!!#IFDEF SPINDOWN
!!!                  ! Bt2t1_dn = Bdtau1_dn * Bt2t1_dn
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Bdtau1_dn,ndim,Bt2t1_dn,ndim,czero,Btmp,ndim)  ! Btmp =  Bdtau1_dn * Bt2t1_dn
!!!                  Bt2t1_dn = Btmp
!!!
!!!                  UR_dn(:,:) = Ust_dn(:,:,n)
!!!                  DRvec_dn(:) = Dst_dn(:,n)
!!!                  VR_dn(:,:) = Vst_dn(:,:,n)
!!!
!!!                  ! gt0dn = gttdn * Bt2t1_dn * UDV
!!!
!!!                  ! -> Bt2t1_dn * U
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Bt2t1_dn,ndim,UR_dn,ndim,czero,Atmp,ndim)  ! Atmp =  Bt2t1_dn * UR_dn
!!!
!!!                  ! -> Bt2t1_dn * UD = U1 * D1 * V1
!!!                  call s_z_x_diag_d(ndim,Atmp,DRvec_dn,Btmp) ! Btmp = Atmp * D = Bt2t1_dn * UD
!!!                  call s_svd_zg(ndim, ndim, ndim, Btmp, Umat1, Dvec1, Vmat1)
!!!
!!!                  ! -> then gt0dn = gttdn * U1 * D1 * V1 * V
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,grdn,ndim,Umat1,ndim,czero,Vtmp,ndim)  ! Vtmp = gttdn * Umat1
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,Vmat1,ndim,VR_dn,ndim,czero,Btmp,ndim)  ! Btmp = V1 * V
!!!                  call s_v_d_u( ndim, Vtmp, Dvec1, Btmp, gt0dn )
!!!#IFDEF DYNERROR
!!!                  call zgemm('n','n',ndim,ndim,ndim,cone,grdn,ndim,Atmp,ndim,czero,Vtmp,ndim)  ! Vtmp = gttdn * Atmp
!!!                  call s_v_d_u( ndim, Vtmp, DRvec_dn, VR_dn, Btmp )  ! Btmp = gt0dn_tmp
!!!                  call s_compare_max_z( ndim, gt0dn, Btmp, xmax_dyn_tmp )
!!!                  if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
!!!#IFDEF TEST
!!!                  write(fout,*)
!!!                  write(fout, '(a,e16.8)') ' gt0dn, gt0dn_tmp, xmax_dyn_tmp = ',  xmax_dyn_tmp
!!!#ENDIF
!!!#ENDIF
!!!
!!!#ENDIF
              end if

#IFDEF TEST_LEVEL3
              write(fout,*)
              write(fout, '(a)') ' gt0up(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') gt0up(i,:)
              end do

              write(fout,*)
              write(fout, '(a)') ' g0tup(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') g0tup(i,:)
              end do
#IFDEF SPINDOWN
              write(fout,*)
              write(fout, '(a)') ' gt0dn(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') gt0dn(i,:)
              end do

              write(fout,*)
              write(fout, '(a)') ' g0tdn(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') g0tdn(i,:)
              end do
#ENDIF
              write(fout,*)
              write(fout, '(a)') ' gttup(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grup(i,:)
              end do
#IFDEF SPINDOWN
              write(fout,*)
              write(fout, '(a)') ' gttdn(:,:) = '
              do i = 1, ndim
                  write(fout,'(4(2e12.4))') grdn(i,:)
              end do
#ENDIF
#ENDIF
              call obsert(nt,gt0up,gt0dn,g0tup,g0tdn,grup,grdn,g00up,g00dn)

              END IF  ! IF ( .false. ) THEN

          end if ! if(ltau) then
