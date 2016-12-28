          ! dyn
          ! g00up g00dn
          ! gt0up gt0dn
          if( ltau .and. lmeasure_dyn .and. (.not.lupdate) ) then
              if( iwrap_nt(nt) .gt. 0 ) then
                  ! at stablization point, already have gt0,g0t from green_tau
              else
                  ! B(nt1,nt2) with nt1 >= nt2
                  nt1 = nt
                  nt2 = nt
                  ! G(t',0) = B(t',t) * G(t,0)
                  call Bmat_tau_R( nt1, nt2, gt0up, gt0dn)

                  ! G(0,t') = G(0,t) * B(t',t)^-1
                  call Bmatinv_tau_L( nt1, nt2, g0tup, g0tdn)
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

          end if ! if( ltau .and. lmeasure_dyn .and. (.not.lupdate) ) then
