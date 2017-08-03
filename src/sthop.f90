subroutine sthop

  use blockc
  implicit none
  
  ! local
#IFDEF BREAKUP_T
  integer, parameter :: nch = 4
  integer, parameter :: nch2 = 2
  integer, parameter :: nch3 = 4
  complex (dp) ::  zx, zy, zz, zzm, zx3, zy3, z0, z1
  complex(dp) :: hlp2(nch,nch), hlp1(nch,nch)
  complex(dp) :: hlpc2(nch2,nch2), hlpc1(nch2,nch2)
  complex(dp) :: hlpd3(nch3,nch3), hlpd1(nch3,nch3)
  real(dp) :: wc(nch)
  real(dp) :: wc2(nch2)
  real(dp) :: wc3(nch3)
  integer :: nf, i_1, ist, i1, i2, i3, i4, i5, i6, i7, i8, m, n, i, j, nf_tmp
  complex(dp), external :: zthp

  ! checkboard breakup for hopping term
  do nf = 1,2
     do i_1 = 1,lq/4
        ist = i_1 + (nf - 1)*lq/4
        i1 = lthf(i_1,nf)
        i2 = nnlist(i1,1)
        i3 = nnlist(i1,5)
        i4 = nnlist(i1,2)
        i5 = nnlist(i1,8)
        i6 = nnlist(i2,1)
        i7 = nnlist(i3,5)
        i8 = nnlist(i4,2)
        
        hlp2 = dcmplx(0.d0,0.d0)
        
        n = i1
        zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
        zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
        zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
        hlp2(1,2) =        zx
        hlp2(2,1) = dconjg(zx)
        hlp2(1,4) =        zy
        hlp2(4,1) = dconjg(zy)
        hlp2(1,3) =        zz
        hlp2(3,1) = dconjg(zz)
        !hopping_tmp(1,ist)=zx
        !hopping_tmp(2,ist)=zy

        n = i2
        zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
        zz = -dcmplx(rt2,0.d0)* zthp(n,-1,1,xmag,flux_x,flux_y)
        hlp2(2,3) =        zy
        hlp2(3,2) = dconjg(zy)
        hlp2(2,4) =        zz
        hlp2(4,2) = dconjg(zz)
        !hopping_tmp(3,ist)=zy


        n = i4
        zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
        hlp2(4,3) =        zx
        hlp2(3,4) = dconjg(zx)
        !hopping_tmp(4,ist)=zx

#IFDEF TEST
        if(irank.eq.0)then

        write(fout,*)
        write(fout,'(a)') 'hlp2 ='
        do i = 1, nch
            write(fout, "(4('(',2f7.3,' )'))" ) hlp2(:,i)
        end do

        end if
#ENDIF

        ! add the chemical potential term
        do i = 1, nch
            !hlp2(i,i) = hlp2(i,i) - dcmplx(0.25d0*mu,0.d0) ! 0.25 because every site count 4 times
            hlp2(i,i) = hlp2(i,i) - dcmplx(mu/6.d0,0.d0) ! 1/6 because every site count 6 times
            if(mod(i,2).eq.1) then
                !hlp2(i,i) = hlp2(i,i) - dcmplx(0.25d0*muA,0.d0)
                hlp2(i,i) = hlp2(i,i) - dcmplx(muA/6.d0,0.d0)
            else
                !hlp2(i,i) = hlp2(i,i) - dcmplx(0.25d0*muB,0.d0)
                hlp2(i,i) = hlp2(i,i) - dcmplx(muB/6.d0,0.d0)
            end if
        end do

        call s_eig_he(nch,nch,hlp2,wc,hlp1)

        do i = 1,nch
           do j = 1,nch
              z0 = dcmplx(0.d0,0.d0)
              z1 = dcmplx(0.d0,0.d0)
              do m = 1,nch
                 z0 = z0 +  hlp1(i,m) * dcmplx(dexp(-dtau *wc(m)),0.d0) * dconjg(hlp1(j,m))
                 z1 = z1 +  hlp1(i,m) * dcmplx(dexp( dtau *wc(m)),0.d0) * dconjg(hlp1(j,m))
              enddo
              urt  (ist,i,j) = z0
              urtm1(ist,i,j) = z1
           enddo
        enddo
     enddo
  enddo

  !For 3rd_nearest hopping family
  do nf = 1,8
     do i_1 = 1,lq/16
        ist = i_1 + (nf - 1)*lq/16
        i1 = lthf3(i_1,nf)
        i2 = nnlist(i1,1)
        i3 = nnlist(i1,5)
        i4 = nnlist(i1,2)
        i5 = nnlist(i1,8)
        i6 = nnlist(i2,1)
        i7 = nnlist(i3,5)
        i8 = nnlist(i4,2)
        
        hlpd3 = dcmplx(0.d0,0.d0)
        
        n = i1
        zx = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
        zy = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
        hlpd3(1,2) =        zx
        hlpd3(2,1) = dconjg(zx)
        hlpd3(1,4) =        zy
        hlpd3(4,1) = dconjg(zy)
                                                               
        n = i6
        zy = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
        hlpd3(2,3) =        zy
        hlpd3(3,2) = dconjg(zy)
                                                               
        n = i8
        zx = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
        hlpd3(4,3) =        zx
        hlpd3(3,4) = dconjg(zx)

        do i = 1, nch3
            hlpd3(i,i) = hlpd3(i,i) - dcmplx(mu/6.d0,0.d0) ! 1/6 because every site count 6 times
            if(mod(i,2).eq.1) then
                hlpd3(i,i) = hlpd3(i,i) - dcmplx(muA/6.d0,0.d0)
            else
                hlpd3(i,i) = hlpd3(i,i) - dcmplx(muB/6.d0,0.d0)
            end if
        end do

        call s_eig_he(nch3,nch3,hlpd3,wc3,hlpd1)
                                                                                               
        do i = 1,nch3
           do j = 1,nch3
              z0 = dcmplx(0.d0,0.d0)
              z1 = dcmplx(0.d0,0.d0)
              do m = 1,nch
                 z0 = z0 +  hlpd1(i,m) * dcmplx(dexp(-dtau *wc3(m)),0.d0) * dconjg(hlpd1(j,m))
                 z1 = z1 +  hlpd1(i,m) * dcmplx(dexp( dtau *wc3(m)),0.d0) * dconjg(hlpd1(j,m))
              end do
              urtd  (ist,i,j) = z0
              urtdm1(ist,i,j) = z1
           end do
        end do
     enddo
  enddo


  !For the next_nearest hopping family
  do nf = 1, 4
     do i_1 = 1, lq/4
        nf_tmp = (nf-1)/2 + 1
        ist = i_1 + (nf - 1)*lq/4
        i1 = lthf2(i_1,nf_tmp)
        i2 = nnlist(i1,5)
        i3 = nnlist(i1,1)
        i4 = nnlist(i3,6)
        hlpc2 = dcmplx(0.d0,0.d0)
        
        if ( mod(nf,2) .eq. 0 ) then
             n = i1
             zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
             hlpc2(1,2) =        zz
             hlpc2(2,1) = dconjg(zz)
             !hopping_tmp(1,ist)=zx
             !hopping_tmp(2,ist)=zy
        else
             n = i3
             zz = -dcmplx(rt2,0.d0)* zthp(n,-1,1,xmag,flux_x,flux_y)
             hlpc2(1,2) =        zz
             hlpc2(2,1) = dconjg(zz)
             !hopping_tmp(1,ist)=zx
             !hopping_tmp(2,ist)=zy
        end if

#IFDEF TEST
        if(irank.eq.0)then

        write(fout,*)
        write(fout,'(a)') 'hlpc2 ='
        do i = 1, nch2
            write(fout, "(2('(',2f7.3,' )'))" ) hlpc2(:,i)
        end do

        end if
#ENDIF

        ! add the chemical potential term
        do i = 1, nch2
            !hlpc2(i,i) = hlpc2(i,i) - dcmplx(0.25d0*mu,0.d0) ! 0.25 because every site count 4 times
            hlpc2(i,i) = hlpc2(i,i) - dcmplx(mu/6.d0,0.d0) ! 1/6 because every site count 6 times
            if(mod(i,2).eq.1) then
                !hlpc2(i,i) = hlpc2(i,i) - dcmplx(0.25d0*muA,0.d0)
                hlpc2(i,i) = hlpc2(i,i) - dcmplx(muA/6.d0,0.d0)
            else
                !hlpc2(i,i) = hlpc2(i,i) - dcmplx(0.25d0*muB,0.d0)
                hlpc2(i,i) = hlpc2(i,i) - dcmplx(muB/6.d0,0.d0)
            end if
        end do
                                                                                                        
        call s_eig_he(nch2,nch2,hlpc2,wc2,hlpc1)
                                                                                                        
        do i = 1,nch2
           do j = 1,nch2
              z0 = dcmplx(0.d0,0.d0)
              z1 = dcmplx(0.d0,0.d0)
              do m = 1,nch2
                 z0 = z0 +  hlpc1(i,m) * dcmplx(dexp(-dtau *wc2(m)),0.d0) * dconjg(hlpc1(j,m))
                 z1 = z1 +  hlpc1(i,m) * dcmplx(dexp( dtau *wc2(m)),0.d0) * dconjg(hlpc1(j,m))          
                 urtc  (ist,i,j) = z0
                 urtcm1(ist,i,j) = z1
              enddo
           enddo
        enddo
     enddo
  enddo

  !calculate hopping_tmp for spin up
  do n = 1,lq                                                
     zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
     zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
     zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
     zzm = -dcmplx(rt2,0.d0)* zthp(n,1,-1,xmag,flux_x,flux_y)
     zx3 = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
     zy3 = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
     hopping_tmp(1,n)=zx  
     hopping_tmp(2,n)=zy
     hopping_tmp(3,n)=zz
     hopping_tmp(4,n)=zzm
     hopping_tmp(5,n)=zx3
     hopping_tmp(6,n)=zy3
  end do                 

#IFDEF SPINDOWN
  do nf = 1,2
     do i_1 = 1,lq/4
        ist = i_1 + (nf - 1)*lq/4
        i1 = lthf(i_1,nf)
        i2 = nnlist(i1,1)
        i3 = nnlist(i1,5)
        i4 = nnlist(i1,2)
        i5 = nnlist(i1,8)
        i6 = nnlist(i2,1)
        i7 = nnlist(i3,5)
        i8 = nnlist(i4,2)
        
        hlp2 = dcmplx(0.d0,0.d0)
        
        n = i1
        zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
        zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
        zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
        hlp2(1,2) =        zx
        hlp2(2,1) = dconjg(zx)
        hlp2(1,4) =        zy
        hlp2(4,1) = dconjg(zy)
        hlp2(1,3) =        zz
        hlp2(3,1) = dconjg(zz)
        !hopping_tmp_dn(1,ist)=zx
        !hopping_tmp_dn(2,ist)=zy

        n = i2
        zy = -dcmplx(rt,0.d0)*zthp(n,0,1,xmag,flux_x,flux_y)
        zz = -dcmplx(rt2,0.d0)* zthp(n,-1,1,xmag,flux_x,flux_y)
        hlp2(2,3) =        zy
        hlp2(3,2) = dconjg(zy)
        hlp2(2,4) =        zz
        hlp2(4,2) = dconjg(zz)
        !hopping_tmp_dn(3,ist)=zy


        n = i4
        zx = -dcmplx(rt,0.d0)*zthp(n,1,0,xmag,flux_x,flux_y)
        hlp2(4,3) =        zx
        hlp2(3,4) = dconjg(zx)
        !hopping_tmp_dn(4,ist)=zx

        !3rd nearest hopping
        n = i1
        zx = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
        zy = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
        hlpd3(1,2) =        zx
        hlpd3(2,1) = dconjg(zx)
        hlpd3(1,4) =        zy
        hlpd3(4,1) = dconjg(zy)
                                                               
        n = i6
        zy = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
        hlpd3(2,3) =        zy
        hlpd3(3,2) = dconjg(zy)
                                                               
        n = i8
        zx = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
        hlpd3(4,3) =        zx
        hlpd3(3,4) = dconjg(zx)

#IFDEF TEST
        if(irank.eq.0)then

        write(fout,*)
        write(fout,'(a)') 'hlp2 ='
        do i = 1, nch
            write(fout, "(4('(',2f7.3,' )'))" ) hlp2(:,i)
        end do

        end if
#ENDIF

        ! add the chemical potential term
        do i = 1, nch
            !hlp2(i,i) = hlp2(i,i) - dcmplx(0.25d0*mu,0.d0) ! 0.25 because every site count 4 times
            hlp2(i,i) = hlp2(i,i) - dcmplx(mu/6.d0,0.d0) ! 1/6 because every site count 6 times
            if(mod(i,2).eq.1) then
                !hlp2(i,i) = hlp2(i,i) - dcmplx(0.25d0*muA,0.d0)
                hlp2(i,i) = hlp2(i,i) - dcmplx(muA/6.d0,0.d0)
            else
                !hlp2(i,i) = hlp2(i,i) - dcmplx(0.25d0*muB,0.d0)
                hlp2(i,i) = hlp2(i,i) - dcmplx(muB/6.d0,0.d0)
            end if
        end do

        call s_eig_he(nch,nch,hlp2,wc,hlp1)

        do i = 1,nch
           do j = 1,nch
              z0 = dcmplx(0.d0,0.d0)
              z1 = dcmplx(0.d0,0.d0)
              do m = 1,nch
                 z0 = z0 +  hlp1(i,m) * dcmplx(dexp(-dtau *wc(m)),0.d0) * dconjg(hlp1(j,m))
                 z1 = z1 +  hlp1(i,m) * dcmplx(dexp( dtau *wc(m)),0.d0) * dconjg(hlp1(j,m))
              enddo
              urt_dn  (ist,i,j) = z0
              urtm1_dn(ist,i,j) = z1
           enddo
        enddo
     enddo
  enddo

  !For 3rd_nearest hopping family
  do nf = 1,8
     do i_1 = 1,lq/16
        ist = i_1 + (nf - 1)*lq/16
        i1 = lthf3(i_1,nf)
        i2 = nnlist(i1,1)
        i3 = nnlist(i1,5)
        i4 = nnlist(i1,2)
        i5 = nnlist(i1,8)
        i6 = nnlist(i2,1)
        i7 = nnlist(i3,5)
        i8 = nnlist(i4,2)
        
        hlpd3 = dcmplx(0.d0,0.d0)
        
        n = i1
        zx = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
        zy = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
        hlpd3(1,2) =        zx
        hlpd3(2,1) = dconjg(zx)
        hlpd3(1,4) =        zy
        hlpd3(4,1) = dconjg(zy)
                                                               
        n = i6
        zy = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
        hlpd3(2,3) =        zy
        hlpd3(3,2) = dconjg(zy)
                                                               
        n = i8
        zx = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
        hlpd3(4,3) =        zx
        hlpd3(3,4) = dconjg(zx)
                                                                                                  
        do i = 1, nch3
            hlpd3(i,i) = hlpd3(i,i) - dcmplx(mu/6.d0,0.d0) ! 1/6 because every site count 6 times
            if(mod(i,2).eq.1) then
                hlpd3(i,i) = hlpd3(i,i) - dcmplx(muA/6.d0,0.d0)
            else
                hlpd3(i,i) = hlpd3(i,i) - dcmplx(muB/6.d0,0.d0)
            end if
        end do
                                                                                                  
        call s_eig_he(nch3,nch3,hlpd3,wc3,hlpd1)
                                                                                               
        do i = 1,nch3
           do j = 1,nch3
              z0 = dcmplx(0.d0,0.d0)
              z1 = dcmplx(0.d0,0.d0)
              do m = 1,nch
                 z0 = z0 +  hlpd1(i,m) * dcmplx(dexp(-dtau *wc3(m)),0.d0) * dconjg(hlpd1(j,m))
                 z1 = z1 +  hlpd1(i,m) * dcmplx(dexp( dtau *wc3(m)),0.d0) * dconjg(hlpd1(j,m))
              end do
              urtd_dn  (ist,i,j) = z0
              urtdm1_dn(ist,i,j) = z1
           end do
        end do
     enddo
  enddo

  !For the next_nearest hopping family
  do nf = 1, 4
     do i_1 = 1, lq/4
        nf_tmp = (nf-1)/2 + 1
        ist = i_1 + (nf - 1)*lq/4
        i1 = lthf2(i_1,nf_tmp)
        i2 = nnlist(i1,5)
        i3 = nnlist(i1,1)
        i4 = nnlist(i3,6)
        hlpc2 = dcmplx(0.d0,0.d0)
        
        if ( mod(nf,2) .eq. 0 ) then
             n = i1
             zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
             hlpc2(1,2) =        zz
             hlpc2(2,1) = dconjg(zz)
             !hopping_tmp(1,ist)=zx
             !hopping_tmp(2,ist)=zy
        else
             n = i3
             zz = -dcmplx(rt2,0.d0)* zthp(n,-1,1,xmag,flux_x,flux_y)
             hlpc2(1,2) =        zz
             hlpc2(2,1) = dconjg(zz)
             !hopping_tmp(1,ist)=zx
             !hopping_tmp(2,ist)=zy
        end if

#IFDEF TEST
        if(irank.eq.0)then

        write(fout,*)
        write(fout,'(a)') 'hlpc2 ='
        do i = 1, nch2
            write(fout, "(2('(',2f7.3,' )'))" ) hlpc2(:,i)
        end do

        end if
#ENDIF

        ! add the chemical potential term
        do i = 1, nch2
            !hlpc2(i,i) = hlpc2(i,i) - dcmplx(0.25d0*mu,0.d0) ! 0.25 because every site count 4 times
            hlpc2(i,i) = hlpc2(i,i) - dcmplx(mu/6.d0,0.d0) ! 1/6 because every site count 6 times
            if(mod(i,2).eq.1) then
                !hlpc2(i,i) = hlpc2(i,i) - dcmplx(0.25d0*muA,0.d0)
                hlpc2(i,i) = hlpc2(i,i) - dcmplx(muA/6.d0,0.d0)
            else
                !hlpc2(i,i) = hlpc2(i,i) - dcmplx(0.25d0*muB,0.d0)
                hlpc2(i,i) = hlpc2(i,i) - dcmplx(muB/6.d0,0.d0)
            end if
        end do
                                                                                              
        call s_eig_he(nch2,nch2,hlpc2,wc2,hlpc1)
                                                                                              
        do i = 1,nch2
           do j = 1,nch2
              z0 = dcmplx(0.d0,0.d0)
              z1 = dcmplx(0.d0,0.d0)
              do m = 1,nch2
                 z0 = z0 +  hlpc1(i,m) * dcmplx(dexp(-dtau *wc2(m)),0.d0) * dconjg(hlpc1(j,m))
                 z1 = z1 +  hlpc1(i,m) * dcmplx(dexp( dtau *wc2(m)),0.d0) * dconjg(hlpc1(j,m)) 
                 urtc_dn  (ist,i,j) = z0
                 urtcm1_dn(ist,i,j) = z1
              enddo
           enddo
        enddo
     enddo
  enddo

  !calculate hopping_tmp for spin down
  do n = 1,lq                                                
     zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
     zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
     zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
     zzm = -dcmplx(rt2,0.d0)* zthp(n,1,-1,xmag,flux_x,flux_y)
     zx3 = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
     zy3 = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
     hopping_tmp_dn(1,n)=zx  
     hopping_tmp_dn(2,n)=zy
     hopping_tmp_dn(3,n)=zz
     hopping_tmp_dn(4,n)=zzm
     hopping_tmp_dn(5,n)=zx3
     hopping_tmp_dn(6,n)=zy3
  end do

#ENDIF

#ELSE
  ! local
  integer :: i_1, ist, i1, i2, i3, i4, i5, i6, i7, i8, n, i, j, m, nf
  complex(dp) :: zx, zy, zz, zzm, zx3, zy3, z0, z1
  complex(dp), dimension(:,:), allocatable :: hvec, hmat
  real(dp), dimension(:), allocatable :: heig
  real(dp) :: en_free
  complex(dp), external :: zthp

  allocate( hmat(ndim,ndim) )
  allocate( hvec(ndim,ndim) )
  allocate( heig(ndim) )

  hmat(:,:) = dcmplx(0.d0,0.d0)
  IF( l .gt. 1 ) THEN
      do n = 1,lq

            i1 = n
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            i5 = nnlist(i1,8)
            i6 = nnlist(i2,1)
            i7 = nnlist(i3,5)
            i8 = nnlist(i4,2)
            
            zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
            zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
            zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
            zzm = -dcmplx(rt2,0.d0)* zthp(n,1,-1,xmag,flux_x,flux_y)
            zx3 = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
            zy3 = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
            hmat(i1,i2) = hmat(i1,i2) +        zx
            hmat(i2,i1) = hmat(i2,i1) + dconjg(zx)
            hmat(i1,i3) = hmat(i1,i3) +        zz
            hmat(i3,i1) = hmat(i3,i1) + dconjg(zz)
            hmat(i1,i4) = hmat(i1,i4) +        zy
            hmat(i4,i1) = hmat(i4,i1) + dconjg(zy)
            hmat(i1,i5) = hmat(i1,i5) +        zzm
            hmat(i5,i1) = hmat(i5,i1) + dconjg(zzm)
            hmat(i1,i6) = hmat(i1,i6) +        zx3
            hmat(i6,i1) = hmat(i6,i1) + dconjg(zx3)
            hmat(i1,i8) = hmat(i1,i8) +        zy3
            hmat(i8,i1) = hmat(i8,i1) + dconjg(zy3)
            hopping_tmp(1,n)=zx      
            hopping_tmp(2,n)=zy
            hopping_tmp(3,n)=zz
            hopping_tmp(4,n)=zzm
            hopping_tmp(5,n)=zx3
            hopping_tmp(6,n)=zy3

      end do
  ELSE
      hmat(1,1) = dcmplx(-4.d0*rt, 0.d0 )
  END IF
#IFDEF TEST
  write(fout,*)
  write(fout, '(a)') ' hmat = '
  do i = 1, ndim
      write(fout,'(40(2f8.4))') hmat(i,:)
  end do
#ENDIF
  ! add the chemical potential term
  do i = 1, ndim
      hmat(i,i) = hmat(i,i) - dcmplx(mu,0.d0)
      if ( mod( list(i,1)+list(i,2), 2 ) .eq. 0 ) then
          hmat(i,i) = hmat(i,i) - dcmplx(muA,0.d0)
      else
          hmat(i,i) = hmat(i,i) - dcmplx(muB,0.d0)
      end if
  end do

  call s_eig_he(ndim,ndim,hmat,heig,hvec)
  do i = 1,ndim
     do j = 1,ndim
        z0 = dcmplx(0.d0,0.d0)
        z1 = dcmplx(0.d0,0.d0)
        do m = 1,ndim
           z0 = z0 +  hvec(i,m) * dcmplx(dexp(-dtau *heig(m)),0.d0) * dconjg(hvec(j,m))
           z1 = z1 +  hvec(i,m) * dcmplx(dexp( dtau *heig(m)),0.d0) * dconjg(hvec(j,m))
        enddo
        urt  (i,j) = z0
        urtm1(i,j) = z1
     enddo
  enddo

  ! urt
#IFDEF TEST

  write(fout,*)
  write(fout,'(a)') 'urt(:,:) = '
  do i = 1, ndim
      write(fout,'(40(2f9.5))') urt(:,i)
  end do

#ENDIF

  if( irank.eq.0 ) then
      write(fout,*)
      write(fout,'(a)') ' heig(:) = '
      do i = 1, ndim
          write(fout,'(e16.8)') heig(i)
      end do

      en_free = 0.d0
      do i = 1, ndim/2
          en_free = en_free + heig(i)
      end do
      write(fout,*)
      write(fout,'(a,e16.8)') ' half-filling en_free = ', 2.d0*en_free   ! 2 spin flavor
  end if



#IFDEF SPINDOWN
  hmat(:,:) = dcmplx(0.d0,0.d0)
  IF( l .gt. 1 ) THEN
      do n = 1, lq 

         i1 = n
         i2 = nnlist(i1,1)
         i3 = nnlist(i1,5)
         i4 = nnlist(i1,2)
         i5 = nnlist(i1,8)
         i6 = nnlist(i2,1)
         i7 = nnlist(i3,5)
         i8 = nnlist(i4,2)
         
         zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
         zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
         zz = -dcmplx(rt2,0.d0)* zthp(n,1,1,xmag,flux_x,flux_y)
         zzm = -dcmplx(rt2,0.d0)* zthp(n,1,-1,xmag,flux_x,flux_y)
         zx3 = -dcmplx(rt3,0.d0)* zthp(n,2,0,xmag,flux_x,flux_y)
         zy3 = -dcmplx(rt3,0.d0)* zthp(n,0,2,xmag,flux_x,flux_y)
         hmat(i1,i2) = hmat(i1,i2) +        zx
         hmat(i2,i1) = hmat(i2,i1) + dconjg(zx)
         hmat(i1,i3) = hmat(i1,i3) +        zz
         hmat(i3,i1) = hmat(i3,i1) + dconjg(zz)
         hmat(i1,i4) = hmat(i1,i4) +        zy
         hmat(i4,i1) = hmat(i4,i1) + dconjg(zy)
         hmat(i1,i5) = hmat(i1,i5) +        zzm
         hmat(i5,i1) = hmat(i5,i1) + dconjg(zzm)
         hmat(i1,i6) = hmat(i1,i6) +        zx3
         hmat(i6,i1) = hmat(i6,i1) + dconjg(zx3)
         hmat(i1,i8) = hmat(i1,i8) +        zy3
         hmat(i8,i1) = hmat(i8,i1) + dconjg(zy3)
         hopping_tmp_dn(1,n)=zx      
         hopping_tmp_dn(2,n)=zy
         hopping_tmp_dn(3,n)=zz
         hopping_tmp_dn(4,n)=zzm
         hopping_tmp_dn(5,n)=zx3 
         hopping_tmp_dn(6,n)=zy3

      end do
  ELSE
      hmat(1,1) = dcmplx(-4.d0*rt, 0.d0 )
  END IF
#IFDEF TEST
  write(fout,*)
  write(fout, '(a)') ' hmat = '
  do i = 1, ndim
      write(fout,'(40(2f8.4))') hmat(i,:)
  end do
#ENDIF
  ! add the chemical potential term
  do i = 1, ndim
      hmat(i,i) = hmat(i,i) - dcmplx(mu,0.d0)
      if ( mod( list(i,1)+list(i,2), 2 ) .eq. 0 ) then
          hmat(i,i) = hmat(i,i) - dcmplx(muA,0.d0)
      else
          hmat(i,i) = hmat(i,i) - dcmplx(muB,0.d0)
      end if
  end do

  call s_eig_he(ndim,ndim,hmat,heig,hvec)
  do i = 1,ndim
     do j = 1,ndim
        z0 = dcmplx(0.d0,0.d0)
        z1 = dcmplx(0.d0,0.d0)
        do m = 1,ndim
           z0 = z0 +  hvec(i,m) * dcmplx(dexp(-dtau *heig(m)),0.d0) * dconjg(hvec(j,m))
           z1 = z1 +  hvec(i,m) * dcmplx(dexp( dtau *heig(m)),0.d0) * dconjg(hvec(j,m))
        enddo
        urt_dn  (i,j) = z0
        urtm1_dn(i,j) = z1
     enddo
  enddo

  ! urt_dn
#IFDEF TEST

  write(fout,*)
  write(fout,'(a)') 'urt_dn(:,:) = '
  do i = 1, ndim
      write(fout,'(40(2f9.5))') urt_dn(:,i)
  end do

#ENDIF

  if( irank.eq.0 ) then
      write(fout,*)
      write(fout,'(a)') ' heig(:) = '
      do i = 1, ndim
          write(fout,'(e16.8)') heig(i)
      end do

      en_free = 0.d0
      do i = 1, ndim/2
          en_free = en_free + heig(i)
      end do
      write(fout,*)
      write(fout,'(a,e16.8)') ' half-filling en_free = ', 2.d0*en_free   ! 2 spin flavor
  end if
#ENDIF

  deallocate( heig )
  deallocate( hvec )

#ENDIF

end subroutine sthop
