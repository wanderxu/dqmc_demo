subroutine sthop

  use blockc
  implicit none
  
  ! local
#IFDEF BREAKUP_T
  integer, parameter :: nch = 4
  complex (dp) ::  zx, zy, z0, z1
  complex(dp) :: hlp2(nch,nch), hlp1(nch,nch)
  real(dp) :: wc(nch)
  integer :: nf, i_1, ist, i1, i2, i3, i4, m, n, i, j
  complex(dp), external :: zthp

  ! checkboard breakup for hopping term
  do nf = 1,2
     do i_1 = 1,lq/4
        ist = i_1 + (nf - 1)*lq/4
        i1 = lthf(i_1,nf)
        i2 = nnlist(i1,1)
        i3 = nnlist(i1,5)
        i4 = nnlist(i1,2)
        
        hlp2 = dcmplx(0.d0,0.d0)
        
        n = i1	
        zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
        zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
        hlp2(1,2) =        zx
        hlp2(2,1) = dconjg(zx)
        hlp2(1,4) =        zy
        hlp2(4,1) = dconjg(zy)
        hopping_tmp(1,ist)=zx
        hopping_tmp(2,ist)=zy

        n = i2
        zy = -dcmplx(rt,0.d0)*zthp(n,0,1,xmag,flux_x,flux_y)
        hlp2(2,3) =        zy
        hlp2(3,2) = dconjg(zy)
        hopping_tmp(3,ist)=zy


        n = i4
        zx = -dcmplx(rt,0.d0)*zthp(n,1,0,xmag,flux_x,flux_y)
        hlp2(4,3) =        zx
        hlp2(3,4) = dconjg(zx)
        hopping_tmp(4,ist)=zx

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
            hlp2(i,i) = hlp2(i,i) - dcmplx(0.5d0*mu,0.d0) ! 0.5 because every site count 2 times
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
#IFDEF SPINDOWN
  do nf = 1,2
     do i_1 = 1,lq/4
        ist = i_1 + (nf - 1)*lq/4
        i1 = lthf(i_1,nf)
        i2 = nnlist(i1,1)
        i3 = nnlist(i1,5)
        i4 = nnlist(i1,2)
        
        hlp2 = dcmplx(0.d0,0.d0)
        
        n = i1	
        zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
        zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
        hlp2(1,2) =        zx
        hlp2(2,1) = dconjg(zx)
        hlp2(1,4) =        zy
        hlp2(4,1) = dconjg(zy)
        hopping_tmp_dn(1,ist)=zx
        hopping_tmp_dn(2,ist)=zy

        n = i2
        zy = -dcmplx(rt,0.d0)*zthp(n,0,1,xmag,flux_x,flux_y)
        hlp2(2,3) =        zy
        hlp2(3,2) = dconjg(zy)
        hopping_tmp_dn(3,ist)=zy


        n = i4
        zx = -dcmplx(rt,0.d0)*zthp(n,1,0,xmag,flux_x,flux_y)
        hlp2(4,3) =        zx
        hlp2(3,4) = dconjg(zx)
        hopping_tmp_dn(4,ist)=zx

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
            hlp2(i,i) = hlp2(i,i) - dcmplx(0.5d0*mu,0.d0) ! 0.5 because every site count 2 times
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
#ENDIF

#ELSE
  ! local
  integer :: i_1, ist, i1, i2, i3, i4, n, i, j, m, nf
  complex(dp) :: zx, zy, z0, z1
  complex(dp), dimension(:,:), allocatable :: hvec, hmat
  real(dp), dimension(:), allocatable :: heig
  real(dp) :: en_free
  complex(dp), external :: zthp

  allocate( hmat(ndim,ndim) )
  allocate( hvec(ndim,ndim) )
  allocate( heig(ndim) )

  hmat(:,:) = dcmplx(0.d0,0.d0)
  IF( l .gt. 1 ) THEN
      do nf = 1,2
         do i_1 = 1,lq/4
            ist = i_1 + (nf - 1)*lq/4
            i1 = lthf(i_1,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            
            
            n = i1	
            zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
            zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
            hmat(i1,i2) = hmat(i1,i2) +        zx
            hmat(i2,i1) = hmat(i2,i1) + dconjg(zx)
            hmat(i1,i4) = hmat(i1,i4) +        zy
            hmat(i4,i1) = hmat(i4,i1) + dconjg(zy)
            hopping_tmp(1,ist)=zx
            hopping_tmp(2,ist)=zy

            n = i2
            zy = -dcmplx(rt,0.d0)*zthp(n,0,1,xmag,flux_x,flux_y)
            hmat(i2,i3) = hmat(i2,i3) +        zy
            hmat(i3,i2) = hmat(i3,i2) + dconjg(zy)
            hopping_tmp(3,ist)=zy

            n = i4
            zx = -dcmplx(rt,0.d0)*zthp(n,1,0,xmag,flux_x,flux_y)
            hmat(i4,i3) = hmat(i4,i3) +        zx
            hmat(i3,i4) = hmat(i3,i4) + dconjg(zx)
            hopping_tmp(4,ist)=zx
        end do
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
      do nf = 1,2
         do i_1 = 1,lq/4
            ist = i_1 + (nf - 1)*lq/4
            i1 = lthf(i_1,nf)
            i2 = nnlist(i1,1)
            i3 = nnlist(i1,5)
            i4 = nnlist(i1,2)
            
            
            n = i1	
            zx = -dcmplx(rt,0.d0)* zthp(n,1,0,xmag,flux_x,flux_y)
            zy = -dcmplx(rt,0.d0)* zthp(n,0,1,xmag,flux_x,flux_y)
            hmat(i1,i2) = hmat(i1,i2) +        zx
            hmat(i2,i1) = hmat(i2,i1) + dconjg(zx)
            hmat(i1,i4) = hmat(i1,i4) +        zy
            hmat(i4,i1) = hmat(i4,i1) + dconjg(zy)
            hopping_tmp_dn(1,ist)=zx
            hopping_tmp_dn(2,ist)=zy

            n = i2
            zy = -dcmplx(rt,0.d0)*zthp(n,0,1,xmag,flux_x,flux_y)
            hmat(i2,i3) = hmat(i2,i3) +        zy
            hmat(i3,i2) = hmat(i3,i2) + dconjg(zy)
            hopping_tmp_dn(3,ist)=zy

            n = i4
            zx = -dcmplx(rt,0.d0)*zthp(n,1,0,xmag,flux_x,flux_y)
            hmat(i4,i3) = hmat(i4,i3) +        zx
            hmat(i3,i4) = hmat(i3,i4) + dconjg(zx)
            hopping_tmp_dn(4,ist)=zx
        end do
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
