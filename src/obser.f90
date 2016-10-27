module obser
  use blockc

  complex(dp), allocatable, dimension(:), save :: spinz
  complex(dp), allocatable, dimension(:,:), save :: spin_corrlt

  integer, allocatable, dimension(:), save :: isingzz_corrlt
  integer, allocatable, dimension(:,:), save :: isingzztau_corrlt

  complex(dp), save :: obs_bin(10)

  complex(dp), allocatable, dimension(:,:), save :: gtau_up, gtau_dn
  complex(dp), allocatable, dimension(:,:), save :: chiszsz, chijxjx, jttbin, j00bin

  contains

  subroutine allocate_obs
    implicit none
    allocate( spin_corrlt(lq,lq) )
    allocate( spinz(lq) )
    allocate( isingzz_corrlt(lq) )
    if(lsstau) allocate( isingzztau_corrlt(lq,ltrot) )
    if(ltau) then
        allocate( gtau_up(ndim,ltrot) )
#IFDEF SPINDOWN
        allocate( gtau_dn(ndim,ltrot) )
#ENDIF
        allocate( chiszsz(ndim,ltrot) )
        allocate( chijxjx(ndim,ltrot) )
        allocate( jttbin(ndim,ltrot) )
        allocate( j00bin(ndim,ltrot) )
    end if
  end subroutine allocate_obs

  subroutine deallocate_obs
    implicit none
    if(ltau) then
        deallocate( j00bin )
        deallocate( jttbin )
        deallocate( chijxjx )
        deallocate( chiszsz )
#IFDEF SPINDOWN
        deallocate( gtau_dn )
#ENDIF
        deallocate( gtau_up )
    end if
    if(lsstau) deallocate( isingzztau_corrlt )
    deallocate( isingzz_corrlt )
    deallocate( spinz )
    deallocate( spin_corrlt )
  end subroutine deallocate_obs

  subroutine obser_init
    implicit none
    nobs = 0
    spinz(:) = czero
    spin_corrlt(:,:) = czero
    isingzz_corrlt(:) = 0
    if(lsstau) isingzztau_corrlt(:,:) = 0
    obs_bin(:) = czero
    if(ltau) then
        gtau_up(:,:) = czero
#IFDEF SPINDOWN
        gtau_dn(:,:) = czero
#ENDIF
        chiszsz(:,:) = czero
        chijxjx(:,:) = czero
        jttbin(:,:) = czero
        j00bin(:,:) = czero
    end if

  end subroutine obser_init

  subroutine obser_equaltime(nt)
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    implicit none
    integer,intent(in) :: nt

    ! local 
    integer :: i, j, i0, i1, i2, i3, i4, i_1, n, nf, ist, ijs, imj
    real(dp) :: xi, xj
    complex(dp) :: szsz_tmp, zx, zy, zkint, zecoup, zne, zejs
    real(dp) :: rne_up, rne_dn, ehx
    integer :: order_branch
    complex(dp), external :: zthp

    nobs = nobs + 1

    !grup (i,j) = < c_i c^+_j >
    !grupc (i,j) = < c^+_i c_j >

    ! get grupc
!$OMP PARALLEL &
!$OMP PRIVATE ( i, j )
!$OMP DO
    do i = 1, ndim
        do j = 1, ndim
            grupc(j,i) = - grup(i,j)
        end do
        grupc(i,i) = grupc(i,i) + 1
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! get grdn and grdnc
#IFDEF SPINDOWN
!$OMP PARALLEL &
!$OMP PRIVATE ( i, j )
!$OMP DO
    do i = 1, ndim
        do j = 1, ndim
            grdnc(j,i) = - grdn(i,j)
        end do
        grdnc(i,i) = grdnc(i,i) + 1
    end do
!$OMP END DO
!$OMP END PARALLEL
#ELSE
!$OMP PARALLEL &
!$OMP PRIVATE ( i, xi, j, xj )
!$OMP DO
    do i = 1,ndim
       xi = 1.d0
       if (orblist(i) == 2 ) xi = -1.d0
       do j = 1,ndim
          xj = 1.d0
          if (orblist(j) == 2 ) xj= -1.d0
          grdn (i,j) = dcmplx(xi*xj, 0.d0) * dconjg ( grup (i,j) )
          grdnc(i,j) = dcmplx(xi*xj, 0.d0) * dconjg ( grupc(i,j) )
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
#ENDIF

    ! zne
    zne = czero
    do i = 1, ndim
        zne = zne + grupc(i,i) + grdnc(i,i)
    end do
    obs_bin(1) = obs_bin(1) + zne

    ! order_branch
    order_branch = 0
    do i = 1, ndim
        order_branch = order_branch + nsigl_u(i,nt)
    end do

    ! measure rne_up, rne_dn
    rne_up = zero
    rne_dn = zero
    do i = 1, ndim
        if( order_branch .gt. 0 ) then
            rne_up = rne_up + dble(grupc(i,i))
            rne_dn = rne_dn + dble(grdnc(i,i))
        else
            rne_up = rne_up + dble(grdnc(i,i))
            rne_dn = rne_dn + dble(grupc(i,i))
        end if
    end do
    obs_bin(2) = obs_bin(2) + dcmplx(rne_up, rne_dn)

    ! measure kinetic energy
    zkint = czero
    IF ( l .gt. 1 ) THEN
    do nf = 1,2
!$OMP PARALLEL &
!$OMP PRIVATE ( i_1, ist, i1, i2, i3, i4, zx, zy )
!$OMP DO REDUCTION ( + : zkint )
       do i_1 = 1,lq/4
          ist = i_1 + (nf - 1)*lq/4
          i1 = lthf(i_1,nf)
          i2 = nnlist(i1,1)
          i3 = nnlist(i1,5)
          i4 = nnlist(i1,2)
          
          zx = hopping_tmp(1,ist)
          zy = hopping_tmp(2,ist)
          zkint = zkint  +        zx *  grupc(i1,i2) + &
                           dconjg(zx) * grupc(i2,i1)
          zkint = zkint  +        zy *  grupc(i1,i4) + &
                           dconjg(zy) * grupc(i4,i1) 
#IFDEF SPINDOWN
          zx = hopping_tmp_dn(1,ist)
          zy = hopping_tmp_dn(2,ist)
          zkint = zkint  +        zx *  grdnc(i1,i2) + &
                           dconjg(zx) * grdnc(i2,i1)
          zkint = zkint  +        zy *  grdnc(i1,i4) + &
                           dconjg(zy) * grdnc(i4,i1) 
#ENDIF

          zy = hopping_tmp(3,ist)
          zkint = zkint  +        zy *  grupc(i2,i3) + &
                           dconjg(zy) * grupc(i3,i2) 
#IFDEF SPINDOWN
          zy = hopping_tmp_dn(3,ist)
          zkint = zkint  +        zy *  grdnc(i2,i3) + &
                           dconjg(zy) * grdnc(i3,i2)
#ENDIF

          zx = hopping_tmp(4,ist)
          zkint = zkint  +        zx *  grupc(i4,i3) + &
                           dconjg(zx) * grupc(i3,i4) 
#IFDEF SPINDOWN
          zx = hopping_tmp_dn(4,ist)
          zkint = zkint  +        zx *  grdnc(i4,i3) + &
                           dconjg(zx) * grdnc(i3,i4) 
#ENDIF
       end do
!$OMP END DO
!$OMP END PARALLEL
    end do
    ELSE
        zkint = zkint + dcmplx(-4.d0*rt,0.d0) * ( grupc(1,1) + grdnc(1,1) )
    ENDIF
    obs_bin(3) = obs_bin(3) + zkint + dconjg(zkint)  ! layer 1 and layer 2

    ! zecoup
    zecoup = czero
    do i = 1, lq
        zecoup = zecoup + ( grupc(i,i) - grdnc(i,i) ) * nsigl_u(i,nt)
    end do
    zecoup = zecoup*dcmplx(-rhub*0.5d0, 0.d0)  ! note 0.5 for fermion spin 1/2
    obs_bin(4) = obs_bin(4) + zecoup + dconjg(zecoup) ! layer 1 and layer 2

    ! zejs
    ijs = 0
    do nf = 1,2
!$OMP PARALLEL &
!$OMP PRIVATE ( i_1, i0, i1, i2, i3, i4 )
!$OMP DO REDUCTION ( + : ijs )
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
!$OMP END DO
!$OMP END PARALLEL
    end do
    zejs = dcmplx( dble(ijs)*(-js), 0.d0 )
    obs_bin(5) = obs_bin(5) + zejs

    ! ehx
    ehx = zero
    do i = 1, lq
        if ( nsigl_u(i,mod(nt,ltrot)+1) .eq. nsigl_u(i,nt) ) then
            ehx = ehx + tanhdth
        else
            ehx = ehx + cothdth
        end if
    end do
    obs_bin(6) = obs_bin(6) + dcmplx( -hx*ehx, 0.d0 )

    ! uncomment if needed
    !!!!! measure S(pi,pi)
    !!!!! S(pi,pi) = 1/L^2 \sum_ij <  1/2 (ni_up -ni_dn) *  1/2 ( nj_up - nj_dn ) >  exp( i * (pi,pi) * r(i-j) )
    !!!!!   ( ni_up -ni_dn ) * ( nj_up - nj_dn ) 
    !!!!! = ni_up * nj_up + ni_dn * nj_dn - ni_dn * nj_up - ni_up * nj_dn
    !!!!! = grupc(i,i)*grupc(j,j) + grupc(i,j)*grup(i,j) + grdnc(i,i)*grdnc(j,j) + grdnc(i,j)*grdn(i,j) - grdnc(i,i)*grupc(j,j) - grupc(i,i)*grdnc(j,j)
    !!!!!                                                  
    !!!!do i = 1, lq
    !!!!    do j = 1, lq
    !!!!        szsz_tmp = grupc(i,i)*grupc(j,j) + grupc(i,j)*grup(i,j) + &
    !!!!                   grdnc(i,i)*grdnc(j,j) + grdnc(i,j)*grdn(i,j) - &
    !!!!                   grdnc(i,i)*grupc(j,j) - grupc(i,i)*grdnc(j,j)
    !!!!        spin_corrlt(j,i) = spin_corrlt(j,i) + szsz_tmp
    !!!!    end do
    !!!!end do

    !!!!! the spinz, for calculate the backgroud
    !!!!! ni_up - ni_dn
    !!!!do i = 1, lq
    !!!!    spinz(i) = spinz(i) +  grupc(i,i) - grdnc(i,i)
    !!!!end do

    ! measure isingzz_corrlt
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i, imj )
!$OMP DO REDUCTION ( + : isingzz_corrlt )
    do j = 1, lq
        do i = 1, lq
            imj = latt_imj(i,j)
            isingzz_corrlt(imj) = isingzz_corrlt(imj) + nsigl_u(i,nt) * nsigl_u(j,nt)
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! isingzztau_corrlt
    ! < s_i(tau) s_j(0) >
    if( lsstau ) then
!$OMP PARALLEL &
!$OMP PRIVATE ( n, j, i, imj )
!$OMP DO
    do n = 1, ltrot
        do j = 1, lq
            do i = 1, lq
                imj = latt_imj(i,j)
                isingzztau_corrlt(imj,n) = isingzztau_corrlt(imj,n) + nsigl_u(i,mod(n+nt-2,ltrot)+1) * nsigl_u(j,nt)
                !isingzztau_corrlt(imj,n) = isingzztau_corrlt(imj,n) + nsigl_u(i,n) * nsigl_u(j,1)
            end do
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    end if

  end subroutine obser_equaltime

  subroutine obsert(nt, grt0_up, grt0_dn, gr0t_up, gr0t_dn, grtt_up, grtt_dn, gr00_up, gr00_dn)
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    implicit none
    integer, intent(in) :: nt
    complex(dp), dimension(ndim,ndim), intent(in) :: grt0_up, grt0_dn, gr0t_up, gr0t_dn, grtt_up, grtt_dn, gr00_up, gr00_dn
    complex(dp) :: ztmp, ztmp1, ztmp2

    ! local 
    integer :: i, j, imj, it, order_branch, iax, imx, jax, jmx

!$OMP PARALLEL &
!$OMP PRIVATE ( j, jax, jmx, i, imj, iax, imx, ztmp, ztmp1, ztmp2 )
!$OMP DO REDUCTION ( + : gtau_up, gtau_dn, chiszsz, chijxjx)
        do j = 1, lq
            jax = nnlist(j,1)
            jmx = nnlist(j,3)
            do i = 1, lq
                imj  = latt_imj(i,j)
                iax = nnlist(i,1) ! i+x
                imx = nnlist(i,3) ! i-x

                gtau_up(imj,nt) = gtau_up(imj,nt) + grt0_up(i,j)
                gtau_dn(imj,nt) = gtau_dn(imj,nt) + grt0_dn(i,j)

                ztmp = ( grtt_up(i,i) - grtt_dn(i,i) ) * ( gr00_up(j,j) - gr00_dn(j,j) ) * chalf - &
                       ( gr0t_up(j,i)*grt0_up(i,j) + gr0t_dn(j,i)*grt0_dn(i,j) ) * cquarter
                chiszsz(imj,nt) = chiszsz(imj,nt) + ztmp + dconjg(ztmp)

                ztmp = ( gr0t_up(jax,i) - gr0t_up(jmx,i) ) * ( grt0_up(iax,j) - grt0_up(imx,j) ) +  &
                       ( gr0t_dn(jax,i) - gr0t_dn(jmx,i) ) * ( grt0_dn(iax,j) - grt0_dn(imx,j) )
                ztmp1 = grtt_up(iax,i) - grtt_up(imx,i) + grtt_dn(iax,i) - grtt_dn(imx,i)
                ztmp2 = gr00_up(jax,j) - gr00_up(jmx,j) + gr00_dn(jax,j) - gr00_dn(jmx,j)
                chijxjx(imj,nt) = chijxjx(imj,nt) + ztmp + dconjg(ztmp) - ( ztmp1 + dconjg(ztmp1) ) * ( ztmp2 + dconjg(ztmp2) )

            end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!!!      do i = 1, lq
!!!          iax = nnlist(i,1) ! i+x
!!!          imx = nnlist(i,3) ! i-x
!!!          ztmp1 = - grtt_up(iax,i) + grtt_up(imx,i) - grtt_dn(iax,i) + grtt_dn(imx,i)
!!!          ztmp2 = - gr00_up(iax,i) + gr00_up(imx,i) - gr00_dn(iax,i) + gr00_dn(imx,i)
!!!          jttbin(i,nt) = jttbin(i,nt) + ztmp1 + dconjg(ztmp1)
!!!          j00bin(i,nt) = j00bin(i,nt) + ztmp2 + dconjg(ztmp2)
!!!      end do
  end subroutine obsert

end module obser
