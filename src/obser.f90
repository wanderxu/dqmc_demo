module obser
  use blockc
  complex(dp), save :: obs_bin(10), pair_bin(19), high_pair_bin(4)
  complex(dp), allocatable, dimension(:,:), save :: gtau_up, gtau_dn
  complex(dp), allocatable, dimension(:,:), save :: chiszsz, chijxjx
  complex(dp), allocatable, dimension(:,:), save :: chijxjxaa, chijxjxab, chijxjxba, chijxjxbb

  contains

  subroutine allocate_obs
    implicit none
    if(ltau) then
        allocate( gtau_up(ndim,ltrot) )
#IFDEF SPINDOWN
        allocate( gtau_dn(ndim,ltrot) )
#ENDIF
        allocate( chiszsz(ndim,ltrot) )
        allocate( chijxjx(ndim,ltrot) )
        allocate( chijxjxaa(ndim,ltrot) )
        allocate( chijxjxab(ndim,ltrot) )
        allocate( chijxjxba(ndim,ltrot) )
        allocate( chijxjxbb(ndim,ltrot) )
    end if
  end subroutine allocate_obs

  subroutine deallocate_obs
    implicit none
    if(ltau) then
        deallocate( chijxjxbb )
        deallocate( chijxjxba )
        deallocate( chijxjxab )
        deallocate( chijxjxaa )
        deallocate( chijxjx )
        deallocate( chiszsz )
#IFDEF SPINDOWN
        deallocate( gtau_dn )
#ENDIF
        deallocate( gtau_up )
    end if
  end subroutine deallocate_obs

  subroutine obser_init
    implicit none
    nobs = 0
    obs_bin(:) = czero
    pair_bin(:) = czero
    high_pair_bin(:) = czero
    if(ltau) then
        gtau_up(:,:) = czero
#IFDEF SPINDOWN
        gtau_dn(:,:) = czero
#ENDIF
        chiszsz(:,:) = czero
        chijxjx(:,:) = czero
        chijxjxbb(:,:) = czero
        chijxjxba(:,:) = czero
        chijxjxab(:,:) = czero
        chijxjxaa(:,:) = czero
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

    integer :: iax, imx, iay, imy, jax, jmx, jay, jmy
    complex(dp) :: Cotss, Cost0s, Cost1s, Cnnt1ss, Cnnt1sd, Cnnt0ss, Cnnt0sd, Cnnst0s, Cnnst0d, Cnnst1s, Cnnst1d, &
                   Cnnt1t0px, Cnnt1t0pxipy, Cnnt0t0px, Cnnt0t0pxipy, Cnnt1t1px, Cnnt1t1pxipy, Cnnt0t1px, Cnnt0t1pxipy 

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
    obs_bin(2) = obs_bin(2) + zkint + dconjg(zkint)  ! layer 1 and layer 2

    ! zecoup
    zecoup = czero
    do i = 1, lq
        zecoup = zecoup + ( grupc(i,i) - grdnc(i,i) ) * nsigl_u(i,nt)
    end do
    zecoup = zecoup*dcmplx(-rhub*0.5d0, 0.d0)  ! note 0.5 for fermion spin 1/2
    obs_bin(3) = obs_bin(3) + zecoup + dconjg(zecoup) ! layer 1 and layer 2

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
    obs_bin(4) = obs_bin(4) + zejs

    ! ehx
    ehx = zero
    do i = 1, lq
        if ( nsigl_u(i,mod(nt,ltrot)+1) .eq. nsigl_u(i,nt) ) then
            ehx = ehx + tanhdth
        else
            ehx = ehx + cothdth
        end if
    end do
    obs_bin(5) = obs_bin(5) + dcmplx( -hx*ehx, 0.d0 )

    ! pairing
    Cotss = czero
    Cost0s = czero
    Cost1s = czero
    Cnnt1ss = czero
    Cnnt1sd = czero
    Cnnt0ss = czero
    Cnnt0sd = czero
    Cnnst0s = czero
    Cnnst0d = czero
    Cnnst1s = czero
    Cnnst1d = czero
    Cnnt1t0px = czero
    Cnnt1t0pxipy = czero
    Cnnt0t0px = czero
    Cnnt0t0pxipy = czero
    Cnnt1t1px = czero
    Cnnt1t1pxipy = czero
    Cnnt0t1px = czero
    Cnnt0t1pxipy = czero
    do j = 1, lq
        jax = nnlist(j,1)
        jmx = nnlist(j,3)
        jay = nnlist(j,2)
        jmy = nnlist(j,4)
        do i = 1, lq
            iax = nnlist(i,1)
            imx = nnlist(i,3)
            iay = nnlist(i,2)
            imy = nnlist(i,4)
            Cotss = Cotss + grdnc(i,j)*grupc(i,j) + dconjg(grdnc(i,j)*grupc(i,j))
            Cost0s = Cost0s + dconjg(grupc(i,j))*grupc(i,j) + dconjg(grdnc(i,j))*grdnc(i,j)
            Cost1s = Cost1s + dconjg(grdnc(i,j))*grupc(i,j) + dconjg(grupc(i,j))*grdnc(i,j)
          Cnnt1ss = Cnnt1ss + grdnc(iax,jax)*grupc(i,j) + dconjg(grupc(iax,jax)*grdnc(i,j)) + &
                              grdnc(iax,jmx)*grupc(i,j) + dconjg(grupc(iax,jmx)*grdnc(i,j)) + &
                              grdnc(iax,jay)*grupc(i,j) + dconjg(grupc(iax,jay)*grdnc(i,j)) + &
                              grdnc(iax,jmy)*grupc(i,j) + dconjg(grupc(iax,jmy)*grdnc(i,j)) + &
                              grdnc(imx,jax)*grupc(i,j) + dconjg(grupc(imx,jax)*grdnc(i,j)) + &
                              grdnc(imx,jmx)*grupc(i,j) + dconjg(grupc(imx,jmx)*grdnc(i,j)) + &
                              grdnc(imx,jay)*grupc(i,j) + dconjg(grupc(imx,jay)*grdnc(i,j)) + &
                              grdnc(imx,jmy)*grupc(i,j) + dconjg(grupc(imx,jmy)*grdnc(i,j)) + &
                              grdnc(iay,jax)*grupc(i,j) + dconjg(grupc(iay,jax)*grdnc(i,j)) + &
                              grdnc(iay,jmx)*grupc(i,j) + dconjg(grupc(iay,jmx)*grdnc(i,j)) + &
                              grdnc(iay,jay)*grupc(i,j) + dconjg(grupc(iay,jay)*grdnc(i,j)) + &
                              grdnc(iay,jmy)*grupc(i,j) + dconjg(grupc(iay,jmy)*grdnc(i,j)) + &
                              grdnc(imy,jax)*grupc(i,j) + dconjg(grupc(imy,jax)*grdnc(i,j)) + &
                              grdnc(imy,jmx)*grupc(i,j) + dconjg(grupc(imy,jmx)*grdnc(i,j)) + &
                              grdnc(imy,jay)*grupc(i,j) + dconjg(grupc(imy,jay)*grdnc(i,j)) + &
                              grdnc(imy,jmy)*grupc(i,j) + dconjg(grupc(imy,jmy)*grdnc(i,j))
          Cnnt1sd = Cnnt1sd + (grdnc(iax,jax)*grupc(i,j) + dconjg(grupc(iax,jax)*grdnc(i,j)))  &
                            + (grdnc(iax,jmx)*grupc(i,j) + dconjg(grupc(iax,jmx)*grdnc(i,j)))  &
                            - (grdnc(iax,jay)*grupc(i,j) + dconjg(grupc(iax,jay)*grdnc(i,j)))  &
                            - (grdnc(iax,jmy)*grupc(i,j) + dconjg(grupc(iax,jmy)*grdnc(i,j)))  &
                            + (grdnc(imx,jax)*grupc(i,j) + dconjg(grupc(imx,jax)*grdnc(i,j)))  &
                            + (grdnc(imx,jmx)*grupc(i,j) + dconjg(grupc(imx,jmx)*grdnc(i,j)))  &
                            - (grdnc(imx,jay)*grupc(i,j) + dconjg(grupc(imx,jay)*grdnc(i,j)))  &
                            - (grdnc(imx,jmy)*grupc(i,j) + dconjg(grupc(imx,jmy)*grdnc(i,j)))  &
                            - (grdnc(iay,jax)*grupc(i,j) + dconjg(grupc(iay,jax)*grdnc(i,j)))  &
                            - (grdnc(iay,jmx)*grupc(i,j) + dconjg(grupc(iay,jmx)*grdnc(i,j)))  &
                            + (grdnc(iay,jay)*grupc(i,j) + dconjg(grupc(iay,jay)*grdnc(i,j)))  &
                            + (grdnc(iay,jmy)*grupc(i,j) + dconjg(grupc(iay,jmy)*grdnc(i,j)))  &
                            - (grdnc(imy,jax)*grupc(i,j) + dconjg(grupc(imy,jax)*grdnc(i,j)))  &
                            - (grdnc(imy,jmx)*grupc(i,j) + dconjg(grupc(imy,jmx)*grdnc(i,j)))  &
                            + (grdnc(imy,jay)*grupc(i,j) + dconjg(grupc(imy,jay)*grdnc(i,j)))  &
                            + (grdnc(imy,jmy)*grupc(i,j) + dconjg(grupc(imy,jmy)*grdnc(i,j)))
          Cnnt0ss = Cnnt0ss + dconjg(grupc(iax,jax))*grupc(i,j) + grdnc(iax,jax)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(iax,jmx))*grupc(i,j) + grdnc(iax,jmx)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(iax,jay))*grupc(i,j) + grdnc(iax,jay)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(iax,jmy))*grupc(i,j) + grdnc(iax,jmy)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imx,jax))*grupc(i,j) + grdnc(imx,jax)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imx,jmx))*grupc(i,j) + grdnc(imx,jmx)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imx,jay))*grupc(i,j) + grdnc(imx,jay)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imx,jmy))*grupc(i,j) + grdnc(imx,jmy)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(iay,jax))*grupc(i,j) + grdnc(iay,jax)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(iay,jmx))*grupc(i,j) + grdnc(iay,jmx)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(iay,jay))*grupc(i,j) + grdnc(iay,jay)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(iay,jmy))*grupc(i,j) + grdnc(iay,jmy)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imy,jax))*grupc(i,j) + grdnc(imy,jax)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imy,jmx))*grupc(i,j) + grdnc(imy,jmx)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imy,jay))*grupc(i,j) + grdnc(imy,jay)*dconjg(grdnc(i,j)) + &
                              dconjg(grupc(imy,jmy))*grupc(i,j) + grdnc(imy,jmy)*dconjg(grdnc(i,j))
          Cnnt0sd = Cnnt0sd + (dconjg(grupc(iax,jax))*grupc(i,j) + grdnc(iax,jax)*dconjg(grdnc(i,j))) &
                            + (dconjg(grupc(iax,jmx))*grupc(i,j) + grdnc(iax,jmx)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(iax,jay))*grupc(i,j) + grdnc(iax,jay)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(iax,jmy))*grupc(i,j) + grdnc(iax,jmy)*dconjg(grdnc(i,j))) &
                            + (dconjg(grupc(imx,jax))*grupc(i,j) + grdnc(imx,jax)*dconjg(grdnc(i,j))) &
                            + (dconjg(grupc(imx,jmx))*grupc(i,j) + grdnc(imx,jmx)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(imx,jay))*grupc(i,j) + grdnc(imx,jay)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(imx,jmy))*grupc(i,j) + grdnc(imx,jmy)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(iay,jax))*grupc(i,j) + grdnc(iay,jax)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(iay,jmx))*grupc(i,j) + grdnc(iay,jmx)*dconjg(grdnc(i,j))) &
                            + (dconjg(grupc(iay,jay))*grupc(i,j) + grdnc(iay,jay)*dconjg(grdnc(i,j))) &
                            + (dconjg(grupc(iay,jmy))*grupc(i,j) + grdnc(iay,jmy)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(imy,jax))*grupc(i,j) + grdnc(imy,jax)*dconjg(grdnc(i,j))) &
                            - (dconjg(grupc(imy,jmx))*grupc(i,j) + grdnc(imy,jmx)*dconjg(grdnc(i,j))) &
                            + (dconjg(grupc(imy,jay))*grupc(i,j) + grdnc(imy,jay)*dconjg(grdnc(i,j))) &
                            + (dconjg(grupc(imy,jmy))*grupc(i,j) + grdnc(imy,jmy)*dconjg(grdnc(i,j)))
            Cnnst0s = Cnnst0s + dconjg(grupc(iax,jax))*grupc(i,j) + dconjg(grdnc(iax,jax))*grdnc(i,j) + &
                                dconjg(grupc(iax,jmx))*grupc(i,j) + dconjg(grdnc(iax,jmx))*grdnc(i,j) + &
                                dconjg(grupc(iax,jay))*grupc(i,j) + dconjg(grdnc(iax,jay))*grdnc(i,j) + &
                                dconjg(grupc(iax,jmy))*grupc(i,j) + dconjg(grdnc(iax,jmy))*grdnc(i,j) + &
                                dconjg(grupc(imx,jax))*grupc(i,j) + dconjg(grdnc(imx,jax))*grdnc(i,j) + &
                                dconjg(grupc(imx,jmx))*grupc(i,j) + dconjg(grdnc(imx,jmx))*grdnc(i,j) + &
                                dconjg(grupc(imx,jay))*grupc(i,j) + dconjg(grdnc(imx,jay))*grdnc(i,j) + &
                                dconjg(grupc(imx,jmy))*grupc(i,j) + dconjg(grdnc(imx,jmy))*grdnc(i,j) + &
                                dconjg(grupc(iay,jax))*grupc(i,j) + dconjg(grdnc(iay,jax))*grdnc(i,j) + &
                                dconjg(grupc(iay,jmx))*grupc(i,j) + dconjg(grdnc(iay,jmx))*grdnc(i,j) + &
                                dconjg(grupc(iay,jay))*grupc(i,j) + dconjg(grdnc(iay,jay))*grdnc(i,j) + &
                                dconjg(grupc(iay,jmy))*grupc(i,j) + dconjg(grdnc(iay,jmy))*grdnc(i,j) + &
                                dconjg(grupc(imy,jax))*grupc(i,j) + dconjg(grdnc(imy,jax))*grdnc(i,j) + &
                                dconjg(grupc(imy,jmx))*grupc(i,j) + dconjg(grdnc(imy,jmx))*grdnc(i,j) + &
                                dconjg(grupc(imy,jay))*grupc(i,j) + dconjg(grdnc(imy,jay))*grdnc(i,j) + &
                                dconjg(grupc(imy,jmy))*grupc(i,j) + dconjg(grdnc(imy,jmy))*grdnc(i,j)
            Cnnst0d = Cnnst0d + (dconjg(grupc(iax,jax))*grupc(i,j) + dconjg(grdnc(iax,jax))*grdnc(i,j))  &
                              + (dconjg(grupc(iax,jmx))*grupc(i,j) + dconjg(grdnc(iax,jmx))*grdnc(i,j))  &
                              - (dconjg(grupc(iax,jay))*grupc(i,j) + dconjg(grdnc(iax,jay))*grdnc(i,j))  &
                              - (dconjg(grupc(iax,jmy))*grupc(i,j) + dconjg(grdnc(iax,jmy))*grdnc(i,j))  &
                              + (dconjg(grupc(imx,jax))*grupc(i,j) + dconjg(grdnc(imx,jax))*grdnc(i,j))  &
                              + (dconjg(grupc(imx,jmx))*grupc(i,j) + dconjg(grdnc(imx,jmx))*grdnc(i,j))  &
                              - (dconjg(grupc(imx,jay))*grupc(i,j) + dconjg(grdnc(imx,jay))*grdnc(i,j))  &
                              - (dconjg(grupc(imx,jmy))*grupc(i,j) + dconjg(grdnc(imx,jmy))*grdnc(i,j))  &
                              - (dconjg(grupc(iay,jax))*grupc(i,j) + dconjg(grdnc(iay,jax))*grdnc(i,j))  &
                              - (dconjg(grupc(iay,jmx))*grupc(i,j) + dconjg(grdnc(iay,jmx))*grdnc(i,j))  &
                              + (dconjg(grupc(iay,jay))*grupc(i,j) + dconjg(grdnc(iay,jay))*grdnc(i,j))  &
                              + (dconjg(grupc(iay,jmy))*grupc(i,j) + dconjg(grdnc(iay,jmy))*grdnc(i,j))  &
                              - (dconjg(grupc(imy,jax))*grupc(i,j) + dconjg(grdnc(imy,jax))*grdnc(i,j))  &
                              - (dconjg(grupc(imy,jmx))*grupc(i,j) + dconjg(grdnc(imy,jmx))*grdnc(i,j))  &
                              + (dconjg(grupc(imy,jay))*grupc(i,j) + dconjg(grdnc(imy,jay))*grdnc(i,j))  &
                              + (dconjg(grupc(imy,jmy))*grupc(i,j) + dconjg(grdnc(imy,jmy))*grdnc(i,j))
            Cnnst1s = Cnnst1s + dconjg(grdnc(iax,jax))*grupc(i,j) + dconjg(grupc(iax,jax))*grdnc(i,j) + &
                                dconjg(grdnc(iax,jmx))*grupc(i,j) + dconjg(grupc(iax,jmx))*grdnc(i,j) + &
                                dconjg(grdnc(iax,jay))*grupc(i,j) + dconjg(grupc(iax,jay))*grdnc(i,j) + &
                                dconjg(grdnc(iax,jmy))*grupc(i,j) + dconjg(grupc(iax,jmy))*grdnc(i,j) + &
                                dconjg(grdnc(imx,jax))*grupc(i,j) + dconjg(grupc(imx,jax))*grdnc(i,j) + &
                                dconjg(grdnc(imx,jmx))*grupc(i,j) + dconjg(grupc(imx,jmx))*grdnc(i,j) + &
                                dconjg(grdnc(imx,jay))*grupc(i,j) + dconjg(grupc(imx,jay))*grdnc(i,j) + &
                                dconjg(grdnc(imx,jmy))*grupc(i,j) + dconjg(grupc(imx,jmy))*grdnc(i,j) + &
                                dconjg(grdnc(iay,jax))*grupc(i,j) + dconjg(grupc(iay,jax))*grdnc(i,j) + &
                                dconjg(grdnc(iay,jmx))*grupc(i,j) + dconjg(grupc(iay,jmx))*grdnc(i,j) + &
                                dconjg(grdnc(iay,jay))*grupc(i,j) + dconjg(grupc(iay,jay))*grdnc(i,j) + &
                                dconjg(grdnc(iay,jmy))*grupc(i,j) + dconjg(grupc(iay,jmy))*grdnc(i,j) + &
                                dconjg(grdnc(imy,jax))*grupc(i,j) + dconjg(grupc(imy,jax))*grdnc(i,j) + &
                                dconjg(grdnc(imy,jmx))*grupc(i,j) + dconjg(grupc(imy,jmx))*grdnc(i,j) + &
                                dconjg(grdnc(imy,jay))*grupc(i,j) + dconjg(grupc(imy,jay))*grdnc(i,j) + &
                                dconjg(grdnc(imy,jmy))*grupc(i,j) + dconjg(grupc(imy,jmy))*grdnc(i,j)
            Cnnst1d = Cnnst1d + (dconjg(grdnc(iax,jax))*grupc(i,j) + dconjg(grupc(iax,jax))*grdnc(i,j))  &
                              + (dconjg(grdnc(iax,jmx))*grupc(i,j) + dconjg(grupc(iax,jmx))*grdnc(i,j))  &
                              - (dconjg(grdnc(iax,jay))*grupc(i,j) + dconjg(grupc(iax,jay))*grdnc(i,j))  &
                              - (dconjg(grdnc(iax,jmy))*grupc(i,j) + dconjg(grupc(iax,jmy))*grdnc(i,j))  &
                              + (dconjg(grdnc(imx,jax))*grupc(i,j) + dconjg(grupc(imx,jax))*grdnc(i,j))  &
                              + (dconjg(grdnc(imx,jmx))*grupc(i,j) + dconjg(grupc(imx,jmx))*grdnc(i,j))  &
                              - (dconjg(grdnc(imx,jay))*grupc(i,j) + dconjg(grupc(imx,jay))*grdnc(i,j))  &
                              - (dconjg(grdnc(imx,jmy))*grupc(i,j) + dconjg(grupc(imx,jmy))*grdnc(i,j))  &
                              - (dconjg(grdnc(iay,jax))*grupc(i,j) + dconjg(grupc(iay,jax))*grdnc(i,j))  &
                              - (dconjg(grdnc(iay,jmx))*grupc(i,j) + dconjg(grupc(iay,jmx))*grdnc(i,j))  &
                              + (dconjg(grdnc(iay,jay))*grupc(i,j) + dconjg(grupc(iay,jay))*grdnc(i,j))  &
                              + (dconjg(grdnc(iay,jmy))*grupc(i,j) + dconjg(grupc(iay,jmy))*grdnc(i,j))  &
                              - (dconjg(grdnc(imy,jax))*grupc(i,j) + dconjg(grupc(imy,jax))*grdnc(i,j))  &
                              - (dconjg(grdnc(imy,jmx))*grupc(i,j) + dconjg(grupc(imy,jmx))*grdnc(i,j))  &
                              + (dconjg(grdnc(imy,jay))*grupc(i,j) + dconjg(grupc(imy,jay))*grdnc(i,j))  &
                              + (dconjg(grdnc(imy,jmy))*grupc(i,j) + dconjg(grupc(imy,jmy))*grdnc(i,j))
          Cnnt1t0px = Cnnt1t0px + (grdnc(iax,jax)*grupc(i,j) + dconjg(grupc(iax,jax)*grdnc(i,j))) &
                                - (grdnc(iax,jmx)*grupc(i,j) + dconjg(grupc(iax,jmx)*grdnc(i,j))) &
                                - (grdnc(imx,jax)*grupc(i,j) + dconjg(grupc(imx,jax)*grdnc(i,j))) &
                                + (grdnc(imx,jmx)*grupc(i,j) + dconjg(grupc(imx,jmx)*grdnc(i,j)))
    Cnnt1t0pxipy = Cnnt1t0pxipy +     (grdnc(iax,jax)*grupc(i,j) + dconjg(grupc(iax,jax)*grdnc(i,j))) &
                                -     (grdnc(iax,jmx)*grupc(i,j) + dconjg(grupc(iax,jmx)*grdnc(i,j))) &
                                + czi*(grdnc(iax,jay)*grupc(i,j) + dconjg(grupc(iax,jay)*grdnc(i,j))) &
                                - czi*(grdnc(iax,jmy)*grupc(i,j) + dconjg(grupc(iax,jmy)*grdnc(i,j))) &
                                -     (grdnc(imx,jax)*grupc(i,j) + dconjg(grupc(imx,jax)*grdnc(i,j))) &
                                +     (grdnc(imx,jmx)*grupc(i,j) + dconjg(grupc(imx,jmx)*grdnc(i,j))) &
                                - czi*(grdnc(imx,jay)*grupc(i,j) + dconjg(grupc(imx,jay)*grdnc(i,j))) &
                                + czi*(grdnc(imx,jmy)*grupc(i,j) + dconjg(grupc(imx,jmy)*grdnc(i,j))) &
                                - czi*(grdnc(iay,jax)*grupc(i,j) + dconjg(grupc(iay,jax)*grdnc(i,j))) &
                                + czi*(grdnc(iay,jmx)*grupc(i,j) + dconjg(grupc(iay,jmx)*grdnc(i,j))) &
                                +     (grdnc(iay,jay)*grupc(i,j) + dconjg(grupc(iay,jay)*grdnc(i,j))) &
                                -     (grdnc(iay,jmy)*grupc(i,j) + dconjg(grupc(iay,jmy)*grdnc(i,j))) &
                                + czi*(grdnc(imy,jax)*grupc(i,j) + dconjg(grupc(imy,jax)*grdnc(i,j))) &
                                - czi*(grdnc(imy,jmx)*grupc(i,j) + dconjg(grupc(imy,jmx)*grdnc(i,j))) &
                                -     (grdnc(imy,jay)*grupc(i,j) + dconjg(grupc(imy,jay)*grdnc(i,j))) &
                                +     (grdnc(imy,jmy)*grupc(i,j) + dconjg(grupc(imy,jmy)*grdnc(i,j)))
          Cnnt0t0px = Cnnt0t0px + (dconjg(grupc(iax,jax))*grupc(i,j) + grdnc(iax,jax)*dconjg(grdnc(i,j))) &
                                - (dconjg(grupc(iax,jmx))*grupc(i,j) + grdnc(iax,jmx)*dconjg(grdnc(i,j))) &
                                - (dconjg(grupc(imx,jax))*grupc(i,j) + grdnc(imx,jax)*dconjg(grdnc(i,j))) &
                                + (dconjg(grupc(imx,jmx))*grupc(i,j) + grdnc(imx,jmx)*dconjg(grdnc(i,j)))
    Cnnt0t0pxipy = Cnnt0t0pxipy +     (dconjg(grupc(iax,jax))*grupc(i,j) + grdnc(iax,jax)*dconjg(grdnc(i,j))) &
                                -     (dconjg(grupc(iax,jmx))*grupc(i,j) + grdnc(iax,jmx)*dconjg(grdnc(i,j))) &
                                + czi*(dconjg(grupc(iax,jay))*grupc(i,j) + grdnc(iax,jay)*dconjg(grdnc(i,j))) &
                                - czi*(dconjg(grupc(iax,jmy))*grupc(i,j) + grdnc(iax,jmy)*dconjg(grdnc(i,j))) &
                                -     (dconjg(grupc(imx,jax))*grupc(i,j) + grdnc(imx,jax)*dconjg(grdnc(i,j))) &
                                +     (dconjg(grupc(imx,jmx))*grupc(i,j) + grdnc(imx,jmx)*dconjg(grdnc(i,j))) &
                                - czi*(dconjg(grupc(imx,jay))*grupc(i,j) + grdnc(imx,jay)*dconjg(grdnc(i,j))) &
                                + czi*(dconjg(grupc(imx,jmy))*grupc(i,j) + grdnc(imx,jmy)*dconjg(grdnc(i,j))) &
                                - czi*(dconjg(grupc(iay,jax))*grupc(i,j) + grdnc(iay,jax)*dconjg(grdnc(i,j))) &
                                + czi*(dconjg(grupc(iay,jmx))*grupc(i,j) + grdnc(iay,jmx)*dconjg(grdnc(i,j))) &
                                +     (dconjg(grupc(iay,jay))*grupc(i,j) + grdnc(iay,jay)*dconjg(grdnc(i,j))) &
                                -     (dconjg(grupc(iay,jmy))*grupc(i,j) + grdnc(iay,jmy)*dconjg(grdnc(i,j))) &
                                + czi*(dconjg(grupc(imy,jax))*grupc(i,j) + grdnc(imy,jax)*dconjg(grdnc(i,j))) &
                                - czi*(dconjg(grupc(imy,jmx))*grupc(i,j) + grdnc(imy,jmx)*dconjg(grdnc(i,j))) &
                                -     (dconjg(grupc(imy,jay))*grupc(i,j) + grdnc(imy,jay)*dconjg(grdnc(i,j))) &
                                +     (dconjg(grupc(imy,jmy))*grupc(i,j) + grdnc(imy,jmy)*dconjg(grdnc(i,j)))
    Cnnt1t1px = Cnnt1t1px + dcmplx(2.d0*dble( grupc(iax,jax)*grupc(i,j) - grupc(iax,j)*grupc(i,jax) + grdnc(iax,jax)*grdnc(i,j) - grdnc(iax,j)*grdnc(i,jax) ),0.d0) &
                          - dcmplx(2.d0*dble( grupc(iax,jmx)*grupc(i,j) - grupc(iax,j)*grupc(i,jmx) + grdnc(iax,jmx)*grdnc(i,j) - grdnc(iax,j)*grdnc(i,jmx) ),0.d0) &
                          - dcmplx(2.d0*dble( grupc(imx,jax)*grupc(i,j) - grupc(imx,j)*grupc(i,jax) + grdnc(imx,jax)*grdnc(i,j) - grdnc(imx,j)*grdnc(i,jax) ),0.d0) &
                          + dcmplx(2.d0*dble( grupc(imx,jmx)*grupc(i,j) - grupc(imx,j)*grupc(i,jmx) + grdnc(imx,jmx)*grdnc(i,j) - grdnc(imx,j)*grdnc(i,jmx) ),0.d0)
Cnnt1t1pxipy = Cnnt1t1pxipy +     dcmplx(2.d0*dble( grupc(iax,jax)*grupc(i,j) - grupc(iax,j)*grupc(i,jax) + grdnc(iax,jax)*grdnc(i,j) - grdnc(iax,j)*grdnc(i,jax) ),0.d0) &
                            -     dcmplx(2.d0*dble( grupc(iax,jmx)*grupc(i,j) - grupc(iax,j)*grupc(i,jmx) + grdnc(iax,jmx)*grdnc(i,j) - grdnc(iax,j)*grdnc(i,jmx) ),0.d0) &
                            + czi*dcmplx(2.d0*dble( grupc(iax,jay)*grupc(i,j) - grupc(iax,j)*grupc(i,jay) + grdnc(iax,jay)*grdnc(i,j) - grdnc(iax,j)*grdnc(i,jay) ),0.d0) &
                            - czi*dcmplx(2.d0*dble( grupc(iax,jmy)*grupc(i,j) - grupc(iax,j)*grupc(i,jmy) + grdnc(iax,jmy)*grdnc(i,j) - grdnc(iax,j)*grdnc(i,jmy) ),0.d0) &
                            -     dcmplx(2.d0*dble( grupc(imx,jax)*grupc(i,j) - grupc(imx,j)*grupc(i,jax) + grdnc(imx,jax)*grdnc(i,j) - grdnc(imx,j)*grdnc(i,jax) ),0.d0) &
                            +     dcmplx(2.d0*dble( grupc(imx,jmx)*grupc(i,j) - grupc(imx,j)*grupc(i,jmx) + grdnc(imx,jmx)*grdnc(i,j) - grdnc(imx,j)*grdnc(i,jmx) ),0.d0) &
                            - czi*dcmplx(2.d0*dble( grupc(imx,jay)*grupc(i,j) - grupc(imx,j)*grupc(i,jay) + grdnc(imx,jay)*grdnc(i,j) - grdnc(imx,j)*grdnc(i,jay) ),0.d0) &
                            + czi*dcmplx(2.d0*dble( grupc(imx,jmy)*grupc(i,j) - grupc(imx,j)*grupc(i,jmy) + grdnc(imx,jmy)*grdnc(i,j) - grdnc(imx,j)*grdnc(i,jmy) ),0.d0) &
                            - czi*dcmplx(2.d0*dble( grupc(iay,jax)*grupc(i,j) - grupc(iay,j)*grupc(i,jax) + grdnc(iay,jax)*grdnc(i,j) - grdnc(iay,j)*grdnc(i,jax) ),0.d0) &
                            + czi*dcmplx(2.d0*dble( grupc(iay,jmx)*grupc(i,j) - grupc(iay,j)*grupc(i,jmx) + grdnc(iay,jmx)*grdnc(i,j) - grdnc(iay,j)*grdnc(i,jmx) ),0.d0) &
                            +     dcmplx(2.d0*dble( grupc(iay,jay)*grupc(i,j) - grupc(iay,j)*grupc(i,jay) + grdnc(iay,jay)*grdnc(i,j) - grdnc(iay,j)*grdnc(i,jay) ),0.d0) &
                            -     dcmplx(2.d0*dble( grupc(iay,jmy)*grupc(i,j) - grupc(iay,j)*grupc(i,jmy) + grdnc(iay,jmy)*grdnc(i,j) - grdnc(iay,j)*grdnc(i,jmy) ),0.d0) &
                            + czi*dcmplx(2.d0*dble( grupc(imy,jax)*grupc(i,j) - grupc(imy,j)*grupc(i,jax) + grdnc(imy,jax)*grdnc(i,j) - grdnc(imy,j)*grdnc(i,jax) ),0.d0) &
                            - czi*dcmplx(2.d0*dble( grupc(imy,jmx)*grupc(i,j) - grupc(imy,j)*grupc(i,jmx) + grdnc(imy,jmx)*grdnc(i,j) - grdnc(imy,j)*grdnc(i,jmx) ),0.d0) &
                            -     dcmplx(2.d0*dble( grupc(imy,jay)*grupc(i,j) - grupc(imy,j)*grupc(i,jay) + grdnc(imy,jay)*grdnc(i,j) - grdnc(imy,j)*grdnc(i,jay) ),0.d0) &
                            +     dcmplx(2.d0*dble( grupc(imy,jmy)*grupc(i,j) - grupc(imy,j)*grupc(i,jmy) + grdnc(imy,jmy)*grdnc(i,j) - grdnc(imy,j)*grdnc(i,jmy) ),0.d0)
          Cnnt0t1px = Cnnt0t1px + (dconjg(grdnc(iax,jax))*grupc(i,j) + dconjg(grupc(iax,jax))*grdnc(i,j)) &
                                - (dconjg(grdnc(iax,jmx))*grupc(i,j) + dconjg(grupc(iax,jmx))*grdnc(i,j)) &
                                - (dconjg(grdnc(imx,jax))*grupc(i,j) + dconjg(grupc(imx,jax))*grdnc(i,j)) &
                                + (dconjg(grdnc(imx,jmx))*grupc(i,j) + dconjg(grupc(imx,jmx))*grdnc(i,j))
    Cnnt0t1pxipy = Cnnt0t1pxipy +     (dconjg(grdnc(iax,jax))*grupc(i,j) + dconjg(grupc(iax,jax))*grdnc(i,j)) &
                                -     (dconjg(grdnc(iax,jmx))*grupc(i,j) + dconjg(grupc(iax,jmx))*grdnc(i,j)) &
                                + czi*(dconjg(grdnc(iax,jay))*grupc(i,j) + dconjg(grupc(iax,jay))*grdnc(i,j)) &
                                - czi*(dconjg(grdnc(iax,jmy))*grupc(i,j) + dconjg(grupc(iax,jmy))*grdnc(i,j)) &
                                -     (dconjg(grdnc(imx,jax))*grupc(i,j) + dconjg(grupc(imx,jax))*grdnc(i,j)) &
                                +     (dconjg(grdnc(imx,jmx))*grupc(i,j) + dconjg(grupc(imx,jmx))*grdnc(i,j)) &
                                - czi*(dconjg(grdnc(imx,jay))*grupc(i,j) + dconjg(grupc(imx,jay))*grdnc(i,j)) &
                                + czi*(dconjg(grdnc(imx,jmy))*grupc(i,j) + dconjg(grupc(imx,jmy))*grdnc(i,j)) &
                                - czi*(dconjg(grdnc(iay,jax))*grupc(i,j) + dconjg(grupc(iay,jax))*grdnc(i,j)) &
                                + czi*(dconjg(grdnc(iay,jmx))*grupc(i,j) + dconjg(grupc(iay,jmx))*grdnc(i,j)) &
                                +     (dconjg(grdnc(iay,jay))*grupc(i,j) + dconjg(grupc(iay,jay))*grdnc(i,j)) &
                                -     (dconjg(grdnc(iay,jmy))*grupc(i,j) + dconjg(grupc(iay,jmy))*grdnc(i,j)) &
                                + czi*(dconjg(grdnc(imy,jax))*grupc(i,j) + dconjg(grupc(imy,jax))*grdnc(i,j)) &
                                - czi*(dconjg(grdnc(imy,jmx))*grupc(i,j) + dconjg(grupc(imy,jmx))*grdnc(i,j)) &
                                -     (dconjg(grdnc(imy,jay))*grupc(i,j) + dconjg(grupc(imy,jay))*grdnc(i,j)) &
                                +     (dconjg(grdnc(imy,jmy))*grupc(i,j) + dconjg(grupc(imy,jmy))*grdnc(i,j))
        end do
    end do
    pair_bin( 1) = pair_bin( 1) + Cotss         / dcmplx( 2.d0, 0.d0 ) 
    pair_bin( 2) = pair_bin( 2) + Cost0s        / dcmplx( 2.d0, 0.d0 ) 
    pair_bin( 3) = pair_bin( 3) + Cost1s        / dcmplx( 2.d0, 0.d0 ) 

    pair_bin( 4) = pair_bin( 4) + Cnnt1ss       / dcmplx( 4.d0, 0.d0 ) 
    pair_bin( 5) = pair_bin( 5) + Cnnt1sd       / dcmplx( 4.d0, 0.d0 ) 
    pair_bin( 6) = pair_bin( 6) + Cnnt0ss       / dcmplx( 4.d0, 0.d0 ) 
    pair_bin( 7) = pair_bin( 7) + Cnnt0sd       / dcmplx( 4.d0, 0.d0 ) 
    pair_bin( 8) = pair_bin( 8) + Cnnst0s       / dcmplx( 4.d0, 0.d0 ) 
    pair_bin( 9) = pair_bin( 9) + Cnnst0d       / dcmplx( 4.d0, 0.d0 ) 
    pair_bin(10) = pair_bin(10) + Cnnst1s       / dcmplx( 4.d0, 0.d0 ) 
    pair_bin(11) = pair_bin(11) + Cnnst1d       / dcmplx( 4.d0, 0.d0 ) 

    pair_bin(12) = pair_bin(12) + Cnnt1t0px     / dcmplx( 2.d0, 0.d0 ) 
    pair_bin(13) = pair_bin(13) + Cnnt1t0pxipy  / dcmplx( 4.d0, 0.d0 ) 
    pair_bin(14) = pair_bin(14) + Cnnt0t0px     / dcmplx( 2.d0, 0.d0 ) 
    pair_bin(15) = pair_bin(15) + Cnnt0t0pxipy  / dcmplx( 4.d0, 0.d0 ) 
    pair_bin(16) = pair_bin(16) + Cnnt1t1px     / dcmplx( 8.d0, 0.d0 ) 
    pair_bin(17) = pair_bin(17) + Cnnt1t1pxipy  / dcmplx(16.d0, 0.d0 ) 
    pair_bin(18) = pair_bin(18) + Cnnt0t1px     / dcmplx( 2.d0, 0.d0 ) 
    pair_bin(19) = pair_bin(19) + Cnnt0t1pxipy  / dcmplx( 4.d0, 0.d0 ) 

    ! high order pairing
    do j = 1, lq
        do i = 1, lq

            ! charge 4e
            high_pair_bin(1) = high_pair_bin(1) + grdnc(i,j)*grupc(i,j) * dconjg(grdnc(i,j)*grupc(i,j))

            ! spin nematic
            high_pair_bin(2) = high_pair_bin(2) + grdn(i,j)*grupc(i,j) * dconjg(grdn(i,j)*grupc(i,j))

            ! triplet a or b
            high_pair_bin(3) = high_pair_bin(3) + dconjg(grupc(i,j))*grupc(i,j)

            ! triplet a or b
            high_pair_bin(4) = high_pair_bin(4) + grdnc(i,j)*dconjg(grdnc(i,j))
        end do
    end do

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
!$OMP DO REDUCTION ( + : gtau_up, gtau_dn, chiszsz, chijxjx, chijxjxaa, chijxjxab, chijxjxba, chijxjxbb)
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

                ztmp = ( gr0t_up(jax,i) - gr0t_up(jmx,i) ) * ( grt0_up(iax,j) - grt0_up(imx,j) ) - &
                       ( grtt_up(iax,i) - grtt_up(imx,i) ) * ( gr00_up(jax,j) - gr00_up(jmx,j) )
                chijxjxaa(imj,nt) = chijxjxaa(imj,nt) + ztmp + dconjg(ztmp)

                ztmp = ( gr0t_dn(jax,i) - gr0t_dn(jmx,i) ) * ( grt0_dn(iax,j) - grt0_dn(imx,j) ) - &
                       ( grtt_dn(iax,i) - grtt_dn(imx,i) ) * ( gr00_dn(jax,j) - gr00_dn(jmx,j) )
                chijxjxbb(imj,nt) = chijxjxbb(imj,nt) + ztmp + dconjg(ztmp)

                ztmp = - ( grtt_up(iax,i) - grtt_up(imx,i) ) * ( gr00_dn(jax,j) - gr00_dn(jmx,j) )
                chijxjxab(imj,nt) = chijxjxab(imj,nt) + ztmp + dconjg(ztmp)

                ztmp = - ( grtt_dn(iax,i) - grtt_dn(imx,i) ) * ( gr00_up(jax,j) - gr00_up(jmx,j) )
                chijxjxba(imj,nt) = chijxjxba(imj,nt) + ztmp + dconjg(ztmp)

            end do
        end do
!$OMP END DO
!$OMP END PARALLEL
  end subroutine obsert

end module obser
