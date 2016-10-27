subroutine sli
  use blockc

  implicit none

  ! local
  integer :: ncount, nx, ny, n, i, j, iq, ix, iy, nk, imj_nx, imj_ny, imj
  real(dp) :: ri(2), qvec(2), rr
  real(dp) :: fd(4)
  integer, external :: npbc

  real(dp), dimension(:,:), allocatable :: rlist

  list = 0; invlist = 0; nnlist = 0

  ncount = 0
  do nx = 1,l
  do ny = 1,l
     ncount = ncount + 1
     list(ncount,1) = nx
     list(ncount,2) = ny
     invlist(nx,ny) = ncount
  enddo
  enddo

  do i = 1, ndim
      orblist(i) = (i-1)/lq + 1
  end do

  do n = 1,lq
     nx = list(n,1)
     ny = list(n,2)
     nnlist(n,0) = invlist( nx , ny )
     nnlist(n,1) = invlist( npbc(nx+1,l) , ny )
     nnlist(n,2) = invlist( nx , npbc(ny+1,l) )
     nnlist(n,3) = invlist( npbc(nx-1,l) , ny )
     nnlist(n,4) = invlist( nx , npbc(ny-1,l) )
     nnlist(n,5) = invlist( npbc(nx+1,l) , npbc(ny+1,l) )
     nnlist(n,6) = invlist( npbc(nx-1,l) , npbc(ny+1,l) )
     nnlist(n,7) = invlist( npbc(nx-1,l) , npbc(ny-1,l) )
     nnlist(n,8) = invlist( npbc(nx+1,l) , npbc(ny-1,l) )
  enddo
  
  fd(1) =  1.d0
  fd(2) = -1.d0
  fd(3) =  1.d0
  fd(4) = -1.d0


  !zkron = dcmplx(0.d0,0.d0)
  !do i = 1,lq
  !   zkron(i,i) = dcmplx(1.d0,0.d0)
  !enddo
  
  list_plaq = 0
  do i = 1,lq
     ix = list(i,1)
     iy = list(i,2)
     list_plaq(i,1) = invlist(     ix      ,iy           ) 
     list_plaq(i,2) = invlist(npbc(ix+1,l) ,iy           ) 
     list_plaq(i,3) = invlist(npbc(ix+1,l) ,npbc(iy+1,l) )
     list_plaq(i,4) = invlist(     ix      ,npbc(iy+1,l) )
     list_plaq(i,5) = invlist(     ix      ,iy           ) 
  enddo

  ! latt_imj
  do j = 1, lq
      do i = 1, lq
          imj_nx = npbc( list(i,1) - list(j,1), l )
          imj_ny = npbc( list(i,2) - list(j,2), l )
          latt_imj(i,j) = invlist( imj_nx, imj_ny )
      end do
  end do

  ! latt_listk
  nk = 0
  do j = 0, l-1
      do i = 0, l-1
          nk = nk+1
          listk(nk,1) = i-l/2
          listk(nk,2) = j-l/2
      end do
  end do
  if( nk .ne. lq ) then
      stop " Error in set listk "
  end if

  ! set zexpiwt, zexpiqr
  do n = 0, nuse
      do i = 0, ltrot-1
          zexpiwt(i,n) = exp( dcmplx( 0.d0, wn(n)*dble(i)*dtau ) )
      end do
  end do

  do iq = 1, lq
      qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
      do i = 1, lq
          ri = dble(list(i,1))*a1_p + dble(list(i,2))*a2_p
          zexpiqr(i,iq) = exp( dcmplx( 0.d0,  qvec(1)*ri(1) + qvec(2)*ri(2) ) )
      end do
  end do

  ! the degerate of imj
  imjdeg(:) = 0
  do j = 1, lq
      do i = 1, lq
          imj = latt_imj(i,j)
          imjdeg(imj) = imjdeg(imj) + 1
      end do
  end do

  ! the distance list
  allocate( rlist(2,lq) )
  do i = 1, lq
      ix = list(i,1)
      iy = list(i,2)
      if( ix > (l+1)/2 ) ix = ix - l
      if( iy > (l+1)/2 ) iy = iy - l
      ri = dble(ix)*a1_p + dble(iy)*a2_p
      rr = dsqrt( ri(1)*ri(1) + ri(2)*ri(2) )
      rlist(2,i) = dble(i)
      rlist(1,i) = rr
      !write(fout, '(i4,f16.8)')i, rr
  end do
  call s_heapsort(lq,2,rlist)

  do i = 1, lq
      distance_index(i) = anint(rlist(2,i))
      distance_len(i) = rlist(1,i)
  end do

  irre_distance_len(:) = zero
  irre_distance_deg(:) = 0
  equ_distance(1) = 1
  rr = distance_len(1)
  irre_distance_len(1) = rr
  irre_distance_deg(1) = 1
  j = 1
  do i = 2, lq
      if (rr < distance_len(i) ) then
          rr = distance_len(i)
          j = j + 1
          equ_distance(i) = j
          irre_distance_len(j) = distance_len(i)
          irre_distance_deg(j) = 1
      else
          equ_distance(i) = j
          irre_distance_deg(j) = irre_distance_deg(j) + 1
      end if
  end do
  num_equ_distance = j

  deallocate( rlist )

end subroutine sli
