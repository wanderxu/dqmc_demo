subroutine generate_neighbor(zmax,nnmax,nntable)
  implicit none
  integer, intent(in) :: zmax, nnmax
  integer, dimension(2,zmax,nnmax), intent(out) :: nntable

  ! local 
  integer :: i, ix, iy, nlent, nn, ni
  
  real(8), dimension(:,:), allocatable :: len_table
  real(8), dimension(:), allocatable :: len_nn

  nntable = 0

  nlent = (2*nnmax+1)**2
  allocate( len_table(3,nlent) )
  allocate( len_nn(nnmax) )

  i = 0
  do iy = -nnmax, nnmax
      do ix = -nnmax, nnmax
          i = i + 1
          len_table(1,i) = dsqrt( dble(iy*iy)+dble(ix*ix) )
          len_table(2,i) = dble(ix)
          len_table(3,i) = dble(iy)
      end do
  end do
  call s_heapsort(nlent,3,len_table)

  nn = 1
  ni = 1
  i = 2  ! 1 is itself, 2 is the start
  len_nn(nn) = len_table(1,i) 
  nntable(1,ni,nn) = nint( len_table(2,i) )
  nntable(2,ni,nn) = nint( len_table(3,i) )
  do i = 3, nlent
      if( len_table(1,i) .gt. len_nn(nn) ) then
          nn = nn + 1
          ni = 1
          if( nn .gt. nnmax ) exit
          len_nn(nn) = len_table(1,i)
          nntable(1,ni,nn) = nint( len_table(2,i) )
          nntable(2,ni,nn) = nint( len_table(3,i) )
      else
          ni = ni + 1
          if(ni>zmax) stop ' you should increase zmax '
          nntable(1,ni,nn) = nint( len_table(2,i) )
          nntable(2,ni,nn) = nint( len_table(3,i) )
      end if
  end do

  deallocate( len_nn, len_table )

end subroutine
