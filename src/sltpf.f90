subroutine sltpf

  use blockc
  implicit none
  
  ! local
  logical :: ltest
  integer :: nf, nc, ix, iy, i, i0, nc1, nc2
  
  nf = 1
  nc = 0
  do ix = 1,l	
  do iy = 1,l
     if (mod(ix,2).ne.0 ) then
        nc = nc + 1
        ltpf(nc,nf) = invlist(ix,iy)
     endif
  enddo
  enddo
  ! write(6,*) 'length fam1: ', nc
  
  nf = 2
  nc = 0
  do ix = 1,l	
  do iy = 1,l
     if (mod(ix,2).eq.0 ) then
        nc = nc + 1
        ltpf(nc,nf) = invlist(ix,iy)
     endif
  enddo
  enddo
  ! write(6,*) 'length fam2: ', nc
  
  nf = 3
  nc = 0
  do ix = 1,l	
  do iy = 1,l
     if (mod(iy,2).ne.0 ) then
        nc = nc + 1
        ltpf(nc,nf) = invlist(ix,iy)
     endif
  enddo
  enddo
  ! write(6,*) 'length fam3: ', nc
  
  nf = 4
  nc = 0
  do ix = 1,l	
  do iy = 1,l
     if (mod(iy,2).eq.0 ) then
        nc = nc + 1
        ltpf(nc,nf) = invlist(ix,iy)
     endif
  enddo
  enddo
  
  
  ltest = .false.
  !ltest = .true.
  if (ltest) then
     do nf = 1,4
     do i  = 1,lfam
        i0 = ltpf(i,nf)
        write(6,*) nf, list(i0,1), list(i0,2)
     enddo
     enddo
  endif
  
  ! for particle hopping.
  nc1 = 0
  nc2 = 0
  do i = 1,lq
     ix = list(i,1)
     iy = list(i,2)
     if (mod(ix,2).eq.0 .and. mod(iy,2).eq.0 ) then
        nc1 = nc1 + 1
        lthf(nc1,1) = i
     endif
     if (mod(ix,2).ne.0 .and. mod(iy,2).ne.0 ) then
        nc2 = nc2 + 1
        lthf(nc2,2) = i
     endif
  enddo
  if( l .gt. 1 ) then
      if (nc1.ne.lq/4 .or. nc2  .ne. lq/4 ) then
         write(6,*) 'error 1'
         stop
      endif
  end if
end subroutine sltpf
