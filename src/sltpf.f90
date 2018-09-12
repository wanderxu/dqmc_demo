subroutine sltpf

  use blockc
  implicit none
  
  ! local
  logical :: ltest
  integer :: nf, nc, ix, iy, i, i0, nc1, nc2, nc3, nc4
  integer :: nc5, nc6, nc7, nc8, nc9, nc10, nc11, nc12

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
  nc3 = 0
  nc4 = 0
  nc5 = 0
  nc6 = 0
  nc7 = 0
  nc8 = 0
  nc9 = 0
  nc10 = 0
  nc11 = 0
  nc12 = 0
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
#ifdef BREAKUP_T
     !3rd nearest
     !=====================================================!
     if (mod(ix-3,4).eq.0 .and. mod(iy-3,4).eq.0 ) then
        nc5 = nc5 + 1
        lthf3(nc5,1) = i
     endif
     if (mod(ix-1,4).eq.0 .and. mod(iy-1,4).eq.0 ) then
        nc6 = nc6 + 1
        lthf3(nc6,2) = i
     endif
     if (mod(ix-4,4).eq.0 .and. mod(iy-3,4).eq.0 ) then
        nc7 = nc7 + 1
        lthf3(nc7,3) = i
     endif
     if (mod(ix-2,4).eq.0 .and. mod(iy-1,4).eq.0 ) then
        nc8 = nc8 + 1
        lthf3(nc8,4) = i
     endif
     if (mod(ix-3,4).eq.0 .and. mod(iy-4,4).eq.0 ) then
        nc9 = nc9 + 1
        lthf3(nc9,5) = i
     endif
     if (mod(ix-1,4).eq.0 .and. mod(iy-2,4).eq.0 ) then
        nc10 = nc10 + 1
        lthf3(nc10,6) = i
     endif
     if (mod(ix-4,4).eq.0 .and. mod(iy-4,4).eq.0 ) then
        nc11 = nc11 + 1
        lthf3(nc11,7) = i
     endif
     if (mod(ix-2,4).eq.0 .and. mod(iy-2,4).eq.0 ) then
        nc12 = nc12 + 1
        lthf3(nc12,8) = i
     endif
     !=====================================================!
#endif
     if (mod(ix,2).eq.0 .and. mod(iy,2).ne.0 ) then
        nc3 = nc3 + 1
        lthf2(nc3,1) = i
     endif
     if (mod(ix,2).ne.0 .and. mod(iy,2).eq.0 ) then
        nc4 = nc4 + 1
        lthf2(nc4,2) = i
     endif
  enddo
  if( l .gt. 1 ) then
      if (nc1.ne.lq/4 .or. nc2  .ne. lq/4 ) then
         write(6,*) 'error 1'
         stop
      endif
  end if
#ifdef BREAKUP_T
  if( mod(l,4) .ne. 0 ) then
      write(6,*) 'error 2'
      stop
  endif
#endif
end subroutine sltpf
