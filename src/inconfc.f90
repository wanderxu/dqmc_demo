subroutine inconfc

    use spring
    use blockc
    implicit none

	include 'mpif.h'
    ! local
	integer  status(mpi_status_size)
    integer, dimension(:,:),   allocatable ::  itmpu
    integer, dimension(:,:,:), allocatable ::  itmpk, itmpj
    integer, dimension(:), allocatable :: b2int
    integer :: iseed0, i, nn, n, nf, nt, iit, ibt, icount, nbits2int, eof, n_re
    real(dp) :: x

	call mpi_comm_size(mpi_comm_world,isize,ierr)
	call mpi_comm_rank(mpi_comm_world,irank,ierr)

    allocate (itmpu(lq,ltrot), itmpk(lq,2,ltrot), itmpj(lq,2,ltrot))

    allocate( nsigl_u(lq,ltrot) )
    allocate( nsigl_k(lq,2,ltrot) )
    allocate( nsigl_j(lq,2,ltrot) )
        
    if (irank .eq. 0 ) then
       !!open (unit=30,file='confin',     status='unknown')
       open (unit=30,file='confin', status='unknown', form='unformatted', access='sequential')
    endif

	if ( irank.eq.0 )  then
	   read(30) iseed0
	   if (iseed0.eq.0) then
          ! start from scratch
          lwarnup = .true.
          write(fout,'(a)') ' start from scratch, need warnup '
	      do n = 1,isize - 1
             ! setup node i and send data.
	         do nt = 1,ltrot
                 do i  = 1,lq
		             !!!do nn = 1,2
		             !!!   x = spring_sfmt_stream()

		             !!!   nsigl_k(i,nn,nt) = 1
		             !!!   if (x.gt.0.5) nsigl_k(i,nn,nt) = -1
		             !!!   x = spring_sfmt_stream()
		             !!!   nsigl_j(i,nn,nt) = 1
		             !!!   if (x.gt.0.5) nsigl_j(i,nn,nt) = -1
		             !!!enddo
		             x = spring_sfmt_stream()
		             nsigl_u(i,nt) = 1
		             if (x.gt.0.5) nsigl_u(i,nt) = -1
                 enddo
             enddo
	         call mpi_send(nsigl_u,lq*ltrot,mpi_integer, n, n+512,mpi_comm_world,ierr)
	         !!!call mpi_send(nsigl_k, 2*lq*ltrot,mpi_integer, n, n+1024,mpi_comm_world,ierr)
	         !!!call mpi_send(nsigl_j, 2*lq*ltrot,mpi_integer, n, n+1536,mpi_comm_world,ierr)

	      enddo
          !	set node zero.
	      do nt = 1,ltrot
              do i  = 1,lq
		          !!!do nn = 1,2
		          !!!   x = spring_sfmt_stream()
		          !!!   nsigl_k(i,nn,nt) = 1
		          !!!   if (x.gt.0.5) nsigl_k(i,nn,nt) = -1
		          !!!   x = spring_sfmt_stream()
		          !!!   nsigl_j(i,nn,nt) = 1
		          !!!   if (x.gt.0.5) nsigl_j(i,nn,nt) = -1
		          !!!enddo
		          x = spring_sfmt_stream()
		          nsigl_u(i,nt) = 1
		          if (x.gt.0.5) nsigl_u(i,nt) = -1
              enddo
          enddo
	   else
!!#IFDEF GEN_CONFC_LEARNING
           read(30) weight_track
!!#ENDIF
          ! read all confins from node 0. 
          lwarnup = .false.
          write(fout,'(a)') ' start from old conf, do not need warnup '
          !	setup node 0
          nbits2int = ltrot*lq/32
          if(mod(ltrot*lq,32).ne.0) nbits2int = nbits2int + 1
          allocate( b2int( nbits2int ) )

          do i = 1, nbits2int
              read(30,IOSTAT=eof) b2int(i)
          end do
          icount = -1
          do nt = 1,ltrot
             do i  = 1,lq
                  icount = icount + 1
                  iit = icount / 32 + 1
                  ibt = mod(icount,32)
                  nsigl_u(i,nt) = ibits( b2int(iit), ibt, 1 ) * 2 - 1
             enddo
          enddo

          do n = 1,isize - 1
             if(eof.lt.0) exit 
             do i = 1, nbits2int
                 read(30,IOSTAT=eof) b2int(i)
                 if(eof.lt.0) exit
             end do
             if(eof.lt.0) exit
             icount = -1
             do nt = 1,ltrot
                do i  = 1,lq
                     icount = icount + 1
                     iit = icount / 32 + 1
                     ibt = mod(icount,32)
                     itmpu(i,nt) = ibits( b2int(iit), ibt, 1 ) * 2 - 1
                enddo
             enddo

             call mpi_send(itmpu,lq*ltrot,mpi_integer, n,  n+512,mpi_comm_world,ierr)
             !!!call mpi_send(itmpk, 2*lq*ltrot,mpi_integer, n, n+1024,mpi_comm_world,ierr)
             !!!call mpi_send(itmpj, 2*lq*ltrot,mpi_integer, n, n+1536,mpi_comm_world,ierr)

          enddo

          ! if we do not have enough configurations, we have to copy configurations from master process
          if( eof .lt. 0 ) then
          do n_re = n, isize-1
              itmpu(:,:) = nsigl_u(:,:)
              call mpi_send(itmpu,lq*ltrot,mpi_integer, n_re,  n_re+512,mpi_comm_world,ierr)
          end do
          end if


	   endif
	else
	   call mpi_recv(nsigl_u, lq*ltrot, mpi_integer,0,  irank + 512,  mpi_comm_world,status,ierr)
	   !!!call mpi_recv(nsigl_k, 2*lq*ltrot, mpi_integer,0, irank + 1024, mpi_comm_world,status,ierr)
	   !!!call mpi_recv(nsigl_j, 2*lq*ltrot, mpi_integer,0,  irank + 1536, mpi_comm_world,status,ierr)
	endif

    call mpi_bcast( lwarnup, 1, mpi_logical, 0, mpi_comm_world, ierr )

    if (irank .eq. 0 ) then
      if(allocated(b2int)) deallocate(b2int)
       close(30)
    endif

    deallocate (itmpu, itmpk, itmpj )

end subroutine inconfc
