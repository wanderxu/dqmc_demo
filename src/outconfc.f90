subroutine outconfc

    use blockc
    implicit none

	include 'mpif.h'

    ! local
    integer, dimension(:), allocatable :: b2int
    integer, dimension(:,:), allocatable ::  itmpu
	integer  status(mpi_status_size)
    integer :: i, n, nf, nt, iit, ibt, icount, itmp, nbits2int 

	call mpi_comm_size(mpi_comm_world,isize,ierr)
	call mpi_comm_rank(mpi_comm_world,irank,ierr)

    allocate ( itmpu(lq,ltrot) )
	
    if (irank.eq.0) then
       !!!open (unit=35, file='confout', status='unknown')
       open (unit=35,file='confout', status='unknown', form='unformatted', access='sequential')
    endif

	if ( irank.ne.0 )  then
	   call mpi_send(nsigl_u,lq*ltrot,mpi_integer, 0, irank+512,mpi_comm_world,ierr)

	   !!!call mpi_send(nsigl_k, 2*lq*ltrot,mpi_integer, 0, irank+1024,mpi_comm_world,ierr)

	   !!!call mpi_send(nsigl_j, 2*lq*ltrot,mpi_integer, 0, irank+1536,mpi_comm_world,ierr)
	endif
	if (irank.eq.0)  then

       write(35) 1
#IFDEF GEN_CONFC_LEARNING
       write(35) weight_track
#ENDIF

       nbits2int = ltrot*lq/32
       if(mod(ltrot*lq,32).ne.0) nbits2int = nbits2int + 1
       allocate( b2int( nbits2int ) )
       b2int = 0

       icount = -1
       do nt = 1,ltrot
          do i  = 1,lq
                icount = icount + 1
                iit = icount / 32 + 1
                ibt = mod(icount,32)
                if( nsigl_u(i,nt) .eq. 1 ) b2int(iit) = ibset( b2int(iit), ibt )
          enddo
       enddo

       do i = 1, nbits2int
           write(35) b2int(i)
       end do

       do n = 1,isize - 1
	      call mpi_recv(itmpu,lq*ltrot, mpi_integer,n, n+512, mpi_comm_world,status,ierr)
          b2int = 0
          icount = -1
          do nt = 1,ltrot
             do i  = 1,lq
                   icount = icount + 1
                   iit = icount / 32 + 1
                   ibt = mod(icount,32)
                   if( itmpu(i,nt) .eq. 1 ) b2int(iit) = ibset( b2int(iit), ibt )
             enddo
          enddo

          do i = 1, nbits2int
              write(35) b2int(i)
          end do

       enddo
    endif


    if (irank.eq.0) then
       if( allocated(b2int) ) deallocate(b2int)
       close(35)
    endif

    deallocate ( itmpu )

end subroutine outconfc
