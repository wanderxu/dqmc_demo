program main

    implicit none
    integer :: l, ltrot, lq, nnimax, nntmax, zmax
    real(8) :: weight_track

    integer, dimension(:,:),   allocatable ::  nsigl_u
    integer, dimension(:), allocatable :: b2int, jjicorr, jjtcorr
    integer, dimension(:,:,:), allocatable :: nntable
    integer :: iseed0, i, nn, n, nf, nt, iit, ibt, icount, nbits2int, eof, n_re

    integer :: nx, ny, jx, jy, nc, ni, j, itmp, ntj

    ! read in parameters
    open (unit=40,file='in.para',status='unknown')
    read(40,*) l, ltrot
    close(40)
    lq = l*l
    nnimax = 2  ! spatial
    nntmax = 1  ! tempeoral
    zmax = 16
#IFDEF TEST
    write(*,'(a,i6)')  ' l = ', l
    write(*,'(a,i6)')  ' lq = ', lq
    write(*,'(a,i6)')  ' ltrot = ', ltrot
    write(*,'(a,i6)')  ' nnimax = ', nnimax
    write(*,'(a,i6)')  ' zmax = ', zmax
#ENDIF

    ! allocate data
    allocate( nsigl_u(lq,ltrot) )
    allocate( nntable(2, zmax, nnimax) )
    allocate( jjicorr(nnimax) )
    allocate( jjtcorr(nntmax) )
    nbits2int = ltrot*lq/32
    if(mod(ltrot*lq,32).ne.0) nbits2int = nbits2int + 1
    allocate( b2int( nbits2int ) )

    ! generate neightbor table
    call generate_neighbor(zmax,nnimax,nntable)
#IFDEF TEST
    write(*,*)
    write(*,'(a)') '  nn    ni     nntable(:,ni,i) '
    do i = 1, nnimax
        do ni = 1, zmax
            write(*,'(4i6)') i, ni, nntable(:,ni,i)
        end do
    end do
#ENDIF
        
    open (unit=30,file='confout.bin', status='unknown', form='unformatted', access='sequential')
    nc = 0
    do
        !!! read weight
        read(30,IOSTAT=eof) itmp
        read(30,IOSTAT=eof) weight_track
        if(eof.lt.0) exit 

        !!! read configuration
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

#IFDEF TEST
        do nt = 1, ltrot
            write(*,'(20i4)') nsigl_u(:,nt)
        end do
#ENDIF

        !!! count number of configuration
        nc = nc + 1

        !!! calculate interaction energy
        ! spatial
        jjicorr(:) = 0
        do nn = 1, nnimax
            do nt = 1, ltrot
            do nx = 1,l
            do ny = 1,l
                i = (nx-1)*l + ny
                do ni = 1, zmax
                    if( nntable(1,ni,nn) .eq.0 .and. nntable(2,ni,nn) .eq. 0 ) exit
                    jx = nx + nntable(1,ni,nn)
                    jy = ny + nntable(2,ni,nn)
                    if(jx.gt.l) jx = mod(jx,l)
                    if(jx.lt.1) jx = mod(jx,l)+l
                    if(jy.gt.l) jy = mod(jy,l)
                    if(jy.lt.1) jy = mod(jy,l)+l
                    j = (jx-1)*l + jy
                    jjicorr(nn) = jjicorr(nn) + nsigl_u(i,nt) * nsigl_u(j,nt)
                end do
            end do
            end do
            end do
        end do

        ! temporal
        jjtcorr(:) = 0
        do nn = 1, nntmax
            do i = 1, lq
                do nt = 1, ltrot
                    ntj = nt+nn
                    if(ntj>ltrot) ntj = mod(ntj,ltrot)
                    jjtcorr(nn) = jjtcorr(nn) + nsigl_u(i,nt) * nsigl_u(i,ntj)
                end do
            end do
        end do

        ! output
        write(*,'(e16.8,10i20)') weight_track, jjicorr(:)/2, jjtcorr(:)
    end do
#IFDEF TEST
    write(*,'(a,i6)')  ' total configuration = ', nc
#ENDIF

    if(allocated(b2int)) deallocate(b2int)
    close(30)

    deallocate (nsigl_u )
    deallocate( jjicorr )
    deallocate( jjtcorr )
    deallocate( nntable )

end program main
