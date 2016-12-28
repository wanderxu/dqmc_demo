program main
! calculate Ising spin correlation
! read confout.bin
! for each bin, average over tau of s(i,t) to get s(i)
! MC average s(i)*s(j) to get jjcorrR(i)
! output jjcorrR(i) in X, Y, and XY direction
! output chi = \sum_{i} jjcorrR(i) / L^2
    implicit none
    integer :: l, ltrot, lq, nnimax, nntmax, nnimax_hyb, nntmax_hyb, zmax
    real(8) :: weight_track

    integer, dimension(:,:), allocatable :: list, invlist
    integer, dimension(:,:), allocatable ::  nsigl_u, latt_imj
    integer, dimension(:), allocatable :: b2int, nsiglR, jjcorr_R
    real(8), dimension(:), allocatable :: jjcorr_X, jjcorr_XY, jjcorr_Y
    integer :: i, nn, n, nf, nt, iit, ibt, icount, nbits2int, eof, n_re, nn_t, nn_i
    integer :: ncount, imj_nx, imj_ny, imj

    integer :: nx, ny, jx, jy, nc, ni, j, itmp, ntj
    real(8) :: chi
    
    integer, external :: npbc

    ! nnimax:  spatial
    ! nntmax:  tempeoral
    ! nnimax_hyb
    ! nntmax_hyb
    ! read in parameters
    open (unit=40,file='in.para',status='unknown')
    read(40,*) l, ltrot
    close(40)
    lq = l*l
#IFDEF TEST
    write(*,'(a,i6)')  ' l = ', l
    write(*,'(a,i6)')  ' lq = ', lq
    write(*,'(a,i6)')  ' ltrot  = ', ltrot
#ENDIF

    ! allocate data
    allocate( list(lq,2) )
    allocate( invlist(l,l) )
    allocate( latt_imj(lq,lq) )
    allocate( nsigl_u(lq,ltrot) )
    allocate( nsiglR(lq) )
    allocate( jjcorr_R(lq) )
    allocate( jjcorr_X(l), jjcorr_Y(l), jjcorr_XY(l) )
    nbits2int = ltrot*lq/32
    if(mod(ltrot*lq,32).ne.0) nbits2int = nbits2int + 1
    allocate( b2int( nbits2int ) )

    !! set list
    ! list, invlist
    ncount = 0
    do nx = 1,l
    do ny = 1,l
       ncount = ncount + 1
       list(ncount,1) = nx
       list(ncount,2) = ny
       invlist(nx,ny) = ncount
    enddo
    enddo
    ! latt_imj
    do j = 1, lq
        do i = 1, lq
            imj_nx = npbc( list(i,1) - list(j,1), l )
            imj_ny = npbc( list(i,2) - list(j,2), l )
            latt_imj(i,j) = invlist( imj_nx, imj_ny )
        end do
    end do

    open( unit=1001, file='jjcorrx.bin', status='unknown' )
    open( unit=1002, file='jjcorry.bin', status='unknown' )
    open( unit=1003, file='jjcorrxy.bin', status='unknown')
    open( unit=1004, file='chi.bin', status='unknown')
        
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

        !! first average over time
        nsiglR(:) = 0
        do nt = 1, ltrot
            do i = 1, lq
                nsiglR(i) = nsiglR(i) + nsigl_u(i,nt)
            end do
        end do

        !!! calculate spin-spin interaction
        jjcorr_R(:) = 0
        do j = 1, lq
            do i = 1, lq
                imj = latt_imj(i,j)
                jjcorr_R(imj) = jjcorr_R(imj) + nsiglR(i)*nsiglR(j)
            end do
        end do

        jjcorr_Y(1:l) = dble( jjcorr_R(lq:lq-l+1:-1) ) / dble( lq*ltrot*ltrot )
        do i = 1, l
            j = lq - (i-1)*l
            jjcorr_X(i) = dble( jjcorr_R(j) ) / dble( lq*ltrot*ltrot )

            j = lq - (i-1)*l - i + 1
            jjcorr_XY(i) = dble( jjcorr_R(j) ) / dble( lq*ltrot*ltrot )
        end do

        chi = 0.d0
        do imj = 1, lq
            chi = chi + dble(jjcorr_R(imj))
        end do
            
        ! output
        write(1001,'(50e16.8)') jjcorr_X(1:l/2)
        write(1002,'(50e16.8)') jjcorr_Y(1:l/2)
        write(1003,'(50e16.8)') jjcorr_XY(1:l/2)
        write(1004,'(e16.8)') dble(chi)/dble(lq*lq)/dble(ltrot*ltrot)
    end do

    if(allocated(b2int)) deallocate(b2int)
    close(30)
    close(1001)
    close(1002)
    close(1003)
    close(1004)

    deallocate( jjcorr_XY, jjcorr_Y, jjcorr_X )
    deallocate( jjcorr_R )
    deallocate( nsiglR )
    deallocate( nsigl_u )
    deallocate( latt_imj )
    deallocate( invlist )
    deallocate( list )

end program main

integer function npbc(nr,l)
  implicit none
  integer, intent(in) :: nr
  integer, intent(in) :: l
  npbc = nr
  if (nr.gt.l) npbc = nr - l
  if (nr.lt.1) npbc = nr + l
end function npbc
