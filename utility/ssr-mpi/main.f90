program main
! calculate Ising spin correlation
! read jjcorrR.bin
! MC average jjcorrR
! output jjcorrR(|i-j|) in X, Y, and XY direction
! output chi = \sum_{i} jjcorrR(i) / L^2

    implicit none
    integer :: l, ltrot, lq, nnimax, nntmax, nnimax_hyb, nntmax_hyb, zmax
    real(8) :: weight_track

    integer, dimension(:), allocatable :: jjcorr_R
    real(8), dimension(:), allocatable :: jjcorr_X, jjcorr_XY, jjcorr_Y
    integer :: i, nn, n, nf, nt, iit, ibt, eof, n_re, nn_t, nn_i
    integer :: ncount, imj_nx, imj_ny, imj

    integer :: nx, ny, jx, jy, nc, ni, j, itmp, ntj
    real(8) :: chi
    
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
    allocate( jjcorr_R(lq) )
    allocate( jjcorr_X(l), jjcorr_Y(l), jjcorr_XY(l) )

    open( unit=1001, file='jjcorrx_mpi.bin', status='unknown' )
    open( unit=1002, file='jjcorry_mpi.bin', status='unknown' )
    open( unit=1003, file='jjcorrxy_mpi.bin', status='unknown')
    open( unit=1004, file='chi_mpi.bin', status='unknown')
        
    open (unit=30,file='jjcorrR.bin',status='unknown')
    nc = 0
    do
        !!! read configuration
        do i = 1, lq
            read(30,*,IOSTAT=eof) jjcorr_R(i)
        end do
        if(eof.lt.0) exit 
#IFDEF TEST
        do i = 1, lq
            write(*,'(e16.8)') jjcorr_R(i)
        end do
#ENDIF

        !!! count number of configuration
        nc = nc + 1

        jjcorr_Y(1:l) = jjcorr_R(lq:lq-l+1:-1) / dble( lq )
        do i = 1, l
            j = lq - (i-1)*l
            jjcorr_X(i) = jjcorr_R(j) / dble( lq )

            j = lq - (i-1)*l - i + 1
            jjcorr_XY(i) = jjcorr_R(j) / dble( lq )
        end do

        chi = 0.d0
        do imj = 1, lq
            chi = chi + jjcorr_R(imj)
        end do
            
        ! output
        write(1001,'(50e16.8)') jjcorr_X(1:l/2)
        write(1002,'(50e16.8)') jjcorr_Y(1:l/2)
        write(1003,'(50e16.8)') jjcorr_XY(1:l/2)
        write(1004,'(e16.8)') chi/dble(lq*lq)
    end do

    close(30)
    close(1001)
    close(1002)
    close(1003)
    close(1004)

    deallocate( jjcorr_XY, jjcorr_Y, jjcorr_X )
    deallocate( jjcorr_R )

end program main
