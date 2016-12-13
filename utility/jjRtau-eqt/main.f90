program main
! calculate Ising spin correlation
! read jjcorrRtau.bin
! MC average jjcorrRtau(i,t=0)  to get chi(i), and also error for chi(i)
! output chi(i) and chierr(i) in X direction

    implicit none
    integer :: l, ltrot, lq, nnimax, nntmax, nnimax_hyb, nntmax_hyb, zmax
    real(8) :: weight_track

    real(8), dimension(:,:), allocatable :: jjcorr_Rtau
    real(8), dimension(:), allocatable :: jjcorr_X, jjcorr_X2, jjcorr_Xerr
    integer :: i, nn, n, nf, nt, iit, ibt, eof, n_re, nn_t, nn_i
    integer :: ncount, imj_nx, imj_ny, imj

    integer :: nx, ny, jx, jy, nc, ni, j, itmp, ntj, nb
    real(8), dimension(:), allocatable :: chir, chir2, chirerr

    integer, parameter :: nskip = 5
    
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
    allocate( jjcorr_Rtau(lq,ltrot) )
    allocate( chir(lq) )
    allocate( chirerr(lq) )
    allocate( chir2(lq) )
    allocate( jjcorr_X(l) )
    allocate( jjcorr_X2(l) )
    allocate( jjcorr_Xerr(l) )

    open( unit=1001, file='jjRtau_sx.dat', status='unknown' )
        
    open (unit=30,file='jjcorrRtau.bin',status='unknown')
    nc = 0
    chir(:) = 0.d0
    chir2(:) = 0.d0
    do
        !!! read configuration
        do nt = 1, ltrot/2+1
        do i = 1, lq
            read(30,*,IOSTAT=eof) jjcorr_Rtau(i,nt)
        end do
        end do
        if(eof.lt.0) exit 
        !!! count number of configuration
        nc = nc + 1
    end do
    close(30)

    write(*,'(a,i6)') " total number of bins = ", nc


    write(*,'(a,i4,a)') " skip ", nskip , " bins"
    open (unit=30,file='jjcorrRtau.bin',status='unknown')
    do nb = 1, nskip
        !!! read configuration
        do nt = 1, ltrot/2+1
        do i = 1, lq
            read(30,*,IOSTAT=eof) jjcorr_Rtau(i,nt)
        end do
        end do
    end do

    ! read bins
    nc = nc - nskip
    do nb = 1, nc
        !!! read configuration
        do nt = 1, ltrot/2+1
        do i = 1, lq
            read(30,*,IOSTAT=eof) jjcorr_Rtau(i,nt)
        end do
        end do
        if(eof.lt.0) exit 

        do i= 1, lq
            chir(i) = chir(i) + jjcorr_Rtau(i,1) /dble(ltrot)
            chir2(i) = chir2(i) + jjcorr_Rtau(i,1)*jjcorr_Rtau(i,1)/dble(ltrot)/dble(ltrot)
        end do
    end do

    do i = 1, l
        j = lq - i + 1
        jjcorr_X(i) = chir(j) / dble(nc)
        jjcorr_X2(i) = chir2(j)/ dble(nc)
    end do
    do i = 1, l
        jjcorr_Xerr(i) = dsqrt( dabs( jjcorr_X2(i) - jjcorr_X(i)*jjcorr_X(i) )/ dble(nc*20*28-1) )
        jjcorr_X(i) = jjcorr_X(i) / dble(lq)
        jjcorr_Xerr(i) = jjcorr_Xerr(i) / dble(lq)
    end do

    do i = 1, l/2
        write(1001, '(i6,2e16.8)') i-1, jjcorr_X(i), jjcorr_Xerr(i)
    end do

    close(30)
    close(1001)

    deallocate( jjcorr_Xerr )
    deallocate( jjcorr_X2 )
    deallocate( jjcorr_X )
    deallocate( chir2, chirerr, chir )
    deallocate( jjcorr_Rtau )

end program main
