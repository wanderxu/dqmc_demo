program main
! calculate Ising spin correlation
! read jjcorrRtau.bin
! MC average \sum_i jjcorrRtau(i,t) /L^2  to get chitau(t), and also error for chitau(t)
! output chitau(t) and chitauerr(t)

    implicit none
    integer :: l, ltrot, lq, nnimax, nntmax, nnimax_hyb, nntmax_hyb, zmax
    real(8) :: weight_track

    real(8), dimension(:,:), allocatable :: jjcorr_Rtau
    real(8), dimension(:), allocatable :: jjcorr_X, jjcorr_XY, jjcorr_Y
    integer :: i, nn, n, nf, nt, iit, ibt, eof, n_re, nn_t, nn_i
    integer :: ncount, imj_nx, imj_ny, imj

    integer :: nx, ny, jx, jy, nc, ni, j, itmp, ntj, nb
    real(8), dimension(:), allocatable :: chitau, chitau2, chitauerr

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
    allocate( chitau(ltrot+1) )
    allocate( chitauerr(ltrot+1) )
    allocate( chitau2(ltrot+1) )

    open( unit=1001, file='chitau.dat', status='unknown' )
        
    open (unit=30,file='jjcorrRtau.bin',status='unknown')
    nc = 0
    chitau(:) = 0.d0
    chitau2(:) = 0.d0
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

        do nt = 1, ltrot/2+1
            do i= 1, lq
                chitau(nt) = chitau(nt) + jjcorr_Rtau(i,nt) /dble(ltrot)
                chitau2(nt) = chitau2(nt) + jjcorr_Rtau(i,nt)*jjcorr_Rtau(i,nt)/dble(ltrot)/dble(ltrot)
            end do
        end do
    end do

    do nt = 1, ltrot/2+1
        chitau(nt) = chitau(nt)/dble(nc)
        chitau2(nt) = chitau2(nt)/dble(nc)
        chitauerr(nt) = dsqrt( dabs( chitau2(nt) - chitau(nt)*chitau(nt) )/ dble(nc*20*28-1) )
        chitau(nt) = chitau(nt) / dble(lq)
        chitauerr(nt) = chitauerr(nt) / dble(lq)
    end do
    do nt = ltrot/2+2, ltrot+1
        chitau(nt) = chitau( ltrot+2-nt )
        chitauerr(nt) = chitauerr( ltrot+2-nt )
    end do

    do nt = 1, ltrot+1
        write(1001, '(3e16.8)') dble(nt-1)*0.05d0, chitau(nt), chitauerr(nt)
    end do

    close(30)
    close(1001)

    deallocate( chitau2, chitauerr, chitau )
    deallocate( jjcorr_Rtau )

end program main
