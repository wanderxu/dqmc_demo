module mod_neighbor
    use blockc
    implicit none
    integer, save :: nnimax, nntmax, nnimax_hyb, nntmax_hyb, zmax, num_nei
    integer, dimension(:,:,:,:), allocatable, save :: nei_cord
    real(dp), dimension(:,:,:), allocatable, save :: nei_Jeff
    real(dp), dimension(:), allocatable, save :: Jeff_i, Jeff_t, Jeff_it

    contains

    subroutine deallocate_neighbor
        implicit none
        deallocate( nei_Jeff )
        deallocate( nei_cord )
    end subroutine

    subroutine set_neighbor
    
        implicit none
        integer, dimension(:,:,:), allocatable :: nntable
        integer :: iseed0, i, nn, n, nf, nt, iit, ibt, icount, nbits2int, eof, n_re, nn_t, nn_i
    
        integer :: nx, ny, jx, jy, nc, ni, j, itmp, ntj
    
        zmax = 16
        ! read in parameters
        open (unit=40,file='Heff.para',status='unknown')
        read(40,*) nnimax, nntmax, nnimax_hyb, nntmax_hyb
        do nn = 1, nnimax
            read(40,*) Jeff_i(nn)
        end do
        do nn = 1, nntmax
            read(40,*) Jeff_t(nn)
        end do
        do nn_i = 1, nnimax_hyb
            do nn_t = 1, nntmax_hyb
                nn = (nn_i-1)*nntmax_hyb + nn_t
                read(40,*) Jeff_it(nn)
            end do
        end do
        close(40)
    
        ! allocate data
        allocate( nntable(2, zmax, nnimax) )
    
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

        !!! count number of neighbor
        ! count spatial and temporal
        inn = 0
        do nn = 1, nnimax ! spatial part
            do ni = 1, zmax
                if( nntable(1,ni,nn) .eq.0 .and. nntable(2,ni,nn) .eq. 0 ) exit
                inn = inn + 1
            end do
        end do
        nei_num = inn + nntmax*2

        ! count spatial temporal hybird
        inn = 0
        do nn = 1, nnimax_hyb
            do ni = 1, zmax ! spatial part
                if( nntable(1,ni,nn) .eq.0 .and. nntable(2,ni,nn) .eq. 0 ) exit
                inn = inn + 1
            end do
        end do
        nei_num = nei_num + inn*nntmax*2

        allocate( nei_cord(2,num_nei,lq,ltrot) )
        allocate( nei_Jeff(num_nei,lq,ltrot) )
        
        !!!  set nei_cord, nei_Jeff
        ! spatial
        nncount = 0
        do nt = 1, ltrot
        do nx = 1,l
        do ny = 1,l
            i = (nx-1)*l + ny
            inn = nncount
            do nn = 1, nnimax
                do ni = 1, zmax
                    if( nntable(1,ni,nn) .eq.0 .and. nntable(2,ni,nn) .eq. 0 ) exit
                    jx = nx + nntable(1,ni,nn)
                    jy = ny + nntable(2,ni,nn)
                    if(jx.gt.l) jx = mod(jx,l)
                    if(jx.lt.1) jx = mod(jx,l)+l
                    if(jy.gt.l) jy = mod(jy,l)
                    if(jy.lt.1) jy = mod(jy,l)+l
                    j = (jx-1)*l + jy
                    inn = inn + 1
                    nei_cord(1,inn,i,nt) = j
                    nei_cord(2,inn,i,nt) = nt
                    nei_Jeff(inn,i,nt) = Jeff_i(nn)
                end do
            end do
        end do
        end do
        end do
        nncount = inn
    
        ! temporal
        do i = 1, lq
            do nt = 1, ltrot
                do nn = 1, nntmax
                    ntj = nt+nn
                    if(ntj>ltrot) ntj = mod(ntj,ltrot)
                    inn = nncount + 1
                    nei_cord(1,inn,i,nt) = i
                    nei_cord(2,inn,i,nt) = ntj
                    nei_Jeff(inn,i,nt) = Jeff_t(nn)

                    ntj = nt-nn
                    if(ntj<1) ntj = mod(ntj,ltrot) + ltrot
                    inn = nncount + 2
                    nei_cord(1,inn,i,nt) = i
                    nei_cord(2,inn,i,nt) = ntj
                    nei_Jeff(inn,i,nt) = Jeff_t(nn)
                end do
            end do
        end do
        nncount = inn
    
        ! spatial temporal hybrid
        do nt = 1, ltrot
        do nx = 1,l
        do ny = 1,l
            i = (nx-1)*l + ny
            inn = nncount
            do nn_i = 1, nnimax_hyb
                do ni = 1, zmax
                    if( nntable(1,ni,nn_i) .eq.0 .and. nntable(2,ni,nn_i) .eq. 0 ) exit
                    jx = nx + nntable(1,ni,nn_i)
                    jy = ny + nntable(2,ni,nn_i)
                    if(jx.gt.l) jx = mod(jx,l)
                    if(jx.lt.1) jx = mod(jx,l)+l
                    if(jy.gt.l) jy = mod(jy,l)
                    if(jy.lt.1) jy = mod(jy,l)+l
                    j = (jx-1)*l + jy
                    do nn_t = 1, nntmax_hyb
    
                        nn = (nn_i-1)*nntmax_hyb + nn_t
    
                        ntj = nt + nn_t
                        if(ntj>ltrot) ntj = mod(ntj, ltrot)
                        inn = inn + 1
                        nei_cord(1,inn,i,nt) = j
                        nei_cord(2,inn,i,nt) = ntj
                        nei_Jeff(inn,i,nt) = Jeff_it(nn)
    
                        ntj = nt - nn_t
                        if(ntj<1) ntj = mod(ntj, ltrot) + ltrot
                        inn = inn + 1
                        nei_cord(1,inn,i,nt) = j
                        nei_cord(2,inn,i,nt) = ntj
                        nei_Jeff(inn,i,nt) = Jeff_it(nn)
                    end do
                end do
            end do
        end do
        end do
        end do

        if( nncount .ne. num_nei ) stop 'wrong in count num_nei'
    
        deallocate( nntable )
    
    end subroutine set_neighbor
end module mod_neighbor
