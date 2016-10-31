module mod_cumulate
    use blockc
    implicit none
    integer, save :: nnimax, nntmax, nnimax_hyb, nntmax_hyb, zmax, num_nei
    integer, dimension(:,:,:,:), allocatable, save :: nei_cord
    real(dp), dimension(:,:,:), allocatable, save :: nei_Jeff
    real(dp), dimension(:), allocatable, save :: Jeff_i, Jeff_t, Jeff_it
    real(dp), dimension(:,:), allocatable, save :: heff

    integer, save :: ncumulate

    contains

    subroutine deallocate_cumulate
        implicit none
        deallocate( heff )
        deallocate( nei_Jeff )
        deallocate( nei_cord )
        deallocate( Jeff_it )
        deallocate( Jeff_t )
        deallocate( Jeff_i )
    end subroutine deallocate_cumulate

    subroutine initial_heff
      implicit none
      integer :: inn, j, ntj, i, nt
      allocate( heff(lq,ltrot) )
      heff(:,:) = 0.d0
      do nt = 1, ltrot
          do i = 1, lq
              do inn = 1, num_nei
                  j = nei_cord(1,inn,i,nt)
                  ntj = nei_cord(2,inn,i,nt)
                  heff(i,nt) = heff(i,nt) + nei_Jeff(inn,i,nt)*nsigl_u(j,ntj)
              end do
          end do
      end do
#IFDEF TEST
      write(*,'(a)') ' nsigl_u = '
      do nt = 1, ltrot
          write(*,'(40i7)') nsigl_u(:,nt)
      end do
      write(*,'(a)') ' after intial_heff, heff = '
      do nt = 1, ltrot
          write(*,'(40f7.3)') heff(:,nt)
      end do
#ENDIF
    end subroutine initial_heff

    subroutine set_neighbor
    
        implicit none
        integer, dimension(:,:,:), allocatable :: nntable
        integer :: iseed0, i, nn, n, nf, nt, iit, ibt, icount, nbits2int, eof, n_re, nn_t, nn_i
    
        integer :: nx, ny, jx, jy, nc, ni, j, itmp, ntj, inn, nncount

    
        zmax = 16
        ! read in parameters
        open (unit=40,file='Heff.para',status='unknown')
        read(40,*) nnimax, nntmax, nnimax_hyb, nntmax_hyb

        allocate( Jeff_i( nnimax ) )
        allocate( Jeff_t( nntmax ) )
        allocate( Jeff_it( nnimax_hyb*nntmax_hyb ) )

        read(40,*) ncumulate
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
        num_nei = inn + nntmax*2

        ! count spatial temporal hybird
        inn = 0
        do nn = 1, nnimax_hyb
            do ni = 1, zmax ! spatial part
                if( nntable(1,ni,nn) .eq.0 .and. nntable(2,ni,nn) .eq. 0 ) exit
                inn = inn + 1
            end do
        end do
        num_nei = num_nei + inn*nntmax_hyb*2
#IFDEF TEST
        write(*,'(a,i8)') ' first count, num_nei = ', num_nei
#ENDIF
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
#IFDEF TEST
        write(*,'(a,i8)') ' after spatial, nncount = ', nncount
#ENDIF
    
        ! temporal
        do i = 1, lq
            do nt = 1, ltrot
                inn = nncount
                do nn = 1, nntmax
                    ntj = nt+nn
                    if(ntj>ltrot) ntj = mod(ntj,ltrot)
                    inn = inn + 1
                    nei_cord(1,inn,i,nt) = i
                    nei_cord(2,inn,i,nt) = ntj
                    nei_Jeff(inn,i,nt) = Jeff_t(nn)

                    ntj = nt-nn
                    if(ntj<1) ntj = mod(ntj,ltrot) + ltrot
                    inn = inn + 1
                    nei_cord(1,inn,i,nt) = i
                    nei_cord(2,inn,i,nt) = ntj
                    nei_Jeff(inn,i,nt) = Jeff_t(nn)
                end do
            end do
        end do
        nncount = inn
#IFDEF TEST
        write(*,'(a,i8)') ' after temporal, nncount = ', nncount
#ENDIF
    
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
        nncount = inn
#IFDEF TEST
        write(*,'(a,i8)') ' after spatial-temporal-hybrid, nncount = ', nncount
        write(*,'(a,i8)') ' recheck count, num_nei = ', nncount
        write(*,*)

        write(*,'(a)') ' nei_cord = '
        do nt = 1, ltrot
            do i = 1, lq
                write(*,"(6('(',i2,',',i2,')  '))") nei_cord(:,1,i,nt), &
                                                nei_cord(:,2,i,nt), &
                                                nei_cord(:,3,i,nt), &
                                                nei_cord(:,4,i,nt), &
                                                nei_cord(:,5,i,nt), &
                                                nei_cord(:,6,i,nt)
            end do
        end do

        write(*,'(a)') ' nei_Jeff = '
        do nt = 1, ltrot
            do i = 1, lq
                write(*,'(40f8.4)') nei_Jeff(:,i,nt)
            end do
        end do
#ENDIF
        if( nncount .ne. num_nei ) stop 'wrong in count num_nei'
    
        deallocate( nntable )
    
    end subroutine set_neighbor
end module mod_cumulate
