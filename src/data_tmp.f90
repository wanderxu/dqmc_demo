module data_tmp
  use blockc, only: ndim, dp
  implicit none
  complex(dp), dimension(:), allocatable :: v1, v2, v3, v4, v5, v6, v7, v8
  complex(dp), dimension(:), allocatable :: vec1, vec2, vhlp1, uhlp1, vhlp2, uhlp2, u1, u2
  complex(dp), dimension(:,:), allocatable :: Atmp, Btmp, Vtmp, vvtmp, uutmp, dvvtmp, dvvdtmp, grtmp, gt0tmp, g0ttmp

  contains

  subroutine allocate_data_tmp
    implicit none
    allocate( v1(ndim), v2(ndim), v3(ndim), v4(ndim), v5(ndim), v6(ndim), v7(ndim), v8(ndim) ) ! 1
    allocate( vec1(ndim), vec2(ndim), vhlp1(ndim), uhlp1(ndim), vhlp2(ndim), uhlp2(ndim), u1(ndim), u2(ndim) )
    allocate( Atmp(ndim,ndim) )    ! 2
    allocate( Btmp(ndim,ndim) )    ! 3
    allocate( Vtmp(ndim,ndim) )    ! 4
    allocate( vvtmp(ndim,ndim) )   ! 5
    allocate( uutmp(ndim,ndim) )   ! 5
    allocate( dvvtmp(ndim,ndim) )  ! 6
    allocate( dvvdtmp(ndim,ndim) ) ! 7
    allocate( grtmp(ndim,ndim) )   ! 8
    allocate( gt0tmp(ndim,ndim) )   ! 8
    allocate( g0ttmp(ndim,ndim) )   ! 8
  end subroutine allocate_data_tmp

  subroutine deallocate_data_tmp
    implicit none
    deallocate( g0ttmp )   ! 8
    deallocate( gt0tmp )   ! 8
    deallocate( grtmp )   ! 8
    deallocate( dvvdtmp ) ! 7
    deallocate( dvvtmp )  ! 6
    deallocate( uutmp )   ! 5
    deallocate( vvtmp )   ! 5
    deallocate( Vtmp )    ! 4
    deallocate( Btmp )    ! 3
    deallocate( Atmp )    ! 2
    deallocate( u2, u1, uhlp2, vhlp2, uhlp1, vhlp1, vec2, vec1 )
    deallocate( v8, v7, v6, v5, v4, v3, v2, v1 ) ! 1
  end subroutine deallocate_data_tmp

end module data_tmp
