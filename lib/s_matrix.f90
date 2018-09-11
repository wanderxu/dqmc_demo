  !!! collections of subroutines for matrix operations
  !!! written by Xiao Yan Xu
  !!! email: wanderxu@gmail.com
  !  
  ! subroutine s_identity_z(n, A) ! build complex identity matrix
  ! subroutine s_logdet_z(ndim, zmat, zlogdet) ! logdet of a complex matrix
  ! subroutine s_eig_he(ldim, ndim, amat, eval, evec) ! computes all eigenvalues and eigenvectors of a complex Hermitian matrix
  ! subroutine s_svd_zg(m, n, min_mn, amat, umat, svec, vmat) ! UDV decomposition of a complex matrix
  ! subroutine s_svd_zg_logdetU(m, n, min_mn, amat, umat, svec, vmat, zlogdetU) ! UDV decomposition of a complex matrix, also return logdet(U)
  ! subroutine s_compare_max_z( ndim, Amat, Bmat, max_diff ) ! elements maxdiff of two complex matrices
  ! subroutine s_compare_max_d( ndim, Amat, Bmat, max_diff ) ! elements maxdiff of two real matrices
  ! subroutine s_z_x_diag_d( ndim, Amat, dvec, Bmat ) ! product of a complex matrix and a real diagonal matrix
  ! subroutine s_diag_d_x_z( ndim, dvec, Amat, Bmat ) ! product of a diagonal matrix and a complex matrix
  ! subroutine s_diag_dvd( ndim, dvecr, Amat, dvecl, Bmat ) ! product of a real diagonal matrix, a complex matrix and another real diagonal matrix
  ! subroutine s_v_invd_u( ndim, vmat, dvec, umat, zmat ) ! product of a complex matrix, inverse a real diagonal matrix and another  complex matrix
  ! subroutine s_v_d_u( ndim, vmat, dvec, umat, zmat ) ! product of a complex matrix, a real diagonal matrix and another complex matrix
  ! subroutine s_dvec_min_max(ndim,dvec,dmax,dmin) ! decomposition of a real diagonal matrix as two diagonal matrices, one with elements >= 1, the other <=1
  ! subroutine s_zmcpt(rdim,cdim,amat,jpvt) ! make matrix amat to column norm decreasing
  ! subroutine s_zmrpt(rdim,cdim,amat,jpvt) ! make matrix amat to row norm decreasing
  ! subroutine s_zpm(rdim,cdim,jpv,amat) ! perform permutation of a complex matrix from left side
  ! subroutine s_pp_R(ndim,jpv1,jpv2) ! merge two right side permutations
  ! subroutine s_pp_L(ndim,jpv1,jpv2) ! merge two left side permutations
  ! subroutine s_zmp(rdim,cdim,amat,jpv) ! perform permutation of a complex matrix from right side
  ! subroutine s_dppd_z(ndim, dvecr, jpvrin, jpvlin, dvecl, amat) ! product of a real diagonal matrix, two permutations and another real diagonal matrix
  ! subroutine s_zgeQRPT(rdim, cdim, amat, qmat, rmat, jpvt) ! QR decomposition with pivoting
  ! subroutine s_zgeQRPT_logdetQ(rdim, cdim, amat, qmat, rmat, jpvt, zlogdet) ! QR decomposition with pivoting, also return logdet(Q)
  ! subroutine s_invdiag_d_x_zr( ndim, dvec, Amat, Bmat ) ! product of inverse a real diagonal and a complex matrix
  ! subroutine s_zgeQR(rdim, cdim, amat, qmat, rmat ) ! QR decomposition
  ! subroutine s_zgeQR_logdetQ(rdim, cdim, amat, qmat, rmat, zlogdet ) ! QR decomposition, also return logdet(Q)
  ! subroutine s_invqr_z(ndim, amat ) ! inverse of a complex matrix, using QR decomposition
  ! subroutine s_inv_logdet_qr_z(ndim, amat, zlogdet ) ! inverse of a complex matrix, using QR decomposition, also return logdet of this matrix
  ! subroutine s_inv_det_qr_z(ndim, amat, zdet ) ! inverse of a complex matrix, using QR decomposition, also return det of this matrix
  ! subroutine s_invlu_z(ndim, zmat) ! inverse of a complex matrix, using LU decomposition
  ! subroutine s_inv_logdet_lu_z(ndim, zmat, zlogdet) ! inverse of a complex matrix, using LU decomposition, also return logdet of this matrix
  ! subroutine s_adfac_z( rdim, cdim, amat, dvec ) ! factorize a complex matrix into product of a dense matrix and a real diagonal matrix with norm of couumn
  ! subroutine s_dafac_z( rdim, cdim, dvec, amat ) ! factorize a complex matrix into prodcut of a real diagonal matrix with norm of row and a dense matrix

  subroutine s_identity_z(n, A)
  !! build complex identity matrix
     use constants, only : dp, czero, cone
     implicit none
     integer, intent(in)      :: n  ! size of matrix
     complex(dp), intent(out) :: A(n,n)

     integer :: i
     A = czero
     do i=1,n
         A(i,i) = cone
     end do
     return
  end subroutine s_identity_z

  subroutine s_logdet_z(ndim, zmat, zlogdet)
  !! log of determinant of a complex matrix
     use constants, only : dp, czero, pi
     implicit none
     integer, intent(in)        :: ndim
     complex(dp), intent(out)   :: zlogdet

     ! object matrix, on entry, it contains the original matrix, on exit,
     ! it is destroyed and replaced with the L and U matrix
     complex(dp), intent(in) :: zmat(ndim,ndim)

     integer :: i
     real(dp) :: im_zlogdet ! tmp for imaginary part
     integer :: ierror
     ! working arrays for lapack subroutines
     integer, allocatable :: ipiv(:)
     complex(dp), allocatable, dimension(:,:) :: zmat_tmp

     allocate( zmat_tmp(ndim,ndim) )
     allocate(ipiv(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         write(*,'(a)') 'error in s_logdet_z : can not allocate enough memory'
         stop
     end if
     zmat_tmp = zmat

     ! computes the LU factorization of a general m-by-n matrix, need lapack package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat_tmp, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         write(*,'(a)') 'error in s_logdet_z : error in lapack subroutine zgetrf'
         stop
     end if

     ! calculate determinant
     zlogdet = czero
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zlogdet = zlogdet + log( +zmat_tmp(i,i) )
         else
             zlogdet = zlogdet + log( -zmat_tmp(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop
     im_zlogdet = dimag(zlogdet)
     i=nint(im_zlogdet/2.d0/pi)
     im_zlogdet = im_zlogdet - 2.d0*pi*dble(i)
     zlogdet = dcmplx( dble(zlogdet), im_zlogdet )

     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(zmat_tmp) ) deallocate(zmat_tmp)
     return
  end subroutine s_logdet_z

  subroutine s_eig_he(ldim, ndim, amat, eval, evec)
  !! computes all eigenvalues and eigenvectors of a complex Hermitian matrix
     use constants, only : dp, zero
     implicit none
     integer, intent(in)  :: ldim ! number of rows of the matrix amat
     integer, intent(in)  :: ndim ! number of columns of of the matrix amat
     complex(dp), intent(in)  :: amat(ldim,ndim) ! original complex Hermitian matrix to compute eigenvals and eigenvectors
     real(dp), intent(out)    :: eval(ndim) ! if info = 0, the eigenvalues in ascending order
     complex(dp), intent(out) :: evec(ldim,ndim) ! if info = 0, orthonormal eigenvectors of the matrix

     integer :: istat ! status flag
     integer :: info ! return information from subroutine zheev

     ! the length of the array work and rwork
     ! lwork >= max(1,2*ndim-1), lrwork >= max(1,3*ndim-2)
     integer :: lwork
     integer :: lrwork

     ! workspace array
     real(dp), allocatable    :: rwork(:)
     complex(dp), allocatable :: work(:)

     ! initialize lwork (lrwork)
     lwork = 2 * ndim - 1
     lrwork = 3 * ndim - 2

     ! allocate memory
     allocate(work(lwork),   stat=istat)
     allocate(rwork(lrwork), stat=istat)
     if ( istat /= 0 ) then
         write(*,'(a)') 'error in s_eig_he : can not allocate enough memory'
         stop
     endif ! back if ( istat /= 0 ) block

     ! initialize output arrays
     eval = zero
     evec = amat

     ! call the computational subroutine: zheev
     call ZHEEV('V', 'U', ndim, evec, ldim, eval, work, lwork, rwork, info)

     ! check the status
     if ( info /= 0 ) then
         write(*, '(a)') 'error in s_eig_he : error in lapack subroutine zheev'
         stop
     endif ! back if ( info /= 0 ) block

     ! dealloate memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)
     return
  end subroutine s_eig_he

  subroutine s_svd_zg(m, n, min_mn, amat, umat, svec, vmat)
  !! s_svd_zg: perform the singular values decomposition for a general
  !! complex(dp) m-by-n matrix A, where A = U * SIGMA * conjugate-transpose(V),
  !! return its left vectors, right vectors, and singular values
     use constants, only : dp
     implicit none
     integer, intent(in)        :: m ! number of rows of A matrix
     integer, intent(in)        :: n ! number of columns of A matrix
     integer, intent(in)        :: min_mn ! minimal value of m and n
     complex(dp), intent(inout) :: amat(m,n) ! A matrix
     complex(dp), intent(out)   :: umat(m,min_mn) ! left vectors of svd, U
     real(dp), intent(out)      :: svec(min_mn) ! singular values of svd, SIGMA
     complex(dp), intent(out)   :: vmat(min_mn,n) ! right vectors of svd, conjugate-transpose(V)

     integer :: istat ! status flag
     integer :: info ! return information from zgesvd
     integer :: lwork ! length of work array, lwork >= max(1, 2 * min_mn + max(m,n))
     ! workspace arrays
     complex(dp), allocatable :: work(:)
     real(dp), allocatable :: rwork(:)

     allocate(rwork(5*min_mn), stat=istat) ! allocate memory

     ! lwork query, to get optimal lwork
     lwork = -1
     allocate(work(1),  stat=istat)
     call ZGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, rwork, info)
     lwork = nint(dble(work(1)))
     deallocate(work)

     ! with optimal lwork, do SVD
     allocate(work(lwork),  stat=istat)
     if ( istat /= 0 ) then
         write(*,'(a)') 'error in s_svd_zg : can not allocate enough memory'
         stop
     endif ! back if ( istat /= 0 ) block

     ! call the computational subroutine: zgesvd
     call ZGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, rwork, info)

     ! check the status
     if ( info /= 0 ) then
         write(*,'(a)') 'error in s_svd_zg : error in lapack subroutine zgesvd'
         stop
     endif ! back if ( info /= 0 ) block

     ! deallocate the memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)

     return
  end subroutine s_svd_zg

  subroutine s_svd_zg_logdetU(m, n, min_mn, amat, umat, svec, vmat, zlogdetU)
  !! s_svd_zg_logdetU: perform the singular values decomposition for a general
  !! complex(dp) m-by-n matrix A, where A = U * SIGMA * conjugate-transpose(V),
  !! return its left vectors, right vectors, and singular values
  !! also return the log determinant of Umat by calling s_logdet_z, thus cost another N^3, not efficient
  !! we'd better calculate logdetU during svd, not implement yet.
  !! note the subroutine is only suiteable for m=min(m,n) case.
     use constants, only : dp
 
     implicit none

     integer, intent(in)        :: m ! number of rows of A matrix
     integer, intent(in)        :: n ! number of columns of A matrix
     integer, intent(in)        :: min_mn ! minimal value of m and n
     complex(dp), intent(inout) :: amat(m,n) ! A matrix
     complex(dp), intent(out)   :: umat(m,min_mn) ! left vectors of svd, U
     real(dp), intent(out)      :: svec(min_mn) ! singular values of svd, SIGMA
     complex(dp), intent(out)   :: vmat(min_mn,n) ! right vectors of svd, conjugate-transpose(V)
     complex(dp), intent(out) :: zlogdetU ! log(det(umat))

     integer :: istat ! status flag
     integer :: info ! return information from zgesvd
     integer :: lwork ! length of work array, lwork >= max(1, 2 * min_mn + max(m,n))
     ! workspace arrays
     complex(dp), allocatable :: work(:)
     real(dp), allocatable :: rwork(:)

     ! check m and n
     if ( m > n ) then
         write(*,'(a)') 'error in s_svd_zg_logdetU : m should less equal than n'
         stop
     end if

     ! allocate memory
     allocate(rwork(5*min_mn), stat=istat)

     ! lwork query, to get optimal lwork
     lwork = -1
     allocate(work(1),  stat=istat)
     call ZGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, rwork, info)
     lwork = nint(dble(work(1)))
     deallocate(work)

     ! with optimal lwork, do SVD
     allocate(work(lwork),  stat=istat)
     if ( istat /= 0 ) then
         write(*,'(a)') 'error in s_svd_zg_logdetU : can not allocate enough memory'
         stop
     endif ! back if ( istat /= 0 ) block

     ! call the computational subroutine: zgesvd
     call ZGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, rwork, info)

     ! check the status
     if ( info /= 0 ) then
         write(*,'(a)') 'error in s_svd_zg_logdetU : error in lapack subroutine zgesvd'
         stop
     endif ! back if ( info /= 0 ) block

     ! deallocate the memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)

     call s_logdet_z(m, umat, zlogdetU)

     return
  end subroutine s_svd_zg_logdetU

  subroutine s_compare_max_z( ndim, Amat, Bmat, max_diff )
  !! maximum elements of |Amat - Bmat|, complex number version
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only : dp, zero
    implicit none

    integer, intent(in) :: ndim
    complex(dp), dimension(ndim,ndim), intent(in) :: Amat, Bmat
    real(dp), intent(out) :: max_diff
    
    ! local
    integer :: i, j
    real(dp) :: max_tmp

    max_diff = zero
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i, max_tmp  )
!$OMP DO REDUCTION ( MAX : max_diff )
    do j = 1, ndim
        do i = 1, ndim
            max_tmp = abs( Amat(i,j) - Bmat(i,j) )
            if( max_tmp .gt. max_diff ) then
                max_diff = max_tmp
            end if
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
  end subroutine s_compare_max_z

  subroutine s_compare_max_d( ndim, Amat, Bmat, max_diff )
  !! maximum elements of |Amat - Bmat|, real number version
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only : dp, zero
    implicit none
    integer, intent(in) :: ndim
    real(dp), dimension(ndim,ndim), intent(in) :: Amat, Bmat
    real(dp), intent(out) :: max_diff
    
    ! local
    integer :: i, j
    real(dp) :: max_tmp

    max_diff = zero
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i, max_tmp  )
!$OMP DO REDUCTION ( MAX : max_diff )
    do j = 1, ndim
        do i = 1, ndim
            max_tmp = dabs( Amat(i,j) - Bmat(i,j) )
            if( max_tmp .gt. max_diff ) then
                max_diff = max_tmp
            end if
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
  end subroutine s_compare_max_d

  subroutine s_z_x_diag_d( ndim, Amat, dvec, Bmat )
  !! matrix product of a complex matrix Amat and a real diagonal matrix (diagonal elements are stored in a vector dvec)
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only: dp, zero
    implicit none
    integer, intent(in) :: ndim
    complex(dp), dimension(ndim,ndim), intent(in) :: Amat
    real(dp), dimension(ndim), intent(in) :: dvec
    complex(dp), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
!$OMP PARALLEL &
!$OMP PRIVATE ( i, j )
!$OMP DO
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dcmplx( dvec(i), zero ) * Amat(j,i)
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL

  end subroutine s_z_x_diag_d

  subroutine s_diag_d_x_z( ndim, dvec, Amat, Bmat )
  !! matrix product of a real diagonal matrix (diagonal elements are stored in a vector dvec) and a complex matrix Amat
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only: dp, zero
    implicit none
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in) :: dvec
    complex(dp), dimension(ndim,ndim), intent(in) :: Amat
    complex(dp), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
!$OMP PARALLEL &
!$OMP PRIVATE ( i, j )
!$OMP DO
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dcmplx( dvec(j), zero ) * Amat(j,i)
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL

  end subroutine s_diag_d_x_z

  subroutine s_diag_dvd( ndim, dvecr, Amat, dvecl, Bmat )
  !! matrix product of a real diagonal matrix (diagonal elements are stored in a vector dvecr), a complex matrix Amat
  !! and another real diagonal matrix (diagonal elements are stored in a vector dvecr)
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only: dp, zero
    implicit none
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in) :: dvecr, dvecl
    complex(dp), dimension(ndim,ndim), intent(in) :: Amat
    complex(dp), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
!$OMP PARALLEL &
!$OMP PRIVATE ( i, j )
!$OMP DO
    do i = 1, ndim
        do j = 1, ndim
            Bmat(j,i) = dcmplx( dvecr(j)*dvecl(i), zero ) * Amat(j,i)
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL

  end subroutine s_diag_dvd

  subroutine s_v_invd_u( ndim, vmat, dvec, umat, zmat )
  !! matrix product of a complex matrix vmat, inverse of a real diagonal matrix (diagonal elements are stored in a vector dvecr)
  !! and another complex matrix umat
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only: dp, czero, cone
    implicit none
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in) :: dvec
    complex(dp), dimension(ndim,ndim), intent(in) :: vmat, umat
    complex(dp), dimension(ndim,ndim), intent(out) :: zmat
    
    ! local
    integer :: i, j
    complex(dp), allocatable, dimension(:,:) :: v_invd
    allocate( v_invd(ndim,ndim) )
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
    do j = 1, ndim
        do i = 1, ndim
            v_invd(i,j) = vmat(i,j) / dvec(j)
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    call zgemm('n','n',ndim,ndim,ndim,cone,v_invd,ndim,umat,ndim,czero,zmat,ndim)
    deallocate( v_invd )
  end subroutine s_v_invd_u

  subroutine s_v_d_u( ndim, vmat, dvec, umat, zmat )
  !! matrix product of a complex matrix vmat, a real diagonal matrix (diagonal elements are stored in a vector dvecr)
  !! and another complex matrix umat
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only: dp, czero, cone
    implicit none
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in) :: dvec
    complex(dp), dimension(ndim,ndim), intent(in) :: vmat, umat
    complex(dp), dimension(ndim,ndim), intent(out) :: zmat
    
    ! local
    integer :: i, j
    complex(dp), allocatable, dimension(:,:) :: v_d
    allocate( v_d(ndim,ndim) )
!$OMP PARALLEL &
!$OMP PRIVATE ( j, i )
!$OMP DO
    do j = 1, ndim
        do i = 1, ndim
            v_d(i,j) = vmat(i,j) * dvec(j)
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    call zgemm('n','n',ndim,ndim,ndim,cone,v_d,ndim,umat,ndim,czero,zmat,ndim)
    deallocate( v_d )
  end subroutine s_v_d_u

  subroutine s_dvec_min_max(ndim,dvec,dmax,dmin)
  !! decomposition of a real diagonal matrix (diagonal elements are stored in a vector dvec)
  !! as product of two real diagonal matrices (diagonal elements are stored in vectors dmax and dmin), with
  !! one contains diagonal elements greater than one (less than one elements are replaced with one)
  !! the other contains diagonal elements less than one (greater than one elements are replaced with one)
    use constants, only: dp, one
    implicit none
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in) :: dvec
    real(dp), dimension(ndim), intent(out) :: dmax, dmin

    ! local
    integer :: i

    do i = 1, ndim
        if( abs(dvec(i)) .gt. one ) then
            dmax(i) = abs( dvec(i) )
            dmin(i) = sign(one,dvec(i))
        else
            dmax(i) = 1
            dmin(i) = dvec(i)
        end if
    end do
  end subroutine s_dvec_min_max

  subroutine s_zmcpt(rdim,cdim,amat,jpvt)
  ! make matrix amat to column norm decreasing
  ! amat = amat*P
  ! output amat and PT
    use constants, only: dp
    implicit none
    integer, intent(in) :: rdim, cdim
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    integer, dimension(cdim), intent(out) :: jpvt

    ! external
    integer, external :: idamax
    real(dp), external :: dznrm2

    ! local
    integer :: k, pvt, itemp
    real(dp), allocatable, dimension(:) :: vn1

    allocate( vn1(cdim) )

    do k = 1, cdim
        vn1(k) = dznrm2(rdim,amat(1,k),1)
        jpvt(k) = k
    end do

    do k = 1, cdim
        pvt = (k-1) + idamax(cdim-k+1, vn1(k), 1)
        if( pvt .ne. k ) then
            call zswap( rdim, amat(1,pvt), 1, amat(1,k), 1 )
            itemp = jpvt( pvt )
            jpvt( pvt ) = jpvt( k )
            jpvt( k ) = itemp
            vn1(pvt) = vn1(k)
        end if
    end do
    !write(*,*) ' after s_zmcpt, jpvt = '
    !write(*,*) jpvt(:)
    deallocate(vn1)
  end subroutine s_zmcpt

  subroutine s_zmrpt(rdim,cdim,amat,jpvt)
  ! make matrix amat to row norm decreasing
  ! amat = P*amat
  ! output amat and PT
    use constants, only: dp
    implicit none
    integer, intent(in) :: rdim, cdim
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    integer, dimension(rdim), intent(out) :: jpvt

    ! external
    integer, external :: idamax
    real(dp), external :: dznrm2

    ! local
    integer :: k, pvt, itemp
    real(dp), allocatable, dimension(:) :: vn1


    allocate( vn1(rdim) )

    do k = 1, rdim
        vn1(k) = dznrm2(cdim,amat(k,1),rdim)
        jpvt(k) = k
    end do

    do k = 1, rdim
        pvt = (k-1) + idamax(rdim-k+1, vn1(k), 1)
        if( pvt .ne. k ) then
            call zswap( cdim, amat(pvt,1), rdim, amat(k,1), rdim )
            itemp = jpvt( pvt )
            jpvt( pvt ) = jpvt( k )
            jpvt( k ) = itemp
            vn1(pvt) = vn1(k)
        end if
    end do
    deallocate(vn1)
  end subroutine s_zmrpt

  subroutine s_zpm(rdim,cdim,jpv,amat)
  ! amat = P * amat, where P is permutation
  !       note P = (a1,a2,a3,....,an) and is defined with the right order when perform column permutation, namely, A*P.
  !                ( 1, 2, 3,...., n)
  !       if we want to use it to perform row permutation, we should first reverse the permutation order, which is
  !       equivalent to do the transpose of P, P^T = ( 1, 2, 3,......, n) = ( b1,b2,b3,......,bn)
  !                                                  (a1,a2,a3,......,an)   (  1, 2, 3,......, n)
  !       then perform P^T on amat with row permutation
    use constants, only: dp
    implicit none
    integer, intent(in) :: rdim, cdim
    integer, dimension(rdim), intent(inout) :: jpv
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat

    ! local
    integer :: k, pv, itemp, npt
    integer, dimension(:), allocatable :: jpvt

    allocate(jpvt(rdim))

    do k = 1, rdim
        jpvt(k) = k
    end do

    !write(*,*) ' in zpm, input jpv = '
    !write(*,*) jpv(:)
    do
        npt = 0
        do k = 1, rdim
            pv = jpv(k)
            if( pv .ne. k ) then 
                npt = npt + 1
                ! jpv
                itemp = jpv( pv )
                jpv( pv ) = jpv( k )
                jpv( k ) = itemp
                ! jpvt
                itemp = jpvt( pv )
                jpvt( pv ) = jpvt( k )
                jpvt( k ) = itemp
            end if
        end do
        if( npt.eq.0) exit
    end do
    !write(*,*) ' in zpm, jpv^T = '
    !write(*,*) jpvt(:)

    do
        npt = 0
        do k = 1, rdim
            pv = jpvt(k)
            if( pv .ne. k ) then
                npt = npt + 1
                call zswap( cdim, amat(pv,1), rdim, amat(k,1), rdim )
                itemp = jpvt( pv )
                jpvt( pv) = jpvt( k )
                jpvt( k ) = itemp
            end if
        end do
        if( npt.eq.0) exit
    end do
    !write(*,*) ' after zpm, jpvt = '
    !write(*,*) jpvt(:)

    deallocate(jpvt)

  end subroutine s_zpm

  subroutine s_pp_R(ndim,jpv1,jpv2)
  !! jpv1, jpv2 are column pivoting permutation
  !! perform * jpv1 * jpv2, results store in jpv1
  !! algorithm: push jpv1^T into jpv2, results store in jpv1
    implicit none
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(inout) :: jpv1
    integer, dimension(ndim), intent(in) :: jpv2

    ! local 
    integer :: k, pv, itemp, npt
    integer, dimension(:), allocatable :: jpvt

    allocate(jpvt(ndim) )

    !write(*,'(a,16i6)') ' in s_pp_R, input jpv1 = ', jpv1(:)
    !write(*,'(a,16i6)') '            input jpv2 = ', jpv2(:)
    !! transpose on jpv1, store in jpvt
    do k = 1, ndim
        jpvt(k) = k
    end do
    do
        npt = 0
        do k = 1, ndim
            pv = jpv1(k)
            if( pv .ne. k ) then 
                npt = npt + 1
                ! jpv1
                itemp = jpv1( pv )
                jpv1( pv ) = jpv1( k )
                jpv1( k ) = itemp
                ! jpvt
                itemp = jpvt( pv )
                jpvt( pv ) = jpvt( k )
                jpvt( k ) = itemp
            end if
        end do
        if( npt.eq.0) exit
    end do
    !! push jpv1^T into jpv2
    jpv1(:) = jpv2(:) ! copy jpv2 to jpv1, after push jpv1^T into it, directly return jpv1
    do
        npt = 0
        do k = 1, ndim
            pv = jpvt(k)
            if( pv .ne. k ) then
                npt = npt + 1
                ! jpv1
                itemp = jpv1( pv )
                jpv1( pv ) = jpv1( k )
                jpv1( k ) = itemp
                ! jpvt
                itemp = jpvt( pv )
                jpvt( pv) = jpvt( k )
                jpvt( k ) = itemp
            end if
        end do
        if( npt.eq.0) exit
    end do
    !write(*,'(a,16i6)') '           jpv1 * jpv2 = ', jpv1(:)
    deallocate(jpvt)
  end subroutine s_pp_R

  subroutine s_pp_L(ndim,jpv1,jpv2)
  !! jpv1, jpv2 are row pivoting permutation
  !! perform jpv1 * jpv2 * , results store in jpv2
  !! algorithm: push jpv2^T into jpv1, results sotre in jpv2
    implicit none
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(inout) :: jpv2
    integer, dimension(ndim), intent(in) :: jpv1

    ! local 
    integer :: k, pv, itemp, npt
    integer, dimension(:), allocatable :: jpvt

    allocate(jpvt(ndim))

    !write(*,'(a,4i6)') ' in s_pp_L, input jpv1 = ', jpv1(:)
    !write(*,'(a,4i6)') '            input jpv2 = ', jpv2(:)

    ! perform jpv2^T
    do k = 1, ndim
        jpvt(k) = k
    end do
    do
        npt = 0
        do k = 1, ndim
            pv = jpv2(k)
            if( pv .ne. k ) then
                npt = npt + 1
                ! jpvt
                itemp = jpvt( pv )
                jpvt( pv ) = jpvt( k )
                jpvt( k ) = itemp
                ! jpv2
                itemp = jpv2( pv )
                jpv2( pv) = jpv2( k )
                jpv2( k ) = itemp
            end if
        end do
        if( npt.eq.0) exit
    end do
    !! push jpv2^T into jpv1
    jpv2(:) = jpv1(:) ! copy jpv1 to jpv2, after push jpv2^T into it, directly return jpv2
    do
        npt = 0
        do k = 1, ndim
            pv = jpvt(k)
            if( pv .ne. k ) then
                npt = npt + 1
                ! jpv2
                itemp = jpv2( pv )
                jpv2( pv ) = jpv2( k )
                jpv2( k ) = itemp
                ! jpvt
                itemp = jpvt( pv )
                jpvt( pv) = jpvt( k )
                jpvt( k ) = itemp
            end if
        end do
        if( npt.eq.0) exit
    end do
    !write(*,'(a,4i6)') '           jpv1 * jpv2 = ', jpv2(:)
    deallocate(jpvt)
  end subroutine s_pp_L

  subroutine s_zmp(rdim,cdim,amat,jpv)
  ! amat = amat * P, where P is permutation
  ! note after the permutation, P will become (1234....)
    use constants, only: dp
    implicit none
    integer, intent(in) :: rdim, cdim
    integer, dimension(cdim), intent(inout) :: jpv
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat

    ! local
    integer :: k, pv, itemp, npt

    do 
        npt = 0
        do k = 1, cdim
            pv = jpv(k)
            if( pv .ne. k ) then
                npt = npt + 1
                call zswap( rdim, amat(1,pv), 1, amat(1,k), 1)
                itemp = jpv( pv )
                jpv( pv ) = jpv( k )
                jpv( k ) = itemp
            end if
        end do
        if( npt.eq.0) exit
    end do
    !write(*,*) ' after zmp, jpv = '
    !write(*,*) jpv(:)
  end subroutine s_zmp

  subroutine s_dppd_z(ndim, dvecr, jpvrin, jpvlin, dvecl, amat)
  !! amat = (dvecr*jpvr)*(jpvl*dvecl)
    use constants, only: dp, czero
    implicit none
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(in) :: jpvrin, jpvlin
    real(dp), dimension(ndim), intent(in) :: dvecr, dvecl
    complex(dp), dimension(ndim,ndim), intent(out) :: amat

    ! local
    integer :: i, j, k
    integer, dimension(:), allocatable :: jpvr, jpvl
    allocate(jpvr(ndim),jpvl(ndim))
    jpvr(:) = jpvrin(:)
    jpvl(:) = jpvlin(:)
    amat(:,:) = czero
    do i = 1, ndim
        k = jpvr(i)
        do j = 1, ndim
            if( jpvl(j) .eq. k ) then
                amat(i,j) = dcmplx( dvecr(i)*dvecl(j), 0.d0 )
            end if
        end do
    end do
    deallocate(jpvl,jpvr)
  end subroutine s_dppd_z

  subroutine s_zgeQRPT(rdim, cdim, amat, qmat, rmat, jpvt)
    use constants, only : dp, zero, czero, one
    implicit none
    integer, intent(in) :: rdim, cdim
    integer, dimension(cdim), intent(out) :: jpvt
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    complex(dp), dimension(min(rdim,cdim),cdim), intent(out) :: rmat
    complex(dp), dimension(min(rdim,cdim), min(rdim,cdim) ), intent(out) :: qmat

    ! local
    integer  :: i, ierror, lwork, mindim
    real(dp), dimension(:), allocatable :: rwork
    complex(dp), dimension(:), allocatable :: work
    complex(dp), dimension(:), allocatable :: tau

    mindim = min(rdim,cdim)

    !! perform QR factorization
    lwork=-1
    allocate( tau(mindim) )
    allocate( rwork(2*cdim) )
    allocate( work(1) )
    ! perform lwork query
    jpvt(:) = 0
    call zgeqp3(rdim, cdim, amat, rdim, jpvt, tau, work, lwork, rwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zgeqp3, lwork =', lwork
    deallocate(work)

    ! perform QR factorization
    allocate(work(lwork),stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT, the place above zgeqp3 : can not allocate enough memory'
        stop
    endif
    jpvt(:) = 0
    call zgeqp3(rdim, cdim, amat, rdim, jpvt, tau, work, lwork, rwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT : error in lapack subroutine zgeqp3'
        stop
    endif
    deallocate(work)
    deallocate(rwork)

    !! get R
    rmat(:,:) = czero
    do i = 1, mindim
        rmat(1:i,i) = amat(1:i,i)
    end do
    do i = mindim+1, cdim
        rmat(:,i) = amat(:,i)
    end do

    !! get Q
    ! get reflectors, stored in qmat
    call s_identity_z(mindim,qmat)
    do i = 1, mindim-1
        qmat(i+1:mindim,i) = amat(i+1:mindim,i)
    end do
    lwork = -1
    allocate(work(1))
    ! perform lwork query
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in dungqr, lwork =', lwork
    deallocate(work)

    ! get Q
    allocate(work(lwork), stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT, the place above zungqr : can not allocate enough memory'
        stop
    endif
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT : error in lapack subroutine zungqr'
        stop
    endif
    deallocate(work)
    deallocate(tau)

  end subroutine s_zgeQRPT

  subroutine s_zgeQRPT_logdetQ(rdim, cdim, amat, qmat, rmat, jpvt, zlogdet)
  ! perform QR decomposition with pivoting, and calculate log det of Q
    use constants, only : dp, zero, czero, one, cone
    implicit none
    integer, intent(in) :: rdim, cdim
    integer, dimension(cdim), intent(out) :: jpvt
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    complex(dp), dimension(min(rdim,cdim),cdim), intent(out) :: rmat
    complex(dp), dimension(min(rdim,cdim), min(rdim,cdim) ), intent(out) :: qmat
    complex(dp),intent(out) :: zlogdet

    ! local
    integer  :: i, ierror, lwork, mindim
    real(dp), dimension(:), allocatable :: rwork
    complex(dp), dimension(:), allocatable :: work
    complex(dp), dimension(:), allocatable :: tau

    mindim = min(rdim,cdim)

    !! perform QR factorization
    lwork=-1
    allocate( tau(mindim) )
    allocate( rwork(2*cdim) )
    allocate( work(1) )
    ! perform lwork query
    jpvt(:) = 0
    call zgeqp3(rdim, cdim, amat, rdim, jpvt, tau, work, lwork, rwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zgeqp3, lwork =', lwork
    deallocate(work)

    ! perform QR factorization
    allocate(work(lwork),stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT_logdetQ, the place above zgeqp3 : can not allocate enough memory'
        stop
    endif
    jpvt(:) = 0
    call zgeqp3(rdim, cdim, amat, rdim, jpvt, tau, work, lwork, rwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT_logdetQ : error in lapack subroutine zgeqp3'
        stop
    endif
    deallocate(work)
    deallocate(rwork)

    !! get determinant of Q
    zlogdet = czero
    do i = 1, mindim
        if(abs(tau(i)).ne.0.d0) then
            zlogdet = zlogdet + log( cone - tau(i)*dcmplx(2.d0*dble(tau(i)), 0.d0 ) / (tau(i)*dconjg(tau(i))) )
        end if
    end do

    !! get R
    rmat(:,:) = czero
    do i = 1, mindim
        rmat(1:i,i) = amat(1:i,i)
    end do
    do i = mindim+1, cdim
        rmat(:,i) = amat(:,i)
    end do

    !! get Q
    ! get reflectors, stored in qmat
    call s_identity_z(mindim,qmat)
    do i = 1, mindim-1
        qmat(i+1:mindim,i) = amat(i+1:mindim,i)
    end do
    lwork = -1
    allocate(work(1))
    ! perform lwork query
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in dungqr, lwork =', lwork
    deallocate(work)

    ! get Q
    allocate(work(lwork), stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT_logdetQ, the place above zungqr : can not allocate enough memory'
        stop
    endif
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQRPT_logdetQ : error in lapack subroutine zungqr'
        stop
    endif
    deallocate(work)
    deallocate(tau)

  end subroutine s_zgeQRPT_logdetQ

  subroutine s_invdiag_d_x_zr( ndim, dvec, Amat, Bmat )
#IFDEF _OPENMP
    USE OMP_LIB
#ENDIF
    use constants, only: dp, czero, zero
    implicit none
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in) :: dvec
    complex(dp), dimension(ndim,ndim), intent(in) :: Amat
    complex(dp), dimension(ndim,ndim), intent(out) :: Bmat
    
    ! local
    integer :: i, j
    Bmat=czero
!$OMP PARALLEL &
!$OMP PRIVATE ( i, j )
!$OMP DO SCHEDULE(DYNAMIC)
    do i = 1, ndim
        do j = 1, i
            Bmat(j,i) =  Amat(j,i) / dcmplx( dvec(j), zero )
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
  end subroutine s_invdiag_d_x_zr

  subroutine s_zgeQR(rdim, cdim, amat, qmat, rmat )
  ! perform QR decomposition
    use constants, only : dp, zero, czero, one
    implicit none
    integer, intent(in) :: rdim, cdim
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    complex(dp), dimension(min(rdim,cdim),cdim), intent(out) :: rmat
    complex(dp), dimension(min(rdim,cdim), min(rdim,cdim) ), intent(out) :: qmat

    ! local
    integer  :: i, ierror, lwork, mindim
    complex(dp), dimension(:), allocatable :: work
    complex(dp), dimension(:), allocatable :: tau

    mindim = min(rdim,cdim)

    !! perform QR factorization
    lwork=-1
    allocate( tau(mindim) )
    allocate( work(1) )
    ! perform lwork query
    call zgeqrf(rdim, cdim, amat, rdim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zgeqrf, lwork =', lwork
    deallocate(work)

    ! perform QR factorization
    allocate(work(lwork),stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR, the place above zgeqrf : can not allocate enough memory'
        stop
    endif
    call zgeqrf(rdim, cdim, amat, rdim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR : error in lapack subroutine zgeqrf'
        stop
    endif
    deallocate(work)

    !! get R
    rmat(:,:) = czero
    do i = 1, mindim
        rmat(1:i,i) = amat(1:i,i)
    end do
    do i = mindim+1, cdim
        rmat(:,i) = amat(:,i)
    end do

    !! get Q
    ! get reflectors, stored in qmat
    call s_identity_z(mindim,qmat)
    do i = 1, mindim-1
        qmat(i+1:mindim,i) = amat(i+1:mindim,i)
    end do
    lwork = -1
    allocate(work(1))
    ! perform lwork query
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zungqr, lwork =', lwork
    deallocate(work)

    ! get Q
    allocate(work(lwork), stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR, the place above zungqr : can not allocate enough memory'
        stop
    endif
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR : error in lapack subroutine zungqr'
        stop
    endif
    deallocate(work)
    deallocate(tau)

  end subroutine s_zgeQR

  subroutine s_zgeQR_logdetQ(rdim, cdim, amat, qmat, rmat, zlogdet )
    use constants, only : dp, zero, czero, one, cone
    implicit none
    integer, intent(in) :: rdim, cdim
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    complex(dp), dimension(min(rdim,cdim),cdim), intent(out) :: rmat
    complex(dp), dimension(min(rdim,cdim), min(rdim,cdim) ), intent(out) :: qmat
    complex(dp), intent(out) :: zlogdet

    ! local
    integer  :: i, ierror, lwork, mindim
    complex(dp), dimension(:), allocatable :: work
    complex(dp), dimension(:), allocatable :: tau

    mindim = min(rdim,cdim)

    !! perform QR factorization
    lwork=-1
    allocate( tau(mindim) )
    allocate( work(1) )
    ! perform lwork query
    call zgeqrf(rdim, cdim, amat, rdim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zgeqrf, lwork =', lwork
    deallocate(work)

    ! perform QR factorization
    allocate(work(lwork),stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR_logdetQ, the place above zgeqrf : can not allocate enough memory'
        stop
    endif
    call zgeqrf(rdim, cdim, amat, rdim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR_logdetQ : error in lapack subroutine zgeqrf'
        stop
    endif
    deallocate(work)

    !! get determinant of Q
    zlogdet = czero
    do i = 1, mindim
        if(abs(tau(i)).ne.0.d0) then
            zlogdet = zlogdet + log( cone - tau(i)*dcmplx(2.d0*dble(tau(i)), 0.d0 ) / (tau(i)*dconjg(tau(i))) )
        end if
    end do

    !! get R
    rmat(:,:) = czero
    do i = 1, mindim
        rmat(1:i,i) = amat(1:i,i)
    end do
    do i = mindim+1, cdim
        rmat(:,i) = amat(:,i)
    end do

    !! get Q
    ! get reflectors, stored in qmat
    call s_identity_z(mindim,qmat)
    do i = 1, mindim-1
        qmat(i+1:mindim,i) = amat(i+1:mindim,i)
    end do
    lwork = -1
    allocate(work(1))
    ! perform lwork query
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zungqr, lwork =', lwork
    deallocate(work)

    ! get Q
    allocate(work(lwork), stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR_logdetQ, the place above zungqr : can not allocate enough memory'
        stop
    endif
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_zgeQR_logdetQ : error in lapack subroutine zungqr'
        stop
    endif
    deallocate(work)
    deallocate(tau)

  end subroutine s_zgeQR_logdetQ

  subroutine s_invqr_z(ndim, amat )
  ! inverse of a complex matrix, using QR decomposition
    use constants, only : dp, zero, czero, one, cone
    implicit none
    integer, intent(in) :: ndim
    complex(dp), dimension(ndim,ndim), intent(inout) :: amat

    ! local
    integer  :: i, ierror, lwork
    complex(dp), dimension(:), allocatable :: work
    complex(dp), dimension(:), allocatable :: tau
    complex(dp), dimension(:,:), allocatable :: qmat
    complex(dp), dimension(:,:), allocatable :: rmat

    allocate( qmat(ndim,ndim) )
    allocate( rmat(ndim,ndim) )

    !! perform QR factorization
    lwork=-1
    allocate( tau(ndim) )
    allocate( work(1) )
    ! perform lwork query
    call zgeqrf(ndim, ndim, amat, ndim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zgeqrf, lwork =', lwork
    deallocate(work)

    ! perform QR factorization
    allocate(work(lwork),stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_invqr_z, the place above zgeqrf : can not allocate enough memory'
        stop
    endif
    call zgeqrf(ndim, ndim, amat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_invqr_z : error in lapack subroutine zgeqrf'
        stop
    endif
    deallocate(work)

    !! get Q
    ! get reflectors, stored in qmat
    call s_identity_z(ndim,qmat)
    do i = 1, ndim-1
        qmat(i+1:ndim,i) = amat(i+1:ndim,i)
    end do
    lwork = -1
    allocate(work(1))
    ! perform lwork query
    call zungqr(ndim, ndim, ndim, qmat, ndim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zungqr, lwork =', lwork
    deallocate(work)

    ! get Q
    allocate(work(lwork), stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_invqr_z, the place above zungqr : can not allocate enough memory'
        stop
    endif
    call zungqr(ndim, ndim, ndim, qmat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_invqr_z : error in lapack subroutine zungqr'
        stop
    endif
    deallocate(work)
    deallocate(tau)

    !! get inverse of R
    call ztrtri('U','N',ndim,amat,ndim,ierror)
    rmat(:,:) = czero
    do i = 1, ndim
        rmat(1:i,i) = amat(1:i,i)
    end do

    !! get A^-1 = R^-1*Q^H
    call zgemm('n','c',ndim,ndim,ndim,cone,rmat,ndim,qmat,ndim,czero,amat,ndim)

    deallocate(rmat)
    deallocate(qmat)
  end subroutine s_invqr_z

  subroutine s_inv_logdet_qr_z(ndim, amat, zlogdet )
    use constants, only : dp, zero, czero, one, cone, pi
    implicit none
    integer, intent(in) :: ndim
    complex(dp), dimension(ndim,ndim), intent(inout) :: amat
    complex(dp), intent(out) :: zlogdet

    ! local
    integer  :: i, ierror, lwork
    real(dp) :: im_zlogdet
    complex(dp) :: zlogtmp
    complex(dp), dimension(:), allocatable :: work
    complex(dp), dimension(:), allocatable :: tau
    complex(dp), dimension(:,:), allocatable :: qmat
    complex(dp), dimension(:,:), allocatable :: rmat

    allocate( qmat(ndim,ndim) )
    allocate( rmat(ndim,ndim) )

    !! perform QR factorization
    lwork=-1
    allocate( tau(ndim) )
    allocate( work(1) )
    ! perform lwork query
    call zgeqrf(ndim, ndim, amat, ndim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zgeqrf, lwork =', lwork
    deallocate(work)

    ! perform QR factorization
    allocate(work(lwork),stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_logdet_qr_z, the place above zgeqrf : can not allocate enough memory'
        stop
    endif
    call zgeqrf(ndim, ndim, amat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_logdet_qr_z : error in lapack subroutine zgeqrf'
        stop
    endif
    deallocate(work)

    !! get Q
    ! get reflectors, stored in qmat
    call s_identity_z(ndim,qmat)
    do i = 1, ndim-1
        qmat(i+1:ndim,i) = amat(i+1:ndim,i)
    end do
    lwork = -1
    allocate(work(1))
    ! perform lwork query
    call zungqr(ndim, ndim, ndim, qmat, ndim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zungqr, lwork =', lwork
    deallocate(work)

    !! get determinant of Q
    zlogtmp = czero
    do i = 1, ndim
        if(abs(tau(i)).ne.0.d0) then
            zlogtmp = zlogtmp + log( cone - tau(i)*dcmplx(2.d0*dble(tau(i)), 0.d0 ) / (tau(i)*dconjg(tau(i))) )
        end if
    end do

    ! get Q
    allocate(work(lwork), stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_logdet_qr_z, the place above zungqr : can not allocate enough memory'
        stop
    endif
    call zungqr(ndim, ndim, ndim, qmat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_logdet_qr_z : error in lapack subroutine zungqr'
        stop
    endif
    deallocate(work)
    deallocate(tau)

    zlogdet = czero
    !! get determinant of R
    do i = 1, ndim
        zlogdet = zlogdet + log(amat(i,i))
    end do

    !! get total determinant
    zlogdet = zlogdet + zlogtmp
    im_zlogdet = dimag(zlogdet)
    i=nint(im_zlogdet/2.d0/pi)
    im_zlogdet = im_zlogdet - 2.d0*pi*dble(i)
    zlogdet = dcmplx( dble(zlogdet), im_zlogdet )

    !! get inverse of R
    call ztrtri('U','N',ndim,amat,ndim,ierror)
    rmat(:,:) = czero
    do i = 1, ndim
        rmat(1:i,i) = amat(1:i,i)
    end do

    !! get A^-1 = R^-1*Q^H
    call zgemm('n','c',ndim,ndim,ndim,cone,rmat,ndim,qmat,ndim,czero,amat,ndim)

    deallocate(rmat)
    deallocate(qmat)
  end subroutine s_inv_logdet_qr_z

  subroutine s_inv_det_qr_z(ndim, amat, zdet )
    use constants, only : dp, zero, czero, one, cone, pi
    implicit none
    integer, intent(in) :: ndim
    complex(dp), dimension(ndim,ndim), intent(inout) :: amat
    complex(dp), intent(out) :: zdet

    ! local
    integer  :: i, ierror, lwork
    complex(dp) :: ztmp
    complex(dp), dimension(:), allocatable :: work
    complex(dp), dimension(:), allocatable :: tau
    complex(dp), dimension(:,:), allocatable :: qmat
    complex(dp), dimension(:,:), allocatable :: rmat

    allocate( qmat(ndim,ndim) )
    allocate( rmat(ndim,ndim) )

    !! perform QR factorization
    lwork=-1
    allocate( tau(ndim) )
    allocate( work(1) )
    ! perform lwork query
    call zgeqrf(ndim, ndim, amat, ndim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zgeqrf, lwork =', lwork
    deallocate(work)

    ! perform QR factorization
    allocate(work(lwork),stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_det_qr_z, the place above zgeqrf : can not allocate enough memory'
        stop
    endif
    call zgeqrf(ndim, ndim, amat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_det_qr_z : error in lapack subroutine zgeqrf'
        stop
    endif
    deallocate(work)

    !! get Q
    ! get reflectors, stored in qmat
    call s_identity_z(ndim,qmat)
    do i = 1, ndim-1
        qmat(i+1:ndim,i) = amat(i+1:ndim,i)
    end do
    lwork = -1
    allocate(work(1))
    ! perform lwork query
    call zungqr(ndim, ndim, ndim, qmat, ndim, tau, work, lwork, ierror)
    lwork = nint(dble(work(1)))
    !write(*,*) 'in zungqr, lwork =', lwork
    deallocate(work)

    !! get determinant of Q
    ztmp = cone
    do i = 1, ndim
        if(abs(tau(i)).ne.0.d0) then
            ztmp = ztmp * ( cone - tau(i)*dcmplx(2.d0*dble(tau(i)), 0.d0 ) / (tau(i)*dconjg(tau(i))) )
        end if
    end do

    ! get Q
    allocate(work(lwork), stat=ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_det_qr_z, the place above zungqr : can not allocate enough memory'
        stop
    endif
    call zungqr(ndim, ndim, ndim, qmat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        write(*,'(a)') 'error in s_inv_det_qr_z : error in lapack subroutine zungqr'
        stop
    endif
    deallocate(work)
    deallocate(tau)

    zdet = cone
    !! get determinant of R
    do i = 1, ndim
        zdet = zdet * amat(i,i)
    end do

    !! get total determinant
    zdet = zdet * ztmp

    !! get inverse of R
    call ztrtri('U','N',ndim,amat,ndim,ierror)
    rmat(:,:) = czero
    do i = 1, ndim
        rmat(1:i,i) = amat(1:i,i)
    end do

    !! get A^-1 = R^-1*Q^H
    call zgemm('n','c',ndim,ndim,ndim,cone,rmat,ndim,qmat,ndim,czero,amat,ndim)

    deallocate(rmat)
    deallocate(qmat)
  end subroutine s_inv_det_qr_z

  subroutine s_invlu_z(ndim, zmat)
  !! inverse of a complex matrix, using LU decomposition
     use constants, only : dp, czero, cone

     implicit none
     integer, intent(in)        :: ndim ! dimension of zmat matrix
     ! object matrix, on entry, it contains the original matrix, on exit,
     ! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

     integer     :: ierror ! error flag
     integer :: i

     ! working arrays for lapack subroutines
     integer, allocatable     :: ipiv(:)
     complex(dp), allocatable :: work(:)

     complex(dp), allocatable, dimension(:,:) :: umat, lmat

     integer :: k, pv

     allocate(umat(ndim,ndim))
     allocate(lmat(ndim,ndim))

     ! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     allocate(work(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         write(*,'(a)') 'error in s_invlu_z : can not allocate enough memory'
         stop
     end if

     ! computes the LU factorization of a general m-by-n matrix, need lapack
     ! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         write(*,'(a)') 'error in s_invlu_z : error in lapack subroutine zgetrf'
         stop
     end if

     !! get inverse of U
     call ztrtri('U','N',ndim,zmat,ndim,ierror)
     umat(:,:) = czero
     do i = 1, ndim
         umat(1:i,i) = zmat(1:i,i)
     end do

     !! get inverse of L
     call ztrtri('L','U',ndim,zmat,ndim,ierror)
     call s_identity_z(ndim,lmat)
     do i = 1, ndim-1
         lmat(i+1:ndim,i) = zmat(i+1:ndim,i)
     end do

     !! L*PT, this part is get from lapack: zgetri
     do k = ndim-1, 1, -1
         pv = ipiv(k)
         if( pv .ne. k ) then
             call zswap( ndim, lmat(1,k), 1, lmat(1,pv), 1)
         end if
     end do

     !! A^-1 = U^-1 * L^-1 * PT
     call zgemm('n','n',ndim,ndim,ndim,cone,umat,ndim,lmat,ndim,czero,zmat,ndim)

     ! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)
     deallocate(lmat)
     deallocate(umat)

     return
  end subroutine s_invlu_z

  subroutine s_inv_logdet_lu_z(ndim, zmat, zlogdet)
  !! inverse of a complex matrix, using LU decomposition, also calculate log det of this matrix
     use constants, only : dp, czero, cone, pi

     implicit none

     ! external arguments
     ! dimension of zmat matrix
     integer, intent(in)        :: ndim

     ! object matrix, on entry, it contains the original matrix, on exit,
     ! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)
     complex(dp), intent(out) :: zlogdet

     ! local variables
     ! error flag
     integer     :: ierror
     integer :: i

     ! tmp for imaginary part
     real(dp) :: im_zlogdet

     ! working arrays for lapack subroutines
     integer, allocatable     :: ipiv(:)
     complex(dp), allocatable :: work(:)

     complex(dp), allocatable, dimension(:,:) :: umat, lmat

     integer :: k, pv

     allocate(umat(ndim,ndim))
     allocate(lmat(ndim,ndim))

     ! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     allocate(work(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         write(*,'(a)') 'error in s_inv_logdet_lu_z : can not allocate enough memory'
         stop
     end if

     ! computes the LU factorization of a general m-by-n matrix, need lapack
     ! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         write(*,'(a)') 'error in s_inv_logdet_lu_z : error in lapack subroutine zgetrf'
         stop
     end if

     ! calculate determinant
     zlogdet = czero
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zlogdet = zlogdet + log( +zmat(i,i) )
         else
             zlogdet = zlogdet + log( -zmat(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop
     im_zlogdet = dimag(zlogdet)
     i=nint(im_zlogdet/2.d0/pi)
     im_zlogdet = im_zlogdet - 2.d0*pi*dble(i)
     zlogdet = dcmplx( dble(zlogdet), im_zlogdet )

     !! get inverse of U
     call ztrtri('U','N',ndim,zmat,ndim,ierror)
     umat(:,:) = czero
     do i = 1, ndim
         umat(1:i,i) = zmat(1:i,i)
     end do

     !! get inverse of L
     call ztrtri('L','U',ndim,zmat,ndim,ierror)
     call s_identity_z(ndim,lmat)
     do i = 1, ndim-1
         lmat(i+1:ndim,i) = zmat(i+1:ndim,i)
     end do

     !! L*PT, this part is get from lapack: zgetri
     do k = ndim-1, 1, -1
         pv = ipiv(k)
         if( pv .ne. k ) then
             call zswap( ndim, lmat(1,k), 1, lmat(1,pv), 1)
         end if
     end do

     !! A^-1 = U^-1 * L^-1 * PT
     call zgemm('n','n',ndim,ndim,ndim,cone,umat,ndim,lmat,ndim,czero,zmat,ndim)

     ! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)
     deallocate(lmat)
     deallocate(umat)

     return
  end subroutine s_inv_logdet_lu_z

  subroutine s_adfac_z( rdim, cdim, amat, dvec )
  !! factorize matrix amat to a dense matrix * diagonal matrix
  !! with i-th element of diagonal matrix is norm of i-th column of matrix amat
    use constants, only:dp
    implicit none
    integer, intent(in) :: rdim, cdim
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    real(dp), dimension(cdim), intent(out) :: dvec

    ! local
    integer :: i

    ! external
    real(dp), external :: dznrm2

    do i = 1, cdim
        dvec(i) = dznrm2(rdim,amat(1,i),1)
        amat(:,i) = amat(:,i)/dcmplx(dvec(i),0.d0)
    end do
  end subroutine s_adfac_z

  subroutine s_dafac_z( rdim, cdim, dvec, amat )
  !! factorize matrix amat to a diagonal matrix * dense matrix
  !! with i-th element of diagonal matrix is norm of i-th row of matrix amat
    use constants, only:dp
    implicit none
    integer, intent(in) :: rdim, cdim
    complex(dp), dimension(rdim,cdim), intent(inout) :: amat
    real(dp), dimension(rdim), intent(out) :: dvec

    ! local
    integer :: i

    ! external
    real(dp), external :: dznrm2

    do i = 1, rdim
        dvec(i) = dznrm2(cdim,amat(i,1),rdim)
        amat(i,:) = amat(i,:)/dcmplx(dvec(i),0.d0)
    end do
  end subroutine s_dafac_z
