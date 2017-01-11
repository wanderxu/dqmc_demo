!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_zeros_i
!!!           s_zeros_d
!!!           s_zeros_z
!!!           s_ones_i
!!!           s_ones_d
!!!           s_ones_z
!!!           s_any_i
!!!           s_any_d
!!!           s_any_z
!!!           s_eye_i
!!!           s_eye_d
!!!           s_eye_z
!!!           s_identity_i
!!!           s_identity_d
!!!           s_identity_z
!!!           s_diag_i
!!!           s_diag_d
!!!           s_diag_z

!             s_diag_dz


!!!           s_trace_d
!!!           s_trace_z
!!!           s_det_d
!!!           s_det_z
!!!           s_inv_d
!!!           s_inv_z
!!!           s_eig_dg
!!!           s_eig_zg
!!!           s_eigvals_dg
!!!           s_eigvals_zg
!!!           s_eig_sy
!!!           s_eig_he
!!!           s_eigvals_sy
!!!           s_eigvals_he
!!!           s_solve_dg
!!!           s_solve_zg
!!!           s_solve_sy
!!!           s_solve_he
!!!           s_svd_dg
!!!           s_svd_zg

!             s_compare_max_z
!             s_compare_max_d
!             s_z_x_diag_d
!             s_diag_d_x_z


!!! source  : s_matrix.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/26/2014 by li huang
!!!           11/10/2014 by li huang
!!!           11/28/2014 by yilin wang
!!! purpose : these subroutines are used to encapsulate some important and
!!!           frequently used linear algebra operations.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! 1. build constants (0) matrix
!! -----------------------------
!!
!! subroutine s_zeros_i(...)
!! subroutine s_zeros_d(...)
!! subroutine s_zeros_z(...)
!!
!! 2. build constants (1) matrix
!! -----------------------------
!!
!! subroutine s_ones_i(...)
!! subroutine s_ones_d(...)
!! subroutine s_ones_z(...)
!!
!! 3. build constants (any values) matrix
!! --------------------------------------
!!
!! subroutine s_any_i(...)
!! subroutine s_any_d(...)
!! subroutine s_any_z(...)
!!
!! 4. build diagonal matrix
!! ------------------------
!!
!! subroutine s_eye_i(...)
!! subroutine s_eye_d(...)
!! subroutine s_eye_z(...)
!!
!! 5. build identity matrix
!! ------------------------
!!
!! subroutine s_identity_i(...)
!! subroutine s_identity_d(...)
!! subroutine s_identity_z(...)
!!
!! 6. build diagonal matrix from vector
!! ------------------------------------
!!
!! subroutine s_diag_i(...)
!! subroutine s_diag_d(...)
!! subroutine s_diag_z(...)
!!
!! 7. calculate trace for matrix
!! -----------------------------
!!
!! subroutine s_trace_d(...)
!! subroutine s_trace_z(...)
!!
!! 8. calculate determinant for matrix
!! -----------------------------------
!!
!! subroutine s_det_d(...)
!! subroutine s_det_z(...)
!!
!! 9. calculate matrix inversion
!! -----------------------------
!!
!! subroutine s_inv_d(...)
!! subroutine s_inv_z(...)
!!
!! 10. general eigensystem problem
!! -------------------------------
!!
!! subroutine s_eig_dg(...)
!! subroutine s_eig_zg(...)
!! subroutine s_eigvals_dg(...)
!! subroutine s_eigvals_zg(...)
!!
!! 11. symmetric eigensystem problem
!! ---------------------------------
!!
!! subroutine s_eig_sy(...)
!! subroutine s_eig_he(...)
!! subroutine s_eigvals_sy(...)
!! subroutine s_eigvals_he(...)
!!
!! 12. linear equation solver
!! --------------------------
!!
!! subroutine s_solve_dg(...)
!! subroutine s_solve_zg(...)
!! subroutine s_solve_sy(...)
!! subroutine s_solve_he(...)
!!
!! 13. general singular value decomposition
!! ----------------------------------------
!!
!! subroutine s_svd_dg(...)
!! subroutine s_svd_zg(...)
!!
!! Note: _i means integer version, _d real(dp) version, and _z complex(dp)
!! version. _dg means real(dp) general version, _zg complex(dp) general
!! version, _sy real(dp) symmetric version, _he complex(dp) Hermitian version.
!!
!!

!!========================================================================
!!>>> matrix construction: build zeros/ones/any matrix                 <<<
!!========================================================================

!!>>> s_zeros_i: build an integer matrix with all elements are zero
  subroutine s_zeros_i(n, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input/output matrix
     integer, intent(out) :: A(n,n)

     A = 0

     return
  end subroutine s_zeros_i

!!>>> s_zeros_d: build a real(dp) matrix with all elements are zero
  subroutine s_zeros_d(n, A)
     use constants, only : dp, zero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input/output matrix
     real(dp), intent(out) :: A(n,n)

     A = zero

     return
  end subroutine s_zeros_d

!!>>> s_zeros_z: build a complex(dp) matrix with all elements are zero
  subroutine s_zeros_z(n, A)
     use constants, only : dp, czero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

     A = czero

     return
  end subroutine s_zeros_z

!!>>> s_ones_i: build an integer matrix with all elements are one
  subroutine s_ones_i(n, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input/output matrix
     integer, intent(out) :: A(n,n)

     A = 1

     return
  end subroutine s_ones_i

!!>>> s_ones_d: build a real(dp) matrix with all elements are one
  subroutine s_ones_d(n, A)
     use constants, only : dp, one

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input/output matrix
     real(dp), intent(out) :: A(n,n)

     A = one

     return
  end subroutine s_ones_d

!!>>> s_ones_z: build a complex(dp) matrix with all elements are one
  subroutine s_ones_z(n, A)
     use constants, only : dp, cone

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

     A = cone

     return
  end subroutine s_ones_z

!!>>> s_any_i: build an integer matrix with all elements are given by i
  subroutine s_any_i(n, i, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! value of matrix element
     integer, intent(in)  :: i

! input/output matrix
     integer, intent(out) :: A(n,n)

     A = i

     return
  end subroutine s_any_i

!!>>> s_any_d: build a real(dp) matrix with all elements are given by d
  subroutine s_any_d(n, d, A)
     use constants, only : dp

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! value of matrix element
     real(dp), intent(in)  :: d

! input/output matrix
     real(dp), intent(out) :: A(n,n)

     A = d

     return
  end subroutine s_any_d

!!>>> s_any_z: build a complex(dp) matrix with all elements are given by z
  subroutine s_any_z(n, z, A)
     use constants, only : dp

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! value of matrix element
     complex(dp), intent(in)  :: z

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

     A = z

     return
  end subroutine s_any_z

!!========================================================================
!!>>> matrix construction: build diagonal matrix                       <<<
!!========================================================================

!!>>> s_eye_i: build integer matrix with ones on the diagonal and zeros elsewhere.
  subroutine s_eye_i(n, k, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! index of the diagonal: 0 refers to the main diagonal, a positive value
! refers to an upper diagonal, and a negative value to a lower diagonal.
     integer, intent(in)  :: k

! input/output matrix
     integer, intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = 0
     do i=1,n
         if ( i - k < 1 .or. i - k > n ) CYCLE
         A(i,i-k) = 1
     enddo ! over i={1,n} loop

     return
  end subroutine s_eye_i

!!>>> s_eye_d: build real(dp) matrix with ones on the diagonal and zeros elsewhere.
  subroutine s_eye_d(n, k, A)
     use constants, only : dp, zero, one

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! index of the diagonal: 0 refers to the main diagonal, a positive value
! refers to an upper diagonal, and a negative value to a lower diagonal.
     integer, intent(in)   :: k

! input/output matrix
     real(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = zero
     do i=1,n
         if ( i - k < 1 .or. i - k > n ) CYCLE
         A(i,i-k) = one
     enddo ! over i={1,n} loop

     return
  end subroutine s_eye_d

!!>>> s_eye_z: build complex(dp) matrix with ones on the diagonal and zeros elsewhere.
  subroutine s_eye_z(n, k, A)
     use constants, only : dp, czero, cone

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! index of the diagonal: 0 refers to the main diagonal, a positive value
! refers to an upper diagonal, and a negative value to a lower diagonal.
     integer, intent(in)      :: k

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = czero
     do i=1,n
         if ( i - k < 1 .or. i - k > n ) CYCLE
         A(i,i-k) = cone
     enddo ! over i={1,n} loop

     return
  end subroutine s_eye_z

!!>>> s_identity_i: build integer identity matrix
  subroutine s_identity_i(n, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input/output matrix
     integer, intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = 0
     do i=1,n
         A(i,i) = 1
     enddo ! over i={1,n} loop

     return
  end subroutine s_identity_i

!!>>> s_identity_d: build real(dp) identity matrix
  subroutine s_identity_d(n, A)
     use constants, only : dp, zero, one

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input/output matrix
     real(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = zero
     do i=1,n
         A(i,i) = one
     enddo ! over i={1,n} loop

     return
  end subroutine s_identity_d

!!>>> s_identity_z: build complex(dp) identity matrix
  subroutine s_identity_z(n, A)
     use constants, only : dp, czero, cone

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = czero
     do i=1,n
         A(i,i) = cone
     enddo ! over i={1,n} loop

     return
  end subroutine s_identity_z

!!>>> s_diag_i: build integer diagonal matrix from a vector
  subroutine s_diag_i(n, v, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input integer vector
     integer, intent(in)  :: v(n)

! output integer diagonal matrix
     integer, intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = 0
     do i=1,n
         A(i,i) = v(i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_diag_i

!!>>> s_diag_d: build real(dp) diagonal matrix from a vector
  subroutine s_diag_d(n, v, A)
     use constants, only : dp, zero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input real(dp) vector
     real(dp), intent(in)  :: v(n)

! output real(dp) diagonal matrix
     real(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = zero
     do i=1,n
         A(i,i) = v(i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_diag_d

!!>>> s_diag_z: build complex(dp) diagonal matrix from a vector
  subroutine s_diag_z(n, v, A)
     use constants, only : dp, czero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input complex(dp) vector
     complex(dp), intent(in)  :: v(n)

! output complex(dp) diagonal matrix
     complex(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = czero
     do i=1,n
         A(i,i) = v(i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_diag_z

!!>>> s_diag_dz: build complex(dp) diagonal matrix from a vector
  subroutine s_diag_dz(n, v, A)
     use constants, only : dp, czero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input complex(dp) vector
     real(dp), intent(in)  :: v(n)

! output complex(dp) diagonal matrix
     complex(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = czero
     do i=1,n
         A(i,i) = dcmplx( v(i), 0.d0 )
     enddo ! over i={1,n} loop

     return
  end subroutine s_diag_dz

!!========================================================================
!!>>> matrix query: return matrix's trace or determinant               <<<
!!========================================================================

!!>>> s_trace_d: return trace for a real(dp) array
  subroutine s_trace_d(n, A, tr)
     use constants, only : dp, zero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! output matrix's trace
     real(dp), intent(out) :: tr

! input real(dp) matrix
     real(dp), intent(in)  :: A(n,n)

! local variables
! loop index
     integer :: i

     tr = zero
     do i=1,n
         tr = tr + A(i,i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_trace_d

!!>>> s_trace_z: return trace for a complex(dp) array
  subroutine s_trace_z(n, A, tr)
     use constants, only : dp, czero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! output matrix's trace
     complex(dp), intent(out) :: tr

! input complex(dp) matrix
     complex(dp), intent(in)  :: A(n,n)

! local variables
! loop index
     integer :: i

     tr = czero
     do i=1,n
         tr = tr + A(i,i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_trace_z

!!>>> s_det_d: calculate the determinant of a real(dp) matrix
  subroutine s_det_d(ndim, dmat, ddet)
     use constants, only : dp, one, cone

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in)     :: ndim

! determinant of dmat matrix
     real(dp), intent(out)   :: ddet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! loop index
     integer  :: i

! error flag
     integer  :: ierror

! size of working array work
     integer  :: lwork

! used to calculate determinant
     complex(dp) :: cres

! working arrays for lapack subroutines: dgetrf
     integer, allocatable  :: ipiv(:)

! working arrays for lapack subroutines: dgeev
     real(dp), allocatable :: work(:)

! real and imaginary parts of the computed eigenvalues
     real(dp), allocatable :: wi(:)
     real(dp), allocatable :: wr(:)

! left and right eigenvectors
     real(dp), allocatable :: vl(:,:)
     real(dp), allocatable :: vr(:,:)

! dummy arrays, used to save dmat
     real(dp), allocatable :: amat(:,:)

! setup lwork
     lwork = 4 * ndim

! allocate memory
     allocate(ipiv(ndim),      stat=ierror)
     allocate(work(lwork),     stat=ierror)
     allocate(wi(ndim),        stat=ierror)
     allocate(wr(ndim),        stat=ierror)
     allocate(vl(ndim,ndim),   stat=ierror)
     allocate(vr(ndim,ndim),   stat=ierror)
     allocate(amat(ndim,ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_d','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! copy dmat to amat at first
     amat = dmat

!-------------------------------------------------------------------------
! method A: preferred method
!-------------------------------------------------------------------------
! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call DGETRF(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_exception('s_det_d','error in lapack subroutine dgetrf')
     endif ! back if ( ierror /= 0 ) block

! calculate determinant
     ddet = one
     do i=1,ndim
         if ( ipiv(i) == i ) then
             ddet = ddet * ( +dmat(i,i) )
         else
             ddet = ddet * ( -dmat(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop

! everything is ok!
     if ( ierror == 0 ) RETURN

!-------------------------------------------------------------------------
! method B: as a backup
!-------------------------------------------------------------------------
! diagonalize amat to obtain its eigenvalues: wr and wi
     call DGEEV('N', 'N', ndim, amat, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_d','error in lapack subroutine dgeev')
     endif ! back if ( ierror /= 0 ) block

! evaluate the final determinant
     cres = cone
     do i=1,ndim
         cres = cres * dcmplx( wr(i), wi(i) )
     enddo ! over i={1,ndim} loop
     ddet = cres

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)
     if ( allocated(wi  ) ) deallocate(wi  )
     if ( allocated(wr  ) ) deallocate(wr  )
     if ( allocated(vl  ) ) deallocate(vl  )
     if ( allocated(vr  ) ) deallocate(vr  )
     if ( allocated(amat) ) deallocate(amat)

     return
  end subroutine s_det_d

!!>>> s_det_z: calculate the determinant of a complex(dp) matrix
  subroutine s_det_z(ndim, zmat, zdet)
     use constants, only : dp, cone

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in)        :: ndim

! determinant of zmat matrix
     complex(dp), intent(out)   :: zdet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

! error flag
     integer :: ierror

! working arrays for lapack subroutines
     integer, allocatable :: ipiv(:)

! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_z','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_z','error in lapack subroutine zgetrf')
     endif ! back if ( ierror /= 0 ) block

! calculate determinant
     zdet = cone
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zdet = zdet * ( +zmat(i,i) )
         else
             zdet = zdet * ( -zmat(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_det_z

!!>>> s_logdet_z: calculate the log of determinant of a complex(dp) matrix
  subroutine s_logdet_z(ndim, zmat, zlogdet)
     use constants, only : dp, czero, pi

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in)        :: ndim

! determinant of zmat matrix
     complex(dp), intent(out)   :: zlogdet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

     ! tmp for imaginary part
     real(dp) :: im_zlogdet

! error flag
     integer :: ierror

! working arrays for lapack subroutines
     integer, allocatable :: ipiv(:)

! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_logdet_z','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_logdet_z','error in lapack subroutine zgetrf')
     endif ! back if ( ierror /= 0 ) block

! calculate determinant
     zlogdet = czero
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zlogdet = zlogdet + log( +zmat(i,i) )
         else
             zlogdet = zlogdet + log( -zmat(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop
     im_zlogdet = aimag(zlogdet)
     i=nint(im_zlogdet/2.d0/pi)
     im_zlogdet = im_zlogdet - 2.d0*pi*dble(i)
     zlogdet = dcmplx( real(zlogdet), im_zlogdet )

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_logdet_z

!!>>> s_logdet_z: calculate the log of determinant of a complex(dp) matrix
  subroutine s_logdet_zc2(ndim, zmat, zlogdet)
     use constants, only : dp, czero, pi

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in)        :: ndim

! determinant of zmat matrix
     complex(dp), intent(out)   :: zlogdet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

     ! tmp for imaginary part
     real(dp) :: im_zlogdet

! error flag
     integer :: ierror

! working arrays for lapack subroutines
     integer, allocatable :: ipiv(:), jpiv(:)

! allocate memory
     allocate(ipiv(ndim), jpiv(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_logdet_zc2','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetc2 subroutine
     call ZGETC2(ndim, zmat, ndim, ipiv, jpiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_logdet_zc2','error in lapack subroutine zgetrf')
     endif ! back if ( ierror /= 0 ) block

! calculate determinant
     zlogdet = czero
     do i=1,ndim
         if ( ( (ipiv(i) .eq. i) .and. (jpiv(i) .eq. i) ) .and. ( (ipiv(i) .ne. i) .and. (jpiv(i) .ne. i) ) ) then
             zlogdet = zlogdet + log( +zmat(i,i) )
         else
             zlogdet = zlogdet + log( -zmat(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop
     im_zlogdet = aimag(zlogdet)
     i=nint(im_zlogdet/2.d0/pi)
     im_zlogdet = im_zlogdet - 2.d0*pi*dble(i)
     zlogdet = dcmplx( real(zlogdet), im_zlogdet )

! deallocate memory
     if ( allocated(jpiv) ) deallocate(jpiv)
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_logdet_zc2

!!========================================================================
!!>>> matrix manipulation: calculate matrix's inversion                <<<
!!========================================================================

!!>>> s_inv_d: invert real(dp) matrix using lapack subroutines
  subroutine s_inv_d(ndim, dmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in)     :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! error flag
     integer  :: ierror

! working arrays for lapack subroutines
     integer, allocatable  :: ipiv(:)
     real(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     allocate(work(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_d','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call DGETRF(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_d','error in lapack subroutine dgetrf')
     endif ! back if ( ierror /= 0 ) block

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call DGETRI(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_d','error in lapack subroutine dgetri')
     endif ! back if ( ierror /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_inv_d

!!>>> s_inv_z: invert complex(dp) matrix using lapack subroutines
  subroutine s_inv_z(ndim, zmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in)        :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! error flag
     integer     :: ierror

! working arrays for lapack subroutines
     integer, allocatable     :: ipiv(:)
     complex(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     allocate(work(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_z','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_z','error in lapack subroutine zgetrf')
     endif ! back if ( ierror /= 0 ) block

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call ZGETRI(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_z','error in lapack subroutine zgetri')
     endif ! back if ( ierror /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_inv_z

!!========================================================================
!!>>> matrix manipulation: solve eigenvalues and eigenvectors problem  <<<
!!========================================================================

!!>>> s_eig_dg: diagonalize a general real(dp) matrix and return eigenvalues
!!>>> and eigenvectors
  subroutine s_eig_dg(ldim, ndim, amat, eval, evec)
     use constants, only : dp, zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original general real(dp) matrix to compute eigenvals and eigenvectors
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     real(dp), intent(out) :: evec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dgeev
     integer :: info

! the length of the array work, lwork >= max(1,4*ndim)
     integer :: lwork

! workspace array
     real(dp), allocatable :: work(:)

! auxiliary real(dp) matrix: real and imaginary parts of eigenvalues
     real(dp), allocatable :: wr(:)
     real(dp), allocatable :: wi(:)

! auxiliary real(dp) matrix: left and right eigenvectors
     real(dp), allocatable :: vr(:,:)
     real(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 4 * ndim

! allocate memory
     allocate(work(lwork),   stat=istat)
     allocate(wr(ndim),      stat=istat)
     allocate(wi(ndim),      stat=istat)
     allocate(vr(ndim,ndim), stat=istat)
     allocate(vl(ndim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eig_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dgeev
     call DGEEV('N', 'V', ndim, evec, ldim, wr, wi, vl, ndim, vr, ndim, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_dg','error in lapack subroutine dgeev')
     endif ! back if ( info /= 0 ) block

! copy eigenvalues and eigenvectors
     eval(1:ndim) = wr(1:ndim)
     evec(1:ndim,1:ndim) = vr(1:ndim,1:ndim)

! dealloate memory for workspace array
     if ( allocated(work) ) deallocate(work)
     if ( allocated(wr  ) ) deallocate(wr  )
     if ( allocated(wi  ) ) deallocate(wi  )
     if ( allocated(vr  ) ) deallocate(vr  )
     if ( allocated(vl  ) ) deallocate(vl  )

     return
  end subroutine s_eig_dg

!!>>> s_eig_zg: diagonalize a general complex(dp) matrix and return eigenvalues
!!>>> and eigenvectors
  subroutine s_eig_zg(ldim, ndim, zmat, zeig, zvec)
     use constants, only : dp, czero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)      :: ldim

! the order of the matrix amat
     integer, intent(in)      :: ndim

! original general complex(dp) matrix to compute eigenvals and eigenvectors
     complex(dp), intent(in)  :: zmat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     complex(dp), intent(out) :: zeig(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     complex(dp), intent(out) :: zvec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zgeev
     integer :: info

! the length of the array work, lwork >= max(1,2*ndim)
     integer :: lwork

! workspace array
     complex(dp), allocatable :: work(:)

! auxiliary real(dp) matrix
     complex(dp), allocatable :: rwork(:)

! auxiliary complex(dp) matrix: left and right eigenvectors
     complex(dp), allocatable :: vr(:,:)
     complex(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 2 * ndim

! allocate memory
     allocate(work(lwork),   stat=istat)
     allocate(rwork(lwork),  stat=istat)
     allocate(vr(ndim,ndim), stat=istat)
     allocate(vl(ndim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eig_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     zeig = czero
     zvec = zmat

! call the computational subroutine: zgeev
     call ZGEEV('N', 'V', ndim, zvec, ldim, zeig, vl, ndim, vr, ndim, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_zg','error in lapack subroutine zgeev')
     endif ! back if ( info /= 0 ) block

! copy eigenvectors
     zvec = vr

! dealloate memory for workspace array
     if ( allocated(work ) )  deallocate(work )
     if ( allocated(rwork) )  deallocate(rwork)
     if ( allocated(vr   ) )  deallocate(vr   )
     if ( allocated(vl   ) )  deallocate(vl   )

     return
  end subroutine s_eig_zg

!!>>> s_eigvals_dg: diagonalize a general real(dp) matrix and return eigenvalues only
  subroutine s_eigvals_dg(ldim, ndim, amat, eval)
     use constants, only : dp, zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original general real(dp) matrix to compute eigenvals
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dgeev
     integer :: info

! the length of the array work, lwork >= max(1,4*ndim)
     integer :: lwork

! workspace array, used to store amat
     real(dp), allocatable :: evec(:,:)

! workspace array
     real(dp), allocatable :: work(:)

! auxiliary real(dp) matrix: real and imaginary parts of eigenvalues
     real(dp), allocatable :: wr(:)
     real(dp), allocatable :: wi(:)

! auxiliary real(dp) matrix: left and right eigenvectors
     real(dp), allocatable :: vr(:,:)
     real(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 4 * ndim

! allocate memory
     allocate(evec(ldim,ndim), stat=istat)
     allocate(work(lwork),     stat=istat)
     allocate(wr(ndim),        stat=istat)
     allocate(wi(ndim),        stat=istat)
     allocate(vr(ndim,ndim),   stat=istat)
     allocate(vl(ndim,ndim),   stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dgeev
     call DGEEV('N', 'N', ndim, evec, ldim, wr, wi, vl, ndim, vr, ndim, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_dg','error in lapack subroutine dgeev')
     endif ! back if ( info /= 0 ) block

! copy eigenvalues
     eval(1:ndim) = wr(1:ndim)

! dealloate memory for workspace array
     if ( allocated(evec) ) deallocate(evec)
     if ( allocated(work) ) deallocate(work)
     if ( allocated(wr  ) ) deallocate(wr  )
     if ( allocated(wi  ) ) deallocate(wi  )
     if ( allocated(vr  ) ) deallocate(vr  )
     if ( allocated(vl  ) ) deallocate(vl  )

     return
  end subroutine s_eigvals_dg

!!>>> s_eigvals_zg: diagonalize a general complex(dp) matrix and return eigenvalues only
  subroutine s_eigvals_zg(ldim, ndim, zmat, zeig)
     use constants, only : dp, czero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)      :: ldim

! the order of the matrix amat
     integer, intent(in)      :: ndim

! original general complex(dp) matrix to compute eigenvals
     complex(dp), intent(in)  :: zmat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     complex(dp), intent(out) :: zeig(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zgeev
     integer :: info

! the length of the array work, lwork >= max(1,2*ndim)
     integer :: lwork

! workspace array, used to store amat
     complex(dp), allocatable :: zvec(:,:)

! workspace array
     complex(dp), allocatable :: work(:)

! auxiliary real(dp) matrix
     complex(dp), allocatable :: rwork(:)

! auxiliary complex(dp) matrix: left and right eigenvectors
     complex(dp), allocatable :: vr(:,:)
     complex(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 2 * ndim

! allocate memory
     allocate(zvec(ldim,ndim), stat=istat)
     allocate(work(lwork),     stat=istat)
     allocate(rwork(lwork),    stat=istat)
     allocate(vr(ndim,ndim),   stat=istat)
     allocate(vl(ndim,ndim),   stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     zeig = czero
     zvec = zmat

! call the computational subroutine: zgeev
     call ZGEEV('N', 'N', ndim, zvec, ldim, zeig, vl, ndim, vr, ndim, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_zg','error in lapack subroutine zgeev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(zvec ) )  deallocate(zvec )
     if ( allocated(work ) )  deallocate(work )
     if ( allocated(rwork) )  deallocate(rwork)
     if ( allocated(vr   ) )  deallocate(vr   )
     if ( allocated(vl   ) )  deallocate(vl   )

     return
  end subroutine s_eigvals_zg

!!>>> s_eig_sy: computes all eigenvalues and eigenvectors of real symmetric matrix
  subroutine s_eig_sy(ldim, ndim, amat, eval, evec)
     use constants, only : dp, zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original real symmetric matrix to compute eigenvals and eigenvectors
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     real(dp), intent(out) :: evec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,3*ndim-1)
     integer :: lwork

! workspace array
     real(dp), allocatable :: work(:)

! initialize lwork
     lwork = 3 * ndim - 1

! allocate memory
     allocate(work(lwork), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eig_sy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dsyev
     call DSYEV('V', 'U', ndim, evec, ldim, eval, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_sy','error in lapack subroutine dsyev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_eig_sy

!!>>> s_eig_he: computes all eigenvalues and eigenvectors of complex Hermitian matrix
  subroutine s_eig_he(ldim, ndim, amat, eval, evec)
     use constants, only : dp, zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)      :: ldim

! the order of the matrix amat
     integer, intent(in)      :: ndim

! original complex Hermitian matrix to compute eigenvals and eigenvectors
     complex(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out)    :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     complex(dp), intent(out) :: evec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zheev
     integer :: info

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
         call s_print_error('s_eig_he','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: zheev
     call ZHEEV('V', 'U', ndim, evec, ldim, eval, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_he','error in lapack subroutine zheev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)

     return
  end subroutine s_eig_he

!!>>> s_eigvals_sy: computes all eigenvalues of real symmetric matrix
  subroutine s_eigvals_sy(ldim, ndim, amat, eval)
     use constants, only : dp, zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original real symmetric matrix to compute eigenvals
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,3*ndim-1)
     integer :: lwork

! workspace array
     real(dp), allocatable :: work(:)

! workspace array, used to store amat
     real(dp), allocatable :: evec(:,:)

! initialize lwork
     lwork = 3 * ndim - 1

! allocate memory
     allocate(work(lwork),     stat=istat)
     allocate(evec(ldim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_sy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dsyev
     call DSYEV('N', 'U', ndim, evec, ldim, eval, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_sy','error in lapack subroutine dsyev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work) ) deallocate(work)
     if ( allocated(evec) ) deallocate(evec)

     return
  end subroutine s_eigvals_sy

!!>>> s_eigvals_he: computes all eigenvalues of complex Hermitian matrix
  subroutine s_eigvals_he(ldim, ndim, amat, eval)
     use constants, only : dp, zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)     :: ldim

! the order of the matrix amat
     integer, intent(in)     :: ndim

! original complex Hermitian matrix to compute eigenvals and eigenvectors
     complex(dp), intent(in) :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out)   :: eval(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zheev
     integer :: info

! the length of the array work and rwork
! lwork >= max(1,2*ndim-1), lrwork >= max(1,3*ndim-2)
     integer :: lwork
     integer :: lrwork

! workspace array
     real(dp), allocatable    :: rwork(:)
     complex(dp), allocatable :: work(:)

! workspace array, used to store amat
     complex(dp), allocatable :: evec(:,:)

! initialize lwork (lrwork)
     lwork = 2 * ndim - 1
     lrwork = 3 * ndim - 2

! allocate memory
     allocate(work(lwork),     stat=istat)
     allocate(rwork(lrwork),   stat=istat)
     allocate(evec(ldim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_he','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: zheev
     call ZHEEV('N', 'U', ndim, evec, ldim, eval, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_he','error in lapack subroutine zheev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)
     if ( allocated(evec ) ) deallocate(evec )

     return
  end subroutine s_eigvals_he

!!========================================================================
!!>>> matrix manipulation: solve linear equations                      <<<
!!========================================================================

!!>>> s_solve_dg: solve linear system AX = B, real(dp) general version
  subroutine s_solve_dg(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)     :: n

! the number of right-hand sides
     integer, intent(in)     :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     real(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     real(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine dgesv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable :: ipiv(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: dgesv
     call DGESV(n, nrhs, A, n, ipiv, B, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_dg','error in lapack subroutine dgesv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_solve_dg

!!>>> s_solve_zg: solve linear system AX = B, complex(dp) general version
  subroutine s_solve_zg(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)        :: n

! the number of right-hand sides
     integer, intent(in)        :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     complex(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     complex(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine zgesv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable :: ipiv(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: zgesv
     call ZGESV(n, nrhs, A, n, ipiv, B, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_zg','error in lapack subroutine zgesv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_solve_zg

!!>>> s_solve_sy: solve linear system AX = B, real(dp) symmetric version
  subroutine s_solve_sy(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)     :: n

! the number of right-hand sides
     integer, intent(in)     :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     real(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     real(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine dsysv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable  :: ipiv(:)

! workspace array, its dimension is at least max(1, lwork) and lwork >= 1
     real(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     allocate(work(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_sy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: dsysv
     call DSYSV('U', n, nrhs, A, n, ipiv, B, n, work, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_sy','error in lapack subroutine dsysv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_solve_sy

!!>>> s_solve_he: solve linear system AX = B, complex(dp) Hermitian version
  subroutine s_solve_he(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)        :: n

! the number of right-hand sides
     integer, intent(in)        :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     complex(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     complex(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine zhesv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable     :: ipiv(:)

! workspace array, its dimension is at least max(1, lwork) and lwork >= 1
     complex(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     allocate(work(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_he','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: zhesv
     call ZHESV('U', n, nrhs, A, n, ipiv, B, n, work, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_he','error in lapack subroutine zhesv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_solve_he

!!========================================================================
!!>>> matrix manipulation: singular values decomposition               <<<
!!========================================================================

!!>>> s_svd_dg: perform the singular values decomposition for a general
!!>>> real(dp) m-by-n matrix A, where A = U * SIGMA * transpose(V), return
!!>>> its left vectors, right vectors, and singular values
  subroutine s_svd_dg(m, n, min_mn, amat, umat, svec, vmat)
     use constants, only : dp
 
     implicit none

! external arguments
! number of rows of A matrix
     integer, intent(in)     :: m

! number of columns of A matrix
     integer, intent(in)     :: n

! minimal value of m and n
     integer, intent(in)     :: min_mn

! A matrix
     real(dp), intent(inout) :: amat(m,n)

! left vectors of svd, U
     real(dp), intent(out)   :: umat(m,min_mn)

! singular values of svd, SIGMA
     real(dp), intent(out)   :: svec(min_mn)

! right vectors of svd, transpose(V)
     real(dp), intent(out)   :: vmat(min_mn,n)

! local variables
! status flag
     integer :: istat

! return information from dgesvd
     integer :: info

! length of work array, lwork >= max(1, 3 * min_mn + max(m,n), 5 * min_mn)
     integer :: lwork 

! workspace array
     real(dp), allocatable :: work(:)

! initialize lwrok
     lwork = max(1, 3 * min_mn + max(m,n), 5 * min_mn)

! allocate memory
     allocate(work(lwork), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_svd_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: dgesvd
     call DGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_svd_dg','error in lapack subroutine dgesvd')
     endif ! back if ( info /= 0 ) block

! deallocate the memory for workspace array
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_svd_dg

!!>>> s_svd_zg: perform the singular values decomposition for a general
!!>>> complex(dp) m-by-n matrix A, where A = U * SIGMA * conjugate-transpose(V),
!!>>> return its left vectors, right vectors, and singular values
  subroutine s_svd_zg(m, n, min_mn, amat, umat, svec, vmat)
     use constants, only : dp
 
     implicit none

! external arguments
! number of rows of A matrix
     integer, intent(in)        :: m

! number of columns of A matrix
     integer, intent(in)        :: n

! minimal value of m and n
     integer, intent(in)        :: min_mn

! A matrix
     complex(dp), intent(inout) :: amat(m,n)

! left vectors of svd, U
     complex(dp), intent(out)   :: umat(m,min_mn)

! singular values of svd, SIGMA
     real(dp), intent(out)      :: svec(min_mn)

! right vectors of svd, conjugate-transpose(V)
     complex(dp), intent(out)   :: vmat(min_mn,n)

! local variables
! status flag
     integer :: istat

! return information from zgesvd
     integer :: info

! length of work array, lwork >= max(1, 2 * min_mn + max(m,n))
     integer :: lwork 

! workspace arrays
     complex(dp), allocatable :: work(:)
     real(dp), allocatable :: rwork(:)

!!!!! initialize lwrok
!!!!     lwork = max(1, 2 * min_mn + max(m,n))

! allocate memory
     allocate(rwork(5*min_mn), stat=istat)

     ! lwork query
     lwork = -1
     allocate(work(1),  stat=istat)
     call ZGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, rwork, info)
     lwork = nint(dble(work(1)))
     deallocate(work)

     ! with optimal lwork, do SVD
     allocate(work(lwork),  stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_svd_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: zgesvd
     call ZGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_svd_zg','error in lapack subroutine zgesvd')
     endif ! back if ( info /= 0 ) block

! deallocate the memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)

     return
  end subroutine s_svd_zg

  ! following is added by Xiao Yan Xu (wanderxu@gmail.com)
  subroutine s_compare_max_z( ndim, Amat, Bmat, max_diff )
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
        call s_print_error('s_zgeQRPT, above zgeqp3','can not allocate enough memory')
    endif
    jpvt(:) = 0
    call zgeqp3(rdim, cdim, amat, rdim, jpvt, tau, work, lwork, rwork, ierror)
    if ( ierror /= 0 ) then
        call s_print_error('s_zgeQRPT, ','error in lapack subroutine zgeqp3')
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
        call s_print_error('s_zgeQRPT, above zungqr','can not allocate enough memory')
    endif
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        call s_print_error('s_zgeQRPT, ','error in lapack subroutine zungqr')
    endif
    deallocate(work)
    deallocate(tau)

  end subroutine s_zgeQRPT

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
        call s_print_error('s_zgeQR, above zgeqrf','can not allocate enough memory')
    endif
    call zgeqrf(rdim, cdim, amat, rdim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        call s_print_error('s_zgeQR, ','error in lapack subroutine zgeqrf')
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
        call s_print_error('s_zgeQR, above zungqr','can not allocate enough memory')
    endif
    call zungqr(mindim, mindim, mindim, qmat, mindim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        call s_print_error('s_zgeQR, ','error in lapack subroutine zungqr')
    endif
    deallocate(work)
    deallocate(tau)

  end subroutine s_zgeQR

  subroutine s_invqr_z(ndim, amat )
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
        call s_print_error('s_zgeQR, above zgeqrf','can not allocate enough memory')
    endif
    call zgeqrf(ndim, ndim, amat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        call s_print_error('s_zgeQR, ','error in lapack subroutine zgeqrf')
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
        call s_print_error('s_zgeQR, above zungqr','can not allocate enough memory')
    endif
    call zungqr(ndim, ndim, ndim, qmat, ndim, tau, work, lwork, ierror)
    if ( ierror /= 0 ) then
        call s_print_error('s_zgeQR, ','error in lapack subroutine zungqr')
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

  !!>>> s_invlu_z: invert complex(dp) matrix using lapack subroutines
  subroutine s_invlu_z(ndim, zmat)
     use constants, only : dp, czero, cone

     implicit none

     ! external arguments
     ! dimension of zmat matrix
     integer, intent(in)        :: ndim

     ! object matrix, on entry, it contains the original matrix, on exit,
     ! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

     ! local variables
     ! error flag
     integer     :: ierror
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
         call s_print_error('s_inv_z','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

     ! computes the LU factorization of a general m-by-n matrix, need lapack
     ! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_z','error in lapack subroutine zgetrf')
     endif ! back if ( ierror /= 0 ) block

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

  subroutine s_adfac_z( rdim, cdim, amat, dvec )
  !! factorize matrix amat to a dense matrix * diagonal matrix
  !! with i-th element of diagnoal matrix is norm of i-th column of matrix amat
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
  !! with i-th element of diagnoal matrix is norm of i-th row of matrix amat
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
