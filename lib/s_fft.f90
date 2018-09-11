  subroutine discrete_fft2d( lx, ly, ndim, X )
  ! 2D complex to complex discrete fourier transform
  Use MKL_DFTI
  Use constants, only : dp
  implicit none
  integer, intent(in) :: lx, ly, ndim
  Complex(dp), intent(inout) ::  X(ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
  Integer :: Status, LN(2)
  ! if the data is in a 2d matrix X_2D, use equivalence(X_2D,X) before calling this subroutine
  
  LN(1) = lx; LN(2) = ly
   
  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
            DFTI_COMPLEX, 2, LN)
  Status = DftiCommitDescriptor( My_Desc1_Handle)
  Status = DftiComputeForward( My_Desc1_Handle, X)
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by X
  
  end subroutine discrete_fft2d

  subroutine discrete_fft3d( lx, ly, lz, ndim, X )
  ! 3D complex to complex discrete fourier transform
  Use MKL_DFTI
  Use constants, only : dp
  implicit none
  integer, intent(in) :: lx, ly, lz, ndim
  Complex(dp), intent(inout) ::  X(ndim)
  type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
  Integer :: Status, LN(3)
  ! if the data is in a 3d matrix X_3D, use equivalence(X_3D,X) before calling this subroutine
  
  LN(1) = lx; LN(2) = ly; LN(3)=lz
   
  ! Perform a complex to complex transform
  Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,&
            DFTI_COMPLEX, 3, LN)
  Status = DftiCommitDescriptor( My_Desc1_Handle)
  Status = DftiComputeForward( My_Desc1_Handle, X)
  Status = DftiFreeDescriptor(My_Desc1_Handle)
  ! result is given by X
  
  end subroutine discrete_fft3d
