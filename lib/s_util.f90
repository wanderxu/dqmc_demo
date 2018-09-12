  subroutine s_heapsort(n,width,brxyz)
     use constants, only : dp
! modified from BSTATE code by Xiaoyan Xu(wanderxu@gmail.com) 2013.9
!########################################################################
! Former discription
!C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
!C-----CALL HPSORT(KG,IG1,IG2,IG3,GX,GY,GZ,GR)
!      subroutine hpsort(kn,n,bxyz,br)
!c                           @(#)hpsort.f 9.1 97/05/08 14:48:33 
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!C     THIS SUBROUTINE SORTS THE ELEMENTS OF ARRAY BR IN ASCENDING
!C     ORDER. THE CODE IS BASED UPON HEAPSORT ALGORITHM WHICH HAS
!C     N*LOG(N) SCALING BEHAVIOUR, PARTICULARLY FAVOURABLE FOR LARGE
!C     VECTORS BR
!C                              STEFAN BL"UGEL, ISSP, NOV 1989
!C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!########################################################################


! description
! This subroutine sorts the elements of array brxyz(width,n) in ascending order
! based on brxyz(1,:).

     IMPLICIT LOGICAL(A-Z)
     integer n, width
     real(dp), dimension(width, n) :: brxyz
     integer :: l, i, ii, iheap
     real(dp) :: brr
     real(dp), dimension(width-1) :: bxyz
!     integer :: kn, n,l,i,ii,iheap
!     real*8    bxyz(kn,3),br(*),brr,bxx,byy,bzz
!C
!C
!C=====> BUILD-UP OF THE HEAP
!C
!C-----> LOOP OVER ALL HIERACHICAL LEVELS OF THE HEAP TREE
!C
      DO 10 L = N/2 , 1 , -1
!         brr = br(l)
         brr = brxyz(1,l)
         bxyz = brxyz(2:width,l)
!         bxx = bxyz(l,1)
!         byy = bxyz(l,2)
!         bzz = bxyz(l,3)
         I   = L
         II  = L + L
!C
!C-----> GO DOWN ALL THE HIERACHICAL LEVELS OF THE HEAP TREE
!C
  20    IF ( II .LE. N ) THEN
!C
!C-----> COMPARE NEIGHBOURING ELEMENTS
           IF ( II .LT. N ) THEN
              IF ( BRXYZ(1,II) .LT. BRXYZ(1,II+1) ) II = II + 1
           ENDIF
!C
!C-----> COMPARE THE ELEMENTS OF TWO HIRACHICAL LEVELS
!C       PROMOTE LARGER ONE, DEMOTE SMALLER ONE
           IF ( BRR .LT. BRXYZ(1,II) ) THEN
              BRXYZ(1,I) = BRXYZ(1,II)
              brxyz(2:width,i) = brxyz(2:width,ii)
!              bxyz(i,1) = bxyz(ii,1)
!              bxyz(i,2) = bxyz(ii,2)
!              bxyz(i,3) = bxyz(ii,3)
              I     = II
              II    = II + II
           ELSE
!C
!C-----> THIS PART OF THE TREE IS ORDERED , STOP BY PUTTING II=N+1
              II    = N + 1
           END IF
                                                GO TO 20
        END IF
!C-----> PUT ELEMENTS IN THE PROPER SLOT
        brxyz(1,i) = brr
        brxyz(2:width,i) = bxyz(:)
!        bxyz(i,1) = bxx
!        bxyz(i,2) = byy
!        bxyz(i,3) = bzz
  10  continue
!C
!C=====> NOW COLLECT ALL ELEMENTS FROM THE HEAP
!C
      DO 30 IHEAP = N , 1 , -1
!C
         brr = brxyz(1,iheap)
         bxyz = brxyz(2:width,iheap)
!         bxx = bxyz(iheap,1)
!         byy = bxyz(iheap,2)
!         bzz = bxyz(iheap,3)
!C
!C-----> THE FIRST ELEMENT IS ALWAYS THE LARGEST
              brxyz(1,iheap) = brxyz(1,1)
              brxyz(2:width,iheap) = brxyz(2:width,1)
!              bxyz(iheap,1) = bxyz(1,1)
!              bxyz(iheap,2) = bxyz(1,2)
!              bxyz(iheap,3) = bxyz(1,3)
!C-----> NOW LARGEST ELEMENT OF ALL BR(I) WITH 1<=I<=IHEAP IS STORED
!C
         I   = 1
         II  = 2
!C
!C-----> NOW GENERATE LARGEST ELEMENT OF BR(I) WITH 1<=I<=IHEAP-1
!C
  40    IF ( II .LE. IHEAP - 1 ) THEN
!C
!C-----> COMPARE NEIGHBOURING ELEMENTS
           IF ( II .LT. IHEAP - 1 ) THEN
              IF ( BRXYZ(1,II) .LT. BRXYZ(1,II+1) ) II = II + 1
           ENDIF
!C
!C-----> PROMOTE EACH ELEMENT OF THE TWIG OF BR UNTIL BRR > BR(I)
           IF ( BRR .LT. BRXYZ(1,II) ) THEN
              brxyz(1,i) = brxyz(1,ii)
              brxyz(2:width,i) = brxyz(2:width,ii)
!              bxyz(i,1) = bxyz(ii,1)
!              bxyz(i,2) = bxyz(ii,2)
!              bxyz(i,3) = bxyz(ii,3)
              I     = II
              II    = II + II
           ELSE
!C
!C-----> THIS PART OF THE TREE IS PROMOTED , STOP BY PUTTING II=IHEAP+1
              II    = IHEAP + 1
           END IF
                                                GO TO 40
        END IF
!C-----> PUT ELEMENTS IN THE PROPER SLOT
        brxyz(1,i) = brr
        brxyz(2:width,i) = bxyz(:)
!        bxyz(i,1) = bxx
!        bxyz(i,2) = byy
!        bxyz(i,3) = bzz
  30  continue
  end subroutine s_heapsort


real(dp) function exp_numeric_d8(x)
!! Written by Xiao Yan Xu (wanderxu@gmail.com)
!! calculate exponentail within an approximation method.
!! the relative error is controled by 'drange' in the following code
!! 
!! algorithm:
!!    first reduce x to rsex = x/2^nlog2 to make |rsex| <= drange
!!    then use the second order pade approximate equation: exp(z) ~ ( (z+3)^2+3 ) / ( (z-3)^2+3 )
!!    to calculate exp(rsex) and exp(x) will be (exp(rsex)^(2^nlog2)
  use constants, only : dp
  implicit none
  real(dp), intent(in) :: x
  real(dp), parameter :: ee = 2.7182818284590452353602874713527d0
  real(dp), parameter :: drange = 0.03d0 ! the relativ error is about 10^-8 for drange=0.03
  
  ! local
  integer :: nlog2, ilog2
  real(dp) :: rsex

  nlog2 = 0
  rsex = x
  do while ( abs(rsex)>drange )  ! control precision
      rsex = rsex / 2.d0
      nlog2 = nlog2 + 1
  end do

  exp_numeric_d8 = ( (rsex+3.d0)*(rsex+3.d0) + 3.d0 ) / ( (rsex-3.d0)*(rsex-3.d0) + 3.d0 )
  do ilog2 = 1, nlog2
      exp_numeric_d8 = exp_numeric_d8*exp_numeric_d8
  end do

end function exp_numeric_d8
