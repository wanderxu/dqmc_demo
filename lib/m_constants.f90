  module constants
     implicit none
     integer, public, parameter :: sp    = kind(1.0) ! single precision
     integer, public, parameter :: dp    = kind(1.0d0) ! double precision
     real(dp), public, parameter :: pi   = 3.141592653589793238462643383279_dp
     real(dp), public, parameter :: zero = 0.0_dp
     real(dp), public, parameter :: one  = 1.0_dp
     real(dp), public, parameter :: two  = 2.0_dp
     real(dp), public, parameter :: half = 0.5_dp
     real(dp), public, parameter :: eps6 = 1.0E-6
     real(dp), public, parameter :: eps8 = 1.0E-8
     real(dp), public, parameter :: epst = 1.0E-10
     real(dp), public, parameter :: epss = 1.0E-12
     complex(dp), public, parameter :: czi   = dcmplx(0.0_dp, 1.0_dp)
     complex(dp), public, parameter :: cone  = dcmplx(1.0_dp, 0.0_dp)
     complex(dp), public, parameter :: czero = dcmplx(0.0_dp, 0.0_dp)
  end module constants
