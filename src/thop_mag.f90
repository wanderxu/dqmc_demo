! do not forgett to declare zthp as complex in calling programm.
function zthp(i,nax,nay,xmag,flux_x,flux_y)
  use blockc, only:dp, pi, l, lq, list, dimer
  implicit none
  complex(dp) :: zthp
  integer, intent(in) :: i, nax, nay
  real(dp), intent(in) :: xmag, flux_x, flux_y

  !  local 
  real(dp) :: x, x1, xmag1

  ! xmag is magnetic xmag per plaquette.
  ! flux  is the twisting of boundary condition in  x-direction.
  ! both xmag and flux are in units of flux quantum.
  
  !write(6,*) 'lq in thop_mag: ', lq


  !     uses landau gauge to compute the matix element 
  !     c^{dagger}_i c_j exp(2 pi i / phi_0 \int_i^j a dl),  j = i + (nax,nay)
  !     a(x) = -b(x_2,0,0) with bondary conditions.
  !     i_x, i_y in [1,l]

  xmag1  = xmag/dble(lq)   ! for thermerdynamic limit
  
  x = -2.0*pi * xmag1 * dble( nax ) * dble(list(i,2) )
  
  x1 = 0.0
  if ( list(i,2) .eq. l  .and. nay .gt. 0)  then
     x1 = 2.0*pi * xmag1 * dble(l) *dble(list(i,1)) 
  endif
  
  if ( list(i,2) .eq. 1  .and. nay .lt. 0)  then
     x1 = -2.0*pi * xmag1 * dble(l) *dble(list(i,1)) 
  endif
  
  zthp = exp( dcmplx(0.d0, x + x1) )
  
  !     flux.
  if (nax.eq.1 .and. nay.eq.0) then
     zthp  = zthp*exp(dcmplx(0.d0, 2.d0*pi*flux_x/dble(l))) 
  endif
  if (nax.eq.0 .and. nay.eq.1) then
     zthp  = zthp*exp(dcmplx(0.d0, 2.d0*pi*flux_y/dble(l))) 
  endif
  
  !     dimerization.
  if (nax.eq.1 .and. nay.eq.0) then
     if ( mod(list(i,1),2).eq.0 ) then
        !             write(6,*) 'dimerize: ', list(nc,1)
        zthp = zthp*dcmplx(1.d0-dimer,0.d0)
     else
        zthp = zthp*dcmplx(1.d0+dimer,0.d0)
     endif
  endif
  if (nax.eq.0 .and. nay.eq.1) then
     zthp = zthp*dcmplx(1.d0-dimer,0.d0)
  endif
  
end function zthp
