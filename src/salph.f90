subroutine salph
  ! calculate dtau and alpha

  use blockc
  implicit none

  ! local 
  integer ::  m, n, id, it, inn_st, i
  real(dp) :: dth, sumjp, sumj, gssb
  integer :: sb_tau  , sb_taum1, sb_taua1, sbd1_tau, sbd2_tau, sbd3_tau, sbd4_tau, sb_taup
  if(irank.eq.0) then
      write(fout,*)
      write(fout,'(a)') ' In salph: '
  end if

  xsigp2(-1) = dcmplx( dexp(  dtau * rj * ( -1.d0 ) ), 0.d0 )  ! for positive eigenvalue,  -1.d0 is the field 
  xsigp2( 1) = dcmplx( dexp(  dtau * rj * (  1.d0 ) ), 0.d0 )  ! for positive eigenvalue,   1.d0 is the field
  xsigm2(-1) = dcmplx( dexp( -dtau * rj * ( -1.d0 ) ), 0.d0 )  ! for minux eigenvalue   ,  -1.d0 is the field
  xsigm2( 1) = dcmplx( dexp( -dtau * rj * (  1.d0 ) ), 0.d0 )  ! for minux eigenvalue   ,   1.d0 is the field

  if(irank.eq.0) then
      write(fout, '(a,2e16.8)') ' xsigp2(-1) = ', xsigp2(-1)
      write(fout, '(a,2e16.8)') ' xsigp2( 1) = ', xsigp2( 1)
      write(fout, '(a,2e16.8)') ' xsigm2(-1) = ', xsigm2(-1)
      write(fout, '(a,2e16.8)') ' xsigm2( 1) = ', xsigm2( 1)
  end if

  dellp2(-1, 1) = ( xsigp2( 1) / xsigp2(-1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from -1 to  1
  dellp2( 1, 1) = ( xsigp2(-1) / xsigp2( 1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from  1 to -1
  dellm2(-1, 1) = ( xsigm2( 1) / xsigm2(-1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from -1 to  1
  dellm2( 1, 1) = ( xsigm2(-1) / xsigm2( 1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from  1 to -1

  if(irank.eq.0) then
      write(fout, '(a,2e16.8)') ' dellp2(-1, 1) = ', dellp2(-1, 1)
      write(fout, '(a,2e16.8)') ' dellp2( 1, 1) = ', dellp2( 1, 1)
      write(fout, '(a,2e16.8)') ' dellm2(-1, 1) = ', dellm2(-1, 1)
      write(fout, '(a,2e16.8)') ' dellm2( 1, 1) = ', dellm2( 1, 1)
  end if
  
  !	for flipping of pair hopping-spin.
  nflipl(-1, 1) =  1
  nflipl( 1, 1) = -1
  
  ! pair-hopping. kinetic.
  ur_k(1,1) =   dcmplx( 1.d0/dsqrt(2.d0),0.d0 )
  ur_k(1,2) =   dcmplx( 1.d0/dsqrt(2.d0),0.d0 )
  ur_k(2,1) =   dcmplx( 1.d0/dsqrt(2.d0),0.d0 )
  ur_k(2,2) =   dcmplx(-1.d0/dsqrt(2.d0),0.d0 )
  
  do m = 1,2
  do n = 1,2
     urt_k(m,n) = dconjg(ur_k(n,m))
  enddo
  enddo
  
  ! pair-hopping.  current.
  ur_j(1,1) =    dcmplx( 1.d0/dsqrt(2.d0),0.d0 )
  ur_j(1,2) =    dcmplx( 1.d0/dsqrt(2.d0),0.d0 )
  ur_j(2,1) =    dcmplx( 0.d0,-1.d0/dsqrt(2.d0))
  ur_j(2,2) =    dcmplx( 0.d0, 1.d0/dsqrt(2.d0))
  
  do m = 1,2
  do n = 1,2
     urt_j(m,n) = dconjg(ur_j(n,m))
  enddo
  enddo

  ! onsite
  ! spin up, positive coupling
  xsigma_u_up(-1) = dcmplx( dexp(  dtau * rhub * ( -0.5d0 ) ), 0.d0 )  ! -1.d0 is the field, 0.5 is from fermion spin s=1/2 
  xsigma_u_up( 1) = dcmplx( dexp(  dtau * rhub * (  0.5d0 ) ), 0.d0 )  !  1.d0 is the field, 0.5 is from fermion spin s=1/2
  ! spin down, negative coupling
  xsigma_u_dn(-1) = dcmplx( dexp( -dtau * rhub * ( -0.5d0 ) ), 0.d0 )  ! -1.d0 is the field, 0.5 is from fermion spin s=1/2
  xsigma_u_dn( 1) = dcmplx( dexp( -dtau * rhub * (  0.5d0 ) ), 0.d0 )  !  1.d0 is the field, 0.5 is from fermion spin s=1/2

  delta_u_up(-1,1) = ( xsigma_u_up( 1) / xsigma_u_up(-1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from -1 to  1
  delta_u_up( 1,1) = ( xsigma_u_up(-1) / xsigma_u_up( 1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from  1 to -1

  delta_u_dn(-1,1) = ( xsigma_u_dn( 1) / xsigma_u_dn(-1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from -1 to  1
  delta_u_dn( 1,1) = ( xsigma_u_dn(-1) / xsigma_u_dn( 1) ) - dcmplx( 1.d0, 0.d0 )   ! flip from  1 to -1
  

  ! set wsxsz, the isping spin part of monte carlo ratio
  id = 0
  if( hx .lt. 1.e-6 ) then
      hx = 1.e-6
      if ( irank .eq.0 ) then
          write(50,'(a,e16.8)') ' alert!! the hx is too small, here we set it to ',hx,' to avoid nan '
      end if
  end if
  dth = dtau * hx
  tanhdth = dtanh(dth) ! set tanhdth, maybe useful in obser
  cothdth = 1.d0/dtanh(dth) ! set cothdth, maybe useful in obser

  if( irank .eq. 0 ) then
      write(50,*)
      write(50,'(a)') ' the wsxsz for transverse ising spin part mc ratio '
  end if

  do sb_tau   = -1, 1, 2
      sb_taup = nflipl( sb_tau, 1 )
  do sb_taum1 = -1, 1, 2
  do sb_taua1 = -1, 1, 2
  do sbd1_tau = -1, 1, 2
  do sbd2_tau = -1, 1, 2
  do sbd3_tau = -1, 1, 2
  do sbd4_tau = -1, 1, 2
      id = id + 1
      sumjp = dexp ( dble( sb_taup*sbd1_tau + sb_taup*sbd2_tau + sb_taup*sbd3_tau + sb_taup*sbd4_tau ) * dtau * js )
      sumj  = dexp ( dble( sb_tau *sbd1_tau + sb_tau *sbd2_tau + sb_tau *sbd3_tau + sb_tau *sbd4_tau ) * dtau * js )
      wsxsz( id ) = ( gssb( sb_taua1, sb_taup, dth ) * sumjp * gssb( sb_taup, sb_taum1, dth ) ) / & 
                    ( gssb( sb_taua1, sb_tau , dth ) * sumj  * gssb( sb_tau , sb_taum1, dth ) )
      !write(50,*)
      if ( irank .eq. 0 ) write(50,'(b7,e16.8)') id-1, wsxsz(id)
      !if ( sb_taua1 .eq. sb_taum1 .and. sb_tau .ne. sb_taua1 ) then
      !    write(50,'(a,e16.8)') ' sumjp = ', sumjp
      !    write(50,'(a,e16.8)') ' sumj  = ', sumj
      !end if
  end do
  end do
  end do
  end do
  end do
  end do
  end do

  ! set wjs
  id = 0
  if( irank .eq. 0 ) then
      write(50,*)
      write(50,'(a)') ' the wjs, the isping part for calculating global update mc ratio '
  end if
  do sb_tau   = -1, 1, 2
      sb_taup = nflipl( sb_tau, 1 )
  do sbd1_tau = -1, 1, 2
  do sbd2_tau = -1, 1, 2
  do sbd3_tau = -1, 1, 2
  do sbd4_tau = -1, 1, 2
      id = id + 1
      sumjp = dexp ( dble( sb_taup*sbd1_tau + sb_taup*sbd2_tau + sb_taup*sbd3_tau + sb_taup*sbd4_tau ) * dtau * js )
      sumj  = dexp ( dble( sb_tau *sbd1_tau + sb_tau *sbd2_tau + sb_tau *sbd3_tau + sb_tau *sbd4_tau ) * dtau * js )
      wjs( id ) = sumjp / sumj
      if ( irank .eq. 0 ) write(50,'(b5,e16.8)') id-1, wjs(id)
  end do
  end do
  end do
  end do
  end do


  ! set stbonds_neib ( space time neighbor ) for space time global update
  if( lstglobal ) then
      do it = 1, ltrot
      do i = 1, lq
          do inn_st = 1, 4
              stbonds_neib(1, inn_st, i, it ) = nnlist(i,inn_st)
              stbonds_neib(2, inn_st, i, it ) = it
          end do

          inn_st = 5  ! +t direction
          stbonds_neib(1, inn_st, i, it ) = i
          stbonds_neib(2, inn_st, i, it ) = mod(it,ltrot) + 1

          inn_st = 6  ! -t direction
          stbonds_neib(1, inn_st, i, it ) = i
          stbonds_neib(2, inn_st, i, it ) = mod(it+ltrot-2,ltrot) + 1

      end do
      end do

      ratio_nn_st(1:4) = 1.d0 - dexp(-2.d0*dtau*abs(js))
      ratio_nn_st(5:6) = 1.d0 - tanhdth

      if ( irank .eq. 0 ) then
          write( fout, * )
          write( fout, '(a)') '  i    ratio_nn_st '
          do i = 1, 6
              write( fout, '(i4, e16.8)' ) i, ratio_nn_st(i)
          end do
          write( fout, * )
      end if

   end if

end subroutine salph

function gssb( s1, s2, ralpha )
  integer :: s1, s2
  real(8) :: ralpha
  real(8) :: gssb
  if ( s1 .eq. s2 ) then
      gssb = dcosh( ralpha )
  else
      gssb = dsinh( ralpha )
  end if
  !write(50,'(a,e16.8)') ' gssb = ', gssb
  
  return
end function gssb
