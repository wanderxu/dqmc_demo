subroutine upgradej(ntau,nf,green_up,green_dn)

  use spring
  use blockc
  use data_tmp
  
  implicit none
  
  !arguments
  integer,intent(in) :: ntau, nf
  complex(dp), dimension(ndim,ndim) :: green_up, green_dn
  !complex(dp) ::  phasej
      
  !local
  complex(dp) ::  g44up, g55up, g45up, g54up, ratioup, ratiotot, del44, del55, z1, z2, z3, z4
  integer ::  nf1, nn, nrflip, i, i1, i4, i5, nl, nl1, nl2, nfb, id, ntm1, nta1
  real(dp) :: accm, ratio_re, ratio_re_abs, random, weight
  
  if (nf.gt.4) then ! current.
     nf1 = nf -4
  else              ! kinetic.
     nf1 = nf
  endif
  if (nf1.eq.1) nn = 1
  if (nf1.eq.2) nn = 1
  if (nf1.eq.3) nn = 2
  if (nf1.eq.4) nn = 2
  
  accm = 0.d0
  do i = 1,lq
     i1 = i
     i4 = i1
     i5 = i1 + lq
     nrflip = 1
     if (nf.gt.4) then ! current.
        del44   =  dellp2( nsigl_j(i1,nn,ntau), nrflip )
        del55   =  dellm2( nsigl_j(i1,nn,ntau), nrflip )
     else              ! kenitic
        del44   =  dellp2( nsigl_k(i1,nn,ntau), nrflip )
        del55   =  dellm2( nsigl_k(i1,nn,ntau), nrflip )
     endif
     g44up = dcmplx(0.d0,0.d0)
     g45up = dcmplx(0.d0,0.d0)
     g54up = dcmplx(0.d0,0.d0)
     g55up = dcmplx(0.d0,0.d0)
  
     g44up =   del44 * ( cone - green_up(i4,i4) )
     g45up = - del44 * green_up( i4, i5 )
     g54up = - del55 * green_up( i5, i4 )
     g55up =   del55 * ( cone - green_up(i5,i5) )
  
     ratioup = (dcmplx(1.d0,0.d0) + g44up) * (dcmplx(1.d0,0.d0) + g55up) - g45up*g54up

#IFDEF TEST
     write(fout,'(a,2e16.8)') 'in upgradej, ratioup = ', ratioup
#ENDIF

     ! get Ising part ratio
     id = 0
     ntm1 = ntau - 1
     if ( ntm1 .lt. 1 ) ntm1 = ntm1 + ltrot
     nta1 = ntau + 1
     if ( nta1 .gt. ltrot ) nta1 = nta1 - ltrot
     if( nsigl_k( i, nf, ntau ) .eq. 1 ) id = ibset( id, 6 )
     if( nsigl_k( i, nf, ntm1 ) .eq. 1 ) id = ibset( id, 5 )
     if( nsigl_k( i, nf, nta1 ) .eq. 1 ) id = ibset( id, 4 )
     do nfb = 1, 4
         if( nsigl_k( nnlist(i1,nfb), nn, ntau ) .eq. 1 ) id = ibset( id, 4-nfb )
     end do
     id = id + 1
  
     ! total ratio
     ratiotot = ratioup*dconjg(ratioup) * wsxsz(id)
     ratio_re = dble( ratiotot )

#IFDEF TEST
     write(fout,'(a,e16.8)') 'in upgradej, ratio_re = ', ratio_re
#ENDIF
     
     ratio_re_abs = ratio_re
     if (ratio_re .lt. 0.d0 )  ratio_re_abs = - ratio_re 
     ! write(6,*) 'upgrade phasej: ', z
  
     random = spring_sfmt_stream()
     if (ratio_re_abs.gt.random) then
  
        ! write(50,*) 'accepted'
        ! upgrade the inverse.
  
        accm = accm + 1.d0
        weight = dsqrt(dble(ratiotot*dconjg(ratiotot)))
        !phasej =  phasej*ratiotot/dcmplx(weight,0.d0)

        z1 = cone / ( cone + g44up )
        z2 = g45up * z1
        z3 = g54up * z1
        z4 = cone + g55up - g45up*g54up*z1
        z4 = cone / z4

        ! v1(:) = ( 1 - G ) (i4, :)
        ! v2(:) = ( 1 - G ) (i5, :)
        do nl = 1, ndim
            u1(nl) = - del44 * green_up(i4,nl)
            u2(nl) = - del55 * green_up(i5,nl)
        end do
        u1(i4) = del44 + u1(i4)
        u2(i5) = del55 + u2(i5)

        do nl = 1, ndim
            uhlp1(nl) = green_up(nl,i4)
            vhlp1(nl) = z1 * u1(nl)
        end do

        do nl =1, ndim
            uhlp2(nl) = green_up(nl,i5) - green_up(nl,i4) * z2
            vhlp2(nl) = z4 * ( u2(nl) - u1(nl) * z3 )
        end do
  
        do nl2 = 1,ndim
        do nl1 = 1,ndim
           green_up(nl1,nl2) = green_up(nl1,nl2) - uhlp1(nl1)*vhlp1(nl2) - uhlp2(nl1)*vhlp2(nl2)
        enddo
        enddo
            !	      flip:
        if (nf.gt.4) then !                current.
  	     nsigl_j(i1,nn,ntau) =    nflipl(nsigl_j(i1,nn,ntau), nrflip)
        else !                kenitic
  	     nsigl_k(i1,nn,ntau) =    nflipl(nsigl_k(i1,nn,ntau), nrflip)
        endif
        
     endif
     
  enddo
  main_obs(2) = main_obs(2) + dcmplx( accm, dble(lq) )
end subroutine upgradej
