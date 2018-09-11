module m_variance
  interface errcalc
     module procedure errcalc, errcalc_c
  end interface
  interface errcalcj
     module procedure errcalc_j, errcalc_j_rebin,  errcalc_js, errcalc_js_rebin, &
          &   errcalc_j_c, errcalc_js_c
  end interface
  interface cov
     module procedure covj, covjs, covjs_c
  end interface
  interface cov_err
     module procedure cov_err
  end interface
  interface intergrate_f
     module procedure inter_f
  end interface
  interface intergrate
     module procedure inter_qmc
  end interface
  interface fit
     module procedure fit
  end interface
  interface auto_cor
     module procedure  auto_cor
  end interface
  interface bootstrap
     module procedure  bootstrap
  end interface

  

  contains 

    subroutine errcalc(en,xm,xerr)
    ! calculates error on the input vector en.  just the standard deviation.
    ! sigma^2 = sum_i ( x_i - x_mean )^2 / (n-1) = n/(n-1) * ( xsq_mean - x_mean^2 )

      implicit none
      real (kind=8), dimension(:) ::  en
      real (kind=8)               ::  xm, xerr, xsq
      integer                     ::  np, nt

      np = size(en)
      
      xm  = 0.d0
      xsq = 0.d0
      do nt = 1,np
         xm  = xm  + en(nt)
         xsq = xsq + en(nt)**2
      enddo
      xm    = xm /dble(np)
      xsq   = xsq/dble(np)

      if( np .gt. 1 ) then
        xerr  = (xsq - xm**2) * dble(np) /dble(np-1)
      else
        xerr = 0.d0
      end if

      if (xerr.gt.0.d0) then
         xerr = sqrt(xerr)
      else
         xerr = 0.d0 
      endif
      
      return
    end subroutine errcalc


    subroutine errcalc_c(en,zm,zerr)
    ! calculates error on the input vector en.  just the standard deviation.

      implicit none
      complex (kind=8), dimension(:) ::  en
      complex (kind=8)               ::  zm, zerr
      integer                     ::  np, nt

      ! local 
      real (kind=8), dimension(:), allocatable :: rhelp
      real (kind=8) :: xm, xerr

      np = size(en)
      allocate (rhelp(np))
      
      do nt = 1,np
         rhelp(nt) = dble(en(nt))
      enddo
      call errcalc(rhelp, xm, xerr)
      zm   =  dcmplx(xm  , 0.d0)
      zerr =  dcmplx(xerr, 0.d0)

      do nt = 1,np
         rhelp(nt) = aimag(en(nt))
      enddo
      call errcalc(rhelp, xm, xerr)
      zm   =  zm   + dcmplx( 0.d0, xm   )
      zerr =  zerr + dcmplx( 0.d0, xerr )

      return
    end subroutine errcalc_c

    subroutine errcalc_j(en,xm,xerr)
      !   calculates jacknife error on the input vector en.  mean and  variance.
      !   the input are the bins.
      
      implicit none

      real (kind=8), dimension(:) ::  en
      real (kind=8)               ::  xm, xerr, x
      real (kind=8), dimension(:), allocatable ::  en1
      integer     :: np, n, n1

      np = size(en)
      allocate (en1(np))
      
      ! build the jackknife averages and send to errcalc.

      do n = 1,np
         x = 0.d0
         do n1 = 1,np
            if (n1.ne.n) x = x + en(n1)
         enddo
         en1(n) = x / dble(np -1)
      enddo
      call errcalc(en1,xm,xerr)
      deallocate  ( en1 )

      return
    end subroutine errcalc_j


    subroutine errcalc_j_c(en,zm,zerr)
      ! calculates jacknife error on the input vector en.  mean and  variance.
      ! the input are the bins.
      
      implicit none

      complex (kind=8), dimension(:) ::  en
      complex (kind=8)               ::  zm, zerr, z
      complex (kind=8), dimension(:), allocatable ::  en1
      integer     :: np, n, n1

      np = size(en)
      allocate (en1(np))
      
      ! build the jackknife averages and send to errcalc.

      do n = 1,np
         z = cmplx(0.d0, 0.d0)
         do n1 = 1,np
            if (n1.ne.n) z = z + en(n1)
         enddo
         en1(n) = z / cmplx(dble(np -1) , 0.d0)
      enddo
      call errcalc(en1,zm,zerr)
      deallocate  ( en1 )

      return
    end subroutine errcalc_j_c


    subroutine errcalc_j_rebin(en,xm,xerr,nrebin)
      ! calculates jacknife error on the input vector en with rebinning.  mean and  variance.
      ! the input are the bins.

      implicit none

      real (kind=8), dimension(:) ::  en
      real (kind=8)               ::  xm, xerr, x
      real (kind=8), dimension(:), allocatable ::  en1
      integer :: nrebin, nc, n, nb, np1, np

      np = size(en)
      np1 = np/nrebin
      allocate (en1(np1))
      
      ! rebin
      nc = 0
      do n = 1,np1
         x = 0.d0
         do nb = 1,nrebin
            nc = nc + 1
            x = x + en(nc)
         enddo
         x = x/dble(nrebin)
         en1(n) = x
      enddo
      call errcalc_j(en1,xm,xerr)
      xerr = xerr / sqrt( dble( np1 ) ) ! error = standard variance / sqrt(nbin)

      deallocate(en1)
      return
    end subroutine errcalc_j_rebin

    subroutine errcalc_js(en,si,xm,xerr)
      ! calculates error on the input vector en.  just the variance.
      ! the input are the bins

      implicit none

      real (kind=8), dimension(:) ::  en, si
      real (kind=8)               ::  xm, xerr, x,xs
      real (kind=8), dimension(:), allocatable :: en1
      integer                     ::  n, n1, np, np1

      np = size(en)
      np1= size(si)
      if (np1.ne.np) then
         write(6,*) 'error in errcalc_js'
         stop
      endif
      allocate (en1(np))
      
      ! build the jackknife averages and send to errcalc

      do n = 1,np
         x  = 0.d0
         xs = 0.d0
         do n1 = 1,np
            if (n1.ne.n)  x  = x  + en(n1)
            if (n1.ne.n)  xs = xs + si(n1)
         enddo
         en1(n) = x / xs
      enddo
      call errcalc(en1,xm,xerr)
      deallocate  ( en1 )

      return
    end subroutine errcalc_js

    subroutine errcalc_js_c(en,si,xm,xerr)
      ! calculates error on the input vector en.  just the variance.
      ! the input are the bins

      implicit none

      complex (kind=8), dimension(:) ::  en, si
      complex (kind=8)               ::  xm, xerr, x,xs
      complex (kind=8), dimension(:), allocatable :: en1
      integer                     ::  n, n1, np, np1

      np = size(en)
      np1= size(si)
      if (np1.ne.np) then
         write(6,*) 'error in errcalc_js'
         stop
      endif
      allocate (en1(np))
      
      ! build the jackknife averages and send to errcalc

      do n = 1,np
         x  = cmplx(0.d0,0.d0)
         xs = cmplx(0.d0,0.d0)
         do n1 = 1,np
            if (n1.ne.n)  x  = x  + en(n1)
            if (n1.ne.n)  xs = xs + si(n1)
         enddo
         en1(n) = x / xs
      enddo
      call errcalc(en1,xm,xerr)
      deallocate  ( en1 )

      return
    end subroutine errcalc_js_c



    subroutine errcalc_js_rebin(en,si,xm,xerr,nrebin)
      ! calculates jacknife error on the input vector en with rebinning.  mean and  variance.
      ! the input are the bins.
      
      implicit none

      real (kind=8), dimension(:) ::  en, si
      real (kind=8)               ::  xm, xerr, x, y
      real (kind=8), dimension(:), allocatable ::  en1, si1
      integer :: nrebin, nc, n, nb, np, np1

      np = size(en)
      np1 = np/nrebin
      allocate (en1(np1))
      allocate (si1(np1))
      
      ! rebin
      nc = 0
      do n = 1,np1
         x = 0.d0; y = 0.d0
         do nb = 1,nrebin
            nc = nc + 1
            x = x + en(nc)
            y = y + si(nc)
         enddo
         x = x/dble(nrebin)
         y = y/dble(nrebin)
         en1(n) = x
         si1(n) = y
      enddo
      call errcalc_js(en1,si1,xm,xerr)
      xerr = xerr / sqrt( dble( np1 ) ) ! error = standard variance / sqrt(nbin)

      deallocate (en1,si1)

      return
    end subroutine errcalc_js_rebin

    subroutine inter_qmc(gr, sign1, dtau, res, err) 
      
      implicit none
      ! given gr(times, bins)  and sign1(bins) calculates the integral and error
      ! the sign is the same for all times.
      real (kind=8), dimension(:,:) ::  gr
      real (kind=8), dimension(:)   ::   sign1

      !local
      real (kind=8), dimension(:  ), allocatable  ::  hlp
      real (kind=8), dimension(:,:), allocatable  ::  hlp1
      real (kind=8)                 ::  x, xm, xerr, y, err, res, dtau
      integer :: nt, nt1, nb, nb1, ntdm, ndata
      
      ntdm  = size(gr,1)
      ndata = size(gr,2)


      allocate( hlp(ndata), hlp1(ntdm,ndata) )
      do nt = 1,ntdm
         do nb= 1, ndata
            x = 0.d0
            y = 0.d0
            do nb1 = 1,ndata
               if (nb1.ne.nb) then
                  x = x + gr(nt,nb1)
                  y = y + sign1(nb1)
               endif
            enddo
            hlp1(nt,nb) = x/y
         enddo
      enddo
      
      do nb = 1,ndata
         x = 0.d0
         do nt = 1,ntdm-1
            x = x + (hlp1(nt,nb) + hlp1(nt+1,nb))*0.5d0
         enddo
         hlp (nb   )  = x * dtau
      enddo
      
      call errcalc(hlp, res, err) 
      err = err*dble(ndata)

      deallocate( hlp, hlp1 )

      return
    end subroutine inter_qmc
    
    real (kind=8) function inter_f(a,b,n,f)
      ! integrates the function f from a to b  using n points.
      
      implicit none

      integer  :: n, i
      real (kind=8) ::  a, b, res, x, x1
      real (kind=8), external :: f

      real (kind=8) ::  del

      del = (b-a)/dble(n)
      inter_f = 0.d0
      do i = 0, n-1
         x  = a + dble(i  )*del
         x1 = a + dble(i+1)*del 
         inter_f = inter_f + ( f(x) + f(x1) )*0.5d0  
      enddo
      inter_f = inter_f*del
    end function inter_f

    ! least square fits:
    subroutine fit(xdata,fdata,error,ares,chsq,f,aerr)

      implicit none

      real(kind=8), dimension(:) ::  xdata, fdata, error, ares
      real(kind=8), dimension(:), optional :: aerr
      real(kind=8)               ::  chsq,  x
      real(kind=8), dimension(:,:),  allocatable :: a, u,v,vinv,v1
      real(kind=8), dimension(:  ),  allocatable :: b,d
      real(kind=8), external  :: f
      integer                     ::  ndata, nbasis, i, m, m1, ncon, n

      ndata = size(xdata)
      nbasis= size(ares)

      !write(6,*) 'ndata, nbasis: ',ndata, nbasis 
      allocate (a(ndata,nbasis))
      allocate (u(ndata,nbasis))
      allocate (d(nbasis))
      allocate (v   (nbasis,nbasis))
      allocate (v1  (nbasis,nbasis))
      allocate (vinv(nbasis,nbasis))
      allocate (b(ndata))
      
      a = 0.d0
      u = 0.d0
      d = 0.d0
      v = 0.d0
      vinv = 0.d0
      v1 = 0.d0
      b = 0.d0
      ncon = 1
      do m = 1,nbasis
         do i = 1,ndata
            a(i,m) = f(m,xdata(i))/error(i)
         enddo
      enddo
      do i = 1,ndata
         b(i) = fdata(i)/error(i)
      enddo
      !write(6,*) a
      !!!call udv(a,u,d,v,ncon)
      call s_svd_dg(ndata,nbasis,nbasis,u,d,v)
      write(6,*) ' WARNING, we assume nbasis < ndata in subroutine fit of m_variance module '
      do m = 1,nbasis
         do i = 1,nbasis
            v1(i,m) = v(m,i)
         enddo
      enddo
      x = 0.d0
      vinv = v1
      call s_inv_d(nbasis,vinv)
      !!!call inv(v1,vinv,x)

      do m1 = 1,nbasis
         x = 0.d0
         do m = 1,nbasis
            do i = 1,ndata
               x = x + b(i)*u(i,m)*vinv(m,m1)/d(m)
            enddo
         enddo
         ares(m1) = x
      enddo

        if (present(aerr)) then
           do m1 = 1,nbasis
              x = 0.d0
              do m = 1,nbasis
                 x = x + (vinv(m1,m)/d(m))**2
              enddo
              aerr(m1) = sqrt(x)
           enddo
        end if

      chsq = 0.d0
      do n = 1,ndata
         x = 0.d0
         do m = 1,nbasis
            x = x + ares(m)*f(m,xdata(n))
         enddo
         chsq = chsq + (fdata(n) - x)**2/error(n)**2
      enddo
      chsq = chsq/dble(ndata)
      
      deallocate (a)
      deallocate (u)
      deallocate (d)
      deallocate (v)
      deallocate (v1)
      deallocate (vinv)
      deallocate (b)
      
    end subroutine fit

    subroutine covj(gr, xcov, xmean) 
      
      implicit none
      !given gr(times, bins)  calculates the mean and the covariance.
      real (kind=8), dimension(:,:) ::  gr, xcov
      real (kind=8), dimension(:)   ::  xmean

      !local
      real (kind=8), dimension(:  ), allocatable  ::  hlp
      real (kind=8), dimension(:,:), allocatable  ::  hlp1
      real (kind=8)                 ::  x, xm, xerr
      integer :: nt, nt1, nb, nb1, ntdm, ndata
      
      ntdm  = size(gr,1)
      ndata = size(gr,2)

      if ( (size(xcov,1).ne.size(xcov,2) ) .or. (size(xcov,1).ne.ntdm) ) then
         write(6,*) 'error in cov'
         stop
      endif

      allocate( hlp(ndata), hlp1(ntdm,ndata) )
      do nt = 1,ntdm
         do nb= 1, ndata
            x = 0.0
            do nb1 = 1,ndata
               if (nb1.ne.nb) then
                  x = x + gr(nt,nb1)
               endif
            enddo
            hlp1(nt,nb) = x/dble(ndata-1)
            hlp (nb   ) = x/dble(ndata-1)
         enddo
         call errcalc(hlp,xm ,xerr)
         xmean(nt) = xm
      enddo

      
      do nt = 1,ntdm
         do nt1= 1,ntdm
            x = 0.0
            do nb = 1,ndata
               x = x +  hlp1(nt,nb)*hlp1(nt1,nb)
            enddo
            x = x/dble(ndata)
            xcov(nt,nt1)  = ( x - xmean(nt)*xmean(nt1) )*dble(ndata)
         enddo
      enddo
      

      deallocate( hlp, hlp1 )

      return
    end subroutine covj


    subroutine covjs(gr, sign1, xcov, xmean) 
      
      implicit none
      ! given gr(times, bins)  and sign1(bins) calculates the mean and the covariance.
      ! the sign is the same for all times.
      real (kind=8), dimension(:,:) ::  gr, xcov
      real (kind=8), dimension(:)   ::  xmean, sign1

      !local
      real (kind=8), dimension(:  ), allocatable  ::  hlp
      real (kind=8), dimension(:,:), allocatable  ::  hlp1
      real (kind=8)                 ::  x, xm, xerr, y
      integer :: nt, nt1, nb, nb1, ntdm, ndata
      
      ntdm  = size(gr,1)
      ndata = size(gr,2)

      if ( (size(xcov,1).ne.size(xcov,2) ) .or. (size(xcov,1).ne.ntdm) ) then
         write(6,*) 'error in cov'
         stop
      endif

      allocate( hlp(ndata), hlp1(ntdm,ndata) )
      do nt = 1,ntdm
         do nb= 1, ndata
            x = 0.d0
            y = 0.d0
            do nb1 = 1,ndata
               if (nb1.ne.nb) then
                  x = x + gr(nt,nb1)
                  y = y + sign1(nb1)
               endif
            enddo
            hlp1(nt,nb) = x/y
            hlp (nb   ) = x/y
         enddo
         call errcalc(hlp,xm ,xerr)
         xmean(nt) = xm
      enddo

      
      do nt = 1,ntdm
         do nt1= 1,ntdm
            x = 0.0
            do nb = 1,ndata
               x = x +  hlp1(nt,nb)*hlp1(nt1,nb)
            enddo
            x = x/dble(ndata)
            xcov(nt,nt1)  = ( x - xmean(nt)*xmean(nt1) )*dble(ndata)
         enddo
      enddo
      

      deallocate( hlp, hlp1 )

      return
    end subroutine covjs




    subroutine covjs_c(gr, sign1, xcov, xmean) 
      
      implicit none
      ! given gr(times, bins)  and sign1(bins) calculates the mean and the covariance.
      ! the sign is the same for all times. 
      complex (kind=8), dimension(:,:) ::  gr, xcov
      complex (kind=8), dimension(:)   ::  xmean
      real    (kind=8), dimension(:)   ::  sign1


      !local
      real (kind=8), dimension(:  ), allocatable  ::  hlp, xmean_r
      real (kind=8), dimension(:,:), allocatable  ::  hlp1
      real (kind=8)                 ::  x, xm, xerr, y
      integer :: nt, nt1, nb, nb1, ntdm, ndata, nth
      complex (kind=8) :: z

      ntdm  = size(gr,1)
      ndata = size(gr,2)


      !write(6,*) 'errors.f90 ', ntdm, ndata
      if ( (size(xcov,1).ne.size(xcov,2) ) .or. (size(xcov,1).ne.ntdm) ) then
         write(6,*) 'error in cov'
         stop
      endif


      
      allocate( hlp(ndata), hlp1(ntdm,ndata), xmean_r(ntdm) )
      xmean = cmplx(0.d0,0.d0)
      xcov  = cmplx(0.d0,0.d0)
      
      do nth = 1,2
         z = cmplx(1.0, 0.0) 
         if (nth .eq. 2 ) z = cmplx( 0.0, -1.0 )
         do nt = 1,ntdm
            do nb= 1, ndata
               x = 0.d0
               y = 0.d0
               do nb1 = 1,ndata
                  if (nb1.ne.nb) then
                     x = x + dble ( z*gr(nt,nb1) ) 
                     y = y + sign1(nb1)
                  endif
               enddo
               hlp1(nt,nb) = x/y
               hlp (nb   ) = x/y
            enddo
            call errcalc(hlp,xm ,xerr)
            xmean(nt) = xmean(nt) + conjg(z)*cmplx(xm,0.d0)
            xmean_r(nt) = xm
            !if (nth.eq.2) write(6,*) xm
         enddo
         
         
         do nt = 1,ntdm
            do nt1= 1,ntdm
               x = 0.0
               do nb = 1,ndata
                  x = x +  hlp1(nt,nb)*hlp1(nt1,nb)
               enddo
               x = x/dble(ndata)
               xcov(nt,nt1)  = xcov(nt,nt1) + conjg(z)* &
                    &    cmplx( ( x - xmean_r(nt)*xmean_r(nt1) )*dble(ndata) , 0.d0 )
            enddo
         enddo
      enddo

      deallocate( hlp, hlp1, xmean_r )
      
      return
    end subroutine covjs_c
    



    subroutine cov_err(xmean, xcov) 
      !  given mean and cov, diagonalizes the cov and produces a new data set within 
      !  the errorbars

      use spring, only: spring_sfmt_gaussian
      
      implicit none
      ! parameters
      real (kind=8), dimension(:,:) :: xcov
      real (kind=8), dimension(:)   :: xmean 
      
      integer :: ntau, i, m
      real (kind = 8) :: x, rang

      real (kind=8), dimension(:,:),  allocatable ::  uc
      real (kind=8), dimension(:),    allocatable ::  xmean_1, sig_1

      ntau = size(xmean,1) 
      allocate (uc(ntau,ntau), xmean_1(ntau), sig_1(ntau) ) 

      !!!call diag(xcov,uc,sig_1)
      call s_eig_dg(ntau,ntau,xcov,sig_1,uc)
      
      do i = 1,ntau
         x = 0.d0
         do m = 1,ntau
            x  = x + uc(m,i)* xmean(m)
         enddo
         xmean_1(i) = x
      enddo
      do i = 1,ntau
     if (sig_1(i).lt.0.d0) then 
         write(6,*) 'error in cov_err', sig_1(i)
     endif
         xmean_1(i) = xmean_1(i) + spring_sfmt_gaussian(0.d0,sqrt(abs(sig_1(i))))
      enddo
      do i = 1,ntau
         x = 0.d0
         do m = 1,ntau
            x  = x + uc(i,m)*xmean_1(m)
         enddo
         xmean(i) = x
      enddo

      deallocate (uc, xmean_1, sig_1) 


    end subroutine cov_err
    
    subroutine  auto_cor(data,res)

      implicit none
      
      real (kind=8),  dimension(:)  :: data,res

      !local 
      integer  :: nb, nt, ntau, nt1
      real (kind=8) :: x1, x2, x3
     
      nb = size(data)
      nt = size(res) 
      if (nb.lt.nt) then
         write(6,*) 'error in autocor'
         stop
      end if
      
      do ntau = 1,  nt
         x1 = 0.0
         x2 = 0.0
         x3 = 0.0
         do nt1 = 1, nb - ntau
            x1 = x1 + data(nt1)*data(nt1 + ntau) 
            x2 = x2 + data(nt1)*data(nt1)
            x3 = x3 + data(nt1)
         enddo
         x1 = x1 / dble(nb - ntau)
         x2 = x2 / dble(nb - ntau)
         x3 = x3 / dble(nb - ntau)
         
         res(ntau)  = ( x1 - x3**2)/(x2 - x3**2)
         
      enddo
      
    end subroutine auto_cor


    subroutine bootstrap(en,xm,xerr,nboot)

      use spring, only: spring_sfmt_string

      implicit none
      real (kind=8), dimension(:) ::  en
      real (kind=8)               ::  xm, xerr,  x, ranf
      integer                     ::  np, nt, nboot, nb, i

      np = size(en)
      
      ! build the bootstrap samples
      
      xm   = 0.d0
      xerr = 0.d0
      do nb = 1,nboot
         x = 0.d0
         do nt = 1, np
            i = nint( dble(np)* spring_sfmt_string() + 0.5 )
            if (i.eq.0 .or. i.gt.np ) then
               write(6,*) 'error in bootstrap'
               stop
            endif
            x = x + en(i)
         enddo
         x = x/dble(np)
         xm   = xm + x
         xerr = xerr + x*x
      enddo

      xm   = xm  /dble(nboot)
      xerr = xerr/dble(nboot)

      x = xerr - xm*xm
      xerr = 0.d0
      if (x.gt.0.d0) xerr = sqrt(x)

    end subroutine bootstrap

end module m_variance
