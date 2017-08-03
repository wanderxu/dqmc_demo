       program cov_tau
         
         use m_variance
         use constants,only: pi
         
         implicit real (kind=8) (a-g,o-z)
         implicit integer (h-n)
          

         complex (kind=8), dimension(:,:,:)  , allocatable :: g_bins
         complex (kind=8), dimension(:,:)    , allocatable :: z_mat, xcov
         complex (kind=8), dimension(:)      , allocatable :: xmean
         complex (kind=8)  :: z, zm, zerr
         real    (kind=8), dimension(:,:),  allocatable  :: x_k
         real    (kind=8), dimension(:)  ,  allocatable  :: sign
         character (16) :: file_k

         real(kind=8) :: xktmp(2)
         
         open ( unit=20, file='paramC_sets', status='unknown' )
            read(20,*) beta, ltrot, ntdm, l, flux_x, flux_y, nbin_skip
         close(20)

         !!!! WARNNING 
         !!!flux_x = -flux_x
         !!!flux_y = -flux_y
         

         norb = 1
         lq   = l*l
         dtau = beta/dble(ltrot)
         ! determine the number of bins. 
         open ( unit=10, file="intau", status="unknown" ) 
         nbins = 0
         do
            do n = 1,lq
               read(10,*,end=10)  x, y
               do nt = 1,ntdm+1
                  do no = 1,norb
                     do no1 = 1,norb
                        read(10,*) z
                     enddo
                  enddo
               enddo
            enddo
            nbins = nbins + 1
         enddo
10       continue
         write(6,*) "# of bins: ", nbins
         close(10) 

         nbins = nbins - nbin_skip


         ! allocate  space
         allocate ( g_bins(0:lq, ntdm+1, nbins)  )
         allocate ( x_k(lq,2), z_mat(norb,norb) )

         open ( unit=10, file="intau", status="unknown" ) 

         ! skip nskip bins
         do nb = 1, nbin_skip
            do n = 1,lq
               read(10,*)  x, y
               do nt = 1,ntdm+1
                  do no = 1,norb
                     do no1 = 1,norb
                        read(10,*) z
                     enddo
                  enddo
               enddo
            enddo
         enddo

         ! read-in the bins.
         do nb = 1,nbins
            do nk = 1,lq
               read(10,*)  x_k(nk,1), x_k(nk,2)
               do nt = 1,ntdm+1
                  do no = 1,norb
                     do no1 = 1,norb
                        read(10,*) z_mat(no,no1)
                     enddo
                  enddo
                  !g_bins(nk,nt,nb)  = ( z_mat(1,1) + z_mat(2,2) ) /cmplx(2.d0,0.d0)
                  g_bins(nk,nt,nb)  = z_mat(1,1)
               enddo
            enddo
         enddo
         close(10)

         do nt = 1,ntdm + 1 
            do nb = 1,nbins
               z = cmplx(0.d0,0.d0)
               do nk = 1,lq
                  z = z + g_bins(nk,nt,nb)
               enddo
               z = z/cmplx(lq,0.d0)
               g_bins ( 0,nt,nb ) = z
            enddo
         enddo
         
         
         allocate(sign(nbins))
         allocate(xcov(ntdm+1,ntdm+1), xmean(ntdm+1))
         sign = 1.d0
         do nk = 1,lq
            call cov(g_bins (nk,:,:), sign, xcov, xmean ) 

            ! itself
            xktmp(1) = x_k(nk,1)+flux_x*2.d0*pi/dble(l)
            xktmp(2) = x_k(nk,2)+flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

            ! reflection x=0, (-x,y)
            xktmp(1) = -x_k(nk,1)-flux_x*2.d0*pi/dble(l)
            xktmp(2) =  x_k(nk,2)+flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

            ! reflection y=0, (x,-y)
            xktmp(1) =  x_k(nk,1)+flux_x*2.d0*pi/dble(l)
            xktmp(2) = -x_k(nk,2)-flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

            ! reflection x=y, (y,x)
            xktmp(2) =  x_k(nk,1)+flux_x*2.d0*pi/dble(l)
            xktmp(1) =  x_k(nk,2)+flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

            ! reflection x=-y, (-y,-x)
            xktmp(2) = -x_k(nk,1)-flux_x*2.d0*pi/dble(l)
            xktmp(1) = -x_k(nk,2)-flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

            ! reflection x=0 and then reflection x=y, (y,-x)
            xktmp(2) = -x_k(nk,1)-flux_x*2.d0*pi/dble(l)
            xktmp(1) =  x_k(nk,2)+flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

            ! reflection y=0 and then reflection x=y, (-y,x)
            xktmp(2) =  x_k(nk,1)+flux_x*2.d0*pi/dble(l)
            xktmp(1) = -x_k(nk,2)-flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

            ! rotate pi, (-x,-y)
            xktmp(1) = -x_k(nk,1)-flux_x*2.d0*pi/dble(l)
            xktmp(2) = -x_k(nk,2)-flux_y*2.d0*pi/dble(l)
            if(xktmp(1) .gt. (pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)-2.d0*pi
            else if ( xktmp(1) .lt. (-pi-0.00001d0) ) then
                xktmp(1) = xktmp(1)+2.d0*pi
            end if
            if(xktmp(2) .gt. (pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)-2.d0*pi
            else if ( xktmp(2) .lt. (-pi-0.00001d0) ) then
                xktmp(2) = xktmp(2)+2.d0*pi
            end if
            nfile = nint((xktmp(2)+pi)*8.d0*dble(l)/(2.d0*pi))*8*l + nint((xktmp(1)+pi)*8.d0*dble(l)/(2.d0*pi)) + 1
           !nfile = nint((xktmp(2)+pi)*2.d0*dble(l)/(2.d0*pi))*2*l + nint((xktmp(1)+pi)*2.d0*dble(l)/(2.d0*pi)) + 1
            call files(nfile, file_k)
            open (unit=33,file=file_k,status="unknown")
            write(33,"(f12.6,2x,f12.6)") xktmp(1), xktmp(2)
            write(33,*) ntdm + 1 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
            close(33)

         enddo

       end program cov_tau
