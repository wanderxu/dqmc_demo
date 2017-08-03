       program cov_tau
         
         use m_variance
         
         implicit real (kind=8) (a-g,o-z)
         implicit integer (h-n)
          

         complex (kind=8), dimension(:,:,:)  , allocatable :: g_bins
         complex (kind=8), dimension(:,:)    , allocatable :: z_mat, xcov
         complex (kind=8), dimension(:)      , allocatable :: xmean
         complex (kind=8)  :: z, zm, zerr
         real    (kind=8), dimension(:,:),  allocatable  :: x_k
         real    (kind=8), dimension(:)  ,  allocatable  :: sign
         character (16) :: file_k
         
         open ( unit=20, file='paramC_sets', status='unknown' )
            read(20,*) beta, ltrot, ntdm, l
         close(20)
         

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


         ! allocate  space
         allocate ( g_bins(0:lq, ntdm+1, nbins)  )
         allocate ( x_k(lq,2), z_mat(norb,norb) )

         ! read-in the bins.
         open ( unit=10, file="intau", status="unknown" ) 
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
         do nk = 0,lq
            call files(nk, file_k  )
            open (unit=33,file=file_k,status="unknown")
            if (nk.gt.0) write(33,"(f12.6,2x,f12.6)") x_k(nk,1), x_k(nk,2)
            write(33,*) ntdm + 1 
            call cov(g_bins (nk,:,:), sign, xcov, xmean ) 
            if( ntdm .eq. 0 ) then
                write(33,"(f14.7,2x,f16.8,2x,f16.8)") beta/2.d0,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
            else
                do nt = 1,ntdm+1
                   write(33,"(f14.7,2x,f16.8,2x,f16.8)") dble(nt)*dtau,  dble(xmean(nt)), sqrt(abs(dble(xcov(nt,nt)))) 
                enddo
            end if
!!            do nt = 1,ntdm+1
!!               do nt1 = 1,ntdm+1
!!                  write(33,*) dble(xcov(nt,nt1))
!!               enddo
!!            enddo
         enddo
         close(33)


       end program cov_tau
