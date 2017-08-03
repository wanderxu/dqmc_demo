       Program Cov_eq
         
         Use Errors
         Use MyMats
	 Use Lattices_v3

         Implicit Real (KIND=8) (A-G,O-Z)
         Implicit Integer (H-N)

	 TYPE(Lattice) :: Latt
         Complex (Kind=8),  Dimension(:,:,:,:),  Allocatable :: g_bins_k, g_bins_r
         Complex (Kind=8)  :: Z, ZM, ZERR
         Real    (Kind=8), Dimension(:,:),   allocatable  :: X_K
	 Real(Kind=8), Dimension(2) :: a1_p, a2_p, L1_p, L2_p
	 Real(Kind=8), DIMENSION(2) :: xk_p, ir_p
         
         OPEN ( UNIT=20, FILE='paramC_sets', STATUS='UNKNOWN' )
            READ(20,*) BETA, LTROT, NWRAP, RHUB, XL, RT3, NBIN, NSWEEP, LTAU, NTDM
            READ(20,*) L ,  TwistX             
         Close(20)

         Norb = 2
         LQ   = L*L
         ! Determine the number of bins. 
         Open ( Unit=10, File="ineq", status="unknown" ) 
         nbins = 0
         do
            do n = 1,LQ
               read(10,*,End=10)  X, Y
               do no = 1,Norb
                  do no1 = 1,Norb
                     read(10,*) Z
                  enddo
               enddo
            enddo
            nbins = nbins + 1
         enddo
10       continue
         Write(6,*) "# of bins: ", Nbins
         Close(10) 


         ! Allocate  space
         Allocate ( g_bins_k(LQ, norb,norb, Nbins), g_bins_r(LQ, norb, norb, Nbins)  )
         Allocate ( X_K(LQ,2) )

         ! Read-in the bins.
         Open ( Unit=10, File="ineq", status="unknown" ) 
         do nb = 1,Nbins
            do nk = 1,LQ
               read(10,*)  X_K(nk,1), X_K(nk,2)
               do no = 1,Norb
                  do no1 = 1,Norb
                     read(10,*) g_bins_k(nk,no,no1,nb)
                  enddo
               enddo
            enddo
         enddo
         close(10)
         
!  initialize lattice
	a1_p(1) =  1.0 ; a1_p(2) = 0.d0
	a2_p(1) =  0.5 ; a2_p(2) = sqrt(3.0)/2.d0
	L1_p    =  dble(L)*a1_p
	L2_p    =  dble(L)*a2_p
	PI      =  acos(-1.d0)
	CALL Make_Lattice(L1_p,L2_p,a1_p,a2_p,Latt)
	NDIM = 2*LQ
!  inverse FT to real space
	DO nb = 1,Nbins
		DO nk = 1,LQ
		xk_p =  dble(Latt%listk(nk,1))*Latt%b1_p + dble(Latt%listk(nk,2))*Latt%b2_p
		DO no = 1,norb
		DO no1 = 1,norb
		DO imj = 1,LQ
		ir_p =  dble(Latt%list(imj,1))*Latt%a1_p + dble(Latt%list(imj,2))*Latt%a2_p
		g_bins_r(imj,no,no1,nb) = g_bins_r(imj,no,no1,nb) + &
        &    exp( cmplx( 0.d0,-Iscalar(xk_p, ir_p)) ) * g_bins_k(nk,no,no1,nb)
		END DO
		END DO
		END DO
		END DO
	END DO
	g_bins_r = g_bins_r/cmplx(LQ,0.d0)


         Open (Unit=33,File="equalJ" ,status="unknown")
         Open (Unit=34,File="equalJ0",status="unknown")
         Zero = 1.E-10
         Do nk = 1,LQ
            write(33,"(F14.7,2x,F14.7)") X_K(nk,1), X_K(nk,2)
!!          ir_p = dble(Latt%list(nk,1))*Latt%a1_p + dble(Latt%list(nk,2))*Latt%a2_p
!!                      WRITE(33,'(F15.8,2X,F15.8,2X,I0)') ir_p(1), ir_p(2)
!!            X = abs(X_K(nk,1)) + abs(X_K(nk,2) )
            do no = 1,Norb
               do no1 = 1,Norb
                  call ERRCALCJ(g_bins_k(nk,no,no1,:), ZM, ZERR ) 
!!                  call ERRCALCJ(g_bins_r(nk,no,no1,:), ZM, ZERR ) 
                  Write(33,"(I3,2x,I3,F14.7,2x,F14.7,F14.7,2x,F14.7)") &
                       & no,no1, dble(ZM),dble(Zerr), aimag(ZM), Aimag(Zerr)

!!               call ERRCALCJ(g_bins_k(nk,1,1,:), ZM, ZERR )
!!               Write(33,"(F14.7,2x,F14.7,1x,F14.7,2x,F14.7,F14.7,2x,F14.7)") &
!!                       & X_K(nk,1), X_K(nk,2), dble(ZM),dble(Zerr), aimag(ZM), Aimag(Zerr)

!!                  if (X < Zero )  &
!!                       &  Write(34,"(I3,2x,I3,F14.7,2x,F14.7,F14.7,2x,F14.7)") &
!!                       &         no,no1, dble(ZM),dble(Zerr), aimag(ZM), Aimag(Zerr)
            
                   enddo
            enddo
            
         enddo
         close(33)
         close(34)
         
	DEALLOCATE(g_bins_k, g_bins_r, X_K)
       end Program Cov_eq
