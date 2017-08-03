      Program enerJ
        Use ERRORS
	IMPLICIT REAL (KIND=8) (A-G,O-Z)
        
        REAL    (KIND=8), DIMENSION(:,:), ALLOCATABLE :: OBS
        REAL    (KIND=8), DIMENSION(:), ALLOCATABLE   :: EN, SIGN
        REAL    (KIND=8) :: XM, XERR

        
        ! Count the number of bins
        Open (Unit=10, File="ener1", status="unknown") 
        !Open (Unit=12, File="ener_hist", status="unknown") 
        nbins = 0
        do
            read(10,*,End=10)  X1, X2, X3, X4, X5, X6, X7, X8
            nbins = nbins + 1
            !write(6,*) nbins
            !write(12,"(I4,2x,F14.8,2x,F14.7,2x,F14.7)") nbins, X6, X2, X3
         enddo
10       continue
         Write(6,*) "# of bins: ", Nbins
         Close(10) 
         !Close(12)
         
         NP = NBINS
         NOBS = 8 
         
	ALLOCATE(OBS(NP,NOBS))
        ISEED = 99244        
!	Error on energy

        Open (Unit=25, File="statdat1", status="unknown") 
        read(25,*) NST, NS1, NS2, NSTEP
        Close(25)
	OPEN (UNIT=20, FILE='ener1', STATUS='old')
	NC = 0
	DO N = 1,NP
	   IF (N.GE.NST) THEN
	      NC = NC + 1
	      READ(20,*) OBS(NC,1), OBS(NC,2), OBS(NC,3), OBS(NC,4), OBS(NC,5), OBS(NC,6), OBS(NC,7), OBS(NC,8)
	   ELSE
	      READ(20,*) X1, X2, X3, X4, X5, X6, X7, X8
	   ENDIF
	ENDDO
        CLOSE(20)
2100	FORMAT(I6,2X,F16.8)

	OPEN (UNIT=21, FILE='enerJ', STATUS='unknown')
	WRITE(21,*) 'Effective number of bins, and bins: ', NC, NP
	NP_EFF = NC
        ALLOCATE (EN(NP_EFF), SIGN(NP_EFF))
	DO IOBS = 1,NOBS
           WRITE(21,*)
           DO I = 1,NP_EFF
              EN  (I) = OBS(I,IOBS)
              SIGN(I) = OBS(I,7)
           ENDDO
           IF (IOBS.EQ.1) WRITE(21,*) ' energy        '
           IF (IOBS.EQ.2) WRITE(21,*) ' rho           '
           IF (IOBS.EQ.3) WRITE(21,*) ' kint          '
           IF (IOBS.EQ.4) WRITE(21,*) ' ecoup         '
           IF (IOBS.EQ.5) WRITE(21,*) ' eJs           '
           IF (IOBS.EQ.6) WRITE(21,*) ' ehx           '
           IF (IOBS.EQ.7) WRITE(21,*) ' phase         '
		   IF (IOBS.EQ.8) WRITE(21,*) ' nflip         '
           DO NBIN = NS1, NS2, NSTEP
              if (NBIN.gt.0) then
                 IF (IOBS.EQ.7 ) then 
                    CALL ERRCALCJ(EN,XM,XERR,NBIN)
                 else
                    CALL ERRCALCJ(EN,SIGN,XM,XERR,NBIN)
                 endif
                 WRITE(21,2001) IOBS, XM,  XERR
                 ! Test
                 ! NBOOT = 40
                 ! CALL BOOTSTRAP( EN,XM_BS,XERR_BS,NBOOT,ISEED)
                 ! WRITE(21,2001) IOBS, XM_BS,  XERR_BS
              endif
           ENDDO
	ENDDO

	CLOSE(21)
2001    FORMAT('OBS : ', I4,4x,F12.6,2X, F12.6)

        DEALLOCATE (EN,SIGN,OBS)

	STOP
      END Program enerJ
