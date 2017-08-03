!=======================================================================
   PROGRAM rebin
!=======================================================================
!         
!   USE precdef
   IMPLICIT NONE
   
   Integer, Parameter :: &
        long = selected_int_kind(9),            &
        single = kind(1.0), 			&
        double = kind(1.0D0)
   REAL(kind=double), DIMENSION(:,:), ALLOCATABLE  :: kvec, obs_ener1
   REAL(kind=double) :: r1,r2,mn,c1,c2,c3,c4,c5,c6
   COMPLEX(kind=double), DIMENSION(:,:,:,:), ALLOCATABLE :: obs
   COMPLEX(kind=double) :: z1
   INTEGER :: I,J,no,no1,nk,nb,norb,nbins,binskip,LQ,L,ifile, &
   	&	nfile,nbinSize,binSize,bs
   CHARACTER(len=64) :: doSkip,fname
   CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: file

!  Read in lattice size, here I assume the lattice size is the first parameter in the second line of paramC_sets
   OPEN(UNIT=5,FILE='paramC_sets',STATUS='OLD') 
   READ(5,*) i
   READ(5,*) L
   CLOSE(5)
   LQ   = L*L
   norb = 0

!  Read in the filenames of the files to rebin
!  Note always put ener1 as the first file to read: ifile = 1.
   OPEN(UNIT=5,FILE='rebin.in',STATUS='OLD') 
   READ(5,*) nfile
   ALLOCATE(file(1:nfile))
   DO i = 1,nfile
		READ(5,*) file(i)
   END DO
   CLOSE(5)

   DO ifile = 1,nfile   
!     Determine the number of bins and orbitals
		WRITE(*,'(A)') trim(file(ifile))
		WRITE(*,'(3X,A)') 'Analyzing ...'

      OPEN(UNIT=5,FILE=trim(file(ifile)),STATUS='OLD') 
      READ(5,*,iostat=i) doSkip,binskip
      IF ((i /= 0).OR.(trim(doSkip) /= 'skip')) THEN
         REWIND(5)
         binskip = 0
      END IF
      
   IF (ifile == 1) THEN
      WRITE(*,'(3X,A)') 'Analyzing ener1'
      nbins = 0
      DO
		READ(5,*,iostat=i) c1,c2,c3
		READ(5,*,iostat=i) c4,c5,c6
!!		IF (i /= 0) EXIT
!!	  END DO
      IF(i /= 0) EXIT
      nbins = nbins + 1
      END DO   
      
   ELSE 
      DO J = 1,binskip
         DO nk = 1,LQ
            IF ((nk == 1).AND.(binskip == 1)) THEN
               READ(5,*,iostat=i) r1,r2,norb
               IF (i < 0) THEN
                  EXIT
               ELSEIF ((i > 0).OR.(norb == 0))THEN
                  norb = 2
                  BACKSPACE(5)
               END IF
            ELSE
               READ(5,*,iostat=i) r1,r2
               IF (i < 0) EXIT
            END IF
            DO no = 1,norb
               DO no1 = 1,norb
                  READ(5,*,iostat=i) z1
                  IF (i /= 0) EXIT
               END DO
            END DO
         END DO
         IF (i /= 0) EXIT
      END DO
      nbins = 0
      DO
         DO nk = 1,LQ
            IF ((nk == 1).AND.(nbins == 0)) THEN
               READ(5,*,iostat=i) r1,r2,norb
               IF (i < 0) THEN
                  EXIT
               ELSEIF ((i > 0).OR.(norb == 0))THEN
                  norb = 2
                  BACKSPACE(5)
               END IF
            ELSE
               READ(5,*,iostat=i) r1,r2
               IF (i < 0) EXIT
            END IF
            DO no = 1,norb
               DO no1 = 1,norb
                  READ(5,*,iostat=i) z1
                  IF (i /= 0) EXIT
               END DO
            END DO
         END DO
         IF (i /= 0) EXIT
         nbins = nbins+1
      END DO
    END IF
		WRITE(*,'(3X,A,I0,A)') 'There are ',nbins,' valid bin(s)'
      ALLOCATE(obs_ener1(6,Nbins))
      obs_ener1 = 0.0D0
      ALLOCATE(obs(LQ,norb,norb,Nbins),kvec(LQ,2))
      obs = 0.0D0
      
!     Read in data
   IF (ifile == 1) THEN
     REWIND(5)
     DO nb = 1, Nbins
	  READ(5,*) obs_ener1(1,nb), obs_ener1(2,nb), obs_ener1(3,nb)
	  READ(5,*) obs_ener1(4,nb), obs_ener1(5,nb), obs_ener1(6,nb)
     END DO
   ELSE
     REWIND(5)
      DO J = 1,binskip
         IF (J == 1) READ(5,*,iostat=i) doSkip,binskip
         DO nk = 1,LQ
            READ(5,*,iostat=i) r1,r2
            IF (i /= 0) EXIT
            DO no = 1,norb
               DO no1 = 1,norb
                  READ(5,*,iostat=i) z1
                  IF (i /= 0) EXIT
               END DO
            END DO
         END DO
         IF (i /= 0) EXIT
      END DO
      r1 = 0.0D0
      mn = 0.0D0
      DO nb = 1,Nbins
         DO nk = 1,LQ
            READ(5,*) kvec(nk,1),kvec(nk,2)
            DO no = 1,norb
               DO no1 = 1,norb
                  READ(5,*) obs(nk,no,no1,nb)
               END DO
            END DO
         END DO
      END DO
   END IF
   CLOSE(5)
		
!     Determine the number of possible rebinning sizes, rebin and write to disk
		nbinSize = int(log(dble(nbins))/log(2.0d0))
		IF (nbinSize == 0) WRITE(*,'(3X,A)') 'There are not enough bins to rebin!'
		
		DO bs = 1,nbinSize
			binSize = 2**bs
			WRITE(fname,'(2A,I0)') trim(file(ifile)),'_',binSize
			WRITE(*,'(3X,A,A18,A,I4,A,I4,A)') 'Writing ',trim(fname),' (re-binsize ',binsize,', ',int(nbins/binsize),' bin(s))'

		   OPEN(UNIT=5,FILE=trim(fname),STATUS='UNKNOWN')
	IF (ifile == 1) THEN
	 DO nb = 1,int(nbins/binsize)
	    WRITE(5,*) sum(obs_ener1(1,nb:nb+binSize-1)/dble(binSize)), sum(obs_ener1(2,nb:nb+binSize-1)/dble(binSize)), sum(obs_ener1(3,nb:nb+binSize-1)/dble(binSize))
	    WRITE(5,*) sum(obs_ener1(4,nb:nb+binSize-1)/dble(binSize)), sum(obs_ener1(5,nb:nb+binSize-1)/dble(binSize)), sum(obs_ener1(6,nb:nb+binSize-1)/dble(binSize))
	 END DO
	ELSE
         DO nb = 1,int(nbins/binsize)
	         DO nk = 1,LQ
   	         WRITE(5,*) kvec(nk,1),kvec(nk,2),norb
	   	   	DO no = 1,norb
   	   	   	DO no1 = 1,norb
            	      WRITE(5,*) sum(obs(nk,no,no1,nb:nb+binSize-1)/dble(binSize))
	   	      	END DO
	   	   	END DO
		END DO
		
	END DO
	END IF
   		CLOSE(5)				
		END DO
		
   DEALLOCATE(obs_ener1,obs,kvec)
   END DO

   END PROGRAM rebin
