subroutine preq
#IFDEF _OPENMP
  USE OMP_LIB
#ENDIF
  use blockc
  use obser
  implicit none

  include 'mpif.h'

  integer :: i, j, imj
  complex(dp) :: expiqr, spipi, spipi0
  real(dp) :: sq_ising_00

  complex(dp), allocatable, dimension(:) :: mpi_c1
  complex(dp), allocatable, dimension(:,:) :: mpi_c2
  real(dp), allocatable, dimension(:) :: mpi_r1
  real(dp), allocatable, dimension(:,:) :: mpi_r2
  integer, allocatable, dimension(:) :: mpi_i1
  integer, allocatable, dimension(:,:) :: mpi_i2
  complex(dp) :: mpi_obs_bin(10)

  real(dp) :: qvec(2), rij(2)
  integer :: iq, n, itau
  complex(dp) :: sq_ising_qwn_tmp, zexpiwtqr
  real(dp), allocatable, dimension(:) :: isingtau0r_crr


  !allocate(mpi_c1(lq))
  !allocate(mpi_c2(lq,lq))
  allocate(mpi_r1(lq))
  if(lsstau) allocate(mpi_r2(lq,ltrot))
  allocate(mpi_i1(lq) )
  if(lsstau) allocate(mpi_i2(lq,ltrot) )

  !!!!call mpi_reduce( spinz, mpi_c1, lq, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  !!!!spinz(:) = mpi_c1(:)

  !!!!call mpi_reduce(spin_corrlt,mpi_c2,lq*lq,mpi_complex16,mpi_sum, 0,mpi_comm_world,ierr)
  !!!!spin_corrlt(:,:) = mpi_c2(:,:)

  call mpi_reduce(isingzz_corrlt,mpi_i1,lq,mpi_integer,mpi_sum, 0,mpi_comm_world,ierr)
  isingzz_corrlt(:) = mpi_i1(:)

  if(lsstau) then
  call mpi_reduce( isingzztau_corrlt, mpi_i2, lq*ltrot, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr )
  isingzztau_corrlt(:,:) = mpi_i2(:,:)
  end if

  call mpi_reduce( obs_bin, mpi_obs_bin, 10, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
  obs_bin(:) = mpi_obs_bin(:)

  if( irank .eq. 0 ) then
      !!!!!! calculate  the S(pi,pi)
      !!!!!! S(pi,pi) = 1/L^2 \sum_ij <  1/2 (ni_up -ni_dn) *  1/2 ( nj_up - nj_dn ) >  exp( i * (pi,pi) * r(i-j) )
      !!!!!! normalize
      !!!!!spin_corrlt(:,:) = spin_corrlt(:,:) / dcmplx( isize * nobs )
      !!!!!spinz(:) = spinz(:) / dcmplx( isize * nobs )
      !!!!!spipi0 = czero
      !!!!!spipi = czero
      !!!!!do i = 1, lq
      !!!!!    do j = 1, lq
      !!!!!        expiqr = exp( dcmplx( 0.d0, ( pi*dble(list(j,1)-list(i,1)) + pi*dble( list(j,2)-list(i,2) ) ) ) )
      !!!!!        spipi  = spipi + spin_corrlt(j,i) * expiqr
      !!!!!        spipi0 = spipi0 + spinz(j)*spinz(i) * expiqr
      !!!!!    end do
      !!!!!end do

      !!!!!spipi  = spipi  / dcmplx( dble( 4 * lq * lq), 0.d0 )  ! 1/4 comes from 1/2 spin
      !!!!!spipi0 = spipi0 / dcmplx( dble( 4 * lq * lq), 0.d0 )  ! 1/4 comes from 1/2 spin

      !!!!!open (unit=40,file='spin_corrlt.bin',status='unknown', action="write", position="append")

      !!!!!write(40, '(3(2e16.8,2x))') spipi, spipi0, spipi-spipi0

      !!!!!close(40)


      ! calcuate the Sising(0,0)
      mpi_r1(:) = dble(isingzz_corrlt(:)) / dble( isize * nobs )
      sq_ising_00 = 0.d0
      do imj = 1, lq
          sq_ising_00 = sq_ising_00 + mpi_r1(imj)
      end do

      sq_ising_00 = sq_ising_00 / dble( lq * lq )
      !!!sq_ising_00 = sq_ising_00 / dble( lq * lq * ltrot )

      open (unit=80,file='isingzz_corrlt.bin',status='unknown', action="write", position="append")

      write(80, '(e16.8)') sq_ising_00

      close(80)

      ! calculate Sising(q,iwn)
      if( lsstau ) then
          mpi_r2(:,:) = dble(isingzztau_corrlt(:,:)) / dble( isize * nobs )
          open (unit=85,file='isingzztau_corrlt.bin',status='unknown', action="write", position="append")
          do iq = 1, lq
              qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
              do n = 0, nuse
                  sq_ising_qwn_tmp = czero
!$OMP PARALLEL &
!$OMP PRIVATE ( itau, imj, zexpiwtqr )
!$OMP DO REDUCTION ( + : sq_ising_qwn_tmp )
                  do itau = 1, ltrot
                      do imj = 1, lq
                          zexpiwtqr = zexpiwt(itau-1,n) / zexpiqr(imj,iq)
                          sq_ising_qwn_tmp = sq_ising_qwn_tmp + dcmplx( mpi_r2(imj,itau), 0.d0 ) * zexpiwtqr
                      end do
                  end do
!$OMP END DO
!$OMP END PARALLEL
                  !sq_ising_qwn(n,iq) = sq_ising_qwn_tmp
                  write(85, '(2f16.8,f16.8,2e16.8)') qvec(:), wn(n), sq_ising_qwn_tmp*dcmplx(dtau/dble(lq),0.d0)
              end do
          end do
          close(85)

          open (unit=89,file='sstaur_corrlt.bin',status='unknown', action="write", position="append")
          do itau = 1, ltrot
              do imj = 1, lq
                  write(89,'(e22.12)') mpi_r2(imj,itau)
              end do
          end do
          close(89)

          if( lsstau0r ) then
              open (unit=195,file='isingzztau0r_corrlt.bin',status='unknown', action="write", position="append")
              allocate( isingtau0r_crr(num_equ_distance) )
              isingtau0r_crr(:) = zero
              do i = 1, lq
                  imj = distance_index(i) ! imj no.
                  j = equ_distance(i)
                  isingtau0r_crr(j) = isingtau0r_crr(j) + mpi_r2(imj,1)/ dble( imjdeg(imj) )
              end do
              do i = 1, num_equ_distance
                  write( 195, '(2f16.8)') irre_distance_len(i), isingtau0r_crr(i)/dble( irre_distance_deg(i) )
              end do
              deallocate( isingtau0r_crr )

              !!!do i = 1, lq
              !!!    imj = distance_index(i) ! imj no.
              !!!    write( 195, '(f16.8)') mpi_r2(imj,1)
              !!!end do

              close(195)

          end if
      end if

      ! calculate obs_bin
      obs_bin(:) = obs_bin(:) / dcmplx( isize * nobs )
      obs_bin(8) = obs_bin(9) / (obs_bin(8)**2)
      open (unit=90,file='ener1.bin',status='unknown', action="write", position="append")
      write(90, '(8(2e16.8,4x))') obs_bin(1)/dcmplx(dble(ndim),0.d0), obs_bin(2), obs_bin(8), obs_bin(3)/dcmplx(dble(ndim),0.d0), obs_bin(4), obs_bin(5), obs_bin(6), obs_bin(7)
      close(90)

  end if

  call mpi_barrier( mpi_comm_world, ierr )

  if(lsstau) deallocate(mpi_i2)
  deallocate(mpi_i1)
  if(lsstau) deallocate(mpi_r2)
  deallocate(mpi_r1)
  !deallocate(mpi_c2)
  !deallocate(mpi_c1)

end subroutine preq
