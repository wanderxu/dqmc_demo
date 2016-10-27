program cal_ssq
  implicit none

  integer, parameter :: dp = 8
  real(dp), parameter :: pi = dacos(-1.d0)
  integer :: l, lq, ltrot
  integer :: a1_p(2), a2_p(2)
  real(dp) :: b1_p(2), b2_p(2)
  integer, allocatable, dimension(:,:) :: listk, list

  ! local
  integer :: i, it, iq, j, eof, nk, nbin, ncount, nx, ny
  real(dp), allocatable, dimension(:,:) :: sstaur
  complex(dp), allocatable, dimension(:,:) :: zexpiqr
  complex(dp), allocatable, dimension(:) :: ssqtmp, ssq, ssq2
  character(len=80) :: arg, arg_mat(0:2)
  real(dp) :: qvec(2), ri(2), rtmp
  
  ! read the command argument
  i = 0
  do
      call get_command_argument(i,arg)
      if( len_trim(arg) == 0 ) exit
      arg_mat(i) = trim(arg)
      i = i + 1
  end do
  read(arg_mat(1),*) ltrot
  read(arg_mat(2),*) lq

  l = nint( sqrt(dble(lq)) )
  nk = (l+1)*(l+1)
  a1_p(1) = 1 ; a1_p(2) =  0
  a2_p(1) = 0 ; a2_p(2) =  1
  b1_p(1) = 2.d0*pi/dble(l) ; b1_p(2) = 0.d0
  b2_p(1) = 0.d0            ; b2_p(2) = 2.d0*pi/dble(l)

  allocate( sstaur(lq,ltrot) )
  allocate( ssqtmp(nk), ssq(nk), ssq2(nk) )
  allocate( zexpiqr(lq,nk) )
  allocate( listk(nk,2), list(lq,2) )

  ncount = 0
  do nx = 1,l
  do ny = 1,l
     ncount = ncount + 1
     list(ncount,1) = nx
     list(ncount,2) = ny
  enddo
  enddo

  nk = 0
  do j = 0, l
      do i = 0, l
          nk = nk+1
          listk(nk,1) = i-l/2
          listk(nk,2) = j-l/2
      end do
  end do

  do iq = 1, nk
      qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
      do i = 1, lq
          ri = dble(list(i,1))*a1_p + dble(list(i,2))*a2_p
          zexpiqr(i,iq) = exp( dcmplx( 0.d0,  qvec(1)*ri(1) + qvec(2)*ri(2) ) )
      end do
  end do

  ssq(:) = 0.d0
  ssq2(:) = 0.d0

  open( unit=555, file='sstaur_corrlt.bin', status='unknown' )

      eof = 0
      nbin = 0
      do 
          ssqtmp(:) = 0.d0
          do it = 1, ltrot
              do i = 1, lq
                  read(555,*,IOSTAT=eof) rtmp
                  if(eof.lt.0) go to 1001
                  do iq = 1, nk
                      ssqtmp(iq) = ssqtmp(iq) + rtmp*zexpiqr(i,iq)
                  end do
              end do
          end do
          nbin = nbin + 1
          ssqtmp(:) = ssqtmp(:) / dcmplx( dble(lq*ltrot), 0.d0 )
          ssq(:) = ssq(:) + ssqtmp(:)
          do iq = 1, nk
              ssq2(iq) = ssq2(iq) + dcmplx( real(ssqtmp(iq))*real(ssqtmp(iq)), aimag(ssqtmp(iq))*aimag(ssqtmp(iq)) )
          end do
      end do
1001 continue 
  
  ssq(:) = ssq(:) / nbin
  ssq2(:) = ssq2(:) / nbin
  do iq = 1, nk
      ssq2(iq) = ssq2(iq) - dcmplx( real(ssq(iq))*real(ssq(iq)), aimag(ssq(iq))*aimag(ssq(iq)) )
      qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
      write(*,'(2f16.8, 4e16.8)') qvec, ssq(iq), ssq2(iq)
      if(mod(iq,l+1)==0) write(*,*)
  end do

  close(555)
  deallocate( list, listk )
  deallocate( zexpiqr )
  deallocate( ssq2, ssq, ssqtmp )
  deallocate( sstaur )
end program cal_ssq
