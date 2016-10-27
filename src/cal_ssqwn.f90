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
  complex(dp), allocatable, dimension(:,:) :: zexpiqr, zexpiwt
  complex(dp), allocatable, dimension(:) :: ssqtmp, ssq, ssq2
  real(dp), allocatable, dimension(:) :: wn
  character(len=80) :: arg, arg_mat(0:2)
  real(dp) :: qvec(2), ri(2), rtmp, dtau
  integer :: nuse
  integer :: n, imj, itau
  complex(dp) :: sq_ising_qwn_tmp, zexpiwtqr 
  
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
  read(arg_mat(3),*) nuse

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
  do j = 0, l-1
      do i = 0, l-1
          nk = nk+1
          listk(nk,1) = i-l/2
          listk(nk,2) = j-l/2
      end do
  end do

  allocate( wn(0:nuse) )
  dtau = 0.05d0
  do i = 0, nuse
      wn(i) = 2.d0*dble(i)*pi/(ltrot*0.05d0)
  end do

  allocate( zexpiwt(0:ltrot-1,0:nuse) )
  ! set zexpiwt, zexpiqr
  do n = 0, nuse
      do i = 0, ltrot-1
          zexpiwt(i,n) = exp( dcmplx( 0.d0, wn(n)*dble(i)*dtau ) )
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
  open( unit=85,file='isingzztau_corrlt_morewn.bin',status='unknown')

      eof = 0
      nbin = 0
      do 
          do itau = 1, ltrot
              do imj = 1, lq
                  read(555,*,IOSTAT=eof) sstaur(imj,itau)
                  if(eof.lt.0) go to 1001
              end do
          end do

          nbin = nbin + 1

          do iq = 1, lq
              qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
              do n = 0, nuse
                  sq_ising_qwn_tmp = dcmplx(0.d0,0.d0)
                  do itau = 1, ltrot
                      do imj = 1, lq
                          zexpiwtqr = zexpiwt(itau-1,n) / zexpiqr(imj,iq)
                          sq_ising_qwn_tmp = sq_ising_qwn_tmp + dcmplx( sstaur(imj,itau), 0.d0 ) * zexpiwtqr
                      end do
                  end do
                  !sq_ising_qwn(n,iq) = sq_ising_qwn_tmp
                  write(85, '(2f16.8,f16.8,2e16.8)') qvec(:), wn(n), sq_ising_qwn_tmp*dcmplx(dtau/dble(lq),0.d0)
              end do
          end do

      end do
1001 continue 

  close(85)
  close(555)
  deallocate( zexpiwt )
  deallocate( wn )
  deallocate( list, listk )
  deallocate( zexpiqr )
  deallocate( ssq2, ssq, ssqtmp )
  deallocate( sstaur )
end program cal_ssq
