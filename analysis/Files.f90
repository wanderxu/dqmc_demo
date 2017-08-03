       Subroutine Files(I1, File_k  )

          Implicit Real  (KIND=8) (A-G,O-Z)
          Implicit Integer (H-N)

          Integer :: I1, Ntype  ! = 1 Files for Green function.
                                ! = 2 Files for self-energy.
          Character (16) :: File_k

          Character (1), Dimension(:), Allocatable :: C1
          Character (8)  :: Res 


          I = I1
          nc = 0
          i_tens = 10
          do 
             nc = nc + 1
             if (I .ge. i_tens) then
                i_tens = i_tens * 10 
             else
                exit
             endif
          enddo

          !write(6,*) 'nc is: ', nc 
          !write(6,*) 'I is: ', I
          Allocate(C1(nc))
          
          do n = nc,1, -1
             J = mod(I,10)
             Open (Unit=10, File='tmp',status='unknown')
             write(10,2001) J
             close(10)
             Open(Unit=10,file='tmp',status='unknown')
             read(10,*) C1(n)
             close(10)
             I = I/10
          enddo
2001      format(I1)
          Open (Unit=10, File='tmp',status='unknown')
          write(10,*) C1
          close(10)
          Open (Unit=10, File='tmp',status='unknown')
          read(10,*) Res
          close(10)

          File_k   = 'g_k'//Res 

        end Subroutine Files
