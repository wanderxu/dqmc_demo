integer function npbc(nr,l)
  implicit none
  integer, intent(in) :: nr
  integer, intent(in) :: l
  npbc = nr
  if (nr.gt.l) npbc = nr - l
  if (nr.lt.1) npbc = nr + l
end function npbc
