function yesno(flag)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Assign y or n to logical value
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  character(len=1) :: yesno     ! y or n function
  logical          :: flag      ! distribution
!
! *********************** y or n assignment ****************************
!
! yesno: y or n function
!
  if (flag) then
    yesno = 'y'
  else
    yesno = 'n'
  endif
  return
end function yesno
! Copyright A.J. Koning 2021
