function truefalse(flag)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Assign T or F to integer value
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  logical :: truefalse ! true or false function
  integer :: flag      ! input integer
!
! *********************** y or n assignment ****************************
!
! yesno: y or n function
!
  if (flag <= 0) then
    truefalse = .false.
  else
    truefalse = .true.
  endif
  return
end function truefalse
! Copyright A.J. Koning 2021
