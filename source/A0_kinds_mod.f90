module A0_kinds_mod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Definition of single and double precision variables
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
end module A0_kinds_mod
! Copyright A.J. Koning 2021
