function match(Eex)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Matching function
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl       ! single precision kind
!   dbl       ! double precision kind
! Variables for level density
!   E0save    ! E0 value saved for matching routine
!   EL        ! lower matching level energy
!   EP        ! higher matching level energy
!   logrho    ! logarithm of level density
!   NLo       ! lowest discrete level for temperature matching
!   NP        ! highest discrete level for temperature matching
!   temprho   ! temperature
!
! *** Declaration of local data
!
  implicit none
  integer   :: i       ! counter
  real(sgl) :: dEx     ! excitation energy bin for population arrays
  real(sgl) :: Eex     ! excitation energy
  real(sgl) :: logrhof ! log of value for total level density
  real(sgl) :: match   ! matching function
  real(sgl) :: temp    ! nuclear temperature
  real(dbl) :: factor1 ! help variable
  real(dbl) :: factor2 ! help variable
  real(dbl) :: rhof    ! value for total level density
  real(dbl) :: term    ! help variable
!
! *********************** Matching function ****************************
!
! match       : matching function
! pol1        : subroutine for interpolation of first order
!
  match = 0.
  dEx = 0.1
  i = max(int(Eex / dEx), 1)
  call pol1(i * dEx, (i + 1) * dEx, temprho(i), temprho(i + 1), Eex, temp)
  if (temp > 0.) then
    call pol1(i * dEx, (i + 1) * dEx, logrho(i), logrho(i + 1), Eex, logrhof)
    rhof = exp(dble(logrhof))
    if (E0save == 1.e-20) then
      factor1 = exp( - Eex / temp)
      factor2 = exp(dble(EP / temp))
      if (EL /= 0.) factor2 = factor2 - exp(EL / temp)
      term = real(min(temp * rhof * factor1 * factor2, 1.d30))
      match = term + NLo - NP
    else
      match = Eex - temp * log(temp * rhof) - E0save
    endif
  endif
  return
end function match
! Copyright A.J. Koning 2021
