function quasideuteron(Egamma)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Quasi-deuteron model of Chadwick and Oblozinsky
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
!   sgl        ! single precision kind
! Variables for main input
!   Atarget    ! mass number of target nucleus
!   Ntarget    ! neutron number of target nucleus
!   Ztarget    ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: Egamma        ! gamma energy
  real(sgl) :: fpauli        ! Pauli-blocking function of Chadwick
  real(sgl) :: freedeut      ! free deuteron cross section
  real(sgl) :: levinger      ! Levinger parameter
  real(sgl) :: quasideuteron ! Quasi-deuteron function of Chadwick and Oblozinsky
!
! ****** Calculate quasi-deuteron photo-absorption cross section *******
!
  levinger = 6.5
  if (Egamma > 2.224) then
    freedeut = 61.2 / (Egamma **3) * (Egamma - 2.224) **1.5
  else
    freedeut = 0.
  endif
  if (Egamma <= 140.) then
    if (Egamma >= 20.) then
      fpauli = 8.3714e-2 - 9.8343e-3 * Egamma + 4.1222e-4 * Egamma * Egamma - &
 &      3.4762e-6 * (Egamma **3) + 9.3537e-9 * (Egamma **4)
    else
      fpauli = exp( - 73.3 / Egamma)
    endif
  else
    fpauli = exp( - 24.2348 / Egamma)
  endif
  quasideuteron = levinger * Ntarget * Ztarget / real(Atarget) * freedeut * fpauli
  return
end function quasideuteron
! Copyright A.J. Koning 2021
