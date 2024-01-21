subroutine sumrules
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Giant resonance sum rules
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
!   sgl         ! single precision kind
! All global variables
!   numlev2     ! maximum number of levels
! Variables for main input
!   Atarget     ! mass number of target nucleus
!   k0          ! index of incident particle
! Variables for nuclides
!   Nindex      ! neutron number index for residual nucleus
!   Zindex      ! charge number index for residual nucleus
! Constants
!   onethird    ! 1 / 3
!   twothird    ! 2 / 3
! Variables for deformation parameters
!   betagr      ! deformation parameter for giant resonance
!   deform      ! deformation parameter
!   deftype     ! deformation length (D) or parameter (B)
!   Egrcoll     ! energy of giant resonance
!   Ggrcoll     ! width of giant resonance
! Variables for levels
!   edis        ! energy of level
!   jdis        ! spin of level
!
! *** Declaration of local data
!
  implicit none
  integer   :: i       ! counter
  integer   :: Ncomp   ! neutron number index for compound nucleus
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: Zcomp   ! proton number index for compound nucleus
  integer   :: Zix     ! charge number index for residual nucleus
  real(sgl) :: betasq  ! square of deformation parameter
  real(sgl) :: deform1 ! deformation parameter
  real(sgl) :: edis1   ! energy of level
  real(sgl) :: jdis1   ! spin of level
  real(sgl) :: Sgmr    ! sum rule for giant monopole resonance
  real(sgl) :: Sgor    ! sum rule for giant octupole resonance
  real(sgl) :: Sgqr    ! sum rule for giant quadrupole resonance
  real(sgl) :: Sheor   ! sum rule for high energy octupole resonance
  real(sgl) :: Sleor   ! sum rule for low energy octupole resonance
!
! ******************* Strength of giant resonances *********************
!
! Using sum rules for giant resonances, the collective strength in the continuum is determined.
!
  Sgmr = 23.*Atarget**(-5./3.)
  Egrcoll(0, 1) = 18.7 - 0.025 * Atarget
  Ggrcoll(0, 1) = 3.
  Sgqr = 575. * Atarget **( - 5. / 3.)
  Egrcoll(2, 1) = 65. * Atarget **( - onethird)
  Ggrcoll(2, 1) = 85. * Atarget **( - twothird)
  Sgor = 1208. * Atarget **( - 5. / 3.)
  Sleor = 0.3 * Sgor
  Egrcoll(3, 1) = 31. * Atarget **( - onethird)
  Ggrcoll(3, 1) = 5.
  Sheor = 0.7 * Sgor
  Egrcoll(3, 2) = 115. * Atarget **( - onethird)
  Ggrcoll(3, 2) = 9.3 - Atarget / 48.
!
! Subtract the deformation parameters of the collective low lying states from the sum.
!
  Zcomp = 0
  Ncomp = 0
  Zix = Zindex(Zcomp, Ncomp, k0)
  Nix = Nindex(Zcomp, Ncomp, k0)
  do i = 1, numlev2
    deform1=deform(Zix,Nix,i)
    if (deform1 /= 0.) then
      if (deftype(Zix, Nix) == 'D') deform1=deform1/(1.24*Atarget**onethird)
      betasq = deform1 * deform1
      edis1 = edis(Zix, Nix, i)
      jdis1 = int(jdis(Zix, Nix, i))
      if (jdis1 == 0) Sgmr = Sgmr - betasq * edis1
      if (jdis1 == 2) Sgqr = Sgqr - betasq * edis1
      if (jdis1 == 3) Sleor = Sleor - betasq * edis1
    endif
  enddo
!
! Determine final GR deformation parameters.
!
  if (Sgmr > 0.) betagr(0, 1) = sqrt(Sgmr / Egrcoll(0, 1))
  if (Sgqr > 0.) betagr(2, 1) = sqrt(Sgqr / Egrcoll(2, 1))
  if (Sleor > 0.) betagr(3, 1) = sqrt(Sleor / Egrcoll(3, 1))
  betagr(3, 2) = sqrt(Sheor / Egrcoll(3, 2))
  return
end subroutine sumrules
! Copyright A.J. Koning 2021
