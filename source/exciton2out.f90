subroutine exciton2out
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of two-component exciton model parameters
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
!   sgl           ! single precision kind
! Variables for preequilibrium
!   M2constant    ! constant for matrix element in exciton model
!   Rnunu         ! ratio for two - component matrix element
!   Rnupi         ! ratio for two - component matrix element
!   Rpinu         ! ratio for two - component matrix element
!   Rpipi         ! ratio for two - component matrix element
! Variables for main input
!   Ainit         ! mass number of initial compound nucleus
! Constants
!   hbar          ! Planck's constant / 2.pi in MeV.s
!   parname       ! name of particle
! Variables for preequilibrium initialization
!   maxpar        ! maximal particle number
! Variables for exciton model
!   M2nunu        ! square of neutron - neutron matrix element
!   M2nupi        ! square of neutron - proton matrix element
!   M2pinu        ! square of proton - neutron matrix element
!   M2pipi        ! square of proton - proton matrix element
!   wemispart2    ! two - component emission rate per par
!   wemistot2     ! total two - component emission rate p
! Variables for preequilibrium
!   Ecomp         ! total energy of composite system
!   p0            ! initial particle number
!   pnu0          ! initial neutron number
!   ppi0          ! initial proton number
!   Spre          ! time - integrated strength of two - component exciton state
!
! *** Declaration of local data
!
  implicit none
  integer   :: h            ! help variable
  integer   :: hnu          ! neutron hole number
  integer   :: hpi          ! proton hole number
  integer   :: n            ! exciton number
  integer   :: p            ! particle number
  integer   :: pnu          ! neutron particle number
  integer   :: ppi          ! proton particle number
  integer   :: type         ! particle type
  real(sgl) :: lambdanupi   ! neutron-proton transition rate for n --> n
  real(sgl) :: lambdanuplus ! neutron transition rate for n --> n+2
  real(sgl) :: lambdapinu   ! proton-neutron transition rate for n --> n
  real(sgl) :: lambdapiplus ! proton transition rate for n --> n+2
!
! ************************ Exciton model *******************************
!
  write(*, '(/" ++++++++++ TWO-COMPONENT EXCITON MODEL ++++++++++")')
!
! 1. Output of matrix element
!
! matrix    : subroutine for matrix element for exciton model
!
  write(*, '(/" 1. Matrix element for E= ", f8.3/)') Ecomp
  write(*, '(" Constant for matrix element : ", f7.3)') M2constant
  write(*, '(" p-p ratio for matrix element: ", f7.3)') Rpipi
  write(*, '(" n-n ratio for matrix element: ", f7.3)') Rnunu
  write(*, '(" p-n ratio for matrix element: ", f7.3)') Rpinu
  write(*, '(" n-p ratio for matrix element: ", f7.3/)') Rnupi
  write(*, '(" p(p) h(p) p(n) h(n)     M2pipi      M2nunu      M2pinu      M2nupi"/)')
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          h = hpi + hnu
          n = p + h
          call matrix(Ainit, n)
          write(*, '(1x, 4(i2, 3x), 4es12.5)') ppi, hpi, pnu, hnu, M2pipi, M2nunu, M2pinu, M2nupi
        endif
      enddo
    enddo
  enddo
!
! 2. Output of emission rates or escape widths
!
  write(*, '(/" 2. Emission rates or escape widths"/)')
  write(*, '(" A. Emission rates ( /sec)"/)')
  write(*, '(" p(p) h(p) p(n) h(n)", 4x, 7(a8, 4x), "Total"/)') (parname(type), type = 0, 6)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(*, '(1x, 4(i2, 3x), 8es12.5)') ppi, hpi, pnu, hnu, (wemispart2(type, ppi, hpi, pnu, hnu), type = 0, 6), &
 &        wemistot2(ppi, hpi, pnu, hnu)
        endif
      enddo
    enddo
  enddo
  write(*, '(/" B. Escape widths (MeV)"/)')
  write(*, '(" p(p) h(p) p(n) h(n)", 4x, 7(a8, 4x), "Total"/)') (parname(type), type = 0, 6)
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(*, '(1x, 4(i2, 3x), 8es12.5)') ppi, hpi, pnu, hnu, &
 &        (wemispart2(type, ppi, hpi, pnu, hnu) * hbar, type = 0, 6), wemistot2(ppi, hpi, pnu, hnu) * hbar
        endif
      enddo
    enddo
  enddo
!
! 3. Output of transition rates or damping widths and total widths
!
  write(*, '(/" 3. Internal transition rates or damping widths, total widths"/)')
  write(*, '(" A. Internal transition rates ( /sec)"/)')
  write(*, '(" p(p) h(p) p(n) h(n)     lambdapiplus   lambdanuplus    lambdapinu     lambdanupi"/)')
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(*, '(1x, 4(i2, 3x), 4es15.6)') ppi, hpi, pnu, hnu, lambdapiplus(0, 0, ppi, hpi, pnu, hnu), &
 &          lambdanuplus(0, 0, ppi, hpi, pnu, hnu), lambdapinu(0, 0, ppi, hpi, pnu, hnu), lambdanupi(0, 0, ppi, hpi, pnu, hnu)
        endif
      enddo
    enddo
  enddo
  write(*, '(/" B. Damping widths (MeV)"/)')
  write(*, '(" p(p) h(p) p(n) h(n)     gammapiplus    gammanuplus    gammapinu      gammanupi"/)')
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(*, '(1x, 4(i2, 3x), 4es15.6)') ppi, hpi, pnu, hnu, lambdapiplus(0, 0, ppi, hpi, pnu, hnu) * hbar, &
 &          lambdanuplus(0, 0, ppi, hpi, pnu, hnu) * hbar, lambdapinu(0, 0, ppi, hpi, pnu, hnu) * hbar, &
 &          lambdanupi(0, 0, ppi, hpi, pnu, hnu) * hbar
        endif
      enddo
    enddo
  enddo
  write(*, '(/" C. Total widths (MeV)"/)')
  write(*, '(" p(p) h(p) p(n) h(n)      gammatot"/)')
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(*, '(1x, 4(i2, 3x), es15.6)') ppi, hpi, pnu, hnu, hbar * (lambdapiplus(0, 0, ppi, hpi, pnu, hnu) + &
 &          lambdanuplus(0, 0, ppi, hpi, pnu, hnu) + lambdapinu(0, 0, ppi, hpi, pnu, hnu) + &
 &          lambdanupi(0, 0, ppi, hpi, pnu, hnu) + wemistot2(ppi, hpi, pnu, hnu))
        endif
      enddo
    enddo
  enddo
!
! 4. Output of lifetimes of exciton states
!
  write(*, '(/" 4. Lifetimes")')
  write(*, '(" p(p) h(p) p(n) h(n)      Strength"/)')
  do p = p0, maxpar
    do ppi = ppi0, maxpar
      hpi = ppi - ppi0
      do pnu = pnu0, maxpar
        if (ppi + pnu == p) then
          hnu = pnu - pnu0
          write(*, '(1x, 4(i2, 3x), es15.6)') ppi, hpi, pnu, hnu, Spre(ppi, hpi, pnu, hnu)
        endif
      enddo
    enddo
  enddo
  return
end subroutine exciton2out
! Copyright A.J. Koning 2021
