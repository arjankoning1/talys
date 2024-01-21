function phdens2(Zix, Nix, ppi, hpi, pnu, hnu, gsp, gsn, Eex, Ewell, surfwell)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Two-component particle-hole state density
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
!   phmodel       ! particle - hole state density model
! Variables for level density
!   edens         ! energy grid for tabulated level density
! Variables for particle-hole state densities
!   Ephdensmax    ! maximum energy on p - h state density
!   nenphdens     ! number of energies for p - h state density
!   phexist2      ! flag for existence of p - h state density
!   phtable2      ! p - h state density
! Variables for preequilibrium initialization
!   Apauli2       ! two - component Pauli blocking correction
!   nfac          ! n!
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell   ! flag for surface effects in finite well
  integer   :: h          ! help variable
  integer   :: hnu        ! neutron hole number
  integer   :: hpi        ! proton hole number
  integer   :: n1         ! number of coordinate grid points
  integer   :: nex2       ! counter
  integer   :: Nix        ! neutron number index for residual nucleus
  integer   :: nnu        ! neutron exciton number
  integer   :: npi        ! proton exciton number
  integer   :: p          ! particle number
  integer   :: pnu        ! neutron particle number
  integer   :: ppi        ! proton particle number
  integer   :: Zix        ! charge number index for residual nucleus
  real(sgl) :: Ap         ! two-component Pauli blocking correction factor
  real(sgl) :: eb         ! help variable
  real(sgl) :: ee         ! energy
  real(sgl) :: Eex        ! excitation energy
  real(sgl) :: Ewell      ! depth of potential well
  real(sgl) :: fac1       ! help variable
  real(sgl) :: factor     ! multiplication factor
  real(sgl) :: factorn    ! help variable
  real(sgl) :: factorp    ! help variable
  real(sgl) :: finitewell ! correction function for finite well depth
  real(sgl) :: gsn        ! single-particle neutron level density parameter
  real(sgl) :: gsp        ! single-particle proton level density parameter
  real(sgl) :: ldb        ! level density
  real(sgl) :: lde        ! level density
  real(sgl) :: ldtab      ! tabulated level density
  real(sgl) :: lldb       ! log of level density
  real(sgl) :: llde       ! log of level density
  real(sgl) :: phdens2    ! function for two-component particle-hole state density
!
! 1. State density of Kalbach, Phys. Rev. C33, 818 (1986).
!
! finitewell: function for correction for finite well depth
!
! The finite depth of the hole is included.
! If the uncorrected state density is required, it should be specified by Ewell=0.
! (which thus actually means Ewell=infinity).
!
  phdens2 = 0.
  if (ppi < 0 .or. hpi < 0 .or. pnu < 0 .or. hnu < 0) return
  if (ppi + hpi + pnu + hnu == 0) return
  if (phmodel == 1 .or. .not. phexist2(Zix, Nix, ppi, hpi, pnu, hnu)) then
    Ap = Apauli2(ppi, hpi, pnu, hnu)
    factorn = (pnu * pnu + hnu * hnu + pnu + hnu) / (4. * gsn)
    factorp = (ppi * ppi + hpi * hpi + ppi + hpi) / (4. * gsp)
    if (Ap + factorn + factorp >= Eex) return
    h = hpi + hnu
    p = ppi + pnu
    npi = ppi + hpi
    nnu = pnu + hnu
    n1 = npi + nnu - 1
    fac1 = nfac(ppi) * nfac(hpi) * nfac(pnu) * nfac(hnu) * nfac(n1)
    factor = gsp **npi * gsn **nnu / fac1
    phdens2 = factor * (Eex - Ap) **n1
    phdens2 = phdens2 * finitewell(p, h, Eex, Ewell, surfwell)
  else
!
! 2. Tabulated particle-hole state densities
!
! locate       : subroutine to find value in ordered table
!
    if (Eex <= 0.) return
    if (Eex <= Ephdensmax) then
      call locate(edens, 0, nenphdens, Eex, nex2)
      eb = edens(nex2)
      ee = edens(nex2 + 1)
      ldb = phtable2(Zix, Nix, ppi, hpi, pnu, hnu, nex2)
      lde = phtable2(Zix, Nix, ppi, hpi, pnu, hnu, nex2 + 1)
    else
      eb = edens(nenphdens - 1)
      ee = edens(nenphdens)
      ldb = phtable2(Zix, Nix, ppi, hpi, pnu, hnu, nenphdens - 1)
      lde = phtable2(Zix, Nix, ppi, hpi, pnu, hnu, nenphdens)
    endif
    if (ldb > 1..and.lde > 1.) then
      lldb = log(ldb)
      llde = log(lde)
      ldtab = exp(lldb + (Eex - eb) / (ee - eb) * (llde - lldb))
    else
      ldtab = ldb + (Eex - eb) / (ee - eb) * (lde - ldb)
    endif
    phdens2 = ldtab
  endif
  if (phdens2 < 1.e-10) phdens2 = 0.
  return
end function phdens2
! Copyright A.J. Koning 2021
