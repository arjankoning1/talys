function phdens(Zix, Nix, p, h, gs, Eex, Ewell, surfwell)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Particle-hole state density
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
!   edens         ! energy grid for tabulated level densities
! Variables for particle-hole state densities
!   Ephdensmax    ! maximum energy on p - h state density table
!   nenphdens     ! number of energies for p - h state density grid
!   phexist1      ! flag for existence of p - h state
!   phtable1      ! p - h state density
! Variables for preequilibrium initialization
!   Apauli        ! Pauli blocking correction
!   nfac          ! n!
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell   ! flag for surface effects in finite well
  integer   :: h          ! help variable
  integer   :: n1         ! number of coordinate grid points
  integer   :: nex2       ! counter
  integer   :: Nix        ! neutron number index for residual nucleus
  integer   :: p          ! particle number
  integer   :: Zix        ! charge number index for residual nucleus
  real(sgl) :: Ap         ! two-component Pauli blocking correction factor
  real(sgl) :: eb         ! help variable
  real(sgl) :: ee         ! energy
  real(sgl) :: Eex        ! excitation energy
  real(sgl) :: Ewell      ! depth of potential well
  real(sgl) :: fac1       ! help variable
  real(sgl) :: factor     ! multiplication factor
  real(sgl) :: finitewell ! correction function for finite well depth
  real(sgl) :: gs         ! single-particle level density parameter
  real(sgl) :: ldb        ! level density
  real(sgl) :: lde        ! level density
  real(sgl) :: ldtab      ! tabulated level density
  real(sgl) :: lldb       ! log of level density
  real(sgl) :: llde       ! log of level density
  real(sgl) :: phdens     ! function for particle-hole state density
!
! ***** State density of Betak and Dobes, Z. Phys. A279 (1976) 319. ****
!
! finitewell: function for correction for finite well depth
!
! In general, the finite depth of the hole is included.
! If the uncorrected state density is required, it should be specified by Ewell=0.
! (which thus actually means Ewell=infinity).
!
  phdens = 0.
  if (p < 0 .or. h < 0) return
  if (p + h == 0) return
  if (phmodel == 1 .or. .not. phexist1(Zix, Nix, p, h)) then
    Ap = Apauli(p, h)
    factor = (p * p + h * h + p + h) / (4. * gs)
    if (Ap + factor >= Eex) return
    n1 = p + h - 1
    fac1 = nfac(p) * nfac(h) * nfac(n1)
    factor = gs **(p + h) / fac1
    phdens = factor * (Eex - Ap) **n1
    phdens = phdens * finitewell(p, h, Eex, Ewell, surfwell)
  else

! 2. Tabulated particle-hole state densities
!
! locate       : subroutine to find value in ordered table
!
    if (Eex <= 0.) return
    if (Eex <= Ephdensmax) then
      call locate(edens, 0, nenphdens, Eex, nex2)
      eb = edens(nex2)
      ee = edens(nex2 + 1)
      ldb = phtable1(Zix, Nix, p, h, nex2)
      lde = phtable1(Zix, Nix, p, h, nex2 + 1)
    else
      eb = edens(nenphdens - 1)
      ee = edens(nenphdens)
      ldb = phtable1(Zix, Nix, p, h, nenphdens - 1)
      lde = phtable1(Zix, Nix, p, h, nenphdens)
    endif
    if (ldb > 1..and.lde > 1.) then
      lldb = log(ldb)
      llde = log(lde)
      ldtab = exp(lldb + (Eex - eb) / (ee - eb) * (llde - lldb))
    else
      ldtab = ldb + (Eex - eb) / (ee - eb) * (lde - ldb)
    endif
    phdens = ldtab
  endif
  if (phdens < 1.e-10) phdens = 0.
  return
end function phdens
! Copyright A.J. Koning 2021
