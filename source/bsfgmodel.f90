function bsfgmodel(Zix, Nix, ald, Eex, P, ibar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Back-shifted Fermi gas level density formula
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               dbl  , &             ! double precision kind
               sgl                  ! single precision kind
!
! *** Declaration of local data
!
  integer   :: ibar      ! fission barrier
  integer   :: Nix       ! neutron number index for residual nucleus
  integer   :: Zix       ! charge number index for residual nucleus
  real(sgl) :: ald       ! level density parameter
  real(sgl) :: an        ! neutron level density parameter
  real(sgl) :: ap        ! proton level density parameter
  real(sgl) :: Eex       ! excitation energy
  real(sgl) :: P         ! pairing energy
  real(sgl) :: sigma     ! help variable
  real(sgl) :: spincut   ! spin cutoff factor
  real(sgl) :: T2        ! square of temperature
  real(sgl) :: U         ! excitation energy minus pairing energy
  real(dbl) :: bsfgmodel ! level density
  real(dbl) :: deninv    ! help variable
  real(dbl) :: expo      ! help variable
  real(dbl) :: fermi     ! function for Fermi gas level density formula
  real(dbl) :: invfermi  ! help variable
  real(dbl) :: term      ! help variable
!
! *********************** Level density formula ************************
!
! fermi    : function for Fermi gas level density formula
!
! Back-shifted Fermi gas
!
! We apply the method of Grossjean and Feldmeier, as implemented by Demetriou and Goriely, to avoid the unphysical divergence
! near zero energy.
! The contribution given by the 1./term goes rapidly to zero with increasing excitation energy.
!
  sigma = sqrt(spincut(Zix, Nix, ald, Eex, ibar, 0))
  an = 0.5 * ald
  ap = 0.5 * ald
  term = exp(1.) / 24. * (an + ap) **2 / sqrt(an * ap) / sigma
  U = Eex - P
  if (U > 0.) then
    invfermi = 1. / fermi(Zix, Nix, ald, Eex, P, ibar)
    T2 = U / ald
    expo = 4. * an * ap * T2
    if (expo < 80.) then
      deninv = invfermi + 1. / (term * exp(expo))
    else
      deninv = invfermi
    endif
    bsfgmodel = 1. / deninv
  else
    bsfgmodel = term
  endif
  return
end function bsfgmodel
! Copyright A.J. Koning 2021
