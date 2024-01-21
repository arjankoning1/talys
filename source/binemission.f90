subroutine binemission
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Compound emission cross sections for binary reaction
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
!   sgl             ! single precision kind
! All global variables
!   numex           ! maximum number of excitation energies
! Variables for basic reaction
!   flagchannels    ! flag for exclusive channels calculation
! Variables for energy grid
!   deltaE          ! energy bin around outgoing energies
!   ebegin          ! first energy point of energy grid
!   Ebottom         ! bottom of outgoing energy bin
!   egrid           ! outgoing energy grid
!   Etop            ! top of outgoing energy bin
! Variables for energies
!   eend            ! last energy point of energy grid
!   eoutdis         ! outgoing energy of discrete state reaction
!   nendisc         ! last discrete bin
!   speceps         ! limit for cross section spectra
! Variables for excitation energy grid
!   deltaEx         ! excitation energy bin for population arrays
!   Ex              ! excitation energy
!   maxex           ! maximum excitation energy bin for residual nucleus
! Variables for binary emission spectra
!   binemis         ! emission spectra from initial compound nucleus
!   binnorm         ! normalization factor for binary spectra
!   xsemis          ! cross section for emission from compound nucleus
! Variables for compound nucleus from target
!   Exinc           ! excitation energy of entrance bin
! Variables for nuclides
!   Nindex          ! neutron number index for residual nucleus
!   parskip         ! logical to skip outgoing particle
!   Zindex          ! charge number index for residual nucleus
! Variables for level density
!   Nlast           ! last discrete level
! Variables for masses
!   S               ! separation energy
! Variables for incident channel
!   contrib         ! contribution to emission spectrum
!   xscompcont      ! compound cross section for continuum
!   xsgr            ! total smoothed giant resonance cross section
!   xspreeq         ! preeq. cross section per particle typ and outgoing energye
!
! *** Declaration of local data
!
  implicit none
  integer   :: na               ! help variable
  integer   :: nb               ! help variable
  integer   :: nen              ! energy counter
  integer   :: nenbeg           ! help variable
  integer   :: nenend           ! help variable
  integer   :: nexout           ! energy index for outgoing energy
  integer   :: Nix              ! neutron number index for residual nucleus
  integer   :: NL               ! last discrete level
  integer   :: type             ! particle type
  integer   :: Zix              ! charge number index for residual nucleus
  real(sgl) :: dE               ! help variable
  real(sgl) :: dEx              ! excitation energy bin for population arrays
  real(sgl) :: Ea               ! start energy of local adjustment
  real(sgl) :: Eb               ! end energy of local adjustment
  real(sgl) :: emax             ! maximal emission energy within bin decay
  real(sgl) :: emin             ! minimal emission energy
  real(sgl) :: emissum          ! integrated binary emission spectrum
  real(sgl) :: Eo(0:numex)      ! outgoing energy grid based on excitation energy
  real(sgl) :: Eout             ! outgoing energy
  real(sgl) :: Exout            ! excitation energy
  real(sgl) :: frac             ! help variable
  real(sgl) :: fracbot          ! help variable
  real(sgl) :: fractop          ! help variable
  real(sgl) :: SS               ! separation energy
  real(sgl) :: xs               ! help variable
  real(sgl) :: xsa              ! help variable
  real(sgl) :: xsb              ! help variable
  real(sgl) :: xsMeV(0:numex+1) ! cross section on outgoing energy grid
!
! ******************** Calculation of binary spectra *******************
!
! The decay from the primary compound nucleus to residual nuclei is converted from the excitation energy grid to emission energies.
! Decay from the primary compound nucleus to discrete states is already taken into account and the emission energies for these
! transitions are exactly known.
! The emission energies Eo correspond exactly with the excitation energy grid points.
! The contribution per MeV, xsMeV, is constructed.
!
  do type = 0, 6
    if (parskip(type)) cycle
    binnorm(type) = 0.
    if (ebegin(type) >= eend(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    NL = Nlast(Zix, Nix, 0)
    SS = S(0, 0, type)
    if (maxex(Zix, Nix) <= NL) cycle
    do nexout = NL + 1, maxex(Zix, Nix)
      dEx = deltaEx(Zix, Nix, nexout)
      Exout = Ex(Zix, Nix, nexout)
      Eo(nexout) = Exinc - SS - Exout
      xsMeV(nexout) = contrib(type, nexout) / dEx
    enddo
!
! To avoid unphysical interpolations, the contribution for the last discrete level is temporarily set to that of the first
! continuum bin.
!
    Eo(NL) = Exinc - SS - Ex(Zix, Nix, NL)
    xsMeV(NL) = xsMeV(NL + 1)
    xsMeV(maxex(Zix, Nix) + 1) = 0.
    do nen = ebegin(type), eend(type)
      Eout = egrid(nen)
!
! Transitions to discrete states are skipped.
! For the lowest emission energies we ensure that the limit is zero for zero outgoing energy.
!
      if (Eout >= Eo(NL)) cycle
      if (Eout <= Eo(maxex(Zix, Nix))) then
        nexout = maxex(Zix, Nix)
        dE = Eout / Eo(nexout)
        xs = xsMeV(nexout) * dE
        xsemis(type, nen) = xs
        cycle
      endif
!
! To obtain the emission energies belonging to the excitation energy grid points, we use first order interpolation.
!
! locate : subroutine to find value in ordered table
! pol1   : subroutine for interpolation of first order
!
      call locate(Eo, NL, maxex(Zix, Nix), Eout, nexout)
      na = nexout
      nb = nexout + 1
      Ea = Eo(na)
      Eb = Eo(nb)
      xsa = xsMeV(na)
      xsb = xsMeV(nb)
      call pol1(Ea, Eb, xsa, xsb, Eout, xs)
      if (xs < speceps) xs = 0.
      xsemis(type, nen) = xs
    enddo
!
! Normalisation of binary spectra to binary cross sections
!
! Due to interpolation errors, there is always a small difference between the binary continuum cross section and the integral of the
! spectrum over the continuum.
! The spectrum is accordingly normalized.
! The fraction of the first continuum bin is taken into account.
!
    emissum = 0.
    do nen = ebegin(type), nendisc(type)
      emissum = emissum + xsemis(type, nen) * deltaE(nen)
    enddo
    if (eoutdis(type, NL) > 0.) then
      frac = Etop(nendisc(type)) - eoutdis(type, NL)
      emissum = emissum - xsemis(type, nendisc(type)) * frac
    endif
    if (emissum /= 0.) binnorm(type) = xscompcont(type) / emissum
    do nen = ebegin(type), eend(type)
      xsemis(type, nen) = binnorm(type) * xsemis(type, nen)
    enddo
!
! Exclusive channels.
! The formula for exclusive reaction channels shows that we need to store the emission spectrum per excitation energy bin for
! each outgoing particle type. This is done in array binemis.
! First the bottom (emin) and top (emax) of the excitation energy bin are determined.
! Then the begin and end points for the spectrum are located.
! The preequilibrium contribution is also added.
!
    if (flagchannels) then
      do nexout = NL + 1, maxex(Zix, Nix)
        dEx = deltaEx(Zix, Nix, nexout)
        emin = Eo(nexout) - 0.5 * dEx
        emax = Eo(nexout) + 0.5 * dEx
        call locate(Ebottom, ebegin(type), eend(type), emin, nenbeg)
        call locate(Ebottom, ebegin(type), eend(type), emax, nenend)
        nenbeg = max(nenbeg, 1)
        if (nenend < nenbeg) cycle
        do nen = nenbeg, nenend
          if (egrid(nen) >= Eo(NL)) cycle
          binemis(type, nexout, nen) = xsemis(type, nen) + xspreeq(type, nen) + xsgr(type, nen)
        enddo
!
! The end points of the bin are taken into account and this is corrected in binemis.
! There are also cases where the excitation energy bin falls completely in the emission energy bin (nenbeg=nenend).
!
        fracbot = (emin - Ebottom(nenbeg)) / deltaE(nenbeg)
        fractop = (Etop(nenend) - emax) / deltaE(nenend)
        if (nexout == NL + 1) fractop = 0.
        if (nenend == nenbeg) then
          binemis(type, nexout, nenbeg) = binemis(type, nexout, nenbeg) * (1. - fracbot - fractop)
        else
          binemis(type, nexout, nenbeg) = binemis(type, nexout, nenbeg) * (1. - fracbot)
          binemis(type, nexout, nenend) = binemis(type, nexout, nenend) * (1. - fractop)
        endif
      enddo
    endif
  enddo
  return
end subroutine binemission
! Copyright A.J. Koning 2021
