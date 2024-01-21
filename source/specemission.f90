subroutine specemission(Zcomp, Ncomp, nex, idorg, type, nexout)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Exclusive emission cross sections for continuum
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
!   numen       ! maximum number of outgoing energies
! Variables for energy grid
!   deltaE      ! energy bin around outgoing energies
!   ebegin      ! first energy point of energy grid
!   Ebottom     ! bottom of outgoing energy bin
!   egrid       ! outgoing energy grid
!   Etop        ! top of outgoing energy bin
! Variables for energies
!   eend        ! last energy point of energy grid
!   speceps     ! limit for cross section spectra
! Variables for excitation energy grid
!   deltaEx     ! excitation energy bin for population arrays
!   Ex          ! excitation energy
!   maxex       ! maximum excitation energy bin for residual nucleus
!   nexmax      ! maximum excitation energy bin for residual nucleus
! Variables for exclusive channels
!   specemis    ! exclusive emission contribution
!   xsexcl      ! exclusive cross section per excitation energy
! Variables for multiple emission
!   feedexcl    ! feeding terms from compound excitation ene
!   popexcl     ! population cross section of bin just before decay
! Variables for compound nucleus from target
!   dExinc      ! excitation energy bin for mother nucleus
!   Exinc       ! excitation energy of entrance bin
! Variables for nuclides
!   Nindex      ! neutron number index for residual nucleus
!   Zindex      ! charge number index for residual nucleus
! Variables for level density
!   Nlast       ! last discrete level
! Variables for masses
!   S           ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer   :: idorg         ! identifier for previous channel
  integer   :: Ncomp         ! neutron number index for compound nucleus
  integer   :: nen           ! energy counter
  integer   :: nenbeg        ! help variable
  integer   :: nenend        ! help variable
  integer   :: nex           ! excitation energy bin of compound nucleus
  integer   :: nexout        ! energy index for outgoing energy
  integer   :: nexout2       ! energy index for outgoing energy
  integer   :: Nix           ! neutron number index for residual nucleus
  integer   :: NL            ! last discrete level
  integer   :: type          ! particle type
  integer   :: Zcomp         ! proton number index for compound nucleus
  integer   :: Zix           ! charge number index for residual nucleus
  real(sgl) :: dE            ! help variable
  real(sgl) :: dEhalf        ! half of emission bin
  real(sgl) :: dEx           ! excitation energy bin for population arrays
  real(sgl) :: edist         ! help variable
  real(sgl) :: emax          ! maximal emission energy within bin decay
  real(sgl) :: emin          ! minimal emission energy
  real(sgl) :: Eout          ! outgoing energy
  real(sgl) :: Ex0min        ! lower boundary of entrance bin
  real(sgl) :: Ex0plus       ! upper boundary of entrance bin
  real(sgl) :: Ex1min        ! lower boundary of residual bin
  real(sgl) :: Ex1plus       ! upper boundary of residual bin
  real(sgl) :: Exm           ! maximal attainable energy
  real(sgl) :: Exmin         ! help variable
  real(sgl) :: Exout         ! excitation energy
  real(sgl) :: factor        ! multiplication factor
  real(sgl) :: feed0         ! help variable
  real(sgl) :: feedm         ! help variable
  real(sgl) :: feedp         ! help variable
  real(sgl) :: fracbot       ! help variable
  real(sgl) :: fractop       ! help variable
  real(sgl) :: SS            ! separation energy
  real(sgl) :: sumweight     ! help variable
  real(sgl) :: term          ! help variable
  real(sgl) :: term1         ! help variable
  real(sgl) :: term2         ! help variable
  real(sgl) :: weight(numen) ! weight of emission energy bin
  real(sgl) :: xso00         ! cross section on excitation energy grid
  real(sgl) :: xso0m         ! help variable
  real(sgl) :: xso0p         ! help variable
  real(sgl) :: xsom0         ! help variable
  real(sgl) :: xsomp         ! help variable
  real(sgl) :: xsop0         ! help variable
  real(sgl) :: xsopm         ! help variable
!
! ****************** Smearing into emission spectra ********************
!
!
! The decay from compound to residual nuclei is converted from the excitation energy grid to emission energies.
! The spectrum is obtained by spreading the decay over the mother bin and, in the case of continuum-continuum transitions, the
! residual bin.
!
! 1. Loop over excitation energy bins of compound nucleus.
!
! For each mother excitation energy bin, determine the highest possible excitation energy bin nexmax for the residual nuclei.
! As reference, we take the top of the mother bin.
! The maximal residual excitation energy Exm is obtained by subtracting the separation energy from this.
! The bin in which this Exm falls is denoted by nexmax.
!
  do nen = 0, numen
    specemis(nen) = 0.
  enddo
  if (popexcl(Zcomp, Ncomp, nex) <= speceps) return
  Exinc = Ex(Zcomp, Ncomp, nex)
  dExinc = deltaEx(Zcomp, Ncomp, nex)
  Ex0plus = Exinc + 0.5 * dExinc
  Ex0min = Exinc - 0.5 * dExinc
  SS = S(Zcomp, Ncomp, type)
  Exm = Ex0plus - SS
  if (type > 1) Exm = Exm - egrid(ebegin(type))
  Zix = Zindex(Zcomp, Ncomp, type)
  Nix = Nindex(Zcomp, Ncomp, type)
  do nexout2 = maxex(Zix, Nix), 0, - 1
    dEx = deltaEx(Zix, Nix, nexout2)
    Exmin = Ex(Zix, Nix, nexout2) - 0.5 * dEx
    if (Exmin < Exm) then
      nexmax(type) = nexout2
!
! 2. Determine the widths over which the decay must be spread.
!
      Exout = Ex(Zix, Nix, nexout)
      dEx = deltaEx(Zix, Nix, nexout)
      NL = Nlast(Zix, Nix, 0)
!
! Decay from continuum to continuum.
! For most residual continuum bins, no special care needs to be taken and the emission energy Eout that characterizes the
! transition is simply the average between the highest energetic transition that is possible (emax, from the top of the mother bin
! to the bottom of the residual bin) and the lowest (emin).
! However, the highest residual bin (nexout=nexmax) is characterized by different energies (Ex1plus is the maximal residual
! excitation energy and Eout is again the average between emin and emax).
!
      if (nexout > NL) then
        Ex1min = Exout - 0.5 * dEx
        if (nexout == nexmax(type) .and. type >= 1) then
          Ex1plus = Ex0plus - SS
          Exout = 0.5 * (Ex1plus + Ex1min)
        else
          Ex1plus = Exout + 0.5 * dEx
        endif
        emax = Ex0plus - SS - Ex1min
        emin = max(Ex0min - SS - Ex1plus, 0.)
        Eout = 0.5 * (emin + emax)
      else
!
! Decay from continuum to discrete.
! The lowest possible mother excitation bin can not entirely decay to the discrete state.
! For the residual discrete state, it is checked whether the mother excitation bin is such a boundary case.
! This is done by adding the article separation energy to the excitation energy of the residual discrete state.
!
        Exm = Exout + SS
        if (Exm <= Ex0plus .and. Exm > Ex0min) then
          Eout = 0.5 * (Ex0plus + Exm) - SS - Exout
          emin = 0.
        else
          Eout = Exinc - SS - Exout
          emin = Ex0min - SS - Exout
        endif
        emax = Ex0plus - SS - Exout
      endif
!
! 3. Redistribution of decay from population on emission energy grid.
!
! After the bottom (emin) and top (emax) of the excitation energy bin have been determined.
! Then the begin and end points for the spectrum are located.
!
! locate       : subroutine to find value in ordered table
!
      call locate(Ebottom, ebegin(type), eend(type), emin, nenbeg)
      call locate(Ebottom, ebegin(type), eend(type), emax, nenend)
      nenbeg = max(nenbeg, 1)
      if (nenend < nenbeg) return
      dE = emax - emin
      dEhalf = 0.5 * dE
!
! For the interpolation, we first determine the feeding from adjacent mother bins.
! Note that we interpolate the whole product of terms from the exclusive cross section that depend on excitation energy.
!
      feed0 = 0.
      if (popexcl(Zcomp, Ncomp, nex) /= 0.) feed0 = xsexcl(idorg, nex) / popexcl(Zcomp, Ncomp, nex)
      xso00 = feedexcl(Zcomp, Ncomp, type, nex, nexout) * feed0
      if (xso00 == 0.) return
      feedm = 0.
      if (popexcl(Zcomp, Ncomp, nex - 1) /= 0.) feedm = xsexcl(idorg, nex - 1) / popexcl(Zcomp, Ncomp, nex - 1)
      xsom0 = feedexcl(Zcomp, Ncomp, type, nex - 1, nexout) * feedm
      feedp = 0.
      if (popexcl(Zcomp, Ncomp, nex + 1) /= 0.) feedp = xsexcl(idorg, nex + 1) / popexcl(Zcomp, Ncomp, nex + 1)
      xsop0 = feedexcl(Zcomp, Ncomp, type, nex + 1, nexout) * feedp
!
! Decay from continuum to continuum.
! We need to interpolate between adjacent residual bins.
! They are also determined.
! For the lowest emission energies we ensure that the limit is zero for zero outgoing energy.
!
      if (nexout > NL) then
        xso0m = feedexcl(Zcomp, Ncomp, type, nex, nexout - 1) * feed0
        xsopm = feedexcl(Zcomp, Ncomp, type, nex + 1, nexout - 1) * feedp
        xsomp = feedexcl(Zcomp, Ncomp, type, nex - 1, nexout + 1) * feedm
        xso0p = feedexcl(Zcomp, Ncomp, type, nex, nexout + 1) * feed0
        do nen = nenbeg, nenend
!
! A. Emission energy lower than the average emission energy.
!    Interpolate with lower bins.
!
          if (egrid(nen) <= Eout) then
            edist = Eout - egrid(nen)
            if (emin <= 0.0001) then
              term = xso00 * (1. - edist / dEhalf)
            else
              term1 = xso00 + 0.5 * edist / dE * (xsom0 - xso00)
              term2 = xso0p + 0.5 * edist / dE * (xsomp - xso0p)
              term = term1 + 0.5 * edist / dE * (term2 - term1)
            endif
          else
!
! B. Emission energy higher than the average emission energy.
!    Interpolate with higher bins.
!
            edist = egrid(nen) - Eout
            if (xsop0 == 0.) then
              term1 = xso00 * (1. - 0.5 * edist / dEhalf)
              term2 = 0.
            else
              term1 = xso00 + 0.5 * edist / dE * (xsop0 - xso00)
              term2 = xso0m + 0.5 * edist / dE * (xsopm - xso0m)
            endif
            term = term1 + 0.5 * edist / dE * (term2 - term1)
          endif
!
! Each emission energy that is possible within the transition from bin to bin gets a weight.
!
          term = max(term, 0.)
          weight(nen) = term * deltaE(nen)
        enddo
      else
!
! Decay from continuum to discrete.
! We only need to interpolate between adjacent mother bins.
!
        do nen = nenbeg, nenend
          if (egrid(nen) <= Eout) then
            edist = Eout - egrid(nen)
            if (Exm <= Ex0plus .and. Exm > Ex0min) then
              term = xso00 * (1. - edist / dEhalf)
            else
              term = xso00 + edist / dExinc * (xsom0 - xso00)
            endif
          else
            edist = egrid(nen) - Eout
            if (xsop0 == 0.) then
              term = xso00 * (1. - edist / dEhalf)
            else
              term = xso00 + edist / dExinc * (xsop0 - xso00)
            endif
          endif
          term = max(term, 0.)
          weight(nen) = term * deltaE(nen)
        enddo
      endif
!
! The end points of the bin are taken into account and this is corrected in weight.
! There are also cases where the excitation energy bin falls completely in the emission energy bin (nenbeg=nenend).
!
      fracbot = (emin - Ebottom(nenbeg)) / deltaE(nenbeg)
      fractop = (Etop(nenend) - emax) / deltaE(nenend)
      if (nexout == NL + 1) fractop = 0.
      if (nenend == nenbeg) then
        weight(nenbeg) = weight(nenbeg) * (1. - fracbot - fractop)
      else
        weight(nenbeg) = weight(nenbeg) * (1. - fracbot)
        weight(nenend) = weight(nenend) * (1. - fractop)
      endif
!
! The weights are summed to ensure that we get exact flux conservation by normalizing with sumweight.
!
      sumweight = 0.
      do nen = nenbeg, nenend
        sumweight = sumweight + weight(nen)
      enddo
      if (sumweight == 0.) return
!
! The contribution for mother bin --> residual bin is added to the spectrum.
!
      do nen = nenbeg, nenend
        factor = xso00 / sumweight
        specemis(nen) = factor * weight(nen) / deltaE(nen)
      enddo
      exit
    endif
  enddo
  return
end subroutine specemission
! Copyright A.J. Koning 2021
