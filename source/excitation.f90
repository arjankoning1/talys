subroutine excitation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Excitation energy population
!
! Author    : Arjan Koning and Stephane Hilaire
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
! All global variables
!   numex         ! maximum number of excitation energies
!   numJ          ! maximum J - value
!   numpop        ! number of population bins
! Variables for input energies
!   EdistE        ! excitation energy of population distribution
!   flagpopMeV    ! flag to use initial population per MeV instead of histogram
!   npopE         ! number of energies for population distribution
!   npopJ         ! number of spins for population distribution
!   PdistE        ! population distribution, spin - independent
!   PdistJP       ! population distribution per spin and parity
! Variables for excitation energy grid
!   deltaEx       ! excitation energy bin for population arrays
!   Ex            ! excitation energy
!   maxex         ! maximum excitation energy bin for residual nucleus
!   maxJ          ! maximal J - value
! Variables for multiple emission
!   feedexcl      ! feeding terms from compound excitation ene
!   popexcl       ! population cross section of bin just before decay
!   xsinitpop     ! initial population cross section
! Variables for incident channel
!   xspop         ! population cross section
!   xspopex       ! population cross section summed over spin and parity
!   xspopexP      ! population cross section per parity
!   xspopnuc      ! population cross section per nucleus
!   xspopnucP     ! population cross section per nucleus per parity
! Variables for nuclides
!   AA            ! mass number of residual nucleus
! Constants
!   pardis        ! parity distribution
! Variables for levels
!   jdis          ! spin of level
!   parlev        ! parity of level
! Variables for level density
!   Nlast         ! last discrete level
!
! *** Declaration of local data
!
  implicit none
  integer   :: A                            ! mass number of target nucleus
  integer   :: J                            ! spin of level
  integer   :: Ncomp                        ! neutron number index for compound nucleus
  integer   :: nen                          ! energy counter
  integer   :: nex                          ! excitation energy bin of compound nucleus
  integer   :: nex0                         ! base number for discrete level
  integer   :: NL                           ! last discrete level
  integer   :: odd                          ! odd (1) or even (0) nucleus
  integer   :: parity                       ! parity
  integer   :: Zcomp                        ! proton number index for compound nucleus
  real(sgl) :: ald                          ! level density parameter
  real(sgl) :: dE                           ! help variable
  real(sgl) :: dElow                        ! energy increment at lower bound
  real(sgl) :: dEup                         ! energy increment at upper bound
  real(sgl) :: dEx                          ! excitation energy bin for population arrays
  real(sgl) :: Ea                           ! start energy of local adjustment
  real(sgl) :: Eb                           ! end energy of local adjustment
  real(sgl) :: Edistmax                     ! excitation energy with maximum population
  real(sgl) :: Edistlow(0:numpop)           ! lower bound of energy bin of population distribution
  real(sgl) :: Edistup(0:numpop)            ! upper bound of energy bin of population distribution
  real(sgl) :: Edlow                        ! help variable
  real(sgl) :: Edup                         ! help variable
  real(sgl) :: Eex                          ! excitation energy
  real(sgl) :: Eexlow(0:numex)              ! lower bound excitation energy
  real(sgl) :: Eexup(0:numex)               ! upper bound excitation energy
  real(sgl) :: Eexmax                       ! maximal excitation energy
  real(sgl) :: Eexmin                       ! minimal excitation energy
  real(sgl) :: Exlow                        ! lower bound excitation energy
  real(sgl) :: Exup                         ! upper bound excitation energy
  real(sgl) :: factor                       ! help variable
  real(sgl) :: frac                         ! help variable
  real(sgl) :: ignatyuk                     ! function for energy dependent level density parameter a
  real(sgl) :: normJ                        ! normalization for spin distribution
  real(sgl) :: Pa                           ! probability distribution value
  real(sgl) :: Pb                           ! probability distribution value
  real(sgl) :: Pex(0:numex)                 ! probability at excitation energy
  real(sgl) :: PexJP(0:numex, 0:numJ, -1:1) ! probability at excitation energy, J and P
  real(sgl) :: Prob                         ! probability
  real(sgl) :: Probex                       ! probability
  real(sgl) :: Rspin                        ! residual spin
  real(sgl) :: spindis                      ! Wigner spin distribution
  real(sgl) :: sumJP                        ! help variable
  real(sgl) :: sumPex                       ! help variable
  real(sgl) :: xsinputpop                   ! total population given in input
  real(sgl) :: spincut                      ! spin cutoff factor
  real(sgl) :: sc                           ! spin cutoff factor
!
! ******************** Fill energy bins with population ****************
!
! locate   : subroutine to find value in ordered table
! ignatyuk : function for energy dependent level density parameter a residual excitation energy bin
!
! Write initial population by summing the bins of the input energy grid.
!
  Zcomp = 0
  Ncomp = 0
  if ( .not. flagpopMeV) then
    xsinputpop = 0.
    if (npopJ > 0) then
      do parity = - 1, 1, 2
        do J = 0, npopJ - 1
          do nen = 1, npopE
            xsinputpop = xsinputpop + PdistJP(nen, J, parity)
          enddo
        enddo
      enddo
    else
      Edistmax = 10.
      do nen = 1, npopE
        xsinputpop = xsinputpop + PdistE(nen)
        if (PdistE(nen) > PdistE(nen-1)) Edistmax = EdistE(nen)
      enddo
      ald = ignatyuk(Zcomp, Ncomp, Edistmax, 0)
      sc = spincut(Zcomp, Ncomp, ald, Edistmax, 0, 1)
    endif
    write(*, '(/" Total population of input excitation energy grid:", es12.5/)') xsinputpop
    if (xsinputpop == 0.) return
  endif
!
! Set population bins on basis of input
!
  if (EdistE(1) == 0.) then
    do nen = 1, npopE
      EdistE(nen - 1) = EdistE(nen)
      PdistE(nen - 1) = PdistE(nen)
      do parity = - 1, 1, 2
        do J = 0, npopJ - 1
          PdistJP(nen - 1, J, parity) = PdistJP(nen, J, parity)
        enddo
      enddo
    enddo
    npopE = npopE - 1
  endif
  if (npopJ > 0) then
    do nen = 1, npopE
      PdistE(nen) = 0.
      do parity = - 1, 1, 2
        do J = 0, npopJ - 1
          PdistE(nen) = PdistE(nen) + PdistJP(nen, J, parity)
        enddo
      enddo
    enddo
  endif
  do nen = 0, npopE
    Edistlow(nen) = 0.5 * (EdistE(nen) + EdistE(max(nen - 1, 0)))
    Edistup(nen) = 0.5 * (EdistE(nen) + EdistE(min(npopE, nen + 1)))
  enddo
  Edistlow(0) = 0.
  Edistup(npopE) = EdistE(npopE)
!
! Set boundaries of excitation energy bins
!
  NL = Nlast(Zcomp, Ncomp, 0)
  do nex = 0, maxex(Zcomp, Ncomp)
    Eex = Ex(Zcomp, Ncomp, nex)
    if (nex <= NL) then
      Eexmin = Ex(Zcomp, Ncomp, max(nex - 1, 0))
      Eexlow(nex) = 0.5 * (Eex + Eexmin)
      Eexmax = Ex(Zcomp, Ncomp, min(nex + 1, maxex(Zcomp, Ncomp)))
      Eexup(nex) = 0.5 * (Eex + Eexmax)
    else
      Eexlow(nex) = Eex - 0.5 * deltaEx(Zcomp, Ncomp, nex)
      Eexup(nex) = Eex + 0.5 * deltaEx(Zcomp, Ncomp, nex)
    endif
  enddo
!
! Redistribute bins
!
  A = AA(Zcomp, Ncomp, 0)
  odd = mod(A, 2)
  nex0 = maxex(Zcomp, Ncomp) + 1
  sumPex = 0.
  do nex = 0, maxex(Zcomp, Ncomp)
    Eex = Ex(Zcomp, Ncomp, nex)
    dEx = deltaEx(Zcomp, Ncomp, nex)
    if ( .not. flagpopMeV) then
      Pex(nex) = 0.
      do parity = - 1, 1, 2
        do J = 0, maxJ(Zcomp, Ncomp, nex)
          PexJP(nex, J, parity) = 0.
        enddo
      enddo
      Exlow = Eexlow(nex)
      Exup = Eexup(nex)
      do nen = 0, npopE - 1
        Edlow = Edistlow(nen)
        Edup = Edistup(nen)
        if (Edlow >= Exup) cycle
        if (Edup <= Exlow) cycle
        dE = Edup - Edlow
        dElow = max(Exlow - Edlow, 0.)
        dEup = max(Edup - Exup, 0.)
        frac = (dE - dElow - dEup) / dE
        Pex(nex) = Pex(nex) + frac * PdistE(nen)
        if (npopJ > 0 .and. nex > NL) then
          do parity = - 1, 1, 2
            do J = 0, maxJ(Zcomp, Ncomp, nex)
              PexJP(nex, J, parity) = PexJP(nex, J, parity) + frac * PdistJP(nen, J, parity)
            enddo
          enddo
       endif
      enddo
      if (nex <= NL) then
        J = int(jdis(Zcomp, Ncomp, nex))
        parity = parlev(Zcomp, Ncomp, nex)
        PexJP(nex, J, parity) = Pex(nex)
      else
        if (npopJ == 0) then
          do parity = - 1, 1, 2
            normJ = 0.
            do J = 0, maxJ(Zcomp,Ncomp,nex)
              Rspin = real(J) + 0.5 * odd
              normJ = normJ + spindis(sc, Rspin)
            enddo
            do J = 0, maxJ(Zcomp, Ncomp, nex)
              Rspin = real(J) + 0.5 * odd
              PexJP(nex, J, parity) = Pex(nex) * pardis * spindis(sc, Rspin) / normJ
            enddo
          enddo
        endif
      endif
      sumPex = sumPex + Pex(nex)
    else
      call locate(EdistE, 0, npopE, Eex, nen)
      nen = max(0, nen)
      Ea = EdistE(nen)
      Eb = EdistE(nen + 1)
      if (npopJ == 0) then
        Pa = PdistE(nen)
        Pb = PdistE(nen + 1)
        call pol1(Ea, Eb, Pa, Pb, Eex, Probex)
      endif
      do parity = - 1, 1, 2
        do J = 0, maxJ(Zcomp, Ncomp, nex)
          Rspin = real(J) + 0.5 * odd
          if (nex <= NL .and. (jdis(Zcomp, Ncomp, nex) /= J .or. parlev(Zcomp, Ncomp, nex) /= parity)) cycle
          if (npopJ == 0) then
            Prob = Probex * spindis(sc, Rspin) * pardis
          else
            Pa = PdistJP(nen, J, parity)
            Pb = PdistJP(nen + 1, J, parity)
            call pol1(Ea, Eb, Pa, Pb, Eex, Prob)
          endif
          xspop(Zcomp, Ncomp, nex, J, parity) = Prob * dEx
          xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) + xspop(Zcomp, Ncomp, nex, J, parity)
          xspopexP(Zcomp, Ncomp, nex, parity) = xspopexP(Zcomp, Ncomp, nex, parity) + xspop(Zcomp, Ncomp, nex, J, parity)
        enddo
        xspopnucP(Zcomp, Ncomp, parity) = xspopnucP(Zcomp, Ncomp, parity) + xspopexP(Zcomp, Ncomp, nex, parity)
      enddo
      xspopnuc(Zcomp, Ncomp) = xspopnuc(Zcomp, Ncomp) + xspopex(Zcomp, Ncomp, nex)
      feedexcl(Zcomp, Ncomp, 0, nex0, nex) = xspopex(Zcomp, Ncomp, nex)
    endif
  enddo
  if ( .not. flagpopMeV .and. sumPex > 0.) then
    factor = sumPex / xsinputpop
    sumJP = 0.
    do nex = 0, maxex(Zcomp, Ncomp)
      do parity = - 1, 1, 2
        do J = 0, maxJ(Zcomp, Ncomp, nex)
          xspop(Zcomp, Ncomp, nex, J, parity) = PexJP(nex, J, parity) / factor
          sumJP = sumJP + xspop(Zcomp, Ncomp, nex, J, parity)
        enddo
      enddo
    enddo
    do nex = 0, maxex(Zcomp, Ncomp)
      factor = sumJP / xsinputpop
      do parity = - 1, 1, 2
        do J = 0, maxJ(Zcomp, Ncomp, nex)
          xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) + xspop(Zcomp, Ncomp, nex, J, parity) / factor
          xspopexP(Zcomp, Ncomp, nex, parity) = xspopexP(Zcomp, Ncomp, nex, parity) + &
            xspop(Zcomp, Ncomp, nex, J, parity) / factor
        enddo
        xspopnucP(Zcomp, Ncomp, parity) = xspopnucP(Zcomp, Ncomp, parity) + xspopexP(Zcomp, Ncomp, nex, parity)
      enddo
      xspopnuc(Zcomp, Ncomp) = xspopnuc(Zcomp, Ncomp) + xspopex(Zcomp, Ncomp, nex)
      feedexcl(Zcomp, Ncomp, 0, nex0, nex) = xspopex(Zcomp, Ncomp, nex)
    enddo
  endif
!
! Special case for one population energy
!
  if (npopE == 1) then
    Eex = EdistE(1)
    nex = maxex(Zcomp, Ncomp)
    do parity = - 1, 1, 2
      if (npopJ == 0) then
        normJ = 0.
        do J = 0, maxJ(Zcomp,Ncomp,nex)
          Rspin = real(J) + 0.5 * odd
          normJ = normJ + spindis(sc, Rspin)
        enddo
      endif
      do J = 0, maxJ(Zcomp, Ncomp, nex)
        if (npopJ == 0) then
          Rspin = real(J) + 0.5 * odd
          xspop(Zcomp, Ncomp, nex, J, parity) = PdistE(1) * spindis(sc, Rspin) * pardis / normJ
        else
          xspop(Zcomp, Ncomp, nex, J, parity) = PdistJP(1, J, parity)
        endif
        xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) + xspop(Zcomp, Ncomp, nex, J, parity)
        xspopexP(Zcomp, Ncomp, nex, parity) = xspopexP(Zcomp, Ncomp, nex, parity) + xspop(Zcomp, Ncomp, nex, J, parity)
      enddo
      xspopnucP(Zcomp, Ncomp, parity) = xspopnucP(Zcomp, Ncomp, parity) + xspopexP(Zcomp, Ncomp, nex, parity)
    enddo
    xspopnuc(Zcomp, Ncomp) = xspopnuc(Zcomp, Ncomp) + xspopex(Zcomp, Ncomp, nex)
    feedexcl(Zcomp, Ncomp, 0, nex0, nex) = xspopex(Zcomp, Ncomp, nex)
  endif
  xsinitpop = xspopnuc(Zcomp, Ncomp)
  popexcl(Zcomp, Ncomp, nex0) = xsinitpop
  return
end subroutine excitation
! Copyright A.J. Koning 2021
