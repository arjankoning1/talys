function fstrength(Zcomp, Ncomp, Efs, Egamma, irad, l)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Gamma ray strength functions
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
! All global variables
!   numgamqrpa    ! number of energies for QRPA strength function
! Variables for gamma rays
!   egr           ! energy of GR
!   egradjust     ! adjustable factor for energy of GR
!   epr           ! energy of PR
!   epradjust     ! adjustable factor for energy of PR
!   Exlfile       ! tabulated gamma ray strength function
!   flagupbend    ! flag for low-energy correction to photon strength function
!   gamadjust     ! logical for energy - dependent gamma adjustment
!   ggr           ! width of GR
!   ggradjust     ! adjustable factor for width of GR
!   gpr           ! width of PR
!   gpradjust     ! adjustable factor for width of PR
!   sgr           ! strength of GR
!   sgradjust     ! adjustable factor for strength of GR
!   strength      ! E1 strength function model
!   strengthM1    ! model for M1 gamma - ray strength function
!   tpr           ! strength of PR
!   tpradjust     ! adjustable factor for strength of PR
!   upbend        ! properties of the low - energy upbend of given multipolarity
! Variables for main input
!   k0            ! index of incident particle
! Variables for level density
!   alev          ! level density parameter
! Variables for energy grid
!   Einc          ! incident energy in MeV
! Constants
!   twopi         ! 2 * pi
!  Variables for gamma-ray strength functions
!   eqrpa         ! energy grid for QRPA strength function
!   fqrpa         ! tabulated QRPA strength function
!   kgr           ! constant for gamma - ray strength function
!   ngr           ! number of GR
!   nTqrpa        ! number of temperatures for QRPA
!   qrpaexist     ! flag for existence of tabulated QRPA strength func.
!   Tqrpa         ! temperature for QRPA
! Variables for level density
!   delta         ! energy shift
! Variables for masses
!   S             ! separation energy
! Variables for masses
!   beta2         ! deformation parameter
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key                 ! keyword
  integer           :: i                   ! level
  integer           :: irad                ! variable to indicate M(=0) or E(=1) radiation
  integer           :: it                  ! counter for tritons
  integer           :: itemp               ! end of do loop
  integer           :: jt                  ! temperature index
  integer           :: l                   ! multipolarity
  integer           :: Ncomp               ! neutron number index for compound nucleus
  integer           :: nen                 ! energy counter
  integer           :: nT                  ! temperature index
  integer           :: nT0                 ! temperature index
  integer           :: Zcomp               ! proton number index for compound nucleus
  real(sgl)         :: denom               ! help variable
  real(sgl)         :: e                   ! energy
  real(sgl)         :: eb                  ! help variable
  real(sgl)         :: ee                  ! energy
  real(sgl)         :: Efs                 ! fast particle energy for gamma ray strength function
  real(sgl)         :: Egam2               ! help variable
  real(sgl)         :: Egamma              ! gamma energy
  real(sgl)         :: egr1                ! energy of GR
  real(sgl)         :: egr2                ! help variable
  real(sgl)         :: enum                ! enumerator of Lorentzian
  real(sgl)         :: epr1                ! energy of PR
  real(sgl)         :: epr2                ! energy of PR
  real(sgl)         :: Eq(0:numgamqrpa)    ! energy grid for QRPA strength function
  real(sgl)         :: f1                  ! help variable
  real(sgl)         :: f2                  ! help variable
  real(sgl)         :: factor1             ! help variable
  real(sgl)         :: factor2             ! help variable
  real(sgl)         :: fb                  ! help variable
  real(sgl)         :: fe                  ! help variable
  real(sgl)         :: fstrength           ! gamma ray strength function
  real(sgl)         :: gamb                ! tabulated QRPA strength function
  real(sgl)         :: game                ! tabulated QRPA strength function
  real(sgl)         :: ggr1                ! width of GR
  real(sgl)         :: ggr2                ! help variable
  real(sgl)         :: ggredep             ! energy dependent damping width
  real(sgl)         :: ggredep0            ! energy dependent damping width at zero gamma energy
  real(sgl)         :: gpr1                ! width of PR
  real(sgl)         :: gpr2                ! width of PR
  real(sgl)         :: kgr1                ! constant for gamma-ray strength function
  real(sgl)         :: sgr1                ! strength of GR
  real(sgl)         :: Tb                  ! temperature for QRPA
  real(sgl)         :: Te                  ! temperature for QRPA
  real(sgl)         :: Tnuc                ! nuclear temperature
  real(sgl)         :: tpr1                ! strength of PR
  real(sgl)         :: upbendc             ! properties of the low-energy upbend of given multipolarity
  real(sgl)         :: upbende             ! properties of the low-energy upbend of given multipolarity
  real(sgl)         :: upbendf             ! properties of the low-energy upbend of given multipolarity
!
! ************************* Strength functions *************************
!
! fstrength: gamma-ray strength function
! adjust   : subroutine for energy-dependent parameter adjustment
! kgr,kgr1 : constant for gamma-ray strength function
! strength : strength function of Kopecky-Uhl (1) or Brink-Axel (2)
!   or microscopic tables (>=3)
!
! Models for E1 gamma-ray strength function:
!
! 1. Kopecky-Uhl
! 2. Brink-Axel
! 3. Goriely HFbcs tables
! 4. Goriely HFB tables
! 5. Goriely Hybrid model
! 6. Goriely T-dep. HFB Tables
! 7. Goriely T-dep. RMF Tables
! 8. Gogny D1M HFB+QRPA Tables
! 9. IAEA-CRP SMLO 2019 Tables
! 10. BSk27+QRPA 2018 Tables
! 11. D1M-Intra-E1
! 12. Shellmodel-E1
!
  fstrength = 0.
  do i = 1, ngr(Zcomp, Ncomp, irad, l)
   if (gamadjust(Zcomp, Ncomp)) then
      key = 'sgr'
      call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor1)
      key = 'sgradjust'
      call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor2)
      sgr1 = factor1 * factor2 * sgr(Zcomp, Ncomp, irad, l, i) * sgradjust(Zcomp, Ncomp, irad, l, i)
      key = 'ggr'
      call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor1)
      key = 'ggradjust'
      call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor2)
      ggr1 = factor1 * factor2 * ggr(Zcomp, Ncomp, irad, l, i) * ggradjust(Zcomp, Ncomp, irad, l, i)
      key = 'egr'
      call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor1)
      key = 'egradjust'
      call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor2)
      egr1 = factor1 * factor2 * egr(Zcomp, Ncomp, irad, l, i) * egradjust(Zcomp, Ncomp, irad, l, i)
    else
      sgr1 = sgr(Zcomp, Ncomp, irad, l, i)
      egr1 = egr(Zcomp, Ncomp, irad, l, i)
      ggr1 = ggr(Zcomp, Ncomp, irad, l, i)
    endif
    kgr1 = kgr(l)
    egr2 = egr1 **2
    ggr2 = ggr1 **2
    Egam2 = Egamma **2
!
! 1. Kopecky-Uhl generalized Lorentzian.
!
! Efs      : fast particle energy for gamma ray strength function
! qrpaexist: flag for existence of tabulated QRPA strength functions
!
    Tnuc = 0.
    if (strength == 1 .and. l == 1 .and. irad == 1) then
      if (Egamma /= Einc) then
        e = min(Efs, 20.) + S(Zcomp, Ncomp, k0) - delta(Zcomp, Ncomp, 0) - Egamma
        if (e > 0..and.alev(Zcomp, Ncomp) > 0.) Tnuc = sqrt(e / alev(Zcomp, Ncomp))
      endif
      ggredep0 = ggr1 * twopi **2 * Tnuc **2 / egr2
      ggredep = ggredep0 + ggr1 * Egam2 / egr2
      enum = ggredep * Egamma
      denom = (Egam2 - egr2) **2 + Egam2 * ggredep **2
      factor1 = enum / denom
      factor2 = 0.7 * ggredep0 / (egr1 **3)
      fstrength = fstrength + kgr1 * sgr1 * ggr1 * (factor1 + factor2)
    endif
!
! 2. Brink-Axel standard Lorentzian.
!
    if (strength == 2 .or. ((strength == 3 .or. strength == 4 .or. strength >= 6) .and. &
      .not. qrpaexist(Zcomp, Ncomp, 1, 1)) .or. l /= 1 .or. irad /= 1) then
      enum = 0.
      if (Egamma > 0.001) then
        enum = ggr2 * Egamma **(3 - 2 * l)
        denom = (Egam2 - egr2) **2 + Egam2 * ggr2
        fstrength = fstrength + kgr1 * sgr1 * enum / denom
      endif
    endif
!
! 3+4+6+7+8+9+10+11+12. Tabulated QRPA strength functions
!
! locate    : subroutine to find value in ordered table
! numgamqrpa: number of energies for QRPA strength function
! eqrpa,Eq  : energy grid for QRPA strength function
! fqrpa     : tabulated QRPA strength function
! gamb      : tabulated QRPA strength function
! game      : tabulated QRPA strength function
!
    if ((strength == 3 .or. strength == 4 .or. strength >= 6 .or. Exlfile(Zcomp, Ncomp, 1, 1)(1:1) /= ' ') .and. &
      ((qrpaexist(Zcomp, Ncomp, 1, 1) .and. l == 1 .and. irad == 1) .or. &
      (qrpaexist(Zcomp, Ncomp, 0, 1) .and. strengthM1 >= 8 .and. l == 1 .and. irad == 0))) then
      nT0 = nTqrpa
      if (irad /= 1 .or. l /= 1) nT0 = 1
      if ((irad == 1 .and. strength == 11) .or. (irad == 0 .and. strengthM1 == 11)) then
        if (Zcomp == 0 .and. Ncomp == 0 .and. flagupbend) then
          nTqrpa=11
        else
          nTqrpa=1
        endif
        nT0=nTqrpa
      endif
      if (Egamma /= Einc .and. nT0 > 1) then
        e = min(Efs, 20.) + S(Zcomp, Ncomp, k0) - delta(Zcomp, Ncomp, 0) - Egamma
        if (e > 0..and.alev(Zcomp, Ncomp) > 0.) Tnuc = sqrt(e / alev(Zcomp, Ncomp))
        if ((irad == 1 .and. strength == 11) .or. (irad == 0 .and. strengthM1 == 11)) Tnuc = Efs + S(Zcomp, Ncomp, k0)
        nT = nT0
        do it = 1, nT0
          if (Tqrpa(it) > Tnuc) then
            nT = it - 1
            exit
          endif
        enddo
        Tb = Tqrpa(nT)
        if (nT < nT0) then
          Te = Tqrpa(nT + 1)
        else
          Te = Tb
        endif
        itemp = 2
      else
        Tb = 0.
        Te = 0.
        itemp = 1
        nT = 1
      endif
      do nen = 0, numgamqrpa
        Eq(nen) = eqrpa(Zcomp, Ncomp, nen, irad, 1)
      enddo
      if ((Egamma > Tnuc + 0.1 .and. strength.eq.11) .or. (Egamma < Eq(1) .and. strength.ne.11)) then
        fstrength = 0.
      else
        do it = 1, itemp
          jt = nT
          if (it == 2) jt = nT + 1
          if (jt > nT0) jt = nT0
          if (Egamma <= Eq(numgamqrpa)) then
            call locate(Eq, 0, numgamqrpa, Egamma, nen)
            eb = Eq(nen)
            ee = Eq(nen + 1)
            gamb = fqrpa(Zcomp, Ncomp, nen, jt, irad, 1)
            game = fqrpa(Zcomp, Ncomp, nen + 1, jt, irad, 1)
            if (gamb > 0..and.game > 0.) then
              f1 = log10(gamb) + (Egamma - eb) / (ee - eb) * (log10(game) - log10(gamb))
              f2 = 10. **f1
            else
              f2 = gamb + (Egamma - eb) / (ee - eb) * (game - gamb)
            endif
          else
            eb = Eq(numgamqrpa - 1)
            ee = Eq(numgamqrpa)
            gamb = fqrpa(Zcomp, Ncomp, numgamqrpa - 1, jt, irad, 1)
            game = fqrpa(Zcomp, Ncomp, numgamqrpa, jt, irad, 1)
            if (gamb > 0..and.game > 0.) then
              f1 = log10(gamb) + (Egamma - eb) / (ee - eb) * (log10(game) - log10(gamb))
              f2 = 10. **f1
            else
              f2 = gamb + (Egamma - eb) / (ee - eb) * (game - gamb)
            endif
          endif
          if (it == 1) fb = f2
          if (it == 2) fe = f2
        enddo
        if (nT0 > 1 .and. Tb /= Te) then
          if (fb > 0..and.fe > 0.) then
            f1 = log10(fb) + (Tnuc - Tb) / (Te - Tb) * (log10(fe) - log10(fb))
            f2 = 10. **f1
          else
            f2 = fb + (Tnuc - Tb) / (Te - Tb) * (fe - fb)
          endif
        endif
        fstrength = f2
      endif
    endif
!
! 5. Goriely Hybrid model
!
    if (strength == 5 .and. Exlfile(Zcomp, Ncomp, 1, 1)(1:1) == ' ' .and. l == 1 .and. irad == 1) then
      if (Egamma /= Einc) then
        e = min(Efs, 20.) + S(Zcomp, Ncomp, k0) - delta(Zcomp, Ncomp, 0) - Egamma
        if (e > 0..and.alev(Zcomp, Ncomp) > 0.) Tnuc = sqrt(e / alev(Zcomp, Ncomp))
      endif
      if (Egamma > 0.) then
        ggredep = 0.7 * ggr1 * (Egamma / egr1 + twopi **2 * Tnuc **2 / Egamma / egr1)
        enum = ggredep * Egamma
        denom = (Egam2 - egr2) **2 + Egam2 * ggredep * ggr1
        factor1 = enum / denom
        fstrength = fstrength + kgr1 * sgr1 * ggr1 * factor1
      endif
    endif
  enddo
!
! Inclusion of additional extra strength (Pygmy Resonance),
! only if explicitly specified in the input
!
  do i= 1, 2
    tpr1 = tpr(Zcomp, Ncomp, irad, l, i)
    if (Egamma > 0.001 .and. tpr1 > 0.) then
      if (gamadjust(Zcomp, Ncomp)) then
        key = 'tpr'
        call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor1)
        key = 'tpradjust'
        call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor2)
        tpr1 = factor1 * factor2 * tpr(Zcomp, Ncomp, irad, l, i) * tpradjust(Zcomp, Ncomp, irad, l, i)
        key = 'gpr'
        call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor1)
        key = 'gpradjust'
        call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor2)
        gpr1 = factor1 * factor2 * gpr(Zcomp, Ncomp, irad, l, i) * gpradjust(Zcomp, Ncomp, irad, l, i)
        key = 'epr'
        call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor1)
        key = 'epradjust'
        call adjust(Egamma, key, Zcomp, Ncomp, 0, l, factor2)
        epr1 = factor1 * factor2 * epr(Zcomp, Ncomp, irad, l, i) * epradjust(Zcomp, Ncomp, irad, l, i)
      else
        epr1 = epr(Zcomp, Ncomp, irad, l, i)
        gpr1 = gpr(Zcomp, Ncomp, irad, l, i)
      endif
      kgr1 = kgr(l)
      epr2 = epr1 **2
      gpr2 = gpr1 **2
      Egam2 = Egamma **2
      enum = gpr2 * Egamma **(3 - 2 * l)
      denom = (Egam2 - epr2) **2 + Egam2 * gpr2
      fstrength = fstrength + kgr1 * tpr1 * enum / denom
    endif
  enddo
!
! Inclusion of an additional low-E limit of E1 nature
!
  if (flagupbend) then
    upbendc = upbendadjust(Zcomp, Ncomp, irad, l, 1) * upbend(Zcomp, Ncomp, irad, l, 1)
    upbende = upbendadjust(Zcomp, Ncomp, irad, l, 2) * upbend(Zcomp, Ncomp, irad, l, 2)
    upbendf = upbendadjust(Zcomp, Ncomp, irad, l, 3) * upbend(Zcomp, Ncomp, irad, l, 3)
    if ((strengthM1 == 8 .or. strengthM1 == 10) .and. irad == 0 .and. l == 1 .and. Zcomp + Ncomp >= 105) upbendf=0.
    if (irad == 1 .and. l == 1) then
      e = min(Efs, 20.) + S(Zcomp, Ncomp, k0) 
      if (e > 1) fstrength = fstrength + upbendc * e / (1. + exp(Egamma - upbende))
    endif
!
! Inclusion of an additional low-E upbend of M1 nature
!
    if (irad == 0)  fstrength = fstrength + upbendc * exp(-upbende * Egamma) * exp(-upbendf * abs(beta2(Zcomp, Ncomp, 0)))
  endif
  return
end function fstrength
! Copyright A.J. Koning 2021
