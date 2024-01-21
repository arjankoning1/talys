subroutine prodyield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculate production yields
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
!   dbl             ! double precision kind
! All global variables
!   numisom         ! number of isomers
!   numN            ! maximum number of neutrons from initial compound nucleus
!   numtime         ! number of time points
!   numZ            ! maximum number of protons from initial compound nucleus
! Variables for medical isotope production
!   Ibeam           ! beam current in mA for isotope production
!   radiounit       ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi, Ci or kCi
!   rhotarget       ! target material density
!   yieldunit       ! unit for isotope yield: num (number), mug, mg, g, or kg
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   Ninit           ! neutron number of initial compound nucleus
!   Zinit           ! charge number of initial compound nucleus
!   Ztarget         ! charge number of target nucleus
! Constants
!   avogadro        ! Avogadro's number
! Variables for levels
!   Nisomer         ! number of isomers for this nuclide
! Variables for decay data
!   daysec          ! number of seconds in a day
!   hoursec         ! number of seconds in an hour
!   lambda          ! decay rate per isotope
!   minutesec       ! number of seconds in a minute
!   rtyp            ! type of beta decay, beta - : 1 , beta + : 2 (from ENDF format)
!   yearsec         ! number of seconds in a year
! Variables for isotope production
!   activity        ! activity of produced isotope in MBq
!   Niso            ! number of isotopes produced after irradiation
!   Nisorel         ! fraction of number of produced isotopes per element
!   Nisotot         ! number of elemental isotopes produced after irradiation
!   Ntar0           ! number of original target atoms
!   Ntime           ! number of time points
!   prate           ! production rate per isotope
!   Tco             ! cooling time per unit
!   Tgrid           ! time
!   Tir             ! irradiation time per unit
!   Tmax            ! irradiation time with maximal yield
!   Tmaxactivity    ! time of maximum activity of produced isotope in MBq
!   Tp              ! irradiation time with maximal yield per time unit
!   Vtar            ! active target volume
!   yield           ! yield of produced isotope in MBq / (mA.h)
! Variables for existence libraries
!   Yexist          ! flag for existence of yield
!
! *** Declaration of local data
!
  implicit none
  integer   :: A       ! mass number of target nucleus
  integer   :: is      ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: isob    ! counter
  integer   :: it      ! counter for tritons
  integer   :: N       ! neutron number of residual nucleus
  integer   :: Ncool   ! number of cooling time steps
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: Nparent ! N of parent isotope
  integer   :: Z       ! charge number of target nucleus
  integer   :: Zix     ! charge number index for residual nucleus
  integer   :: Zparent ! Z of parent isotope
  real(sgl) :: acmax   ! maximum activity
  real(sgl) :: N0      ! number of isotopes
  real(sgl) :: rfac    ! conversion factor for radioactivity
  real(sgl) :: yfac    ! conversion factor for isotope yield
  real(dbl) :: C1      ! constant
  real(dbl) :: CP1     ! constant
  real(dbl) :: denom   ! help variable
  real(dbl) :: denomD  ! help variable
  real(dbl) :: dT      ! time step
  real(dbl) :: exp1    ! exponent
  real(dbl) :: exp2    ! exponent
  real(dbl) :: expo1   ! exponent
  real(dbl) :: expo2   ! exponent
  real(dbl) :: lamD    ! decay rate per isotope
  real(dbl) :: lamPD   ! decay rate for parent isotope
  real(dbl) :: prate0  ! production rate for all isotopes
  real(dbl) :: pratei  ! production rate per isotope
  real(dbl) :: prateP  ! production rate
  real(dbl) :: t1      ! help variable
  real(dbl) :: t2      ! help variable
  real(dbl) :: Tc      ! function for Coulomb dip
  real(dbl) :: term    ! help variable
  real(dbl) :: Th      ! time in hours
  real(dbl) :: TT      ! help variable
!
! ************ Initial condition for irradiation ***********************
!
  Ntar0 = avogadro / Atarget * rhotarget * Vtar
!
! ******************** Set time grid ***********************************
!
  Ntime = numtime / 2
  Th = Tir / hoursec
  dT = Th / Ntime
  if (dT <= 0.1) dT = 0.1
  if (dT > 0.1 .and. dT <= 0.2) dT = 0.2
  if (dT > 0.2 .and. dT <= 0.5) dT = 0.5
  if (dT > 0.5 .and. dT <= 1.0) dT = 1.0
  do it = 0, Ntime
    Tgrid(it) = it * dT
    if (abs(Tgrid(it) - Th) < 1./hoursec .or. Tgrid(it) > Th) then
      Tgrid(it) = Th
      Ntime = it
      exit
    endif
  enddo
  if (Tco > 0.) then
    Tc = Tco / hoursec
    Ncool = numtime - Ntime
    dT = log(Tc + 1.) / Ncool
    do it = Ntime + 1, numtime
      Tgrid(it) = Th + exp((it - Ntime) * dT) - 1.
    enddo
  endif
!
! ******************** Activity (yield) in MBq *************************
!
  Nisotot = 0.
  Tmax = 0.
  Tmaxactivity = 0
  Niso = 0.
  activity = 0.
  yield = 0.
  Nisorel = 0.
  Tp = 0
  Erp = 0.
  Nenrp = 0
  prate0 = dble(prate( - 1, - 1, - 1))
  do Zix = 0, maxZ
    Z = Zinit - Zix
    do Nix = 0, maxN
      N = Ninit - Nix
      A = Z + N
      do is = - 1, Nisomer(Zix, Nix)
        pratei = dble(prate(Zix, Nix, is))
        if (is >= 0 .and. pratei == 0.) cycle
        if (Z == Ztarget .and. A == Atarget .and. is ==  - 1) then
          Niso(Zix, Nix, is, 0) = dble(Ntar0)
          Nisotot(Zix, 0) = Nisotot(Zix, 0) + Niso(Zix, Nix, is, 0)
        endif
        lamD = dble(lambda(Zix, Nix, is))
        denomD = lamD - prate0
        C1 = dble(Ntar0) * pratei
        do it = 1, numtime
          TT = dble(Tgrid(it) * hoursec)
!
! Depletion of target
!
          if (Z == Ztarget .and. A == Atarget .and. is ==  - 1) then
            if (it <= Ntime) then
              Niso(Zix, Nix, is, it) = dble(Ntar0) * exp( - prate0 * TT)
            else
              Niso(Zix, Nix, is, it) = Niso(Zix, Nix, is, Ntime)
            endif
          else
!
! Production and decay of other isotopes
!
! 1. Production directly from target
!
            if (it <= Ntime) then
              t1 = prate0 * TT
              expo1 = exp( - t1)
              t2 = lamD * TT
              expo2 = exp( - t2)
              exp1 = expo1 / denomD - expo2 / denomD
              Niso(Zix, Nix, is, it) = C1 * exp1
            else
              t2 = lamD * (TT - Tir)
              expo2 = exp( - t2)
              term = Niso(Zix, Nix, is, Ntime) * expo2
              Niso(Zix, Nix, is, it) = term
            endif
!
! 2. Production from decay of other isotope
!
            do isob = - 1, 1
              Zparent = Zix - isob
              Nparent = Nix + isob
              if (Zparent < 0 .or. Nparent < 0) cycle
              if ((isob ==  - 1 .and. rtyp(Zparent, Nparent, - 1) == 1) .or. (isob == 1 .and. rtyp(Zparent, Nparent, - 1) == 2)) &
 &              then
                lamPD = dble(lambda(Zparent, Nparent, is))
                if (it <= Ntime) then
                  prateP = dble(prate(Zparent, Nparent, is))
                  if (prateP > 0..and.lamPD > 0.) then
                    CP1 = dble(Ntar0) * prateP
                    t1 = lamPD * TT
                    expo1 = exp( - t1)
                    denom = lamD - lamPD
                    exp2 = expo1 / denom - expo2 / denom
                    term = lamPD * CP1 / (lamPD - prate0) * (exp1 - exp2)
                    Niso(Zix, Nix, is, it) = Niso(Zix, Nix, is, it) + term
                  endif
                else
!
! Cooling only
!
                  N0 = Niso(Zparent, Nparent, is, Ntime)
                  t1 = lamPD * (TT - Tir)
                  exp1 = exp( - t1)
                  t2 = lamD * (TT - Tir)
                  exp2 = exp( - t2)
                  denom = lamD - lamPD
                  if (denom /= 0.) then
                    term = N0 * lamPD / denom * (exp1 - exp2)
                    Niso(Zix, Nix, is, it) = Niso(Zix, Nix, is, it) + term
                  endif
                endif
              endif
            enddo
            activity(Zix, Nix, is, it) = lamD * Niso(Zix, Nix, is, it) * 1.e-6
            if (it <= Ntime) yield(Zix, Nix, is, it) = &
              max((activity(Zix, Nix, is, it) - activity(Zix, Nix, is, it - 1)) / &
 &              (Ibeam * hoursec * dble(Tgrid(it) - Tgrid(it - 1))),0.)
          endif
          if (Niso(Zix, Nix, is, it) > 0.) Yexist(Zix, Nix, is) = .true.
          Nisotot(Zix, it) = Nisotot(Zix, it) + Niso(Zix, Nix, is, it)
        enddo
        if ( .not. Yexist(Zix, Nix, is)) cycle
        if (lamD > 0..and.prate0 > 0.) then
          Tmax(Zix, Nix, is) = log(lamD / prate0) / (lamD - prate0)
        else
          Tmax(Zix, Nix, is) = 1.e30
        endif
!
! Write irradiation time with maximum yield in years, days, etc.
!
        TT = Tmax(Zix, Nix, is)
        Tp(Zix, Nix, is, 1) = int(TT / yearsec)
        TT = TT - Tp(Zix, Nix, is, 1) * yearsec
        Tp(Zix, Nix, is, 2) = int(TT / daysec)
        TT = TT - Tp(Zix, Nix, is, 2) * daysec
        Tp(Zix, Nix, is, 3) = int(TT / hoursec)
        TT = TT - Tp(Zix, Nix, is, 3) * hoursec
        Tp(Zix, Nix, is, 4) = int(TT / minutesec)
        TT = TT - Tp(Zix, Nix, is, 4) * minutesec
        Tp(Zix, Nix, is, 5) = int(TT)
      enddo
    enddo
    do Nix = 0, maxN
      do is = - 1, Nisomer(Zix, Nix)
        if ( .not. Yexist(Zix, Nix, is)) cycle
        do it = 0, numtime
          if (Nisotot(Zix, it) /= 0.) Nisorel(Zix, Nix, is, it) = Niso(Zix, Nix, is, it) / Nisotot(Zix, it)
        enddo
      enddo
    enddo
  enddo
!
! Transform quantities to user-dependent units
!
  rfac = 1.
  yfac = 1.
  if (radiounit == 'bq') rfac = 1.e6
  if (radiounit == 'kbq') rfac = 1.e3
  if (radiounit == 'gbq') rfac = 1.e-3
  if (radiounit == 'ci') rfac = 1./3.7e4
  if (radiounit == 'kci') rfac = 1./3.7e7
  if (radiounit == 'mci') rfac = 1./3.7e1
  do Zix = 0, maxZ
    Z = Zinit - Zix
    do Nix = 0, maxN
      N = Ninit - Nix
      A = Z + N
      if (yieldunit == 'g') yfac = real(A)/avogadro
      if (yieldunit == 'mug') yfac = real(A)/avogadro*1.e6
      if (yieldunit == 'mg') yfac = real(A)/avogadro*1.e3
      if (yieldunit == 'kg') yfac = real(A)/avogadro*1.e-3
      do is = - 1, Nisomer(Zix, Nix)
        acmax = 0.
        do it = 1, numtime
          activity(Zix, Nix, is, it) = rfac * activity(Zix, Nix, is, it)
          yield(Zix, Nix, is, it) = rfac * yield(Zix, Nix, is, it)
          Niso(Zix, Nix, is, it) = yfac * Niso(Zix, Nix, is, it)
          if (activity(Zix, Nix, is, it) > acmax) then
            Tmaxactivity(Zix, Nix, is) = it
            acmax = activity(Zix, Nix, is, it)
          endif
        enddo
      enddo
    enddo
  enddo
  return
end subroutine prodyield
! Copyright A.J. Koning 2021
