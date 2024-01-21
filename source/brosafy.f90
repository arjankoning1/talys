subroutine brosafy(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fission fragment yields based on Brosa model
!
! Author    : Marieke Duijvestijn
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
!   dbl           ! double precision kind
! All global variables
!   numelem       ! number of elements
!   nummass       ! number of masses
! Variables for fission
!   flagffevap    ! flag for calculation of particle evaporation from fissi
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Variables for files
!   path          ! directory containing files to be read
! Constants
!   amu           ! atomic mass unit in MeV
!   clight        ! speed of light in vacuum in m / s
!   hbar          ! Planck's constant / 2.pi in MeV.s
!   nuc           ! symbol of nucleus
!   parmass       ! mass of particle in a.m.u.
! Variables for mass distribution
!   disa          ! normalised fission fragment mass yield per excitation energy bin
!   disacor       ! normalised fission product mass yield per excitation energy bin
!   disaz         ! normalised fission fragment isotope yield per excitation energy
!   disazcor      ! normalised fission product isotope yield per excitation energy b
!   excfis        ! excitation energy at fission
! Variables for Brosa model
!   bf            ! barrier heigt
!   bfsplin       ! barrier height splin fit parameters
!   hw            ! barrier width
!   hwsplin       ! barrier width splin fit parameters
!   numtemp       ! running variable denoting the numtempth temperature
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist                          ! logical to determine existence
  character(len=6)  :: gschar                          ! structure database ground-state binding energy file
  character(len=8)  :: filen                           ! file name
  character(len=132):: barfile                         ! path of barrier parameter file in structure data base
  character(len=132):: gsfile                          ! full path for structure database file
  character(len=132):: precfile                        ! path with prescission shape info in structure data base
  integer           :: A                               ! mass number of target nucleus
  integer           :: amassar(1:7)                    ! mass array
  integer           :: amassdum                        ! dummy mass variable
  integer           :: amassmax                        ! heaviest isotope for which parameters have been calculated
  integer           :: i                               ! level
  integer           :: iloop                           ! loop counter
  integer           :: index1                          ! denote array elements used to interpolate
  integer           :: index2                          ! denote array elements used to interpolate
  integer           :: istat                           ! logical for file access
  integer           :: k                               ! designator for particle
  integer           :: massdif                         ! mass difference between amassmax and the nucleus studied
  integer           :: Nix                             ! neutron number index for residual nucleus
  integer           :: noff(8)                         ! number of neutrons away from heaviest isotope given  in structure data base
  integer           :: numoff                          ! number of isotopes for which information is present  in structure data base
  integer           :: numtempsl                       ! number of array indices per fission mode
  integer           :: numtempst                       ! number of array indices per fission mode
  integer           :: numtempst2                      ! number of array indices per fission mode
  integer           :: Z                               ! charge number of target nucleus
  integer           :: Zbrosa                          ! charge number within range of Brosa tables
  integer           :: Zix                             ! charge number index for residual nucleus
  real(dbl)         :: sl                              ! SL,ST I, ST II transmission coefficients
  real(dbl)         :: st                              ! denominator of compound nucleus formula
  real(dbl)         :: st2                             ! help variable
  real(dbl)         :: stot                            ! sum of sl, st, and st2
  real(dbl)         :: transm                          ! function for transmission coefficient per fission mode
  real(dbl)         :: trcof                           ! transmission coefficient
  real(sgl)         :: ald                             ! level density parameter
  real(sgl)         :: bar                             ! fission barrier height
  real(sgl)         :: bf_sl(9)                        ! SL outer fission barrier heights
  real(sgl)         :: bf_st(9)                        ! ST I outer fission barrier heights
  real(sgl)         :: bf_st2(9)                       ! ST II outer fission barrier heights
  real(sgl)         :: bfinter(1:9, 1:2)               ! barrier heights used for interpolation
  real(sgl)         :: bfsplin_sl(9)                   ! SL barrier height splin fit parameters
  real(sgl)         :: bfsplin_st(9)                   ! ST I barrier height splin fit parameters
  real(sgl)         :: bfsplin_st2(9)                  ! ST II barrier height splin fit parameters
  real(sgl)         :: BFT                             ! temperature-dependent barrier height
  real(sgl)         :: binddum                         ! dummy binding energy variable
  real(sgl)         :: bindgs(1:9, 1:2)                ! ground-state binding energy
  real(sgl)         :: bindsc                          ! Brosa prescission parameter
  real(sgl)         :: crel                            ! scaling factor for neck curvature
  real(sgl)         :: Edefo                           ! excitation energy at scission
  real(sgl)         :: elsc                            ! Brosa prescission parameter
  real(sgl)         :: ELT                             ! Brosa prescission parameter
  real(sgl)         :: ET                              ! sum of all partial energies
  real(sgl)         :: fraction                        ! fraction of potential energy gain available as internal  excitation energy
  real(sgl)         :: hm                              ! mean mass heavy fragment at scission point used  for interpolation
  real(sgl)         :: HMT                             ! help variable
  real(sgl)         :: hw_sl(9)                        ! width opf barrier
  real(sgl)         :: hw_st(9)                        ! width opf barrier
  real(sgl)         :: hwinter(1:9, 1:2)               ! barrier widths used for interpolation
  real(sgl)         :: ignatyuk                        ! function for energy dependent level density parameter
  real(sgl)         :: mn                              ! neutron mass in MeV
  real(sgl)         :: mp                              ! number of radial grid point
  real(sgl)         :: rmass                           ! help variable
  real(sgl)         :: somtot                          ! help variable
  real(sgl)         :: temps(9)                        ! temperature
  real(sgl)         :: Tmp                             ! temperature
  real(sgl)         :: Tmpmax                          ! maximum temperature
  real(sgl)         :: width                           ! Full width at half maximum of the breakup nucleon  energy distribution, Kal
  logical           :: dont                            ! logical for mass yield calculation
  real(sgl)         :: fmass(nummass)                  ! fission fragment mass yield
  real(sgl)         :: fmass_sl(nummass)               ! fission fragment mass yield for SL
  real(sgl)         :: fmass_st(nummass)               ! fission fragment mass yield for ST I
  real(sgl)         :: fmass_st2(nummass)              ! fission fragment mass yield for ST II
  real(sgl)         :: fmasscor_sl(nummass)            ! corrected fission fragment mass yield for SL
  real(sgl)         :: fmasscor_st(nummass)            ! corrected fission fragment mass yield for ST I
  real(sgl)         :: fmasscor_st2(nummass)           ! corrected fission fragment mass yield for ST II
  real(sgl)         :: fmasscor(nummass)               ! corrected fission fragment mass yield
  real(sgl)         :: fmz(nummass, numelem)           ! fission fragment isotope yield
  real(sgl)         :: fmzcor(nummass, numelem)        ! corrected fission fragment isotope yield
  real(sgl)         :: fmz_sl(nummass, numelem)        ! fission fragment isotope yield for SL
  real(sgl)         :: fmz_st(nummass, numelem)        ! fission fragment isotope yield for ST I
  real(sgl)         :: fmz_st2(nummass, numelem)       ! fission fragment isotope yield for ST II
  real(sgl)         :: fmzcor_sl(nummass, numelem)     ! corrected fission fragment isotope yield for SL
  real(sgl)         :: fmzcor_st(nummass, numelem)     ! corrected fission fragment isotope yield for ST
  real(sgl)         :: fmzcor_st2(nummass, numelem)    ! corrected fission fragment isotope yield for ST II
  real(sgl)         :: hmsplin(9)                      ! mean mass heavy fragment at scission point  splin fit parameters
  real(sgl)         :: elsplin(9)                      ! nucleus half length at scission point splin  fit parameters
  real(sgl)         :: Esplin(9)                       ! prescission energy splin fit parameters
  real(sgl)         :: hmneck_sl(9)                    ! mean mass heavy fragment at scission point for SL
  real(sgl)         :: hmneck_st(9)                    ! mean mass heavy fragment at scission point for ST
  real(sgl)         :: elneck_sl(9)                    ! nucleus half length at scission point for SL
  real(sgl)         :: elneck_st(9)                    ! nucleus half length at scission point for ST
  real(sgl)         :: Eneck_sl(9)                     ! prescission energy for SL
  real(sgl)         :: Eneck_st(9)                     ! prescission energy for ST
  real(sgl)         :: Eneck(9)                        ! prescission energy
  real(sgl)         :: hmneck(9)                       ! mean mass heavy fragment at scission point
  real(sgl)         :: elneck(9)                       ! nucleus half length at scission point
  real(sgl)         :: Eneckinter(9, 2)                ! prescission energy used for interpolation
  real(sgl)         :: hmneckinter(9, 2)               ! mean mass heavy fragment at scission point used  for interpolation
  real(sgl)         :: elneckinter(9, 2)               ! nucleus half length at scission point used for  interpolation
  temps =                  (/0.0, 0.3, 0.6, 0.9, 1.2, 1.6, 2.0, 2.5, 3.0/)
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
!
! Brosa parameters available between Z=72 and Z=96, outside this range
! we adopt the values of the boundary nuclides
!
  Zbrosa = max(72, min(Z, 96))
  mn = parmass(1) * amu
  mp = parmass(2) * amu
!
!   up to Np(Z=93) 7 isotopes per element calculated, above 6.
!
  noff(1) = 0
  noff(2) = 2
  noff(3) = 4
  noff(4) = 6
  noff(5) = 8
  if (Z < 94) then
    noff(6) = 11
    noff(7) = 14
    noff(8) = 0
    numoff = 7
  else
    noff(6) = 10
    noff(7) = 0
    noff(8) = 0
    numoff = 6
  endif
!
!   initialize
!
  do k = 1, 9
    bf_sl(k) = 0.
    bf_st(k) = 0.
    bf_st2(k) = 0.
    bfsplin_sl(k) = 0.
    bfsplin_st(k) = 0.
    bfsplin_st2(k) = 0.
  enddo
  do k = 1, 9
     do i = 1, 2
     bindgs(k, i) = 0.
    enddo
  enddo
!
!   reading ground state binding energies
!
  gschar = trim(nuc(Zbrosa))//'.fis'
  gsfile = trim(path)//'fission/brosa/groundstate/'//gschar
  open (unit = 2, file = gsfile, status = 'old')
!
!   determine position nucleus in array
!
  read(2, '(4x, i4)') amassmax
  rewind(2)
  do k = 1, numoff
    amassar(k) = amassmax - noff(k)
  enddo
  massdif = amassmax - A
  do k = 1, numoff
     if(massdif.GE.noff(k) .and. massdif.LE.noff(k + 1))then
        if(massdif.EQ.noff(k))then
           index1 = k
           index2 = - 1
        else
           if(massdif.EQ.noff(k + 1))then
              index1 = k + 1
              index2 = - 1
           else
              index1 = k
              index2 = k + 1
           endif
        endif
        cycle
     endif
     if(massdif.GT.noff(numoff))then
        index1 = numoff
        index2 = - 1
     else
        if(massdif.LT.noff(1))then
           index1 = 1
           index2 = - 1
        endif
     endif
  enddo
!
!   read in ground state binding energy as function of T
!
  i = 1
  k = 1
  if(index2.EQ. - 1)then
    do
      read(2, '(4x, i4, 15x, f15.5)', iostat = istat) amassdum, binddum
      if (istat == -1) exit
      if(amassdum.EQ.amassar(index1))then
        bindgs(i, 1) = abs(binddum)
        i = i + 1
      endif
    enddo
    close (unit = 2)
  else
    do
      read(2, '(4x, i4, 15x, f15.5)', iostat = istat) amassdum, binddum
      if (istat == -1) exit
      if(amassdum.EQ.amassar(index1))then
        bindgs(i, 1) = abs(binddum)
        i = i + 1
      else
        if(amassdum.EQ.amassar(index2))then
           bindgs(k, 2) = abs(binddum)
           k = k + 1
        endif
      endif
    enddo
    close (unit = 2)
  endif
!
!   superlong (sl), standard 1 (st), standard 2 (st2) loop
!
  do iloop = 1, 3
     do k = 1, 9
        bfsplin(k) = 0.
        hwsplin(k) = 0.
        bf(k) = 0.
        hw(k) = 0.
        do i = 1, 2
           bfinter(k, i) = 0.
           hwinter(k, i) = 0.
        enddo
      enddo
     if(iloop.EQ.1)then
        filen = trim(nuc(Zbrosa))//'.sl '
     else
        if(iloop.EQ.2)then
           filen = trim(nuc(Zbrosa))//'.st '
        else
           filen = trim(nuc(Zbrosa))//'.st2'
        endif
     endif
     barfile = trim(path)//'fission/brosa/barrier/'//filen
     inquire (file = barfile, exist = lexist)
     if (lexist)then
       open (unit = 10, file = barfile, status = 'old')
!
!   if fission modes exists:
!   read and calculate barrier parameters height, width, and position (in calculating the width, the mass at T=0 is used, since
!   the dependence on T is very small)
!
        i = 1
        k = 1
        do
          read(10, '(4x, i4, 15x, 2f15.5)', iostat = istat) amassdum, bar, width
          if (istat == -1) exit
          if(amassdum.EQ.amassar(index1))then
            rmass = (1 / (clight * clight)) * (Z * mp + (A - Z) * mn - bindgs(1, 1))
            bfinter(i, 1) = bar
            hwinter(i, 1) = hbar * sqrt(width * 1.D+30 / rmass)
            i = i + 1
          endif
          if(index2.NE. - 1)then
            if(amassdum.EQ.amassar(index2))then
              rmass = (1 / (clight * clight)) * (Z * mp + (A - Z) * mn - bindgs(1, 2))
              bfinter(k, 2) = bar
              hwinter(k, 2) = hbar * sqrt(width * 1.D+30 / rmass)
              k = k + 1
            endif
          endif
        enddo
        close (unit = 10)
     endif
!
!   interpolate, if necessary, barrier parameters
!
     numtemp = 1
     if(index2.EQ. - 1)then
        do k = 1, 9
           if(bfinter(k, 1).NE.0.)then
              bf(k) = bfinter(k, 1)
              hw(k) = hwinter(k, 1)
              numtemp = k
           endif
        enddo
     else
        if(abs(noff(index1) - noff(index2)).EQ.3)then
           do k = 1, 9
              if(bfinter(k, 1).NE.0 .and. bfinter(k, 2).NE.0.)then
                 massdif = abs(amassar(index1) - A)
                 bf(k) = bfinter(k, 1) + (massdif / 3.) * (bfinter(k, 2) - bfinter(k, 1))
                 hw(k) = hwinter(k, 1) + (massdif / 3.) * (hwinter(k, 2) - hwinter(k, 1))
                 numtemp = k
              endif
           enddo
        else
           do k = 1, 9
              if(bfinter(k, 1).NE.0 .and. bfinter(k, 2).NE.0.)then
                 bf(k) = 0.5 * (bfinter(k, 1) + bfinter(k, 2))
                 hw(k) = 0.5 * (hwinter(k, 1) + hwinter(k, 2))
                 numtemp = k
              endif
           enddo
        endif
     endif
!
!   fit spline to barrier parameters
!
     if(iloop.EQ.2)then
        bf(numtemp + 1) = bf_sl(numtemp + 1)
        hw(numtemp + 1) = hw_sl(numtemp + 1)
        numtemp = numtemp + 1
     endif
     if(iloop.EQ.3)then
        bf(numtemp + 1) = bf_st(numtemp + 1)
        hw(numtemp + 1) = hw_st(numtemp + 1)
        numtemp = numtemp + 1
     endif
     tmpmax = temps(numtemp)
     if(bf(1) == 0)numtemp = 9
     call spline(temps, bf, numtemp, 2.e+30, 2.e+30, bfsplin)
     call spline(temps, hw, numtemp, 2.e+30, 2.e+30, hwsplin)
!
! ignatyuk : function for energy dependent level density parameter
! trans    : subroutine to determine transmission coefficients per fission mode
!
!   calculate transmission coefficient trcof
!
     ald = ignatyuk(Zix, Nix, excfis, 0)
     Tmp = sqrt(excfis / ald)
!
!   check existence fission mode by looking at bf(1)
!
     if(Tmp < tmpmax .and. bf(1).NE.0.)then
        call trans(Zix, Nix, transm)
        trcof = transm
     else
        trcof = 0.d0
     endif
!
!   fill arrays with trcof, bf, bfsplin
!
     if(iloop.EQ.1)then
        sl = trcof
        numtempsl = numtemp
        do k = 1, numtemp
           bf_sl(k) = bf(k)
           hw_sl(k) = hw(k)
           bfsplin_sl(k) = bfsplin(k)
        enddo
     else
        if(iloop.EQ.2)then
           st = trcof
           numtempst = numtemp
           do k = 1, numtemp
              bf_st(k) = bf(k)
              hw_st(k) = hw(k)
              bfsplin_st(k) = bfsplin(k)
           enddo
        else
           st2 = trcof
           numtempst2 = numtemp
           do k = 1, numtemp
              bf_st2(k) = bf(k)
              bfsplin_st2(k) = bfsplin(k)
           enddo
        endif
     endif
!
!   end first loop over fission modes sl, st, st2
!
  enddo
!
!   final transmission coefficients
!
  stot = sl + st + st2
  if(stot.EQ.0.)then
     sl = 1.
     st = 0.
     st2 = 0.
  else
     sl = sl / stot
     st = st / stot
     st2 = st2 / stot
  endif
!
! neck       : subroutine to determine mass and isotope yield per fission mode
!
!   MASS DISTRIBUTION
!   second loop over fission modes starts here
!
  do iloop = 1, 3
     do k = 1, 9
        bf(k) = 0.
        bfsplin(k) = 0.
     enddo
     if(iloop.EQ.1)then
        dont = .false.
        numtemp = numtempsl
        if(sl == 0.)dont = .true.
        do k = 1, numtemp
           bf(k) = bf_sl(k)
           bfsplin(k) = bfsplin_sl(k)
        enddo
        crel = 0.1
        filen = trim(nuc(Zbrosa))//'.sl '
     else
        if(iloop.EQ.2)then
           dont = .false.
           numtemp = numtempst
           if(st == 0.)dont = .true.
           do k = 1, numtemp
              bf(k) = bf_st(k)
              bfsplin(k) = bfsplin_st(k)
           enddo
           crel = 0.3
           filen = trim(nuc(Zbrosa))//'.st '
        else
           dont = .false.
           numtemp = numtempst2
           if(st2 == 0.)dont = .true.
           do k = 1, numtemp
              bf(k) = bf_st2(k)
              bfsplin(k) = bfsplin_st2(k)
           enddo
           crel = 0.3
           filen = trim(nuc(Zbrosa))//'.st2'
        endif
     endif
!
!   initialize arrays for mass and charge yields
!
     do k = 1, nummass
        fmass(k) = 0.
        fmasscor(k) = 0.
        do i = 1, numelem
           fmz(k, i) = 0.
           fmzcor(k, i) = 0.
    enddo
  enddo
!
!   initialize arrays for neck parameters
!
     do k = 1, 9
        hmneck(k) = 0.
        elneck(k) = 0.
        Eneck(k) = 0.
        do i = 1, 2
           hmneckinter(k, i) = 0.
           elneckinter(k, i) = 0.
           Eneckinter(k, i) = 0.
    enddo
  enddo
     do k = 1, numtemp
        if(iloop.EQ.2 .and. k.EQ.numtemp)then
           hmneck(k) = hmneck_sl(k)
           elneck(k) = elneck_sl(k)
           Eneck(k) = Eneck_sl(k)
           cycle
        endif
        if(iloop.EQ.3 .and. k.EQ.numtemp)then
           hmneck(k) = hmneck_st(k)
           elneck(k) = elneck_st(k)
           Eneck(k) = Eneck_st(k)
           cycle
        endif
     enddo
     if(dont)goto 17998
!
!  open file with prescission output and read values heavy fragment mass
!  prescission energy, corresponding nucleus half length, and mass heavy fragment
!
   precfile = trim(path)//'fission/brosa/prescission/'//filen
   open (unit = 10, file = precfile, status = 'old')
     i = 1
     k = 1
     do
       read(10, '(4x, i4, 15x, 3f15.5)', iostat = istat) amassdum, hm, elsc, bindsc
       if (istat == -1) exit
       if(amassdum.EQ.amassar(index1))then
         hmneckinter(i, 1) = hm
         elneckinter(i, 1) = elsc
         Eneckinter(i, 1) = abs(bindsc) - bindgs(i, 1)
         i = i + 1
       endif
       if(index2.NE. - 1)then
         if(amassdum.EQ.amassar(index2))then
           hmneckinter(k, 2) = hm
           elneckinter(k, 2) = elsc
           Eneckinter(k, 2) = abs(bindsc) - bindgs(k, 2)
           k = k + 1
         endif
       endif
     enddo
     close (unit = 10)
!
!   interpolate, if necessary, prescission shape parameters
!
     if(index2.EQ. - 1)then
        do k = 1, 9
           if(hmneckinter(k, 1).NE.0.)then
              hmneck(k) = hmneckinter(k, 1) * A / amassar(index1)
              elneck(k) = elneckinter(k, 1)
              Eneck(k) = Eneckinter(k, 1)
           endif
        enddo
     else
        if(abs(noff(index1) - noff(index2)).EQ.3)then
           do k = 1, 9
              if(hmneckinter(k, 1).NE.0 .and. hmneckinter(k, 2).NE.0.)then
                 massdif = abs(amassar(index1) - A)
                 hmneck(k) = hmneckinter(k, 1) + (massdif / 3.) * (hmneckinter(k, 2) &
                      - hmneckinter(k, 1))
                 elneck(k) = elneckinter(k, 1) + (massdif / 3.) * (elneckinter(k, 2) &
                      - elneckinter(k, 1))
                 Eneck(k) = Eneckinter(k, 1) + (massdif / 3.) * (Eneckinter(k, 2) &
                      - Eneckinter(k, 1))
              endif
           enddo
        else
           do k = 1, 9
              if(hmneckinter(k, 1).NE.0 .and. hmneckinter(k, 2).NE.0.)then
                hmneck(k) = 0.5 * (hmneckinter(k, 1) + hmneckinter(k, 2))
                elneck(k) = 0.5 * (elneckinter(k, 1) + elneckinter(k, 2))
                Eneck(k) = 0.5 * (Eneckinter(k, 1) + Eneckinter(k, 2))
              endif
           enddo
        endif
     endif
!
!   spline fitting to neck input parameters
!
     call spline(temps, hmneck, numtemp, 2.e+30, 2.e+30, hmsplin)
     call spline(temps, elneck, numtemp, 2.e+30, 2.e+30, elsplin)
     call spline(temps, Eneck, numtemp, 2.e+30, 2.e+30, Esplin)
!
!   calculate excitation energy at scission (including an iteration
!   to calculate Edefo and ET consistently, one iteration suffices)
!
     ald = ignatyuk(Zix, Nix, excfis, 0)
     Tmp = sqrt(excfis / ald)
     call splint(temps, Eneck, Esplin, numtemp, Tmp, ET)
     call splint(temps, bf, bfsplin, numtemp, Tmp, BFT)
     fraction = 1.
     if(excfis.GT.BFT)then
        Edefo = fraction * (ET + BFT) + excfis - BFT
        if(ET.LE.0.)Edefo = excfis
        ald = ignatyuk(Zix, Nix, Edefo, 0)
        Tmp = sqrt(Edefo / ald)
        call splint(temps, Eneck, Esplin, numtemp, Tmp, ET)
        Edefo = fraction * (ET + BFT) + excfis - BFT
        if(ET.LE.0.)Edefo = excfis
     else
        Edefo = (ET + excfis) * fraction
        if(ET.LE.0.)Edefo = excfis
        ald = ignatyuk(Zix, Nix, Edefo, 0)
        Tmp = sqrt(Edefo / ald)
        call splint(temps, Eneck, Esplin, numtemp, Tmp, ET)
        Edefo = (ET + excfis) * fraction
        if(ET.LE.0.)Edefo = excfis
     endif
!
!    check if Tmp(Edefo) does not exceed the highest value for which the barrier is defined
!
     ald = ignatyuk(Zix, Nix, Edefo, 0)
     Tmp = sqrt(Edefo / ald)
     if(Tmp.GT.temps(numtemp))then
        ald = ignatyuk(Zix, Nix, excfis, 0)
        Tmp = sqrt(excfis / ald)
        if(Tmp.LT.temps(numtemp))then
           call splint(temps, hmneck, hmsplin, numtemp, Tmp, HMT)
           call splint(temps, elneck, elsplin, numtemp, Tmp, ELT)
           call neck(Z, A, fmass, fmasscor, fmz, fmzcor, HMT , Edefo, ELT, crel)
        endif
     else
        call splint(temps, hmneck, hmsplin, numtemp, Tmp, HMT)
        call splint(temps, elneck, elsplin, numtemp, Tmp, ELT)
        call neck(Z, A, fmass, fmasscor, fmz, fmzcor, HMT , Edefo, ELT, crel)
     endif
17998    continue
     if(iloop.EQ.1)then
        do k = 1, nummass
           fmass_sl(k) = fmass(k) * sl
           fmasscor_sl(k) = fmasscor(k) * sl
           do i = 1, numelem
              fmz_sl(k, i) = fmz(k, i) * sl
              fmzcor_sl(k, i) = fmzcor(k, i) * sl
    enddo
  enddo
        do k = 1, numtemp
           hmneck_sl(k) = hmneck(k)
           elneck_sl(k) = elneck(k)
           Eneck_sl(k) = Eneck(k)
        enddo
     else
        if(iloop.EQ.2)then
           do k = 1, nummass
              fmass_st(k) = fmass(k) * st
              fmasscor_st(k) = fmasscor(k) * st
           do i = 1, numelem
              fmz_st(k, i) = fmz(k, i) * st
              fmzcor_st(k, i) = fmzcor(k, i) * st
    enddo
  enddo
           do k = 1, numtemp
              hmneck_st(k) = hmneck(k)
              elneck_st(k) = elneck(k)
              Eneck_st(k) = Eneck(k)
           enddo
        else
           do k = 1, nummass
              fmass_st2(k) = fmass(k) * st2
              fmasscor_st2(k) = fmasscor(k) * st2
           do i = 1, numelem
              fmz_st2(k, i) = fmz(k, i) * st2
              fmzcor_st2(k, i) = fmzcor(k, i) * st2
    enddo
  enddo
        endif
     endif
!
!   end loop over sl, st and st2 modes
!
  enddo
!
!
  somtot = 1.
  do k = 1, nummass
     disa(k) = (fmass_sl(k) + fmass_st(k) + fmass_st2(k)) / somtot
     if (flagffevap) then
        disacor(k) = (fmasscor_sl(k) + fmasscor_st(k) + fmasscor_st2(k)) / somtot
     else
        disacor(k) = 0.
     endif
     do i = 1, numelem
     disaz(k, i) = (fmz_sl(k, i) + fmz_st(k, i) + fmz_st2(k, i)) / somtot
     if(flagffevap)then
        disazcor(k, i) = (fmzcor_sl(k, i) + fmzcor_st(k, i) + fmzcor_st2(k, i)) / somtot
     else
        disazcor(k, i) = 0.
     endif
    enddo
  enddo
  return
end subroutine brosafy
! Copyright A.J. Koning 2021
