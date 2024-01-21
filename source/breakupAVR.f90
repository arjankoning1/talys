subroutine breakupAVR
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deuteron breakup fractions, Avrigeanu model
!
! Author    : Marilena Avrigeanu
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! All global variables
!   numen        ! maximum number of outgoing energies
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0            ! index of incident particle
!   Ztarget      ! charge number of target nucleus
! Variables for energy grid
!   deltaE       ! energy bin around outgoing energies
!   ebegin       ! first energy point of energy grid
!   egrid        ! outgoing energy grid
!   Einc         ! incident energy in MeV
! Variables for nuclides
!   NN       ! neutron number of residual nucleus
!   ZZ       ! charge number of residual nucleus
! Variables for energies
!   eend         ! last energy point of energy grid
! Variables for incident channel
!   xsreacinc    ! reaction cross section for incident channel
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
! Constants
!   onethird     ! 1 / 3
!   nuc          ! symbol of nucleus
! Variables for preequilibrium
!   breakupmodel   ! model for break - up reaction: 1. Kalbach 2. Avrigeanu
! Variables for preequilibrium
!   breakupexist    ! logical for break up file
!   xsBF         ! nucleon inelastic breakup cross section
!   xsBFnuc      ! inelastic breakup enhancement brought by break-up
!   xsBUnuc      ! nucleon breakup cross section
!   xsEB         ! elastic breakup cross section
!   xspreeqbu    ! preequilibrium cross section per particle type and outgoing energy for break-up
!
! *** Declaration of local data
!
  implicit none
  integer   :: Acomp                ! CN mass number
  integer   :: nen                  ! energy counter
  integer   :: type                 ! particle type
  integer   :: N                    ! neutron number
  integer   :: Ncomp                ! CN neutron number
  integer   :: Z                    ! charge number
  integer   :: Zcomp                ! CN charge number
  real(sgl) :: arg                  ! help variable
  real(sgl) :: arg2                 ! squre of energy
  real(sgl) :: argcm                ! C.M. energy
  real(sgl) :: AT13                 ! A**1/3
  real(sgl) :: BCin                 ! effective incident/outgoing Coulomb barrier
  real(sgl) :: BCout(2)             ! effective incident/outgoing Coulomb barrier
  real(sgl) :: Bdeut                ! deuteon binding energy
  real(sgl) :: BFcont               ! help variable
  real(sgl) :: BFsum                ! help variable
  real(sgl) :: E0n03                ! break-up energy
  real(sgl) :: E0n03CM              ! C.M. break-up energy
  real(sgl) :: ebreakCM(2)          ! C.M. break-up energy
  real(sgl) :: ebreakLS(2)          ! LAB break-up energy
  real(sgl) :: EhighBU              ! help variable
  real(sgl) :: EnoBU                ! help variable
  real(sgl) :: EmaxnCM              ! maximum breakup nucleons energy, C.M. System
  real(sgl) :: EmaxnLS              ! maximum breakup nucleons energy, Lab. System
  real(sgl) :: En                   ! outgoing energy
  real(sgl) :: EnCM                 ! C.M. energy
  real(sgl) :: EnormEB              ! energy for elastic breakup fraction
  real(sgl) :: Eout                 ! outgoing energy
  real(sgl) :: fracBUBF             ! nucleon inelastic breakup fraction
  real(sgl) :: fracBUE              ! elastic breakup fraction
  real(sgl) :: fracBUp              ! nucleon breakup fraction
  real(sgl) :: fracBUT              ! total breakup fraction
  real(sgl) :: Gauss                ! function for Gaussian
  real(sgl) :: RnormEB              ! normalization factor for elastic breakup fraction
  real(sgl) :: Sigr                 ! reaction cross section for incident channel
  real(sgl) :: spec03CM(2, 0:numen) ! nucleon breakup spectrum in Center of Mass System
  real(sgl) :: spec03LS(2, 0:numen) ! nucleon breakup spectrum in Laboratory System
  real(sgl) :: speccorCM(2)         ! correction factors for breakup spectra
  real(sgl) :: speccorLS(2)         ! correction factors for breakup spectra
  real(sgl) :: sumtest(2)           ! check of breakup spectrum
  real(sgl) :: TnormBU              ! normalization factor for total breakup fraction
  real(sgl) :: w03                  ! width
  real(sgl) :: width                ! Full width at half maximum of the breakup nucleon  energy distribution, Kalbach 200
  real(sgl) :: xsaBU                ! breakup cross section
  real(sgl) :: xsBUT                ! TOTAL breakup cross section
!
! ************************** Avrigeanu model ***************************
!
  if ( .not. breakupexist) then
    breakupexist = .true.
    open (unit = 8, file = 'breakup.dat', status = 'unknown')
    write(8, * )"      "
    write(8, *)'***** TALYSREACTBU: breakupmodel k0 Einc call', ' npxsratios', breakupmodel, k0, Einc
    write(8, * ) "                                  "
    call npxsratios
  else
    open (unit = 8, file = 'breakup.dat', status = 'unknown', position = 'append')
  endif
  write(8, '(3x, "  d + ",i3, a2)') Atarget, nuc(Ztarget)
  write(8, *)'  DEUTERON break-up parameterization, M. Avrigeanu', ' and V. Avrigeanu'
  write(8, *)'  Phys. Rev. C95, 024607(2017), eqs 1-5, and', ' Refs. therein'
!
! Calculation of terms independent of emission energy.
!
! xsreacinc, Sigr : reaction cross section for incident channel
! Gauss           : function for Gaussian
!
!   xsBUnuc = xsEB + xsBF,  Eq. 4  PRC95, 024607.
!
! Avoiding to count xsEB twice, TOTAL breakup cross section is:
!   xsBUT = xsEB + 2*xsBF,  Eq. 5, PRC89,044613, or
!   xsBUT = 2*xsBUnuc -xsEB ,  Eq. 5, PRC95,024607.
!
  Sigr = xsreacinc
  Bdeut = 2.225
  BCin = Ztarget / 9.5
  BCout(1) = 0.
  BCout(2) = BCin
  AT13 = Atarget **onethird
!
! maximum breakup nucleon energy: Laboratory System
!
  EmaxnLS = Einc * (Atarget + 1.) / (Atarget + 2.) - Bdeut * (Atarget + 1.) / Atarget
!
! Center of Mass System
!
  EmaxnCM  = Einc * Atarget / (Atarget + 2.) - Bdeut
  if(EmaxnLS <= 0.) then
     EmaxnLS = 0.
     write(8, *)" deuteron energy lower than Bd, no breakup"
     write(* ,*)" deuteron energy lower than Bd, no breakup"
     write(8, *)" incident deuteron energy too low"
     goto 699
  endif
  write(8, * )"    "
  write(8, * )" deuteron energy & maximum outgoing fragments energy"
  write(8, * )"    Einc   EmaxnLS   EmaxnCM  "
  write(8, 889)Einc, EmaxnLS, EmaxnCM
  RnormEB = 1.
  TnormBU = 1.
  xsBUT = 0.
  xsaBU = 0.
  do type = 1, 2
    ebreakLS(type) = 0.
    ebreakCM(type) = 0.
    sumtest(type) = 0.
    speccorCM(type) = 0.
    speccorLS(type) = 0.
    do nen = 0, eend(type)
      spec03LS(type, nen) = 0.
      spec03CM(type, nen) = 0.
      xspreeqbu(type, nen) = 0.
    enddo
  enddo
  call checkBU(EnormEB, RnormEB, TnormBU)
  write(8, * )"Ei_EB_norm   EB_norm ", "  total_BU_norm "
  write(8, 889)EnormEB, RnormEB, TnormBU
  arg = Einc
  arg2 = arg * arg
  argcm = arg * Atarget / (Atarget + 2)
!
!   ********   breakup fractions        ********
!
!   NUCLEON TOTAL breakup fraction, Eq. 1, PRC95,024607.
!
  fracBUp = 0.087 - 0.0066 * Ztarget + 0.00163 * Ztarget * AT13 + 0.0017 * AT13 * arg - 0.000002 * Ztarget * arg2
  fracBUp = TnormBU * fracBUp
  if (fracBUp <= 0.) then
     fracBUp = 0.
     write(8, * )" nucleon total breakup fraction zero"
     goto 699
  endif
  fracBUE = 0.031 - 0.0028 * Ztarget + 0.00051 * Ztarget * AT13 + 0.0005 * AT13 * arg - 0.000001 * Ztarget * arg2
  if(TnormBU /= 1 .or. RnormEB /= 1.) then
    if(Einc < EnormEB) then
      fracBUE = TnormBU * fracBUE
    else
      fracBUE = RnormEB * fracBUp
    endif
  endif
!
!   nucleon INELASTIC breakup fraction, PRC88,014612, Eq. 3
  fracBUBF = fracBUp - fracBUE
!
!  fracBUT: TOTAL breakup fraction, Eq. 5: PRC95,024607 & PRC89,044613
!
  fracBUT = 2.D+00 * fracBUBF + fracBUE
!
!   ********     end    breakup fractions       ********
!
!   ******************************
!    TOTAL BREAKUP cross section
  xsBUT = Sigr * fracBUT
!    crossafterbreakup  xsaBU
  xsaBU = Sigr - xsBUT
!   ******************************
!
!-------------           deuteron fragments   loop         --------
!
  do type = 1, 2
!
!   ******************************
!
!    ELASTIC breakup cross section
    xsEB(type) = Sigr * fracBUE
!    INELASTIC breakup cross section (BF) equal for neutron and proton
    xsBF(type) = Sigr * fracBUBF
!    TOTAL NUCLEON breakup cross section = xsEB+xsBF
    xsBUnuc(type) = Sigr * fracBUp
!
!   ******************************
!
    if(type == 1) then
      write(8, * )"     BU fractions (equal for n and p)        ", "    BU cross sections"
      write(8, * )"    Einc    fracBUp   fracBUE", "   fracBUBF      xsBUnuc    xsEB      xsBF"
      write(8, 889)Einc, fracBUp, fracBUE, fracBUBF, xsBUnuc(1), xsEB(1), xsBF(1)
      write(8, * )" xsBUnuc = xsEB + xsBF=", xsBUnuc(Type), ", Eq. 4,  PRC95, 024607 (2017). "
      write(8, * )" xsBUTOTAL = xsEB + 2*xsBF=", xsBUT, ", Eq. 5, PRC95, 024607 & PRC89, 044613."
    endif
!
!  ebreakLS, ebreakCM  : Centroid energy of the breakup nucleon energy
!   distributions, Kalbach 2003 in Laboratory, and
!   respectively Center of Mass Systems
!   nucleon energy distribution, Kalbach 2003
!
  ebreakCM(type) = (0.5 * (argcm - Bdeut - BCin) + BCout(type))
  ebreakLS(type) = 0.5 * (Atarget + 1.) / (Atarget + 2.) * Einc + &
     0.5 * (Atarget + 1.) / Atarget * ( - Bdeut - BCin + 2. * BCout(type))
  width = 1.15 + 0.12 * Einc - Atarget / 140.
!
  if (width <= 0.) then
    width = 0.
  else
!
    if (ebreakLS(type) <= 0.) ebreakLS(type) = 0.01
    if (ebreakCM(type) <= 0.) ebreakCM(type) = 0.01
!
    E0n03 = ebreakLS(type)
    E0n03CM = ebreakCM(type)
    w03 = width
!
    speccorCM(type) = 0.
    speccorLS(type) = 0.
!
!------------------    breakup nucleon energy   loop   -----------------
!
    do nen = ebegin(type), eend(type)
!
      Eout = egrid(nen) !ine breakup
      En = Eout
      EnCM = Eout
!
!   breakup nucleon spectra in Laboratory System
!
      if (Eout <= EmaxnLS) then
        spec03LS(type, nen) = GAUSS(En, E0n03, w03)
!
!   BREAKUP THRESHOLD:  En>Emax BU
!
      else
        spec03LS(type, nen) = 0.D+00
      endif
!
!   breakup nucleon spectra in Center of Mass
!
      if  (Eout <= EmaxnCM) then
        spec03CM(type, nen) = GAUSS(EnCM, E0n03CM, w03)
!
!   BREAKUP THRESHOLD: EnCM>Emax BU
!
      else
        spec03CM(type, nen) = 0.D+00
      endif
!
      speccorCM(type) = speccorCM(type) + spec03CM(type, nen) * deltaE(nen)
      speccorLS(type) = speccorLS(type) + spec03LS(type, nen) * deltaE(nen)
!
    enddo
!
!-------------------  end   breakup nucleon energy  loop ---------------
!
  endif
!
  enddo
!
!----------------------     end   deuteron fragments      --------------
!
  do type = 1, 2
    if (speccorCM(type) == 0..or.speccorLS(type) == 0.) cycle
    sumtest(type) = 0.D+00
    do nen = ebegin(type), eend(type)
      Eout = egrid(nen)
      spec03LS(type, nen) = xsBUnuc(type) * spec03LS(type, nen) / speccorLS(type)
      spec03CM(type, nen ) = xsBUnuc(type) * spec03CM(type, nen) / speccorCM(type)
        sumtest(type) = sumtest(type) + spec03CM(type, nen) * deltaE(nen)
!
!   ******************************
  xspreeqbu(type, nen) = spec03CM(type, nen)
!   ******************************
    enddo
  enddo
!
! Inelastic Break-up enhancement of Avrigeanu model
! Included smooth disappearance of enhancement above 80 MeV
!
699 do Acomp = 1, 9
    do Zcomp = 0, 5
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxN) cycle
      Z = ZZ(Zcomp, Ncomp, 0)
      N = NN(Zcomp, Ncomp, 0)
      BFsum = 0.
      EhighBU = 80.
      EnoBU = 120.
      if (Einc <= EhighBU) then
        call buenhance(Z, N, BFsum, Einc)
        BFcont = BFsum
      else
        call buenhance(Z, N, BFsum, EhighBU)
        BFcont = BFsum * (min(Einc, EnoBU) - EhighBU) / (EhighBU - EnoBU)
      endif
      xsBFnuc(Zcomp, Ncomp) = BFcont
    enddo
  enddo
  write(8, *) '  end ..ine breakupAVR'
  close (unit = 8)
!
889   FORMAT(2x,f8.3,2x,f8.4,2x,f8.4,2x,f8.5,5x,f8.3,2x,F8.3,2x, F8.3,2x,F8.3,2x,F8.3)
!
  return
end subroutine breakupAVR
! Copyright A.J. Koning 2021
subroutine checkBU(EnormEB, RnormEB, TnormBU)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : check the elastic and nucleon breakup fraction
!
! Author    : Marilena Avrigeanu
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
! Definition of single and double precision variables
!   sgl         ! single precision kind
! Variables for main input
!   Atarget     ! mass number of target nucleus
!   Ztarget     ! charge number of target nucleus
! Constants
!   onethird    ! 1 / 3
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: IRDIM=6001       ! dimension fo breakup subroutine
  integer            :: I                ! counter
  integer            :: Inorm            ! help variable
  real(sgl)          :: arg              ! help variable
  real(sgl)          :: arg2             ! squre of energy
  real(sgl)          :: AT13             ! A**1/3
  real(sgl)          :: bind             ! energy bin
  real(sgl)          :: Edi(IRDIM)       ! energy
  real(sgl)          :: EnormEB          ! energy for elastic breakup fraction
  real(sgl)          :: fBU(IRDIM)       ! total breakup fraction
  real(sgl)          :: fBUBF(IRDIM)     ! nucleon inelastic breakup fraction
  real(sgl)          :: fBUE(0:IRDIM)    ! elastic breakup fraction
  real(sgl)          :: fBUT(IRDIM)      ! nucleon breakup fraction
  real(sgl)          :: FEDrat(IRDIM)    ! ratio for breakup calculation
  real(sgl)          :: RnormEB          ! normalization factor for elastic breakup fraction
  real(sgl)          :: Tnorm(0:IRDIM)   ! normalization factor
  real(sgl)          :: TnormBU          ! normalization factor for total breakup fraction
!
!   ***  Deuteron breakup fractions, Eqs.(1-3), PRC 88,014612(2013)  ***
!
! IRDIM              : dimension fo breakup subroutine
!
!  According to the CDCC predictions for the elastic breakup cross
!  section behavior(PRC 82, 037601 (2010)), for higher energy than the
!  energetic domain (~30 MeV) where the parametrization was obtained,
!  the elastic breakup fraction is mentained constant, through the
!  RnormEB factor,    fBUE=RnormEB*fBUT.
!
!  Mainly for heavy nuclei, A~200, and  for Einc~Coulomb barrier, the
!  normalization constant TnormBu prevents the total breakup cross
!  section to exceed the reaction cross section.
!
  AT13 = Atarget**onethird
  RnormEB = 1.
  TnormBU = 1.
  bind = 0.1
  Inorm = int (80. / 0.1 + 1.1)
!
  do I = 1, IRDIM
     Edi(I) = 0.
     FEDrat(I) = 0.
     Tnorm(I) = 0.
     fBUT(I) = 0.
     fBUE(I) = 0.
     fBUBF(I) = 0.
     fBU(I) = 0.
  enddo
  Tnorm(0) = 0.
  fBUE(0) = 0.
!
!----------------------     deuteron incident energy   -----------------
!
  do i = 1, Inorm
!
  Edi(I) = i * bind
  arg = Edi(I)
  arg2 = arg * arg
!
!   ***    PRC 95,024607 (2017)    parametrization        ***
!
    fBUT(I) = 0.087 - 0.0066 * Ztarget + 0.00163 * Ztarget * AT13 + 0.0017 * AT13 * arg - 0.000002 * Ztarget * arg2
     if(fBUT(I) <= 0.d+00) then
        fBUT(I) = 0.
     else
!
       fBUE(I) = 0.031 - 0.0028 * Ztarget + 0.00051 * Ztarget * AT13 + 0.0005 * AT13 * arg - 0.000001 * Ztarget * arg2
!
       FEDrat(I) = fBUE(I) / fBUT(I)
!
!   elastic breakup normalization for Ed>25 MeV
!
       if(fBUE(I) < fBUE(I - 1)) then
         if((fBUE(I) / fBUE(I - 1)) >= 0.9999) then
           INORM = I - 1
           EnormEB = INORM * bind
           RnormEB = FEDrat(INORM)
         endif
!
         fBUE(I) = RnormEB * fBUT(I)
!
!   end elastic breakup normalization for Ed>25 MeV
!
      endif
!
      fBUBF(i) = fBUT(i) - fBUE(i)
      fBU(i) = 2 * fBUT(i) - fBUE(i)
!
!    total normalization for very heavy nuclei
!
      if(fBU(I) > 0.9D+00) then
        Tnorm(I) = 0.9D+00 / fBU(I)
!
        if (Tnorm(I) < Tnorm(I - 1)) TnormBU = Tnorm(I)
!
!   end total normalization for very heavy nuclei
!
      endif
    endif
  enddo
!
!   ***    end Phys. Rev. C 95, 024607 (2017)  parametrization    ***
!
  return
end subroutine checkBU
! Copyright A.J. Koning 2021
function GAUSS(En, E0n, w0)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Marilena Avrigeanu
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl    ! single precision kind
!
! *** Declaration of local data
!
  implicit none

  real(sgl) :: arg              !
  real(sgl) :: E0n              !
  real(sgl) :: En               !
  real(sgl) :: Gauss            !
  real(sgl) :: term1            !
  real(sgl) :: PI               !
  real(sgl) :: w0               !
  real(sgl) :: XG               !
  PI = ACOS(-1.)
  term1 = 1. / (w0 * sqrt(2.D+00 * pi))
  arg = - (En - E0n) * (En - E0n) / (2.D+00 * w0 * w0)
  XG = term1 * exp( - arg)
  GAUSS = term1 * exp(arg)
  return
end function GAUSS
! Copyright A.J. Koning 2021
