subroutine npxsratios
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Reads from 'structure' breakup n & p c.s. ratios requested
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
!   dbl         ! double precision kind
! Constants
!   nuc          ! symbol of nucleus
! All global variables
!   numenout     ! number of outgoing energies
!   numN         ! maximum number of neutrons from initial compound nucleus
!   numZ         ! maximum number of protons from initial compound nucleus
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   Ninit           ! neutron number of initial compound nucleus
!   Zinit           ! charge number of initial compound nucleus
!   Ztarget      ! charge number of target nucleus
! Variables for preequilibrium
!   ebubin      ! outgoing breakup nucleon energy bin for integration
!   ENHratio    ! breakup nucleons enhancing reaction cross
! Variables for nuclides
!   AA       ! mass number of residual nucleus
!   ZZ       ! charge number of residual nucleus
! Variables for files
!   path            ! directory containing files to be read
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist          ! logical to determine existence
  character(len=3)  :: Zstring         !
  character(len=132):: enhBUnxs        !
  character(len=132):: enhBUpxs        !
  integer           :: A               ! mass number of target nucleus
  integer           :: Acomp           ! mass number index for compound nucleus
  integer           :: ien             !
  integer           :: iep             !
  integer           :: Ioutn           !
  integer           :: Ioutp           !
  integer           :: N               ! neutron number of residual nucleus
  integer           :: Ncomp           ! neutron number index for compound nucleus
  integer           :: Nix             ! neutron number index for residual nucleus
  integer           :: type            ! particle type
  integer           :: Z               ! charge number of target nucleus
  integer           :: Zcomp           ! proton number index for compound nucleus
  integer           :: Zix             ! charge number index for residual nucleus
  real(dbl)         :: ENHN            !
  real(dbl)         :: ENHP            !
  real(dbl)         :: enout           !
  real(dbl)         :: Eoutn           !
  real(dbl)         :: Eoutp           !
  real(dbl)         :: epout           !
  real(dbl)         :: NUCLEUN         !
  real(dbl)         :: NUCLEUNr        !
  real(dbl)         :: NUCLEUP         !
  real(dbl)         :: NUCLEUPr        !
!
! ***  n,p reaction cross sections ratios for breakup enhancement  ***
!
!   Eq. (2), PRC 94,014606 (2016)
!   n + Atarget:  sig(n,Z,A,Enout)/sig_Total(n,Enout)
!   p + Atarget:  sig(p,Z,A,Epout)/sig_Reaction(p,Epout)
!   included in STRUCTURE directory as follows:
!   structure\breakup\neutrons\ENHratioN.da
!   structure\breakup\protons\ENHratioP.da
!
  write(8, *)' '
  write(8, *) '  start npxsratios  for   d + ', Atarget, nuc(Ztarget)
  write(8, *)' '
  write(8, *)' numZ numN numA numenout', numZ, NumN, numZ+numN, numenout
!
  ENHratio = 0.
!
!   Inquire whether breakup enhancing ratios file are presen
!
   Zstring = '000'
   write(Zstring, '(i3.3)') Ztarget
   do type = 1, 2
    if (type == 1) then
!   enhBUnxs=trim(path)//'breakup/neutrons/ENHratioN.dat'
    enhBUnxs = trim(path)//'breakup/neutrons/ENHratioN.'//Zstring
    inquire (file = enhBUnxs, exist = lexist)
    if ( .not. lexist) go to 409
    open (unit = 25, status = 'unknown', file = enhBUnxs)
!   enhBUpxs=trim(path)//'breakup/protons/ENHratioP.dat'
    enhBUpxs = trim(path)//'breakup/protons/ENHratioP.'//Zstring
    inquire (file = enhBUpxs, exist = lexist)
    if ( .not. lexist) go to 409
    open (unit = 26, status = 'unknown', file = enhBUpxs)
    endif
  enddo
!
  ebubin = 0.1
!
! ************     Inelastic breakup enhancing ratios     ************
!
!   in Eq.(2), PRC 94,014606
!   brought by type=1 (breakup neutrons) are those in
!   the area:  Zcomp=1:5
!   brought by type=2 (breakup protons) are those in
!   the area:  Zcomp=0:4
!
! ************   Read breakup neutrons enhancing ratios   ************
!
  write(8, *) 'residual nuclei for BU enhancement calculation'
  do 33 type = 1, 2
   if(type == 1)then
     do 15 Zcomp = 1, 5
     do 13 Acomp = 1, 9
       Ncomp = Acomp - Zcomp
       if (Ncomp < 0 .or. Ncomp > 4) goto 11
       Z = ZZ(Zcomp, Ncomp, 0)
       A = AA(Zcomp, Ncomp, 0)
       N = AA(Zcomp, Ncomp, 0) - ZZ(Zcomp, Ncomp, 0)
       NUCLEUN = Ztarget * 1.0D+9 + Atarget * 1.0D+6 + Z * 1.0D+3 + A
5            read(25, * , end = 409)NUCLEUNr, Eoutn
       Ioutn = int(Eoutn)
       if(NUCLEUNr.EQ.NUCLEUN)then
         write(8, *)' NUCLEUN=', NUCLEUN, ' Z=', Z, ' A=', A, ' N=', N, ' * NUCLEUNr=', NUCLEUNr
         Zix = Zinit - Z
         Nix = Ninit - N
         do ien = 1, Ioutn
           enout = ebubin * ien
           read(25, * )enout, ENHN
           ENHratio(type, Zix, Nix, ien) = ENHN
         enddo
       else
         do ien = 1, Ioutn + 1
           read(25, * )
         enddo
         go to 5
       endif
11           continue
13         rewind (25)
     continue
15  continue
!
   else
     write(8, *)'               ***'
!
! ************   Read breakup protons enhancing ratios   ************
!
     do 29 Zcomp = 0, 4
       do 27 Acomp = 1, 9
         Ncomp = Acomp - Zcomp
         if (Ncomp <= 0 .or. Ncomp > 5) goto 25
         Z = ZZ(Zcomp, Ncomp, 0)
         A = AA(Zcomp, Ncomp, 0)
         N = AA(Zcomp, Ncomp, 0) - ZZ(Zcomp, Ncomp, 0)
         NUCLEUP = Ztarget * 1.0D+9 + Atarget * 1.0D+6 + Z * 1.0D+3 + A
21            read(26, * , end = 409)NUCLEUPr, Eoutp
         Ioutp = int(Eoutp)
         if (NUCLEUPr == NUCLEUP) then
           write(8, *)' NUCLEUP = ', NUCLEUP, ' Z = ', Z, ' A = ', A, ' N = ', N, ' * NUCLEUPr = ', NUCLEUPr
           Zix = Zinit - Z
           Nix = Ninit - N
           do iep = 1, Ioutp
             epout = ebubin * Iep
             read(26, * )epout, ENHP
             ENHratio(type, Zix, Nix, Iep) = ENHP
           enddo
         else
           do ien = 1, Ioutp + 1
             read(26, * )
           enddo
          go to 21
         endif
25       continue
27       rewind(26)
     continue
29 continue
   endif
33 continue
!
      write(8, *)'  '
      write(8, *)' *** Test npxsratios *** '
!
  do type = 1, 2
    do Zcomp = 1, 5
      do Acomp = 1, 9
        Ncomp = Acomp - Zcomp
        if (Ncomp < 0 .or. Ncomp > 4) goto 313
        Z = ZZ(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        N = AA(Zcomp, Ncomp, 0) - ZZ(Zcomp, Ncomp, 0)
        Zix = Zinit - Z
        Nix = Ninit - N
        write(8, *)' type Z A N ENHratio(type,Z,N,320)', type, Z, A, N, ENHratio(type, Zix, Nix, 320)
313     continue
      enddo
    enddo
  enddo
!
  write(8, *)' *** END test npxsratios *** '
  write(8, *)'  '
!
  go to 601
!
409     write(8, *)' '
  write(8, *)'subroutine npxsratios: missing breakup library in TALYS structure or missing xs files from TENDL-2019 '
  write(8, *)'***     NO BUenhancement taken into account     ***'
  write(8, *)'  '
!
  write(*, *)' '
  write(*, *)'subroutine npxsratios: missing breakup library in TALYS structure or missing xs files from TENDL-2019 '
  write(*, *)'***     NO BUenhancement taken into account     ***'
  write(*, *)' '
!
  ENHratio = 0.
!
601     continue
  return
end subroutine npxsratios
! Copyright A.J. Koning 2021
