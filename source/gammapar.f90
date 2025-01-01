subroutine gammapar(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Gamma ray parameters
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
! All global variables
!   numgamqrpa    ! number of energies for QRPA strength function
!   numTqrpa      ! number of temperatures for QRPA strength functions
! Variables for basic reaction
!   flagngfit       ! flag for using fitted (n,g) nuclear model parameters
! Variables for gamma rays
!   egr           ! energy of GR
!   epr           ! energy of PR
!   etable        ! constant to adjust tabulated strength functions
!   etableadjust  ! correction to adjust tabulated strength functions
!   flagpsfglobal ! flag for global photon strength functions only
!   Exlfile       ! tabulated gamma ray strength function
!   ftable        ! constant to adjust tabulated strength functions
!   ftableadjust  ! correction to adjust tabulated strength functions
!   gamadjust     ! logical for energy - dependent gamma adjustment
!   gammax        ! number of l - values for gamma multipolarity
!   ggr           ! width of GR
!   gpr           ! width of PR
!   sgr           ! strength of GR
!   strength      ! E1 strength function model
!   strengthM1    ! model for M1 gamma - ray strength function
!   tpr           ! strength of PR
!   wtable        ! constant to adjust tabulated strength functions
!   wtableadjust  ! correction to adjust tabulated strength functions
! Variables for level density
!   ldmodel         ! level density model
! Variables for main input
!   k0          ! index of incident particle
! Variables for masses
!   beta2         ! deformation parameter
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Constants
!   nuc           ! symbol of nucleus
! Variables for files
!   path          ! directory containing files to be read
! Constants
!   onethird      ! 1 / 3
!   pi            ! pi
!   pi2h2c2       ! 1 / (pi * pi * clight * clight * hbar **2) in mb ** - 1.MeV ** - 2
!  Variables for gamma-ray strength functions
!   eqrpa         ! energy grid for QRPA strength functio
!   fqrpa         ! tabulated QRPA strength functi
!   kgr           ! constant for gamma - ray strength functio
!   ngr           ! number of GR
!   nTqrpa        ! number of temperatures for QRPA
!   qrpaexist     ! flag for existence of tabulated QRPA st
!   Tqrpa         ! temperature for QRPA
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist           ! logical to determine existence
  character(len=6)  :: gamchar          ! help variable
  character(len=132):: key              ! keyword
  character(len=132):: gamfile          ! giant resonance parameter file
  integer           :: A                ! mass number of target nucleus
  integer           :: ia               ! mass number from abundance table
  integer           :: irad             ! variable to indicate M(=0) or E(=1) radiation
  integer           :: istat            ! logical for file access
  integer           :: it               ! counter for tritons
  integer           :: l                ! multipolarity
  integer           :: N                ! neutron number of residual nucleus
  integer           :: nen              ! energy counter
  integer           :: Nix              ! neutron number index for residual nucleus
  integer           :: Z                ! charge number of target nucleus
  integer           :: Zix              ! charge number index for residual nucleus
  real(sgl)         :: dE               ! delta energy
  real(sgl)         :: denom            ! help variable
  real(sgl)         :: dtemp            ! temperature increment
  real(sgl)         :: ee               ! energy
  real(sgl)         :: eg1              ! Gamow energies at T9(min) and T9(max)
  real(sgl)         :: eg2              ! Gamow energies at T9(min) and T9(max)
  real(sgl)         :: egamref          ! help variable
  real(sgl)         :: enum             ! enumerator of Lorentzian
  real(sgl)         :: Eq               ! help variable
  real(sgl)         :: Emid             ! middle energy of GDR
  real(sgl)         :: et               ! help variable
  real(sgl)         :: factor           ! multiplication factor
  real(sgl)         :: fe1(numTqrpa)    ! tabulated QRPA strength function
  real(sgl)         :: fe1t             ! tabulated strength function
  real(sgl)         :: fmax             ! maximum of GDR
  real(sgl)         :: fm1              ! tabulated M1 strength function
  real(sgl)         :: fstrength        ! gamma ray strength function
  real(sgl)         :: ft               ! help variable
  real(sgl)         :: wt               ! help variable
  real(sgl)         :: gg1              ! width of GR
  real(sgl)         :: gg2              ! width of GR
  real(sgl)         :: sg1              ! strength of GR
  real(sgl)         :: sg2              ! strength of GR
  real(sgl)         :: temp             ! nuclear temperature
!
! ***************** Default giant resonance parameters *****************
!
! We use indices to indicate the type of radiation, the multipolarity and the number of the peak (which sometimes is two).
! They are (0(M) or 1(E), l, number), i.e. egr(0,1,1) means a constant for M1-radiation and egr(1,2,1) a constant for E2-radiation.
!
! 1. Read experimental E1 values from GR parameter file
!
! GDR parameters from the table can always be overruled by a value given in the input file.
!
  nTqrpa = 1
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  N = A - Z
  if (.not.flagpsfglobal) then
    gamchar = trim(nuc(Z))//'.gdr'
    gamfile = trim(path)//'gamma/gdr/'//gamchar
    inquire (file = gamfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = gamfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(gamfile, istat)
      do
        read(2, '(4x, i4, 6f8.2)', iostat = istat) ia, eg1, sg1, gg1, eg2, sg2, gg2
        if (istat == -1) exit
        if (A == ia) then
          if (egr(Zix, Nix, 1, 1, 1) == 0.) egr(Zix, Nix, 1, 1, 1) = eg1
          if (sgr(Zix, Nix, 1, 1, 1) == 0.) sgr(Zix, Nix, 1, 1, 1) = sg1
          if (ggr(Zix, Nix, 1, 1, 1) == 0.) ggr(Zix, Nix, 1, 1, 1) = gg1
          if (egr(Zix, Nix, 1, 1, 2) == 0.) egr(Zix, Nix, 1, 1, 2) = eg2
          if (sgr(Zix, Nix, 1, 1, 2) == 0.) sgr(Zix, Nix, 1, 1, 2) = sg2
          if (ggr(Zix, Nix, 1, 1, 2) == 0.) ggr(Zix, Nix, 1, 1, 2) = gg2
          exit
        endif
      enddo
    endif
    close (unit = 2)
  endif
!
! 1. Default GR parameterization compiled by Kopecky in RIPL IAEA-TECDOC-1034, August 1998.
!
! E1-radiation
!
  if (egr(Zix, Nix, 1, 1, 1) == 0.) egr(Zix, Nix, 1, 1, 1) = 31.2 * A **( - onethird) + 20.6 * A **( - (1. / 6.))
  if (ggr(Zix, Nix, 1, 1, 1) == 0.) ggr(Zix, Nix, 1, 1, 1) = 0.026 * (egr(Zix, Nix, 1, 1, 1) **1.91)
  if (sgr(Zix, Nix, 1, 1, 1) == 0.) sgr(Zix, Nix, 1, 1, 1) = 1.2 * 120. * Z * N / (A * pi * ggr(Zix, Nix, 1, 1, 1))
!
! E2-radiation
!
  if (egr(Zix, Nix, 1, 2, 1) == 0.) egr(Zix, Nix, 1, 2, 1) = 63. * A **( - onethird)
  if (ggr(Zix, Nix, 1, 2, 1) == 0.) ggr(Zix, Nix, 1, 2, 1) = 6.11 - 0.012 * A
  if (sgr(Zix, Nix, 1, 2, 1) == 0.) sgr(Zix, Nix, 1, 2, 1) = 1.4e-4 * (Z **2) * egr(Zix, Nix, 1, 2, 1) / &
    (A **onethird * ggr(Zix, Nix, 1, 2, 1))
!
! E3-6 radiation
!
  do l = 3, gammax
    if (egr(Zix, Nix, 1, l, 1) == 0.) egr(Zix, Nix, 1, l, 1) = egr(Zix, Nix, 1, l - 1, 1)
    if (ggr(Zix, Nix, 1, l, 1) == 0.) ggr(Zix, Nix, 1, l, 1) = ggr(Zix, Nix, 1, l - 1, 1)
    if (sgr(Zix, Nix, 1, l, 1) == 0.) sgr(Zix, Nix, 1, l, 1) = sgr(Zix, Nix, 1, l - 1, 1) * 8.e-4
  enddo
!
! Check whether number of giant resonances is two (in which case all parameters must be specified)
!
  do irad = 0, 1
    do l = 1, gammax
      if (egr(Zix, Nix, irad, l, 2) /= 0..and. sgr(Zix, Nix, irad, l, 2) /= 0..and. egr(Zix, Nix, irad, l, 2) /= 0.) &
 &      ngr(Zix, Nix, irad, l) = 2
    enddo
  enddo
  if (strength > 2 .and. strength /= 5) then
!
! ***************** HFbcs or HFB QRPA strength functions ***************
!
! For Goriely's HFbcs or HFB QRPA strength function we overwrite the E1 strength function with tabulated results, if available.
!
    gamchar = trim(nuc(Z))//'.psf'
    if (strength == 3) gamfile = trim(path)//'gamma/hfbcs/'//gamchar
    if (strength == 4) gamfile = trim(path)//'gamma/hfb/'//gamchar
    if (strength == 6) gamfile = trim(path)//'gamma/hfbt/'//gamchar
    if (strength == 7) gamfile = trim(path)//'gamma/rmf/'//gamchar
    if (strength == 8) gamfile = trim(path)//'gamma/gogny/'//gamchar
    if (strength == 9) then
      if (flagpsfglobal) then
        gamfile = trim(path)//'gamma/smlo2019global/'//gamchar
      else
        gamfile = trim(path)//'gamma/smlo2019/'//gamchar
      endif
    endif
    if (strength == 10) gamfile = trim(path)//'gamma/bsk27_E1/'// gamchar
    if (strength == 11) gamfile = trim(path)//'gamma/d1m-intra-e1/'// gamchar
    if (strength == 12) gamfile = trim(path)//'gamma/shellmodel-e1/'// gamchar
    inquire (file = gamfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = gamfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(gamfile, istat)
      do
        read(2, '(10x, i4)', iostat = istat) ia
        if (istat <= -1) then
          call read_error(gamfile, istat, eor = 'continue', eof = 'continue')
          exit
        endif
        read(2, *, iostat = istat)
        if (istat /= 0) exit
        if (ia /= A) then
          do nen = 1, numgamqrpa + 1
            read(2, *, iostat = istat)
            if (istat /= 0) exit
          enddo
          cycle
        endif
        if (strength == 6 .or. strength == 7 .or. strength == 9 .or. strength == 10 .or. strength == 11) nTqrpa = 11
        do nen = 1, numgamqrpa
          read(2, '(f9.3, 20es12.3)', iostat = istat) ee, (fe1(it), it = 1, nTqrpa)
          if (istat /= 0) exit
          if (gamadjust(Zix, Nix)) then
            key = 'etable'
            call adjust(ee, key, Zix, Nix, 0, 0, factor)
            et = etable(Zix, Nix, 1, 1) + etableadjust(Zix, Nix, 1, 1) + factor - 1.
            key = 'ftable'
            call adjust(ee, key, Zix, Nix, 0, 0, factor)
            ft = ftable(Zix, Nix, 1, 1) * ftableadjust(Zix, Nix, 1, 1) + factor - 1.
          else
            et = etable(Zix, Nix, 1, 1) + etableadjust(Zix, Nix, 1, 1)
            ft = ftable(Zix, Nix, 1, 1) * ftableadjust(Zix, Nix, 1, 1)
          endif
          eqrpa(Zix, Nix, nen, 1, 1) = ee + et
          do it = 1, nTqrpa
            fqrpa(Zix, Nix, nen, it, 1, 1) = onethird * pi2h2c2 * fe1(it) * ft
          enddo
        enddo
!
! Avoid having an increasing function for extrapolation purposes
!
        do it = 1, nTqrpa
          if (fqrpa(Zix, Nix, numgamqrpa, it, 0, 1) > fqrpa(Zix, Nix, numgamqrpa - 1, it, 0, 1)) &
 &          fqrpa(Zix, Nix, numgamqrpa, it, 0, 1) = fqrpa(Zix, Nix, numgamqrpa - 1, it, 0, 1)
        enddo
        if (nTqrpa > 1) then
          if (strength == 11) then
            dtemp = 2.
          else
            dtemp = 0.2
          endif
          temp = - dtemp
          do it = 1, nTqrpa
            if (strength == 11) then
              if (it == 7) dtemp = 3.
              if (it >= 8) dtemp = 5.
            endif
            temp = temp + dtemp
            Tqrpa(it) = temp
          enddo
        endif
        qrpaexist(Zix, Nix, 1, 1) = .true.
        exit
      enddo
      close (unit = 2)
    endif
  endif
!
! ############################# M radiation ############################
!
! M1 radiation: strength function is related to that of E1
!
! strengthM1: model for M1 gamma-ray strength function
!
  if (strengthM1 <= 2 .or. strengthM1 == 4) then
    if (egr(Zix, Nix, 0, 1, 1) == 0.) egr(Zix, Nix, 0, 1, 1) = 41. * A **( - onethird)
    if (ggr(Zix, Nix, 0, 1, 1) == 0.) ggr(Zix, Nix, 0, 1, 1) = 4.
    if (sgr(Zix, Nix, 0, 1, 1) == 0.) then
      egamref = 7.
      enum = kgr(1) * egamref * ggr(Zix, Nix, 0, 1, 1) **2
      denom = (egamref **2 - egr(Zix, Nix, 0, 1, 1) **2) **2 + (ggr(Zix, Nix, 0, 1, 1) * egamref) **2
      if (strengthM1 == 1) then
        factor = 1.58e-9 * A **0.47
      else
        factor = fstrength(Zix, Nix, 0., egamref, 1, 1) / (0.0588 * A **0.878)
      endif
      sgr(Zix, Nix, 0, 1, 1) = factor * denom / enum
    endif
  endif
  if (strengthM1 == 3) then
!
! Spin-flip mode
!
    if (egr(Zix, Nix, 0, 1, 1) == 0.) egr(Zix, Nix, 0, 1, 1) = 18. / A **(1. / 6.)
    if (ggr(Zix, Nix, 0, 1, 1) == 0.) ggr(Zix, Nix, 0, 1, 1) = 4.
    if (sgr(Zix, Nix, 0, 1, 1) == 0.) sgr(Zix, Nix, 0, 1, 1) = 0.03 * A **(5. / 6.)
!
! Scissors mode
!
    if (epr(Zix, Nix, 0, 1, 1) == 0.) epr(Zix, Nix, 0, 1, 1) = 5. / A**(0.1)
    if (gpr(Zix, Nix, 0, 1, 1) == 0.) gpr(Zix, Nix, 0, 1, 1) = 1.5
    if (tpr(Zix, Nix, 0, 1, 1) == 0.) tpr(Zix, Nix, 0, 1, 1) =  1.0e-2 * abs(beta2(Zix, Nix, 0)) * A**(0.9)
  endif
!
! Includes Spin-flip from RIPL & Scissors mode from Kawano
!
  if (strengthM1 == 4) then
    if (epr(Zix, Nix, 0, 1, 1) == 0.) epr(Zix, Nix, 0, 1, 1) = 80. * abs(beta2(Zix, Nix, 0)) / A **(1. / 3.)
    if (gpr(Zix, Nix, 0, 1, 1) == 0.) gpr(Zix, Nix, 0, 1, 1) = 1.5
    if (tpr(Zix, Nix, 0, 1, 1) == 0.) tpr(Zix, Nix, 0, 1, 1) = 42.4 * beta2(Zix, Nix, 0) **2 /gpr(Zix, Nix, 0, 1, 1)
  endif
!
! Tabulated M1 strength
!
  if (strengthM1 == 8 .or. strengthM1 == 10 .or. strengthM1 == 11 .or. strengthM1 == 12) then
    gamchar = trim(nuc(Z))//'.psf'
    gamfile = trim(path)//'gamma/gognyM1/'//gamchar
    if (strengthM1 == 10) gamfile = trim(path)//'gamma/bsk27_M1/'//gamchar
    if (strengthM1 == 11 .and. Zix == 0 .and. Nix == 0) then
      gamfile = trim(path)//'gamma/d1m-intra-m1/'//gamchar
      nTqrpa=11
      dtemp = 2.
      temp = - dtemp
      do it = 1, nTqrpa
        if (it == 7) dtemp = 3.
        if (it >= 8) dtemp = 5.
        temp = temp + dtemp
        Tqrpa(it) = temp
      enddo
    endif
    if (strengthM1 == 12) gamfile = trim(path)//'gamma/shellmodel-m1/'//gamchar
    inquire (file = gamfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = gamfile, status = 'old')
      do
        read(2, '(10x, i4)', iostat = istat) ia
        if (istat == -1) exit
        read(2, * )
        if (ia == A) then
          do nen = 1, numgamqrpa
            read(2, '(f9.3, 20es12.3)') ee, fm1
            if (gamadjust(Zix, Nix)) then
              key = 'etable'
              call adjust(ee, key, Zix, Nix, 0, 0, factor)
              et = etable(Zix, Nix, 0, 1) + etableadjust(Zix, Nix, 0, 1) + factor - 1.
              key = 'ftable'
              call adjust(ee, key, Zix, Nix, 0, 0, factor)
              ft = ftable(Zix, Nix, 0, 1) * ftableadjust(Zix, Nix, 0, 1) + factor - 1.
            else
              et = etable(Zix, Nix, 0, 1) + etableadjust(Zix, Nix, 0, 1)
              ft = ftable(Zix, Nix, 0, 1) * ftableadjust(Zix, Nix, 0, 1)
            endif
            eqrpa(Zix, Nix, nen, 0, 1) = ee + et
            fqrpa(Zix, Nix, nen, 1, 0, 1) = onethird * pi2h2c2 * fm1 * ft
          enddo
!
! Avoid having an increasing function for extrapolation purposes
!
          if (fqrpa(Zix, Nix, numgamqrpa, 1, 0, 1) > fqrpa(Zix, Nix, numgamqrpa - 1, 1, 0, 1)) &
 &          fqrpa(Zix, Nix, numgamqrpa, 1, 0, 1) = fqrpa(Zix, Nix, numgamqrpa - 1, 1, 0, 1)
          qrpaexist(Zix, Nix, 0, 1) = .true.
          exit
        else
          do nen = 1, numgamqrpa + 1
            read(2, * )
          enddo
          cycle
        endif
      enddo
      close (unit = 2)
    endif
  endif
!
! Add some scissors mode contribution to spherical QRPA calculation
!
  if (strengthM1 == 10) then
    if (epr(Zix, Nix, 0, 1, 1) == 0.) epr(Zix, Nix, 0, 1, 1) = 5. / A**(0.1)
    if (gpr(Zix, Nix, 0, 1, 1) == 0.) gpr(Zix, Nix, 0, 1, 1) = 1.5
    if (tpr(Zix, Nix, 0, 1, 1) == 0.) tpr(Zix, Nix, 0, 1, 1) = 1.0e-2 * abs(beta2(Zix, Nix, 0)) * A**(0.9)
  endif
!
! M2-6 radiation
!
  do l = 2, gammax
    if (egr(Zix, Nix, 0, l, 1) == 0.) egr(Zix, Nix, 0, l, 1) = egr(Zix, Nix, 0, l - 1, 1)
    if (ggr(Zix, Nix, 0, l, 1) == 0.) ggr(Zix, Nix, 0, l, 1) = ggr(Zix, Nix, 0, l - 1, 1)
    if (sgr(Zix, Nix, 0, l, 1) == 0.) sgr(Zix, Nix, 0, l, 1) = sgr(Zix, Nix, 0, l - 1, 1) * 8.e-4
  enddo
!
! External strength functions
!
  do irad = 0, 1
    do l = 1, gammax
      if (Exlfile(Zix, Nix, irad, l)(1:1) /= ' ') then
        nen = 0
        open (unit = 2, file = Exlfile(Zix, Nix, irad, l), status = 'old')
        do
          read(2, * , iostat = istat) ee, fe1t
          if (istat == -1) exit
          if (istat /= 0) cycle
          if (gamadjust(Zix, Nix)) then
            key = 'etable'
            call adjust(ee, key, Zix, Nix, 0, 0, factor)
            et = etable(Zix, Nix, irad, l) + factor - 1.
            key = 'ftable'
            call adjust(ee, key, Zix, Nix, 0, 0, factor)
            ft = ftable(Zix, Nix, irad, l) + factor - 1.
          else
            et = etable(Zix, Nix, irad, l)
            ft = ftable(Zix, Nix, irad, l)
          endif
          nen = nen + 1
          eqrpa(Zix, Nix, nen, irad, l) = ee + et
          do it = 1, nTqrpa
            fqrpa(Zix, Nix, nen, it, irad, l) = onethird * pi2h2c2 * fe1t * ft
          enddo
          if (nen == numgamqrpa) exit
        enddo
        close (unit = 2)
        if (nen > 0) qrpaexist(Zix, Nix, irad, l) = .true.
        nTqrpa = 1
      endif
    enddo
  enddo
!
! Adjustment of width of tabulated PSF
!
  if (.not.flagpsfglobal) then
    do irad = 0, 1
      do l = 1, gammax
        wt = wtable(Zix, Nix, irad, l) * wtableadjust(Zix, Nix, irad, l)
        Emid = 0.
        fmax = 0.
        do nen = 1, numgamqrpa
          if (fqrpa(Zix, Nix, nen, 1, irad, l) > fmax) then
            fmax = fqrpa(Zix, Nix, nen, 1, irad, l)
            Emid = eqrpa(Zix, Nix ,nen, irad, l)
          endif
        enddo
        do nen = 1, numgamqrpa
          Eq = eqrpa(Zix, Nix, nen, irad, l)
          dE = Eq - Emid
          eqrpa(Zix, Nix, nen, irad, l) = Emid + dE * wt
        enddo
      enddo
    enddo
  endif
  return
end subroutine gammapar
! Copyright A.J. Koning 2021
