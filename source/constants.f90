subroutine constants
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Constants, basic properties of particles and initialization
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
! Constants
!   amu         ! atomic mass unit in MeV
!   avogadro    ! Avogadro's number
!   clight      ! speed of light in vacuum in m / s
!   cparity     ! parity (character)
!   e2          ! square of elementary charge in MeV.fm
!   emass       ! electron mass in MeV / c^2
!   excmass     ! mass excess of particle in a.m.u.
!   hbar        ! Planck's constant / 2.pi in MeV.s
!   isochar     ! symbol of isomer
!   kT          ! energy kT expressed in MeV corresponding to a temperature T9 = 1
!   magic       ! magic numbers
!   nuc         ! symbol of nucleus
!   parA        ! mass number of particle
!   parmass     ! mass of particle in a.m.u.
!   parN        ! neutron number of particle
!   parname     ! name of particle
!   parspin     ! spin of particle
!   parsym      ! symbol of particle
!   parZ        ! charge number of particle
!   pi          ! pi
!   qelem       ! elementary charge in C
!   amu4pi2h2c2 ! amu / (4 * pi * pi * clight * clight * hbar **2) in mb ** - 1.MeV ** - 1
!   amupi2h3c2  ! amu / (pi * pi * clight * clight * hbar **3) in mb ** - 1.MeV ** - 2.s ** - 1
!   deg2rad     ! conversion factor for degrees to radians
!   Emaxtalys   ! maximum acceptable energy for TALYS
!   fislim      ! mass above which nuclide fissions
!   fourpi      ! 4. * pi
!   hbarc       ! hbar.c in MeV.fm
!   iso         ! counter for isotope
!   natstring   ! string extension for file names
!   onethird    ! 1 / 3
!   pardis      ! parity distribution
!   pi2         ! pi **2
!   pi2h2c2     ! 1 / (pi * pi * clight * clight * hbar **2) in mb ** - 1.MeV ** - 2
!   pi2h3c2     ! 1 / (pi * pi * clight * clight * hbar **3) in mb ** - 1.MeV ** - 3.s ** - 1
!   rad2deg     ! conversion factor for radians to degrees
!   sgn         ! sign
!   sqrttwopi   ! sqrt(2. * pi)
!   spin2       ! 2 * spin of particle
!   twopi       ! 2 * pi
!   twopihbar   ! 2 * pi / hbar
!   twothird    ! 2 / 3
! All global variables
!   numl           ! number of l values
! Variables for input energies
!   nin0           ! counter for incident energy
! Variables for reading input lines
!   inline    ! input line
!   nlines    ! number of input lines
! Variables for basic reactions
!   flagffruns ! flag to designate subsequent evap. of fission products
!   flagrpruns ! flag to designate that run is for residual product
!
! *** Declaration of local data
!
  implicit none
  integer :: i ! counter
!
! ****************** General properties of particles *******************
!
!
! Fundamental masses have been taken from physics.nist.gov 
!
! Indices: fission = -1
!   photon  = 0
!   neutron = 1
!   proton  = 2
!   deuteron= 3
!   triton  = 4
!   helium-3= 5
!   alpha   = 6
!
  parname = (/'fission ', 'gamma   ', 'neutron ', 'proton  ', 'deuteron', 'triton  ', 'helium-3', 'alpha   '/)
  parsym =  (/'f', 'g', 'n', 'p', 'd', 't', 'h', 'a'/)
  parZ =    (/ 0, 0, 1, 1, 1, 2, 2 /)
  parN =    (/ 0, 1, 0, 1, 2, 1, 2 /)
  parA =    (/ 0, 1, 1, 2, 3, 3, 4 /)
  parmass = (/ 0., 1.00866491595, 1.007276466621, 2.013553212745, 3.01550071621, 3.014932247175, 4.001506179127 /)
  excmass = (/ 0., 8.66491595e-3, 7.276466621e-3, 1.3553212745e-2, 1.550071621e-2, 1.4932247175e-2, 1.506179127e-3 /)
  parspin = (/ 0., 0.5, 0.5, 1., 0.5, 0.5, 0. /)
!
! ************************ Nuclear symbols *****************************
!
  nuc = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
    'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
    'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
    'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
    'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', &
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'B9', 'C0', 'C1', 'C2', 'C3', 'C4'/)
  magic = (/ 2, 8, 20, 28, 50, 82, 126, 184 /)
  isochar = (/' ', 'g', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v'/)
!
! **************** Set character symbol for parities *******************
!
! The '+' parity will have the value 1 and the '-' parity the value -1
!
  cparity = (/'-', ' ', '+'/)
!
! *********************** Fundamental constants ************************
!
  pi = 3.14159265358979323
  amu = 931.49386
  e2 = 1.439965161
  hbar = 6.5821220e-22
  clight = 2.99792458e8
  kT = 0.086173
  emass = 0.510999
  avogadro = 6.0221367e23
  qelem = 1.6021773e-19
!
! ************************* Derived constants **************************
!
  twopi = 2.*pi
  pi2 = pi * pi
  sqrttwopi = sqrt(2. * pi)
  fourpi = 4. * pi
  deg2rad = pi / 180.
  rad2deg = 180. / pi
  onethird = 1. / 3.
  twothird = 2. * onethird
  twopihbar = twopi / hbar
  hbarc = hbar * clight * 1.e15
  pi2h2c2 = 0.1 / (pi2 * hbarc * hbarc)
  pi2h3c2 = pi2h2c2 / hbar
  amupi2h3c2 = real(amu) * pi2h3c2
  amu4pi2h2c2 = real(amu) * 0.25 * pi2h2c2
!
! ************************** Set signs *********************************
!
! sgn is used in level density and compound nucleus calculations
!
  sgn(0) = 1.
  do i = 2, 2 * numl, 2
    sgn(i) = 1.
    sgn(i - 1) = - 1.
  enddo
!
! As long as we use an equidistant parity distribution, we set it equal to 0.5, instead of calling a subroutine.
! As soon as an analytical non-equidistant parity distribution is used in TALYS, pardis should be changed into
! a function pardis(J,Ex) throughout TALYS.
!
  pardis = 0.5
!
! ************************ Set counter for isotope *********************
!
  iso = 1
  natstring = '    '
  nlines = 0
  inline = ' '
!
! ************ Help variable for Hauser-Feshbach subroutine ************
!
! For alpha-particles and photons, spin2 is set to 1, to prevent division by zero in Hauser-Feshbach routines.
! This is fine, since in these cases the enumerator in the expression in which spin2 appears is always zero.
!
  spin2 = int(2. * parspin)
  spin2(0) = 1
  spin2(6) = 1
!
! *********************** Set default mass for fission *****************
!
  fislim = 215
!
! ***************** Set maximum acceptable energy for TALYS ************
!
  Emaxtalys = 1000.
!
! **************** Set counter for fission fragment evaporation ***********
!
  flagffruns= .false.
  flagrpruns= .false.
  nin0 = 0
  return
end subroutine constants
! Copyright A.J. Koning 2021
