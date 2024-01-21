subroutine incidentgrid(fname, fexist)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Predefined incident energy grids
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
!   sgl          ! single precision kind
! All global variables
!   numenin      ! number of incident energies
! Variables for input energies
!   eninc        ! incident energy in MeV
!   enincmax     ! maximum incident energy
!   enincmin     ! minimum incident energy
!   Ninc         ! number of incident energies
! Constants
!   Emaxtalys    ! maximum acceptable energy for TALYS
!   parsym       ! symbol of particle
!
! *** Declaration of local data
!
  implicit none
  logical           :: fexist        ! flag for energy grid
  character(len=14) :: fname         ! filename
  integer           :: i             ! counter
  integer           :: j             ! counter
  integer           :: k             ! counter
  integer           :: iE1           ! counter of energies
  integer           :: iE2           ! counter for energies
  integer           :: istat         ! logical for file access
  integer           :: ktype         ! particle type
  integer           :: lenE          ! length of energy file
  integer           :: nen           ! energy counter
  integer           :: Nnen          ! total number of energies
  integer           :: pos           ! position
  integer           :: type          ! particle type
  real(sgl)         :: degrid        ! energy increment
  real(sgl)         :: E(numenin)    ! incident energy
  real(sgl)         :: Eeps          ! help variable
  real(sgl)         :: Ein           ! incident energy
  real(sgl)         :: Emax          ! maximal emission energy for particle channel
  real(sgl)         :: Emin          ! minimum energy
  real(sgl)         :: Erest         ! help variable
!
! ************************ Basic incident energy grid ******************
!
! The general incident energy grid we use is:
!
!   1.e-11,2.53e-8,  1.e-6,  1.e-5,  1.e-4 MeV
!   0.001,  0.002,  0.004,  0.007 MeV
!   0.01,   0.02,   0.04,   0.07 MeV
!   0.1-   1 MeV : dE= 0.1 MeV
!   1  -   8 MeV : dE= 0.2 MeV
!   8  -  15 MeV : dE= 0.5 MeV
!   15  -  30 MeV : dE= 1.0 MeV
!   30  -  60 MeV : dE= 2.0 MeV
!   60  -  80 MeV : dE= 5.0 MeV
!   80  - 160 MeV : dE=10.0 MeV
!   160 - 300 MeV : dE=20.0 MeV
!   300 - 600 MeV : dE=50.0 MeV
!   600 -1000 MeV : dE=100.0 MeV
!
! This grid ensures that the excitation functions are sufficiently smooth at low energies, while at higher energies
! a somewhat coarserenergy grid can be used.
! For various pre-defined energy grids (energy filenames), subsets of this grid can be taken.
!
  E(1) = 1.e-11
  E(2) = 2.53e-8
  E(3) = 1.e-6
  E(4) = 1.e-5
  E(5) = 1.e-4
  E(6) = 0.001
  E(7) = 0.002
  E(8) = 0.004
  E(9) = 0.007
  E(10) = 0.01
  E(11) = 0.02
  E(12) = 0.04
  E(13) = 0.07
  E(14) = 0.1
!
! Test for existence of pre-defined energy grid
!
  fexist = .false.
  do type = 0, 6
    if (fname(1:1) == parsym(type)) then
      ktype = type
      do j = 1, 14
        if (fname(j:j) == '-') then
          pos = j
          do k = pos + 2, 10
            if (fname(k:k+4) == '.grid') then
              lenE = k - 1
              read(fname(2:pos - 1), * , iostat = istat) iE1
              if (istat == 0) then
                read(fname(pos + 1:lenE), * , iostat = istat) iE2
                if (istat == 0) then
                  Emin = real(iE1)
                  Emax = real(iE2)
!
! Set basic energy grid
!
                  Ninc = 0
                  nen = 14
                  Ein = 0.1
                  degrid = 0.1
                  do
                    Ein = Ein + degrid
                    Eeps = Ein + 1.e-4
                    if (Ein > Emaxtalys) exit
                    if (nen == numenin) exit
                    nen = nen + 1
                    E(nen) = Ein
                    if (Eeps > 1.) degrid = 0.2
                    if (Eeps > 8.) degrid = 0.5
                    if (Eeps > 15.) degrid = 1.
                    if (Eeps > 30.) degrid = 2.
                    if (Eeps > 60.) degrid = 5.
                    if (Eeps > 80.) degrid = 10.
                    if (Eeps > 160.) degrid = 20.
                    if (Eeps > 300.) degrid = 50.
                    if (Eeps > 600.) degrid = 100.
                  enddo
                  Nnen = nen
!
! Fill energy grid
!
                  i = 0
                  do nen = 1, Nnen
                    Ein = E(nen)
                    Eeps = Ein + 1.e-4
                    if (Eeps < Emin) cycle
                    if (ktype /= 1) then
                      if (Eeps < 1.) cycle
                      Erest = Eeps - real(int(Eeps))
                      if (Erest > 0.1) cycle
                    endif
                    if (Ein <= Emax + 1.e-4) then
                      i = i + 1
                      eninc(i) = Ein
                    else
                      exit
                    endif
                  enddo
                  Ninc = i
                  enincmin = eninc(1)
                  enincmax = eninc(Ninc)
                  if (Ninc > 0) fexist = .true.
                endif
              endif
            endif
          enddo
        endif
      enddo
    endif
  enddo
  return
end subroutine incidentgrid
! Copyright A.J. Koning 2021
