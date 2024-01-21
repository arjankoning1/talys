subroutine decaydata
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Decay data
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
!   sgl          ! single precision kind
! All global variables
!   numelem      ! number of elements
! Variables for numerics
!   maxN         ! maximal number of neutrons away from initial compound nucleus
!   maxZ         ! maximal number of protons away from initial compound nucleus
! Variables for main input
!   Ninit        ! neutron number of initial compound nucleus
!   Zinit        ! charge number of initial compound nucleus
! Variables for files
!   path         ! directory containing files to be read
! Constants
!   nuc          ! symbol of nucleus
! Variables for decay data
!   daysec       ! number of seconds in a day
!   hoursec      ! number of seconds in an hour
!   lambda       ! decay rate per isotope
!   minutesec    ! number of seconds in a minute
!   rtyp         ! type of beta decay, beta - : 1 , beta + : 2 (from ENDF format)
!   Td           ! half life per time unit
!   Thalf        ! half life of nuclide in sec.
!   yearsec      ! number of seconds in a year
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist       ! logical to determine existence
  character(len=8)  :: decaychar    ! help variable
  character(len=132):: decayfile    ! decay data file
  character(len=80) :: string       ! line with parameter value
  integer           :: Abegin       ! first A to be included
  integer           :: Aend         ! last A to be included
  integer           :: i            ! counter
  integer           :: ia           ! mass number from abundance table
  integer           :: iline        ! number of lines
  integer           :: is           ! isotope counter: -1=total, 0=ground state 1=isomer
  integer           :: istat        ! logical for file access
  integer           :: N            ! neutron number of residual nucleus
  integer           :: Nbegin       ! first N to be included
  integer           :: NC           ! number of lines for MF/MT number
  integer           :: Nend         ! maximal neutron number
  integer           :: Nix          ! neutron number index for residual nucleus
  integer           :: Z            ! charge number of target nucleus
  integer           :: Zix          ! charge number index for residual nucleus
  real(sgl)         :: rrt          ! type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
  real(sgl)         :: TT           ! help variable
!
! ******************************** Time units **************************
!
  minutesec = 60.
  hoursec = 60. * minutesec
  daysec = 24. * hoursec
  yearsec = 365.25 * daysec
!
! ************** Find half lives for all involved nuclides *************
!
! Read decay constants from JEFF-3.1.1 radioactive decay data library
!
  do Zix = - 1, maxZ
    Z = Zinit - Zix
    Z = min(Z, numelem)
    Nbegin = Ninit - maxN
    Nend = Ninit
    Abegin = Z + Nbegin
    Aend = Z + Nend
    decaychar = trim(nuc(Z))//'.decay'
    decayfile = trim(path)//'decay/'//decaychar
    inquire (file = decayfile, exist = lexist)
    if (.not.lexist) cycle
    open (unit = 1, file = decayfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(decayfile, istat)
    do
      read(1, '(a80)', iostat = istat) string
      if (istat == -1) exit
      if (string(72:80) /= '1451    5') cycle
      read(string(8:10), '(i3)') ia
      if (ia < Abegin) cycle
      if (ia > Aend) exit
      N = ia - Z
      Nix = Ninit - N
      is = - 1
      if (string(11:11) == 'M') is = 1
      if (string(11:11) == 'N') is = 2
      do
        read(1, '(a80)', iostat = istat) string
        if (istat == -1) exit
        if (string(72:80) == '8457    2') exit
      enddo
      read(string(1:11), * ) Thalf(Zix, Nix, is)
      read(string(45:55), '(i11)') NC
      iline = 1 + (NC - 1) / 6
      do i = 1, iline
        read(1, '(a80)', iostat = istat) string
        if (istat == -1) exit
      enddo
      read(1, '(a80)', iostat = istat) string
      if (istat == -1) exit
      read(1, '(a80)', iostat = istat) string
      if (istat == -1) exit
      read(string(1:11), * ) rrt
      rtyp(Zix, Nix, is) = int(rrt)
      if (Thalf(Zix, Nix, is) == 0.) Thalf(Zix, Nix, is) = 1.e30
!
! Write half life in years, days, etc.
!
      TT = Thalf(Zix, Nix, is)
      Td(Zix, Nix, is, 1) = int(TT / yearsec)
      TT = TT - Td(Zix, Nix, is, 1) * yearsec
      Td(Zix, Nix, is, 2) = int(TT / daysec)
      TT = TT - Td(Zix, Nix, is, 2) * daysec
      Td(Zix, Nix, is, 3) = int(TT / hoursec)
      TT = TT - Td(Zix, Nix, is, 3) * hoursec
      Td(Zix, Nix, is, 4) = int(TT / minutesec)
      TT = TT - Td(Zix, Nix, is, 4) * minutesec
      Td(Zix, Nix, is, 5) = int(TT)
!
! Calculate decay rates.
! Nuclides with a lifetime longer 1.e17 sec are considered stable.
!
      if (Thalf(Zix, Nix, is) > 1.e17) then
        lambda(Zix, Nix, is) = 0.
      else
        lambda(Zix, Nix, is) = log(2.) / Thalf(Zix, Nix, is)
      endif
      if (is ==  - 1) then
        rtyp(Zix, Nix, 0) = rtyp(Zix, Nix, - 1)
        lambda(Zix, Nix, 0) = lambda(Zix, Nix, - 1)
        Thalf(Zix, Nix, 0) = Thalf(Zix, Nix, - 1)
        do i = 1, 5
          Td(Zix, Nix, 0, i) = Td(Zix, Nix, - 1, i)
        enddo
      endif
    enddo
    close (unit = 1)
    do ia = Abegin, Aend
      N = ia - Z
      Nix = Ninit - N
      if (Thalf(Zix, Nix, 1) > 1.e17) then
        rtyp(Zix, Nix, 1) = rtyp(Zix, Nix, - 1)
        lambda(Zix, Nix, 1) = lambda(Zix, Nix, - 1)
        Thalf(Zix, Nix, 1) = Thalf(Zix, Nix, - 1)
        do i = 1, 5
          Td(Zix, Nix, 1, i) = Td(Zix, Nix, - 1, i)
        enddo
      endif
    enddo
  enddo
  return
end subroutine decaydata
! Copyright A.J. Koning 2021
