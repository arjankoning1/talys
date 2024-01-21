subroutine prodres
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Residual production cross sections
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
!   numenrp       ! number of incident energies for residual products
!   numisom       ! number of isomers
! Variables for medical isotope production
!   Eback         ! lower end of energy range in MeV for isotope
!   Ebeam         ! incident energy in MeV for isotope production
! Variables for numerics
!   maxN          ! maximal number of neutrons away from initial compound nucleus
!   maxZ          ! maximal number of protons away from initial compound nucleus
! Variables for main input
!   Atarget       ! mass number of target nucleus
!   Ninit         ! neutron number of initial compound nucleus
!   Zinit         ! charge number of initial compound nucleus
!   Ztarget       ! charge number of target nucleus
! Variables for nuclides
!   strucexist    ! flag to state whether structure info for this nucleus exists
! Constants
!   iso           ! counter for isotope
!   natstring     ! string extension for file names
! Variables for levels
!   Lisomer       ! level number of isomer
!   Nisomer       ! number of isomers for this nuclide
! Variables for isotope production
!   Erp           ! incident energy
!   Nenrp         ! number of incident energies for residual production cross
!   xsrp          ! residual production cross section in mb
! Variables for existence libraries
!   prodexist     ! logical to determine existence of residual production
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: flagpositive    ! flag for existence of non-zero cross sections
  logical           :: lexist          ! logical to determine existence
  character(len=16) :: rpfile          ! file with residual production cross sections
  character(len=18) :: nonfile         ! file with nonelastic cross sections
  character(len=132) :: string          ! line with parameter value
  integer           :: A               ! mass number of target nucleus
  integer           :: iE              ! energy counter
  integer           :: is              ! isotope counter: -1=total, 0=ground state 1=isomer
  integer           :: istat           ! logical for file access
  integer           :: N               ! neutron number of residual nucleus
  integer           :: nen             ! energy counter
  integer           :: Nix             ! neutron number index for residual nucleus
  integer           :: Z               ! charge number of target nucleus
  integer           :: Zix             ! charge number index for residual nucleus
  real(sgl)         :: E               ! incident energy
  real(sgl)         :: xs              ! help variable
!
! **************** Read residual production cross sections *************
!
  do Zix = -1, maxZ
    Z = Zinit - Zix
    do Nix = - 1, maxN
      do is = - 1, numisom
        prodexist(Zix, Nix, is) = .false.
        Nenrp(Zix, Nix, is) = 0
        do nen = 1, numenrp
          Erp(Zix, Nix, is, nen) = 0.
          xsrp(Zix, Nix, is, nen) = 0.
        enddo
      enddo
!
! For convenience in later loops, we store the non-elastic cross section in the (-1,-1) element of xsrp.
! In the next loop, we subtract the inelastic cross section from this, i.e. xsrp will contain all non-elastic cross sections
! other than inelastic.
!
      if (Zix ==  - 1 .and. Nix ==  - 1) then
        is = - 1
        prodexist(Zix, Nix, is) = .true.
        nonfile = 'nonelastic.tot'//natstring(iso)
        inquire (file = nonfile, exist = lexist)
        if (.not. lexist) then
          write(*,'(" TALYS-error: non-elastic cross section file nonelastic.tot does not exist")')
          stop
        endif
        open (unit = 1, file = nonfile, status = 'old')
        iE = 0
        do
          read(1, '(a)', iostat = istat) string
          if (istat == -1) exit
          if (string(1:1) == '#') cycle
          read(string, * ) E, xs
          iE = iE + 1
          if (iE > numenrp) then
            write(*, '(" TALYS-warning: number of energy points exceeds ",i3,", incomplete integration")') numenrp
            exit
          endif
          Erp(Zix, Nix, is, iE) = E
          xsrp(Zix, Nix, is, iE) = xs
        enddo
        close (unit = 1)
        Nenrp(Zix, Nix, is) = iE
        cycle
      endif
      if (Zix ==  - 1) cycle
      if (Nix ==  - 1) cycle
      if ( .not. strucexist(Zix, Nix)) call levels(Zix, Nix)
!
! Residual production cross sections
!
      N = Ninit - Nix
      A = Z + N
      rpfile = 'rp000000.L00'//natstring(iso)
      write(rpfile(3:5), '(i3.3)') Z
      write(rpfile(6:8), '(i3.3)') A
!
! Check for existence of isomeric cross sections.
! The total cross section is considered for is=-1.
!
      do is = - 1, Nisomer(Zix, Nix)
        if (is == -1) rpfile(10:12) = 'tot'
        if (is == 0) rpfile(10:12) = 'L00'
        if (is > 0) write(rpfile(11:12), '(i2.2)') Lisomer(Zix, Nix, is)
        inquire (file = rpfile, exist = lexist)
        if (lexist) then
          prodexist(Zix, Nix, is) = .true.
          flagpositive = .false.
          iE = 0
          open (unit = 1, file = rpfile, status = 'old')
          do
            read(1, '(a)', iostat = istat) string
            if (istat == -1) exit
            if (string(1:1) == '#') cycle
            read(string, * , iostat = istat) E, xs
            if (istat /= 0) call read_error(rpfile, istat)
            if (E >= Eback .and. E <= Ebeam .and. xs > 0.) flagpositive = .true.
            iE = iE + 1
            if (iE > numenrp) then
              write(*, '(" TALYS-warning: number of energy points exceeds ",i3,", incomplete integration")') numenrp
              exit
            endif
            Erp(Zix, Nix, is, iE) = E
            xsrp(Zix, Nix, is, iE) = xs
!
! Subtract inelastic cross section from non-elastic cross section
!
            if (Z == Ztarget .and. A == Atarget .and. is <= 0) &
 &            xsrp( - 1, - 1, - 1, iE) = xsrp( - 1, - 1, - 1, iE) - xsrp(Zix, Nix, is, iE)
          enddo
          close (unit = 1)
          Nenrp(Zix, Nix, is) = iE
          if ( .not. flagpositive) prodexist(Zix, Nix, is) = .false.
        endif
      enddo
    enddo
  enddo
  return
end subroutine prodres
! Copyright A.J. Koning 2021
