subroutine directread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read ECIS results for direct cross section
!
! Author    : Arjan Koning and Eric Bauge
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl              ! single precision kind
!   dbl              ! double precision kind
! All global variables
!   numlev2          ! maximum number of levels
! Variables for numerics
!   nangle           ! number of angles
!   nanglecont       ! number of angles for continuum
! Variables for main input
!   k0               ! index of incident particle
!   Ltarget          ! excited level of target
! Variables for energy grid
!   ecisstatus       ! status of ECIS file
! Variables for energies
!   eninccm          ! center - of - mass incident energy in MeV
!   eoutdis          ! outgoing energy of discrete state reaction
!   flaggiant        ! flag for collective contribution from giant resonances
! Variables for incident channel
!   directad         ! direct angular distribution
!   dleg             ! direct reaction Legendre coefficient
!   dorigin          ! origin of direct cross section (Direct or Preeq)
!   xscollconttot    ! total collective cross section in the continuum
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsdirdiscsum     ! total direct cross section
!   xsdirdisctot     ! direct cross section summed over discrete states
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   parskip          ! logical to skip outgoing particle
!   Zindex           ! charge number index for residual nucleus
! Variables for giant resonances
!   grcollad         ! giant resonance angular distribution
!   xsgrcoll         ! giant resonance cross section
! Constants
!   parA             ! mass number of particle
! Variables for deformation parameters
!   betagr           ! deformation parameter for giant resonance
!   deform           ! deformation parameter
!   Egrcoll          ! energy of giant resonance
! Variables for levels
!   edis             ! energy of level
! Variables for level density
!   Nlast            ! last discrete level
!
! *** Declaration of local data
!
  implicit none
  logical   :: lexist      ! logical to determine existence
  integer   :: i           ! counter
  integer   :: iang        ! running variable for angle
  integer   :: iS          ! counter
  integer   :: istat       ! logical for file access
  integer   :: itype       ! help variable
  integer   :: k           ! designator for particle
  integer   :: l           ! multipolarity
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: NL          ! last discrete level
  integer   :: nleg        ! number of Legendre coefficients
  integer   :: nS          ! number of states
  integer   :: type        ! particle type
  integer   :: Zix         ! charge number index for residual nucleus
  real(sgl) :: levelenergy ! energy of level
  real(sgl) :: xsdwbatot   ! direct DWBA cross section summed over discrete states
  real(dbl) :: ddl         ! direct reaction Legendre coefficient
  real(dbl) :: xs          ! help variable
!
! **************** Read direct, inelastic cross sections ***************
!
  inquire (file = 'ecis.dirang', exist = lexist)
  if ( .not. lexist) return
  open (unit = 8, file = 'ecis.dirang', status = 'unknown')
  open (unit = 9, file = 'ecis.dirleg', status = 'unknown')
  open (unit = 10, file = 'ecis.dirin', status = 'unknown')
  do type = k0, k0
    if (parskip(type)) cycle
    Zix = Zindex(0, 0, type)
    Nix = Nindex(0, 0, type)
    NL = Nlast(Zix, Nix, 0)
!
! 1. Direct collective states
!
    do i = 0, numlev2
      if (i == Ltarget .and. type == k0) cycle
      if (deform(Zix, Nix, i) == 0.) cycle
      if (eoutdis(type, i) <= 0.) cycle
      levelenergy = edis(Zix, Nix, i)
      if (eninccm <= levelenergy + 0.1 * parA(type)) cycle
      read(10, '()')
      read(10, * ) xs
      xsdirdisc(type, i) = real(xs)
      if (i <= NL) dorigin(type, i) = 'Direct'
!
! ******************* Direct reaction Legendre coefficients ************
!
! We read the Legendre coefficients for the direct component of the reaction only, the compound nucleus coefficients are
! calculated by TALYS later on.
!
      read(9, '(55x, i5)') nS
      do iS = 1, nS
        if (iS == 1) then
          read(9, '(5x, i5)')  nleg
          do k = 1, nleg
            read(9, '()')
          enddo
        else
          read(9, '(5x, i5)')  nleg
          do k = 1, nleg
            read(9, '(5x, i5, e20.10)') l, ddl
            if (i <= NL) dleg(type, i, l) = real(ddl)
          enddo
        endif
      enddo
!
! ************************ Angular distributions ***********************
!
! We first skip the elastic angular distribution
!
      read(8, '()')
      read(8, '(12x, i3)') nS
      do iang = 0, nangle
        do k = 1, nS
          read(8, '()')
    enddo
  enddo
!
! Read direct angular distribution
!
      read(8, '(12x, i3)') nS
      do iang = 0, nangle
        do k = 1, nS
          read(8, '(i3, 12x, e12.5)', iostat = istat) itype, xs
          if (istat /= 0) cycle
          if (itype == 0) directad(type, i, iang) = real(xs)
        enddo
      enddo
    enddo
!
! 2. Giant resonance states
!
    if (flaggiant .and. type == k0) then
      do l = 0, 3
        do i = 1, 2
          if (betagr(l, i) == 0.) cycle
          levelenergy = Egrcoll(l, i)
          if (eninccm <= levelenergy + 0.1 * parA(type)) cycle
!
! Giant resonance cross section
!
          read(10, '()')
          read(10, * ) xs
          xsgrcoll(k0, l, i) = real(xs)
!
! Giant resonance angular distribution
!
          read(8, '()')
          read(8, '(12x, i3)') nS
          do iang = 0, nanglecont
            do k = 1, nS
              read(8, '()')
            enddo
          enddo
          read(8, '(12x, i3)') nS
          do iang = 0, nanglecont
            do k = 1, nS
              read(8, '(i3, 12x, e12.5)', iostat = istat) itype, xs
              if (istat /= 0) cycle
              if (itype == 0) grcollad(k0, l, i, iang) = real(xs)
            enddo
          enddo
        enddo
      enddo
    endif
!
! ************* Create total direct inelastic cross section ************
!
    xsdwbatot = 0.
    do i = 0, numlev2
      if (i == 0 .and. type == k0) cycle
      if (deform(Zix, Nix, i) /= 0.) then
        if (i <= NL) then
          xsdwbatot = xsdwbatot + xsdirdisc(type, i)
        else
          xscollconttot(type) = xscollconttot(type) + xsdirdisc(type, i)
        endif
      endif
    enddo
    xsdirdisctot(type) = xsdirdisctot(type) + xsdwbatot
    xsdirdiscsum = xsdirdiscsum + xsdwbatot
  enddo
  close (unit = 8, status = ecisstatus)
  close (unit = 9, status = ecisstatus)
  close (unit = 10, status = ecisstatus)
  open (unit = 3, file = 'ecis.dircs', status = 'unknown')
  close (unit = 3, status = ecisstatus)
  return
end subroutine directread
! Copyright A.J. Koning 2021
