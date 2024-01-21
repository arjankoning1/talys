subroutine totalxs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Total cross sections
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
! Variables for basic reaction
!   flagffruns      ! flag to designate subsequent evaporation of fission products
! Variables for basic reaction
!   flagastro       ! flag for calculation of astrophysics reaction rate
!   flagmassdis     ! flag for calculation of fission fragment mass yields
!   flagchannels    ! flag for exclusive channels calculation
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
! Variables for input energies
!   flaginitpop     ! flag for initial population distribution
!   nin             ! counter for incident energy
! Variables for energy grid
!   Einc     ! incident energy in MeV
! Variables for fission
!   flagfission     ! flag for fission
!   fymodel         ! fission yield model, 1: Brosa 2: GEF 3: GEF + TALYS
! Variables for total cross sections
!   xsexclcont      ! exclusive single channel cross section for continuum
!   xsexclusive     ! exclusive single channel cross section
!   xsfistot        ! total fission cross section
!   xsfistot0       ! total fission cross section
! Variables for incident channel
!   multiplicity    ! particle multiplicity
!   xsparticle      ! total particle production cross section
! Variables for energies
!   idchannel       ! identifier for exclusive channel
! Variables for exclusive channels
!   idnum           ! counter for exclusive channel
!   xschannel       ! channel cross section
! Variables for multiple emission
!   xsinitpop       ! initial population cross section
!   xsfeed          ! cross section from compound to residual nucleus
! Variables for binary reactions
!   xsconttot       ! total cross section for continuum
!   xsdisctot       ! total cross section summed over discrete states
!   xsnonel         ! non - elastic cross
! Variables for nuclides
!   parskip         ! logical to skip outgoing particle
! Variables for mass distribution
!   nubar          ! average nu
! Variables for astro
!   xsastrofis      ! astrophysical fission cross section
!
! *** Declaration of local data
!
  implicit none
  integer :: idc       ! help variable
  integer :: ident     ! exclusive channel identifier
  integer :: Ncomp     ! neutron number index for compound nucleus
  integer :: type      ! particle type
  integer :: Zcomp     ! proton number index for compound nucleus
  real    :: nubarWahl
!
! *********************** Specific cross sections **********************
!
  if (flagchannels) then
    do type = 0, 6
      if (parskip(type)) cycle
      if (type == 0) then
        xsexclusive(0) = xschannel(0)
      else
        ident = 10 **(6 - type)
        do idc = 0, idnum
          if (idchannel(idc) == ident) then
            xsexclusive(type) = xschannel(idc)
            exit
          endif
        enddo
      endif
      xsdisctot(type) = min(xsexclusive(type), xsdisctot(type))
      if (xsconttot(type) == 0.) then
        xsexclcont(type) = 0.
      else
        xsexclcont(type) = max(xsexclusive(type) - xsdisctot(type), 0.)
      endif
    enddo
  endif
!
! *************** Total particle production cross sections *************
!
  do type = 0, 6
    if (parskip(type)) cycle
    xsparticle(type) = 0.
    do Zcomp = 0, maxZ
      do Ncomp = 0, maxN
        xsparticle(type) = xsparticle(type) + xsfeed(Zcomp, Ncomp, type)
      enddo
    enddo
    if (flaginitpop) then
      if (xsinitpop /= 0.) multiplicity(type) = xsparticle(type) / xsinitpop
    else
      if (xsnonel /= 0.) multiplicity(type) = xsparticle(type) / xsnonel
    endif
  enddo
!
! ******************* Total fission cross sections ********************
!
  xsfistot = 0.
  if (flagfission) then
    do Zcomp = 0, maxZ
      do Ncomp = 0, maxN
        xsfistot = xsfistot + xsfeed(Zcomp, Ncomp, - 1)
      enddo
    enddo
    if ( .not. flagffruns) xsfistot0 = xsfistot
    if (flagastro) xsastrofis(nin) = xsfistot
    if (.not. (flagmassdis .and. fymodel >= 3)) nubar(1) = nubarWahl(Einc)
  endif
  return
end subroutine totalxs
! Copyright A.J. Koning 2021
