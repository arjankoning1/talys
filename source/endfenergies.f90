subroutine endfenergies
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Energy grid for ENDF-6 file
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
!   dbl              ! double precision kind
! All global variables
!   numen6           ! number of energies for ENDF6 energy grid
! Variables for basic parameters
!   eninclow         ! minimal incident energy for nuclear model calculations
! Variables for input energies
!   eninc            ! incident energy in MeV
!   enincmax         ! maximum incident energy
!   Estop            ! incident energy above which TALYS stops
!   Ninc           ! number of incident energies
! Variables for main input
!   k0               ! index of incident particle
!   Ltarget          ! excited level of target
! Variables for basic reaction
!   ompenergyfile    ! file with energies for OMP calculation (ENDF files)
! Variables for energies
!   Ethrexcl         ! threshold incident energy for exclusive channel
!   Ethresh          ! threshold incident energy for residual nucleus
!   idchannel        ! identifier for exclusive channel
! Variables for exclusive channels
!   idnum            ! counter for exclusive channel
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   Zindex           ! charge number index for residual nucleus
! Constants
!   Emaxtalys        ! maximum acceptable energy for TALYS
! Variables for level density
!   Nlast            ! last discrete level
! Variables for ENDF data
!   e6               ! energies of ENDF - 6 energy grid in MeV
!   nen6             ! total number of energies
! Error handling
!   range_integer_error     ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  integer   :: i      ! counter
  integer   :: idc    ! help variable
  integer   :: istat  ! logical for file access
  integer   :: k      ! designator for particle
  integer   :: n      ! exciton number
  integer   :: nen    ! energy counter
  integer   :: nen0   ! energy counter
  integer   :: nen2   ! energy counter
  integer   :: nenint ! energy index
  integer   :: nensub ! intermediate number of incident energies
  integer   :: nex    ! excitation energy bin of compound nucleus
  integer   :: Nix    ! neutron number index for residual nucleus
  integer   :: type   ! particle type
  integer   :: Zix    ! charge number index for residual nucleus
  real(dbl) :: degrid ! energy increment
  real(dbl) :: e6tmp  ! help variable
  real(dbl) :: ee     ! energy
  real(dbl) :: Eeps   ! help variable
  real(dbl) :: Ein    ! incident energy
  real(dbl) :: emax   ! maximal emission energy within bin decay
  real(dbl) :: ratio  ! if 0<ratio<1 then x is between xl1 and xl2
!
! ************************ Basic ENDF-6 energy grid ********************
!
! The basic ENDF-6 energy grid we use is:
!
!   0.001 -   0.01 MeV  : dE= 0.001 MeV
!   0.01  -   0.1 MeV   : dE= 0.01 MeV
!   0.1   -   1 MeV     : dE= 0.05 MeV
!   1     -   4 MeV     : dE= 0.1 MeV
!   4     -   7 MeV     : dE= 0.2 MeV
!   7     -  20 MeV     : dE= 0.5 MeV
!   20     - 100 MeV    : dE= 1.0 MeV
!   100     - 200 MeV   : dE= 2.0 MeV
!   200     - 300 MeV   : dE= 5.0 MeV
!   above 300 MeV       : dE=10.0 MeV
!
! This grid ensures that the total, elastic and reaction cross section are calculated on a sufficiently precise energy grid.
! However, this grid can be overwritten using the 'ompenergyfile'.
!
  e6(1) = eninclow
  nen = 1
  if (ompenergyfile(1:1) /= ' ') then
    open (unit = 2, file = ompenergyfile, status = 'old')
    do
      read(2, * , iostat = istat) Ein
      if (istat == -1) exit
      if (istat /= 0) call read_error(ompenergyfile, istat)
      if (Ein <= eninclow) cycle
      nen = nen + 1
      call range_integer_error(ompenergyfile, nen, 1, numen6)
      e6(nen) = Ein
    enddo
    close (unit = 2)
!
! Sort incident energies in ascending order and remove double points
!
    do i = 1, nen
      do k = 1, i
        if (e6(i) >= e6(k)) cycle
        e6tmp = e6(i)
        e6(i) = e6(k)
        e6(k) = e6tmp
      enddo
    enddo
    n = nen
    do i = 1, nen - 1
      if (e6(i) == e6(i + 1)) then
        do k = i + 1, nen
          e6(k) = e6(k + 1)
        enddo
        n = n - 1
      endif
    enddo
    nen = n
  else
    Ein = 0.
    degrid = 0.001
    do
      Ein = Ein + degrid
      if (Ein > e6(1)) then
        nen = nen + 1
        e6(nen) = Ein
      endif
      Eeps = Ein + 1.e-4
      if (Eeps > enincmax) exit
      if (Eeps > 0.01) degrid = 0.01
      if (Eeps > 0.1) degrid = 0.05
      if (Eeps > 1.) degrid = 0.1
      if (Eeps > 4.) degrid = 0.2
      if (Eeps > 8.) degrid = 0.5
      if (Eeps > 20.) degrid = 1.
      if (Eeps > 100.) degrid = 2.
      if (Eeps > 200.) degrid = 5.
      if (Eeps > 300.) degrid = 10.
      if (Eeps > Emaxtalys) exit
    enddo
  endif
!
! Add possible denser energy point from original energy grid
!
  nensub = nen
  if (k0 == 1) then
Loop1:    do nen0 = 1, Ninc - 1
      do nen = 1, nensub
        ratio = e6(nen) / eninc(nen0)
        if (ratio >= 0.999 .and. ratio <= 1.001) cycle Loop1
      enddo
      do nen = 1, nensub
        if (eninc(nen0) > e6(nen) .and. eninc(nen0) <= e6(nen + 1)) then
          nenint = nen + 1
          nensub = nensub + 1
          call range_integer_error(ompenergyfile, nensub, 1, numen6)
          do nen2 = nensub, nenint + 1, - 1
            e6(nen2) = e6(nen2 - 1)
          enddo
          e6(nenint) = eninc(nen0)
          exit
        endif
      enddo
    enddo Loop1
    nen = nensub
!
! *************** Add partial thresholds to energy grid ****************
!
    emax = min(enincmax, 20.)
    do type = 1, 6
      Zix = Zindex(0, 0, type)
      Nix = Nindex(0, 0, type)
      do nex = 0, Nlast(Zix, Nix, 0)
        if (type == k0 .and. nex == Ltarget) cycle
        ee = Ethresh(Zix, Nix, nex)
        if (ee > eninclow .and. ee <= emax) then
          nen = nen + 1
          call range_integer_error(ompenergyfile, nen, 1, numen6)
          e6(nen) = ee
        endif
      enddo
    enddo
    do idc = 0, idnum
      if (idchannel(idc) == 100000) cycle
      if (idchannel(idc) == 10000) cycle
      if (idchannel(idc) == 1000) cycle
      if (idchannel(idc) == 100) cycle
      if (idchannel(idc) == 10) cycle
      if (idchannel(idc) == 1) cycle
      ee = Ethrexcl(idc, 0)
      if (ee > eninclow .and. ee <= emax) then
        nen = nen + 1
        call range_integer_error(ompenergyfile, nen, 1, numen6)
        e6(nen) = ee
      endif
    enddo
  endif
  nen6 = nen
!
! *************************** Sort energies ****************************
!
  do nen = 1, nen6
    do nen2 = nen, nen6
      if (e6(nen) <= e6(nen2)) cycle
      e6tmp = e6(nen)
      e6(nen) = e6(nen2)
      e6(nen2) = e6tmp
    enddo
  enddo
!
! Remove incident energies above energy given by Estop
!
  do i = 1, nen
    if (e6(i) > Estop) then
      nen6 = i - 1
      exit
    endif
  enddo
  return
end subroutine endfenergies
! Copyright A.J. Koning 2021
