subroutine dwbaout(itype, type, nen1, nen2)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of ECIS results for DWBA for MSD
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
! Variables for output
!   flagddx       ! flag for output of double - differential cross sections
! Variables for numerics
!   nanglecont    ! number of angles for continuum
! Variables for energy grid
!   anglecont     ! angle in degrees for continuum
! Constants
!   parsym        ! symbol of particle
! Variables for MSD
!   Emsdin        ! incident MSD energy
!   Emsdout       ! outgoing MSD energy
!   Exmsd         ! excitation energy for MSD energy grid
!   maxJmsd       ! maximal spin for MSD calculation
!   xsdw          ! DWBA angular distribution as a function of incident energy, outgoing energy, ang. mom. and angle
!   xsdwin        ! DWBA cross section as a function of incident energy, outgoing energy and angular momentum
!
! *** Declaration of local data
!
  implicit none
  character(len=72) :: ofor1    ! help variable
  character(len=72) :: ofor2    ! help variable
  integer           :: iang     ! running variable for angle
  integer           :: itype    ! help variable
  integer           :: J        ! spin of level
  integer           :: nen1     ! energy counter
  integer           :: nen2     ! energy counter
  integer           :: type     ! particle type
!
! ********************** Write DWBA cross sections *********************
!
  ofor1='(/,"  Angle",3x,n(" J =",i2,7x))'
  ofor2='(/,10x,n(" J =",i2,7x))'
  write(ofor1(17:17), '(i1)') min(maxJmsd+1, 8)
  write(ofor2(8:8), '(i1)') min(maxJmsd+1, 8)
  write(*, '(/, " (", a1, ",", a1, ") DWBA cross sections (mb/Sr) for E-in=", f6.2, " MeV, Ex=", f6.2, " MeV, E-out=", f6.2, &
 &  " MeV")') parsym(itype), parsym(type), Emsdin, Exmsd, Emsdout
  if (flagddx) then
    write( * , fmt = ofor1) (J, J = 0, min(maxJmsd, 7))
    write(*, '()')
    do iang = 0, nanglecont
      write(*, '(1x, f5.1, 8es13.5)') anglecont(iang), (xsdw(nen1, nen2, J, iang, 0), J = 0, min(maxJmsd, 7))
    enddo
  endif
  write( * , fmt = ofor2) (J, J = 0, min(maxJmsd, 7))
  write(*, '(/, " Angle", 8es13.5)') (xsdwin(nen1, nen2, J, 0), J = 0, min(maxJmsd, 7))
  write(*, '(" integr.")')
  return
end subroutine dwbaout
! Copyright A.J. Koning 2021
