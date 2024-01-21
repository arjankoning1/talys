subroutine phdensitytable(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Tabulated particle-hole state densities
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
!   numconf      ! number of particle - hole combinations
!   numdens      ! number of energy points for tabulated level densities
! Variables for preequilibrium
!   flag2comp    ! flag for two - component pre - equilibrium mod
! Variables for nuclides
!   AA           ! mass number of residual nucleus
!   ZZ           ! charge number of residual nucleus
! Constants
!   nuc          ! symbol of nucleus
! Variables for files
!   path         ! directory containing files to be read
! Variables for particle-hole state densities
!   hhtable      ! hole number from table
!   hnutable     ! neutron hole number from table
!   hpitable     ! proton hole number from table
!   nenphdens    ! number of energies for p - h stat
!   Nphconf1     ! number of 1 - component p - h confi
!   Nphconf2     ! number of 2 - component p - h confi
!   phexist1     ! flag for existence of p - h state
!   phexist2     ! flag for existence of p - h state
!   phtable1     ! p - h state density
!   phtable2     ! p - h state density
!   pnutable     ! neutron particle number from ta
!   ppitable     ! proton particle number from tab
!   pptable      ! particle number from table
! Variables for level density
!   edens        ! energy grid for tabulated level den
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist                        ! logical to determine existence
  character(len=5)  :: denchar                       ! string for level density file
  character(len=132):: denfile                       ! level density parameter file
  integer           :: A                             ! mass number of target nucleus
  integer           :: i                             ! counter
  integer           :: ia                            ! mass number from abundance table
  integer           :: istat                         ! logical for file access
  integer           :: nex                           ! excitation energy bin of compound nucleus
  integer           :: nex2                          ! counter
  integer           :: Nix                           ! neutron number index for residual nucleus
  integer           :: Z                             ! charge number of target nucleus
  integer           :: Zix                           ! charge number index for residual nucleus
  real(sgl)         :: phden(0:numdens, numconf)     ! particle-hole density
  real(sgl)         :: phden2(0:numdens, numconf)    ! particle-hole density
  real(sgl)         :: sum                           ! help variable
!
! ******** Tabulated particle-hole level densities from Hilaire ********
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  denchar = trim(nuc(Z))//'.ld'
  denfile = trim(path)//'density/ph/'//denchar
!
! Check existence of file and read data from the tables.
!
  inquire (file = denfile, exist = lexist)
  if ( .not. lexist) return
  open (unit = 2, file = denfile, status = 'old', iostat = istat)
  if (istat /= 0) call read_error(denfile, istat)
  do
    read(2, '(/31x, i3)', iostat = istat) ia
    if (istat == -1) return
    if (istat /= 0) call read_error(denfile, istat)
    if (A == ia) then
      if (flag2comp) then
        do i = 1, Nphconf2
          phexist2(Zix, Nix, ppitable(i), hpitable(i), pnutable(i), hnutable(i)) = .true.
        enddo
      else
        do i = 1, Nphconf1
          phexist1(Zix, Nix, pptable(i), hhtable(i)) = .true.
        enddo
      endif
      read(2, '(///)')
      do nex = 1, nenphdens
        read(2, '(7x, 86(e9.2))', iostat = istat) (phden2(nex, i), i = 1, Nphconf2), (phden(nex, i), i = 1, Nphconf1)
        if (istat /= 0) call read_error(denfile, istat, ival = nex, xval = real(phden2(nex, i)))
      enddo
      read(2, '()')
      if (flag2comp) then
        do nex = 1, nenphdens
          do i = 1, Nphconf2
            sum = 0.
            do nex2 = max(nex - 2, 1), min(nex + 2, nenphdens)
              sum = sum + phden2(nex2, i) * 0.5 * (edens(min(nex2 + 1, nenphdens)) - edens(max(nex2 - 1, 1)))
            enddo
            phtable2(Zix, Nix, ppitable(i), hpitable(i), pnutable(i), hnutable(i), nex) = &
 &            sum / (edens(min(nex + 2, nenphdens)) - edens(max(nex - 2, 1)))
          enddo
        enddo
      else
        do nex = 1, nenphdens
          do i = 1, Nphconf1
            sum = 0.
            do nex2 = max(nex - 2, 1), min(nex + 2, nenphdens)
              sum = sum + phden(nex2, i) * 0.5 * (edens(min(nex2 + 1, nenphdens)) - edens(max(nex2 - 1, 1)))
            enddo
            phtable1(Zix, Nix, pptable(i), hhtable(i), nex) = sum / (edens(min(nex + 2, nenphdens)) - edens(max(nex - 2, 1)))
          enddo
        enddo
      endif
      exit
    else
      do nex = 1, nenphdens + 5
        read(2, '()')
      enddo
    endif
  enddo
  return
end subroutine phdensitytable
! Copyright A.J. Koning 2021
