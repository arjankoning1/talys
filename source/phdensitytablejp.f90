subroutine phdensitytablejp(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Tabulated spin- and parity-dependent particle-hole level
!
! Author    : Stephane Goriely and Arjan Koning
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
!   sgl             ! single precision kind
!   dbl             ! double precision kind
! All global variables
!   numJ            ! maximum J - value
! Variables for gamma rays
!   ldmodelracap    ! level density model for direct radiative capture
! Variables for main input
!   k0              ! index of incident particle
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   nuc             ! symbol of nucleus
! Variables for files
!   path            ! directory containing files to be read
! Variables for particle-hole state densities
!   nenphdens       ! number of energies for p - h stat
! Variables for direct capture initialization
!   edensphjp       ! energy grid of ph spin - and parity - depe
!   phdensjp        ! ph spin - and parity - dependent level den
!   phdenstot       ! total ph level density from table
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist           ! logical to determine existence
  character(len=5)  :: denchar          ! string for level density file
  character(len=6)  :: phdir            ! directory for particle-hole density
  character(len=132):: denfile          ! level density parameter file
  integer           :: A                ! mass number of target nucleus
  integer           :: ia               ! mass number from abundance table
  integer           :: istat            ! logical for file access
  integer           :: J                ! spin of level
  integer           :: nex              ! excitation energy bin of compound nucleus
  integer           :: Nix              ! neutron number index for residual nucleus
  integer           :: parity           ! parity
  integer           :: Z                ! charge number of target nucleus
  integer           :: Zix              ! charge number index for residual nucleus
  real(sgl)         :: ephjpgrid        ! energy for level density grid
  real(dbl)         :: ld2j1(0:numJ)    ! spin dependent level density
  real(dbl)         :: ldtot            ! total level density
!
! *********** Tabulated level densities from Goriely *******************
!
  if (ldmodelracap /= 1) return
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  if (k0 == 1) then
    phdir = 'ph0011'
  else
    ldmodelracap = 2
    return
  endif
  denchar = trim(nuc(Z))//'.ph'
  denfile = trim(path)//'density/phjp/'//phdir//'/'//denchar
!
! Check existence of file and read data from the tables.
!
  inquire (file = denfile, exist = lexist)
  if (lexist) then
    open (unit = 2, status = 'old', file = denfile)
    do parity = 1, - 1, - 2
      do
        read(2, '(/31x, i3//)', iostat = istat) ia
        if (istat /= 0) exit
        if (A == ia) then
          do nex = 1, nenphdens
            read(2, '(1x, f6.2, 17x, e9.2, 9x, 30e9.2)', iostat = istat) ephjpgrid, ldtot, (ld2j1(J), J = 0, 29)
            if (istat /= 0) call read_error(denfile, istat, ival = Z, xval = ephjpgrid)
!
! Determination of the mass-asymmetric enhancement factor for fission barrier
!
            phdenstot(Zix, Nix, nex) = phdenstot(Zix, Nix, nex) + ldtot
            edensphjp(Zix, Nix, nex) = ephjpgrid
            do J = 0, 29
              phdensjp(Zix, Nix, nex, J, parity) = ld2j1(J)
            enddo
          enddo
          read(2, '()')
          exit
        else
          do nex = 1, nenphdens + 1
            read(2, '()')
          enddo
        endif
      enddo
    enddo
    close (unit = 2)
    if (A /= ia) then
      write(*, '("Input ph file not available:", a)') trim(denfile)
      write(*, '(" For A=", i3, " --> Change of ldmodelracap=1 to ldmodelracap=2")') A
      ldmodelracap = 2
    endif
  else
    write(*, '("Input ph file not available:", a)') trim(denfile)
    write(*, '("  Change of ldmodelracap=1 to ldmodelracap=2")')
    ldmodelracap = 2
  endif
  return
end subroutine phdensitytablejp
! Copyright A.J. Koning 2021
