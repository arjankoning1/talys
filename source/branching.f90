subroutine branching
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Best set of branching ratios
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
!   sgl            ! single precision kind
! All global variables
!   numlev         ! maximum number of discrete levels
!   numN           ! maximum number of neutrons from initial compound nucleus
!   numZ           ! maximum number of protons from initial compound nucleus
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   Ninit          ! neutron number of initial compound nucleus
!   Zinit          ! charge number of initial compound nucleus
!   Ztarget        ! charge number of target nucleus
! Variables for files
!   path           ! directory containing files to be read
! Constants
!   nuc            ! symbol of nucleus
! Variables for levels
!   branchlevel    ! level to which branching takes place
!   branchratio    ! gamma - ray branching ratio to level
!   nbranch        ! number of branching levels
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   range_index_error    ! Test if index is out of range
!   range_real_error    ! Test if real variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist        ! logical to determine existence
  character(len=10) :: branchchar    ! part of branching ratio file
  character(len=132):: line          ! input line
  character(len=132):: word(40)      ! words on input line
  character(len=132):: branchfile    ! branching ratio file
  integer           :: ia            ! mass number from abundance table
  integer           :: ilev0         ! counter for level
  integer           :: ilev1         ! counter for level
  integer           :: istat         ! logical for file access
  integer           :: iword         ! word counter
  integer           :: iz            ! charge number of residual nucleus
  integer           :: k             ! designator for particle
  integer           :: nbr           ! number of branching levels
  integer           :: Nix           ! neutron number index for residual nucleus
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: bra           ! gamma-ray branching ratio to level
  real(sgl)         :: sum           ! help variable
!
! ******************** Read best branching ratios **********************
!
  branchchar = '000.branch'
  write(branchchar(1:3), '(i3.3)') Atarget
  branchfile = trim(path)//'levels/branch/'// trim(nuc(Ztarget)) //branchchar
  inquire (file = branchfile, exist = lexist)
  if (lexist) then
    open (unit = 2, file = branchfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(branchfile, istat)
    do
      read(2, '(a80)', iostat = istat) line
      if (istat == -1) exit
      if (istat /= 0) call read_error(branchfile, istat)
      call getkeywords(line, word)
      read(word(2), * , iostat = istat) iz
      if (istat /= 0) call read_error(line, istat)
      read(word(3), * , iostat = istat) ia
      if (istat /= 0) call read_error(line, istat)
      read(word(4), * , iostat = istat) ilev0
      if (istat /= 0) call read_error(line, istat)
      read(word(5), * , iostat = istat) nbr
      if (istat /= 0) call read_error(line, istat)
      Zix = Zinit - iz
      Nix = Ninit - ia + iz
      if (Zix < 0 .or. Zix > numZ .or. Nix < 0 .or. Nix > numN) cycle
      call range_index_error(branchfile, 'Level', ilev0, 0, numlev)
      call range_index_error(branchfile, '# branchings', nbr, 0, numlev)
      if (nbranch(Zix, Nix, ilev0) /= 0) cycle
      iword = 5
      sum = 0.
      do k = 1, nbr
        iword = iword + 1
        read(word(iword), * , iostat = istat) ilev1
        if (istat /= 0) call read_error(line, istat)
        call range_integer_error(line, ilev1, 0, numlev)
        iword = iword + 1
        read(word(iword), * , iostat = istat) bra
        if (istat /= 0) call read_error(line, istat)
        call range_real_error(line, bra, 0., 1.)
        branchlevel(Zix, Nix, ilev0, k) = ilev1
        branchratio(Zix, Nix, ilev0, k) = bra
        sum = sum + bra
      enddo
      if (sum > 0.) then
        do k = 1, nbr
          branchratio(Zix, Nix, ilev0, k) = branchratio(Zix, Nix, ilev0, k) / sum
        enddo
      endif
      nbranch(Zix, Nix, ilev0) = nbr
    enddo
    close (unit = 2)
  endif
  return
end subroutine branching
! Copyright A.J. Koning 2021
