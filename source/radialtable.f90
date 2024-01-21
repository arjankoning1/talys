subroutine radialtable(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Tabulated radial matter densities
!
! Author    : Arjan Koning and Eric Bauge
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
!   numjlm         ! maximum number of radial points
! Variables for OMP
!   alphaomp       ! alpha optical model
!   flagjlm        ! flag for using semi - microscopic JLM OMP
!   radialfile     ! radial matter density file
!   radialmodel    ! model for radial matter densities (JLM OMP only)
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Constants
!   nuc            ! symbol of nucleus
! Variables for files
!   path           ! directory containing files to be read
! Variables for optical model
!   jlmexist       ! flag for existence of tabulated radial matter density
! Variables for JLM
!   radjlm         ! radial points for JLM potential
!   rhojlmn        ! density for neutrons
!   rhojlmp        ! density for protons
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist     ! logical to determine existence
  character(len=6)  :: radchar    ! help variable
  character(len=132):: radfile    ! radial matter density file
  integer           :: A          ! mass number of target nucleus
  integer           :: i          ! counter
  integer           :: ia         ! mass number from abundance table
  integer           :: istat      ! logical for file access
  integer           :: j          ! counter
  integer           :: Nix        ! neutron number index for residual nucleus
  integer           :: nradjlm    ! number of radial points
  integer           :: Z          ! charge number of target nucleus
  integer           :: Zix        ! charge number index for residual nucleus
  real(sgl)         :: an         ! neutron level density parameter
  real(sgl)         :: ap         ! proton level density parameter
  real(sgl)         :: dr         ! integration bin width
  real(sgl)         :: expo       ! help variable
  real(sgl)         :: h          ! help variable
  real(sgl)         :: rn0        ! JLM density
  real(sgl)         :: rn1        ! number of evaporated neutrons from light fragment
  real(sgl)         :: rp0        ! JLM density
  real(sgl)         :: rp1        ! JLM density
!
! ****************** Read radial matter densities **********************
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  if (radialfile(Zix)(1:1) /= ' ') then
    radfile = radialfile(Zix)
  else
    radchar = trim(nuc(Z))//'.rad'
    if (radialmodel == 1) then
      radfile = trim(path)//'optical/jlm/bskg3/'//radchar
    else
      radfile = trim(path)//'optical/jlm/d1m/'//radchar
    endif
    inquire (file = radfile, exist = lexist)
    if ( .not. lexist) return
  endif
  open (unit = 2, file = radfile, status = 'old', iostat = istat)
  if (istat /= 0) call read_error(radfile, istat)
  do
    read(2, '(4x, i4, i5, f7.3)', iostat = istat) ia, nradjlm, h
    if (istat == -1) then
      close (unit = 2)
      return
    endif
    if (istat /= 0) call read_error(radfile, istat)
    if (A == ia) then
      read(2, * )
      nradjlm = nradjlm - 1
      do i = 1, min(nradjlm, numjlm)
        read(2, * ) radjlm(Zix, Nix, i), (rhojlmp(Zix, Nix, i, j), j = 1, 5), (rhojlmn(Zix, Nix, i, j), j = 1, 5)
      enddo
      exit
    else
      do i = 1, nradjlm
        read(2, * )
      enddo
      cycle
    endif
  enddo
  close (unit = 2)
!
! Extrapolate for untabulated values
!
  if (nradjlm < numjlm) then
    do j = 1, 5
      rp0 = rhojlmp(Zix, Nix, nradjlm, j)
      rp1 = rhojlmp(Zix, Nix, nradjlm - 1, j)
      rn0 = rhojlmn(Zix, Nix, nradjlm, j)
      rn1 = rhojlmn(Zix, Nix, nradjlm - 1, j)
      ap = 0.
      if (rp0 * rp1 > 0.) ap = - log(rp0 / rp1) / h
      an = 0.
      if (rn0 * rn1 > 0.) an = - log(rn0 / rn1) / h
      do i = nradjlm + 1, numjlm
         dr = h * (i - nradjlm)
         if (j == 1) radjlm(Zix, Nix, i) = radjlm(Zix, Nix, nradjlm) + dr
         expo = ap * dr
         if (abs(expo) <= 80. .and. expo > 0.) rhojlmp(Zix, Nix, i, j) = rp0 * exp( - expo)
         expo = an * dr
         if (abs(expo) <= 80. .and. expo > 0.) rhojlmn(Zix, Nix, i, j) = rn0 * exp( - expo)
      enddo
    enddo
  endif
!
! Set JLM flags
!
  if (flagjlm) then
    jlmexist(Zix, Nix, 1) = .true.
    jlmexist(Zix, Nix, 2) = .true.
  endif
  if (alphaomp >= 3 .and. alphaomp <= 5) jlmexist(Zix, Nix, 6) = .true.
  return
end subroutine radialtable
! Copyright A.J. Koning 2021
