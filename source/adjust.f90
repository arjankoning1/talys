subroutine adjust(E, key, Zix, Nix, type, ibar, factor)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Energy-dependent parameter adjustment
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
! Definition of single and double precision variables
!   sgl           ! single precision kind
! All global variables
!   numenadj      ! maximum number of energies for adjustable parameters
! Variables for adjustment
!   adjustfile    ! file for local adjustment
!   adjustix      ! local adjustment index
!   adjustkey     ! keyword for local adjustment
!   adjustpar     ! local adjustment parameters
!   Dadjust       ! tabulated depth of local adjustment
!   Eadjust       ! tabulated energy of local adjustment
!   Nadjust       ! number of adjustable parameters
!   nenadjust     ! number of tabulated energies of local adjustment
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key                 ! keyword
  integer           :: i                   ! counter
  integer           :: ibar                ! fission barrier
  integer           :: j                   ! counter
  integer           :: N                   ! neutron number of residual nucleus
  integer           :: nen                 ! energy counter
  integer           :: Nix                 ! neutron number index for residual nucleus
  integer           :: npow                ! power of distribution
  integer           :: type                ! particle type
  integer           :: Zix                 ! charge number index for residual nucleus
  real(sgl)         :: C                   ! integration constant
  real(sgl)         :: D                   ! depth of local adjustment
  real(sgl)         :: Da                  ! depth of local adjustment
  real(sgl)         :: Dadj(0:numenadj)    ! depth of local adjustment
  real(sgl)         :: Db                  ! depth of local adjustment
  real(sgl)         :: E                   ! incident energy
  real(sgl)         :: Ea                  ! start energy of local adjustment
  real(sgl)         :: Eadj(0:numenadj)    ! energy of local adjustment
  real(sgl)         :: Eb                  ! end energy of local adjustment
  real(sgl)         :: Em                  ! intermediate energy of local adjustment
  real(sgl)         :: factor              ! multiplication factor
!
! ************************* OMP Wv function ****************************
!
  factor = 1.
  Eadj(0) = 0.
  Dadj(0) = 1.
  do i = 1, Nadjust
    if (trim(key) == trim(adjustkey(i)) .and. Zix == adjustix(i, 1) .and. Nix == adjustix(i, 2) .and. &
      type == adjustix(i, 3) .and. ibar == adjustix(i, 4)) then
      if (adjustfile(i)(1:1) /= ' ') then
        N = nenadjust(i)
        do j = 1, nenadjust(i)
          Eadj(j) = Eadjust(i, j)
          Dadj(j) = Dadjust(i, j)
        enddo
        if (E >= Eadj(1) .and. E <= Eadj(N)) then
          call locate(Eadj, 1, N, E, nen)
          Ea = Eadj(nen)
          Eb = Eadj(nen + 1)
          Da = Dadj(nen)
          Db = Dadj(nen + 1)
          call pol1(Ea, Eb, Da, Db, E, factor)
        endif
      else
        Ea = adjustpar(i, 1)
        Eb = adjustpar(i, 2)
        Em = adjustpar(i, 3)
        D = adjustpar(i, 4)
        if (E <= Ea .or. E >= Eb) cycle
        npow = 4
        Db = (D - 1.) * (1. + 0.5 **npow)
        if (E <= Em) then
          C = Ea + 0.5 * (Em - Ea)
          factor = 1. + Db * ((E - Ea) **npow) / ((E - Ea) **npow + (C - Ea) **npow)
        else
          C = Em + 0.5 * (Eb - Em)
          factor = 1. + Db * ((Eb - E) **npow) / ((Eb - E) **npow + (C - Eb) **npow)
        endif
      endif
    endif
  enddo
  return
end subroutine adjust
! Copyright A.J. Koning 2021
