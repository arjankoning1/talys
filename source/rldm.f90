subroutine rldm(iz, ia, il, egs, esp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Saddle point energies, rotating gs energy
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl    ! single precision kind
! Variables for fission parameters
!   x1b    ! parameter for RLDM
!   x1h    ! parameter for RLDM
!   x2b    ! parameter for RLDM
!   x2h    ! parameter for RLDM
!   x3b    ! parameter for RLDM
!   x3h    ! parameter for RLDM
!
! *** Declaration of local data
!
  implicit none
  integer   :: ia                ! mass number from abundance table
  integer   :: il                ! angular momentum
  integer   :: ix                ! help variable
  integer   :: iy                ! help variable
  integer   :: iz                ! charge number of residual nucleus
  real(sgl) :: amass             ! mass number of residual nucleus
  real(sgl) :: b1                ! help variable
  real(sgl) :: b2                ! beta2
  real(sgl) :: bfrldm            ! help variable
  real(sgl) :: bx                ! help variable
  real(sgl) :: by                ! help variable
  real(sgl) :: cx                ! help variable
  real(sgl) :: cy                ! help variable
  real(sgl) :: dx                ! increment
  real(sgl) :: dy                ! increment
  real(sgl) :: egs               ! rotating ground state energy
  real(sgl) :: ero               ! rotational energy
  real(sgl) :: eso               ! surface energy
  real(sgl) :: esp               ! saddle point energy
  real(sgl) :: h1                ! help variable
  real(sgl) :: h2                ! help variable
  real(sgl) :: hf                ! help variable
  real(sgl) :: ll                ! angular momentum
  real(sgl) :: neut              ! neutron number of residual nucleus
  real(sgl) :: paren             ! help variable
  real(sgl) :: x                 ! help variable
  real(sgl) :: y                 ! coordinates of the point to test
  real(sgl) :: zchar             ! charge number of residual nucleus
!
! ******************** Rotating Liquid Drop Model *********************
!
  zchar = real(iz)
  amass = real(ia)
  ll = real(il)
  neut = amass - zchar
  paren = 1. - 1.7826 * ((neut - zchar) / amass) **2
  eso = 17.9439 * paren * amass **0.666666
  x = 0.019655 * zchar * (zchar / amass) / paren
  ero = 34.548 * ll * ll / amass **1.666666
  y = 1.9254 * ll * ll / (paren * amass **2.3333333)
  ix = int(20. * x + 1.)
  cx = real(ix)
  bx = int(20. * x + 1.)
  dx = real(bx) - cx
  if(x - .25)5, 5, 30
    5 by = 10. * y + 1.
  if(by - 9.)15, 15, 10
   10 by = 9.
   15 if(by - 1.)20, 20, 25
   20 by = 1.
   25 iy = int(by)
  cy = iy
  dy = by - cy
  h1 = (x1h(ix + 1, iy) - x1h(ix, iy)) * dx + x1h(ix, iy)
  h2 = (x1h(ix + 1, iy + 1) - x1h(ix, iy + 1)) * dx + x1h(ix, iy + 1)
  hf = (h2 - h1) * dy + h1
  b2 = (x1b(ix + 1, iy + 1 ) - x1b(ix, iy + 1)) * dx + x1b(ix, iy + 1)
  b1 = (x1b(ix + 1, iy) - x1b(ix, iy)) * dx + x1b(ix, iy)
  bfrldm = (b2 - b1) * dy + b1
  goto 95
   30 if(x - .5)35, 35, 60
   35 by = 20. * y + 1.
  if(by - 10.)45, 45, 40
   40 by = 10.
   45 if(by - 1.)50, 50, 55
   50 by = 1.
   55 ix = ix - 5
  iy = int(by)
  cy = iy
  dy = by - cy
  h1 = (x2h(ix + 1, iy) - x2h(ix, iy)) * dx + x2h(ix, iy)
  h2 = (x2h(ix + 1, iy + 1) - x2h(ix, iy + 1)) * dx + x2h(ix, iy + 1)
  hf = (h2 - h1) * dy + h1
  b1 = (x2b(ix + 1, iy) - x2b(ix, iy)) * dx + x2b(ix, iy)
  b2 = (x2b(ix + 1, iy + 1 ) - x2b(ix, iy + 1)) * dx + x2b(ix, iy + 1)
  bfrldm = (b2 - b1) * dy + b1
  goto 95
   60 if(x - .95)70, 70, 65
   65 x = .95
   70 ix = int(20. * x + 1.)
  ix = ix - 10
  by = 100. * y + 1.
  if(by - 19.)80, 80, 75
   75 by = 19.
   80 if(by - 1.)85, 85, 90
   85 by = 1.
   90 iy = int(by)
  cy = iy
  dy = by - cy
  ix = min(ix, 9)
  iy = min(iy, 19)
  h1 = (x3h(ix + 1, iy) - x3h(ix, iy)) * dx + x3h(ix, iy)
  h2 = (x3h(ix + 1, iy + 1) - x3h(ix, iy + 1)) * dx + x3h(ix, iy + 1)
  hf = (h2 - h1) * dy + h1
  b1 = (x3b(ix + 1, iy) - x3b(ix, iy)) * dx + x3b(ix, iy)
  b2 = (x3b(ix + 1, iy + 1 ) - x3b(ix, iy + 1)) * dx + x3b(ix, iy + 1)
  bfrldm = (b2 - b1) * dy + b1
   95 egs = ero + hf * eso
  esp = egs + bfrldm * eso
  return
end subroutine rldm
! Copyright A.J. Koning 2021
