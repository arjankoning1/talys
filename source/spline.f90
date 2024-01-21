subroutine spline(x, y, nk, yp1, ypn, y2)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Spline fit
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl               ! single precision kind
!
! *** Declaration of local data
!
  integer   :: nk     ! help variable
  real(sgl) :: x(nk)  ! help variable
  real(sgl) :: y(nk)  ! coordinates of the point to test
  real(sgl) :: y2(nk) ! variable for final GOE calculation
  real(sgl) :: yp1    ! y value
  real(sgl) :: ypn    ! y value
  integer   :: i      ! level
  integer   :: k      ! designator for particle
  real(sgl) :: psp    ! help variable
  real(sgl) :: qn     ! help variable
  real(sgl) :: sig    ! help variable
  real(sgl) :: un     ! help variable
  real(sgl) :: u(500) ! double folding potential
!
! **********************************************************************
!
  if (yp1 > 0.99e30) then
    y2(1) = 0.
    u(1) = 0.
  else
    y2(1) = - 0.5
    u(1) = (3. / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
  endif
  do i = 2, nk - 1
    sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
    psp = sig * y2(i - 1) + 2.
    y2(i) = (sig - 1.) / psp
    u(i) = (6. * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - (y(i) - y(i - 1)) / &
      (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / psp
  enddo
  if (ypn > .99e30) then
    qn = 0.
    un = 0.
  else
    qn = 0.5
    un = (3. / (x(nk) - x(nk - 1))) * (ypn - (y(nk) - y(nk - 1)) / (x(nk) - x(nk - 1)))
  endif
  y2(nk) = (un - qn * u(nk - 1)) / (qn * y2(nk - 1) + 1.)
  do k = nk - 1, 1, - 1
    y2(k) = y2(k) * y2(k + 1) + u(k)
  enddo
  return
end subroutine spline
! Copyright A.J. Koning 2021
