function intri(x, y, x1, y1, x2, y2, x3, y3)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Test if (x,y) belongs to the triangle defined by the points
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl                 ! single precision kind
!
! *** Declaration of local data
!
  logical   :: intri    ! function to test if (x,y) belongs to the triangle
  real(sgl) :: sideline ! function to indicate if (x,y) is on one side of the segment
  real(sgl) :: sign1    ! help variable in talys.cmb
  real(sgl) :: sign2    ! help variable in talys.cmb
  real(sgl) :: sign3    ! help variable in talys.cmb
  real(sgl) :: signg1   ! help variable
  real(sgl) :: signg2   ! help variable
  real(sgl) :: signg3   ! help variable
  real(sgl) :: x        ! help variable
  real(sgl) :: x1       ! coordinates of intersection points inside the bin
  real(sgl) :: x2       ! coordinates of the 2nd summit of the triangle
  real(sgl) :: x3       ! coordinates of the 3rd summit of the triangle
  real(sgl) :: xg       ! center of gravity of the intersection points
  real(sgl) :: y        ! coordinates of the point to test
  real(sgl) :: y1       ! variable for final GOE calculation
  real(sgl) :: y2       ! variable for final GOE calculation
  real(sgl) :: y3       ! variable for final GOE calculation
  real(sgl) :: yg       ! center of gravity of the intersection points
!
! ********************************* Method *****************************
!
! We test if the center of gravity (xg,yg) of the three summits and the point to test are on the same side of each of the
! three segments defining the triangle
!
  intri = .false.
!
! Center of gravity
!
  xg = (x1 + x2 + x3) / 3.0
  yg = (y1 + y2 + y3) / 3.0
!
! Sideline for each point
!
! first segment (x1,y1) --> (x2,y2)
!
! sideline: function to indicate if (x,y) is on one side of the segment
!
  sign1 = sideline(x, y, x1, y1, x2, y2)
  signg1 = sideline(xg, yg, x1, y1, x2, y2)
!
! second segment (x2,y2) --> (x3,y3)
!
  sign2 = sideline(x, y, x2, y2, x3, y3)
  signg2 = sideline(xg, yg, x2, y2, x3, y3)
!
! third segment (x3,y3) --> (x1,y1)
!
  sign3 = sideline(x, y, x3, y3, x1, y1)
  signg3 = sideline(xg, yg, x3, y3, x1, y1)
!
! Final test
!
  if ((sign1 == signg1) .or. (sign1 == 0.0)) then
    if ((sign2 == signg2) .or. (sign2 == 0.0)) then
      if ((sign3 == signg3) .or. (sign3 == 0.0)) then
        intri = .true.
      endif
    endif
  endif
  return
end function intri
! Copyright A.J. Koning 2021
