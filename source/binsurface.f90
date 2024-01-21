subroutine binsurface(x1, y1, x2, y2, x3, y3, xl, xu, yl, yu, pi, twopi, epsx, epsy, surfloc, surfbin)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of the fraction of the bidimensional grid
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit  none
  real(sgl) :: epsx        ! x-accuracy (used to distinguish 2 points)
  real(sgl) :: epsy        ! y-accuracy (used to distinguish 2 points)
  real(sgl) :: pi          ! pi
  real(sgl) :: twopi       ! 2.*pi
  real(sgl) :: x1          ! coordinates of intersection points inside the bin
  real(sgl) :: x2          ! coordinates of the 2nd summit of the triangle
  real(sgl) :: x3          ! coordinates of the 3rd summit of the triangle
  real(sgl) :: xl          ! coordinates defining a bin
  real(sgl) :: xsurf(19)   ! arrays containing the coordinates of the points  defining the covered surface in the bin.
  real(sgl) :: xu          ! upper x-bin coordinates
  real(sgl) :: y1          ! variable for final GOE calculation
  real(sgl) :: y2          ! variable for final GOE calculation
  real(sgl) :: y3          ! variable for final GOE calculation
  real(sgl) :: yl          ! lower and upper y-bin coordinates
  real(sgl) :: ysurf(19)   ! arrays containing the coordinates of the points  defining the covered surface in the bin.
  real(sgl) :: yu          ! lower and upper y-bin coordinates
  integer   :: isurf       ! number of points limiting the surface in the bin
  real(sgl) :: surfloc     ! surface to calculate
  real(sgl) :: surfbin     ! LAB DDX area contribution
  real(sgl) :: xr(4)       ! help array for storing segment intersection points  with the bin sides
  real(sgl) :: yr(4)       ! help array for storing segment intersection points  with the bin sides
  integer   :: nr          ! number of radial grid point
  integer   :: ir          ! loop counter
  real(sgl) :: xr1         ! help variable
  real(sgl) :: yr1         ! help variable
  real(sgl) :: xpts(12)    ! help array for storing all intersection points
  real(sgl) :: ypts(12)    ! help array for storing all intersection points
  integer   :: ipts        ! loop counter
  integer   :: npts        ! total number of intersection points
  real(sgl) :: xdiff       ! help variable
  real(sgl) :: ydiff       ! difference in y
  real(sgl) :: xr2         ! help variable
  real(sgl) :: yr2         ! help variable
  real(sgl) :: xr3         ! help variable
  real(sgl) :: yr3         ! help variable
  real(sgl) :: xg          ! center of gravity of the intersection points
  real(sgl) :: yg          ! center of gravity of the intersection points
  real(sgl) :: weight      ! weight of emission energy bin
  real(sgl) :: norm1       ! norm of vector (x-x1,y-y1)
  real(sgl) :: norm2       ! norm of vector (x2-x1,y2-y1)
  real(sgl) :: norm12      ! distance
  real(sgl) :: angsurf(19) ! angular surface
  real(sgl) :: cos12       ! cosine
  real(sgl) :: ang12       ! angle
  real(sgl) :: deter12     ! determinant
  real(sgl) :: ang1        ! angle
  real(sgl) :: ang2        ! product of angular cross section
  logical   :: belongs     ! logical function to test if one value is between two others
  logical   :: intri       ! function to test if (x,y) belongs to the triangle
  logical   :: invect      ! function to test if (x,y) belongs to the segment
!
! Determine the bin summits included in the considered surface
!
  isurf = 0
  if (intri(xl, yl, x1, y1, x2, y2, x3, y3)) then
    isurf = isurf + 1
    xsurf(isurf) = xl
    ysurf(isurf) = yl
  endif
  if (intri(xl, yu, x1, y1, x2, y2, x3, y3)) then
    isurf = isurf + 1
    xsurf(isurf) = xl
    ysurf(isurf) = yu
  endif
  if (intri(xu, yl, x1, y1, x2, y2, x3, y3)) then
    isurf = isurf + 1
    xsurf(isurf) = xu
    ysurf(isurf) = yl
  endif
  if (intri(xu, yu, x1, y1, x2, y2, x3, y3)) then
    isurf = isurf + 1
    xsurf(isurf) = xu
    ysurf(isurf) = yu
  endif
  if (isurf == 4) then
    surfloc = surfbin
    return
  endif
!
! Segment coordinates definition
!
! intri: function to test if (x,y) belongs to the triangle
! invect: function to test if (x,y) belongs to the segment
!
  ipts = 0
  call intersection(x1, y1, x2, y2, xl, yl, xu, yu, xr, yr, nr)
  do ir = 1, nr
    xr1 = xr(ir)
    yr1 = yr(ir)
    if (invect(xr1, yr1, x1, y1, x2, y2)) then
      ipts = ipts + 1
      xpts(ipts) = xr1
      ypts(ipts) = yr1
    endif
  enddo
  call intersection(x2, y2, x3, y3, xl, yl, xu, yu, xr, yr, nr)
  do ir = 1, nr
    xr1 = xr(ir)
    yr1 = yr(ir)
    if (invect(xr1, yr1, x2, y2, x3, y3)) then
      ipts = ipts + 1
      xpts(ipts) = xr1
      ypts(ipts) = yr1
    endif
  enddo
  call intersection(x3, y3, x1, y1, xl, yl, xu, yu, xr, yr, nr)
  do ir = 1, nr
    xr1 = xr(ir)
    yr1 = yr(ir)
    if (invect(xr1, yr1, x3, y3, x1, y1)) then
      ipts = ipts + 1
      xpts(ipts) = xr1
      ypts(ipts) = yr1
    endif
  enddo
!
! Check if the triangle summits are inside that bin
!
  if (belongs(x1, xu, xl) .and. belongs(y1, yu, yl)) then
    ipts = ipts + 1
    xpts(ipts) = x1
    ypts(ipts) = y1
  endif
  if (belongs(x2, xu, xl) .and. belongs(y2, yu, yl)) then
    ipts = ipts + 1
    xpts(ipts) = x2
    ypts(ipts) = y2
  endif
  if (belongs(x3, xu, xl) .and. belongs(y3, yu, yl)) then
    ipts = ipts + 1
    xpts(ipts) = x3
    ypts(ipts) = y3
  endif
  npts = ipts
!
! intersection analysis : eliminate redundancies
!
! eliminate redundant points
!
  if (isurf == 0 .and. npts > 0) then
    isurf = 1
    xsurf(1) = xpts(1)
    ysurf(1) = ypts(1)
  endif
Loop1:  do ipts = 1, npts
    xr1 = xpts(ipts)
    yr1 = ypts(ipts)
    do ir = 1, isurf
      xdiff = abs(xr1 - xsurf(ir))
      ydiff = abs(yr1 - ysurf(ir))
      if ((xdiff <= epsx) .and. (ydiff <= epsy)) cycle Loop1
    enddo
    isurf = isurf + 1
    xsurf(isurf) = xr1
    ysurf(isurf) = yr1
  enddo Loop1
!
!
!
  if (isurf <= 2) then
    surfloc = 0.
    return
  endif
  if (isurf == 3) then
    xr1 = xsurf(1)
    yr1 = ysurf(1)
    xr2 = xsurf(2)
    yr2 = ysurf(2)
    xr3 = xsurf(3)
    yr3 = ysurf(3)
    surfloc = abs(0.5 * ((xr3 - xr2) * (yr1 - yr2) - (xr1 - xr2) * (yr3 - yr2)))
    return
  endif
!
! If we have more than 4 points defining the covered surface we have to order these points using (xg,yg) as origin and
! xsurf(1),ysurf(1) as (1,0) point
!
  xg = 0.
  yg = 0.
  do ipts = 1, isurf
    xg = xg + xsurf(ipts)
    yg = yg + ysurf(ipts)
  enddo
  weight = 1. / float(isurf)
  xg = xg * weight
  yg = yg * weight
  xr1 = xsurf(1) - xg
  yr1 = ysurf(1) - yg
  norm1 = xr1 * xr1 + yr1 * yr1
  angsurf(1) = 0.
  do ipts = 2, isurf
    xr2 = xsurf(ipts) - xg
    yr2 = ysurf(ipts) - yg
    norm2 = xr2 * xr2 + yr2 * yr2
    norm12 = norm1 * norm2
    if (norm12 == 0.) then
      cos12 = 1.
    else
      cos12 = (xr1 * xr2 + yr1 * yr2) / (sqrt(norm12))
    endif
    if (cos12 <=  - 1.0) then
      ang12 = pi
      goto 20
    endif
    if (cos12 >= 1.0) then
      ang12 = 0.
      goto 20
    endif
    ang12 = acos(cos12)
 20     deter12 = xr1 * yr2 - xr2 * yr1
    if (deter12 < 0.) ang12 = twopi - ang12
    angsurf(ipts) = ang12
  enddo
  do ipts = 1, isurf - 1
    do ir = ipts + 1, isurf
      xr1 = xsurf(ipts)
      yr1 = ysurf(ipts)
      ang1 = angsurf(ipts)
      xr2 = xsurf(ir)
      yr2 = ysurf(ir)
      ang2 = angsurf(ir)
      if (ang2 > ang1) cycle
      xsurf(ipts) = xr2
      ysurf(ipts) = yr2
      angsurf(ipts) = ang2
      xsurf(ir) = xr1
      ysurf(ir) = yr1
      angsurf(ir) = ang1
    enddo
  enddo
!
! Determine the covered surface by summing all the sub-triangles
!
  xr1 = xsurf(1)
  yr1 = ysurf(1)
  xr2 = xsurf(isurf)
  yr2 = ysurf(isurf)
  surfloc = abs(0.5 * ((xg - xr2) * (yr1 - yr2) - (xr1 - xr2) * (yg - yr2)))
  do ipts = 1, isurf - 1
    xr1 = xsurf(ipts)
    yr1 = ysurf(ipts)
    xr2 = xsurf(ipts + 1)
    yr2 = ysurf(ipts + 1)
    surfloc = surfloc + abs(0.5 * ((xg - xr2) * (yr1 - yr2) - (xr1 - xr2) * (yg - yr2)))
  enddo
  return
end subroutine binsurface
! Copyright A.J. Koning 2021
