function sform(x, y)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Form factor for the Coulomb interaction energy between two
!
! Author    : Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl              ! single precision kind
!
! *** Declaration of local data
!
  integer   :: k     ! designator for particle
  integer   :: m     ! counter
  real(sgl) :: cnm   ! help variable
  real(sgl) :: sform ! function for Form factor for the Coulomb interaction energy
  real(sgl) :: sm    ! help variable
  real(sgl) :: sn    ! help variable
  real(sgl) :: x     ! help variable
  real(sgl) :: y     ! coordinates of the point to test
  real(sgl) :: z1loc ! help variable
  real(sgl) :: z2loc ! help variable
  real(sgl) :: zm    ! help variable
  real(sgl) :: zn    ! help variable
!
! **********************************************************************
!
  if (abs(x) > abs(y)) then
!
! ph. quentin, j. de physique 30 (1969) 497.
!
! sform : function for Form factor for the Coulomb interaction energy between two spheroids
!
    z1loc = x
    z2loc = y
  else
    z1loc = y
    z2loc = x
  endif
  if (z1loc == 0.) then
    sform = 1.
    return
  endif
  z1loc = sign(z1loc **2, z1loc)
  z2loc = sign(z2loc **2, z2loc)
  sform = 0.
  zn = 1. / z1loc
  do k = 0, 20
    zn = zn * z1loc
    sn = 0.75 / ((k + 0.5) * (k + 1.5)) * zn
    zm = 1.
    cnm = 1.
    do m = 1, 15
      zm = zm * z2loc
      cnm = cnm * (k + m) * (k + m - 0.5) / (m * (m - 0.5))
      sm = 0.5625 / ((k + 0.5) * (k + 1.5) * (m + 0.5) * (m + 1.5)) * cnm * zn * zm
      sn = sn + sm
    enddo
    sform = sform + sn
  enddo
  return
end function sform
! Copyright A.J. Koning 2021
