function tripathi(zproj, aproj, iz, ia, e)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Semi-empirical reaction cross section of Tripathi et al.
!
! Author    : Arjan Koning (adapted from R.K. Tripathi)
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod, only: & ! Definition of single and double precision variables
               sgl      ! single precision kind
!
! *** Declaration of local data
!
  implicit none
  integer   :: aproj       ! mass number of particle
  integer   :: ia          ! mass number from abundance table
  integer   :: iz          ! charge number of residual nucleus
  integer   :: xat         ! mass number of projectile
  integer   :: xzt         ! charge number of projectile
  integer   :: zproj       ! charge number of particle
  real(sgl) :: bcm         ! relativistic factor
  real(sgl) :: beta        ! density-dependenceterm of the DDM3Y interaction
  real(sgl) :: bigb        ! help variable
  real(sgl) :: bigr        ! help variable
  real(sgl) :: ce          ! amplitude of the density-independent term of the DDM3Y interaction
  real(sgl) :: const       ! constant
  real(sgl) :: const1      ! constant
  real(sgl) :: delta       ! energy shift
  real(sgl) :: dens        ! total level density
  real(sgl) :: e           ! energy
  real(sgl) :: ecm         ! energy in C.M. frame
  real(sgl) :: ecmp        ! energy in C.M. frame of particle
  real(sgl) :: ecmt        ! energy in C.M. frame of target
  real(sgl) :: expo        ! help variable
  real(sgl) :: expo1       ! exponent
  real(sgl) :: fourthird   ! 4/3
  real(sgl) :: gcm         ! help variable
  real(sgl) :: onethird    ! 1/3
  real(sgl) :: pi          ! pi
  real(sgl) :: plab        ! momentum in LAB framce
  real(sgl) :: radius      ! radius function
  real(sgl) :: rela        ! relativistic energy
  real(sgl) :: rp          ! radii of projectile and target
  real(sgl) :: rt          ! radii of projectile and target
  real(sgl) :: sl          ! SL,ST I, ST II transmission coefficients
  real(sgl) :: t1          ! help variable
  real(sgl) :: term1       ! help variable
  real(sgl) :: tripathi    ! function for semi-empirical reaction cross section of  Tripathi et al.
  real(sgl) :: twxsec      ! cross section
  real(sgl) :: vp          ! optical model parameters for protons
  real(sgl) :: vt          ! help variable
  real(sgl) :: x1          ! coordinates of intersection points inside the bin
  real(sgl) :: xabs        ! cross section
  real(sgl) :: xm          ! help variable
  external radius
!
! ******************* Reaction cross section calculation ***************
!
! radius        : radius function
! all other var.: ask Tripathi
!
  pi =      3.14159265358979323
  onethird = 1. / 3.
  fourthird = 4. / 3.
  tripathi = 0.
  rp = radius(real(aproj))
  rt = radius(real(ia))
  vp = fourthird * pi * rp **3
  vt = fourthird * pi * rt **3
  dens = 0.5 * ((aproj / vp) + (ia / vt))
  const = 1.75 * dens / 8.824728e-02
  if (iz < zproj) then
    xzt = zproj
    xat = aproj
    zproj = iz
    aproj = ia
    zproj = xzt
    aproj = xat
  endif
  if (aproj == 1) const = 2.05
  if (zproj == 2 .and. aproj == 4) const1 = 2.77 - ia * 8.0e-03 + (ia * ia) * 1.8e-05
  if (zproj == 3) const = const * onethird
  t1 = 40.
  if (zproj == 0) then
    if (ia >= 11 .and. ia < 40) t1 = 30.
    if (iz == 14) t1 = 35.
    if (iz == 26) t1 = 30.
  endif
  gcm = (aproj * (1. + e / 938.) + ia) / sqrt(aproj **2 + ia **2 + 2. * aproj * (e + 938.) * ia / 938.)
!
! A. Koning: Safety
!
  if (gcm <= 1.) return
  bcm = sqrt(1. - 1. / gcm **2)
  plab = aproj * sqrt(2. * 938. * e + e * e)
  ecmp = gcm * (e + 938.) * aproj - bcm * gcm * plab - aproj * 938.
  ecmt = gcm * 938. * ia - ia * 938.
  rela = ecmp + ecmt
  ecm = rela
  bigr = rp + rt + 1.2 * (aproj **onethird + ia **onethird) / (ecm **onethird)
  bigb = 1.44 * zproj * iz / bigr
  if (zproj == 1 .and. ia > 56) bigb = 0.90 * bigb
  if (aproj > 56 .and. iz == 1) bigb = 0.90 * bigb
  if (aproj == 1 .and. ia == 12) bigb = 3.5 * bigb
  if (aproj == 1) then
    if (ia <= 16 .and. ia >= 13) bigb = (ia / 7.) * bigb
    if (iz == 12) bigb = 1.8 * bigb
    if (iz == 14) bigb = 1.4 * bigb
    if (iz == 20) bigb = 1.3 * bigb
  endif
  if (aproj == 1 .and. ia < 4) bigb = 21.0 * bigb
  if (aproj < 4 .and. ia == 1) bigb = 21.0 * bigb
  if (aproj == 1 .and. ia == 4) bigb = 27.0 * bigb
  if (aproj == 4 .and. ia == 1) bigb = 27.0 * bigb
  if (zproj == 0 .or. iz == 0) bigb = 0.0
  xm = 1.
  if (zproj == 0) then
    if (ia < 200.) then
      x1 = 2.83 - 3.1e-02 * ia + 1.7e-04 * ia * ia
      if (x1 <= 1) x1 = 1.
      sl = 1.0
      if (ia == 12) sl = 1.6
      if (ia < 12) sl = 0.6
      xm = (1 - x1 * exp( - e / (sl * x1)))
    else
      xm = (1. - 0.3 * exp( - (e - 1.) / 15.)) * (1. - exp( - (e - 0.9)))
    endif
  endif
  if (zproj == 2 .and. aproj == 4) const = const1 - 0.8 / (1. + exp((250. - e) / 75.))
  expo = min((e - 20) / 10., 80.)
  if (zproj == 1 .and. aproj == 1) then
    if (ia > 45) t1 = 40. + ia / 3.
    if (ia < 4) t1 = 55.
    const = 2.05 - 0.05 / (1. + exp((250. - e) / 75.))
    if (ia < 4) const = 1.7
    if (iz == 12) then
      t1 = 40.
      const = 2.05 - 3.0 / (1. + exp(expo))
    endif
    if (iz == 14) then
      t1 = 40.
      const = 2.05 - 1.75 / (1. + exp(expo))
    endif
    if (iz == 18) then
      t1 = 40.
      const = 2.05 - 2.0 / (1. + exp(expo))
    endif
    if (iz == 20) then
      t1 = 40.
      expo1 = min((e - 40) / 10., 80.)
      const = 2.05 - 1.0 / (1. + exp(expo1))
    endif
  endif
  if (zproj == 0 .and. aproj == 1) then
    const = 2. * (0.134457 / dens)
    if (ia > 140 .and. ia < 200) const = const - 1.5 * (ia - 2. * iz) / ia
    if (ia < 60) const = const - 1.5 * (ia - 2. * iz) / ia
    if (ia <= 40) const = const + 0.25 / (1. + exp( - (170. - e) / 100.))
    if (iz > 82) const = const - real(iz) / (ia - iz)
    if (iz >= 82) then
      expo1 = min((e - 20) / 20., 80.)
      const = const - 2.0 / (1. + exp(expo))
    endif
    if (iz <= 20 .and. iz >= 10) const = const - 1.0 / (1. + exp(expo))
  endif
  ce = const * (1. - exp( - e / t1)) - 0.292 * exp( - e / 792) * cos(0.229 * e **0.453)
  term1 = (ia * aproj) **onethird / (ia **onethird + aproj **onethird)
  delta = 1.615 * term1 - 0.873 * ce
  delta = delta + 0.140 * term1 / ecm **onethird
  delta = delta + 0.794 * (ia - 2. * iz) * zproj / (ia * aproj)
  delta = - delta
  beta = 1.
  twxsec = 10. * pi * 1.26 * 1.26 * beta * (0.873 * aproj **onethird + 0.873 * ia **onethird - delta) **2
  xabs = twxsec * (1. - bigb / ecm) * xm
  if (xabs < 0) xabs = 0.
  tripathi = xabs
  return
end function tripathi
! Copyright A.J. Koning 2021
