subroutine labsurface(Zcomp, Ncomp, type, ityp, x1, cos1, sin1, x2, cos2, sin2, x3, cos3, sin3, totsurf, limx, limy, sparte, spartr)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Calculation of the way the surface triangle defined by its summits (x1,y1),(x2,y2) and (x3,y3)
! |         given bidimensional grid depending on ityp.
!
! Author    : Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl              ! single precision kind
! All global variables
!   numangcont       ! maximum number of angles for continuum
!   numangrec        ! maximum number of recoil angles
!   numen2           ! maximum number of outgoing energies
!   numenrec         ! maximum number of recoil energies
! Variables for numerics
!   maxenrec         ! number of recoil energies
!   nanglecont       ! number of angles for continuum
!   nanglerec        ! number of recoil angles
! Variables for nuclides
!   Nindex           ! neutron number index for residual nucleus
!   Zindex           ! charge number index for residual nucleus
! Constants
!   pi               ! pi
!   sgn              ! sign
!   twopi            ! 2 * pi
! Variables for recoil initialization
!   cosangcontmax    ! cosine of maximum of angular bin
!   cosangcontmin    ! cosine of minimum of angular bin
!   cosrecmax        ! array for cosine of upper bin angle values
!   cosrecmin        ! array for cosine of lower bin angle values
! Variables for recoil
!   angcm            ! CM angle with respect to the LAB
!   areaejlab        ! Total surface of LAB ddx bins
!   areareclab       ! Total surface of LAB ddx bins
!   Eejlabmax        ! maximum energy of ejectile lab bin
!   Eejlabmin        ! minimum energy of ejectile lab bin
!   Erecmax          ! minimal energy limit of recoil bin
!   Erecmin          ! minimal energy limit of recoil bin
!   iejlab           ! number of ejectile lab bins
!
! *** Declaration of local data
!
  implicit none
  integer   :: ityp                                ! =0 means ejectile, =1 means recoil
  integer   :: Ncomp                               ! neutron number index for compound nucleus
  integer   :: Nix                                 ! neutron number index for residual nucleus
  integer   :: type                                ! particle type
  integer   :: Zcomp                               ! proton number index for compound nucleus
  integer   :: Zix                                 ! charge number index for residual nucleus
  real(sgl) :: cos1                                ! cos of angle
  real(sgl) :: cos2                                ! cos of angle
  real(sgl) :: cos3                                ! cos of angle
  real(sgl) :: sin1                                ! sin of angle
  real(sgl) :: sin2                                ! sin of angle
  real(sgl) :: sin3                                ! sin of angle
  real(sgl) :: totsurf                             ! triangle surface
  real(sgl) :: x1                                  ! coordinates of intersection points inside the bin
  real(sgl) :: x2                                  ! coordinates of the 2nd summit of the triangle
  real(sgl) :: x3                                  ! coordinates of the 3rd summit of the triangle
  integer   :: limx(2)                             ! x limits covered by the triangle
  integer   :: limy(2)                             ! y limits covered by the triangle
  real(sgl) :: sparte(numen2, 0:2*numangcont+1)    ! array containing the partial surface (filled if ityp=0)
  real(sgl) :: spartr(0:numenrec, 0:2*numangrec+1) ! array containing the partial surface (filled if ityp=1)
  logical   :: belongs                             ! logical function to test if one value is between two others
  integer   :: recoilmode                          ! type of recoil calculation (1=most exact,2=simplest)
  real(sgl) :: angle1                              ! angle
  real(sgl) :: angle2                              ! angle
  real(sgl) :: angle3                              ! angle
  real(sgl) :: yy1                                 ! help variable
  real(sgl) :: yy2                                 ! help variable
  real(sgl) :: yy3                                 ! help variable
  real(sgl) :: yc1                                 ! cosine of yy1
  real(sgl) :: yc2                                 ! cosine of yy2
  real(sgl) :: yc3                                 ! cosine of yy3
  integer   :: is1                                 ! help variable
  integer   :: is2                                 ! help variable
  integer   :: is3                                 ! help variable
  integer   :: nt1                                 ! help variable
  integer   :: nt2                                 ! help variable
  integer   :: nt3                                 ! help variable
  real(sgl) :: add1                                ! help variable
  real(sgl) :: add2                                ! help variable
  real(sgl) :: add3                                ! help variable
  integer   :: ixmax                               ! maximum x-loop index
  integer   :: ix                                  ! help variable
  real(sgl) :: xl                                  ! coordinates defining a bin
  real(sgl) :: xu                                  ! upper x-bin coordinates
  real(sgl) :: yl                                  ! lower and upper y-bin coordinates
  real(sgl) :: yu                                  ! lower and upper y-bin coordinates
  real(sgl) :: xl1                                 ! 1st limit
  real(sgl) :: xu1                                 ! upper x-bin coordinates
  real(sgl) :: xl2                                 ! 2nd limit
  real(sgl) :: xu2                                 ! upper x-bin coordinates
  real(sgl) :: xl3                                 ! lower x-bin coordinates
  real(sgl) :: xu3                                 ! upper x-bin coordinates
  integer   :: ix1                                 ! index locating the summits in the x-axis grid
  integer   :: ix2                                 ! index locating the summits in the x-axis grid
  integer   :: ix3                                 ! index locating the summits in the x-axis grid
  integer   :: iymax                               ! maximum y-loop index
  integer   :: iy                                  ! help variable
  real(sgl) :: yl1                                 ! lower y-bin coordinates
  real(sgl) :: yu1                                 ! upper y-bin coordinates
  real(sgl) :: yl2                                 ! lower y-bin coordinates
  real(sgl) :: yu2                                 ! upper y-bin coordinates
  real(sgl) :: yl3                                 ! lower y-bin coordinates
  real(sgl) :: yu3                                 ! upper y-bin coordinates
  integer   :: iy1                                 ! index locating the summits in the y-axis grid
  integer   :: iy2                                 ! index locating the summits in the y-axis grid
  integer   :: iy3                                 ! index locating the summits in the y-axis grid
  integer   :: ixmin                               ! minimum x-loop index
  integer   :: iymin                               ! minimum y-loop index
  integer   :: iysup                               ! maximum y index
  integer   :: numbinscovered                      ! total number of covered bins
  real(sgl) :: oneovernumbins                      ! help variable to improve speed
  integer   :: iystore                             ! help variable
  real(sgl) :: surfloc                             ! surface to calculate
  real(sgl) :: sumsurf                             ! sum of partial surfaces
  real(sgl) :: epsx                                ! x-accuracy (used to distinguish 2 points)
  integer   :: nturn                               ! general case
  real(sgl) :: epsy                                ! y-accuracy (used to distinguish 2 points)
  real(sgl) :: surfbin                             ! LAB DDX area contribution
  real(sgl) :: renorm                              ! renormalisation factor
!
! ************************** Initializations ***************************
!
! Basic initialisations
!
  recoilmode = 1
  Zix = Zindex(Zcomp, Ncomp, type)
  Nix = Nindex(Zcomp, Ncomp, type)
!
! We here calculate the lab angle resulting from the coupling of the angle of the considered CM with respect to the beam axis
! (angcm) with the emission angles in this CM.
!
  ix1 = 1
  ix2 = 1
  ix3 = 1
  iy1 = 0
  iy2 = 0
  iy3 = 0
  cos1 = max(cos1, - 1.)
  cos1 = min(cos1, 1.)
  angle1 = acos(cos1)
  if (sin1 < 0.) angle1 = twopi - angle1
  cos2 = max(cos2, - 1.)
  cos2 = min(cos2, 1.)
  angle2 = acos(cos2)
  if (sin2 < 0.) angle2 = twopi - angle2
  cos3 = max(cos3, - 1.)
  cos3 = min(cos3, 1.)
  angle3 = acos(cos3)
  if (sin3 < 0.) angle3 = twopi - angle3
  yy1 = angle1 + angcm
  yy2 = angle2 + angcm
  yy3 = angle3 + angcm
  yc1 = cos(yy1)
  yc2 = cos(yy2)
  yc3 = cos(yy3)
  is1 = int(yy1 / pi)
  is2 = int(yy2 / pi)
  is3 = int(yy3 / pi)
  nt1 = int(yy1 / twopi)
  nt2 = int(yy2 / twopi)
  nt3 = int(yy3 / twopi)
  add1 = 2. * is1
  add2 = 2. * is2
  add3 = 2. * is3
  yc1 = - add1 + sgn(is1) * yc1
  yc2 = - add2 + sgn(is2) * yc2
  yc3 = - add3 + sgn(is3) * yc3
!
! Reference area calculated
!
  totsurf = 0.5 * ((x1 - x3) * (yc2 - yc3) - (x2 - x3) * (yc1 - yc3))
  totsurf = abs(totsurf)
!
! if the lab area is weak, we use the simplest approximation
!
  if (totsurf < 1.e-14) then
    totsurf = 1.e-14
    recoilmode = 2
  endif
!
! ****************** Ejectile surface calculation **********************
!
  if (ityp /= 0) then
!
! determine x and y loop limits
!
! xloop
!
    ixmax = iejlab(type)
    if (x1 > Eejlabmax(type, ixmax)) then
      Eejlabmax(type, ixmax) = x1
      ix1 = ixmax
    else
      do ix = 1, ixmax
        xl1 = Eejlabmin(type, ix)
        xu1 = Eejlabmax(type, ix)
        if (belongs(x1, xl1, xu1)) then
          ix1 = ix
          exit
        endif
      enddo
    endif
    if (x2 > Eejlabmax(type, ixmax)) then
      Eejlabmax(type, ixmax) = x2
      ix2 = ixmax
    else
      do ix = 1, ixmax
        xl2 = Eejlabmin(type, ix)
        xu2 = Eejlabmax(type, ix)
        if (belongs(x2, xl2, xu2)) then
          ix2 = ix
          exit
        endif
      enddo
    endif
    if (x3 > Eejlabmax(type, ixmax)) then
      Eejlabmax(type, ixmax) = x3
      ix3 = ixmax
     else
      do ix = 1, ixmax
        xl3 = Eejlabmin(type, ix)
        xu3 = Eejlabmax(type, ix)
        if (belongs(x3, xl3, xu3)) then
          ix3 = ix
          exit
        endif
      enddo
    endif
!
! yloop
!
    iymax = 2 * nanglecont + 1
    do iy = 0, iymax
      yl1 = - 4. * nt1
      yu1 = - 4. * nt1
      if (iy <= nanglecont) then
          yl1 = yl1 + cosangcontmin(iy)
          yu1 = yu1 + cosangcontmax(iy)
        else
          yl1 = yl1 - cosangcontmin(iy) - 2.
          yu1 = yu1 - cosangcontmax(iy) - 2.
      endif
      if (belongs(yc1, yl1, yu1)) then
        iy1 = iy + nt1 * (iymax + 1)
        exit
      endif
    enddo
    do iy = 0, iymax
      yl2 = - 4. * nt2
      yu2 = - 4. * nt2
      if (iy <= nanglecont) then
          yl2 = yl2 + cosangcontmin(iy)
          yu2 = yu2 + cosangcontmax(iy)
        else
          yl2 = yl2 - cosangcontmin(iy) - 2.
          yu2 = yu2 - cosangcontmax(iy) - 2.
      endif
      if (belongs(yc2, yl2, yu2)) then
        iy2 = iy + nt2 * (iymax + 1)
        exit
      endif
    enddo
    do iy = 0, iymax
      yl3 = - 4. * nt3
      yu3 = - 4. * nt3
      if (iy <= nanglecont) then
          yl3 = yl3 + cosangcontmin(iy)
          yu3 = yu3 + cosangcontmax(iy)
        else
          yl3 = yl3 - cosangcontmin(iy) - 2.
          yu3 = yu3 - cosangcontmax(iy) - 2.
      endif
      if (belongs(yc3, yl3, yu3)) then
        iy3 = iy + nt3 * (iymax + 1)
        exit
      endif
    enddo
    ixmin = min(ix1, ix2)
    ixmin = min(ixmin, ix3)
    ixmax = max(ix1, ix2)
    ixmax = max(ixmax, ix3)
    iymin = min(iy1, iy2)
    iymin = min(iymin, iy3)
    iysup = max(iy1, iy2)
    iysup = max(iysup, iy3)
    limx(1) = ixmin
    limx(2) = ixmax
    limy(1) = iymin
    limy(2) = iysup
!
! SPECIAL CASES
!
! if the three summits are inside the same bin
!
    if ((ixmin == ixmax) .and. (iymin == iysup)) then
      iymin = mod(iymin, iymax + 1)
      sparte(max(ixmin, 1), iymin) = totsurf
      return
    endif
!
! simplest approximation
!
    if (recoilmode == 2) then
      numbinscovered = (ixmax - ixmin + 1) * (iysup - iymin + 1)
      oneovernumbins = 1. / float(numbinscovered)
      surfloc = totsurf * oneovernumbins
      do ix = ixmin, ixmax
        do iy = iymin, iysup
          iystore = mod(iy, iymax + 1)
          sparte(max(ix, 1), iystore) = surfloc
        enddo
      enddo
      return
    endif
!
! GENERAL CASE
!
    sumsurf = 0.
    do ix = ixmin, ixmax
      xl = Eejlabmin(type, ix)
      xu = Eejlabmax(type, ix)
      epsx = (xu - xl) / 100000.
      do iy = iymin, iysup
        iystore = mod(iy, iymax + 1)
        nturn = iy / (2 * nanglecont + 2)
        yl = - 4. * nturn
        yu = - 4. * nturn
        if (iystore <= nanglecont) then
          yl = yl + cosangcontmin(iystore)
          yu = yu + cosangcontmax(iystore)
        else
          yl = yl - cosangcontmin(iystore) - 2.
          yu = yu - cosangcontmax(iystore) - 2.
        endif
        epsy = abs(yu - yl) / 100000.
        surfbin = areaejlab(type, ix, iystore)
        call binsurface(x1, yc1, x2, yc2, x3, yc3, xl, xu, yl, yu, pi, twopi, epsx, epsy, surfloc, surfbin)
        sparte(max(ix, 1), iystore) = surfloc
        sumsurf = sumsurf + surfloc
      enddo
    enddo
!
! Renormalisation of ejectiles
!
    if (totsurf /= sumsurf) then
      if (sumsurf /= 0.) then
        renorm = totsurf / sumsurf
        sumsurf = 0.
        do ix = ixmin, ixmax
          do iy = iymin, iysup
            iystore = mod(iy, iymax + 1)
            sparte(max(ix, 1), iystore) = sparte(max(ix, 1), iystore) * renorm
            sumsurf = sumsurf + sparte(max(ix, 1), iystore)
          enddo
        enddo
      else
        numbinscovered = (ixmax - ixmin + 1) * (iysup - iymin + 1)
        oneovernumbins = 1. / float(numbinscovered)
        surfloc = totsurf * oneovernumbins
        sumsurf = 0.
        do ix = ixmin, ixmax
          do iy = iymin, iysup
            iystore = mod(iy, iymax + 1)
            sparte(max(ix, 1), iystore) = surfloc
            sumsurf = sumsurf + surfloc
          enddo
        enddo
      endif
    endif
    return
  else
!
! ********************* Recoil surface calculation *********************
!
! determine x and y loop limits
!
! xloop
!
    if (x1 > Erecmax(Zix, Nix, maxenrec)) then
      Erecmax(Zix, Nix, maxenrec) = x1
      ix1 = maxenrec
    else
      do ix = 0, maxenrec
        xl1 = Erecmin(Zix, Nix, ix)
        xu1 = Erecmax(Zix, Nix, ix)
        if (belongs(x1, xl1, xu1)) then
          ix1 = ix
          exit
        endif
      enddo
    endif
    if (x2 > Erecmax(Zix, Nix, maxenrec)) then
      Erecmax(Zix, Nix, maxenrec) = x2
      ix2 = maxenrec
    else
      do ix = 0, maxenrec
        xl2 = Erecmin(Zix, Nix, ix)
        xu2 = Erecmax(Zix, Nix, ix)
        if (belongs(x2, xl2, xu2)) then
          ix2 = ix
          exit
        endif
      enddo
    endif
    if (x3 > Erecmax(Zix, Nix, maxenrec)) then
      Erecmax(Zix, Nix, maxenrec) = x3
      ix3 = maxenrec
    else
      do ix = 0, maxenrec
        xl3 = Erecmin(Zix, Nix, ix)
        xu3 = Erecmax(Zix, Nix, ix)
        if (belongs(x3, xl3, xu3)) then
          ix3 = ix
          exit
        endif
      enddo
    endif
!
! yloop
!
    iymax = 2 * nanglerec + 1
    do iy = 0, iymax
      yl1 = - 4. * nt1
      yu1 = - 4. * nt1
      if (iy <= nanglerec) then
          yl1 = yl1 + cosrecmin(iy)
          yu1 = yu1 + cosrecmax(iy)
        else
          yl1 = yl1 - cosrecmin(iy) - 2.
          yu1 = yu1 - cosrecmax(iy) - 2.
      endif
      if (belongs(yc1, yl1, yu1)) then
        iy1 = iy + nt1 * (iymax + 1)
        exit
      endif
    enddo
    do iy = 0, iymax
      yl2 = - 4. * nt2
      yu2 = - 4. * nt2
      if (iy <= nanglerec) then
          yl2 = yl2 + cosrecmin(iy)
          yu2 = yu2 + cosrecmax(iy)
        else
          yl2 = yl2 - cosrecmin(iy) - 2.
          yu2 = yu2 - cosrecmax(iy) - 2.
      endif
      if (belongs(yc2, yl2, yu2)) then
        iy2 = iy + nt2 * (iymax + 1)
        exit
      endif
    enddo
    do iy = 0, iymax
      yl3 = - 4. * nt3
      yu3 = - 4. * nt3
      if (iy <= nanglerec) then
          yl3 = yl3 + cosrecmin(iy)
          yu3 = yu3 + cosrecmax(iy)
        else
          yl3 = yl3 - cosrecmin(iy) - 2.
          yu3 = yu3 - cosrecmax(iy) - 2.
      endif
      if (belongs(yc3, yl3, yu3)) then
        iy3 = iy + nt3 * (iymax + 1)
        exit
      endif
    enddo
    ixmin = min(ix1, ix2)
    ixmin = min(ixmin, ix3)
    ixmax = max(ix1, ix2)
    ixmax = max(ixmax, ix3)
    iymin = min(iy1, iy2)
    iymin = min(iymin, iy3)
    iysup = max(iy1, iy2)
    iysup = max(iysup, iy3)
    limx(1) = ixmin
    limx(2) = ixmax
    limy(1) = iymin
    limy(2) = iysup
!
! SPECIAL CASES
!
! if the three summits are inside the same bin
!
    if ((ixmin == ixmax) .and. (iymin == iysup)) then
      iymin = mod(iymin, iymax + 1)
      spartr(ixmin, iymin) = totsurf
      return
    endif
!
! simplest approximation
!
    if (recoilmode == 2) then
      numbinscovered = (ixmax - ixmin + 1) * (iysup - iymin + 1)
      oneovernumbins = 1. / float(numbinscovered)
      surfloc = totsurf * oneovernumbins
      do ix = ixmin, ixmax
        do iy = iymin, iysup
          iystore = mod(iy, iymax + 1)
          spartr(ix, iystore) = surfloc
        enddo
      enddo
      return
    endif
!
! GENERAL CASE
!
    sumsurf = 0.
    do ix = ixmin, ixmax
      xl = Erecmin(Zix, Nix, ix)
      xu = Erecmax(Zix, Nix, ix)
      epsx = (xu - xl) / 100000.
      do iy = iymin, iysup
        iystore = mod(iy, iymax + 1)
        nturn = iy / (2 * nanglerec + 2)
        yl = - 4. * nturn
        yu = - 4. * nturn
        if (iystore <= nanglerec) then
          yl = yl + cosrecmin(iystore)
          yu = yu + cosrecmax(iystore)
        else
          yl = yl - cosrecmin(iystore) - 2.
          yu = yu - cosrecmax(iystore) - 2.
        endif
        epsy = abs(yu - yl) / 100000.
        surfbin = areareclab(Zix, Nix, ix, iystore)
        call binsurface(x1, yc1, x2, yc2, x3, yc3, xl, xu, yl, yu, pi, twopi, epsx, epsy, surfloc, surfbin)
        spartr(ix, iystore) = surfloc
        sumsurf = sumsurf + surfloc
      enddo
    enddo
!
! Renormalisation of recoils
!
    if (totsurf /= sumsurf) then
      if (sumsurf /= 0.) then
        renorm = totsurf / sumsurf
        sumsurf = 0.
        do ix = ixmin, ixmax
          do iy = iymin, iysup
            iystore = mod(iy, iymax + 1)
            spartr(ix, iystore) = spartr(ix, iystore) * renorm
            sumsurf = sumsurf + spartr(ix, iystore)
          enddo
        enddo
      else
        numbinscovered = (ixmax - ixmin + 1) * (iysup - iymin + 1)
        oneovernumbins = 1. / float(numbinscovered)
        surfloc = totsurf * oneovernumbins
        sumsurf = 0.
        do ix = ixmin, ixmax
          do iy = iymin, iysup
            iystore = mod(iy, iymax + 1)
            spartr(ix, iystore) = surfloc
            sumsurf = sumsurf + surfloc
          enddo
        enddo
      endif
    endif
  endif
  return
end subroutine labsurface
! Copyright A.J. Koning 2021
