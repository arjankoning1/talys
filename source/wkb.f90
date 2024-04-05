subroutine wkb(Z, A, Zix, Nix, nbar)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of WKB approximation for fission
!
! Author    : Roberto Capote, Arjan Koning, Stephane Goriely, Guillaume Scamps
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
!   numbar        ! number of fission barriers
! Variables for fission
!   bdamp         ! fission partial damping parameter
!   fbarrier      ! height of fission barrier
!   flagfisout    ! flag for output of fission information
!   flagfispartdamp   ! flag for fission partial damping
!   fwidth        ! width of fission barrier
! Variables for WKB
!   betafis       ! fission path width
!   iiextr        ! WKB counter
!   nbinswkb      ! integration step for WKB calculation
!   nextr         ! WKB counter
!   Twkb          ! transmission coefficient of WKB potential
!   Twkbdir       ! transmission coefficient of WKB potential
!   Twkbphase     ! transmission coefficient of WKB potential
!   Twkbtrans     ! transmission coefficient of WKB potential
!   Uwkb          ! energy of WKB potential
!   vfis          ! adjustable factor for fission path height
!   Vheight       ! height of WKB potential
!   Vpos          ! position of WKB potential
!   Vwidth        ! width of WKB potential
!
! *** Declaration of local data
!
  implicit none
  integer          :: A               ! mass number of target nucleus
  integer          :: Find_Extrem     ! function to find extrema
  integer          :: i               ! counter
  integer          :: ii              ! counter
  integer          :: j               ! counter
  integer          :: n1              ! number of coordinate grid points
  integer          :: nbar            ! number of fission barriers
  integer          :: nbarr           ! number of fission barriers
  integer          :: Nix             ! neutron number index for residual nucleus
  integer          :: nsmooth         ! smoothing parameter
  integer          :: Z               ! charge number of target nucleus
  integer          :: Zix             ! charge number index for residual nucleus
  real(sgl)        :: centr           ! center of fitted parabola
  real(sgl)        :: dE              ! help variable
  real(sgl)        :: height          ! barrier height
  real(sgl)        :: phase(2*numbar) ! phase
  real(sgl)        :: rmiu            ! parameter for WKB
  real(sgl)        :: tdir            ! WKB variable
  real(sgl)        :: tff(2*numbar)   ! fission transmission coefficient
  real(sgl)        :: ucentr          ! energy of center
  real(sgl)        :: uexc            ! excitation energy
  real(sgl)        :: uexc0           ! excitation energy
  real(sgl)        :: uexc1           ! excitation energy
  real(sgl)        :: uheight         ! height of parabola
  real(sgl)        :: uwidth          ! height of parabola
  real(sgl)        :: width           ! Full width at half maximum of the breakup nucleon  energy distribution, Kalbach 2003
  real(sgl)        :: Vwell           ! help variable
  real(sgl)        :: Vwell2          ! help variable
  real(sgl)        :: Vtop            ! help variable
  real(sgl)        :: Vtop2           ! help variable
  real(sgl)        :: b1              ! help variable
  real(sgl)        :: b2              ! help variable
  character(len=9) :: filewkb         ! WKB file
!
! ******* Finding maxima and minima of the deformation energy curve ****
!
! Find_Extrem: function to find extrema
!
  nsmooth = 5
  rmiu = 0.054 * A **(5. / 3.)
! initializes iiextr() and nextr
  nextr = Find_Extrem(nsmooth)
  nextr = min(nextr, 2 * numbar - 1)
  nbarr = nextr / 2 + 1
  nbar = nbarr
!   nwell = nextr/2
!   Fitting parabola
  filewkb = '         '
  write(filewkb, '("wkb", 2i3.3)') Z, A
  if (flagfisout) then
    open (unit = 22, file = filewkb, status = 'replace')
    write(22, '("Z=", i3, " A=", i3)') Z, A
  endif
  do j = 1, 2 * numbar
    Vheight(j) = 0.
    Vwidth(j) = 0.
  enddo
  do j = 1, nextr
    CALL ParabFit(iiextr(j), 3, rmiu, betafis, vfis, centr,  height,  width, ucentr, uheight, uwidth)
    if(width.LT.0.05d0) CYCLE ! Skipping very narrow peaks
    if (flagfisout) then
      write(22, *) ' Def: ', betafis(iiextr(j)), ' (', centr, ' +/- ', ucentr, ')'
      write(22, *) ' Height :', vfis(iiextr(j)), ' (', height, ' +/- ', uheight, ')'
      write(22, *) ' Width :', width, ' +/- ', uwidth
      write(22, *) '*******************************************'
    endif
!--------------------------------------------------
!   Initializes parabola's parameters
!
!   The real height of the barrier is used (not the fitted
!   parabola height)
!
    Vheight(j) = vfis(iiextr(j))
    Vwidth(j) = width
    if (mod(j, 2) == 1) then
      ii = (j + 1) / 2
      fbarrier(Zix, Nix, ii) = vfis(iiextr(j))
      fwidth(Zix, Nix, ii) = width
    endif
!   The real position of the barrier is used (not the fitted
!   parabola centr)
    Vpos(j) = betafis(iiextr(j))
!--------------------------------------------------
  enddo
  Vpos(nextr + 1) = 100.d0
!
!   PAUSE
!==================================================================
!   Calculating the limit for the energy loop
!   (Depends on barriers' height and depth)
!
!   ustep = 0.01d0
  uexc0 = 1.d0
  uexc1 = 7.d0
  uexc0 = Vheight(2)
  if(nextr == 5) uexc0 = min(uexc0, Vheight(4))
  uexc1 =     max(Vheight(1), Vheight(3))
  if(nextr == 5) uexc1 =  max(uexc1, Vheight(5))
!==================================================================
  if (flagfisout) write(22, '(1x, A5, 2x, 18(A10, 1x))')  ' Uexc', '   Tdir   ', &
 &  'TaTb/Ta+Tb', '    Ta    ', '    Tb    ', '    Tc    '
!
!   ENERGY LOOP
!
  b1 = bdamp(Zix, Nix, 1)
  b2 = bdamp(Zix, Nix, 2)
  n1 = 3 * nbinswkb / 4
  uexc = 0.
  Uwkb(Zix, Nix, 0) = 0.
  do j = 1, nbar
    Twkb(Zix, Nix, 0, j) = 0.
    Twkbdir(Zix, Nix, 0, j) = 0.
    Twkbtrans(Zix, Nix, 0, j) = 0.
    Twkbphase(Zix, Nix, 0, j) = 0.
  enddo
  do i = 1, nbinswkb
    if (i <= n1) then
      dE = uexc1 / n1
    else
      dE = 0.2
    endif
    uexc = uexc + dE
    CALL WKBFIS(Z, A - Z, uexc, tff, phase, tdir)
    Uwkb(Zix, Nix, i) = uexc
    if (flagfispartdamp) then
      if (nbar > 1) then
        Vwell = 0.
        Vwell2 = 0.
        if (iiextr(2) > 0) Vwell = vfis(iiextr(2))
        if (iiextr(4) > 0) Vwell2 = vfis(iiextr(4))
        Vtop = min(Vheight(1), Vheight(3))
        Vtop2 = min(Vheight(3), Vheight(5))
        Twkbdir(Zix, Nix, i, 1) = tdir
        if (uexc < Vwell) then
          Twkbtrans(Zix, Nix, i, 1) = 0.
        else if (uexc > Vtop) then
          Twkbtrans(Zix, Nix, i, 1) = 1.
        else
          Twkbtrans(Zix, Nix, i, 1) = (uexc**2 - Vwell**2) / ( (Vtop**2 - Vwell**2) * exp(-(uexc - Vtop) / b1) )
        endif
        if (uexc < Vwell2) then
          Twkbtrans(Zix, Nix, i, 2)=0.
        else if (uexc > Vtop2) then
          Twkbtrans(Zix, Nix, i, 2) = 1.
        else
          Twkbtrans(Zix, Nix, i, 2) = (uexc**2 - Vwell2**2)/ ( (Vtop2**2 -Vwell2**2) * exp(-(uexc - Vtop2) / b2) )
        endif
      else
        Twkbdir(Zix,Nix,i,1)=0.
        Twkbtrans(Zix,Nix,i,1)=1.
        Twkbtrans(Zix,Nix,i,2)=1.
      endif
    endif
    do j = 1, nbar
      Twkb(Zix, Nix, i, j) = tff(2 * j - 1)
      if (flagfispartdamp) Twkbphase(Zix, Nix, i, j) = phase(2 * j)
      if (flagfisout) write(22, '("i=", i2, " E=", f8.3, " j=", i2, " T=", es12.3)') i, uexc, j, Twkb(Zix, Nix, i, j)
    enddo
  enddo
  if (flagfisout) close (22)
end subroutine wkb
! Copyright A.J. Koning 2021
subroutine wkbfis(iz, in, uexcit, tff, phase, tdir)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
! Definition of single and double precision variables
!   sgl           ! single precision kind
! All global variables
!   numbar        ! number of fission barriers
! Variables for fission
!   flagfisout    ! flag for output of fission information
! Constants
!   pi            ! pi
! Variables for WKB
!   iiextr        ! WKB counter
!   nextr         ! WKB counter
!   Vheight       ! height of WKB potential
!   Vwidth        ! width of WKB potential
!
! *** Declaration of local data
!
  implicit none
!
!   Roberto Capote, IAEA/NDS
!   rcapotenoy@yahoo.com
!
!   Mihaela Sin, University of Bucharest
!
!   For theory behind see
!   M. Sin, R. Capote, A. Ventura et al, Phys.Rev.C74, 014608 (2006)
!
!   WKBFIS calculates WKB momentum integrals by Gauss-Legendre method
!   for any given one-dimensional barrier and a given excitation
!   energy.
!
!   * uexcit is the excitation energy
!   * tfis is the fission transmission coefficient (see eq.(4) PRC)
!   * tff(extrema), phase(extrema) are the transmission coefficients
!   and phase integrals for all extrema
!   tff(1)= Ta, tff(3)= Tb, tff(5)= Tc (for triple humped barrier)
!
!   IMPLICIT none

  integer   :: in                         ! counter for neutrons
  integer   :: iz                         ! charge number of residual nucleus
  real(sgl) :: phase(2*numbar)            ! phase
  real(sgl) :: tdir                       ! WKB variable
  real(sgl) :: tdirv(2*numbar)            ! WKB variable
  real(sgl) :: tff(2*numbar)              ! fission transmission coefficient
  real(sgl) :: uexcit                     ! excitation energy
  real(sgl) :: dmom                       ! momentum integral
  real(sgl) :: abserr                     ! absolute error
  real(sgl) :: rmiu                       ! parameter for WKB
  real(sgl) :: epsa                       ! help variable
  real(sgl) :: epsb                       ! help variable
  real(sgl) :: dummy                      ! help variable
  real(sgl) :: eps                        ! help variable
  real(sgl) :: deps                       ! help variable
  real(sgl) :: vdef                       ! deformed potential
  real(sgl) :: phasecal                   ! help variable
  integer   :: j                          ! counter
  integer   :: ieps                       ! counter
  real(sgl) :: Fmoment                    ! function for integrand for Gauss-Legendre integration
  real(sgl) :: FmomentParab               ! function for integrand for Gauss-Legendre integration
  real(sgl) :: GaussLegendre41            ! function for Gauss Legendre integration
  real(sgl) :: FindIntersect              ! function to find intersection
  include "wkb.cmb"
  EXTERNAL FmomentParab, Fmoment
  Uexc   = uexcit
  rmiu = 0.054d0 * (iz + in) **(5.d0 / 3.d0)
!   Below is an square root of (MIU divided by 2)
  smiu = sqrt(0.5d0 * rmiu)
  do k = 1, 2 * numbar ! barriers and wells
    phase(k) = 0.d0
    tff(k)   = 0.d0
    tdirv(k) = 0.d0
  enddo
!-----------------------------------------------------------------------
!   Momentum integrals are calculated
  if (flagfisout) then
    write(22, * )
    write(22, *) ' Excitation Energy : ', uexcit
  endif
  do k = 1, nextr ! barriers and wells
      if(mod(k, 2) == 1) then
!   BARRIERS
        if(uexcit >= Vheight(k)) then
!
!   For excitation energies above barriers the transmission
!   coefficient is calculated by Hill-Wheeler formula
!
        if (Vwidth(k) > 0) then
          dmom = pi * (Vheight(k) - uexcit) / Vwidth(k)
        else
          dmom = - 50.
        endif
        phase(k)   = min(dmom, 50.)
        tff(k) = 1.d0 / (1.d0 + DEXP(2.d0 * dmom))
        if (flagfisout) write(22, '(1x, A6, I2, A10, d10.3, A3, d10.3, A15)') &
 &      ' BARR ', k, '  Mom.Int=', dmom, ' T=', tff(k), ' (Hill-Wheeler)'

        else
!
!   For excitation energies below barriers the transmission
!   coefficient is calculated by WKB method
!
! GaussLegendre41: function for Gauss Legendre integration
! Fmoment: function for integrand for Gauss-Legendre integration
! FmomentParab: function for integrand for Gauss-Legendre integration
!
        epsa = FindIntersect(uexcit, iiextr(k - 1), iiextr(k), .false.)
        epsb = FindIntersect(uexcit, iiextr(k), iiextr(k + 1), .false.)
!
!   Calculating phase integral for real shape
!
        dmom = GaussLegendre41(Fmoment, epsa, epsb, abserr)
        if(flagfisout .and. dmom > 0.d0 .and. abserr > dmom * 0.03) then
          write(22, *) ' WARNING: For extremum ', k, ' phase integral is not accurate (', &
 &        abserr/dmom*100.d0, ' %)'
        endif
        phase(k) = min(dmom, 50.)
        tff(k) = 1.d0 / (1.d0 + DEXP(2.d0 * phase(k)))
!
        if (flagfisout) then
     write(22, '(1x, A6, I2, A10, f7.4, A4, f7.4, 1x, A9, d10.3, A3, d10.3)') &
 &      ' BARR ', k, ' at beta2 ', epsa, ' to ', epsb, ' Mom.Int=', phase(k), ' T=', tff(k)
     deps = (epsb - epsa) / 50.
     eps = epsa - deps
     phasecal = 0.
     do ieps = 1, 50
       eps = eps + deps
       if (mod(ieps, 5) == 0) write(22, '(10x, "eps=", f6.3, " Vdef-E=", f8.3)') eps, &
 &       vdef(eps) - Uexc
       if (vdef(eps) >= Uexc .and. vdef(eps + deps) >= Uexc) phasecal = phasecal + 2. * smiu * ((vdef(eps) - Uexc) **0.5 + &
         (vdef(eps + deps) - Uexc) **0.5) / 2. * deps
     enddo
     if (phasecal < 30.) write(22, '(10x, " Kcal=", f10.4, " Tcal=", es12.4)') &
 &     phasecal, 1. / (1. + exp(2. * phasecal))
     endif
        endif
      else
!   WELLS
        if(uexcit.LE.Vheight(k)) then
!
!   Excitation energies below the well
!
        phase(k) = 0.d0
        else
!
!   For excitation energies above the well the transmission
!   coefficient is calculated by WKB method
!
        epsa = FindIntersect(uexcit, iiextr(k - 1), iiextr(k), .true.)
        epsb = FindIntersect(uexcit, iiextr(k), iiextr(k + 1), .true.)
!
!   Calculating phase integral for real shape
!
        dmom = GaussLegendre41(Fmoment, epsa, epsb, abserr)
     if (flagfisout .and. dmom > 0.d0 .and. abserr.gT.dmom * 0.03) write(22, *) ' WARNING: For extremum ', k, &
 &        ' phase integral is not accurate (', abserr/dmom*100.d0, ' %)'
        phase(k) = min(dmom, 50.)
        if (flagfisout) write(22, '(1x, A6, I2, A10, f7.4, A4, f7.4, 1x, A9, d10.3, A3, d10.3)') &
 &      ' WELL ', k, ' at beta2 ', epsa, ' to ', epsb, ' Mom.Int=', phase(k)
        endif
      endif
  enddo
!
!   Fission transmission for double/triple humped barrier
!
!   Iteration over barriers
!
  if (nextr > 0) tdirv(nextr) = tff(nextr)
  do k = nextr - 2, 1, - 2
     dmom = (1.d0 - tff(k)) * (1.d0 - tdirv(k + 2))
     if(k > 1) then
       tdirv(k) = tff(k) * tdirv(k + 2) / (1.d0 + dmom)
     else
       tdirv(1) = tff(1) * tdirv(3) / (1.d0 + dmom)
     endif
  enddo
    tdir = tdirv(1)
    dummy = 1.d0 / (1.d0 / tff(1) + 1.d0 / tff(3))
    if(nextr > 3) dummy = 1.d0 / (1.d0 / tff(1) + 1.d0 / tff(3) + 1.d0 / tff(5))

  if (flagfisout) write(22, '(1x, f5.2, 2x, 21(d10.3, 1x))') &
 &   uexc, tdir, dummy, (tff(j), j = 1, nextr), (phase(j), j = 1, nextr)
  RETURN
end subroutine wkbfis
! Copyright A.J. Koning 2021
function VdefParab(EPS)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
! Definition of single and double precision variables
!   sgl        ! single precision kind
! Variables for WKB
!   Vheight    ! height of WKB potential
!   Vpos       ! position of WKB potential
!   Vwidth     ! width of WKB potential
!
! *** Declaration of local data
!
  implicit none
!
!   This function calculates parabolic shape of the deformation energy
!
!   Called by gaussian integration
!
  real(sgl) :: EPS                  ! deformation
  real(sgl) :: VdefParab            ! function for parabolic shape
  include "wkb.cmb"
  VdefParab = Vheight(K) + (-1)**K*(SMIU*Vwidth(K)*(EPS-Vpos(K)))**2
end function VdefParab
! Copyright A.J. Koning 2021
function Vdef(EPS)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
! Definition of single and double precision variables
!   sgl        ! single precision kind
! Variables for WKB
!   betafis    ! fission path width
!   nbeta      ! number of beta values
!   vfis       ! adjustable factor for fission path height
!
! *** Declaration of local data
!
  implicit none
!
!   This function calculates real shape of the deformation energy
!   by linear interpolation to obtain the value of the barrier Vdef
!   at deformation EPS (needed to integrte this function)
!
!   Called by gaussian integration
!
  real(sgl) :: EPS             ! deformation
  real(sgl) :: Vdef            ! deformed potential
  integer   :: idef            ! counter for deformation
  real(sgl) :: vi              ! potential
  real(sgl) :: vip             ! potential
  real(sgl) :: ei              ! energy
  real(sgl) :: eip             !
!
!   The four lines below is a very simple (and inefficient) search
!   Should be replaced by efficient search routine to locate the
!   element idef of the array betafis()
      idef=1
      do while (EPS.GT.betafis(idef) .and. idef.LE.nbeta)
        idef = idef + 1
      enddo
      if (idef.ne.1) idef = idef - 1
      vi  = vfis(idef)
      vip = vfis(idef+1)
      ei  = betafis(idef)
      eip = betafis(idef+1)
      if(ei.eq.eip) then
!   Special case treated here to avoid division by zero
!   We assume that in this case vi = vip
        Vdef = vi
        return
      endif
      Vdef = vi + (EPS-ei)/(eip-ei)*(vip-vi)
end function Vdef
! Copyright A.J. Koning 2021
function FindIntersect(uexc, ja, jb, iswell)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
! Definition of single and double precision variables
!   sgl        ! single precision kind
! Variables for WKB
!   betafis    ! fission path width
!   vfis       ! adjustable factor for fission path height
!
! *** Declaration of local data
!
  implicit none
!
!   Calculates solutions (beta2(i)) of the equation
!   V(beta2) = Excitation Energy
!
!   If solution not found, then assign first or last point
!   of the interval depending on whether we are solving the equation
!   in the right(left) side of the barrier(well) case
!
! FindIntersect: function to find intersection
!
  integer   :: ja                       ! J value
  integer   :: jb                       ! J value
  real(sgl) :: FindIntersect            ! function to find intersection
  real(sgl) :: slope                    ! slope
  real(sgl) :: uexc                     ! excitation energy
  logical   :: iswell                   ! logical
  integer   :: j                        ! counter
  integer   :: is0                      ! index
  integer   :: is1                      ! help variable
!
! slope: slope
!
      is0 = -1
      IF(ABS(uexc-vfis(ja)).EQ.uexc-vfis(ja)) is0 = 1
      DO j=ja,jb
!   checking for a sign change in Uexc-Vdef
       is1 = -1
       IF(ABS(uexc-vfis(j)).EQ.uexc-vfis(j)) is1 = 1
       IF(is1.EQ.is0) CYCLE
!   Sign of (Uexc-vfis(j)) changed, calculating the precise value
!   of the deformation EPS at which Uexc = Vdef
       FindIntersect = betafis(j-1) + (betafis(j)-betafis(j-1))*(uexc-vfis(j-1))/(vfis(j)-vfis(j-1))
       RETURN
      ENDDO
!
!   Below is the analysis if intersection not found in [ja,jb]
!
      slope = vfis(jb) - vfis(ja)
      IF(iswell) then
!   WELLS
        IF(slope.ge.0) then
          FindIntersect = betafis(jb) ! ascending
        ELSE
          FindIntersect = betafis(ja) ! descending
        ENDIF
      ELSE
!   BARRIERS
        IF(slope.ge.0) then
          FindIntersect = betafis(ja) ! ascending
        ELSE
          FindIntersect = betafis(jb) ! descending
        ENDIF
      ENDIF
      RETURN
end function FindIntersect
! Copyright A.J. Koning 2021
function Find_Extrem (Nsmooth)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Purpose of this module
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
  use A0_talys_mod
! All global variables
!   numbar    ! number of fission barriers
! Variables for WKB
!   iiextr    ! WKB counter
!   nbeta     ! number of beta values
!   vfis      ! adjustable factor for fission path height
!
! *** Declaration of local data
!
  implicit none
!
!   Find all extremums of the smoothed deformation energy curve
!
!   Nsmooth - Number of smoothing points
!
  integer :: Find_Extrem              ! function to find extrema
  integer :: Nsmooth                  ! smoothing parameter
  logical :: logmax                   ! logical
  logical :: logmin                   ! logical
  integer :: j                        ! counter
  integer :: k                        ! designator for particle
  integer :: iext                     ! counter
  iext = 0
  Find_Extrem = 0
  iiextr(0) = 1
!C----------------------------------------------------------------------
!   determination of the minima   and maxima
  do j = nsmooth + 1 , nbeta - nsmooth

! Modification SCAMPS to prevent problem (ex Fm324):
    if (vfis(j) == vfis(j-1)) cycle
    logmax = .true.
    do k = j - nsmooth, j + nsmooth
      if (k == j) cycle
      if (vfis(k) > vfis(j)) logmax = .false.
    enddo
    if (logmax) then
      iext = iext + 1
      iiextr(iext) = j
    endif
    if (iext == 2 * numbar - 1) exit
    logmin = .true.
    do k = j - nsmooth, j + nsmooth
      if (k == j .or. k < 1) cycle
      if (vfis(k) < vfis(j)) logmin = .false.
    enddo
    if (logmin) then
      iext = iext + 1
      iiextr(iext) = j
    endif
    if (iext == 2 * numbar - 1) exit
  enddo
  Find_Extrem = iext
  iiextr(iext + 1) = nbeta
  return
end function Find_Extrem
! Copyright A.J. Koning 2021
