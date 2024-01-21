subroutine urr
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Unresolved resonance range parameters
!
! Author    : Arjan Koning and Gilles Noguere
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! Variables for OMP
!   Rprime         ! potential scattering radius
! Variables for OMP
!   RprimeU        ! potential scattering radius
! Variables for compound reactions
!   flagurrnjoy    ! normalization of URR parameters with NJOY method
!   lurr           ! maximal orbital angular momentum for URR
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Variables for binary reactions
!   xselastot      ! total elastic cross section (shape + compound)
! Variables for compound nucleus from target
!   JmaxU          ! maximal total angular momentum
!   JminU          ! minimal total angular momentum
!   lmaxU          ! maximal orbital angular momentum
!   lminU          ! minimal orbital angular momentum
!   nulj           ! (l, j) number of degrees of freedom for URR calculation
!   Turrlj         ! transmission coefficient for URR calculation
!   Turrljinc      ! incident channel (l, j) transmission coefficient for URR ca
!   xsbinarylj     ! cross section from initial compound to residual nucleus
! Variables for nuclides
!   targetspin2    ! 2 * target spin
! Constants
!   twopi          ! 2 * pi
! Variables for resonance parameters
!   Dlj            ! mean resonance spacing per J, l value
! Variables for URR
!   Rprime0        ! potential scattering radius
!   sigurrc        ! (l, j) capture cross section for URR
!   sigurrf        ! (l, j) fission cross section for URR
!   sigurrs        ! (l, j) scattering cross section for URR
!   spot           ! potential scattering contribution
!   strengthl      ! l neutron strength function
!   strengthlj     ! (l, j) neutron strength function
!   urrwidth       ! channel width in URR
!   xsurrN         ! URR cross section
!   xsurrT         ! URR cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: i     ! counter
  integer   :: istop ! integer to stop
  integer   :: J     ! spin of level
  integer   :: l     ! multipolarity
  integer   :: type  ! particle type
  real(sgl) :: err   ! error
  real(sgl) :: gJ    ! spin factor
  real(sgl) :: Rsig  ! factor for URR
  real(sgl) :: Rsig0 ! factor for URR
  real(sgl) :: Rsig1 ! factor for URR
  real(sgl) :: sum   ! help variable
  real(sgl) :: sumgJ ! sum over J values
!
! Channel widths and neutron strength function
!
! dtheory     : subroutine for theoretical calculation of average neutron spacings
! strengthfunc: function for (l,j) neutron strength function for URR
!
  call dtheory(0, 0, Einc)
  do l = lminU, min(lmaxU, lurr)
    sum = 0.
    sumgJ = 0.
    do J = JminU(l), JmaxU(l)
      if (nulj(0, l, J) > 0) Turrljinc(l, J) = Turrljinc(l, J) / real(nulj(0, l, J))
      call strengthfunc(Turrljinc(l, J), l, J)
      gJ = (2 * J + 1.) / (2. * (targetspin2 + 1.))
      sumgJ = sumgJ + gJ
      sum = sum + gJ * strengthlj(l, J)
      do type = - 1, 6
        urrwidth(type, l, J) = Dlj(l, J) * Turrlj(type, l, J) / twopi
      enddo
      urrwidth(3, l, J) = Dlj(l, J) * strengthlj(l, J) * 1.e-4
    enddo
    if (sumgJ > 0.) strengthl(l) = sum / sumgJ
  enddo
!
! Cross section after correction with NJOY formalism (lmaxU=2)
!
  if (RprimeU == 0.) RprimeU = Rprime
  if (flagurrnjoy) then
    Rprime0 = 0.1 * RprimeU
    err = 0.005
    Rsig = 0.
    Rsig0 = 1.
    Rsig1 = 1.
    istop = 0
    do while (Rsig > (1. + err) .or. Rsig < (1. - err))
      istop = istop + 1
      do i = 1, 4
        xsurrN(i) = 0.
        xsurrT(i) = 0.
      enddo
      do l = lminU, min(2, lurr)
        do J = JminU(l), JmaxU(l)
          urrwidth(1, l, J) = urrwidth(1, l, J) / Rsig0
          urrwidth( - 1, l, J) = urrwidth( - 1, l, J) * Rsig1
          call csunr2(Rprime0, l, J)
          xsurrN(3) = xsurrN(3) + sigurrf(l, J) * 1000.
          xsurrN(4) = xsurrN(4) + sigurrc(l, J) * 1000.
          xsurrT(1) = xsurrT(1) + xsbinarylj(1, l, J)
          xsurrT(3) = xsurrT(3) + xsbinarylj( - 1, l, J)
          xsurrT(4) = xsurrT(4) + xsbinarylj(0, l, J)
        enddo
      enddo
      if (xsurrN(4) > 0.) then
        Rsig0 = xsurrT(4) / xsurrN(4)
      else
        Rsig0 = 1.
      endif
      if (xsurrN(3) > 0.) then
        Rsig1 = xsurrT(3) / xsurrN(3)
      else
        Rsig1 = 1.
      endif
      if (xsurrT(4) > xsurrT(3)) then
        if (xsurrT(1) > 0.) then
          Rsig = Rsig0
          Rsig1 = 1.
        else
          Rsig = 1.
          Rsig1 = 1.
        endif
      else
        if (Rsig1 > (1. + err) .or. Rsig1 < (1. - err)) then
          Rsig = Rsig1
          Rsig0 = 1.
        else
          if (xsurrT(1) > 0.) then
            Rsig = Rsig0
            Rsig1 = 1.
          else
            Rsig = 1.
            Rsig1 = 1.
          endif
        endif
      endif
      if (istop > 100) Rsig = 1.
    enddo
  endif
!
! Final cross sections
!
  do i = 1, 4
    xsurrN(i) = 0.
    xsurrT(i) = 0.
  enddo
  do l = lminU, min(2, lurr)
    do J = JminU(l), JmaxU(l)
      if (flagurrnjoy) then
        call csunr2(Rprime0, l, J)
        xsurrN(1) = xsurrN(1) + xsbinarylj(1, l, J)
        xsurrN(2) = xsurrN(2) + sigurrs(l, J) * 1000.
        xsurrN(3) = xsurrN(3) + sigurrf(l, J) * 1000.
        xsurrN(4) = xsurrN(4) + sigurrc(l, J) * 1000.
      endif
      xsurrT(1) = xsurrT(1) + xsbinarylj(1, l, J)
      xsurrT(3) = xsurrT(3) + xsbinarylj( - 1, l, J)
      xsurrT(4) = xsurrT(4) + xsbinarylj(0, l, J)
    enddo
    if (flagurrnjoy) xsurrN(2) = xsurrN(2) + spot(l) * 1000.
  enddo
  if (flagurrnjoy) xsurrN(1) = xsurrN(1) + xsurrN(2) + xsurrN(3) + xsurrN(4)
  xsurrT(2) = xselastot
  xsurrT(1) = xsurrT(1) + xsurrT(2) + xsurrT(3) + xsurrT(4)
!
! Output of URR
!
! urrout     : subroutine for output of unresolved resonance parameters in separate files
!
  call urrout
  return
end subroutine urr
