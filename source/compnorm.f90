subroutine compnorm
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Normalization of compound nucleus cross section
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
!   sgl             ! single precision kind
! All global variables
!   numJ            ! maximum J - value
! Variables for output
!   flagcheck       ! flag for output of numerical checks
! Variables for direct reactions
!   flagspher      ! flag to force spherical optical model
! Variables for main input
!   k0              ! index of incident particle
! Variables for astrophysics
!   flagastrogs     ! flag for calculation of astrophysics reaction rate with
! Variables for basic reaction
!   flagastro       ! flag for calculation of astrophysics reaction rate
! Variables for energies
!   Etotal          ! total energy of compound system (target + projectile)
!   flagwidth       ! flag for width fluctuation calculation
!   wavenum         ! wave number
! Variables for incident channel
!   lmaxinc         ! maximal l - value for transm. coeff. for incident channel
!   Tjlinc          ! transm. coeff. as a function of spin and l for inc. channel
!   xscoupled       ! inelastic cross section from coupled channels
!   xsdirdiscsum    ! total direct cross section
!   xsgrsum         ! sum over giant resonance cross sections
!   xspreeqsum      ! total preequilibrium cross section summed over particles
!   xsreacinc       ! reaction cross section for incident channel
! Variables to normalize compound nucleus cross section
!   CNfactor        ! factor for compound nucleus cross section: pi / [ k **2 (2s + 1)(2I + 1) ]
!   CNterm          ! compound nucleus formation cross section per spin
!   J2beg           ! begin of J summation
!   J2end           ! end of J summation
!   pardif          ! difference between target and compound nucleus parity
! Variables for nuclides
!   targetP         ! parity of target
!   targetspin      ! spin of target
!   targetspin2     ! 2 * spin of target
! Constants
!   cparity         ! parity (character)
!   parN            ! neutron number of particle
!   parspin         ! spin of particle
!   parZ            ! charge number of particle
!   pi              ! pi
!   spin2           ! 2 * spin of particle
! Variables for deformation parameters
!   colltype        ! type of collectivity (D, V or R)
! Variables for preequilibrium
!   xsflux          ! cross section flux
!
! *** Declaration of local data
!
  implicit none
  integer   :: J         ! J
  integer   :: J2        ! 2 * J
  integer   :: jj2       ! 2 * j
  integer   :: jj2beg    ! 2 * start of j summation
  integer   :: jj2end    ! 2 * end of j summation
  integer   :: l         ! multipolarity
  integer   :: lbeg      ! start of l summation
  integer   :: lend      ! end of l summation
  integer   :: Nix       ! neutron number index for residual nucleus
  integer   :: parity    ! parity
  integer   :: updown    ! spin index for transmission coefficient
  integer   :: Zix       ! charge number index for residual nucleus
  real(sgl) :: cfratio   ! compound formation ratio
  real(sgl) :: norm      ! normalization factor
  real(sgl) :: pik2      ! pi/(k**2) in mb
  real(sgl) :: xsreacsum ! reaction cross section constructed from transmission  coefficients
!
! There is a small difference between the reaction cross section as calculated by ECIS and the sum over transmission coefficients.
! We therefore normalize the level population accordingly.
!
! *** Set (kinematic) factors for compound cross section calculation ***
!
  xsreacsum = 0.
  pik2 = 10. * pi / (wavenum * wavenum)
  Zix = parZ(k0)
  Nix = parN(k0)
  CNfactor = pik2 / (2. * parspin(k0) + 1.) / (targetspin2 + 1.)
  if (k0 == 0) CNfactor = 0.5 * CNfactor
  J2beg = mod(int(2. * (targetspin + parspin(k0))), 2)
  J2end = int(2 * (lmaxinc + parspin(k0) + targetspin))
  J2end = min(J2end, 2 * numJ)
  if (flagastro .and. .not. flagwidth) J2end = 2 * numJ
  if (flagastro .and. .not. flagwidth .and. .not. flagastrogs) J2end = 2 * numJ
!
! ******************* Loop over incoming channels **********************
!
! In order to get do loops running over integer values, certain quantum numbers are multiplied by 2,
! which can be seen from a 2 present in the corresponding variable names.
! For each loop, the begin and end point is determined from the triangular rule.
!
  do parity = - 1, 1, 2
    pardif = abs(targetP - parity) / 2
!
! Loop over compound nucleus parity
!
    do J2 = J2beg, J2end, 2
!
! Loop over total angular momentum J (J2) of compound nucleus
!
      J = J2 / 2
      CNterm(parity, J) = 0.
      jj2beg = int(abs(J2 - targetspin2))
      jj2end = int(J2 + targetspin2)
!
! Loop over j (jj2) of incident channel
!

      do jj2 = jj2beg, jj2end, 2
        lbeg = int(abs(0.5 * jj2 - parspin(k0)))
        lend = int(0.5 * jj2 + parspin(k0))
!
! Loop over l of incident channel
!
        do l = lbeg, lend
!
! Check parity conservation and make index for transmission coefficient.
! Sum transmission coefficients for incident channel.
! Add partial contribution to the reaction cross section.
!
! A. Incident particles
!
          if (k0 > 0) then
            if (l > lmaxinc .or. mod(l, 2) /= pardif) cycle
            updown = int(jj2 - 2. * l) / spin2(k0)
!
! B. Incident photons
!
! Multipole radiation selection rules
!
          else
            updown = 1
            if (pardif == 0 .and. mod(l, 2) == 1) updown = 0
            if (pardif /= 0 .and. mod(l, 2) == 0) updown = 0
          endif
          CNterm(parity, J) = CNterm(parity, J) + CNfactor * (J2 + 1.) * Tjlinc(updown, l)
        enddo
      enddo
      xsreacsum = xsreacsum + CNterm(parity, J)
    enddo
  enddo
!
! Create the compound nucleus formation cross section xsflux and the associated normalization factors.
!
! For coupled-channels calculations, the transmission coefficients are already depleted by the contribution going into
! the strongly coupled inelastic channels.
!
  cfratio = 0.
  if (colltype(Zix, Nix) /= 'S' .and. .not.flagspher) then
    norm = 1.
    xsflux = xsreacsum - xsdirdiscsum + xscoupled - xspreeqsum - xsgrsum
    if (xsreacsum > 0.) cfratio = xsflux / xsreacsum
  else
    norm = xsreacsum / xsreacinc
    xsflux = xsreacinc - xsdirdiscsum - xspreeqsum - xsgrsum
    if (xsreacinc > 0.) cfratio = xsflux / xsreacinc
  endif
  if (norm > 0.) CNfactor = CNfactor * cfratio / norm
  if (flagcheck) then
    write(*, '(/" ++++++++++ Compound nucleus formation cross", " section ++++++++++"/)')
    write(*,'(" Compound nucleus excitation energy:", f15.5/)') Etotal
    write(*, '(2(" J/Pi cross section "))')
    do J2 = J2beg, J2end, 2
      J = J2/2
      write(*, '(2(f4.1, a1, es12.5, 3x))') (0.5 * J2, cparity(parity), CNterm(parity, J), parity = - 1, 1, 2)
    enddo
    write(*, '(/" ++++++++++ Normalization of reaction cross", " section ++++++++++"/)')
    write(*, '(" Reaction cross section          :", f15.5, " (A)")') xsreacinc
    write(*, '(" Sum over T(j,l)                 :", f15.5, " (B)")') xsreacsum
    write(*, '(" Compound nucleus formation c.s. :", f15.5, " (C)")') xsflux
    write(*, '(" Ratio C/B                       :", f15.5)') cfratio
  endif
  return
end subroutine compnorm
! Copyright A.J. Koning 2021
