subroutine colenhance(Zix, Nix, Eex, ald, ibar, Krot, Kvib, Kcoll)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Collective enhancement
!
! Author    : Marieke Duijvestijn, Arjan Koning, Stephane Hilaire, Stephane Goriely and Pascal Romain
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
! Variables for level density
!   cfermi          ! width of Fermi distribution for damping
!   deltaW          ! shell correction in nuclear mass
!   Exmatch         ! matching point for Ex
!   flagcol         ! flag for collective enhancement of level density
!   flagcolldamp    ! flag for damping of coll. effects in eff. level density (without explicit coll. enh.)
!   Krotconstant    ! normalization constant for rotational enhancement
!   kvibmodel       ! model for vibrational enhancement
!   ldadjust        ! logical for energy - dependent level density adjustment
!   ldmodel         ! level density model
!   pair            ! pairing energy
!   Ufermi          ! energy of Fermi distribution for damping
! Variables for fission
!   axtype          ! type of axiality of barrier
!   fismodel        ! fission model alternative fission model for default barriers
! Variables for OMP
!   beta2           ! deformation parameter
! Variables for nuclides
!   AA              ! mass number of residual nucleus
! Constants
!   onethird        ! 1 / 3
!   pi2             ! pi **2
!   twopi           ! 2 * pi
!   twothird        ! 2 / 3
! Variables for deformation parameters
!   Irigid          ! rigid body value of moment of inertia
! Variables for level density
!   delta           ! energy shift
!
! *** Declaration of local data
!
  implicit none
  character(len=132):: key            ! keyword
  integer           :: A              ! mass number of target nucleus
  integer           :: ibar           ! fission barrier
  integer           :: l              ! multipolarity
  integer           :: Nix            ! neutron number index for residual nucleus
  integer           :: Zix            ! charge number index for residual nucleus
  real(sgl)         :: ald            ! level density parameter
  real(sgl)         :: aldf           ! level density parameter
  real(sgl)         :: aldgs          ! level density parameter
  real(sgl)         :: aldint         ! level density parameter
  real(sgl)         :: avib           ! level density parameter for vibrational model
  real(sgl)         :: Cvib           ! constant for vibrational enhancement
  real(sgl)         :: damper         ! energy damping function
  real(sgl)         :: damprot        ! damping factor for rotational enhancement
  real(sgl)         :: deltaS         ! entropy change
  real(sgl)         :: deltaU         ! excitation energy change
  real(sgl)         :: Eex            ! excitation energy
  real(sgl)         :: expo           ! help variable
  real(sgl)         :: expo0          ! exponent
  real(sgl)         :: factor         ! multiplication factor
  real(sgl)         :: gammavib       ! constant for vibrational enhancement
  real(sgl)         :: ignatyuk       ! function for energy dependent level density parameter a
  real(sgl)         :: Kcoll          ! total collective enhancement
  real(sgl)         :: Kr             ! rotational enhancement factor
  real(sgl)         :: Krot           ! rotational enhancement factor
  real(sgl)         :: Krot0          ! rotational enhancement factor (undamped)
  real(sgl)         :: Kvib           ! vibrational enhancement factor
  real(sgl)         :: Kvib0          ! vibrational enhancement factor (undamped)
  real(sgl)         :: nvib           ! occupation number
  real(sgl)         :: omegavib(3)    ! energy of vibrational excitation
  real(sgl)         :: Pvib           ! total pairing correction
  real(sgl)         :: rotfactor      ! rotational enhancement factor
  real(sgl)         :: spincut        ! spin cutoff factor
  real(sgl)         :: spincutbf      ! spin-cutoff parameter squared (perpendicular projection)  for fission barrier
  real(sgl)         :: spincutgs      ! spin-cutoff parameter squared (perpendicular projection)  for ground state
  real(sgl)         :: Temp           ! nuclear temperature
  real(sgl)         :: term           ! help variable
  real(sgl)         :: Tvib           ! temperature for vibrational model
  real(sgl)         :: U              ! excitation energy minus pairing energy
!
! ************************ Collective enhancement **********************
!
! ignatyuk : function for energy dependent level density parameter a
! adjust   : subroutine for energy-dependent parameter adjustment
!
  Krot = 1.
  Krot0 = 1.
  Kvib = 1.
  Kvib0 = 1.
  Kcoll = 1.
  if ( .not. flagcol(Zix, Nix)) return
  if (ldmodel(Zix, Nix) == 1 .and. Eex < Exmatch(Zix, Nix, ibar)) return
  U = Eex - delta(Zix, Nix, ibar)
  if (U <= 0.) return
  if (ibar == 0) then
    aldgs = ald
  else
    aldgs = ignatyuk(Zix, Nix, Eex, 0)
  endif
  Temp = sqrt(U / aldgs)
  if (ldadjust(Zix, Nix)) then
    key = 'krotconstant'
    call adjust(Eex, key, Zix, Nix, 0, ibar, factor)
    Kr = factor * Krotconstant(Zix, Nix, ibar)
  else
    Kr = Krotconstant(Zix, Nix, ibar)
  endif
!
! 1. Specific collective enhancement for Bruyeres-le-Chatel (Pascal Romain) fission model
!
  if (fismodel == 1 .and. flagcolldamp) then
    if (ibar == 0) return
!
! Tri-axial barrier
!
    if (axtype(Zix, Nix, ibar) == 2) then
      aldint = ald * 8. / 13.5
      rotfactor = Kr * (U / aldint) **(0.25)
    else
!
! Axial and other barriers
!
      rotfactor = max(Kr, 1.)
    endif
!
! Damping of enhancement
!
    damprot = 1. / (1. + exp( - 0.5 * (U - 18.)))
    Krot = rotfactor * (1. - damprot) + damprot
  else
!
! 2. Default calculation
!
! Calculation of Kvib
!
    A = AA(Zix, Nix, 0)
!
! Kvibmodel 1:  Liquid drop
!
    avib = A / 13.
    Pvib = pair(Zix, Nix)
    if (Eex > Pvib) then
      Tvib = sqrt((Eex - Pvib) / avib)
    else
      Tvib = 0.
    endif
    if (Tvib > 0.) then
      if (kvibmodel == 1) then
        expo = min(0.0555 * (A **twothird) * (Tvib **(4. / 3.)), 80.)
        Kvib0 = exp(expo)
      else
!
! Kvibmodel 2: Bose gas
!
! In this case, Kvib is automatically damped at high energies
!
        deltaS = 0.
        deltaU = 0.
        Cvib = 0.0075 * (A **onethird)
        term = A **( - 5. / 6.) / (1. + 0.05 * deltaW(Zix, Nix, ibar))
        omegavib(2) = 65. * term
        omegavib(3) = 100. * term
        do l = 2, 3
          gammavib = Cvib * (omegavib(l) **2 + 4. * pi2 * Tvib **2)
          expo0 = gammavib / (2. * omegavib(l))
          expo = omegavib(l) / Tvib
          if (expo0 <= 80 .and. expo > 0..and.expo <= 80.) then
            nvib = exp( - expo0) / (exp(expo) - 1.)
            deltaS = deltaS + (2 * l + 1) * ((1. + nvib) * log(1. + nvib) - nvib * log(nvib))
            deltaU = deltaU + (2 * l + 1) * omegavib(l) * nvib
          endif
        enddo
        Kvib0 = exp(deltaS - deltaU / Tvib)
      endif
    endif
!
! Calculation of damping function and Krot
!
! damper      : energy damping function
!
! Ground state
!
    damper = 0.
    expo = (U - Ufermi(Zix, Nix, ibar)) / cfermi(Zix, Nix, ibar)
    if (expo <= 80.) damper = 1. / (1. + exp(expo))
    if (ibar == 0) then
      spincutgs = Kr * Irigid(Zix, Nix, 0) * Temp
      Krot0 = max(spincutgs, 1.)
!
! Fission barrier
!
    else
      spincutbf = Krotconstant(Zix, Nix, ibar) * Irigid(Zix, Nix, ibar) * Temp
      if (axtype(Zix, Nix, ibar) == 1) Krot0 = spincutbf
      if (axtype(Zix, Nix, ibar) == 2) Krot0 = 2. * spincutbf
      if (axtype(Zix, Nix, ibar) >= 3) then
        aldf = ignatyuk(Zix, Nix, Eex, ibar)
        term = spincutbf * sqrt(spincut(Zix, Nix, aldf, Eex, ibar, 0) * (1. - twothird * abs(beta2(Zix, Nix, ibar))))
        if (axtype(Zix, Nix, ibar) == 3) Krot0 = 0.5 * sqrt(twopi) * term
        if (axtype(Zix, Nix, ibar) == 4) Krot0 = sqrt(twopi) * term
        if (axtype(Zix, Nix, ibar) == 5) Krot0 = 2. * sqrt(twopi) * term
      endif
      Krot0 = max(Krot0, 1.)
    endif
    Krot = 1. + (Krot0 - 1.) * damper
    if (kvibmodel == 1) then
      Kvib = 1. + (Kvib0 - 1.) * damper
    else
      Kvib = Kvib0
    endif
  endif
  Kcoll = max(Krot * Kvib, 1.)
  return
end subroutine colenhance
! Copyright A.J. Koning 2021
