subroutine compprepare(Zcomp, Ncomp, J2, parity)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Prepare information for initial compound nucleus
!
! Author    : Arjan Koning and Stephane Hilaire
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
!   dbl            ! double precision kind
! All global variables
!   numhill        ! maximum number of Hill - Wheeler points
!   numJ           ! maximum J - value
!   numtrans       ! number of transmission coefficients
! Variables for fission
!   flagfission    ! flag for fission
! Variables for gamma rays
!   gammax         ! number of l - values for gamma multipolarity
! Variables for main input
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for basic reaction
!   flagastro      ! flag for calculation of astrophysics reaction rate
! Variables for astrophysics
!   flagastrogs    ! flag for calculation of astrophysics reaction rate with
! Variables for energies
!   flagcompang    ! flag for compound angular distribution calculation
!   flagwidth      ! flag for width fluctuation calculation
! Variables for excitation energy grid
!   maxex          ! maximum excitation energy bin for residual nucleus
!   maxJ           ! maximal J - value
! Variables for incident channel
!   lmaxinc        ! maximal l - value for transm. coeff. for incident channel
!   Tjlinc         ! transm. coeff. as a function of spin and l for inc. channel
! Variables to normalize compound nucleus cross section
!   pardif         ! difference between target and compound nucleus parity
! Variables to prepare information for initial compound nucleus
!   denomhf        ! denominator for compound nucleus formula
!   enumhf         ! enumerator for compound nucleus formula
!   feed           ! feeding term for compound nucleus
!   fiswidth       ! fission width
!   tnum           ! counter for width fluctuation calculation
!   tNinc          ! counter for width fluctuation calculation
!   transjl        ! array for width fluctuation calculation
! Variables for energy grid, level densities and transmission coefficients
!   lmaxhf         ! maximal l - value for transmission coefficients
!   rho0           ! integrated level density
!   Tgam           ! gamma transmission coefficients
!   Tjlnex         ! transmission coefficients as a function of particle type, energy,
! Variables for fission transmission coefficients
!   rhofisA        ! integrated level density corresponding to tfisA
!   tfis           ! fission transmission coefficient for Hill - Wheeler magnitude
!   tfisA          ! transmission coefficient for Hill - Wheeler magnitude
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   targetP        ! parity of target
!   targetspin     ! spin of target
!   targetspin2    ! 2 * spin of target
!   Zindex         ! charge number index for residual nucleus
! Constants
!   parspin        ! spin of particle
!   spin2          ! 2 * spin of particle
! Variables for levels
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for fission parameters
!   nfisbar        ! number of fission barrier parameters
! Variables for level density
!   Nlast          ! last discrete level
! Variables for astro
!   rhoastrotot    ! total level density for astrophysical case
!   Tastrotot      ! total transmission coefficient for astrophysical case
! Variables for initial compound nucleus
!   wpower         ! power used for rho * (t **wpower)
!
! *** Declaration of local data
!
  implicit none
  integer   :: i           ! counter
  integer   :: ihill       ! counter for Hill-Wheeler magnitude
  integer   :: Ir          ! residual spin
  integer   :: irad        ! variable to indicate M(=0) or E(=1) radiation
  integer   :: Irspin2     ! 2 * residual spin
  integer   :: Irspin2beg  ! 2 * start of residual spin summation
  integer   :: Irspin2end  ! 2 * end of residual spin summation
  integer   :: J           ! spin of level
  integer   :: J2          ! 2 * J
  integer   :: J2maxastro  ! maximum J value for astro
  integer   :: J2res       ! help variable
  integer   :: jj2         ! 2 * j
  integer   :: jj2beg      ! 2 * start of j summation
  integer   :: jj2end      ! 2 * end of j summation
  integer   :: jj2prime    ! 2 * j'
  integer   :: jj2primebeg ! 2 * start of j' summation
  integer   :: jj2primeend ! 2 * end of j' summation
  integer   :: l           ! multipolarity
  integer   :: l2          ! 2 * l
  integer   :: l2beg       ! 2 * start of l summation
  integer   :: l2end       ! 2 * end of l summation
  integer   :: l2maxhf     ! 2 * lmaxhf
  integer   :: l2prime     ! 2 * l'
  integer   :: l2primebeg  ! 2 * start of l summation
  integer   :: l2primeend  ! 2 * end of l summation
  integer   :: lprime      ! 2 * l
  integer   :: modl        ! help variable
  integer   :: Ncomp       ! neutron number index for compound nucleus
  integer   :: nexout      ! energy index for outgoing energy
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: pardif2     ! difference between residual and compound nucleus parity
  integer   :: parity      ! parity
  integer   :: parspin2i   ! 2 * particle spin for incident channel
  integer   :: parspin2o   ! 2 * particle spin for outgoing channel
  integer   :: Pprime      ! parity
  integer   :: Pprimebeg   ! start of residual parity summation
  integer   :: Pprimeend   ! end of residual parity summation
  integer   :: pspin2i     ! 2 * spin of particle (usually) for incident channel
  integer   :: pspin2o     ! 2 * spin of particle (usually) for incident channel
  integer   :: type        ! particle type
  integer   :: updown      ! spin index for transmission coefficient
  integer   :: updown2     ! spin index for transmission coefficient
  integer   :: Zcomp       ! proton number index for compound nucleus
  integer   :: Zix         ! charge number index for residual nucleus
  real(sgl) :: Tinc        ! transmission coefficients as a function of j and l  for the incident channel
  real(sgl) :: Tout        ! transmission coefficients
  real(dbl) :: factor      ! multiplication factor
  real(dbl) :: gamwidth    ! sum over all gamma transmission coefficients
  real(dbl) :: ratio       ! if 0<ratio<1 then x is between xl1 and xl2
  real(dbl) :: rho         ! integrated level density
  real(dbl) :: tfishill    ! help variable
!
! The transmission coefficients are put in arrays for possible width fluctuation calculations.
! Also the total width denomhf appearing in the denominator of the compound nucleus formula is created.
! Note that the complete subroutine is performed inside the loop over J and P in subroutine comptarget.
!
! *************************** Initialization ***************************
!
  J = J2 / 2
  denomhf = 0.
  tnum = 0
  feed = 0.
  parspin2i = int(2. * parspin(k0))
  pspin2i = spin2(k0)
!
! ********************* Loop over incident channels ********************
!
! In order to get do-loops running over integer values, certain quantum numbers are multiplied by 2,
! which can be seen from a 2 present in the corresponding variable names.
! For each loop, the begin and end point is determined from the triangular rule.
!
! We have to distinguish between incident particles and incident photons.
!
! A. Incident particles
!
  if (k0 > 0) then
    jj2beg = abs(J2 - targetspin2)
    jj2end = J2 + targetspin2
!
! Loop over j (jj2) of incident channel
!
    do jj2 = jj2beg, jj2end, 2
      l2beg = abs(jj2 - parspin2i)
      l2end = jj2 + parspin2i
      l2end = min(l2end, 2 * lmaxinc)
!
! Loop over l (l2) of incident channel
!
      do l2 = l2beg, l2end, 2
        l = l2 / 2
!
! Check parity conservation and make index for transmission coefficient.
!
! If the parity of the target nucleus is equal (unequal) to the parity of compound nucleus, i.e. pardif=0(1),
! the l-value must be even (odd).
!
        if (mod(l, 2) /= pardif) cycle
        updown = (jj2 - l2) / pspin2i
        Tinc = Tjlinc(updown, l)
!
! Information needed for width fluctuation calculation.
!
! For the width fluctuation calculation, all transmission coefficients need to be placed in one sequential array.
! Therefore, a counter tnum needs to be followed to keep track of the proper index for the transmission coefficients.
! The order inside the transjl array is:
!
! 1. Incident channel
! 2. Outgoing particle channels
! 3. Gamma channel
! 4. Fission channel (if present)
!
        if (flagwidth) then
          tnum = tnum + 1
          transjl(0, tnum) = 1.
          do i = 1, wpower
            transjl(i, tnum) = 0.
            if (Tinc > 1.e-30 **(1. / i)) transjl(i, tnum) = Tinc **i
          enddo
        else
          if ( .not. flagcompang) feed = feed + Tinc
        endif
      enddo
    enddo
  else
!
! B. Incident photons
!
! Note that for photons, the first index of Tjlinc represents radiation type (M or E)
!
    Tinc = 0.
    if (J2 /= targetspin2 .or. J2 /= 0) then
      do l = 1, gammax
        if (0.5 * J2 < targetspin - l .or. 0.5 * J2 > targetspin + l) cycle
        irad = 1
        if (parity == targetP .and. mod(l, 2) == 1) irad = 0
        if (parity /= targetP .and. mod(l, 2) == 0) irad = 0
        Tinc = Tinc + Tjlinc(irad, l)
      enddo
      feed = feed + Tinc
    endif
  endif
!
! There are two possible types of calculation for the initial compound nucleus.
! If either width fluctuation corrections or compound nucleus angular distributions are wanted, we need to sum explicitly over all
! possible quantum numbers before we calculate the width fluctuation or angular factor.
! If not, the sum over j,l of the transmission coefficients can be lumped into one factor, which decreases the calculation time.
! In the latter case, the partial decay widths are stored in enumhf.
!
  enumhf(k0, Ltarget, int(targetspin), targetP) = feed
  tNinc = tnum
!
! ********************* Loop over outgoing channels ********************
!
! 1. Fission
!
! The fission contribution is calculated and added to the total decay width.
!
  if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) then
    fiswidth = tfis(J, parity)
    denomhf = denomhf + fiswidth
  endif
!
! 2. Gamma and particle channels
!
  gamwidth = 0.
  do type = 0, 6
    if (parskip(type)) cycle
    parspin2o = int(2. * parspin(type))
    pspin2o = spin2(type)
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
!
! This loop is over all discrete levels and continuum bins of the final nucleus.
!
    do nexout = 0, maxex(Zix, Nix)
      l2maxhf = 2 * lmaxhf(type, nexout)
!
! Initialization of summations
!
! For discrete states, the begin and end points of the residual spin/parity summation are both set equal to the residual discrete
! level spin/parity.
!
      if (nexout <= Nlast(Zix, Nix, 0)) then
        Pprimebeg = parlev(Zix, Nix, nexout)
        Pprimeend = Pprimebeg
        Irspin2beg = int(2. * jdis(Zix, Nix, nexout))
        Irspin2end = Irspin2beg
      else
!
! For the continuum, the begin and end points of the residual spin/parity summation are set to the maximally accessible values.
!
        Pprimebeg = - 1
        Pprimeend = 1
        J2res = J2 + parspin2o
        Irspin2beg = mod(J2res, 2)
        Irspin2end = J2res + l2maxhf
        Irspin2end = min(Irspin2end, 2 * maxJ(Zix, Nix, nexout))
      endif
!
! The variable pardif2 is used as an indicator of parity conservation for the outgoing channel.
!
! Loop over residual parity
!
      do Pprime = Pprimebeg, Pprimeend, 2
        pardif2 = abs(parity - Pprime) / 2
!
! Loop over residual spin
!
        do Irspin2 = Irspin2beg, Irspin2end, 2
          Ir = Irspin2 / 2
          enumhf(type, nexout, Ir, Pprime) = 0.
          rho = rho0(type, nexout, Ir, Pprime)
          if (rho < 1.e-20) cycle
          if (flagastro .and. .not. flagastrogs) then
            J2maxastro = int(2 * (lmaxhf(k0, nexout) + parspin(k0) + Ir))
            J2maxastro = min(J2maxastro, numJ)
          endif
          jj2primebeg = abs(J2 - Irspin2)
          jj2primeend = J2 + Irspin2
!
! Loop over j (jj2) of outgoing channel
!
          do jj2prime = jj2primebeg, jj2primeend, 2
            l2primebeg = abs(jj2prime - parspin2o)
            if (type == 0) l2primebeg = max(l2primebeg, 2)
            l2primeend = jj2prime + parspin2o
            l2primeend = min(l2primeend, l2maxhf)
!
! Loop over l (l2) of outgoing channel
!
            do l2prime = l2primebeg, l2primeend, 2
              lprime = l2prime / 2
              modl = mod(lprime, 2)
!
! 1. Photons
!
! We include photons as a special case, with the multipole radiation selection rules (irad=0: M-transition, irad=1: E-transition)
!
              if (type == 0) then
                if (pardif2 == modl) then
                  irad = 1
                else
                  irad = 0
                endif
                Tout = Tgam(nexout, lprime, irad)
              else
!
! 2. Particles
!
! If the parity of the residual nucleus is equal (unequal) to the parity of compound nucleus,i.e. pardif2=0(1),
! the l-value must be even (odd).
!
                if (modl /= pardif2) cycle
                updown2 = (jj2prime - l2prime) / pspin2o
                Tout = Tjlnex(type, nexout, updown2, lprime)
              endif
!
! The contribution is added to the total width.
!
              factor = rho * Tout
              denomhf = denomhf + factor
!
! Information needed for width fluctuation calculation.
! Values for rho*T**i where (i=0,5) are stored.
! The photon contribution is stored in a single gamma width.
!
              if (flagwidth) then
                if (type == 0) then
                  gamwidth = gamwidth + factor
                else
                  if (tnum < numtrans - numhill - 2) then
                    tnum = tnum + 1
                    transjl(0, tnum) = rho
                    do i = 1, wpower
                      transjl(i, tnum) = 0.
                      if (Tout > 1.e-30 **(1. / i)) transjl(i, tnum) = rho * (Tout **i)
                    enddo
                  endif
                endif
              else
!
! If NO width fluctuation corrections or angular distributions are required, the partial decay width is used in
! subroutine comptarget.
!
                if ( .not. flagcompang) enumhf(type, nexout, Ir, Pprime) = &
                  enumhf(type, nexout, Ir, Pprime) + factor
              endif
!
! Astrophysical case
!
              if (flagastro .and. .not. flagastrogs .and. type == k0) then
                if (Tout > 0.) then
                  rhoastrotot = rhoastrotot + rho
                  if (nexout /= Ltarget) Tastrotot = Tastrotot + factor
                endif
                if (type == k0 .and. nexout /= Ltarget) then
                  if (J2 <= J2maxastro) feed = feed + factor
                endif
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
!
! ******** Add fission and gamma transmission coefficients to transmission coefficient array for width fluctuations *******
!
! 1. Fission.
!
  if (flagwidth) then
    if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) then
      do ihill = 1, numhill
        ratio = 0.
        if (tfisA(J, parity, 0) > 0.) ratio = tfisA(J, parity, ihill) / tfisA(J, parity, 0)
        tfishill = ratio * fiswidth
        tnum = tnum + 1
        transjl(0, tnum) = max(rhofisA(J, parity, ihill), 1.d0)
        do i = 1, wpower
          transjl(i, tnum) = 0.
          if (tfishill > 1.e-30 **(1. / i)) transjl(i, tnum) = transjl(0, tnum) * (tfishill **i) / (transjl(0, tnum) **i)
        enddo
      enddo
    endif
!
! 2. Photons.
!
    transjl(0, tnum + 1) = 1.
    do i = 1, wpower
      transjl(i, tnum + 1) = 0.
      if (gamwidth > 1.e-30 **(1. / i) .and. gamwidth < 1.e30 **(1. / i)) transjl(i, tnum + 1) = gamwidth **i
    enddo
  endif
  return
end subroutine compprepare
! Copyright A.J. Koning 2021
