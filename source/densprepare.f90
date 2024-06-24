subroutine densprepare(Zcomp, Ncomp, idfis)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Prepare energy grid, level density and transmission
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
!   sgl            ! single precision kind
!   dbl            ! double precision kind
! All global variables
!   numbinfis      ! maximum number of bins for fission calculation
!   numJ           ! maximum J - value
!   numl           ! number of l values
! Variables for compound reactions
!   flagfullhf     ! flag for full spin dependence of transmission coefficie
! Variables for numerics
!   transeps       ! absolute limit for transmission coefficient
! Variables for fission
!   flagfission    ! flag for fission
! Variables for input energies
!   flaginitpop   ! flag for initial population distribution
! Variables for main input
!   k0             ! index of incident particle
! Variables for gamma rays
!   gammax         ! number of l - values for gamma multipolarity
! Variables for level density
!   ldmodel        ! level density model
! Variables for energy grid
!   ebegin         ! first energy point of energy grid
!   egrid          ! outgoing energy grid
! Variables for energies
!   eend           ! last energy point of energy grid
! Variables for excitation energy grid
!   deltaEx        ! excitation energy bin for population arrays
!   Ex             ! excitation energy
!   maxex          ! maximum excitation energy bin for residual nucleus
!   maxJ           ! maximal J - value
!   nexmax         ! maximum excitation energy bin for residual nucleus
!   rhogrid        ! integrated level density
! Variables for compound nucleus from target
!   dExinc         ! excitation energy bin for mother nucleus
!   Exinc          ! excitation energy of entrance bin
!   Fnorm          ! multiplication factor
! Variables for incident channel
!   lmaxinc        ! maximal l - value for transm. coeff. for incident channel
! Variables for energy grid, level densities and transmission coefficients
!   eintfis        ! excitation energy for fission
!   lmaxhf         ! maximal l - value for transmission coefficients
!   nbintfis       ! number of bins
!   rho0           ! integrated level density
!   rhofis         ! integrated level density corresponding to tfisA
!   Tgam           ! gamma transmission coefficients
!   Tjlnex         ! transmission coefficients as a function of particle type, energy,
!   Tlnex          ! transmission coefficients as a function of particle type, energy
! Variables for inverse channel data
!   Tjl            ! transmission coefficient per particle, energy, spin and l - value
!   Tl             ! transmission coefficients per particle, energy and l - value
!  Variables for gamma-ray strength functions
!   lmax           ! maximal l - value for transmission coefficients
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   primary        ! flag to designate primary (binary) reaction
!   Zindex         ! charge number index for residual nucleus
! Constants
!   twopi          ! 2 * pi
! Variables for levels
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for level density
!   Nlast          ! last discrete level
! Variables for fission parameters
!   fecont         ! start of continuum energy
!   nfisbar        ! number of fission barrier parameters
! Variables for masses
!   S              ! separation energy
! Variables for preequilibrium
!   Ecomp          ! total energy of composite system
!
! *** Declaration of local data
!
  implicit none
  integer           :: A         ! mass number of target nucleus
  integer           :: ibar      ! fission barrier
  integer           :: ibin      ! counter
  integer           :: ibk       ! counter
  integer           :: idfis     ! fission identifier
  integer           :: Ir        ! residual spin
  integer           :: irad      ! variable to indicate M(=0) or E(=1) radiation
  integer           :: J         ! spin of level
  integer           :: l         ! multipolarity
  integer           :: ldmod     ! level density model
  integer           :: na        ! help variable
  integer           :: nb        ! help variable
  integer           :: nc        ! counter
  integer           :: Ncomp     ! neutron number index for compound nucleus
  integer           :: nen       ! energy counter
  integer           :: nexout    ! energy index for outgoing energy
  integer           :: Nix       ! neutron number index for residual nucleus
  integer           :: NL        ! last discrete level
  integer           :: NT        ! Ntop
  integer           :: odd       ! odd (1) or even (0) nucleus
  integer           :: parity    ! parity
  integer           :: Pprime    ! parity
  integer           :: type      ! particle type
  integer           :: updown    ! spin index for transmission coefficient
  integer           :: Zcomp     ! proton number index for compound nucleus
  integer           :: Zix       ! charge number index for residual nucleus
  real(sgl)         :: dEx       ! excitation energy bin for population arrays
  real(sgl)         :: dExhalf   ! half of excitation energy bin for population arrays
  real(sgl)         :: dExmin    ! energy bin
  real(sgl)         :: Ea        ! start energy of local adjustment
  real(sgl)         :: Eb        ! end energy of local adjustment
  real(sgl)         :: Ec        ! help variable for energy
  real(sgl)         :: Efs       ! fast particle energy for gamma ray strength function
  real(sgl)         :: Egamma    ! gamma energy
  real(sgl)         :: elow      ! help variable
  real(sgl)         :: elowest   ! help variable
  real(sgl)         :: emax      ! maximal emission energy within bin decay
  real(sgl)         :: emid      ! help variable
  real(sgl)         :: emin      ! minimal emission energy
  real(sgl)         :: Eout      ! outgoing energy
  real(sgl)         :: Ex0min    ! lower boundary of entrance bin
  real(sgl)         :: Ex0plus   ! upper boundary of entrance bin
  real(sgl)         :: Ex1min    ! lower boundary of residual bin
  real(sgl)         :: Ex1plus   ! upper boundary of residual bin
  real(sgl)         :: exfis     ! help variable
  real(sgl)         :: Exm       ! maximal attainable energy
  real(sgl)         :: Exout     ! excitation energy
  real(sgl)         :: fstrength ! gamma ray strength function
  real(sgl)         :: Rboundary ! factor taking into count first accessible mother bin for  discrete state
  real(sgl)         :: Rodd      ! term to determine integer or half-integer spins
  real(sgl)         :: Rspin     ! residual spin
  real(sgl)         :: SS        ! separation energy
  real(sgl)         :: ta        ! transmission coefficient
  real(sgl)         :: tb        ! transmission coefficient
  real(sgl)         :: tc        ! transmission coefficient
  real(sgl)         :: tint      ! help variable
  real(sgl)         :: weight    ! weight of discrete level
  real(dbl)         :: density   ! level density
!
! ************ Determine energetically allowed transitions *************
!
! Efs       : fast particle energy for gamma ray strength function
!
! The mother excitation energy bin is characterized by the middle Exinc, the top Ex0plus and the bottom Ex0min.
! Discrete levels above Ntop get a weight according to the level density, accounting for missing levels.
!
  Ex0plus = Exinc + 0.5 * dExinc
  Ex0min = Exinc - 0.5 * dExinc
  Efs = Exinc - S(Zcomp, Ncomp, 1)
  do type = 0, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    A = AA(Zcomp, Ncomp, type)
    NL = Nlast(Zix, Nix, 0)
    NT = Ntop(Zix, Nix, 0)
    if (NL > NT) then
      discfactor(Zix, Nix) = (Ncum(Zix, Nix, NL) - NT) / (NL - NT)
      discfactor(Zix, Nix) = min(discfactor(Zix, Nix), 2.)
      discfactor(Zix, Nix) = max(discfactor(Zix, Nix), 0.5)
    else
      discfactor(Zix, Nix) = 1.
    endif
    SS = S(Zcomp, Ncomp, type)
    odd = mod(A, 2)
    Rodd = 0.5 * odd
!
! There are 4 types of decay:
! 1. From discrete state to discrete state: This happens for the primary compound nucleus, which is formed at an energy
!    Etotal and can decay to a discrete state of a residual nucleus (e.g. compound elastic scattering).
! 2. From discrete state to continuum: This happens for the primary compound nucleus, which is formed at an energy
!    Etotal and can decay to the continuum of a residual nucleus. (e.g. continuum inelastic scattering).
! 3. From continuum to discrete state: This happens for multiple emission, where a residual nucleus can be populated in a continuum
!    excitation energy bin which can decay to a discrete state of another residual nucleus.
! 4. From continuum to continuum: This happens for multiple emission, where a residual nucleus can be populated in a continuum
!    excitation energy bin which can decay to a continuum bin of another residual nucleus.
!
! Types 3 and 4 are subject to boundary effects, i.e. they can represent cases where not the entire mother bin can have
! decayed to the residual bin or level, because of the particle separation energy.
! The end points need to be taken care of by a proper normalization.
!
    do nexout = 0, nexmax(type)
      dEx = deltaEx(Zix, Nix, nexout)
      dExhalf = 0.5 * dEx
      Exout = Ex(Zix, Nix, nexout)
      Rboundary = 1.
!
! Types 1 and 2. Decay from the primary compound nucleus.
! No special care needs to be taken.
! The residual bin/level is characterized by an excitation energy Exout.
! For transitions to the continuum, the bins is further characterized by a top Ex1plus and a bottom Ex1min.
!
      if (primary) then
        if (nexout > NL) then
          Ex1min = Exout - dExhalf
          Ex1plus = Exout + dExhalf
        endif
        Eout = Exinc - Exout - SS
      else
!
! Type 3. Decay from continuum to continuum.
! For most residual continuum bins, no special care needs to be taken and the emission energy Eout that characterizes
! the transition is simply the average between the highest energetic transition that is possible (emax, from the top of
! the mother bin to the bottom of the residual bin) and the lowest (emin).
! However, the highest residual bin (nexout=nexmax) is characterized by different energies (Ex1plus is the maximal residual
! excitation energy and Exout is shifted from its original position).
! If some transitions from mother to residual bin are forbidden, the factor Rboundary takes care of the correction.
!
        if (nexout > NL) then
          Ex1min = Exout - dExhalf
          if (nexout == nexmax(type) .and. type >= 1) then
            Ex1plus = Ex0plus - SS
            Exout = 0.5 * (Ex1plus + Ex1min)
          else
            Ex1plus = Exout + dExhalf
          endif
          emax = Ex0plus - SS - Ex1min
          emin = Ex0min - SS - Ex1plus
          Eout = 0.5 * (emin + emax)
          if (emin < 0.) then
            if (Eout > 0.) then
              Rboundary = 1. - 0.5 * (emin / (0.5 * (emax - emin))) **2
            else
              Rboundary = 0.5 * (emax / (0.5 * (emax - emin))) **2
            endif
            emin = 0.
          endif
          Eout = 0.5 * (emin + emax)
        else
!
! Type 4. Decay from continuum to discrete.
! The lowest possible mother excitation bin can not entirely decay to the discrete state.
! For the residual discrete state, it is checked whether the mother excitation bin is such a boundary case.
! This is done by adding the particle separation energy to the excitation energy of the residual discrete state.
! The correction is put in Rboundary.
!
          Exm = Exout + SS
          if (Exm <= Ex0plus .and. Exm > Ex0min) then
            Rboundary = (Ex0plus - Exm) / dExinc
            Eout = 0.5 * (Ex0plus + Exm) - SS - Exout
          else
            Eout = Exinc - SS - Exout
          endif
        endif
      endif
!
! ********************** Determine level densities *********************
!
! The calculation of level densities can be done outside many loops of various quantum numbers performed in other subroutines.
! Therefore, we store the level density as function of residual nucleus (type), excitation energy (nexout),
! spin (Ir) and parity (Pprime) in the array rho0.
!
! For discrete states, the level density is set to Rboundary.
!
      if (nexout <= NL) then
        Pprime = parlev(Zix, Nix, nexout)
        Ir = int(jdis(Zix, Nix, nexout))
        if (nexout > NT) then
          weight = Rboundary * discfactor(Zix, Nix)
        else
          weight = Rboundary
        endif
        rho0(type, nexout, Ir, Pprime) = weight
      else
!
! For decay to the continuum we use a spin and parity dependent level density.
!
        do Pprime = - 1, 1, 2
          do Ir = 0, maxJ(Zix, Nix, nexout)
            rho0(type, nexout, Ir, Pprime) = Rboundary * rhogrid(Zix, Nix, nexout, Ir, Pprime)
          enddo
        enddo
      endif
!
! ************* Interpolation of transmission coefficients *************
!
! 1. Gamma transmission coefficients
!
! adjust   : subroutine for energy-dependent parameter adjustment
! fstrength: gamma ray strength function
!
      if (type == 0) then
        lmaxhf(0, nexout) = gammax
        Egamma = Exinc - Exout
        do l = 0, gammax
          do irad = 0, 1
            Tgam(nexout, l, irad) = 0.
          enddo
        enddo
        if (Egamma <= 0) cycle
        do l = 1, gammax
          do irad = 0, 1
            Tgam(nexout, l, irad) = twopi * (Egamma **(2 * l + 1)) * fstrength(Zcomp, Ncomp, Efs, Egamma, irad, l) * Fnorm(0)
          enddo
        enddo
      else
!
! 2. Particle transmission coefficients
!
! locate    : subroutine to find value in ordered table
! pol2      : subroutine for interpolation of second order
!
        do updown = - 1, 1
          do l = 0, numl
            Tjlnex(type, nexout, updown, l) = 0.
          enddo
        enddo
        do l = 0, numl
          Tlnex(type, nexout, l) = 0.
        enddo
        lmaxhf(type, nexout) = 0
        if (ebegin(type) >= eend(type)) cycle
!
! To get the transmission coefficients on the excitation energy grid, Tjlnex, from those on the emission energy grid, Tjl,
! we use interpolation of the second order.
!
        if (Eout < egrid(ebegin(type))) then
          nen = 0
        else
          call locate(egrid, ebegin(type), eend(type), Eout, nen)
        endif
        if (nen > ebegin(type) + 1 .or. nen >= maxen - 1) then
          na = nen - 1
          nb = nen
          nc = nen + 1
        else
          na = nen
          nb = nen + 1
          nc = nen + 2
        endif
        Ea = egrid(na)
        Eb = egrid(nb)
        Ec = egrid(nc)
        do updown = - 1, 1
          do l = 0, lmax(type, nen)
            ta = Tjl(type, na, updown, l)
            tb = Tjl(type, nb, updown, l)
            tc = Tjl(type, nc, updown, l)
            call pol2(Ea, Eb, Ec, ta, tb, tc, Eout, tint)
            if (tint < transeps) tint = 0.
            Tjlnex(type, nexout, updown, l) = Fnorm(type) * tint
          enddo
        enddo
        if ( .not. flagfullhf) then
          do l = 0, lmax(type, nen)
            ta = Tl(type, na, l)
            tb = Tl(type, nb, l)
            tc = Tl(type, nc, l)
            call pol2(Ea, Eb, Ec, ta, tb, tc, Eout, tint)
            if (tint < transeps) tint = 0.
            Tlnex(type, nexout, l) = Fnorm(type) * tint
          enddo
        endif
!
! The maximal l-values needed in the compound nucleus calculations are determined.
!
        lmaxhf(type, nexout) = lmax(type, nen)
      endif
    enddo
    if (nexmax(type) > 0) lmaxhf(type, nexmax(type)) = lmaxhf(type, nexmax(type) - 1)
  enddo
  if (flaginitpop) then
    lmaxhf(k0, 0) = gammax
  else
    lmaxhf(k0, 0) = lmaxinc
  endif
!
! **** Calculate fission level densities and Hill-Wheeler terms ********
!
  if ((flagfission) .and. (idfis == 1)) then
    A = AA(Zcomp, Ncomp, 0)
    odd = mod(A, 2)
    Rodd = 0.5 * odd
    if (nfisbar(Zcomp, Ncomp) /= 0) then
      ldmod = ldmodel(Zcomp, Ncomp)
      do ibar = 1, nfisbar(Zcomp, Ncomp)
        if (primary) then
          exfis = Exinc - fecont(Zcomp, Ncomp, ibar)
        else
          exfis = Ex(Zcomp, Ncomp, maxex(Zcomp, Ncomp)) - fecont(Zcomp, Ncomp, ibar)
        endif
        nbintfis(ibar) = numbinfis / 2
        if (exfis <= 0.) cycle
        dExmin = 0.01
        dEx = exfis / nbintfis(ibar)
        if (dEx < dExmin) then
          nbintfis(ibar) = max(int(exfis / dExmin), 1)
          dEx = exfis / nbintfis(ibar)
        endif
        ibk = 1
        dExhalf = 0.5 * dEx
        elowest = fecont(Zcomp, Ncomp, ibar)
        do ibin = 0, nbintfis(ibar) - 1
          elow = elowest + ibin * dEx
          emid = elow + dExhalf
          eintfis(ibk, ibar) = elow
          eintfis(ibk + 1, ibar) = emid
          do parity = - 1, 1, 2
            do J = 0, numJ
              Rspin = J + Rodd
              rhofis(ibk, J, parity, ibar) = density(Zcomp, Ncomp, elow, Rspin, parity, ibar, ldmod)
              rhofis(ibk + 1, J, parity, ibar) = density(Zcomp, Ncomp, emid, Rspin, parity, ibar, ldmod)
            enddo
          enddo
          ibk = ibk + 2
        enddo
        nbintfis(ibar) = ibk - 1
        eintfis(nbintfis(ibar), ibar) = exfis + elowest
      enddo
    endif
  endif
  return
end subroutine densprepare
! Copyright A.J. Koning 2021
