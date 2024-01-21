subroutine multipreeq(Zcomp, Ncomp, nex)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Multiple preequilibrium model
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
!   sgl           ! single precision kind
!   dbl           ! double precision kind
! All global variables
!   numex         ! maximum number of excitation energies
!   numNph        ! maximum number of neutrons away from the initial compound nucleus
!   numZph        ! maximum number of protons away from the initial compound nucleus
! Variables for preequilibrium
!   flaggshell    ! flag for energy dependence of single particle level density parameter g
!   g             ! single - particle level density parameter
!   mpreeqmode    ! designator for multiple pre - equilibrium model
! Variables for level density
!   alev          ! level density parameter
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
!   egrid         ! outgoing energy grid
! Variables for energies
!   eend          ! last energy point of energy grid
!   mulpreZN      ! logical for multiple pre - equilibrium per nucleus
! Variables for excitation energy grid
!   deltaEx       ! excitation energy bin for population arrays
!   Ex            ! excitation energy
!   maxJ          ! maximal J - value
!   nexmax        ! maximum excitation energy bin for residual nucleus
! Variables for multiple emission
!   Dmulti        ! depletion factor for multiple preequilibrium
!   mcontrib      ! contribution to emission spectrum
!   mpecontrib    ! contribution to multiple pre - equilibrium emission spectr
!   xsfeed        ! cross section from compound to residual nucleus
!   xsmpe         ! multiple - preequilibrium cross section per energy bin
!   xsmpetot      ! total multiple - preequilibrium cross section
!   xspartial     ! emitted cross section flux per energy bin
!   xspoppreeq    ! preequilibrium population cross section per nucleus
! Variables for compound nucleus from target
!   Exinc         ! excitation energy of entrance bin
!   dExinc        ! excitation energy bin for mother nucleus
! Variables for incident channel
!   popdecay      ! decay from population
!   preeqpop      ! pre - equilibrium population
!   preeqpopex    ! pre - equilibrium population
!   xspop         ! population cross section
!   xspopex       ! population cross section summed over spin and parity
!   xspopnuc      ! population cross section per nucleus
! Variables for inverse channel data
!   Tjl           ! transmission coefficient per particle, energy, spin and l - value
! Variables for nuclides
!   Nindex        ! neutron number index for residual nucleus
!   parskip       ! logical to skip outgoing particle
!   Zindex        ! charge number index for residual nucleus
! Variables for level density
!   Nlast         ! last discrete level
! Variables for masses
!   S             ! separation energy
! Variables for exciton model
!   tauexc        ! lifetime of exciton state
! Variables for preequilibrium initialization
!   Efermi        ! depth of Fermi well
!   maxpar        ! maximal particle number
!   numparx       ! maximum number of particles
!   Rblann        ! Blann's factor
!   RnJ           ! spin distribution for particle - hole states
!   RnJsum        ! (2J + 1) * sum over spin distributions
! Variables for preequilibrium
!   Ecomp         ! total energy of composite system
!   h0            ! initial hole number
!   p0            ! initial particle number
!   wemission     ! emission rate per particle, exciton number and energy
!   xspopph       ! population cross section per particle-hole configuration
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell                      ! flag for surface effects in finite well
  integer   :: h                             ! help variable
  integer   :: ih                            ! hole number
  integer   :: ip                            ! particle number
  integer   :: itype                         ! help variable
  integer   :: J                             ! spin of level
  integer   :: Ncomp                         ! neutron number index for compound nucleus
  integer   :: nen                           ! energy counter
  integer   :: nex                           ! excitation energy bin of compound nucleus
  integer   :: nexout                        ! energy index for outgoing energy
  integer   :: Nix                           ! neutron number index for residual nucleus
  integer   :: p                             ! particle number
  integer   :: parity                        ! parity
  integer   :: type                          ! particle type
  integer   :: Zcomp                         ! proton number index for compound nucleus
  integer   :: Zix                           ! charge number index for residual nucleus
  integer   :: ZNcomp                        ! help variable
  real(sgl) :: dEx                           ! excitation energy bin for population arrays
  real(sgl) :: Eex                           ! excitation energy
  real(sgl) :: Eo                            ! outgoing energy grid based on excitation energy
  real(sgl) :: EoplusS                       ! outgoing energy+S
  real(sgl) :: Exm                           ! maximal attainable energy
  real(sgl) :: Exmin                         ! help variable
  real(sgl) :: factor(2, 0:numex, 0:numparx) ! multiplication factor
  real(sgl) :: gs                            ! single-particle level density parameter
  real(sgl) :: ignatyuk                      ! function for energy dependent level density
  real(sgl) :: Jterm                         ! term dependent on J
  real(sgl) :: omega1p                       ! particle-hole state density for continuum particle
  real(sgl) :: omegap1h                      ! particle-hole state density for residual system
  real(sgl) :: omegaph                       ! particle-hole state density for compound system
  real(sgl) :: Pescape                       ! escape probability
  real(sgl) :: phdens                        ! function for particle-hole state density
  real(sgl) :: proba                         ! probability to find continuum particle
  real(sgl) :: Rfactor                       ! Blann's factor
  real(sgl) :: sumfeed                       ! help variable
  real(sgl) :: term(2, 0:numex)              ! help variable
  real(sgl) :: Tswave                        ! transmission coefficient for s-wave
  real(dbl) :: feedph                        ! feeding term from previous particle-hole calculation
  real(dbl) :: summpe                        ! multiple preequilibrium emission for excitation energy bin
  real(dbl) :: sumph                         ! multiple preequilibrium flux for particle-hole pair
  real(dbl) :: sumterm                       ! help variable
  real(dbl) :: sumtype(2)                    ! help variable
!
! Multiple preequilibrium emission model 2 is adopted from Chadwick and Young, Phys. Rev. C50 (1994) p. 996.
!
! ************ Check presence of multiple pre-equilibrium **************
!
! There must be excited particles and holes present for multiple pre-equilibrium emission to occur.
!
  if (Zcomp == 0 .and. Ncomp == 0) return
  sumfeed = 0.
  factor = 0.
  do ip = 1, maxpar
    do ih = 1, maxpar
      sumfeed = sumfeed + xspopph(Zcomp, Ncomp, nex, ip, ih)
    enddo
  enddo
  if (sumfeed <= 1.e-10) return
!
! ***************** Multiple preequilibrium emission *******************
!
! ignatyuk      : function for energy dependent level density parameter a
!
  Ecomp = Exinc
  gs = g(Zcomp, Ncomp)
  if (flaggshell) gs = gs * ignatyuk(Zcomp, Ncomp, Ecomp, 0) / alev(Zcomp, Ncomp)
  surfwell = .false.
  Rfactor = 0.5
  ZNcomp = Zcomp + Ncomp
!
! For secondary pre-equilibrium emission, determine the type of the first emitted particle.
!
  if (ZNcomp == 1) then
    if (Ncomp == 1) then
      itype = 1
    else
      itype = 2
    endif
  endif
  summpe = 0.
!
! Loops over all possible particle-hole excitations of the mother bin.
!
! emissionrate: subroutine for emission rate
! lifetime    : subroutine for calculation of lifetime of exciton state
!
  do ip = 1, maxpar
    do ih = 1, maxpar
      feedph = xspopph(Zcomp, Ncomp, nex, ip, ih)
      if (feedph <= 1.e-10) cycle
      sumph = 0.
      if (mpreeqmode == 2) then
        omegaph = phdens(Zcomp, Ncomp, ip, ih, gs, Exinc, Efermi, surfwell)
        omegaph = max(omegaph, 1.)
      else
        p0 = ip
        h0 = ih
        do p = p0, maxpar
          h = h0 + p - p0
          if (h > maxpar) cycle
          call emissionrate(Zcomp, Ncomp, p, h)
          call lifetime(Zcomp, Ncomp, p, h)
        enddo
      endif
      sumterm = 0.
      do type = 1, 2
        sumtype(type) = 0.
        do nexout = 0, numex
          term(type, nexout) = 0.
        enddo
        if (parskip(type)) cycle
        Zix = Zindex(Zcomp, Ncomp, type)
        Nix = Nindex(Zcomp, Ncomp, type)
        if (Zix > numZph .or. Nix > numNph) cycle
        if (ZNcomp == 1) Rfactor = Rblann(itype, type, ip)
!
! locate    : subroutine to find value in ordered table
!
        do nexout = Nlast(Zix, Nix, 0) + 1, nexmax(type)
          dEx = deltaEx(Zix, Nix, nexout)
          Eex = Ex(Zix, Nix, nexout)
          Eo = Exinc - Eex - S(Zcomp, Ncomp, type)
          call locate(egrid, ebegin(type), eend(type), Eo, nen)
          if (mpreeqmode == 2) then
            gs = g(Zix, Nix)
            if (flaggshell) gs = g(Zix, Nix) * ignatyuk(Zix, Nix, Eex, 0) / alev(Zix, Nix)
            omegap1h = phdens(Zix, Nix, ip - 1, ih, gs, Eex, Efermi, surfwell)
            EoplusS = Exinc - Eex
            omega1p = phdens(Zix, Nix, 1, 0, gs, EoplusS, Efermi, surfwell)
            proba = omega1p * omegap1h / omegaph / ip * Rfactor
            Tswave = Tjl(type, nen, 1, 0)
            Pescape = proba * Tswave
            if (nexout == nexmax(type)) then
              Exm = Exinc + 0.5 * dExinc - S(Zcomp, Ncomp, type)
              Exmin = Ex(Zix, Nix, nexout) - 0.5 * dEx
              dEx = Exm - Exmin
            endif
            term(type, nexout) = feedph * Pescape * dEx
            sumterm = sumterm + term(type, nexout)
          else
            do p = p0, maxpar
              h = h0 + p - p0
              if (h > maxpar) cycle
              factor(type, nexout, p) = feedph * tauexc(p, h) * wemission(type, p, h, nen) * Rfactor * dEx
              term(type, nexout) = term(type, nexout) + factor(type, nexout, p)
            enddo
            sumterm = sumterm + term(type, nexout)
          endif
          sumtype(type) = sumtype(type) + term(type, nexout)
        enddo
      enddo
!
! Normalization
!
      if (sumterm > feedph) then
        do type = 1, 2
          do nexout = Nlast(Zix, Nix, 0) + 1, nexmax(type)
            term(type, nexout) = term(type, nexout) * feedph / sumterm
            if (mpreeqmode == 1) then
              do p = p0, maxpar
                factor(type, nexout, p) = factor(type, nexout, p) * feedph / sumterm
              enddo
            endif
          enddo
          sumtype(type) = sumtype(type) * feedph / sumterm
        enddo
      endif
!
! Feed new population bins
!
      do type = 1, 2
        if (parskip(type)) cycle
        Zix = Zindex(Zcomp, Ncomp, type)
        Nix = Nindex(Zcomp, Ncomp, type)
        if (Zix > numZph .or. Nix > numNph) cycle
        do nexout = Nlast(Zix, Nix, 0) + 1, nexmax(type)
          if (mpreeqmode == 2) then
            xspopph(Zix, Nix, nexout, ip - 1, ih) = xspopph(Zix, Nix, nexout, ip - 1, ih) + term(type, nexout)
          else
            do p = p0, maxpar
              h = h0 + p - p0
              if (h > maxpar) cycle
              xspopph(Zix, Nix, nexout, p - 1, h) = xspopph(Zix, Nix, nexout, p - 1, h) + factor(type, nexout, p)
            enddo
          endif
          mcontrib(type, nex, nexout) = mcontrib(type, nex, nexout) + term(type, nexout)
          mpecontrib(type, nex, nexout) = mpecontrib(type, nex, nexout) + term(type, nexout)
          xspopex(Zix, Nix, nexout) = xspopex(Zix, Nix, nexout) + term(type, nexout)
          preeqpopex(Zix, Nix, nexout) = preeqpopex(Zix, Nix, nexout) + term(type, nexout)
          do parity = - 1, 1, 2
            do J = 0, maxJ(Zix, Nix, nexout)
              Jterm = 0.5 * (2 * J + 1) * RnJ(2, J) / RnJsum(2) * term(type, nexout)
              xspop(Zix, Nix, nexout, J, parity) = xspop(Zix, Nix, nexout, J, parity) + Jterm
              popdecay(type, nexout, J, parity) = popdecay(type, nexout, J, parity) + Jterm
              preeqpop(Zix, Nix, nexout, J, parity) = preeqpop(Zix, Nix, nexout, J, parity) + Jterm
            enddo
          enddo
        enddo
!
! Add total multiple pre-equilibrium contributions
!
        sumph = sumph + sumtype(type)
        summpe = summpe + sumtype(type)
        xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) + sumtype(type)
        xspoppreeq(Zix, Nix) = xspoppreeq(Zix, Nix) + sumtype(type)
        xspartial(type, nex) = xspartial(type, nex) + sumtype(type)
        xsfeed(Zcomp, Ncomp, type) = xsfeed(Zcomp, Ncomp, type) + sumtype(type)
        xsmpe(type, nex) = xsmpe(type, nex) + sumtype(type)
        xsmpetot(type) = xsmpetot(type) + sumtype(type)
        if (sumtype(type) /= 0.) mulpreZN(Zix, Nix) = .true.
      enddo
!
! Flux that is not emitted during a particular stage, is always transferred to the next stage.
!
      if (mpreeqmode == 2) then
        if (ip <= maxpar - 1 .and. ih <= maxpar - 1) xspopph(Zcomp, Ncomp, nex, ip + 1, ih + 1) = &
          xspopph(Zcomp, Ncomp, nex, ip + 1, ih + 1) + feedph - sumph
      endif
    enddo
  enddo
!
! ************************ Normalization *******************************
!
  Dmulti(nex) = summpe / xspopex(Zcomp, Ncomp, nex)
  xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) - summpe
  preeqpopex(Zcomp, Ncomp, nex) = preeqpopex(Zcomp, Ncomp, nex) - summpe
  return
end subroutine multipreeq
! Copyright A.J. Koning 2021
