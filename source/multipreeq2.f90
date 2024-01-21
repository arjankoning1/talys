subroutine multipreeq2(Zcomp, Ncomp, nex)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Two-component multiple preequilibrium model
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
!   gn            ! single - particle neutron level density parameter
!   gp            ! single - particle proton level density parameter
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
!   dExinc        ! excitation energy bin for mother nucleus
!   Exinc         ! excitation energy of entrance bin
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
! Constants
!   parN          ! neutron number of particle
!   parZ          ! charge number of particle
! Variables for level density
!   Nlast         ! last discrete level
! Variables for masses
!   S             ! separation energy
! Variables for preequilibrium initialization
!   Efermi        ! depth of Fermi well
!   maxpar        ! maximal particle number
!   numparx       ! maximum number of particles
!   RnJ           ! spin distribution for particle - hole states
!   RnJsum        ! (2J + 1) * sum over spin distributions
! Variables for preequilibrium
!   Ecomp         ! total energy of composite system
!   hnu0          ! initial neutron hole number
!   hpi0          ! initial proton hole number
!   p0            ! initial particle number
!   pnu0          ! initial neutron number
!   ppi0          ! initial proton number
!   Spre          ! time - integrated strength of two - component exciton state
!   wemission2    ! two - component emission rate
!   xspopph2      ! population cross section p
!
! *** Declaration of local data
!
  implicit none
  logical   :: surfwell                                 ! flag for surface effects in finite well
  integer   :: h                                        ! help variable
  integer   :: hnu                                      ! neutron hole number
  integer   :: hpi                                      ! proton hole number
  integer   :: ih                                       ! hole number
  integer   :: ihn                                      ! neutron hole number of mother nucleus
  integer   :: ihp                                      ! proton hole number of mother nucleus
  integer   :: ip                                       ! particle number
  integer   :: ipn                                      ! neutron particle number of mother nucleus
  integer   :: ipp                                      ! proton particle number of mother nucleus
  integer   :: J                                        ! spin of level
  integer   :: Ncomp                                    ! neutron number index for compound nucleus
  integer   :: nejec                                    ! neutron number of leading particle
  integer   :: nen                                      ! energy counter
  integer   :: nex                                      ! excitation energy bin of compound nucleus
  integer   :: nexout                                   ! energy index for outgoing energy
  integer   :: Nix                                      ! neutron number index for residual nucleus
  integer   :: p                                        ! particle number
  integer   :: parity                                   ! parity
  integer   :: pnu                                      ! neutron particle number
  integer   :: ppi                                      ! proton particle number
  integer   :: type                                     ! particle type
  integer   :: Zcomp                                    ! proton number index for compound nucleus
  integer   :: zejec                                    ! charge number of leading particle
  integer   :: Zix                                      ! charge number index for residual nucleus
  real(sgl) :: damp                                     ! shell damping factor
  real(sgl) :: dEx                                      ! excitation energy bin for population arrays
  real(sgl) :: Eex                                      ! excitation energy
  real(sgl) :: Eo                                       ! outgoing energy grid based on excitation energy
  real(sgl) :: EoplusS                                  ! outgoing energy+S
  real(sgl) :: Exm                                      ! maximal attainable energy
  real(sgl) :: Exmin                                    ! help variable
  real(sgl) :: factor(2, 0:numex, 0:numparx, 0:numparx) ! multiplication factor
  real(sgl) :: gsn                                      ! single-particle neutron level density parameter
  real(sgl) :: gsp                                      ! single-particle proton level density parameter
  real(sgl) :: ignatyuk                                 ! function for energy dependent level density parameter a
  real(sgl) :: Jterm                                    ! term dependent on J
  real(sgl) :: omega1p                                  ! particle-hole state density for continuum particle
  real(sgl) :: omegap1h                                 ! particle-hole state density for residual system
  real(sgl) :: omegaph                                  ! particle-hole state density for compound system
  real(sgl) :: Pescape                                  ! escape probability
  real(sgl) :: phdens2                                  ! function for two-component particle-hole state density
  real(sgl) :: proba                                    ! probability to find continuum particle
  real(sgl) :: sumfeed                                  ! help variable
  real(sgl) :: term(2, 0:numex)                         ! help variable
  real(sgl) :: Tswave                                   ! transmission coefficient for s-wave
  real(dbl) :: feedph                                   ! feeding term from previous particle-hole calculation
  real(dbl) :: summpe                                   ! multiple preequilibrium emission for excitation energy bin
  real(dbl) :: sumph                                    ! multiple preequilibrium flux for particle-hole pair
  real(dbl) :: sumterm                                  ! help variable
  real(dbl) :: sumtype(2)                               ! help variable
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
  do ipp = 0, maxpar
    do ihp = 0, maxpar
      do ipn = 0, maxpar
        do ihn = 0, maxpar
          sumfeed = sumfeed + xspopph2(Zcomp, Ncomp, nex, ipp, ihp, ipn, ihn)
        enddo
      enddo
    enddo
  enddo
  if (sumfeed <= 1.e-10) return
!
! ***************** Multiple preequilibrium emission *******************
!
! ignatyuk  : function for energy dependent level density parameter a
!
  surfwell = .false.
  Ecomp = Exinc
  gsp = gp(Zcomp, Ncomp)
  gsn = gn(Zcomp, Ncomp)
  if (flaggshell) then
    damp = ignatyuk(Zcomp, Ncomp, Ecomp, 0) / alev(Zcomp, Ncomp)
    gsp = gsp * damp
    gsn = gsn * damp
  endif
  summpe = 0.
!
! phdens2   : function for two-component particle-hole state density
! exchange2 : subroutine for calculation of two-component exchange terms
! lifetime2 : subroutine for calculation of lifetime of two-component exciton state
!
  do ipp = 0, maxpar
    do ihp = 0, maxpar
      do ipn = 0, maxpar
        ip = ipp + ipn
        if (ip == 0 .or. ip > maxpar) cycle
        do ihn = 0, maxpar
          ih = ihp + ihn
          if (ih == 0 .or. ih > maxpar) cycle
          feedph = xspopph2(Zcomp, Ncomp, nex, ipp, ihp, ipn, ihn)
          if (feedph <= 1.e-10) cycle
          sumph = 0.
          if (mpreeqmode == 2) then
            omegaph = phdens2(Zcomp, Ncomp, ipp, ihp, ipn, ihn, gsp, gsn, Exinc, Efermi, surfwell)
          else
            p0 = ip
            ppi0 = ipp
            hpi0 = ihp
            pnu0 = ipn
            hnu0 = ihn
            call exchange2(Zcomp, Ncomp)
            do p = p0, maxpar
              do ppi = ppi0, maxpar
                hpi = hpi0 + ppi - ppi0
                do pnu = pnu0, maxpar
                  hnu = hnu0 + pnu - pnu0
                  h = hpi + hnu
                  if (ppi + pnu == p .and. h <= maxpar) call lifetime2(ppi, hpi, pnu, hnu)
                enddo
              enddo
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
            zejec = parZ(type)
            nejec = parN(type)
            if (mpreeqmode == 2) then
              if (ipp - zejec < 0) cycle
              if (ipn - nejec < 0) cycle
            endif
!
! locate    : subroutine to find value in ordered table
!
            do nexout = Nlast(Zix, Nix, 0) + 1, nexmax(type)
              dEx = deltaEx(Zix, Nix, nexout)
              Eex = Ex(Zix, Nix, nexout)
              Eo = Exinc - Eex - S(Zcomp, Ncomp, type)
              call locate(egrid, ebegin(type), eend(type), Eo, nen)
              if (mpreeqmode == 2) then
                if (omegaph > 0.) then
                  gsp = gp(Zix, Nix)
                  gsn = gn(Zix, Nix)
                  if (flaggshell) then
                    damp = ignatyuk(Zix, Nix, Eex, 0) / alev(Zix, Nix)
                    gsp = gsp * damp
                    gsn = gsn * damp
                  endif
                  omegap1h = phdens2(Zix, Nix, ipp - zejec, ihp, ipn - nejec, ihn, gsp, gsn, Eex, Efermi, surfwell)
                  EoplusS = Exinc - Eex
                  omega1p = phdens2(Zix, Nix, zejec, 0, nejec, 0, gsp, gsn, EoplusS, Efermi, surfwell)
                  proba = omega1p * omegap1h / omegaph / (ipp + ipn)
                  Tswave = Tjl(type, nen, 1, 0)
                  Pescape = proba * Tswave
                  if (nexout == nexmax(type)) then
                    Exm = Exinc + 0.5 * dExinc - S(Zcomp, Ncomp, type)
                    Exmin = Ex(Zix, Nix, nexout) - 0.5 * dEx
                    dEx = Exm - Exmin
                  endif
                  term(type, nexout) = feedph * Pescape * dEx
                  sumterm = sumterm + term(type, nexout)
                endif
              else
                do ppi = ppi0, maxpar
                  if (ppi - zejec < 0) cycle
                  hpi = hpi0 + ppi - ppi0
                  do pnu = pnu0, maxpar
                    hnu = hnu0 + pnu - pnu0
                    if (pnu - nejec < 0) cycle
                    h = hpi + hnu
                    if (h > maxpar) cycle
                    factor(type, nexout, ppi, pnu) = feedph * Spre(ppi, hpi, pnu, hnu) * &
                      wemission2(type, ppi, hpi, pnu, hnu, nen) * dEx
                    term(type, nexout) = term(type, nexout) + factor(type, nexout, ppi, pnu)
                  enddo
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
                  do ppi = ppi0, maxpar
                    do pnu = pnu0, maxpar
                      factor(type, nexout, ppi, pnu) = factor(type, nexout, ppi, pnu) * feedph / sumterm
                    enddo
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
            zejec = parZ(type)
            nejec = parN(type)
            if (mpreeqmode == 2) then
              if (ipp - zejec < 0) cycle
              if (ipn - nejec < 0) cycle
            endif
            do nexout = Nlast(Zix, Nix, 0) + 1, nexmax(type)
              if (mpreeqmode == 2) then
                xspopph2(Zix, Nix, nexout, ipp - zejec, ihp, ipn - nejec, ihn) = &
                  xspopph2(Zix, Nix, nexout, ipp - zejec, ihp, ipn - nejec, ihn) + term(type, nexout)
              else
                do ppi = ppi0, maxpar
                  if (ppi - zejec < 0) cycle
                  hpi = hpi0 + ppi - ppi0
                  do pnu = pnu0, maxpar
                    hnu = hnu0 + pnu - pnu0
                    if (pnu - nejec < 0) cycle
                    h = hpi + hnu
                    if (h > maxpar) cycle
                    xspopph2(Zix, Nix, nexout, ppi - zejec, hpi, pnu - nejec, hnu) &
                      = xspopph2(Zix, Nix, nexout, ppi - zejec, hpi, pnu - nejec, hnu) + factor(type, nexout, ppi, pnu)
                  enddo
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
            if (ip <= maxpar - 1 .and. ih <= maxpar - 1) then
              if (ipp <= maxpar - 1 .and. ihp <= maxpar - 1) xspopph2(Zcomp, Ncomp, nex, ipp + 1, ihp + 1, ipn, ihn) = &
                xspopph2(Zcomp, Ncomp, nex, ipp + 1, ihp + 1, ipn, ihn) + 0.5 * (feedph - sumph)
              if (ipn <= maxpar - 1 .and. ihn <= maxpar - 1) xspopph2(Zcomp, Ncomp, nex, ipp, ihp, ipn + 1, ihn + 1) = &
                xspopph2(Zcomp, Ncomp, nex, ipp, ihp, ipn + 1, ihn + 1) + 0.5 * (feedph - sumph)
            endif
          endif
        enddo
      enddo
    enddo
  enddo
!
! ************************ Normalization *******************************
!
  Dmulti(nex) = summpe / xspopex(Zcomp, Ncomp, nex)
  xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) - summpe
  preeqpopex(Zcomp, Ncomp, nex) = preeqpopex(Zcomp, Ncomp, nex) - summpe
  return
end subroutine multipreeq2
! Copyright A.J. Koning 2021
