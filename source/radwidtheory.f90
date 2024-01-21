subroutine radwidtheory(Zcomp, Ncomp, E)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Theoretical calculation of total radiative width
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
!   numJ          ! maximum J - value
! Variables for gamma rays
!   gammax        ! number of l - values for gamma multipolarity
! Variables for level density
!   ldmodel       ! level density model
! Variables to normalize compound nucleus cross section
!   pardif        ! difference between target and compound nucleus parity
! Variables for nuclides
!   strucexist    ! flag to state whether structure info for this nucleus exists
! Constants
!   parspin       ! spin of particle
! Variables for levels
!   edis          ! energy of level
!   jdis          ! spin of level
!   parlev        ! parity of level
! Variables for resonance parameters
!   D0theo        ! mean s - wave resonance spacing
!   D1theo        ! mean p - wave resonance spacing
!   Eavres        ! average resonance energy
!   gamgamth      ! theoretical total radiative width
!   swaveth       ! theoretical strength function for s - wave
! Variables for level density
!   Nlast         ! last discrete level
! Variables for masses
!   S             ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer   :: irad              ! variable to indicate M(=0) or E(=1) radiation
  integer   :: Irspin2           ! 2 * residual spin
  integer   :: Irspin2beg        ! 2 * start of residual spin summation
  integer   :: Irspin2end        ! 2 * end of residual spin summation
  integer   :: J2                ! 2 * J
  integer   :: J2b               ! 2 * start of J summation
  integer   :: J2e               ! 2 * end of J summation
  integer   :: l                 ! multipolarity
  integer   :: l2                ! 2 * l
  integer   :: l2beg             ! 2 * start of l summation
  integer   :: l2end             ! 2 * end of l summation
  integer   :: ldmod             ! level density model
  integer   :: ll                ! angular momentum
  integer   :: modl              ! help variable
  integer   :: Ncomp             ! neutron number index for compound nucleus
  integer   :: nex               ! excitation energy bin of compound nucleus
  integer   :: nexgam            ! maximum excitation energy bin for gamma normalization
  integer   :: nexout            ! energy index for outgoing energy
  integer   :: Niter             ! counter
  integer   :: NL                ! last discrete level
  integer   :: Pprime            ! parity
  integer   :: Pprimebeg         ! start of residual parity summation
  integer   :: Pprimeend         ! end of residual parity summation
  integer   :: tpar              ! target parity
  integer   :: Zcomp             ! proton number index for compound nucleus
  real(sgl) :: dE1               ! help variable
  real(sgl) :: dE2               ! help variable
  real(sgl) :: dEx               ! excitation energy bin for population arrays
  real(sgl) :: E                 ! incident energy
  real(sgl) :: Egamma            ! gamma energy
  real(sgl) :: Exgam(0:10*numex) ! excitation energy for gamma normalization
  real(sgl) :: Exmid             ! help variable
  real(sgl) :: Exmin             ! help variable
  real(sgl) :: Explus            ! help variable
  real(dbl) :: factor            ! help variable
  real(sgl) :: fstrength         ! gamma ray strength function
  real(sgl) :: Rspin             ! residual spin
  real(sgl) :: Sgamma            ! gamma transmission coefficient/2 pi
  real(sgl) :: Sgamsum           ! sum over gamma-ray transmission coefficients/2 pi  (gamma strength function)
  real(sgl) :: Sn                ! neutron separation energy + incident energy
  real(sgl) :: tspin             ! target spin
  real(dbl) :: density           ! level density
  real(dbl) :: r1log             ! help variable
  real(dbl) :: r2log             ! help variable
  real(dbl) :: r3log             ! help variable
  real(dbl) :: rho               ! integrated level density
  real(dbl) :: rho1              ! help variable
  real(dbl) :: rho2              ! help variable
  real(dbl) :: rho3              ! help variable
!
! We calculate the theoretical total radiative width, for comparison with the experimental value and for possible later
! normalization of the gamma-ray strength functions.
!
! *************** Initialization of excitation energies ****************
!
! If the gamma normalization factor has been given in the input, we do not need to calculate it.
!
  Sn = S(Zcomp, Ncomp, 1) + E
  NL = Nlast(Zcomp, Ncomp, 0)
  do nex = 0, NL
    if (edis(Zcomp, Ncomp, nex) > Sn) then
      nexgam = nex - 1
      goto 30
    endif
    Exgam(nex) = edis(Zcomp, Ncomp, nex)
  enddo
  nexgam = 10 * numex
  dEx = (Sn - Exgam(NL)) / (nexgam - NL)
  do nex = NL + 1, nexgam
    Exgam(nex) = Exgam(NL) + (nex - NL) * dEx
  enddo
!
! ********************** Loop over quantum numbers *********************
!
! Specify the summation boundaries for the neutron channel.
!
30 if ( .not. strucexist(Zcomp, Ncomp + 1)) call levels(Zcomp, Ncomp + 1)
  tspin = jdis(Zcomp, Ncomp + 1, 0)
  do Niter = 0, 100
    do ll = 0, 1
      if (ll == 0) then
        tpar = parlev(Zcomp, Ncomp + 1, 0)
        J2b = int(abs(2. * (tspin - parspin(1))))
        J2e = int(2. * (tspin + parspin(1)))
      else
!
! Change of parity for p-wave
!
        tpar = - parlev(Zcomp, Ncomp + 1, 0)
        J2b = int(2. * (tspin - 1.5))
        if (J2b < 0) J2b = J2b + 2
        if (J2b < 0) J2b = J2b + 2
        J2e = int(2. * (tspin + 1.5))
      endif
      Sgamsum = 0.
      ldmod = ldmodel(Zcomp, Ncomp)
!
! Sum over total angular momentum J (J2) of compound nucleus
!
      do J2 = J2b, J2e, 2
!
!  Sum over outgoing excitation energie
!
        do nexout = 0, nexgam
          Exmid = Exgam(nexout)
          Egamma = Sn - Exmid
          if (Egamma <= 0.) cycle
!
! Set begin and end energies for level density integration
!
          if (nexout > NL) then
            Exmin = max(Exmid - 0.5 * dEx, 0.)
            Explus = Exmid + 0.5 * dEx
            dE1 = Exmid - Exmin
            dE2 = Explus - Exmid
          endif
!
! Initialization of summations
!
! For discrete states, the begin and end points of the residual spin/parity summation are both set equal to the residual discrete
! level spin/parity.
!
          if (nexout <= NL) then
            Pprimebeg = parlev(Zcomp, Ncomp, nexout)
            Pprimeend = Pprimebeg
            Irspin2beg = int(2. * jdis(Zcomp, Ncomp, nexout))
            Irspin2end = Irspin2beg
            rho = 1.
          else
!
! For the continuum, the begin and end points of the residual spin/parity summation are set to the maximally accessible values.
!
            Pprimebeg = - 1
            Pprimeend = 1
            Irspin2beg = mod(J2, 2)
            Irspin2end = J2 + 2 * gammax
            Irspin2end = min(Irspin2end, 2 * numJ)
          endif
!
!  Sum over residual parity
!
          do Pprime = Pprimebeg, Pprimeend, 2
            pardif = abs(tpar - Pprime) / 2
!
! rho is the level density integrated over the bin (Explus - Exmin).
! Note that we use logarithmic integration for more precision.
!
! Sum over residual spin
!
            do Irspin2 = Irspin2beg, Irspin2end, 2
              Rspin = 0.5 * Irspin2
              if (nexout > NL) then
                rho1 = real(density(Zcomp, Ncomp, Exmin, Rspin, Pprime, 0, ldmod)) * (1. + 1.d-10)
                rho2 = real(density(Zcomp, Ncomp, Exmid, Rspin, Pprime, 0, ldmod))
                rho3 = real(density(Zcomp, Ncomp, Explus, Rspin, Pprime, 0, ldmod)) * (1. + 1.d-10)
                r1log = log(rho1)
                r2log = log(rho2)
                r3log = log(rho3)
                if (r2log /= r1log .and. r2log /= r3log) then
                  rho = (rho1 - rho2) / (r1log - r2log) * dE1 + (rho2 - rho3) / (r2log - r3log) * dE2
                else
                  rho = rho2 * (dE1 + dE2)
                endif
              endif
!
! Sum over l of outgoing channel
!
              l2beg = max(abs(J2 - Irspin2), 2)
              l2end = min(J2 + Irspin2, 2 * gammax)
              do l2 = l2beg, l2end, 2
!
! Multipole radiation selection rules (irad=0: M-transition, irad=1: E-transition)
!
                l = l2 / 2
                modl = mod(l, 2)
                if (pardif == modl) then
                  irad = 1
                else
                  irad = 0
                endif
!
! Create sum for normalization
!
                Sgamma = (Egamma **(l2 + 1)) * fstrength(Zcomp, Ncomp, E, Egamma, irad, l)
                Sgamsum = Sgamsum + rho * Sgamma
              enddo
            enddo
          enddo
        enddo
      enddo
!
! ************** Calculation of theoretical values ********************
!
      if (E == Eavres) then
        if (ll == 0) then
          gamgamth(Zcomp, Ncomp, 0) = Sgamsum * D0theo(Zcomp, Ncomp)
          swaveth(Zcomp, Ncomp) = Sgamsum
        else
          gamgamth(Zcomp, Ncomp, 1) = Sgamsum * D1theo(Zcomp, Ncomp)
        endif
      endif
    enddo
    if (flaggnorm .and. gamgam(Zcomp,Ncomp) > 0.) then
      factor = gamgamth(Zcomp,Ncomp,0) / gamgam(Zcomp,Ncomp)
      if (abs(factor-1.) >= 0.0001) then
        ftable(Zcomp,Ncomp,1,1) = ftable(Zcomp,Ncomp,1,1) / factor
        call gammapar(Zcomp, Ncomp)
        cycle
      else
        exit
      endif
    else
      exit
    endif
  enddo
  return
end subroutine radwidtheory
! Copyright A.J. Koning 2021
