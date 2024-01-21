subroutine compound(Zcomp, Ncomp, nex, J2, parity)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Hauser-Feshbach model for multiple emission
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
!   numJ           ! maximum J - value
! Variables for numerics
!   transeps       ! absolute limit for transmission coefficient
! Variables for compound reactions
!   flagfullhf     ! flag for full spin dependence of transmission coefficie
! Variables for output
!   flagpop        ! flag for output of population
! Variables for fission
!   flagfission    ! flag for fission
! Variables for excitation energy grid
!   Exmax          ! maximum excitation energy for residual nucleus
!   maxJ           ! maximal J - value
!   nexmax         ! maximum excitation energy bin for residual nucleus
! Variables for multiple emission
!   fisfeedex      ! fission contribution from excitation energy bin
! Variables for multiple emission
!   Dmulti         ! depletion factor for multiple preequilibrium
!   mcontrib       ! contribution to emission spectrum
!   xsfeed         ! cross section from compound to residual nucleus
!   xspartial      ! emitted cross section flux per energy bin
! Variables for compound nucleus from target
!   dExinc         ! excitation energy bin for mother nucleus
!   Exinc          ! excitation energy of entrance bin
! Variables for incident channel
!   partdecay      ! total decay per particle
!   popdecay       ! decay from population
!   xspop          ! population cross section
!   xspopex        ! population cross section summed over spin and parity
!   xspopexP       ! population cross section per parity
!   xspopnuc       ! population cross section per nucleus
!   xspopnucP      ! population cross section per nucleus per parity
! Variables to prepare information for initial compound nucleus
!   denomhf        ! denominator for compound nucleus formula
!   enumhf         ! enumerator for compound nucleus formula
!   feed           ! feeding term for compound nucleus
! Variables for energy grid, level densities and transmission coefficients
!   lmaxhf         ! maximal l - value for transmission coefficients
!   rho0           ! integrated level density
!   Tgam           ! gamma transmission coefficients
!   Tjlnex         ! transmission coefficients as a function of particle type, energy,
!   Tlnex          ! transmission coefficients as a function of particle type, energy
! Variables for fission transmission coefficients
!   tfis           ! fission transmission coefficient for Hill - Wheeler magnitude
!   tfisdown       ! fission transmission coefficients
!   tfisup         ! fission transmission coefficients
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   Zindex         ! charge number index for residual nucleus
! Constants
!   parspin        ! spin of particle
!   spin2          ! 2  spin of particle
! Variables for levels
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for fission parameters
!   nfisbar        ! number of fission barrier parameters
! Variables for level density
!   Nlast          ! last discrete level
!
! *** Declaration of local data
!
  implicit none
  integer   :: iloop       ! loop counter
  integer   :: Ir          ! residual spin
  integer   :: irad        ! variable to indicate M(=0) or E(=1) radiation
  integer   :: Irspin2     ! 2 * residual spin
  integer   :: Irspin2beg  ! 2 * start of residual spin summation
  integer   :: Irspin2end  ! 2 * end of residual spin summation
  integer   :: J           ! spin of level
  integer   :: J2          ! 2 * J
  integer   :: J2minI2     ! help variable
  integer   :: J2plusI2    ! help variable
  integer   :: J2res       ! help variable
  integer   :: jj2prime    ! 2 * j'
  integer   :: jj2primebeg ! 2 * start of j' summation
  integer   :: jj2primeend ! 2 * end of j' summation
  integer   :: l2maxhf     ! 2 * lmaxhf
  integer   :: l2prime     ! 2 * l'
  integer   :: lb          ! help variable
  integer   :: lprime      ! 2 * l
  integer   :: lprimebeg   ! start of l summation
  integer   :: lprimeend   ! end of l summation
  integer   :: Ncomp       ! neutron number index for compound nucleus
  integer   :: nex         ! excitation energy bin of compound nucleus
  integer   :: nexout      ! energy index for outgoing energy
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: NL          ! last discrete level
  integer   :: pardif2     ! difference between residual and compound nucleus parity
  integer   :: parity      ! parity
  integer   :: parspin2    ! 2 * particle spin
  integer   :: Pprime      ! parity
  integer   :: Pprimebeg   ! start of residual parity summation
  integer   :: Pprimeend   ! end of residual parity summation
  integer   :: pspin2      ! 2 * spin of particle
  integer   :: type        ! particle type
  integer   :: updown2     ! spin index for transmission coefficient
  integer   :: Zcomp       ! proton number index for compound nucleus
  integer   :: Zix         ! charge number index for residual nucleus
  real(sgl) :: dE1         ! help variable
  real(sgl) :: dE2         ! help variable
  real(sgl) :: Exmin       ! help variable
  real(sgl) :: Explus      ! help variable
  real(sgl) :: fisfeed     ! cross section from compound nucleus to fission
  real(sgl) :: s2plus1     ! 2 * particle spin + 1
  real(dbl) :: factor      ! multiplication factor
  real(dbl) :: fiscontr    ! fission contribution
  real(dbl) :: fiscontr1   ! fission contribution
  real(dbl) :: fiscontr2   ! fission contribution
  real(dbl) :: leftover    ! remaining cross section flux after decay (should be close  to zero)
  real(dbl) :: logdtfd     ! help variable
  real(dbl) :: logdtfu     ! help variable
  real(dbl) :: rho         ! integrated level density
  real(dbl) :: sumIP       ! compound contribution summed over residual spin and parity
  real(dbl) :: sumIPE      ! compound contribution summed over residual spin and parity  and energy
  real(dbl) :: tf          ! help variable
  real(dbl) :: tfd         ! help variable
  real(dbl) :: tfu         ! help variable
  real(dbl) :: total       ! help variable
  real(dbl) :: totalrho    ! help variable
!
! **************************** Initialization **************************
!
  J = J2 / 2
  denomhf = 0.
!
! ********************* Loop over outgoing channels ********************
!
! Note that this whole subroutine is performed inside a loop over excitation energy bins nex, compound spin J and parity P.
! The loop over all outgoing particles, energies, residual spins and parities is done twice, and is designated by iloop.
!
! iloop=1: Loop to determine the partial and total widths of the Hauser-Feshbach formula.
!
! iloop=2: Loop to perform the actual compound nucleus calculation, using the ratio partial width/total width.
!
  do iloop = 1, 2
!
! 1. Fission
!
    if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) then
!
! iloop=1: The fission width for the mother excitation energy bin with spin J and parity P is calculated and added to the total
!   decay width of the mother bin. We use a logarithmic integration of the fission transmission coefficients.
!
      if (iloop == 1) then
        tfd = max(tfisdown(J, parity), transeps)
        tf = max(tfis(J, parity), transeps)
        tfu = max(tfisup(J, parity), transeps)
        Explus = min(Exmax(Zcomp, Ncomp), Exinc + 0.5 * dExinc)
        Exmin = max(Exinc - 0.5 * dExinc, 0.)
        dE1 = Exinc - Exmin
        dE2 = Explus - Exinc
        logdtfd = log(tf) - log(tfd)
        logdtfu = log(tfu) - log(tf)
        if (logdtfd == 0.) then
          fiscontr1 = tf * dE1
        else
          fiscontr1 = (tf - tfd) / logdtfd * dE1
        endif
        if (logdtfu == 0.) then
          fiscontr2 = tf * dE2
        else
          fiscontr2 = (tfu - tf) / logdtfu * dE2
        endif
        if (Explus > Exmin) then
          fiscontr = (fiscontr1 + fiscontr2) / (Explus - Exmin)
        else
          fiscontr = 0.
        endif
        if (fiscontr <= 10. * transeps) fiscontr = 0.
        denomhf = fiscontr
      else
!
! iloop=2: The fission contribution for the mother excitation energy bin with spin J and parity P is calculated and added to the
!   various feeding arrays.
!
        fisfeed = 0.
        if (denomhf /= 0.) then
          fisfeed = real(feed * fiscontr)
          xsfeed(Zcomp, Ncomp, - 1) = xsfeed(Zcomp, Ncomp, - 1) + fisfeed
          fisfeedex(Zcomp, Ncomp, nex) = fisfeedex(Zcomp, Ncomp, nex) + fisfeed
        endif
      endif
    endif
!
! Photon and particle channels
!
    if (iloop == 1) then
      do Pprime = - 1, 1, 2
        do Ir = 0, numJ
          do type = 0, 6
            do nexout = 0, nexmax(type)
              enumhf(type, nexout, Ir, Pprime) = 0.
            enddo
          enddo
        enddo
      enddo
    endif
    do type = 0, 6
      if (parskip(type)) cycle
      if (iloop == 1 .and. type == 6 .and. denomhf == 0.) cycle
      parspin2 = int(2. * parspin(type))
      s2plus1 = parspin2 + 1.
      pspin2 = spin2(type)
      Zix = Zindex(Zcomp, Ncomp, type)
      Nix = Nindex(Zcomp, Ncomp, type)
      NL = Nlast(Zix, Nix, 0)
!
! This loop is over all discrete levels and continuum bins of the final nucleus.
!
      sumIPE = 0.
!
! Loop over outgoing excitation energies
!
      do nexout = 0, nexmax(type)
        if (nexout == 0 .and. NL == 0) cycle
!
! Initialization of summations
!
! For discrete states, the begin and end points of the residual spin/parity summation are both set equal to the residual discrete
! level spin/parity.
!
        if (nexout <= NL) then
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
          J2res = J2 + parspin2
          Irspin2beg = mod(J2res, 2)
          Irspin2end = 2 * maxJ(Zix, Nix, nexout)
        endif
        l2maxhf = 2 * lmaxhf(type, nexout)
        if (iloop == 1) then
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
              rho = rho0(type, nexout, Ir, Pprime)
!
! The Hauser-Feshbach formula contains the following triangular relations:
! |J-I| < j < J+I
! |j-s| < l < j+s
! Computationally it is more economic to do the l-summation first.
! This can be done with the following equivalent triangular relations:
! ||J-I|-s| < l < |J+I+s|
! max(|J-I|,|l-s|) < j < min((J+I),l+s)
!
              total = 0.
              J2minI2 = abs(J2 - Irspin2)
              J2plusI2 = J2 + Irspin2
              lprimebeg = abs(J2minI2 - parspin2) / 2
              lprimeend = abs(J2plusI2 + parspin2) / 2
              lprimeend = min(lprimeend, l2maxhf / 2)
!
! Loop over l of outgoing channel
!
! We include photons as a special case, with the multipole radiation selection rules (irad=0: M-transition, irad=1: E-transition)
!
! 1. Photons
!
              if (type == 0) then
                do lprime = lprimebeg, lprimeend
                  irad = 1 - abs(mod(lprime, 2) - pardif2)
                  total = total + Tgam(nexout, lprime, irad)
                enddo
              else
!
! 2. Particles
!
! Loop over j (jj2) of outgoing channel
!
! If the parity of the residual nucleus is equal (unequal) to the parity of compound nucleus,i.e. pardif2=0(1),
! the l-value must be even (odd).
!
! A choice is made between averaged and full spin dependence of the transmission coefficients in the Hauser-Feshbach model (only
! significant in case of very large spin-orbit terms in the OMP).
!
                lb = lprimebeg + abs(mod(lprimebeg, 2) - pardif2)
                if (flagfullhf) then
                  do lprime = lb, lprimeend, 2
                    l2prime = 2 * lprime
                    jj2primebeg = max(J2minI2, abs(l2prime - parspin2))
                    jj2primeend = min(J2plusI2, l2prime + parspin2)
                    do jj2prime = jj2primebeg, jj2primeend, 2
                      updown2 = (jj2prime - l2prime) / pspin2
                      total = total + Tjlnex(type, nexout, updown2, lprime)
                    enddo
                  enddo
                else
                  do lprime = lb, lprimeend, 2
                    total = total + Tlnex(type, nexout, lprime)
                  enddo
                  total = s2plus1 * total
                endif
              endif
!
! iloop=1: The partial decay widths are stored in enumhf and also added to the total decay width denomhf.
!
              totalrho = rho * total
              enumhf(type, nexout, Ir, Pprime) = totalrho
              denomhf = denomhf + totalrho
            enddo
          enddo
        else
!
! ** Populate the outgoing channels using the compound nucleus formula *
!
! iloop=2: Increment of the the population arrays of the residual nuclei and the mcontrib array, which will be used for the
!   interpolation of the compound emission spectra.
!
          sumIP = 0.
          do Pprime = Pprimebeg, Pprimeend, 2
            do Irspin2 = Irspin2beg, Irspin2end, 2
              Ir = Irspin2 / 2
              factor = real(feed * enumhf(type, nexout, Ir, Pprime))
              xspop(Zix, Nix, nexout, Ir, Pprime) = xspop(Zix, Nix, nexout, Ir, Pprime) + factor
              if (flagpop) then
                xspopnucP(Zix, Nix, Pprime) = xspopnucP(Zix, Nix, Pprime) + factor
                xspopexP(Zix, Nix, nexout, Pprime) = xspopexP(Zix, Nix, nexout, Pprime) + factor
                popdecay(type, nexout, Ir, Pprime) = popdecay(type, nexout, Ir, Pprime) + factor
                partdecay(type, Pprime) = partdecay(type, Pprime) + factor
                partdecaytot(type) = partdecaytot(type) + factor
              endif
              sumIP = sumIP + factor
            enddo
          enddo
          xspopex(Zix, Nix, nexout) = xspopex(Zix, Nix, nexout) + sumIP
          xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) - sumIP
          mcontrib(type, nex, nexout) = mcontrib(type, nex, nexout) + sumIP
          sumIPE = sumIPE + sumIP
        endif
      enddo
      if (iloop == 2) then
        xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) + sumIPE
        xspartial(type, nex) = xspartial(type, nex) + sumIPE
        xsfeed(Zcomp, Ncomp, type) = xsfeed(Zcomp, Ncomp, type) + sumIPE
      endif
    enddo
    if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0 .and. iloop == 2) &
  &   xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) - fisfeed
!
! ** Create feeding term for compound nucleus decay in the second loop *
!
    if (iloop == 1) then
      if (denomhf /= 0.) then
        feed = (1. - Dmulti(nex)) * xspop(Zcomp, Ncomp, nex, J, parity) / denomhf
      else
!
! Prevent trapping of cross section in the continuum.
! Any flux that is left in the compound nucleus bin is equally distributed over the discrete states.
!
        NL = Nlast(Zcomp, Ncomp, 0)
        feed = 0.
        leftover = xspop(Zcomp, Ncomp, nex, J, parity) / (NL + 1.)
        xspartial(0, nex) = xspartial(0, nex) + xspop(Zcomp, Ncomp, nex, J, parity)
        do nexout = 0, NL
          xspopex(Zcomp, Ncomp, nexout) = xspopex(Zcomp, Ncomp, nexout) + leftover
          Ir = int(jdis(Zcomp, Ncomp, nexout))
          Pprime = parlev(Zcomp, Ncomp, nexout)
          xspop(Zcomp, Ncomp, nexout, Ir, Pprime) = xspop(Zcomp, Ncomp, nexout, Ir, Pprime) + leftover
          popdecay(0, nexout, Ir, Pprime) = popdecay(0, nexout, Ir, Pprime) + leftover
          mcontrib(0, nex, nexout) = mcontrib(0, nex, nexout) + leftover
        enddo
      endif
    endif
  enddo
  return
end subroutine compound
! Copyright A.J. Koning 2021
