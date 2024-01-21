subroutine astrotarget
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Compound reaction for many target states
!
! Author    : Stephane Hilaire
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
! Variables for output
!   flagcheck      ! flag for output of numerical checks
! Variables for fission
!   flagfission    ! flag for fission
! Variables for compound reactions
!   wmode          ! designator for width fluctuation model
! Variables for main input
!   k0             ! index of incident particle
!   Ltarget        ! excited level of target
! Variables for basic reaction
!   flagastro      ! flag for calculation of astrophysics reaction rate
! Variables for energies
!   flagwidth      ! flag for width fluctuation calculation
! Variables for excitation energy grid
!   maxex          ! maximum excitation energy bin for residual nucleus
!   maxJ           ! maximal J - value
! Variables for incident channel
!   contrib        ! contribution to emission spectrum
!   xsbinary       ! cross section from initial compound to residual nucleus
!   xscompcont     ! compound cross section for continuum
!   xspop          ! population cross section
!   xspopex        ! population cross section summed over spin and parity
!   xspopnuc       ! population cross section per nucleus
! Variables for compound nucleus from target
!   tnumi          ! counter for width fluctuation calculation
!   tnumo          ! counter for width fluctuation calculation
!   Wab            ! width fluctuation factor
! Variables to normalize compound nucleus cross section
!   CNfactor       ! factor for compound nucleus cross section: pi / [ k **2 (2s + 1)(2I + 1) ]
!   J2beg          ! begin of J summation
!   pardif         ! difference between target and compound nucleus parity
! Variables to prepare information for initial compound nucleus
!   denomhf        ! denominator for compound nucleus formula
!   enumhf         ! enumerator for compound nucleus formula
!   feed           ! feeding term for compound nucleus
!   tnum           ! counter for width fluctuation calculation
!   tNinc        ! counter for width fluctuation calculation
! Variables for energy grid, level densities and transmission coefficients
!   lmaxhf         ! maximal l - value for transmission coefficients
!   rho0           ! integrated level density
! Variables for fission transmission coefficients
!   tfis           ! fission transmission coefficient for Hill - Wheeler magnitude
! Variables for nuclides
!   Nindex         ! neutron number index for residual nucleus
!   parskip        ! logical to skip outgoing particle
!   Zindex         ! charge number index for residual nucleus
! Constants
!   cparity        ! parity (character)
!   parspin        ! spin of particle
! Variables for levels
!   jdis           ! spin of level
!   parlev         ! parity of level
! Variables for fission parameters
!   nfisbar        ! number of fission barrier parameters
! Variables for level density
!   Nlast          ! last discrete level
! Variables for astro
!   Tastroinc      ! transmission coefficient for incident channel (Astrophysic
!   Tastroout      ! transmission coefficient for outgoing channel (Astrophysic
!
! *** Declaration of local data
!
  implicit none
  logical   :: elas1       ! logical for elastic channel
  logical   :: elas2       ! logical for elastic channel
  logical   :: elastic     ! designator for elastic channel
  integer   :: ielas       ! designator for elastic channel
  integer   :: Ir          ! residual spin
  integer   :: Irspin2     ! 2 * residual spin
  integer   :: Irspin2beg  ! 2 * start of residual spin summation
  integer   :: Irspin2end  ! 2 * end of residual spin summation
  integer   :: J           ! spin of level
  integer   :: J2          ! 2 * J
  integer   :: J2cnend     ! 2 * end of J summation of CN
  integer   :: J2res       ! help variable
  integer   :: l2maxhf     ! 2 * lmaxhf
  integer   :: Ncomp       ! neutron number index for compound nucleus
  integer   :: nexastro    ! energy index for astrophysics
  integer   :: nexout      ! energy index for outgoing energy
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: Nixtarget   ! neutron number index of target nucleus
  integer   :: NL          ! last discrete level
  integer   :: parity      ! parity
  integer   :: parspin2o   ! 2 * particle spin for outgoing channel
  integer   :: Pbeg        ! begin and end of parity summation
  integer   :: Pend        ! begin and end of parity summation
  integer   :: Pprime      ! parity
  integer   :: Pprimebeg   ! start of residual parity summation
  integer   :: Pprimeend   ! end of residual parity summation
  integer   :: Ptarget     ! parity of target state
  integer   :: spin2beg    ! 2 * start of target spin summation
  integer   :: spin2end    ! 2 * end of target spin summation
  integer   :: spin2target ! 2 * target spin
  integer   :: spintarget  ! target spin
  integer   :: type        ! particle type
  integer   :: Zcomp       ! proton number index for compound nucleus
  integer   :: Zix         ! charge number index for residual nucleus
  integer   :: Zixtarget   ! charge number index of target nucleus
  real(sgl) :: fluxsum     ! check for conservation of flux per P,J,j,l
  real(sgl) :: Tinc        ! transmission coefficients as a function of j and l
  real(sgl) :: Tout        ! transmission coefficients
  real(dbl) :: compterm    ! partial contribution to compound nucleus term
  real(dbl) :: factor1     ! help variable
  real(dbl) :: rho         ! integrated level density
  real(dbl) :: rhoel       ! level density
  real(dbl) :: rhoinc      ! integrated level density for target
  real(dbl) :: suminl      ! sum over inelastic channels
  real(dbl) :: sumIP       ! compound contribution summed over residual spin and parity
  real(dbl) :: sumIPas     ! sumIP for astrophysics
  real(dbl) :: sumIPE      ! compound contribution summed over residual spin and parity  and energy
  real(dbl) :: Wabinelastic! WFC factor
!
! ********************** Astrophysical situation ***********************
!
! The astrophysical calculation may involve also the excited states of the target located above the ground state of the nucleus
! or both below and above when dealing with a isomeric target.
! The temperature during the evolution of the universe or in stars is such that a nucleus can exist in various excited states.
! This would in principle imply that the loops over the quantum number must be performed for all possible inelastic channels.
! However, in practice these loops are reduced by a lumping approximation for lprime and jprime quantum numbers which means that
! we only loop over excited states number,spin and parity of the target nucleus.
!
  Zcomp = 0
  Ncomp = 0
  Zixtarget = Zindex(Zcomp, Ncomp, k0)
  Nixtarget = Nindex(Zcomp, Ncomp, k0)
  if (flagcheck .and. flagwidth) write(*, '(/"Flux check for astrophysical case"/)')
  do nexastro = 0, maxex(Zixtarget, Nixtarget)
    if (nexastro == Ltarget) cycle
    l2maxhf = 2 * lmaxhf(k0, nexastro)
    NL = Nlast(Zixtarget, Nixtarget, 0)
    if (nexastro <= NL) then
      Pbeg = parlev(Zixtarget, Nixtarget, nexastro)
      Pend = Pbeg
      spin2beg = int(2. * jdis(Zixtarget, Nixtarget, nexastro))
      spin2end = spin2beg
    else
      Pbeg = - 1
      Pend = 1
      spin2beg = mod(int(2 * jdis(Zixtarget, Nixtarget, 0)), 2)
      spin2end = 2 * maxJ(Zixtarget, Nixtarget, nexastro)
    endif
!
! Loop over target parity
!
    do Ptarget = Pbeg, Pend, 2
!
! Loop over target spin
!
      do spin2target = spin2beg, spin2end, 2
        spintarget = spin2target / 2
        rhoinc = rho0(k0, nexastro, spintarget, Ptarget)
        if (rhoinc == 0.) cycle
        J2cnend = int(2 * (lmaxhf(k0, nexastro) + parspin(k0) + spintarget))
        J2cnend = min(J2cnend, numJ)
!
! Loop over Compound Nucleus parity
!
        do parity = - 1, 1, 2
          pardif = abs(Ptarget - parity) / 2
!
! tfission    : subroutine for fission transmission coefficients
! astroprepare: subroutine to prepare information for astrophysical compound nucleus calculations
! widthprepare: subroutine for preparation of width fluctuation corrections
!
! Loop over Compound Nucleus total angular momentum J (J2)
!
          do J2 = J2beg, J2cnend, 2
            J = J2 / 2
            if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) call tfission(Zcomp, Ncomp, nexastro, J2, parity)
            call astroprepare(Zcomp, Ncomp, J2, parity, spin2target, Ptarget, nexastro)
            if (denomhf == 0.) cycle
            if (flagwidth) then
              Tinc = Tastroinc(1, nexastro, spintarget, Ptarget)
              rhoinc = Tastroinc(0, nexastro, spintarget, Ptarget)
              tnumi = 1
              if (Tinc < transeps) cycle
              call widthprepare
            else
              if (feed < transeps) cycle
            endif
            fluxsum = 0.
!
! widthfluc    : subroutine for width fluctuation correction pi/[ k**2 (2s+1)(2I+1) ]
!
! 1. Fission channel
!
            if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) then
              if (flagwidth) then
                tnumo = tnum
                if (wmode >= 1) then
                  call widthfluc(0)
                  factor1 = real(Tinc * tfis(J, parity) / denomhf * Wab)
                  fluxsum = fluxsum + factor1
                else
                  factor1 = real(Tinc * tfis(J, parity) / denomhf)
                endif
              else
                factor1 = real(feed * tfis(J, parity) / denomhf)
              endif
              xsbinary( - 1) = xsbinary( - 1) + CNfactor * (J2 + 1.) * factor1
            endif
!
! 2. Gamma and particle channels
!
            tnumo = tnum + 1
            do type = 0, 6
              if (type == 1) tnumo = tNinc
              if (parskip(type)) cycle
              parspin2o = int(2. * parspin(type))
              Zix = Zindex(Zcomp, Ncomp, type)
              Nix = Nindex(Zcomp, Ncomp, type)
              NL = Nlast(Zix, Nix, 0)
              sumIPE = 0.
!
! Loop over over outgoing excitation energies
!
              do nexout = 0, maxex(Zix, Nix)
                l2maxhf = 2 * lmaxhf(type, nexout)
                elas1 = (type == k0 .and. nexout == nexastro)
                if (nexout <= NL) then
                  Pprimebeg = parlev(Zix, Nix, nexout)
                  Pprimeend = Pprimebeg
                  Irspin2beg = int(2. * jdis(Zix, Nix, nexout))
                  Irspin2end = Irspin2beg
                else
                  Pprimebeg = - 1
                  Pprimeend = 1
                  J2res = J2 + parspin2o
                  Irspin2beg = mod(J2res, 2)
                  Irspin2end = J2res + l2maxhf
                  Irspin2end = min(Irspin2end, 2 * maxJ(Zix, Nix, nexout))
                endif
                sumIP = 0.
                sumIPas = 0.
!
! Loop over residual parity
!
                do Pprime = Pprimebeg, Pprimeend, 2
!
! Loop over residual spin
!
                  do Irspin2 = Irspin2beg, Irspin2end, 2
                    Ir = Irspin2 / 2
                    elas2 = (spintarget == Ir .and. Ptarget == Pprime)
                    elastic = (elas1 .and. elas2)
                    if (flagwidth) then
                      rho = Tastroout(0, type, nexout, Ir, Pprime)
                      if (rho < 1.0d-20) cycle
                      if (rho /= 0.) Tout = Tastroout(1, type, nexout, Ir, Pprime) / rho
                      if (type >= 1) tnumo = tnumo + 1
                      if (wmode >= 1) then
                        ielas = 0
                        call widthfluc(ielas)
                        Wabinelastic = Wab
                        if (elastic) then
                          ielas = 1
                          call widthfluc(ielas)
                          rhoel = rho
                        endif
                      endif
                      factor1 = real(Tinc * rho * Tout / denomhf * Wab)
!
! We avoid double counting for the elastic channel, take into account the rho*(rho-1) inelastic contribution of the lumped
!  (l,j) incident channels. For more explanations, please call 911
!
                      if (ielas == 1) then
                         factor1 = factor1 + real((rho - 1.) * Tinc * rho * Tout / denomhf * Wabinelastic)
                         if (rho > 0.) factor1 = factor1 / rho
                         suminl = real((rho - 1.) * Tinc * rho * Tout / denomhf * Wabinelastic)
                         if (rho > 0.) suminl = suminl / rho
                      endif
                      fluxsum = fluxsum + factor1
                    else
                      factor1 = feed * enumhf(type, nexout, Ir, Pprime) / denomhf
                      if (flagastro .and. elastic) then
                        sumIPas = sumIPas + CNfactor * (J2 + 1) * enumhf(type, nexout, Ir, Pprime) **2 / denomhf
                      endif
                    endif
                    compterm = CNfactor * (J2 + 1.) * factor1
                    xspop(Zix, Nix, nexout, Ir, Pprime) = xspop(Zix, Nix, nexout, Ir, Pprime) + compterm
                    sumIP = sumIP + compterm
                  enddo
                enddo
                xspopex(Zix, Nix, nexout) = xspopex(Zix, Nix, nexout) + sumIP
                if (nexout > NL) then
                  xscompcont(type) = xscompcont(type) + sumIP
                  contrib(type, nexout) = contrib(type, nexout) + sumIP
                endif
!
! Compound elastic scattering is excluded from the residual production cross sections
!
                if ( .not. flagastro) then
                  if ( .not. elastic) sumIPE = sumIPE + sumIP
                else
                  if (flagwidth) then
                    if ( .not. elastic) then
                      sumIPE = sumIPE + sumIP
                    else
                      sumIPE = sumIPE + suminl
                    endif
                  else
                    if (elastic) then
                      sumIPE = sumIPE + (sumIP - sumIPas)
                    else
                      sumIPE = sumIPE + sumIP
                    endif
                  endif
                endif
              enddo
              xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) + sumIPE
              xsbinary(type) = xsbinary(type) + sumIPE
            enddo
!
! ****** Check of flux conservation of transmission coefficients *******
!
! This check is included to test the stability of the width fluctuation calculation.
!
            if (flagcheck .and. flagwidth) then
              if (fluxsum == 0.) fluxsum = Tinc
              write(*, '(" Parity=", a1, "  J=", f4.1, " Tinc", es12.5, " Sum over outgoing channels=", es12.5, &
 &              "  Ratio=", f8.5, "  Rho=", 2g14.6)') cparity(parity), 0.5*J2, Tinc, fluxsum, Tinc / fluxsum, rhoinc, rhoel
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  return
end subroutine astrotarget
! Copyright A.J. Koning 2021
