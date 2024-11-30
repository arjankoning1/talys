subroutine comptarget
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Compound reaction for initial compound nucleus
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
!   sgl             ! single precision kind
!   dbl             ! double precision kind
! All global variables
!   numfact         ! number of terms for factorial logarithm
!   numhill         ! maximum number of Hill - Wheeler points
! Variables for numerics
!   transeps        ! absolute limit for transmission coefficient
!   xseps           ! limit for cross sections
! Variables for output
!   flagcheck       ! flag for output of numerical checks
!   flagdecay       ! flag for output of decay of each population bin
!   flagpop         ! flag for output of population
! Variables for compound reactions
!   adjustTJ        ! logical for energy-dependent TJ adjustment
!   ewfc            ! off - set incident energy for width fluctuation
!   flageciscomp    ! flag for compound nucleus calculation by ECIS
!   flagurr         ! flag for output of unresolved resonance parameters
!   lurr            ! maximal orbital angular momentum for URR
!   wmode           ! designator for width fluctuation model
! Variables for input energies
!   flaginitpop     ! flag for initial population distribution
! Variables for main input
!   k0              ! index of incident particle
!   Ltarget         ! excited level of target
!   Ninit           ! neutron number of initial compound nucleus
!   Zinit           ! charge number of initial compound nucleus
! Variables for gamma rays
!   fiso            ! correction factor for isospin forbidden transitions
! Variables for fission
!   flagfisout      ! flag for output of fission information
!   flagfission     ! flag for fission
! Variables for basic reaction
!   flagastro       ! flag for calculation of astrophysics reaction rate
! Variables for astrophysics
!   flagastrogs     ! flag for calculation of astrophysics reaction rate with
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for energies
!   flagcompang     ! flag for compound angular distribution calculation
!   flagwidth       ! flag for width fluctuation calculation
! Variables for excitation energy grid
!   Ex              ! excitation energy
!   fisfeedJP       ! fission contribution from excitation energy bin per J, P
!   maxex           ! maximum excitation energy bin for residual nucleus
!   maxJ            ! maximal J - value
! Variables for compound nucleus from target
!   dExinc          ! excitation energy bin for mother nucleus
!   Fnorm           ! multiplication factor
!   JmaxU           ! maximal total angular momentum
!   JminU           ! minimal total angular momentum
!   lmaxU           ! maximal orbital angular momentum
!   lminU           ! minimal orbital angular momentum
!   nulj            ! (l, j) number of degrees of freedom for URR calculation
!   Purrlj          ! (l, j) parity for URR calculation
!   tnumi           ! counter for width fluctuation calculation
!   tnumo           ! counter for width fluctuation calculation
!   Turrlj          ! transmission coefficient for URR calculation
!   Turrljinc       ! incident channel (l, j) transmission coefficient for URR ca
!   Wab             ! width fluctuation factor
!   xsbinarylj      ! cross section from initial compound to residual nucleus
! Variables for incident channel
!   cleg            ! compound nucleus Legendre coefficient
!   contrib         ! contribution to emission spectrum
!   lmaxinc         ! maximal l - value for transm. coeff. for incident channel
!   partdecay       ! total decay per particle
!   popdecay        ! decay from population
!   Tjlinc          ! transm. coeff. as a function of spin and l for inc. channel
!   xsbinary        ! cross section from initial compound to residual nucleus
!   xscompcont      ! compound cross section for continuum
!   xspop           ! population cross section
!   xspopex         ! population cross section summed over spin and parity
!   xspopexP        ! population cross section per parity
!   xspopnuc        ! population cross section per nucleus
!   xspopnucP       ! population cross section per nucleus per parity
! Variables to normalize compound nucleus cross section
!   CNfactor        ! factor for compound nucleus cross section: pi / [ k **2 (2s + 1)(2I + 1) ]
!   CNterm          ! compound nucleus formation cross section per spin
! Variables to normalize compound nucleus cross section
!   J2beg           ! begin of J summation
!   J2end           ! end of J summation
!   pardif          ! difference between target and compound nucleus parity
! Variables to prepare information for initial compound nucleus
!   denomhf         ! denominator for compound nucleus formula
!   enumhf          ! enumerator for compound nucleus formula
!   feed            ! feeding term for compound nucleus
!   tnum            ! counter for width fluctuation calculation
!   tNinc         ! counter for width fluctuation calculation
! Variables for energy grid, level densities and transmission coefficients
!   lmaxhf          ! maximal l - value for transmission coefficients
!   rho0            ! integrated level density
!   Tgam            ! gamma transmission coefficients
!   Tjlnex          ! transmission coefficients as a function of particle type, energy,
! Variables for fission transmission coefficients
!   tfis            ! fission transmission coefficient for Hill - Wheeler magnitude
!   tfisA           ! transmission coefficient for Hill - Wheeler magnitude
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   Nindex          ! neutron number index for residual nucleus
!   NN              ! neutron number of residual nucleus
!   parskip         ! logical to skip outgoing particle
!   targetP         ! parity of target
!   targetspin      ! spin of target
!   targetspin2     ! 2 * target spin
!   Zindex          ! charge number index for residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   cparity         ! parity (character)
!   fourpi          ! 4. * pi
!   nuc             ! symbol of nucleus
!   parN            ! neutron number of particle
!   parname         ! name of particle
!   parspin         ! spin of particle
!   parZ            ! charge number of particle
!   sgn             ! sign
!   spin2           ! 2 * spin of particle
! Variables for levels
!   jdis            ! spin of level
!   parlev          ! parity of level
! Variables for level density
!   Nlast           ! last discrete level
! Variables for fission parameters
!   nfisbar         ! number of fission barrier parameters
! Variables for astro
!   rhoastrotot     ! total level density for astrophysical case
!   Tastrotot       ! total transmission coefficient for astrophysical case
! Variables for initial compound nucleus
!   logfact         ! factorial logarithm
! Variables for preequilibrium
!   xsflux          ! cross section flux
!
! *** Declaration of local data
!
  implicit none
  character*132 :: key     ! keyword
  character(len=3)  :: massstring !  
  character(len=15) :: col(2*numJ+6)    ! header
  character(len=15) :: un(2*numJ+6)    ! units
  character(len=80) :: quantity   ! quantity
  character(len=6)  :: finalnuclide !
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=80)  :: popfile               ! population file
  character(len=13)  :: Estr
  logical   :: elastic     ! designator for elastic channel
  integer   :: Ares        ! mass number of residual nucleus
  integer   :: ielas       ! designator for elastic channel
  integer   :: ihill       ! counter for Hill-Wheeler magnitude
  integer   :: iphase      ! help variable
  integer   :: Ir          ! residual spin
  integer   :: irad        ! variable to indicate M(=0) or E(=1) radiation
  integer   :: Irspin2     ! 2 * residual spin
  integer   :: Irspin2beg  ! 2 * start of residual spin summation
  integer   :: Irspin2end  ! 2 * end of residual spin summation
  integer   :: J           ! spin of level
  integer   :: J2          ! 2 * J
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
  integer   :: LL          ! counter for l value
  integer   :: LLmax       ! maximal Legendre order
  integer   :: lprime      ! 2 * l
  integer   :: modl        ! help variable
  integer   :: Ncomp       ! neutron number index for compound nucleus
  integer            :: Z
  integer            :: A
  integer            :: Ncol
  integer   :: nex         ! excitation energy bin of compound nucleus
  integer   :: nexout      ! energy index for outgoing energy
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: NL          ! last discrete level
  integer   :: Nres        ! neutron number of residual nucleus
  integer   :: oddres      ! odd (1) or even (0) nucleus
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
  integer   :: Zres        ! charge number of residual nucleus
  real(sgl) :: Ablatt      ! Blatt-Biedenharn A-factor
  real(sgl) :: angfac      ! help variable
  real(sgl) :: clebsch     ! function for Clebsch-Gordan coefficients
  real(sgl) :: clin        ! Clebsch-Gordan coefficients
  real(sgl) :: clout       ! Clebsch-Gordan coefficients
  real(sgl) :: factor      ! help variable
  real(sgl) :: fluxsum     ! check for conservation of flux per P,J,j,l
  real(sgl) :: phase       ! phase
  real(sgl) :: ra1in       ! Racah coefficients
  real(sgl) :: ra1out      ! Racah coefficients
  real(sgl) :: ra2in       ! Racah coefficients
  real(sgl) :: ra2out      ! Racah coefficients
  real(sgl) :: racah       ! function for racah coefficients
  real(sgl) :: ratio       ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: rJ          ! help variable
  real(sgl) :: rjin        ! help variable
  real(sgl) :: rjout       ! help variable
  real(sgl) :: rlin        ! help variable
  real(sgl) :: rLL         ! help variable
  real(sgl) :: rlout       ! help variable
  real(sgl) :: Rspin       ! residual spin
  real(sgl) :: Tinc        ! transmission coefficients as a function of j and l  for the incident channel
  real(sgl) :: Tout        ! transmission coefficients
  real(dbl) :: compterm    ! partial contribution to compound nucleus term
  real(dbl) :: factor1     ! help variable
  real(dbl) :: fisterm     ! fission term
  real(dbl) :: rho         ! integrated level density
  real(dbl) :: sumIP       ! compound contribution summed over residual spin and parity
  real(dbl) :: sumIPas     ! sumIP for astrophysics
  real(dbl) :: sumIPE      ! compound contribution summed over residual spin and parity  and energy
  real(dbl) :: sumjl       ! compound contribution summed over residual j and l
!
! *** Check if initial compound nucleus calculation needs to be done ***
!
! This may occur when direct + pre-equilibrium reactions completely exhaust the reaction cross section for the binary channel, and
! when the initial system is a excitation energy population.
!
  if (flaginitpop) return
  if (xsflux <= xseps) return
!
! *************************** Initialization ***************************
!
! Initially, Wab is set to 1 (no width fluctuations).
!
  Wab = 1.
  parspin2i = int(2. * parspin(k0))
  pspin2i = spin2(k0)
!
! The level densities and transmission coefficients can be prepared before the nested loops over all quantum numbers.
!
! densprepare: subroutine to prepare energy grid, level density and transmission coefficient information for compound nucleus
!
  Zcomp = 0
  Ncomp = 0
  dExinc = 0.
!
! Optional adjustment factors
!
! isotrans: subroutine for correction factors for isospin forbidden transitions
!
  call isotrans(Zinit, Ninit)
!
! Output of initial population 
!
  if (flagpop) then
    Z = ZZ(Zcomp, Ncomp, 0)
    A = AA(Zcomp, Ncomp, 0)
    massstring='   '
    write(massstring,'(i3)') A
    finalnuclide=trim(nuc(Z))//adjustl(massstring)
    reaction='('//parsym(k0)//',x)'
    quantity='initial population'
    Estr=''
    write(Estr,'(es13.6)') Einc
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' of '//trim(finalnuclide)//' at '//Estr//' MeV'
    popfile='initial.pop'
    open (unit = 1, file = popfile, status = 'replace')
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.D0,0.D0,0,0)
    call write_real(2,'E-incident [MeV]',Einc)
    call write_residual(Z,A,finalnuclide)
    call write_real(2,'Excitation energy [MeV]',Exinc)
    quantity='Compound nucleus formation'
    un = ''
    col(1)='J/P'
    col(2)='population'
    un(2) = 'mb'
    Ncol = 2
    call write_datablock(quantity,Ncol,J2end-J2beg+2,col,un)
    do parity = - 1, 1, 2
      do J2 = J2beg, J2end, 2
        J = J2/2
        write(1, '(5x,f4.1, 1x, a1, 4x, es15.6)') 0.5 * J2, cparity(parity), CNterm(parity, J)
      enddo
    enddo
    quantity='Primary compound nucleus decay'
    un = 'mb'
    col(1)='J/P'
    un(1) = ''
    col(2)='population'
    do type = 0, 6
      col(type+3)=parname(type)
    enddo     
    Ncol = 9
    call write_datablock(quantity,Ncol,J2end-J2beg+2,col,un)
  endif
  do type = -1, 6
    if (adjustTJ(Zcomp, Ncomp, type)) then
      key = 'tjadjust'
      call adjust(Einc, key, Zcomp, Ncomp, type, 0, factor)
    else
      factor = 1.
    endif
    Fnorm(type) = factor / fiso(type)
  enddo
  call densprepare(Zcomp, Ncomp, 1)
!
! *** Output of flux conservation check of transmission coefficients ***
!
  if (flagcheck .and. flagwidth) then
    write(*, '(/" ++++++++++ CHECK OF FLUX CONSERVATION", " OF TRANSMISSION COEFFICIENTS ++++++++++",/)')
    if (wmode == 0) write(*, '(" Hauser-Feshbach model"/)')
    if (wmode == 1) write(*, '(" Moldauer model")')
    if (wmode == 2) write(*, '(" HRTW model"/)')
    if (wmode == 3) write(*, '(" GOE model"/)')
  endif
!
! ************** Initialization of transmission coefficients ***********
!
! The transmission coefficients Tjlnex have values interpolated from the Tjl coefficients of the emission energy grid.
! For the incident channel, the Tjl are exactly calculated.
! Hence, we use exact values for the incident channel in the transmission coefficient array.
!
  do updown = - 1, 1
    do l = 0, lmaxinc
      Tjlnex(k0, Ltarget, updown, l) = Tjlinc(updown, l)
    enddo
  enddo
  nex = maxex(Zcomp, Ncomp)
!
! Initialisation of total astrophysical transmission coefficient
!
  if (flagastro) then
    Tastrotot = 0.
    rhoastrotot = 0.
  endif
!
! ************** Loop over incoming and outgoing channels **************
!
! The variable pardif is used as an indicator of parity conservation for the incident channel.
!
! Loop over compound nucleus parity
!
  do parity = - 1, 1, 2
    pardif = abs(targetP - parity) / 2
!
! compprepare : subroutine to prepare information for initial compound nucleus
! widthprepare: subroutine for preparation of width fluctuation corrections
!
! There are two possible types of calculation for the initial compound nucleus.
! If either width fluctuation corrections or compound nucleus angular distributions are wanted, we need to sum explicitly over all
! possible quantum numbers before we calculate the width fluctuation or angular factor.
! If not, the sum over j,l of the transmission coefficients can be lumped into one factor, which decreases the calculation time.
! In the latter case, the partial decay widths enumhf are calculated in subroutine compprepare.
!
! In order to get do-loops running over integer values, certain quantum numbers are multiplied by 2,
! which can be seen from a 2 present in the corresponding variable names.
! For each loop, the begin and end point is determined from the triangular rule.
!
! For every J and P (parity), first the denominator (total width) denomhf for the Hauser-Feshbach formula is constructed in
! subroutine compprepare.
! Also width fluctuation variables that only depend on J and P and not on the other angular momentum quantum numbers are
! calculated in subroutine widthprepare.
! For the width fluctuation calculation, all transmission coefficients need to be placed in one sequential array.
! Therefore, a counter tnum needs to be followed to keep track of the proper index for the transmission coefficients.
!
! Loop over total angular momentum J (J2) of compound nucleus
!
    do J2 = J2beg, J2end, 2
      popdecay=0.
      partdecay = 0.
      partdecaytot = 0.
      J = J2 / 2
      if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) call tfission(Zcomp, Ncomp, nex, J2, parity)
      call compprepare(Zcomp, Ncomp, J2, parity)
      if (denomhf == 0.) cycle
      if (flagwidth) call widthprepare
      tnumi = 0
!
! Initially, we assume that no width fluctuation corrections or angular distributions are calculated.
! This means the various loops over l and j for the incident and outgoing channel do not need to be performed.
! In this case, the terms needed for the Hauser-Feshbach calculation only depend on J and P, and the corresponding widths
! enumhf and denomhf have already been calculated in subroutine compprepare, where the loops over l and j were done.
!
      jj2beg = 1
      jj2end = 1
      l2beg = 1
      l2end = 1
      jj2primebeg = 1
      jj2primeend = 1
      l2primebeg = 1
      l2primeend = 1
!
! On-set of loop over j in the case of width fluctuations or angular distributions.
!
      if (flagwidth .or. flagcompang .or. flagurr) then
        jj2beg = abs(J2 - targetspin2)
        jj2end = J2 + targetspin2
      endif
!
! 130: Sum over j (jj2) of incident channel
!
      do jj2 = jj2beg, jj2end, 2
!
! On-set of loop over l in the case of width fluctuations or angular distributions.
!
        if (flagwidth .or. flagcompang .or. flagurr) then
          l2beg = abs(jj2 - parspin2i)
          l2end = jj2 + parspin2i
          l2end = min(l2end, 2 * lmaxinc)
        endif
!
! Loop over l (l2) of incident channel
!
        do l2 = l2beg, l2end, 2
          l = l2 / 2
!
! Check parity conservation and make index for transmission coefficient for width fluctuation calculation.
!
! If the parity of the target nucleus is equal (unequal) to the parity of compound nucleus, i.e. pardif=0(1),
! the l-value must be even (odd).
!
          if (flagwidth .or. flagcompang .or. flagurr) then
            if (mod(l, 2) /= pardif) cycle
            updown = (jj2 - l2) / pspin2i
            Tinc = Tjlinc(updown, l)
            tnumi = tnumi + 1
          endif
!
! widthfluc    : subroutine for width fluctuation correction
!
          fluxsum = 0.
!
! 1. Fission channel
!
! The fission contribution is calculated and added to the binary fission cross section.
! Also the transmission coefficient index for width fluctuations is increased.
!
          if (flagfission .and. nfisbar(Zcomp, Ncomp) /= 0) then
            if (flagwidth .or. flagcompang .or. flagurr) then
              tnumo = tnum
              if (flagwidth .and. wmode >= 1) then
                tnumo = tnum - numhill
                factor1 = 0.
                if (tfisA(J, parity, 0) > 0.) then
                  do ihill = 1, numhill
                    tnumo = tnumo + 1
                    ratio = tfisA(J, parity, ihill) / tfisA(J, parity, 0)
                    if (ratio == 0.) cycle
                    call widthfluc(0)
                    factor1 = factor1 + real(Tinc * tfis(J, parity) / denomhf * Wab) * ratio
                  enddo
                  fluxsum = fluxsum + factor1
                endif
              else
                factor1 = real(Tinc * tfis(J, parity) / denomhf)
              endif
            else
              factor1 = real(feed * tfis(J, parity) / denomhf)
            endif
            fisterm = CNfactor * (J2 + 1.) * factor1
            xsbinary( - 1) = xsbinary( - 1) + fisterm
            fisfeedJP(0, 0, maxex(0, 0) + 1, J, parity) = fisfeedJP(0, 0, maxex(0, 0) + 1, J, parity) + fisterm
!
! Extract (L,J) dependent parameters for URR (Gilles Noguere)
!
            if (flagurr .and. l <= lurr) then
              Turrlj( - 1, l, J) = factor1 * denomhf / (Tinc * Wab)
              xsbinarylj( - 1, l, J) = xsbinarylj( - 1, l, J) + CNfactor * (J2 + 1.) * real(factor1)
              nulj( - 1, l, J) = 1
            endif
          endif
!
! 2. Gamma and particle channels
!
          tnumo = tnum + 1
!
! Loop over outgoing channels
!
          do type = 0, 6
            if (type == 1) tnumo = tNinc
            if (parskip(type)) cycle
            parspin2o = int(2. * parspin(type))
            pspin2o = spin2(type)
            Zix = Zindex(Zcomp, Ncomp, type)
            Nix = Nindex(Zcomp, Ncomp, type)
            NL = Nlast(Zix, Nix, 0)
!
! This loop is over all discrete levels and continuum bins of the final nucleus.
!
! Loop over outgoing excitation energies
!
            sumIPE = 0.
            do nexout = 0, maxex(Zix, Nix)
              l2maxhf = 2 * lmaxhf(type, nexout)
              elastic = (type == k0 .and. nexout == Ltarget)
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
                J2res = J2 + parspin2o
                Irspin2beg = mod(J2res, 2)
                Irspin2end = J2res + l2maxhf
                Irspin2end = min(Irspin2end, 2 * maxJ(Zix, Nix, nexout))
              endif
              sumIP = 0.
              if (flagastro) sumIPas = 0.
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
                  if (rho < 1.e-20) cycle
!
! On-set of loop over jprime in the case of width fluctuations or angular distributions.
!
                  sumjl = 0.
                  if (flagwidth .or. flagcompang .or. flagurr) then
                    jj2primebeg = abs(J2 - Irspin2)
                    jj2primeend = J2 + Irspin2
                  endif
!
! Loop over j (jj2) of outgoing channel
!

                  do jj2prime = jj2primebeg, jj2primeend, 2
!
! On-set of loop over lprime in the case of width fluctuations or angular distributions.
!
                    if (flagwidth .or. flagcompang .or. flagurr) then
                      l2primebeg = abs(jj2prime - parspin2o)
                      l2primeend = jj2prime + parspin2o
                      l2primeend = min(l2primeend, l2maxhf)
                      if (type == 0) l2primebeg = max(l2primebeg, 2)
                    endif
!
! Loop over l of outgoing channel
!
                    do l2prime = l2primebeg, l2primeend, 2
                      lprime = l2prime / 2
                      modl = mod(lprime, 2)
!
! We include photons as a special case, with the multipole radiation selection rules (irad=0: M-transition, irad=1: E-transition)
!
                      if (flagwidth .or. flagcompang .or. flagurr) then
!
! 1. Photons
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
! ** Populate the outgoing channels using the compound nucleus formula *
!
! We determine the index for the width fluctuation calculation and call the subroutine that calculates the correction factor.
! The contribution of this particular outgoing channel is added to the sum for the incident channel
! (to check for flux conservation).
!
                        if (type >= 1) tnumo = tnumo + 1
                        if (flagwidth .and. wmode >= 1) then
                          ielas = 0
                          if (elastic .and. jj2 == jj2prime .and. l2 == l2prime) ielas = 1
                          call widthfluc(ielas)
                        endif
                        factor1 = real(Tinc * rho * Tout / denomhf * Wab)
                        fluxsum = fluxsum + factor1
                      else
!
! If NO width fluctuation corrections or angular distributions are required, the partial and total decay widths were
! already determined in subroutine compprepare.
! This means we are now in the short loop with l=j=lprime=jprime=1.
!
                        factor1 = feed * enumhf(type, nexout, Ir, Pprime) / denomhf
                        if (flagastro .and. .not. flagastrogs .and. type == k0) &
 &                        sumIPas = sumIPas + CNfactor * (J2 + 1) * enumhf(type, nexout, Ir, Pprime) **2 / denomhf
                      endif
!
! The contribution is added to the total width for the particular incident l and j.
!
                      compterm = CNfactor * (J2 + 1.) * factor1
                      sumjl = sumjl + compterm
!
! ******************* Compound angular distributions *******************
!
! clebsch     : function for Clebsch-Gordan coefficients
! racah       : function for racah coefficients
!
! Compound angular distributions are calculated for discrete states only.
! It can be easily generalized to the continuum, but we postpone that until we find a reason to do so
! (the deviation from isotropy is assumed to be negligible). Note that we are still inside the most inner loop,
! i.e. as indicated in the manual, all quantum numbers are required for a proper calculation of the compound angular distribution.
! The Legendre coefficients are divide by (2LL+1) to bring them on the same level as the direct reaction Legendre
! coefficients that come out of ECIS.
!
                      if (flagcompang .and. nexout <= NL) then
                        Rspin = 0.5 * Irspin2
                        angfac = (J2 + 1) * (jj2prime + 1) * (jj2 + 1) * (l2 + 1) * &
                          (l2prime + 1) / fourpi
                        LLmax = min(l2, l2prime)
                        LLmax = min(LLmax, J2)
                        rJ = 0.5 * J2
                        rlin = real(l)
                        rjin = 0.5 * jj2
                        rlout = real(lprime)
                        rjout = 0.5 * jj2prime
                        iphase = int(abs(Rspin - parspin(type) - targetspin + parspin(k0)) + 0.1)
                        phase = sgn(iphase)
                        do LL = 0, LLmax, 2
                          rLL = real(LL)
                          clin = clebsch(rlin, rlin, rLL, 0., 0., 0., logfact, numfact)
                          ra1in = racah(rJ, rjin, rJ, rjin, targetspin, rLL, logfact, numfact)
                          ra2in = racah(rjin, rjin, rlin, rlin, rLL, parspin(k0), logfact, numfact)
                          clout = clebsch(rlout, rlout, rLL, 0., 0., 0., logfact, numfact)
                          ra1out = racah(rJ, rjout, rJ, rjout, Rspin, rLL, logfact, numfact)
                          ra2out = racah(rjout, rjout, rlout, rlout, rLL, parspin(type), logfact, numfact)
                          Ablatt = phase * angfac * clin * ra1in * ra2in * clout * ra1out * ra2out
                          cleg(type, nexout, LL) = cleg(type, nexout, LL) + Ablatt * compterm / (2 * LL + 1)
                        enddo
                      endif
!
! ************************* end of all loops ***************************
!
! Increment of the the population arrays of the residual nuclei, the binary reaction cross sections and the contrib array,
! which will be used for the interpolation of the compound emission spectra.
!
                    enddo
                  enddo
                  xspop(Zix, Nix, nexout, Ir, Pprime) = xspop(Zix, Nix, nexout, Ir, Pprime) + sumjl
                  sumIP = sumIP + sumjl
                  if (flagpop) then
                    xspopnucP(Zix, Nix, Pprime) = xspopnucP(Zix, Nix, Pprime) + sumjl
                    xspopexP(Zix, Nix, nexout, Pprime) = xspopexP(Zix, Nix, nexout, Pprime) + sumjl
                    popdecay(type, nexout, Ir, Pprime) = popdecay(type, nexout, Ir, Pprime) + sumjl
                    partdecay(type, Pprime) = partdecay(type, Pprime) + sumjl
                    partdecaytot(type) = partdecaytot(type) + sumjl
                  endif
                enddo
              enddo
              xspopex(Zix, Nix, nexout) = xspopex(Zix, Nix, nexout) + sumIP
              if (nexout > NL) then
                xscompcont(type) = xscompcont(type) + sumIP
                contrib(type, nexout) = contrib(type, nexout) + sumIP
              endif
!
! Compound elastic scattering is excluded from the residual production cross sections.
!
              if ( .not. elastic) sumIPE = sumIPE + sumIP
              if (flagastro .and. .not. flagwidth .and. type == k0 .and. sumIP - sumIPas > xseps)  sumIPE = sumIPE - sumIPas
            enddo
            xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) + sumIPE
            xsbinary(type) = xsbinary(type) + sumIPE
!
! Extract (L,J) dependent parameters for URR (advice of Gilles Noguere)
!
            if (flagurr .and. l <= lurr) then
              Turrlj(type, l, J) = sumIPE * denomhf / (CNfactor * (J2 + 1.) * Tinc * Wab)
              xsbinarylj(type, l, J) = xsbinarylj(type, l, J) + real(sumIPE)
              if (type == 0) then
                nulj(type, l, J) = nulj(type, l, J) + 1
              else
                nulj(type, l, J) = 1
              endif
            endif
          enddo
!
! Determine angular momentum range for URR calculation
!
          if (flagurr .and. l <= lurr) then
            Turrljinc(l, J) = Turrljinc(l, J) + Tinc
            Purrlj(l, J) = parity
            if (l <= lminU) lminU = l
            if (l >= lmaxU) lmaxU = l
            if (J <= JminU(l)) JminU(l) = J
            if (J >= JmaxU(l)) JmaxU(l) = J
          endif
!
! ****** Check of flux conservation of transmission coefficients *******
!
! This check is included to test the stability of the width fluctuation calculation.
!
          if (flagcheck .and. flagwidth) then
            if (fluxsum == 0.) fluxsum = Tinc
            write(*, '(" Parity=", a1, "  J=", f4.1, "  j=", f4.1, "  l=", i2, "  T(j,l)=", es12.5, &
 &            "  Sum over outgoing channels=", es12.5, "  Ratio=", f8.5)') cparity(parity), 0.5*J2, 0.5*jj2, l, &
 &            Tinc, fluxsum, Tinc / fluxsum
          endif
        enddo
      enddo
      if (flagpop) then
        rJ = 0.5 * J2
        write(1,'(5x,f4.1,1x,a1,4x,8es15.6)') rJ,cparity(parity),CNterm(parity,J),(partdecaytot(type),type=0,6)
        if (flagdecay) then
          do type = 0, 6
            if (parskip(type)) cycle
            Zix = Zindex(Zcomp, Ncomp, type)
            Nix = Nindex(Zcomp, Ncomp, type)
            Zres = ZZ(Zcomp, Ncomp, type)
            Nres = NN(Zcomp, Ncomp, type)
            Ares = AA(Zcomp, Ncomp, type)
            oddres = mod(Ares, 2)
            do Pprime = - 1, 1, 2
              write(*,'(/" Total Pprime=",i2,":",es10.3," via ",a8," emission"/)') Pprime,partdecay(type,Pprime),parname(type)
              write(*, '(" bin    Ex", 10("    J=", f4.1)/)') (Ir + 0.5 * oddres, Ir = 0, 9)
              do nexout = 0, maxex(Zix, Nix)
                write(*, '(1x, i3, f8.3, 10es10.3)') nexout, Ex(Zix, Nix, nexout), (popdecay(type, nexout, Ir, Pprime), Ir = 0, 9)
              enddo
            enddo
          enddo
        endif
      endif
    enddo
  enddo
  if (flagpop) then
    quantity='Isospin factors to reduce emission'
    un = ''
    col(1)='particle'
    col(2)='isospin factor'
    Ncol = 2
    call write_datablock(quantity,Ncol,7,col,un)
    do type = 0, 6
      write(1, '(3x, a8, 4x, es15.6)') parname(type), fiso(type)
    enddo
    close(1)
    call write_outfile(popfile,flagoutall)
  endif
!
! ************************** Astrophysics ******************************
!
  if (flagastro .and. .not. flagastrogs) then
    if (Tastrotot >= transeps .and. rhoastrotot >= 0..and.flagwidth) call astrotarget
    if (flagwidth) then
      if (maxex(parZ(k0), parN(k0)) + 1 > 3) ewfc = Einc
    endif
  endif
!
! *********** Output of fission transmission coefficients **************
!
! tfissionout: subroutine for output of fission transmission coefficients
!
  if (flagfisout) call tfissionout(Zcomp, Ncomp, nex)
!
! **** ECIS calculation of compound cross sections (reference only) ****
!
! In addition to a calculation by TALYS, a compound nucleus run by ECIS can be requested.
! The results will however not be used for TALYS but are just for comparison.
!
! raynalcomp  : subroutine for ECIS calculation of compound cross sections (reference only)
!
  if (flageciscomp) call raynalcomp
  return
end subroutine comptarget
! Copyright A.J. Koning 2021
