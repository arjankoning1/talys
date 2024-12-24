subroutine thermal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Estimate of thermal cross sections
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
!   sgl              ! single precision kind
! All global variables
!   numisom          ! number of isomers
!   numlev           ! maximum number of discrete levels
!   numNastro        ! maximal number of neutrons away from initial CN for astroph. calcs
!   numZastro        ! maximal number of protons away from initial CN for astroph. calcs
! Variables for compound reactions
!   xsalphatherm     ! thermal (n, a) cross section
!   xscaptherm       ! thermal capture cross section
!   xsptherm         ! thermal (n, p) cross section
! Variables for energies
!   Ninclow        ! number of incident energies below Elow
! Variables for levels
!   Lisomer         ! level number of isomer
! Variables for numerics
!   maxN             ! maximal number of neutrons away from initial compound nucleus
!   maxZ             ! maximal number of protons away from initial compound nucleus
!   xseps            ! limit for cross sections
! Variables for basic reaction
!   flagastro        ! flag for calculation of astrophysics reaction rate
!   flagchannels     ! flag for exclusive channels calculation
! Variables for input energies
!   eninc            ! incident energy in MeV
! Variables for main input
!   k0               ! index of incident particle
! Variables for total cross sections
!   xsexclcont       ! exclusive single channel cross section for continuum
!   xsexclusive      ! exclusive single channel cross section
! Variables for energy grid
!   E1v              ! energy at end of 1 / v region
! Variables for energies
!   Ethresh          ! threshold incident energy for residual nucleus
!   Ethrexcl         ! threshold incident energy for exclusive channel
!   idchannel        ! identifier for exclusive channel
! Variables for exclusive channels
!   exclbranch       ! exclusive channel yield per isomer
!   idnum            ! counter for exclusive channel
!   xschannel        ! channel cross section
!   xsgamchannel     ! gamma channel cross section
!   xsgamdischan     ! discrete gamma channel cross section
!   xsratio          ! ratio of exclusive cross section over residual p
! Variables for multiple emission
!   xsngn            ! total (projectile, gamma - ejectile) cross section
! Variables for binary reactions
!   xscompdisc       ! compound cross section for discrete state
!   xscompel         ! compound elastic cross section
!   xsdisc           ! total cross section for discrete state
!   xsdisctot        ! total cross section summed over discrete states
!   xselastot        ! total elastic cross section (shape + compound)
! Variables for incident channel
!   xsbinary         ! cross section from initial compound to residual nucleus
!   xsbranch         ! branching ratio for isomeric cross section
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsdirdiscsum     ! total direct cross section
!   xselasinc        ! total elastic cross section (neutrons only) for inc. channel
!   xspopex          ! population cross section summed over spin and parity
!   xspopnuc         ! population cross section per nucleus
!   xspreeqsum       ! total preequilibrium cross section summed over particles
! Variables for direct capture
!   xsracape         ! direct capture cross section
! Constants
!   parN             ! neutron number of particle
!   parZ             ! charge number of particle
! Variables for thermal cross sections
!   fexclbranch      ! exclusive channel yield per isomer
!   fxsbinary        ! cross section from initial compound to residual n
!   fxsbranch        ! branching ratio for isomeric cross section
!   fxschaniso       ! channel cross section per isomer
!   fxschannel       ! channel cross section
!   fxscompdisc      ! compound cross section for discrete state
!   fxscompel        ! compound elastic cross section
!   fxscompnonel     ! total compound non - elastic cross section
!   fxsdirdisc       ! direct cross section for discrete state
!   fxsdirdiscsum    ! total direct cross section
!   fxsdisc          ! total cross section for discrete state
!   fxsdisctot       ! total cross section summed over discrete states
!   fxselasinc       ! total elastic cross section (neutrons only) for i
!   fxselastot       ! total elastic cross section (neutrons only) for i
!   fxsexclcont      ! exclusive single channel cross section for contin
!   fxsexclusive     ! exclusive single channel cross section
!   fxsgamchannel    ! gamma channel cross section
!   fxsgamdischan    ! discrete gamma channel cross section
!   fxsngn           ! total (projectile, gamma - ejectile) cross section
!   fxsnonel         ! non - elastic cross section for incident channel
!   fxspopex         ! population cross section summed over spin and par
!   fxspopnuc        ! population cross section per nucleus
!   fxspreeqsum      ! total preequilibrium cross section summed over pa
!   fxsracape        ! direct capture cross section
!   fxsratio         ! ratio of exclusive cross section over residual pr
!   fxsreacinc       ! reaction cross section for incident channel
!   fxstotinc        ! total cross section (neutrons only) for incident
! Variables for level density
!   Nlast            ! last discrete level
! Variables for astro
!   xsastro          ! cross section for astrophysical calculatio
!
! *** Declaration of local data
!
  implicit none
  integer   :: i             ! counter
  integer   :: i1            ! value
  integer   :: i2            ! value
  integer   :: idc           ! help variable
  integer   :: Ncomp         ! neutron number index for compound nucleus
  integer   :: Nix           ! neutron number index for residual nucleus
  integer   :: nen           ! energy counter
  integer   :: Nis           ! number of isotope
  integer   :: nex           ! excitation energy bin of compound nucleus
  integer   :: type          ! particle type
  integer   :: Zix           ! charge number index for residual nucleus
  integer   :: Zcomp         ! proton number index for compound nucleus
  real(sgl) :: branchres     ! branching ratio in resonance range
  real(sgl) :: ctherm        ! constant for 1/v capture cross section function
  real(sgl) :: ealog         ! help variable
  real(sgl) :: elog          ! help variable
  real(sgl) :: Eratio        ! energy ratio
  real(sgl) :: Ereslog       ! logarithm of energy at start of resonance region
  real(sgl) :: Etherm        ! thermal energy
  real(sgl) :: R             ! radius
  real(sgl) :: ratio         ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: ratioalpha    ! ratio thermal/first energy for proton
  real(sgl) :: ratiop        ! ratio thermal/first energy for proton
  real(sgl) :: ratiores      ! ratio start of resonance region/first energy
  real(sgl) :: ratioresalpha ! ratio start of resonance region/first energy
  real(sgl) :: ratioresp     ! ratio start of resonance region/first energy
  real(sgl) :: Rres          ! ratio start of resonance region/first energy
  real(sgl) :: xs            ! help variable
  real(sgl) :: xsa           ! help variable
  real(sgl) :: xsalog        ! help variable
  real(sgl) :: xsalpha1      ! (n,a) cross section
  real(sgl) :: xsalphares    ! (n,a) cross section in resonance region
  real(sgl) :: xscap1        ! capture cross section at first incident energy
  real(sgl) :: xscapres      ! capture cross section in resonance region
  real(sgl) :: xsp1          ! (n,p) cross section
  real(sgl) :: xspres        ! (n,p) cross section in resonance region
  real(sgl) :: xsres         ! cross section at start of resonance region
  real(sgl) :: xsreslog      ! cross section at start of resonance region
!
! *********************** Extrapolate cross sections *******************
!
! For non-threshold channels, the cross sections are extrapolated down to 1.e-5 eV.
! Capture values at thermal energies are used.
! For energies up to 1 eV, the 1/sqrt(E) law is used.
! Between 1 eV and the first energy at which TALYS performs the statistical model calculation, we use logarithmic interpolation.
!
! ctherm       : constant for 1/v capture cross section function
!
  Zcomp = 0
  Ncomp = 0
  Etherm = 2.53e-8
  ctherm = sqrt(Etherm) * xscaptherm(-1)
  xscapres = ctherm / sqrt(E1v)
  xscap1 = xspopnuc(Zcomp, Ncomp)
  if (xscap1 > 0.) then
    ratio = xscaptherm(-1) / xscap1
    ratiores = xscapres / xscap1
  else
    ratio = 1.
    ratiores = 1.
  endif
!
! Protons
!
  ratiop = ratio
  ratioresp = ratiores
  if (xsptherm(-1) /= 0.) then
    ctherm = sqrt(Etherm) * xsptherm(-1)
    xspres = ctherm / sqrt(E1v)
    xsp1 = xspopnuc(1, 0)
    if (xsp1 > 0.) then
      ratiop = xsptherm(-1) / xsp1
      ratioresp = xspres / xsp1
    else
      ratiop = 1.
    endif
  endif
!
! Alpha particles
!
  ratioalpha = ratio
  ratioresalpha = ratiores
  if (xsalphatherm(-1) /= 0.) then
    ctherm = sqrt(Etherm) * xsalphatherm(-1)
    xsalphares = ctherm / sqrt(E1v)
    xsalpha1 = xspopnuc(2, 2)
    if (xsalpha1 > 0.) then
      ratioalpha = xsalphatherm(-1) / xsalpha1
      ratioresalpha = xsalphares / xsalpha1
    else
      ratioalpha = 1.
    endif
  endif
!
! Determine cross sections on low-energy grid
!
  Ereslog = log(E1v)
  do nen = 1, Ninclow
    elog = log(eninc(nen))
    ealog = log(eninc(Ninclow + 1))
    Eratio = sqrt(Etherm) / sqrt(eninc(nen))
    fxsnonel(nen) = 0.
    fxselastot(nen) = 0.
    fxstotinc(nen) = 0.
    fxscompel(nen) = 0.
    fxselasinc(nen) = 0.
    fxsreacinc(nen) = 0.
    fxscompnonel(nen) = 0.
    fxsdirdiscsum(nen) = 0.
    fxspreeqsum(nen) = 0.
    fxsracape(nen) = 0.
!
! Exclusive channel cross sections
!
! pol1         : subroutine for interpolation of first order
!
    if (flagchannels) then
      do idc = 0, idnum
        fxschannel(nen, idc) = 0.
        fxsgamchannel(nen, idc) = 0.
        do i1 = 1, numlev
          do i2 = 0, i1
            fxsgamdischan(nen, idc, i1, i2) = 0.
          enddo
        enddo
        if (eninc(nen) <= Ethrexcl(idc, 0)) cycle
        if (idchannel(idc) == 100000) cycle
        xsa = xschannel(idc)
        if (xsa <= xseps) cycle
        R = ratio
        Rres = ratiores
        branchres = 0.
        Zix = 0
        Nix = 0
        if (idchannel(idc) == 010000) then
          Zix = 1
          Nix = 0
          R = ratiop
          Rres = ratioresp
        endif
        if (idchannel(idc) == 000001) then
          Zix = 2
          Nix = 2
          R = ratioalpha
          Rres = ratioresalpha
        endif
        if (eninc(nen) > E1v) then
          xsres = xsa * Rres
          if (xsres <= 0..or.xsa <= 0.) cycle
          xsalog = log(xsa)
          xsreslog = log(xsres)
          call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
          fxschannel(nen, idc) = exp(xs)
        else
          xs = xsa * R * Eratio
          fxschannel(nen, idc) = xs
        endif
        fxsnonel(nen) = fxsnonel(nen) + fxschannel(nen, idc)
        fxsratio(nen, idc) = xsratio(idc)
        do nex = 0, Nlast(0, 0, 0)
          fexclbranch(nen, idc, nex) = exclbranch(idc, nex)
          Nis = -1
          do i = 0, numisom
            if (Lisomer(Zix, Nix, i) == nex) then
              Nis = i
              exit
            endif
          enddo
          if (Nis == -1) cycle
          branchres = 0.
          if (idchannel(idc) == 0) then
            if (xscaptherm(-1) > 0.) branchres = xscaptherm(Nis) / xscaptherm(-1)
          endif
          if (idchannel(idc) == 010000) then
            if (xsptherm(-1) > 0.) branchres = xsptherm(Nis) / xsptherm(-1)
          endif
          if (idchannel(idc) == 000001) then
            if (xsalphatherm(-1) > 0.) branchres = xsalphatherm(Nis) / xsalphatherm(-1)
          endif
          if (eninc(nen) <= E1v .and. branchres > 0.) fexclbranch(nen, idc, nex) = branchres
          fxschaniso(nen, idc, nex) = 0.
          if (eninc(nen) <= Ethrexcl(idc, nex)) cycle
          if (eninc(nen) > E1v) then
            xsres = xsa * Rres
            if (xsres <= 0. .or. xsa <= 0.) cycle
            xsalog = log(xsa)
            xsreslog = log(xsres)
            call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
            fxschaniso(nen, idc, nex) = exp(xs) * fexclbranch(nen, idc, nex)
          else
            xs = xsa * R * Eratio
            fxschaniso(nen, idc, nex) = xs * fexclbranch(nen, idc, nex)
          endif
        enddo
        xsa = xsgamchannel(idc)
        if (xsa < xseps) cycle
        if (eninc(nen) > E1v) then
          xsres = xsa * Rres
          if (xsres <= 0. .or. xsa <= 0.) cycle
          xsalog = log(xsa)
          xsreslog = log(xsres)
          call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
          fxsgamchannel(nen, idc) = exp(xs)
        else
          xs = xsa * R * Eratio
          fxsgamchannel(nen, idc) = xs
        endif
        do i1 = 1, numlev
          do i2 = 0, i1
            xsa = xsgamdischan(idc, i1, i2)
            if (xsa < xseps) cycle
            if (eninc(nen) > E1v) then
              xsres = xsa * Rres
              if (xsres <= 0..or.xsa <= 0.) cycle
              xsalog = log(xsa)
              xsreslog = log(xsres)
              call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
              fxsgamdischan(nen, idc, i1, i2) = exp(xs)
            else
              xs = xsa * R * Eratio
              fxsgamdischan(nen, idc, i1, i2) = xs
            endif
          enddo
        enddo
      enddo
    endif
!
! Binary cross sections
!
    do type = 0, 6
      fxsbinary(nen, type) = 0.
      if (eninc(nen) <= Ethresh(parZ(type), parN(type), 0)) cycle
      if (type == k0) cycle
      xsa = xsbinary(type)
      if (xsa < xseps) cycle
      R = ratio
      Rres = ratiores
      if (type == 2) then
        R = ratiop
        Rres = ratioresp
      endif
      if (type == 6) then
        R = ratioalpha
        Rres = ratioresalpha
      endif
      if (eninc(nen) > E1v) then
        xsres = xsa * Rres
        if (xsres <= 0..or.xsa <= 0.) cycle
        xsalog = log(xsa)
        xsreslog = log(xsres)
        call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
        fxsbinary(nen, type) = exp(xs)
      else
        xs = xsa * R * Eratio
        fxsbinary(nen, type) = xs
      endif
    enddo
!
! Residual production cross sections
!
    do Zcomp = 0, maxZ
      do Ncomp = 0, maxN
        fxspopnuc(nen, Zcomp, Ncomp) = 0.
        if (eninc(nen) <= Ethresh(Zcomp, Ncomp, 0)) cycle
        if (Zcomp == parZ(k0) .and. Ncomp == parN(k0)) cycle
        xsa = xspopnuc(Zcomp, Ncomp)
        if (xsa < xseps) cycle
        R = ratio
        Rres = ratiores
        if (Zcomp == 1 .and. Ncomp == 0) then
          R = ratiop
          Rres = ratioresp
        endif
        if (Zcomp == 2 .and. Ncomp == 2) then
          R = ratioalpha
          Rres = ratioresalpha
        endif
        if (eninc(nen) > E1v) then
          xsres = xsa * Rres
          if (xsres <= 0..or.xsa <= 0.) cycle
          xsalog = log(xsa)
          xsreslog = log(xsres)
          call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
          fxspopnuc(nen, Zcomp, Ncomp) = exp(xs)
        else
          xs = xsa * R * Eratio
          fxspopnuc(nen, Zcomp, Ncomp) = xs
        endif
        if ( .not. flagchannels) fxsnonel(nen) = fxsnonel(nen) + fxspopnuc(nen, Zcomp, Ncomp)
        do nex = 0, Nlast(Zcomp, Ncomp, 0)
          fxspopex(nen, Zcomp, Ncomp, nex) = 0.
          if (eninc(nen) <= Ethresh(Zcomp, Ncomp, nex)) cycle
          xsa = xspopex(Zcomp, Ncomp, nex)
          if (xsa < xseps) cycle
          if (eninc(nen) > E1v) then
            xsres = xsa * Rres
            if (xsres <= 0..or.xsa <= 0.) cycle
            xsalog = log(xsa)
            xsreslog = log(xsres)
            call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
            fxspopex(nen, Zcomp, Ncomp, nex) = exp(xs)
          else
            xs = xsa * R * Eratio
            fxspopex(nen, Zcomp, Ncomp, nex) = xs
          endif
          fxsbranch(nen, Zcomp, Ncomp, nex) = xsbranch(Zcomp, Ncomp, nex)
        enddo
        if (flagastro .and. Zcomp <= numZastro .and. Ncomp <= numNastro) xsastro(Zcomp, Ncomp, nen) = fxspopnuc(nen, Zcomp, Ncomp)
      enddo
    enddo
!
! Reactions to discrete states
!
    do type = 0, 6
      fxsexclusive(nen, type) = 0.
      fxsdisctot(nen, type) = 0.
      fxsexclcont(nen, type) = 0.
      fxsngn(nen, type) = 0.
      if (eninc(nen) <= Ethresh(parZ(type), parN(type), 0)) cycle
      if (type == k0) cycle
      R = ratio
      Rres = ratiores
      if (type == 2) then
        R = ratiop
        Rres = ratioresp
      endif
      if (type == 6) then
        R = ratioalpha
        Rres = ratioresalpha
      endif
      xsa = xsexclusive(type)
      if (xsa >= xseps) then
        if (eninc(nen) > E1v) then
          xsres = xsa * Rres
          if (xsres > 0..or.xsa > 0.) then
            xsalog = log(xsa)
            xsreslog = log(xsres)
            call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
            fxsexclusive(nen, type) = exp(xs)
          endif
        else
          xs = xsa * R * Eratio
          fxsexclusive(nen, type) = xs
        endif
      endif
      xsa = xsdisctot(type)
      if (xsa >= xseps) then
        if (eninc(nen) > E1v) then
          xsres = xsa * Rres
          if (xsres > 0..or.xsa > 0.) then
            xsalog = log(xsa)
            xsreslog = log(xsres)
            call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
            fxsdisctot(nen, type) = exp(xs)
          endif
        else
          xs = xsa * R * Eratio
          fxsdisctot(nen, type) = xs
        endif
      endif
      xsa = xsexclcont(type)
      if (xsa >= xseps) then
        if (eninc(nen) > E1v) then
          xsres = xsa * Rres
          if (xsres > 0..or.xsa > 0.) then
            xsalog = log(xsa)
            xsreslog = log(xsres)
            call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
            fxsexclcont(nen, type) = exp(xs)
          endif
        else
          xs = xsa * R * Eratio
          fxsexclcont(nen, type) = xs
        endif
      endif
      xsa = xsngn(type)
      if (xsa >= xseps) then
        if (eninc(nen) > E1v) then
          xsres = xsa * Rres
          if (xsres > 0..or.xsa > 0.) then
            xsalog = log(xsa)
            xsreslog = log(xsres)
            call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
            fxsngn(nen, type) = exp(xs)
          endif
        else
          xs = xsa * R * Eratio
          fxsngn(nen, type) = xs
        endif
      endif
      do nex = 0, Nlast(parZ(type), parN(type), 0)
        fxsdisc(nen, type, nex) = 0.
        fxsdirdisc(nen, type, nex) = 0.
        fxscompdisc(nen, type, nex) = 0.
        xsa = xsdisc(type, nex)
        if (xsa >= xseps) then
          if (eninc(nen) > E1v) then
            xsres = xsa * Rres
            if (xsres > 0..or.xsa > 0.) then
              xsalog = log(xsa)
              xsreslog = log(xsres)
              call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
              fxsdisc(nen, type, nex) = exp(xs)
            endif
          else
            xs = xsa * R * Eratio
            fxsdisc(nen, type, nex) = xs
          endif
        endif
        xsa = xsdirdisc(type, nex)
        if (xsa >= xseps) then
          if (eninc(nen) > E1v) then
            xsres = xsa * Rres
            if (xsres > 0..or.xsa > 0.) then
              xsalog = log(xsa)
              xsreslog = log(xsres)
              call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
              fxsdirdisc(nen, type, nex) = exp(xs)
            endif
          else
            xs = xsa * R * Eratio
            fxsdirdisc(nen, type, nex) = xs
          endif
        endif
        xsa = xscompdisc(type, nex)
        if (xsa >= xseps) then
          if (eninc(nen) > E1v) then
            xsres = xsa * Rres
            if (xsres > 0..or.xsa > 0.) then
              xsalog = log(xsa)
              xsreslog = log(xsres)
              call pol1(Ereslog, ealog, xsreslog, xsalog, elog, xs)
              fxscompdisc(nen, type, nex) = exp(xs)
            endif
          else
            xs = xsa * R * Eratio
            fxscompdisc(nen, type, nex) = xs
          endif
        endif
      enddo
    enddo
!
! Total cross sections
!
    fxselastot(nen) = xselastot
    fxstotinc(nen) = fxselastot(nen) + fxsnonel(nen)
    fxscompel(nen) = xscompel
    fxselasinc(nen) = xselasinc
    fxsreacinc(nen) = fxsnonel(nen) + fxscompel(nen)
    fxscompnonel(nen) = fxsnonel(nen)
    fxsdirdiscsum(nen) = xsdirdiscsum
    fxspreeqsum(nen) = xspreeqsum
    fxsracape(nen) = xsracape
  enddo
  return
end subroutine thermal
! Copyright A.J. Koning 2021
