subroutine densitycum(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : CUmulative level density and spin distribution
!
! Author    : Arjan Koning
!
! 2024-06-18: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
!   dbl             ! double precision kind
! Variables for level density
!   D0              ! s - wave resonance spacing in eV
! Variables for level density
!   alev            ! level density parameter
!   alimit          ! asymptotic level density parameter
!   cfermi          ! width of Fermi distribution for damping
!   ctable          ! constant to adjust tabulated level densities
!   deltaW          ! shell correction in nuclear mass
!   E0              ! particle constant of temperature formula
!   Exmatch         ! matching point for Ex
!   filedensity     ! flag for level densities on separate files
!   flagcol         ! flag for collective enhancement of level density
!   gammald         ! gamma - constant for asymptotic level density parameter
!   ldmodel         ! level density model
!   Ncum            ! number of cumulative levels (integral of level density)
!   rhoexp          ! level density of experimental discrete levels
!   Nlow            ! lowest discrete level for temperature matching
!   Ntop            ! highest discrete level for temperature matching
!   pair            ! pairing energy
!   Pshift          ! adjustable pairing shift
!   ptable          ! constant to adjust tabulated level densities
!   T               ! temperature
!   Ufermi          ! energy of Fermi distribution for damping
! Variables for fission
!   flagfission     ! flag for fission
! Variables for masses
!   beta2           ! deformation parameter
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   NN              ! neutron number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   pardis          ! parity distribution
!   nuc             ! symbol of nucleus
! Variables for resonance parameters
!   D0theo          ! mean s - wave resonance spacing
!   D1theo          ! mean p - wave resonance spacing
!   dD0             ! uncertainty in D0
!   Eavres          ! number of resonances
! Variables for levels
!   edis            ! energy of level
!   nlevmax2        ! maximum number of levels
! Variables for level density
!   delta0          ! systematical pairing energy
!   edens           ! energy grid for tabulated level densities
!   nendens         ! number of energies for level density grid
!   Nlast           ! last discrete level
! Variables for fission parameters
!   nfisbar         ! number of fission barrier parameters
! Variables for masses
!   nucmass         ! mass of nucleus
!   S               ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer           :: A             ! mass number of target nucleus
  integer           :: i             ! counter
  integer           :: Nav
  integer           :: ibeg
  integer           :: iend
  integer           :: k             ! counter
  integer           :: Nk            ! counter
  integer           :: ibar          ! fission barrier
  integer           :: J             ! spin of level
  integer           :: ldmod         ! level density model
  integer           :: N             ! neutron number of residual nucleus
  integer           :: nex           ! excitation energy bin of compound nucleus
  integer           :: Nix           ! neutron number index for residual nucleus
  integer           :: NL            ! last discrete level
  integer           :: NT            ! help variable
  integer           :: odd           ! odd (1) or even (0) nucleus
  integer           :: parity        ! parity
  integer           :: Z             ! charge number of target nucleus
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: dEx           ! excitation energy bin for population arrays
  real(sgl)         :: Eex           ! excitation energy
  real(sgl)         :: denom
  real(sgl)         :: x
  real(sgl)         :: pp
  real(sgl)         :: Dexp          ! 
  real(sgl)         :: dDexp         ! 
  real(sgl)         :: Dtheo         ! 
  real(sgl)         :: Dglob
  real(sgl)         :: dDglob
  real(sgl)         :: Ea
  real(sgl)         :: Eb
  real(sgl)         :: Erange
  real(sgl)         :: P             ! pairing energy
  real(sgl)         :: SS            ! separation energy
  real(dbl)         :: chi2sum       ! help variable
  real(dbl)         :: Frmssum       ! help variable
  real(dbl)         :: Ermssum       ! help variable
  real(dbl)         :: avdevsum      ! help variable
  real(dbl)         :: Ri            ! C/E
  real(dbl)         :: dens          ! total level density
  real(dbl)         :: density       ! level density
  real(dbl)         :: densitytot    ! total level density
  real(dbl)         :: densitytotP   ! total level density per parity
  real(dbl)         :: ldtotP        ! total level density per parity
  real(dbl)         :: Rdist(numdens, -1:1, 0:numJ) ! spin distribution
!
! ********************** Level density parameters **********************
!
! aldmatch   : function to determine effective level density parameter
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  SS = S(Zix, Nix, 1)
  P = pair(Zix, Nix)
  ldmod = ldmodel(Zix, Nix)
!
! Cumulative number of discrete levels vs. integrated level density
!
  NL = Nlow(Zix, Nix, 0)
  NT = Ntop(Zix, Nix, 0)
  if (filedensity) then
    Dtheo=D0theo(Zix, Nix)
    Dexp=D0(Zix, Nix)
    dDexp=dD0(Zix, Nix)
    Dglob=D0global(Zix, Nix)
    dDglob=dD0global(Zix, Nix)
    if (dDexp <= 1.e-30) then
      chi2D0(Zix, Nix) = 0.
    else
      chi2D0(Zix, Nix) = (abs(Dtheo - Dexp) / dDexp) **2
    endif
!
! Cumulative probability derived as follows:
!
!          cdf=0.5*(1.+erf(x))
!          pp=2.*(cdf-0.5)
!
    FrmsD0(Zix, Nix) = 0.
    ErmsD0(Zix, Nix) = 0.
    CED0(Zix, Nix) = 0.
    if (Dexp > 1.e-30) then
      CED0(Zix, Nix) = Dtheo / Dexp
      if (Dtheo > 0. .and. Dexp > 0.) then
        if (dDexp > 0.) then
          x = (Dtheo - Dexp)/(dDexp*sqrt(2.))
          if (Dtheo > Dexp) then
            pp = erf(x)
          else
            pp = -erf(x)
          endif
          Ri = 1. + (CED0(Zix, Nix)-1.)*pp
        else
          Ri = CED0(Zix, Nix)
        endif
        if (Ri >= 1.) then
          FrmsD0(Zix, Nix) = Ri
        else
          FrmsD0(Zix, Nix) = 1./Ri
        endif
        ErmsD0(Zix, Nix) = Ri
      endif
    endif
    CGD0(Zix, Nix) = 0.
    if (Dglob > 1.e-30) CGD0(Zix, Nix) = Dtheo / Dglob
    chi2lev(Zix, Nix) = 0.
    Frmslev(Zix, Nix) = 0.
    Ermslev(Zix, Nix) = 0.
    avdevlev(Zix, Nix) = 0.
    k = 0
    chi2sum = 0.
    Frmssum = 0.
    Ermssum = 0.
    avdevsum = 0.
    do i = 0, nlevmax2(Zix, Nix)
      Ncum(Zix, Nix, i) = 0.
      rhoexp(Zix, Nix, i) = 0.
    enddo
    do i = 1, nlevmax2(Zix, Nix)
      if (edis(Zix, Nix, i) == 0.) cycle
      Eex = 0.5 * (edis(Zix, Nix, i) + edis(Zix, Nix, i - 1))
      dEx = edis(Zix, Nix, i) - edis(Zix, Nix, i - 1)
      dens = densitytot(Zix, Nix, Eex, 0, ldmod)
      Ncum(Zix, Nix, i) = Ncum(Zix, Nix, i-1) + dens * dEx
      if (i == NL) Ncum(Zix, Nix, i) = real(NL)
      Nav = 10
      ibeg = max(i - Nav/2,0)
      iend = min(i + Nav/2,nlevmax2(Zix, Nix))
      Nav = iend - ibeg + 1
      Ea = edis(Zix, Nix, ibeg)
      Eb = edis(Zix, Nix, iend)
      Erange = Eb - Ea
      if (Erange > 0.) rhoexp(Zix, Nix, i) = Nav / Erange
      if (i >= NL .and. i <= NT) then
        chi2sum = chi2sum + ((Ncum(Zix, Nix, i) - i) **2 ) / sqrt(real(i)) / max(NT-NL,1)
        Ri = Ncum(Zix, Nix, i) / i
        Frmssum = Frmssum + log(Ri)**2 
        Ermssum = Ermssum + log(Ri) 
        avdevsum = avdevsum + abs(Ncum(Zix, Nix, i) - i) 
      endif
    enddo
    denom = real(NT - NL)
    chi2lev(Zix, Nix) = chi2sum 
    if (denom > 0.) then
      Frmslev(Zix, Nix) = exp(sqrt(Frmssum / denom))
      Ermslev(Zix, Nix) = exp(Ermssum / denom)
      avdevlev(Zix, Nix) = avdevsum / denom
    endif
    Nk = k
  endif
!
! Level densities per parity on separate files, as in tabulated format
!
  if (filedensity) then
!
! Spin distribution
!
    odd = mod(A, 2)   
    do ibar = 0, nfisbar(Zix, Nix)
      do parity = -1, 1, 2
        do nex = 1, nendens(Zix, Nix)
          Eex = edens(nex)
          ldtotP = densitytotP(Zix, Nix, Eex, parity, ibar, ldmod)
          do J = 0, numJ
            Rdist(nex, parity, J) = density(Zix, Nix, Eex, real(J + 0.5 * odd), parity, ibar, ldmod) / ldtotP
          enddo
        enddo
      enddo
    enddo
  endif
  return
end subroutine densitycum
! Copyright A.J. Koning 2024
