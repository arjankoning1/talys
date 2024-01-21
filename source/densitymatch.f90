subroutine densitymatch(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Level density matching solution
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
!   sgl              ! single precision kind
!   dbl              ! double precision kind
! All global variables
!   nummatchT        ! maximum number of energy points for T matching
! Variables for level density
!   deltaW           ! shell correction in nuclear mass
!   E0               ! particle constant of temperature formula
!   E0adjust         ! adjustable factor for E0
!   Exmatch          ! matching point for Ex
!   Exmatchadjust    ! adjustable factor for matching energy
!   flagcol          ! flag for collective enhancement of level density
!   flagctmglob      ! flag for global CTM model (no discrete level info)
!   gammald          ! gamma - constant for asymptotic level density para
!   ldmodel          ! level density model
!   Nlow             ! lowest discrete level for temperature matching
!   Ntop             ! highest discrete level for temperature matching
!   pair             ! pairing energy
!   T                ! temperature
!   Tadjust          ! adjustable factor for temperature
! Variables for nuclides
!   AA               ! mass number of residual nucleus
! Variables for resonance parameters
!   D0theo           ! mean s - wave resonance spacing
!   D1theo           ! mean p - wave resonance spacing
!   Dl               ! mean resonance spacing per l value
! Variables for fission parameters
!   efistrrot        ! energy of rotational transition states
!   nfisbar          ! number of fission barrier parameters
! Variables for levels
!   edis             ! energy of level
! Variables for level density
!   E0save           ! E0 value saved for matching routine
!   EL               ! lower matching level energy
!   EP               ! higher matching level energy
!   Exmemp           ! empirical estimate for matching point for Ex
!   ldexist          ! flag for existence of level density
!   ldparexist       ! flag for existence of tabulated level density
!   logrho           ! logarithm of level density
!   NLo              ! lowest discrete level for temperature matching
!   NP               ! highest discrete level for temperature matching
!   Tmemp            ! empirical estimate for temperature
!   temprho          ! temperature
!
! *** Declaration of local data
!
  implicit none
  character(len=7)  :: key     ! keyword
  character(len=132):: denfile ! level density file
  integer   :: A               ! mass number of target nucleus
  integer   :: Z               ! charge number of target nucleus
  integer   :: ia              ! mass number of target nucleus
  integer   :: iz              ! charge number of target nucleus
  integer   :: i               ! counter
  integer   :: istat           ! file logical
  integer   :: ibar            ! fission barrier
  integer   :: j               ! counter
  integer   :: nEx             ! number of energy points for T matching
  integer   :: Nix             ! neutron number index for residual nucleus
  integer   :: Nstart          ! energy starting point from which T starts to behave  smoothly (for interpolation)
  integer   :: Zix             ! charge number index for residual nucleus
  real(sgl) :: ald             ! level density parameter
  real(sgl) :: dEx             ! excitation energy bin for population arrays
  real(sgl) :: E0m             ! constant of temperature formula
  real(sgl) :: Eex             ! excitation energy
  real(sgl) :: Exend           ! end of possible energy region
  real(sgl) :: Exm             ! maximal attainable energy
  real(sgl) :: ignatyuk        ! function for energy dependent level density parameter a
  real(sgl) :: Kcoll           ! total collective enhancement
  real(sgl) :: Krot            ! rotational enhancement factor
  real(sgl) :: Kvib            ! vibrational enhancement factor
  real(sgl) :: logrholoc(-1:1) ! logarithm of level density
  real(sgl) :: logrhomatch     ! logarithm of level density at matching energy
  real(sgl) :: P               ! pairing energy
  real(sgl) :: Tm              ! nuclear temperature
  real(sgl) :: U               ! excitation energy minus pairing energy
  real(dbl) :: fermi           ! function for Fermi gas level density formula
  real(dbl) :: rhomatch        ! level density at matching point
  real(sgl) :: val             ! help variable
!
! ************* Determine level density matching parameters ************
!
! For non-fissile nuclei, there are no fission barriers and only the loop for ibar=0 is performed,
! i.e. for level densities relative to the ground state.
!
  do ibar = 0, nfisbar(Zix, Nix)
!
! 1. Determine matching levels and energies.
!
    NP = Ntop(Zix, Nix, ibar)
    NLo = Nlow(Zix, Nix, ibar)
!
! A. Parameters for the ground state level density
!
    if (ibar == 0) then
      EL = edis(Zix, Nix, NLo)
      EP = edis(Zix, Nix, NP)
    else
!
! B. Parameters for the fission level density
!
      EL = efistrrot(Zix, Nix, ibar, NLo)
      EP = efistrrot(Zix, Nix, ibar, NP)
    endif
    if (ldmodel(Zix, Nix) == 2 .or. ldmodel(Zix, Nix) == 3 .or. ldexist(Zix, Nix, ibar)) cycle
!
! 2. Solve level density matching problem for Gilbert-Cameron model
!
! ignatyuk  : function for energy dependent level density parameter a
! colenhance: subroutine for collective enhancement
! fermi     : function for Fermi gas level density formula smoothly (for interpolation)
!
! Calculate logarithm of level density and its derivative (temperature) on excitation energy grid.
!
    P = pair(Zix, Nix)
    A = AA(Zix, Nix, 0)
    Exend = 20. + 300. / A
    dEx = 0.1
    nEx = int(Exend / dEx)
    do i = 1, nummatchT
      logrho(i) = 0.
      temprho(i) = 0.
    enddo
    do i = nEx, 1, - 1
      do j = - 1, 1
        Eex = dEx * (i + 0.5 * j)
        U = Eex - P
        if (U > 0.) then
          ald = ignatyuk(Zix, Nix, Eex, ibar)
          call colenhance(Zix, Nix, Eex, ald, ibar, Krot, Kvib, Kcoll)
          logrholoc(j) = real(log(Kcoll * fermi(Zix, Nix, ald, Eex, P, ibar)))
        else
          logrholoc(j) = 0.
        endif
      enddo
      logrho(i) = logrholoc(0)
      if (logrholoc(1) /= logrholoc( - 1)) temprho(i) = dEx / (logrholoc(1) - logrholoc( - 1))
      if (temprho(i) <= 0.1) temprho(i) = temprho(i + 1)
    enddo
    Nstart = 1
    do i = nEx, 1, - 1
      if (i < nEx .and. temprho(i) >= temprho(i + 1)) then
        Nstart = i + 1
        exit
      endif
    enddo
!
! matching   : subroutine to determine matching between temperature and Fermi-gas region
! pol1       : subroutine for interpolation of first order
!
! Light nuclides
!
    if (A <= 18) then
      Z = ZZ(Zix, Nix, 0)
      denfile = trim(path) // 'density/ground/ctm/ctm.light'
      open(unit=1, file=denfile, status='unknown')
      do
        read(1, *, iostat = istat) key, iz, ia, val
        if (istat == -1) exit
        if (Z == iz .and. ia == A) then
          if (trim(key) == 'T') T(Zix, Nix, 0) = val
          if (trim(key) == 'E0') E0(Zix, Nix, 0) = val
          if (trim(key) == 'Exmatch') Exmatch(Zix, Nix, 0) = val
        endif
      enddo
      close(1)
    endif
    do
      Tm = T(Zix, Nix, ibar)
      Exm = Exmatchadjust(Zix, Nix, ibar) * Exmatch(Zix, Nix, ibar)
      E0m = E0(Zix, Nix, ibar)
      E0save = E0m
!
! Empirical estimates needed in case of trouble or no discrete levels.
! In those cases, the empirical formula for T is the starting point.
!
      if (flagcol(Zix, Nix)) then
        Tmemp = - 0.22 + 9.4 / sqrt(max(A * (1. + gammald(Zix, Nix) * deltaW(Zix, Nix, 0)), 1.))
        Exmemp = 2.67 + 253. / A + P
      else
        Tmemp = - 0.25 + 10.2 / sqrt(max(A * (1. + gammald(Zix, Nix) * deltaW(Zix, Nix, 0)), 1.))
        Exmemp = 2.33 + 253. / A + P
      endif
      Tmemp = max(Tmemp, 0.1)
      Exmemp = max(Exmemp, 0.1)
!
! Normal case: CTM parameters are derived from matching automatically, i.e. T an Exmatch are not given in the input file.
!
      if (Tm == 0..and.Exm == 0.) then
!
! A. Discrete levels given
!
        if (ldparexist(Zix, Nix) .and. .not. flagctmglob) then
          call matching(Zix, Nix, Exm, ibar)
!
! If Exm was set to 0 in matching.f, it means that no reasonable value was found.
! In that case we first use an empirical value for T.
!
          if (Exm > 0.) then
            i = int(Exm / dEx)
            call pol1(i * dEx, (i + 1) * dEx, temprho(i), temprho(i + 1), Exm, Tm)
          else
            Tm = Tmemp
            call locate(temprho, Nstart, nEx - 1, Tm, i)
            if (i > 0 .and. i <= nEx - 1) call pol1(temprho(i), temprho(i + 1), i * dEx, (i + 1) * dEx, Tm, Exm)
          endif
        else
!
! B. No discrete levels given and/or global CTM model
!
          Tm = Tmemp
          call locate(temprho, Nstart, nEx - 1, Tm, i)
          if (i > 0) call pol1(temprho(i), temprho(i + 1), i * dEx, (i + 1) * dEx, Tm, Exm)
        endif
!
! If Exm is still unrealistic, we first use an empirical value for T.
!
        ald = ignatyuk(Zix, Nix, Exm, ibar)
        if (Exm <= max(2.25 / ald + P, 0.) + 0.11) Exm = 0.
        if (Exm > 3. * Exmemp) Exm = 0.
      endif
!
! Special case: either T or Exmatch is given in the input.
! If the Exmatch given in the input is unphysical, we choose an empirical value.
!
      if (Exm == 0.) then
        if (Tm == 0.) Tm = Tmemp
        call locate(temprho, Nstart, nEx - 1, Tm, i)
        if (i > 0) call pol1(temprho(i), temprho(i + 1), i * dEx, (i + 1) * dEx, Tm, Exm)
        ald = ignatyuk(Zix, Nix, Exm, ibar)
        if (Exm <= max(2.25 / ald + P, 0.) + 0.11) Exm = Exmemp
        if (Exm == 0.) Exm = Exmemp
        if (Exm > 3. * Exmemp) Exm = Exmemp
      endif
      if (Tm == 0.) then
        ald = ignatyuk(Zix, Nix, Exm, ibar)
        if (Exm <= max(2.25 / ald + P, 0.) + 0.11) Exm = Exmemp
        if (Exm > 3. * Exmemp) Exm = Exmemp
        i = int(Exm / dEx)
        if (i > 0) call pol1(i * dEx, (i + 1) * dEx, temprho(i), temprho(i + 1), Exm, Tm)
      endif
      if (E0m == 1.e-20) then
        i = int(Exm / dEx)
        if (i > 0) call pol1(i * dEx, (i + 1) * dEx, logrho(i), logrho(i + 1), Exm, &
          logrhomatch)
        rhomatch = exp(dble(logrhomatch))
        E0m = Exm - Tm * real(log(dble(Tm) * rhomatch))
      endif
      if (Tm == 0.) Tm = Tmemp
!
! Possible iteration after input defined adjustment
!
      if (T(Zix, Nix, ibar) == 0..and.Tadjust(Zix, Nix, ibar) /= 1.) then
        T(Zix, Nix, ibar) = Tadjust(Zix, Nix, ibar) * Tm
        cycle
      else
        T(Zix, Nix, ibar) = Tm
      endif
      if (E0(Zix, Nix, ibar) == 1.e-20 .and. E0adjust(Zix, Nix, ibar) /= 1.) then
        E0(Zix, Nix, ibar) = E0adjust(Zix, Nix, ibar) * E0m
        cycle
      else
        E0(Zix, Nix, ibar) = E0m
      endif
      Exmatch(Zix, Nix, ibar) = Exm
      exit
    enddo
  enddo
!
! Set theoretical value of D0
!
! dtheory: subroutine for theoretical calculation of average neutron spacings
!
  call dtheory(Zix, Nix, 0.)
  D0theo(Zix, Nix) = Dl(0)
  D1theo(Zix, Nix) = Dl(1)
  return
end subroutine densitymatch
! Copyright A.J. Koning 2021
