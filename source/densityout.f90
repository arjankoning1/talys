subroutine densityout(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of level density
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
!   flagparity      ! flag for non - equal parity distribution
!   gammald         ! gamma - constant for asymptotic level density parameter
!   Krotconstant    ! normalization constant for rotational enhancement
!   ldmodel         ! level density model
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
!   aldcrit         ! critical level density parameter
!   delta0          ! systematical pairing energy
!   Econd           ! condensation energy
!   edens           ! energy grid for tabulated level densities
!   ldexist         ! flag for existence of level density table
!   nendens         ! number of energies for level density grid
!   Nlast           ! last discrete level
!   scutoffdisc     ! spin cutoff factor for discrete level region
!   Tcrit           ! critical temperature
!   Ucrit           ! critical U
! Variables for fission parameters
!   nfisbar         ! number of fission barrier parameters
! Variables for masses
!   nucmass         ! mass of nucleus
!   S               ! separation energy
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=15) :: col(numJ+1)    ! header
  character(len=15) :: un(numJ+1)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=200) :: topline   ! topline
  character(len=8)  :: str(-1:1)     ! input line
  character(len=12) :: ldfile        ! level density file
  character(len=13) :: ldfileout     ! level density file
  character(len=12) :: ldstring      ! string for level density file
  character(len=25) :: model         ! string for level density model
  character(len=30) :: collstring    ! string
  integer           :: A             ! mass number of target nucleus
  integer           :: i             ! counter
  integer           :: Nav
  integer           :: ibeg
  integer           :: iend
  integer           :: i1(numlev2)
  integer           :: k             ! counter
  integer           :: Nk            ! counter
  integer           :: ibar          ! fission barrier
  integer           :: J             ! spin of level
  integer           :: ldmod         ! level density model
  integer           :: N             ! neutron number of residual nucleus
  integer           :: nex           ! excitation energy bin of compound nucleus
  integer           :: Ncol          !
  integer           :: Nix           ! neutron number index for residual nucleus
  integer           :: NL            ! last discrete level
  integer           :: NT            ! help variable
  integer           :: odd           ! odd (1) or even (0) nucleus
  integer           :: parity        ! parity
  integer           :: ploop         ! help variable
  integer           :: Z             ! charge number of target nucleus
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: ald           ! level density parameter
  real(sgl)         :: aldmatch      ! function to determine effective level density parameter
  real(sgl)         :: chi2D0        ! chi-square of D0
  real(sgl)         :: dEx           ! excitation energy bin for population arrays
  real(sgl)         :: Eex           ! excitation energy
  real(sgl)         :: Eex1(numlev2) ! excitation energy
  real(sgl)         :: CE            ! ratio D0theo/Dexp
  real(sgl)         :: CG            ! ratio D0theo/Dglob
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
  real(sgl)         :: FrmsD0        ! 
  real(sgl)         :: ErmsD0        ! 
  real(sgl)         :: rhoexp        ! discrete level density
  real(sgl)         :: x1(numlev2)
  real(sgl)         :: x2(numlev2)
  real(sgl)         :: x3(numlev2)
  real(sgl)         :: x4(numlev2)
  real(sgl)         :: x5(numlev2)
  real(sgl)         :: ignatyuk      ! function for energy dependent level density parameter a
  real(sgl)         :: Kcoll         ! total collective enhancement
  real(sgl)         :: Krot          ! rotational enhancement factor
  real(sgl)         :: Kvib          ! vibrational enhancement factor
  real(sgl)         :: P             ! pairing energy
  real(sgl)         :: sigma         ! help variable
  real(sgl)         :: spincut       ! spin cutoff factor
  real(sgl)         :: SS            ! separation energy
  real(sgl)         :: Tnuc          ! nuclear temperature
  real(dbl)         :: chi2sum       ! help variable
  real(dbl)         :: Frmssum       ! help variable
  real(dbl)         :: Ermssum       ! help variable
  real(dbl)         :: avdevsum      ! help variable
  real(dbl)         :: Ri            ! C/E
  real(dbl)         :: chi2          ! chi-square
  real(dbl)         :: Frms          ! root mean square
  real(dbl)         :: Erms          ! asymmetry
  real(dbl)         :: avdev         ! average deviation
  real(dbl)         :: dens          ! total level density
  real(dbl)         :: density       ! level density
  real(dbl)         :: densitytot    ! total level density
  real(dbl)         :: densitytotP   ! total level density per parity
  real(dbl)         :: ldtot         ! total level density
  real(dbl)         :: ldtotP        ! total level density per parity
  real(dbl)         :: Ncum          ! number of cumulative levels (integral of level density)
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
  ldmod = ldmodel(Zix, Nix)
  write(*, '(/" Level density parameters for Z=", i3, " N=", i3, " (", i3, a2, ") "/)')  Z, N, A, nuc(Z)
  if (ldmod == 1) model = "Constant temperature     "
  if (ldmod == 2) model = "Back-shifted Fermi Gas   "
  if (ldmod == 3) model = "Generalized superfluid   "
  if (ldmod == 4) model = "Goriely Skyrme           "
  if (ldmod == 5) model = "Hilaire-Goriely Skyrme   "
  if (ldmod == 6) model = "Hilaire-Goriely Gogny    "
  if (ldmod == 7) model = "BSKG3                    "
  write(*, '(" Model: ", a25)') model
  if (ldmod >= 4 .and. .not. ldexist(Zix, Nix, 0)) write(*, '(" Tables not available ")')
  if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, 0)) then
    write(*, '(" Collective enhancement: yes"/)')
  else
    write(*, '(" Collective enhancement: no"/)')
  endif
  if (flagfission) then
    write(*, '(21x, " g.s.     Fission barriers ")')
    write(*, '(29x, 3(5x, i1, 4x))') (ibar, ibar = 1, nfisbar(Zix, Nix))
    write(*, '()')
  endif
  write(*, '(" a(Sn)           :", f10.5)') alev(Zix, Nix)
  if (flagfission) write(*, '(" a-effective     :", 10x, 3f10.5)') (aldmatch(Zix, Nix, SS, ibar), ibar = 1, nfisbar(Zix, Nix))
!
! D0
!
  write(*, '(" Experimental D0 :", f18.2, " eV +- ", f15.5)') D0(Zix, Nix), dD0(Zix, Nix)
  write(*, '(" Theoretical D0  :", f18.2, " eV")') D0theo(Zix, Nix)
  write(*, '(" Theoretical D1  :", f18.2, " eV")') D1theo(Zix, Nix)
  write(*, '(" Av. res. energy :", f18.2, " eV")') Eavres * 1.e6
!
! Other parameters
!
! ignatyuk    : function for energy dependent level density parameter a
!
  P = pair(Zix, Nix)
  write(*, '(" Asymptotic a    :", f10.5)') alimit(Zix, Nix)
  write(*, '(" Damping gamma   :", f10.5)') gammald(Zix, Nix)
  write(*, '(" Pairing energy  :", f10.5)') P
  write(*, '(" Shell correction:", 4f10.5)') (deltaW(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
  write(*, '(" Last disc. level:", 4(7x, i3))') (Nlast(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
  write(*, '(" Nlow            :", 4(7x, i3))') (Nlow(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
  write(*, '(" Ntop            :", 4(7x, i3))') (Ntop(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
  if (ldmod == 1) then
    write(*, '(" Matching Ex     :", 4f10.5)') (Exmatch(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
    write(*, '(" Temperature     :", 4f10.5)') (T(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
    write(*, '(" E0              :", 4f10.5)') (E0(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
  endif
  write(*, '(" Adj. pair shift :", 4f10.5)') (Pshift(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
  if (ldmod == 3) then
    write(*, '(" Critical energy :", 4f10.5)') (Ucrit(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
    write(*, '(" Condensation en.:", 4f10.5)') (Econd(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
    write(*, '(" Critical temp.  :", f10.5)') Tcrit(Zix, Nix)
  endif
  write(*, '(" Discrete sigma  :", 4f10.5)') (sqrt(scutoffdisc(Zix, Nix, ibar)), ibar = 0, nfisbar(Zix, Nix))
  write(*, '(" Sigma (Sn)      :", 4f10.5)') (sqrt(spincut(Zix, Nix, ignatyuk(Zix, Nix, SS, ibar), SS, ibar, 0)), &
 &  ibar = 0, nfisbar(Zix, Nix))
  write(*, '(" Rhotot(Sn=", f5.2, "):", 1p, 4e10.3)') SS, (densitytot(Zix, Nix, SS, ibar, ldmod), ibar = 0, nfisbar(Zix, Nix))
  if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, 1)) then
    write(*, '(" beta2           :", f10.5)') beta2(Zix, Nix, 0)
    write(*, '(" Krotconstant    :", 4f10.5)') (Krotconstant(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
    write(*, '(" Ufermi          :", 4f10.5)') (Ufermi(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
    write(*, '(" cfermi          :", 4f10.5)') (cfermi(Zix, Nix, ibar), ibar = 0, nfisbar(Zix, Nix))
  endif
!
! ********************** Total level density ***************************
!
  odd = mod(A, 2)
  do ibar = 0, nfisbar(Zix, Nix)
    if (ibar == 0) then
      write(*, '(/" Level density per parity for ground state")')
    else
      write(*, '(/" Level density per parity for fission barrier", i3)') ibar
    endif
    write(*, '(" (Total level density also per parity)"/)')
    if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, ibar)) then
      write(*, '("    Ex     a    sigma   total ", 9("  JP= ", f4.1), "    Krot      Kvib      Kcoll")') &
 &      (real(J + 0.5 * odd), J = 0, 8)
    else
      write(*, '("    Ex     a    sigma   total ", 9("  JP= ", f4.1)/)') (real(J+0.5*odd), J = 0, 8)
    endif
!
! Tabulated level densities
!
    if (ldmod >= 4 .and. ldexist(Zix, Nix, ibar)) then
      if (ldmod == 4) then
        ploop = 1
      else
        ploop = - 1
      endif
      do parity = 1, ploop, - 2
        if (flagparity .and. parity == 1) write(*, '(/" Positive parity"/)')
        if (parity == -1) write(*, '(/" Negative parity"/)')
        do nex = 1, nendens(Zix, Nix)
          Eex = edens(nex)
          write(*, '(1x, f6.2, 14x, 11es10.3)') Eex, densitytotP(Zix, Nix, Eex, parity, ibar, ldmod), &
 &          (density(Zix, Nix, Eex, real(J + 0.5 * odd), parity, ibar, ldmod), J = 0, 8)
        enddo
      enddo
      write(*, '(/" Normalization:")')
      write(*, '("        ctable=", f10.5)') ctable(Zix, Nix, ibar)
      write(*, '("        ptable=", f10.5)') ptable(Zix, Nix, ibar)
      write(*, '("      s2adjust=", f10.5)') s2adjust(Zix, Nix, ibar)
    else
!
! Analytical level densities
!
! colenhance: subroutine for collective enhancement
!
      do nex = 1, nendens(Zix, Nix)
        Eex = edens(nex)
        ald = ignatyuk(Zix, Nix, Eex, ibar)
        if (ldmod == 3 .and. Eex < Ucrit(Zix, Nix, ibar) - P - Pshift(Zix, Nix, ibar)) ald = aldcrit(Zix, Nix, ibar)
        if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, ibar)) then
          call colenhance(Zix, Nix, Eex, ald, ibar, Krot, Kvib, Kcoll)
          write(collstring, '(3es10.3)') Krot, Kvib, Kcoll
        else
          collstring = ' '
        endif
        write(*, '(1x, f6.2, 2f7.3, 10es10.3, a30)') Eex, ald, sqrt(spincut(Zix, Nix, ald, Eex, ibar, 0)), &
 &        densitytotP(Zix, Nix, Eex, 1, ibar, ldmod), &
 &        (density(Zix, Nix, Eex, real(J + 0.5 * odd), 1, ibar, ldmod), J = 0, 8), collstring
      enddo
    endif
  enddo
!
! Cumulative number of discrete levels vs. integrated level density
!
! Output given in general output file and on separate files.
!
  massstring='   '
  write(massstring,'(i3)') A
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  do ibar = 0, nfisbar(Zix, Nix)
    NL = Nlow(Zix, Nix, ibar)
    NT = Ntop(Zix, Nix, ibar)
    write(*, '(/" Discrete levels versus total level density"/)')
    write(*, '("   Energy Level   N_cumulative"/)')
    if (filedensity) then
      ldfile = 'ld000000.tot'
      if (ibar > 0) write(ldfile(10:12), '(i3.3)') ibar
      write(ldfile(3:8), '(2i3.3)') Z, A
      open (unit = 1, file = ldfile, status = 'replace')
      if (ibar > 0) then
        ldstring = ' Barrier    '
        write(ldstring(10:10), '(i1)') ibar
      else
        ldstring = '            '
      endif
      quantity='level density'
      topline=trim(finalnuclide)//' '//trim(quantity)
      col = ''
      un = ''
      col(1)='E'
      un(1)='MeV'
      col(2)='Level'
      col(3)='N_cumulative'
      col(4)='Total_LD'
      un(4)='MeV^-1'
      col(5)='Exp_LD'
      un(5)='MeV^-1'
      col(6)='a'
      un(6)='MeV^-1'
      col(7)='Sigma'
      call write_header(topline,source,user,date,oformat)
      call write_residual(Z,A,finalnuclide)
      write(1,'("# parameters:")')
      call write_integer(2,'ldmodel keyword',ldmod)
      call write_char(2,'level density model',model)
      if (ldmod <= 3) then
        Ncol=7
        if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, ibar)) then
          call write_char(2,'collective enhancement','y')
        else
          call write_char(2,'collective enhancement','n')
        endif
        call write_real(2,'a(Sn) [MeV^-1]',alev(Zix,Nix))
        call write_real(2,'asymptotic a [MeV^-1]',alimit(Zix,Nix))
        call write_real(2,'shell correction [MeV]',deltaW(Zix,Nix,ibar))
        call write_real(2,'damping gamma',gammald(Zix,Nix))
        call write_real(2,'pairing energy [MeV]',P)
        call write_real(2,'adjusted pairing shift [MeV]',Pshift(Zix, Nix, ibar))
        call write_real(2,'separation energy [MeV]',SS)
        call write_real(2,'discrete spin cutoff parameter',scutoffdisc(Zix, Nix, ibar))
        call write_real(2,'spin cutoff parameter(Sn)',spincut(Zix, Nix, ignatyuk(Zix, Nix, SS, ibar), SS, ibar, 0))
        if (ldmod == 1) then
          call write_real(2,'matching energy [MeV]',Exmatch(Zix, Nix, ibar))
          call write_real(2,'temperature [MeV]',T(Zix, Nix, ibar))
          call write_real(2,'E0 [MeV]',E0(Zix, Nix, ibar))
        endif
        if (ldmod == 3) then
          call write_real(2,'Delta0',delta0(Zix, Nix))
          call write_real(2,'critical a [MeV^-1]',aldcrit(Zix, Nix,ibar))
          call write_real(2,'critical energy [MeV]',Ucrit(Zix, Nix,ibar))
          call write_real(2,'condensation energy [MeV]',Econd(Zix, Nix,ibar))
          call write_real(2,'critical temperature [MeV]',Tcrit(Zix, Nix))
        endif
      else
        Ncol=5
      endif
      call write_integer(2,'number of excited levels',nlevmax2(Zix,Nix))
      call write_integer(2,'Nlow',Nlow(Zix, Nix, ibar))
      call write_integer(2,'Ntop',Ntop(Zix, Nix, ibar))
      call write_real(2,'ctable',ctable(Zix, Nix, ibar))
      call write_real(2,'ptable',ptable(Zix, Nix, ibar))
      write(1,'("# observables:")')
      Dtheo=D0theo(Zix, Nix)
      Dexp=D0(Zix, Nix)
      dDexp=dD0(Zix, Nix)
      Dglob=D0global(Zix, Nix)
      dDglob=dD0global(Zix, Nix)
      if (dDexp <= .1e-30) then
        chi2D0 = 0.
      else
        chi2D0 = (abs(Dtheo - Dexp) / dDexp) **2
      endif
!
! Cumulative probability derived as follows:
!
!          cdf=0.5*(1.+erf(x))
!          pp=2.*(cdf-0.5)
!
      FrmsD0 = 0.
      ErmsD0 = 0.
      CE = 0.
      if (Dexp > 1.e-30) then
        CE = Dtheo / Dexp
        if (Dtheo > 0. .and. Dexp > 0.) then
          if (dDexp > 0.) then
            x = (Dtheo - Dexp)/(dDexp*sqrt(2.))
            if (Dtheo > Dexp) then
              pp = erf(x)
            else
              pp = -erf(x)
            endif
            Ri = 1. + (CE-1.)*pp
          else
            Ri = CE
          endif
          if (Ri >= 1.) then
            FrmsD0 = Ri
          else
            FrmsD0 = 1./Ri
          endif
          ErmsD0 = Ri
        endif
      endif
      CG = 0.
      if (Dglob > 1.e-30) CG = Dtheo / Dglob
      call write_real(2,'experimental D0 [eV]',Dexp)
      call write_real(2,'experimental D0 unc. [eV]',dDexp)
      call write_real(2,'global D0 [eV]',Dglob)
      call write_real(2,'global D0 unc. [eV]',dDglob)
      call write_real(2,'theoretical D0 [eV]',Dtheo)
      call write_real(2,'Chi-2 D0',chi2D0)
      call write_real(2,'C/E D0',CE)
      call write_real(2,'Frms D0',FrmsD0)
      call write_real(2,'Erms D0',ErmsD0)
      call write_real(2,'C/G D0',CG)
    endif
    Ncum = 0
    chi2 = 0.
    Frms = 0.
    Erms = 0.
    avdev = 0.
    k = 0
    Eex1 = 0.
    i1 = 0
    x1 = 0.
    x2 = 0.
    x3 = 0.
    x4 = 0.
    x5 = 0.
    chi2sum = 0.
    Frmssum = 0.
    Ermssum = 0.
    avdevsum = 0.
    do i = 1, nlevmax2(Zix, Nix) - 1
      if (edis(Zix, Nix, i + 1) == 0.) cycle
      Eex = 0.5 * (edis(Zix, Nix, i) + edis(Zix, Nix, i - 1))
      dEx = edis(Zix, Nix, i) - edis(Zix, Nix, i - 1)
      dens = densitytot(Zix, Nix, Eex, ibar, ldmod)
      Ncum = Ncum + dens * dEx
      if (i == NL) Ncum = real(NL)
      rhoexp = 0.
      Nav = 10
      ibeg = max(i - Nav/2,0)
      iend = min(i + Nav/2,nlevmax2(Zix, Nix))
      Nav = iend - ibeg + 1
      Ea = edis(Zix, Nix, ibeg)
      Eb = edis(Zix, Nix, iend)
      Erange = Eb - Ea
      if (Erange > 0.) rhoexp = Nav / Erange
      Eex = edis(Zix, Nix, i)
      if (Ncum < 1.e9) write(*, '(1x, f8.4, i4, f12.3)') Eex, i, Ncum
      if (filedensity) then
        k = k + 1
        if (ldmod <= 3) then
          ald = ignatyuk(Zix, Nix, Eex, ibar)
          if (ldmod == 3 .and. Eex < Ucrit(Zix, Nix, ibar) - P - &
            Pshift(Zix, Nix, ibar)) ald = aldcrit(Zix, Nix, ibar)
          dens = densitytot(Zix, Nix, Eex, ibar, ldmod)
          sigma = sqrt(spincut(Zix, Nix, ald, Eex, ibar, 0))
          x4(k) = ald
          x5(k) = sigma
        endif
        Eex1(k) = Eex
        i1(k) = i
        x1(k) = Ncum
        x2(k) = dens
        x3(k) = rhoexp
      endif
      if (i >= NL .and. i <= NT) then
        chi2sum = chi2sum + ((Ncum - i) **2 ) / sqrt(real(i)) / max(NT-NL,1)
        Ri = Ncum / i
        Frmssum = Frmssum + log(Ri)**2 
        Ermssum = Ermssum + log(Ri) 
        avdevsum = avdevsum + abs(Ncum - i) 
      endif
    enddo
    denom = real(NT - NL)
    chi2 = chi2sum 
    if (denom > 0.) then
      Frms = exp(sqrt(Frmssum / denom))
      Erms = exp(Ermssum / denom)
      avdev = avdevsum / denom
    endif
    Nk = k
    if (filedensity) then
      call write_double(2,'Chi-2 per level',chi2)
      call write_double(2,'Frms per level',Frms)
      call write_double(2,'Erms per level',Erms)
      call write_double(2,'average deviation per level',avdev)
      call write_datablock(quantity,Ncol,Nk,col,un)
      do k = 1, Nk
        if (ldmod <= 3) then
          write(1, '(es15.6, i6, 9x, 5es15.6)') Eex1(k), i1(k), x1(k), x2(k), x3(k), x4(k), x5(k)
        else
          write(1, '(es15.6, i6, 9x, 3es15.6)') Eex1(k), i1(k), x1(k), x2(k), x3(k)
        endif
      enddo
!
! Spin distribution
!
      quantity='spin distribution'
      Ncol = 3
      un=''
      col(1) = 'Spin'
      col(2) = 'R (parity -)'
      col(3) = 'R (parity +)'
      Rdist = 0.
      do parity = -1, 1, 2
        do nex = 1, nendens(Zix, Nix)
          Eex = edens(nex)
          ldtotP = densitytotP(Zix, Nix, Eex, parity, ibar, ldmod)
          do J = 0, numJ
            Rdist(nex, parity, J) = density(Zix, Nix, Eex, real(J + 0.5 * odd), parity, ibar, ldmod) / ldtotP
          enddo
        enddo
      enddo
      do nex = 1, nendens(Zix, Nix)
        write(1,'("# parameters:")')
        call write_real(2,'Excitation energy [MeV]',edens(nex))
        call write_datablock(quantity,Ncol,numJ+1,col,un)
        do J = 0, numJ
          write(1, '(6x,f4.1,5x,2es15.6)') real(J + 0.5 * odd),(Rdist(nex, parity, J), parity = -1, 1, 2)
        enddo
      enddo
      close (1)
    endif
  enddo
!
! Level densities per parity on separate files, as in tabulated format
!
  if (filedensity) then
    ibar = 0
    str(-1) = 'Negative'
    str(1) = 'Positive'
    do parity = 1, -1, -2
      ldfileout = 'nld000000.tab'
      write(ldfileout(4:9), '(2i3.3)') Z, A
      open (unit = 2, status = 'unknown', file = ldfileout)
      write(2, '(20x, 96("*"))')
      write(2, '(20x, "*  Z=", i3, " A=", i3, ": ", a8, "-Parity Spin-dependent Level Density [MeV-1] for ", a2, i3, &
 &      " and ldmodel=", i2, "  *")') Z, A, str(parity), nuc(Z), A, ldmod
      write(2, '(20x, 96("*"))')
      if (mod(A, 2) == 0) then
        write(2, '(" U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT  ", 31("   J=", i2.2, 2x, :))') (J, J = 0, 29)
      else
        write(2, '(" U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT  ", 31("  J=", i2.2, "/2", 1x, :))') (J, J = 1, 59, 2)
      endif
      Ncum = 0.
      dEx = edens(1)
      do nex = 1, nendens(Zix, Nix)
        Eex = edens(nex)
        Tnuc = sqrt(Eex / alev(Zix, Nix))
        if (nex > 1) dEx = Eex - edens(nex - 1)
        dens = densitytotP(Zix, Nix, Eex, parity, ibar, ldmod)
        Ncum = Ncum + dens * dEx
        ldtot = 0.
        do J = 0, numJ
          ldtot = ldtot + (2. * J + 1) * density(Zix, Nix, Eex, real(J + 0.5 * odd), parity, ibar, ldmod)
        enddo
        ldtotP = densitytotP(Zix, Nix, Eex, parity, ibar, ldmod)
        write(2, '(1x, f6.2, f7.3, 1x, 1p, 33e9.2)') Eex, Tnuc, Ncum, ldtotP, ldtot, &
 &        (density(Zix, Nix, Eex, real(J + 0.5 * odd), 1, ibar, ldmod), J = 0, 29)
      enddo
      write(2, * )
    enddo
    close(2)
  endif
  return
end subroutine densityout
! Copyright A.J. Koning 2021
