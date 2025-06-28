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
  character(len=15) :: col(numJ+6)    ! header
  character(len=15) :: un(numJ+6)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=200) :: topline   ! topline
  character(len=8)  :: str(-1:1)     ! input line
  character(len=30) :: ldfile        ! level density file
  character(len=30) :: ldfileout     ! level density file
  character(len=25) :: model         ! string for level density model
  integer           :: A             ! mass number of target nucleus
  integer           :: i             ! counter
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
  integer           :: Z             ! charge number of target nucleus
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: ald           ! level density parameter
  real(sgl)         :: aldmatch      ! function to determine effective level density parameter
  real(sgl)         :: dEx           ! excitation energy bin for population arrays
  real(sgl)         :: Eex           ! excitation energy
  real(sgl)         :: Eex1(numlev2) ! excitation energy
  real(sgl)         :: denom
  real(sgl)         :: x1(numlev2)
  real(sgl)         :: x2(numlev2)
  real(sgl)         :: x3(numlev2)
  real(sgl)         :: x4(numlev2)
  real(sgl)         :: x5(numlev2)
  real(sgl)         :: x6(numlev2)
  real(sgl)         :: x7(numlev2)
  real(sgl)         :: x8(numlev2)
  real(sgl)         :: rJ
  real(sgl)         :: ignatyuk      ! function for energy dependent level density parameter a
  real(sgl)         :: Kcoll         ! total collective enhancement
  real(sgl)         :: Krot          ! rotational enhancement factor
  real(sgl)         :: Kvib          ! vibrational enhancement factor
  real(sgl)         :: P             ! pairing energy
  real(sgl)         :: sigma         ! help variable
  real(sgl)         :: spincut       ! spin cutoff factor
  real(sgl)         :: SS            ! separation energy
  real(sgl)         :: Tnuc          ! nuclear temperature
  real(dbl)         :: dens(0:numJ)  ! level density
  real(dbl)         :: denstot       ! total level density
  real(dbl)         :: denstotP      ! total level density per parity
  real(dbl)         :: density       ! level density
  real(dbl)         :: densitytot    ! total level density
  real(dbl)         :: densitytotP   ! total level density per parity
  real(dbl)         :: ldtot         ! total level density
  real(dbl)         :: ldtotP        ! total level density per parity
  real(dbl)         :: NC            ! cumulative number of levels
  real(dbl)         :: Rdist(numdens, -1:1, 0:numJ) ! spin distribution
!
! ********************** Level density parameters **********************
!
! aldmatch   : function to determine effective level density parameter
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  odd = mod(A, 2)
  SS = S(Zix, Nix, 1)
  P = pair(Zix, Nix)
  massstring='   '
  write(massstring,'(i3)') A
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  ldmod = ldmodel(Zix, Nix)
  if (ldmod == 1) model = "Constant-Temperature     "
  if (ldmod == 2) model = "Back-shifted-Fermi-Gas   "
  if (ldmod == 3) model = "Generalized-Superfluid   "
  if (ldmod == 4) model = "Goriely-Skyrme           "
  if (ldmod == 5) model = "Hilaire-Goriely-Skyrme   "
  if (ldmod == 6) model = "Hilaire-Goriely-Gogny    "
  if (ldmod == 7) model = "BSKG3-Combinatorial      "
  if (ldmod == 8) model = "D1M-QRPA-BE              "
  write(*, '(/" Level densities for Z=", i3, " N=", i3, " (",a,") "/)')  Z, N, trim(finalnuclide)
  write(*, '(" Total level density"/)') 
!
! ********************** Total level density ***************************
!
! Cumulative number of discrete levels vs. integrated level density
!
! Output given in general output file and on separate files.
!
  do ibar = 0, nfisbar(Zix, Nix)
    NL = Nlow(Zix, Nix, ibar)
    NT = Ntop(Zix, Nix, ibar)
    ldfile = 'ld000000.gs'
    if (ibar > 0) write(ldfile(10:12), '("b",i2.2)') ibar
    write(ldfile(3:8), '(2i3.3)') Z, A
    open (unit = 1, file = ldfile, status = 'replace')
    topline=trim(finalnuclide)//' level density'
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
    call write_integer(2,'fission barrier',ibar)
    call write_integer(2,'ldmodel keyword',ldmod)
    call write_char(2,'level density model',model)
    if (ldmod <= 3) then
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
      call write_real(2,'discrete spin cutoff parameter',scutoffdisc(Zix, Nix, ibar))
      call write_real(2,'spin cutoff parameter(Sn)',spincut(Zix, Nix, ignatyuk(Zix, Nix, SS, ibar), SS, ibar, 0))
      if (ldmod == 1) then
        call write_real(2,'matching energy [MeV]',Exmatch(Zix, Nix, ibar))
        call write_real(2,'temperature [MeV]',T(Zix, Nix, ibar))
        call write_real(2,'E0 [MeV]',E0(Zix, Nix, ibar))
      endif
      if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, 1)) then
        call write_real(2,'beta2',beta2(Zix, Nix, 0))
        call write_real(2,'Krotconstant',Krotconstant(Zix, Nix, ibar))
        call write_real(2,'Ufermi',Ufermi(Zix, Nix, ibar))
        call write_real(2,'cfermi',cfermi(Zix, Nix, ibar))
      endif
      if (ibar > 0) call write_real(2,'a-effective [MeV^-1]',aldmatch(Zix, Nix, SS, ibar)) 
      if (ldmod == 3) then
        call write_real(2,'Delta0',delta0(Zix, Nix))
        call write_real(2,'critical a [MeV^-1]',aldcrit(Zix, Nix,ibar))
        call write_real(2,'critical energy [MeV]',Ucrit(Zix, Nix,ibar))
        call write_real(2,'condensation energy [MeV]',Econd(Zix, Nix,ibar))
        call write_real(2,'critical temperature [MeV]',Tcrit(Zix, Nix))
      endif
      Ncol=7
      if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, ibar)) then
        Ncol=10
        col(8)='Krot'
        col(9)='Kvib'
        col(10)='Kcoll'
      endif
    else
      Ncol=5
    endif
    call write_real(2,'separation energy [MeV]',SS)
    call write_double(2,'Rhotot(Sn) [MeV^-1]',densitytot(Zix, Nix, SS, ibar, ldmod))
    call write_integer(2,'number of excited levels',nlevmax2(Zix,Nix))
    call write_integer(2,'Nlow',Nlow(Zix, Nix, ibar))
    call write_integer(2,'Ntop',Ntop(Zix, Nix, ibar))
    call write_real(2,'ctable',ctable(Zix, Nix, ibar))
    call write_real(2,'ptable',ptable(Zix, Nix, ibar))
    if (ibar == 0) then
      call write_real(2,'experimental D0 [eV]',D0(Zix, Nix))
      call write_real(2,'experimental D0 unc. [eV]',dD0(Zix, Nix))
      call write_real(2,'global D0 [eV]',D0global(Zix, Nix))
      call write_real(2,'global D0 unc. [eV]',dD0global(Zix, Nix))
      call write_real(2,'theoretical D0 [eV]',D0theo(Zix, Nix))
      call write_real(2,'Chi-2 D0',chi2D0(Zix, Nix))
      call write_real(2,'C/E D0',CED0(Zix, Nix))
      call write_real(2,'Frms D0',FrmsD0(Zix, Nix))
      call write_real(2,'Erms D0',ErmsD0(Zix, Nix))
      call write_real(2,'C/G D0',CGD0(Zix, Nix))
      call write_real(2,'experimental D1 [eV]',D1r(Zix, Nix))
      call write_real(2,'experimental D1 unc. [eV]',dD1r(Zix, Nix))
      call write_real(2,'theoretical D1 [eV]',D1theo(Zix, Nix))
      k = 0
      Eex1 = 0.
      i1 = 0
      x1 = 0.
      x2 = 0.
      x3 = 0.
      x4 = 0.
      x5 = 0.
      x6 = 0.
      x7 = 0.
      x8 = 0.
      do i = 1, nlevmax2(Zix, Nix)
        if (edis(Zix, Nix, i) == 0.) cycle
        Eex = 0.5 * (edis(Zix, Nix, i) + edis(Zix, Nix, i - 1))
        dEx = edis(Zix, Nix, i) - edis(Zix, Nix, i - 1)
        denstot = densitytot(Zix, Nix, Eex, ibar, ldmod)
        Eex = edis(Zix, Nix, i)
        k = k + 1
        if (ldmod <= 3) then
          ald = ignatyuk(Zix, Nix, Eex, ibar)
          if (ldmod == 3 .and. Eex < Ucrit(Zix, Nix, ibar) - P - Pshift(Zix, Nix, ibar)) ald = aldcrit(Zix, Nix, ibar)
          denstot = densitytot(Zix, Nix, Eex, ibar, ldmod)
          sigma = sqrt(spincut(Zix, Nix, ald, Eex, ibar, 0))
          x4(k) = ald
          x5(k) = sigma
          if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, ibar)) then
            call colenhance(Zix, Nix, Eex, ald, ibar, Krot, Kvib, Kcoll)
            x6(k) = Krot
            x7(k) = Kvib
            x8(k) = Kcoll
          endif
        endif
        Eex1(k) = Eex
        i1(k) = i
        x1(k) = Ncum(Zix, Nix, i)
        x2(k) = denstot
        x3(k) = rhoexp(Zix, Nix, i)
      enddo
      denom = real(NT - NL)
      Nk = k
      call write_double(2,'Chi-2 per level',chi2lev(Zix, Nix))
      call write_double(2,'Frms per level',Frmslev(Zix, Nix))
      call write_double(2,'Erms per level',Ermslev(Zix, Nix))
      call write_double(2,'average deviation per level',avdevlev(Zix, Nix))
      quantity='total level density'
      call write_quantity(quantity)
      call write_datablock(Ncol,Nk,col,un)
      do k = 1, Nk
        if (ldmod <= 3) then
          if (flagcol(Zix, Nix) .and. .not. ldexist(Zix, Nix, ibar)) then
            write(1, '(es15.6, i6, 9x, 8es15.6)') Eex1(k), i1(k), x1(k), x2(k), x3(k), x4(k), x5(k), x6(k), x7(k), x8(k)
          else
            write(1, '(es15.6, i6, 9x, 5es15.6)') Eex1(k), i1(k), x1(k), x2(k), x3(k), x4(k), x5(k)
          endif
        else
          write(1, '(es15.6, i6, 9x, 3es15.6)') Eex1(k), i1(k), x1(k), x2(k), x3(k)
        endif
      enddo
    endif
!
! Level densities per parity on separate files
!
    write(*, '(" Level density per spin and parity"/)') 
    quantity='level density'
    col=''
    col(1)='E'
    col(2)='T'
    col(3)='N_cumulative'
    col(4)='rho_observed'
    col(5)='rho_total'
    do J = 0, numJ
      col(J+6)='rho(J)=     '
      write(col(J+6)(9:12),'(f4.1)') J+0.5*odd 
    enddo      
    un='MeV^-1'
    un(1)='MeV'
    un(2)='MeV'
    un(3)=''
    Ncol=numJ+6
    str(-1) = 'Negative'
    str(1) = 'Positive'
    do parity = 1, -1, -2
      call write_quantity(quantity)
      call write_integer(2,'Parity',parity)
      call write_datablock(Ncol,nendens(Zix,Nix),col,un)
      if (parity == 1) then
        ldfileout = 'nld000000.tab'
        write(ldfileout(4:9), '(2i3.3)') Z, A
        if (ibar > 0) write(ldfileout(14:17), '("_b",i2.2)') ibar
        open (unit = 2, status = 'replace', file = ldfileout)
      endif
      write(2, '(20x, 96("*"))')
      write(2, '(20x, "*  Z=", i3, " A=", i3, ": ", a8, "-Parity Spin-dependent Level Density [MeV-1] for ", a2, i3, &
 &      " and ldmodel=", i2, "  *")') Z, A, str(parity), nuc(Z), A, ldmod
      write(2, '(20x, 96("*"))')
      if (mod(A, 2) == 0) then
        write(2, '(" U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT  ", 31("   J=", i2.2, 2x, :))') (J, J = 0, 29)
      else
        write(2, '(" U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT  ", 31("  J=", i2.2, "/2", 1x, :))') (J, J = 1, 59, 2)
      endif
      Nc = 0.
      dEx = edens(1)
      do nex = 1, nendens(Zix, Nix)
        Eex = edens(nex)
        denstotP = densitytotP(Zix, Nix, Eex, parity, ibar, ldmod)
        if (ldmod > 3) then
          Tnuc = ldtableT(Zix, Nix, nex, parity, ibar)
          Nc = ldtableN(Zix, Nix, nex, parity, ibar)
        else
          Tnuc = sqrt(Eex / alev(Zix, Nix))
          if (nex > 1) then
            dEx = Eex - edens(nex - 1)
          else
            dEx = Eex
          endif
          Nc = Nc + denstotP * dEx
        endif
        ldtot = 0.
        do J = 0, numJ
          rJ = real(J + 0.5 * odd)
          dens(J) = density(Zix, Nix, Eex, rJ, parity, ibar, ldmod)
          ldtot = ldtot + (2. * rJ + 1) * dens(J)
        enddo
        write(2, '(1x, f6.2, f7.3, 1x, 33es9.2)') Eex, Tnuc, Nc, denstotP, ldtot, (dens(J), J = 0, 29)
        write(1, '(47es15.6)') Eex, Tnuc, Nc, denstotP, ldtot, (dens(J), J = 0, numJ)
      enddo
      write(2, * )
    enddo
    close(2)
!
! Spin distribution
!
    write(*, '(/" Level density spin distribution"/)') 
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
      call write_quantity(quantity)
      call write_real(2,'Excitation energy [MeV]',edens(nex))
      call write_datablock(Ncol,numJ+1,col,un)
      do J = 0, numJ
        write(1, '(6x,f4.1,5x,2es15.6)') real(J + 0.5 * odd),(Rdist(nex, parity, J), parity = -1, 1, 2)
      enddo
    enddo
    close(1)
    call write_outfile(ldfile,flagoutall)
  enddo
  return
end subroutine densityout
! Copyright A.J. Koning 2021
