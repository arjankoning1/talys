subroutine gammaout(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of gamma-ray strength functions, transmission
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2022-03-20: Added output for Gamma_gamma
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
! All global variables
!   numen         ! maximum number of energy points
! Variables for level density
!   D0            ! s - wave resonance spacing in eV
! Variables for gamma rays
!   egr           ! energy of GR
!   epr           ! energy of PR
!   etable        ! constant to adjust tabulated strength functions
!   filepsf       ! flag for photon strength functions on separate files
!   flagupbend    ! flag for low-energy upbend of photon strength function
!   ftable        ! constant to adjust tabulated strength functions
!   gamgam        ! total radiative width in eV
!   gammax        ! number of l - values for gamma multipolarity
!   ggr           ! width of GR
!   gpr           ! width of PR
!   sgr           ! strength of GR
!   strength      ! E1 strength function model
!   strengthM1    ! model for M1 gamma - ray strength function
!   tpr           ! strength of PR
!   upbend        ! properties of the low - energy upbend of given multipolarity
!   wtable        ! constant to adjust tabulated strength functions
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
!   egrid         ! outgoing energy grid
!   Einc          ! incident energy in MeV
! Variables for energies
!   eend          ! last energy point of energy grid
! Variables for inverse channel data
!   Tjl           ! transmission coefficient per particle, energy, spin and l - value
!   xsreac        ! reaction cross section
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   NN            ! neutron number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Constants
!   nuc           ! symbol of nucleus
!  Variables for gamma-ray strength functions
!   kgr           ! constant for gamma - ray strength function
!   ngr           ! number of GR
!   nTqrpa        ! number of temperatures for QRPA
! Variables for resonance parameters
!   D0theo        ! mean s - wave resonance spacing
!   D1theo        ! mean p - wave resonance spacing
!   dD0           ! uncertainty in D0
!   dgamgam       ! uncertainty in gamgam
!   gamgamth      ! theoretical total radiative width
!   Eavres        ! number of resonances
!   swaveth       ! theoretical strength function for s - wave
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=7)  :: crossfile !
  character(len=15) :: col(2)    ! header
  character(len=15) :: un(2)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline   ! topline
  character(len=1 ) :: radtype      ! radiation type
  character(len=12) :: psffile      ! photon strength function file
  character(len=25) :: modelM1      ! string for gamma-ray strength function
  character(len=25) :: modelE1      ! string for gamma-ray strength function
  character(len=25) :: psfmodel     ! string for gamma-ray strength function
  integer           :: A            ! mass number of target nucleus
  integer           :: i            ! counter
  integer           :: irad         ! counter
  integer           :: l            ! multipolarity
  integer           :: MM           ! model
  integer           :: N            ! neutron number of residual nucleus
  integer           :: Ncol       ! counter
  integer           :: Ncomp        ! neutron number index for compound nucleus
  integer           :: nen          ! energy counter
  integer           :: nenb         ! begin energy
  integer           :: nene         ! end energy
  integer           :: Npsf         ! number of energy points for PSF
  integer           :: Z            ! charge number of target nucleus
  integer           :: Zcomp        ! proton number index for compound nucleus
  real(sgl)         :: CEgamgam     ! C/E for Gamma_gamma
  real(sgl)         :: e            ! energy
  real(sgl)         :: Epsf(numen)  ! number of energy points
  real(sgl)         :: fstrength    ! gamma-ray strength function
!
! ***** Gamma-ray strength functions and transmission coefficients *****
!
! fstrength : gamma-ray strength function
!
  Z = ZZ(Zcomp, Ncomp, 0)
  N = NN(Zcomp, Ncomp, 0)
  A = AA(Zcomp, Ncomp, 0)
  if (filepsf) then
    i = 0
    do nen = ebegin(0), eend(0)
      i = i + 1
      Epsf(i) = egrid(nen)
    enddo
    Npsf = i
    if (Npsf > 0) then
      if (Epsf(Npsf) < 30.) then
        nenb = int(2. * (Epsf(Npsf) + 0.5))
        nene = 60
        do nen = nenb, nene
          i = i + 1
          if (i <= numen) Epsf(i) = 0.5 * nen
        enddo
        Npsf = min(i, numen)
      endif
    endif
  endif
  write(*, '(/" ########## GAMMA STRENGTH FUNCTIONS, TRANSMISSION", " COEFFICIENTS AND CROSS SECTIONS ##########")')
  write(*, '(/" Gamma-ray information for Z=", i3, " N=", i3, " (", i3, a2, ") "/)') Z, N, A, nuc(Z)
  write(*, '(" S-wave strength function parameters:"/)')
  write(*, '(" Exp. total radiative width=", f10.5, " eV +/-", f8.5, " Theor. total radiative width for l=0:", f15.5, " eV")') &
 &  gamgam(Zcomp, Ncomp), dgamgam(Zcomp, Ncomp), gamgamth(Zcomp, Ncomp, 0)
  write(*, '(53x, " Theor. total radiative width for l=0:", f15.5, " eV")') gamgamth(Zcomp,Ncomp,0)
  write(*, '(53x, " Theor. total radiative width for l=1:", f15.5, " eV")') gamgamth(Zcomp, Ncomp, 1)
  write(*, '(" Exp. D0                   =", f10.2, " eV +/-", f8.2, " Theor. D0                   =", f15.5, " eV")') &
 &  D0(Zcomp, Ncomp), dD0(Zcomp, Ncomp), D0theo(Zcomp, Ncomp)
  write(*, '(53x, " Theor. D1                   =", f15.5, " eV")') D1theo(Zcomp, Ncomp)
  write(*, '(" Theor. S-wave strength f. =", f10.5, "E-4")') 1.e4 * swaveth(Zcomp, Ncomp)
  write(*,'(" Average resonance energy  =", f13.2, " eV")') Eavres * 1.e6
  write(*, '(/" Incident energy: E[MeV]=", f8.3)') Einc
  if (strength == 1) modelE1 = "Kopecky-Uhl              "
  if (strength == 2) modelE1 = "Brink-Axel               "
  if (strength == 3) modelE1 = "Goriely HFbcs tables     "
  if (strength == 4) modelE1 = "Goriely HFB tables       "
  if (strength == 5) modelE1 = "Goriely Hybrid model     "
  if (strength == 6) modelE1 = "Goriely T-dep. HFB Tables"
  if (strength == 7) modelE1 = "Goriely T-dep. RMF Tables"
  if (strength == 8) modelE1 = "Gogny D1M HFB+QRPA Tables"
  if (strength == 9) modelE1 = "IAEA-CRP SMLO 2019 Tables"
  if (strength == 10) modelE1 = "BSk27+QRPA 2018 Tables   "
  write(*, '(/" Gamma-ray strength function model for E1: ", a25)') modelE1
  if (strengthM1 == 1) modelM1 = "RIPL-1                   "
  if (strengthM1 == 2) modelM1 = "RIPL-2                   "
  if (strengthM1 == 3) modelM1 = "IAEA GSF CRP (2018)      "
  if (strengthM1 == 4) modelM1 = "RIPL-2+Scissors Kawano   "
  if (strengthM1 == 8) modelM1 = "Gogny D1M HFB+QRPA Tables"
  if (strengthM1 == 10) modelM1 = "BSk27+QRPA 2018 Tables   "
  write(*, '(/" Gamma-ray strength function model for M1: ", a25)') modelM1
  if (strength == 3 .or. strength == 4 .or. strength >= 6) then
    write(*, '(/" Adjustable parameters for E1: etable=", f10.5, " ftable=", f10.5, " wtable=", f10.5, " number of T=", i3)') &
 &    etable(Zcomp, Ncomp, 1, 1), ftable(Zcomp, Ncomp, 1, 1), wtable(Zcomp, Ncomp, 1, 1), nTqrpa
  endif
  if (strengthM1 == 8 .or. strengthM1 == 10) then
    write(*, '(/" Adjustable parameters for M1: etable=", f10.5, " ftable=", f10.5, " wtable=", f10.5)') &
 &    etable(Zcomp, Ncomp, 0, 1), ftable(Zcomp, Ncomp, 0, 1), wtable(Zcomp, Ncomp, 0, 1)
  endif
  if (flagupbend) then
    write(*,'(/" Inclusion of an E1 upbend C x U*/ (1+exp(E-eta)) with  C=", es10.2, " eta=", f8.2)') &
 &    upbend(Zcomp, Ncomp, 1, 1, 1), upbend(Zcomp, Ncomp, 1, 1, 2)
    write(*,'(/" Inclusion of an M1 upbend C exp(-F*|beta2|) exp(-eta*E) with C=",es10.2," eta=",f8.2," F=",f8.2)') &
 &    upbend(Zcomp, Ncomp, 0, 1, 1), upbend(Zcomp, Ncomp, 0, 1, 2), upbend(Zcomp, Ncomp, 0, 1, 3)
  endif
  do l = 1, gammax
    write(*, '(/" Normalized gamma-ray strength functions and ", "transmission coefficients for l=", i2, /)') l
    write(*, '(" Giant resonance parameters :"/)')
    if (ngr(Zcomp, Ncomp, 1, l) == 2) then
      write(*, '(" sigma0(M", i1, ") =", f8.3, "       sigma0(E", i1, ") =", f8.3, " and ", f8.3, &
 &      "    PR: sigma0(M", i1, ") =", f8.3, "       sigma0(E", i1, ") =", f8.3)') &
 &      l, sgr(Zcomp, Ncomp, 0, l, 1), l, sgr(Zcomp, Ncomp, 1, l, 1), sgr(Zcomp, Ncomp, 1, l, 2), &
        l, tpr(Zcomp, Ncomp, 0, l, 1), l, tpr(Zcomp, Ncomp, 1, l, 1)
    else
      write(*, '(" sigma0(M", i1, ") =", f8.3, "       sigma0(E", i1, ") =", f8.3, &
 &      "    PR: sigma0(M", i1, ") =", f8.3, "       sigma0(E", i1, ") =", f8.3)') &
 &      l, sgr(Zcomp, Ncomp, 0, l, 1), l, sgr(Zcomp, Ncomp, 1, l, 1), &
        l, tpr(Zcomp, Ncomp, 0, l, 1), l, tpr(Zcomp, Ncomp, 1, l, 1)
    endif
    if (ngr(Zcomp, Ncomp, 1, l) == 2) then
      write(*, '("      E(M", i1, ") =", f8.3, "            E(E", i1, ") =", f8.3, " and ", f8.3, &
 &    "    PR:      E(M", i1, ") =", f8.3, "            E(E", i1, ") =", f8.3)') &
 &      l, egr(Zcomp, Ncomp, 0, l, 1), l, egr(Zcomp, Ncomp, 1, l, 1), egr(Zcomp, Ncomp, 1, l, 2), &
        l, epr(Zcomp, Ncomp, 0, l, 1), l, epr(Zcomp, Ncomp, 1, l, 1)
    else
      write(*, '("      E(M", i1, ") =", f8.3, "            E(E", i1, ") =", f8.3, &
 &    "    PR:      E(M", i1, ") =", f8.3, "            E(E", i1, ") =", f8.3)') &
 &      l, egr(Zcomp, Ncomp, 0, l, 1), l, egr(Zcomp, Ncomp, 1, l, 1), &
        l, epr(Zcomp, Ncomp, 0, l, 1), l, epr(Zcomp, Ncomp, 1, l, 1)
    endif
    if (ngr(Zcomp, Ncomp, 1, l) == 2) then
      write(*, '("  gamma(M", i1, ") =", f8.3, "        gamma(E", i1, ") =", f8.3, " and ", f8.3, &
 &      "    PR:  gamma(M", i1, ") =", f8.3, "        gamma(E", i1, ") =", f8.3)') &
 &      l, ggr(Zcomp, Ncomp, 0, l, 1), l, ggr(Zcomp, Ncomp, 1, l, 1), ggr(Zcomp, Ncomp, 1, l, 2), &
        l, gpr(Zcomp, Ncomp, 0, l, 1), l, gpr(Zcomp, Ncomp, 1, l, 1)
    else
      write(*, '("  gamma(M", i1, ") =", f8.3, "        gamma(E", i1, ") =", f8.3, &
 &      "    PR:  gamma(M", i1, ") =", f8.3, "        gamma(E", i1, ") =", f8.3)') &
 &      l, ggr(Zcomp, Ncomp, 0, l, 1), l, ggr(Zcomp, Ncomp, 1, l, 1), &
        l, gpr(Zcomp, Ncomp, 0, l, 1), l, gpr(Zcomp, Ncomp, 1, l, 1)
    endif
    write(*, '("      k(M", i1, ") =", es14.5, "      k(E", i1, ") =", es14.5/)') l, kgr(l), l, kgr(l)
    write(*, '("      E       f(M", i1, ")        f(E", i1, ")", "        T(M", i1, ")        T(E", i1, ")"/)')  l, l, l, l
    do nen = ebegin(0), eend(0)
      e = egrid(nen)
      write(*, '(1x, f8.3, 4es13.5)') e, fstrength(Zcomp, Ncomp, 0., e, 0, l), &
 &      fstrength(Zcomp, Ncomp, 0., e, 1, l), Tjl(0, nen, 0, l), Tjl(0, nen, 1, l)
    enddo
!
! Output on separate files
!
    if (filepsf .and. l == 1) then
      quantity='photon strength function'
      if (gamgam(Zcomp,Ncomp) > 0.) then
        CEgamgam = gamgamth(Zcomp,Ncomp,0) / gamgam(Zcomp,Ncomp)
      else
        CEgamgam = 0.
      endif
      massstring='   '
      write(massstring,'(i3)') A
      finalnuclide=trim(nuc(Z))//adjustl(massstring)
      do irad = 0, 1
        psffile = 'psf000000.E1'
        write(psffile(4:9), '(2i3.3)') Z, A
        if (irad == 0) then
          radtype = 'M'
          MM = strengthM1
          psfmodel = modelM1
        else
          radtype = 'E'
          MM = strength
          psfmodel = modelE1
        endif
        psffile(11:11) = radtype
        open (unit=1, file=psffile, status='replace')
        topline=trim(finalnuclide)//' '//trim(quantity)
        col(1)='E'
        un(1)='MeV'
        col(2)='f('//radtype//'1)'
        un(2)=''
        Ncol=2
        call write_header(topline,source,user,date,oformat)
        call write_residual(Z,A,finalnuclide)
        write(1,'("# parameters:")')
        call write_integer(2,'strength keyword',MM)
        call write_char(2,'PSF model',psfmodel)
        call write_char(2,'radiation type',radtype//'1')
        write(1,'("# observables:")')
        call write_real(2,'experimental Gamma_gamma [eV]',gamgam(Zcomp,Ncomp))
        call write_real(2,'experimental Gamma_gamma unc. [eV]',dgamgam(Zcomp,Ncomp))
        call write_real(2,'theoretical Gamma_gamma [eV]',gamgamth(Zcomp,Ncomp,0))
        call write_real(2,'C/E Gamma_gamma',CEgamgam)
        call write_datablock(quantity,Ncol,Npsf,col,un)
        do nen = 1, Npsf
          e = Epsf(nen)
          write(1,'(2es15.6)') e, fstrength(Zcomp, Ncomp, 0., e, irad, 1)
        enddo
        close (unit=1)
      enddo
    endif
  enddo
!
! **************** Cross sections for inverse channels *****************
!
  quantity='photo absorption cross sections'
  topline=trim(finalnuclide)//' '//trim(quantity)
  crossfile='cross.g'
  col(1)='E'
  un(1)='MeV'
  col(2)='cross section'
  un(2)='mb'
  Ncol = 2
  open (unit=1, file=crossfile, status='replace')
  call write_header(topline,source,user,date,oformat)
  call write_residual(Z,A,finalnuclide)
  write(1,'("# parameters:")')
  call write_char(2,'particle',parname(0))
  Nen =  eend(0) - ebegin(0) + 1
  call write_datablock(quantity,Ncol,Nen,col,un)
  write(*, '(/" Photoabsorption cross sections"/)')
  write(*, '("  E [MeV]   xs [mb]"/)')
  do nen = ebegin(0), eend(0)
    write(1, '(2es15.6)') egrid(nen), xsreac(0, nen)
    write(*, '(2es15.6)') egrid(nen), xsreac(0, nen)
  enddo
  close (unit=1)
  return
end subroutine gammaout
! Copyright A.J. Koning 2021
