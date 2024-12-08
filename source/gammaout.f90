subroutine gammaout(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of gamma-ray strength functions, transmission
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2024-12-08: revised version
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
  character(len=15) :: col(4)    ! header
  character(len=15) :: un(4)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline   ! topline
  character(len=2 ) :: radtype      ! radiation type
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
  real(sgl)         :: xsgamma      ! photo-absorption cross section
  real(sgl)         :: xsgdr        ! photo-absorption cross section from GDR part
  real(sgl)         :: xsqd         ! photo-absorption cross section from QD part
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
  write(*, '(/" ########## GAMMA STRENGTH FUNCTIONS, TRANSMISSION COEFFICIENTS AND CROSS SECTIONS ##########")')
  write(*, '(/" Gamma-ray information for Z=", i3, " N=", i3, " (", i3, a2, ") "/)') Z, N, A, nuc(Z)
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
  if (strength == 11) modelE1 = "D1M-Intra-E1             "
  if (strength == 12) modelE1 = "Shellmodel-E1            "
  if (strengthM1 == 1) modelM1 = "RIPL-1                   "
  if (strengthM1 == 2) modelM1 = "RIPL-2                   "
  if (strengthM1 == 3) modelM1 = "IAEA GSF CRP (2018)      "
  if (strengthM1 == 4) modelM1 = "RIPL-2+Scissors Kawano   "
  if (strengthM1 == 8) modelM1 = "Gogny D1M HFB+QRPA Tables"
  if (strengthM1 == 10) modelM1 = "BSk27+QRPA 2018 Tables   "
  if (strength == 11) modelM1 = "D1M-Intra-M1             "
  if (strength == 12) modelM1 = "Shellmodel-M1            "
  write(*, '(" Gamma-ray strength function model for E1: ", a25)') modelE1
  write(*, '(" Gamma-ray strength function model for M1: ", a25,/)') modelM1
  do l = 1, gammax
!
! Output on separate files
!
    massstring='   '
    write(massstring,'(i3)') A
    finalnuclide=trim(nuc(Z))//adjustl(massstring)
    if (filepsf) then
      if (gamgam(Zcomp,Ncomp) > 0.) then
        CEgamgam = gamgamth(Zcomp,Ncomp,0) / gamgam(Zcomp,Ncomp)
      else
        CEgamgam = 0.
      endif
      do irad = 0, 1
        if (irad == 0) then
          radtype = 'M '
          MM = strengthM1
          psfmodel = modelM1
        else
          radtype = 'E '
          MM = strength
          psfmodel = modelE1
        endif
        write(radtype(2:2),'(i1)') l
        psffile = 'psf000000.'//radtype
        write(psffile(4:9), '(2i3.3)') Z, A
        quantity=radtype//' photon strength function'
        open (unit=1, file=psffile, status='replace')
        topline=trim(finalnuclide)//' '//trim(quantity)
        col(1)='E'
        un(1)='MeV'
        col(2)='f('//radtype//')'
        un(2)=''
        col(3)='T('//radtype//')'
        un(3)=''
        Ncol=3
        call write_header(topline,source,user,date,oformat)
        call write_residual(Z,A,finalnuclide)
        write(1,'("# parameters:")')
        call write_integer(2,'strength keyword',MM)
        call write_char(2,'PSF model',psfmodel)
        if ((irad == 0 .and. strengthM1 <= 4) .or. (irad == 1 .and. (strength <= 2 .or. strength == 5)) .or. l >= 2) then
          call write_real(2,'sigma0 (sgr) [mb]',sgr(Zcomp, Ncomp, irad, l, 1))
          call write_real(2,'E (egr) [MeV]',egr(Zcomp, Ncomp, irad, l, 1))
          call write_real(2,'gamma (ggr) [MeV]',ggr(Zcomp, Ncomp, irad, l, 1))
          if (ngr(Zcomp, Ncomp, irad, l) == 2) then
            call write_real(2,'sigma0_2 (sgr) [mb]',sgr(Zcomp, Ncomp, irad, l, 2))
            call write_real(2,'E_2 (egr) [MeV]',egr(Zcomp, Ncomp, irad, l, 2))
            call write_real(2,'gamma_2 (ggr) [MeV]',ggr(Zcomp, Ncomp, irad, l, 2))
          endif
        endif
        if (flagupbend .and. l == 1) then
          if (irad == 0) then
            call write_real(2,'M1 upbend C',upbendadjust(Zcomp, Ncomp, 0, 1, 1) * upbend(Zcomp, Ncomp, 0, 1, 1))
            call write_real(2,'M1 upbend eta',upbendadjust(Zcomp, Ncomp, 0, 1, 2) * upbend(Zcomp, Ncomp, 0, 1, 2))
            call write_real(2,'M1 upbend F',upbendadjust(Zcomp, Ncomp, 0, 1, 3) * upbend(Zcomp, Ncomp, 0, 1, 3))
          else
            call write_real(2,'E1 upbend C',upbendadjust(Zcomp, Ncomp, 1, 1, 1) * upbend(Zcomp, Ncomp, 1, 1, 1))
            call write_real(2,'E1 upbend eta',upbendadjust(Zcomp, Ncomp, 1, 1, 2) * upbend(Zcomp, Ncomp, 1, 1, 2))
          endif
        endif
        if (l == 1) then
          if (irad == 0) then
            call write_real(2,'pygmy tpr [mb]',tpr(Zcomp, Ncomp, 0, 1, 1))
            call write_real(2,'pygmy epr [MeV]',epr(Zcomp, Ncomp, 0, 1, 1))
            call write_real(2,'pygmy gpr [MeV]',gpr(Zcomp, Ncomp, 0, 1, 1))
          else
            if (tpr(Zcomp, Ncomp, 1, 1, 1) > 0.) then
              call write_real(2,'pygmy tpr [mb]',tpr(Zcomp, Ncomp, 1, 1, 1))
              call write_real(2,'pygmy epr [MeV]',epr(Zcomp, Ncomp, 1, 1, 1))
              call write_real(2,'pygmy gpr [MeV]',gpr(Zcomp, Ncomp, 1, 1, 1))
            endif
          endif
        endif
        call write_real(2,'ftable',ftable(Zcomp, Ncomp, irad, l))
        call write_real(2,'etable',etable(Zcomp, Ncomp, irad, l))
        call write_real(2,'wtable',wtable(Zcomp, Ncomp, irad, l))
        write(1,'("# observables:")')
        call write_real(2,'experimental Gamma_gamma [eV]',gamgam(Zcomp,Ncomp))
        call write_real(2,'experimental Gamma_gamma unc. [eV]',dgamgam(Zcomp,Ncomp))
        call write_real(2,'theoretical Gamma_gamma [eV]',gamgamth(Zcomp,Ncomp,0))
        call write_real(2,'C/E Gamma_gamma',CEgamgam)
        call write_real(2,'average resonance energy [eV]',Eavres * 1.e6)
        call write_real(2,'theoretical S-wave strength function [e-4]',1.e4 * swaveth(Zcomp, Ncomp))
        call write_datablock(quantity,Ncol,Npsf,col,un)
        do nen = 1, Npsf
          e = Epsf(nen)
          write(1,'(3es15.6)') e, fstrength(Zcomp, Ncomp, 0., e, irad, 1),Tjl(0, nen, irad, l)
        enddo
        close (unit=1)
        call write_outfile(psffile,flagoutall)      
      enddo
    endif
  enddo
!
! **************** Cross sections for inverse channels *****************
!
  quantity='photo absorption cross sections'
  topline=trim(finalnuclide)//' '//trim(quantity)
  crossfile='cross.g'
  un='mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='cross section'
  col(3)='GDR'
  col(4)='Quasi-deuteron'
  Ncol = 4
  open (unit=1, file=crossfile, status='replace')
  call write_header(topline,source,user,date,oformat)
  call write_residual(Z,A,finalnuclide)
  write(1,'("# parameters:")')
  call write_char(2,'particle',parname(0))
  Nen =  eend(0) - ebegin(0) + 1
  call write_datablock(quantity,Ncol,Nen,col,un)
  write(*, '(/" Photoabsorption cross sections"/)')
  do nen = ebegin(0), eend(0)
    call gammaxs(Zcomp, Ncomp, egrid(nen), xsgamma, xsgdr, xsqd)
    write(1, '(4es15.6)') egrid(nen), xsgamma, xsgdr, xsqd
  enddo
  close (unit=1)
  call write_outfile(crossfile,flagoutall)      
  return
end subroutine gammaout
! Copyright A.J. Koning 2021
