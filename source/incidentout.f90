subroutine incidentout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Reaction output for incident channel
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
! Variables for OMP
!   Rprime       ! potential scattering radius
!   Sstrength    ! s, p, d, etc - wave strength function
! Variables for basic reaction
!   flagang      ! flag for output of angular distributions
! Variables for numerics
!   nangle       ! number of angles
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   k0           ! index of incident particle
!   Ltarget      ! excited level of target
! Variables for gamma rays
!   fiso         ! correction factor for isospin forbidden transition
!   gammax       ! number of l - values for gamma multipolarity
! Variables for energy grid
!   angle        ! angle in degrees
!   Einc         ! incident energy in MeV
! Variables for incident channel
!   directad     ! direct angular distribution
!   lmaxinc      ! maximal l - value for transm. coeff. for incident channel
!   Tjlinc       ! transm. coeff. as a function of spin and l for inc. channel
!   Tlinc        ! transm. coeff. as a function of l for incident channel
!   xselasinc    ! total elastic cross section (neutrons only) for inc. channel
!   xsreacinc    ! reaction cross section for incident channel
!   xstotinc     ! total cross section (neutrons only) for incident channel
! Constants
!   parname      ! name of particle
!
! *** Declaration of local data
!
  implicit none
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=12) :: Estr
  character(len=15) :: col(5)     ! header
  character(len=15) :: un(5)     ! units
  character(len=80) :: quantity   ! quantity
  integer :: iang              ! running variable for angle
  integer :: l                 ! multipolarity
  integer           :: Ncol        !
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer           :: Nk
  integer :: type              ! particle type
!
! *********** Total cross sections for incident channel ****************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  Estr=''
  write(Estr,'(es12.6)') Einc
  write(*, '(/" Optical model results"/)')
  if (k0 == 1) then
    write(*, '(" Total cross section   :", es11.4, " mb")') xstotinc
  endif
  write(*, '(" Reaction cross section:", es11.4, " mb")') xsreacinc
  if (k0 == 1) then
    write(*, '(" Elastic cross section :", es11.4, " mb")') xselasinc
  endif
!
! For low energy neutrons we give the resonance parameters.
!
! Sstrength: s,p,d,etc-wave strength function
!
  if (k0 == 1 .and. Einc <= 0.1) then
    write(*, '(/" S-wave and P-wave strength functions and potential scattering radius"/)')
    write(*, '("      A      Value"/)')
    write(*, '(" S0:", i4, f8.4, " .e-4")') Atarget, Sstrength(0)*1.e4
    write(*, '(" S1:", i4, f8.4, " .e-4")') Atarget, Sstrength(1)*1.e4
    write(*, '(" R :", i4, f8.4, " fm")') Atarget, Rprime
  endif
  write(*, '(/" Isospin factors to reduce emission "/)')
  do type = 0, 6
    write(*, '(1x, a8, 1x, f8.5)') parname(type), fiso(type)
  enddo
!
! *********** Transmission coefficients for incident channel ***********
!
  if (lmaxinc ==  -1) return
  quantity='transmission coefficients'
  reaction='('//parsym(k0)//',tot)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
  write(*, '(/" Transmission coefficients for incident channel"/)')
  open (unit=1, file='transmission_inc.out', status='replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  call write_reaction(indent,reaction,0.D0,0.D0,0,0)
  call write_real(id2,'E-incident [MeV]',Einc)
  un=''
  col(1)='L'
  if (k0 /= 3 .and. k0 /= 6) then
    col(2)='T(L-1/2,L)'
    col(3)='T(L+1/2,L)'
    col(4)='T(L))'
    Ncol=4
  endif 
  if (k0 == 3) then
    col(2)='T(L-1,L)'
    col(3)='T(L,L)'
    col(4)='T(L+1,L)'
    col(5)='T(L))'
    Ncol=5
  endif
  if (k0 == 6) then
    col(2)='T(L))'
    Ncol=2
  endif 
  if (k0 == 0) then
    Nk=gammax
  else
    Nk=lmaxinc + 1
  endif
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Nk,col,un)
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
  if (k0 == 1 .or. k0 == 2 .or. k0 == 4 .or. k0 == 5) then
    do l = 0, lmaxinc
      write(1, '(i9, 6x, 3es15.6)') l, Tjlinc(-1, l), Tjlinc(1, l), Tlinc(l)
    enddo
  endif
!
! 2. Spin 1 particles: Deuterons
!
  if (k0 == 3) then
    do l = 0, lmaxinc
      write(1, '(i9, 6x, 4es15.6)') l, Tjlinc(-1, l), Tjlinc(0, l), Tjlinc(1, l), Tlinc(l)
    enddo
  endif
!
! 3. Spin 0 particles: Alpha-particles
!
  if (k0 == 6) then
    do l = 0, lmaxinc
      write(1, '(i9, 6x, es15.6)') l, Tjlinc(0, l)
    enddo
  endif
!
! 4. Photons
!
  if (k0 == 0) then
    do l = 1, gammax
      write(1, '(i9, 6x, es15.6)') l, Tjlinc(0, l)
    enddo
  endif
  close (unit = 1)
  call write_outfile('transmission_inc.out',flagoutall)
!
! *********** Shape elastic scattering angular distribution ************
!
  if (flagang .and. k0 > 0) then
    quantity='shape elastic scattering angular distribution'
    reaction='('//parsym(k0)//',el)'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)//' at '//Estr//' MeV'
    open (unit=1, file='shape_el.tot', status='unknown')
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.D0,0.D0,0,0)
    call write_real(id2,'E-incident [MeV]',Einc)
    col(1)='Angle'
    un(1)='deg'
    col(2)='xs'
    un(2)='mb/sr'
    Ncol = 2
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,nangle+1,col,un)
    do iang = 0, nangle
      write(1, '(2es16.5)') angle(iang), directad(k0, Ltarget, iang)
    enddo
    close (unit = 1)
    call write_outfile('shape_el.tot',flagoutall)
  endif
  return
end subroutine incidentout
! Copyright A.J. Koning 2021
