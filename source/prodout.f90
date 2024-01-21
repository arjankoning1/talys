subroutine prodout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of isotope production
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
! All global variables
!   numtime         ! number of time points
! Variables for medical isotope production
!   Area            ! target area in cm^2
!   Eback           ! lower end of energy range in MeV for isotope
!   Ebeam           ! incident energy in MeV for isotope production
!   Ibeam           ! beam current in mA for isotope production
!   radiounit       ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi, Ci or kCi
!   rhotarget       ! target material density
!   Tcool           ! cooling time per unit cooling time unit (y, d, h, m, s)
!   Tirrad          ! irradiation time per unit irradiation time unit (y, d, h, m, s)
!   yieldunit       ! unit for isotope yield: num (number), mug, mg, g, or kg
! Variables for numerics
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxZ            ! maximal number of protons away from initial compound nucleus
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ninit           ! neutron number of initial compound nucleus
!   Zinit           ! charge number of initial compound nucleus
!   Ztarget         ! charge number of target nucleus
! Constants
!   iso             ! counter for isotope
!   isochar         ! symbol of isomer
!   natstring       ! string extension for file names
!   nuc             ! symbol of nucleus
!   parsym          ! symbol of particle
! Variables for decay data
!   lambda          ! decay rate per isotope
!   Td              ! half life per time unit
!   Thalf           ! half life of nuclide in sec.
! Variables for levels
!   Lisomer         ! level number of isomer
!   Nisomer         ! number of isomers for this nuclide
! Variables for isotope production
!   activity        ! activity of produced isotope in MBq
!   heat            ! produced heat
!   Mtar            ! active target mass
!   Niso            ! number of isotopes produced after irradiation
!   Nisorel         ! fraction of number of produced isotopes per element
!   Ntar0           ! number of original target atoms
!   prate           ! production rate per isotope
!   projnum         ! number of incident particles [s^ - 1]
!   targetdx        ! effective thickness of target
!   Tgrid           ! time
!   Tmax            ! irradiation time with maximal yield
!   Tmaxactivity    ! time of maximum activity of produced isotope in MBq
!   Tp              ! irradiation time with maximal yield per time unit
!   Vtar            ! active target volume
!   yield           ! yield of produced isotope in MBq / (mA.h)
! Variables for existence libraries
!   Yexist          ! flag for existence of yield
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: rstr        ! string
  character(len=3)  :: ystr        ! string
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=15) :: Yfile       ! file with production yields
  character(len=35) :: halflife    ! half life
  character(len=35) :: maxprod     ! maximum production
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(5)    ! header
  character(len=15) :: un(5)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: string    ! 
  character(len=132) :: topline    ! topline
  integer           :: Ncol       ! counter
  integer           :: A           ! mass number of target nucleus
  integer           :: is          ! isotope counter: -1=total, 0=ground state 1=isomer
  integer           :: it          ! counter for tritons
  integer           :: k           ! designator for particle
  integer           :: N           ! neutron number of residual nucleus
  integer           :: Nix         ! neutron number index for residual nucleus
  integer           :: Z           ! charge number of target nucleus
  integer           :: Zix         ! charge number index for residual nucleus
!
! ************************* Main output ********************************
!
!   mCi, Ci or kCi
!   mug (micro-gram), mg, g, or kg
!
  rstr = 'MBq'
  ystr = '   '
  if (radiounit == 'bq') rstr = ' Bq'
  if (radiounit == 'kbq') rstr = 'KBq'
  if (radiounit == 'gbq') rstr = 'GBq'
  if (radiounit == 'ci') rstr = ' Ci'
  if (radiounit == 'kci') rstr = 'KCi'
  if (radiounit == 'mci') rstr = 'mCi'
  if (yieldunit == 'g') ystr = '  g'
  if (yieldunit == 'mug') ystr = 'mug'
  if (yieldunit == 'mg') ystr = ' mg'
  if (yieldunit == 'kg') ystr = ' kg'
  write(*, '(/" Summary of isotope production for ", a1, " + ", a/)') ptype0, trim(targetnuclide)
  string=''
  write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tirrad(k), k = 1, 5)
! call write_char(6,'Maximum irradiation time',string)
  write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tcool(k), k = 1, 5)
! call write_char(6,'Cooling time',string)
! write(*, '(" Maximal irradiation time    : ", i3, " years ", i3, &
!&  " days ", i3, " hours ", i3, " minutes ", i3, " seconds ")') (Tirrad(i), i = 1, 5)
! write(*, '(" Cooling time                : ", i3, " years ", i3, &
!&  " days ", i3, " hours ", i3, " minutes ", i3, " seconds ")') (Tcool(i), i = 1, 5)
! call write_real(6,'Beam current [mA]',Ibeam)
! call write_real(6,'E-Beam [MeV]',Ebeam)
! call write_real(6,'E-Back [MeV]',Eback)
! call write_real(6,'Target material density [g/cm^3]',rhotarget)
! call write_real(6,'Target area [cm^2]',Area)
! call write_real(6,'Effective target thickness [cm]',targetdx)
! call write_real(6,'Effective target volume [cm^3]',Vtar)
! call write_real(6,'Effective target mass [g]',Mtar)
! call write_real(6,'Number of target atoms',Ntar0)
! call write_real(6,'Number of incident particles',projnum)
! call write_real(6,'Produced heat in target [kW]',heat)
! call write_real(6,'Total production rate [sec^-1]',prate(-1, -1, -1))
! write(*, '(" Energy range                : ", f8.3, " --> ", f8.3, " MeV")') Ebeam, Eback
! write(*, '(" Beam current                : ", f12.3, " mA")') Ibeam
! write(*, '(" Target material density     : ", f12.3, " g/cm^3")') rhotarget
! write(*, '(" Target area                 : ", f12.3, " cm^2")') Area
! write(*, '(" Effective target thickness  : ", f12.3, " cm")') targetdx
! write(*, '(" Effective target volume     : ", f12.3, " cm^3")') Vtar
! write(*, '(" Effective target mass       : ", f12.3, " g   ")') Mtar
! write(*, '(" Number of target atoms      : ", es12.5)') Ntar0
! write(*, '(" Number of incident particles: ", es12.5, " s^-1")') projnum
! write(*, '(" Produced heat in target     : ", f12.3, " kW")') heat
  write(*, '(/" (Maximum) production and decay rates per isotope"/)')
! write(*, '(" Production rate for all isotopes: ", es12.5, " [s^-1]"/)') prate(-1, -1, -1)
  write(*, '("#  Nuc     Production rate Decay rate     Activity       #isotopes   Yield          Isotopic", &
 &  "                   Half-life               Time of maximum production")')
  write(*, '("#             [s^-1]         [s^-1]         [", a3, "]          [", a3, "]        [", a3, "/mAh]      fraction")') &
 & rstr, ystr, rstr
  do Zix = 0, maxZ
    Z = Zinit - Zix
    do Nix = 0, maxN
      N = Ninit - Nix
      A = Z + N
      do is = - 1, Nisomer(Zix, Nix)
        if ( .not. Yexist(Zix, Nix, is)) cycle
        it = Tmaxactivity(Zix, Nix, is)
        halflife = '                                   '
        if (Thalf(Zix, Nix, is) > 1.e17) then
          write(halflife, '(a13)') '     stable  '
        else
          write(halflife, '(i8, " y ", i3, " d ", i3, " h ", i3, " m ", i3, " s ")') (Td(Zix, Nix, is, k), k = 1, 5)
        endif
        maxprod = '                                   '
        if (Tmax(Zix, Nix, is) > 1.e17) then
          write(maxprod, '(a13)') '     infinite'
        else
          write(maxprod, '(i8, " y ", i3, " d ", i3, " h ", i3, " m ", i3, " s ")') (Tp(Zix, Nix, is, k), k = 1, 5)
        endif
        write(*, '(1x, a2, i4, 1x, a1, 5es15.6, f8.5, 2a35)') &
 &        nuc(Z), A, isochar(is), prate(Zix, Nix, is), lambda(Zix, Nix, is), &
 &        activity(Zix, Nix, is, it), Niso(Zix, Nix, is, it), &
          yield(Zix, Nix, is, it), Nisorel(Zix, Nix, is, it), halflife, maxprod
      enddo
    enddo
  enddo
!
! Output to files per residual product
!
  reaction='('//parsym(k0)//',x)'
  quantity='Isotope production'
  write(*,'(/,"Medical isotope production output files per nuclide:",/)')
  do Zix = 0, maxZ
    Z = Zinit - Zix
    do Nix = 0, maxN
      N = Ninit - Nix
      A = Z + N
      do is = - 1, Nisomer(Zix, Nix)
        if ( .not. Yexist(Zix, Nix, is)) cycle
        Yfile = 'Y000000.tot'//natstring(iso)
        write(Yfile(2:7), '(2i3.3)') Z, A
        if (is >= 0) Yfile(9:11) = 'L00'
        if (is >= 1) write(Yfile(10:11), '(i2.2)') Lisomer(Zix, Nix, is)
        massstring='   '
        write(massstring,'(i3)') A
        finalnuclide=trim(nuc(Z))//trim(adjustl(massstring))//isochar(is)
!       call write_char(6,quantity,Yfile)
        open (unit = 1, file = Yfile, status = 'replace')
        topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
        un = ''
        col(1)='Time'
        un(1)='h'
        col(2)='Activity'
        un(2)=rstr
        col(3)='Isotopes'
        un(3)=ystr
        col(4)='Yield'
        un(4)=rstr
        col(5)='Isotopic fract.'
        Ncol=5
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.d0,0.d0,0,0)
        call write_residual(Z,A,finalnuclide)
        write(1,'("# parameters:")')
        call write_real(2,'Beam current [mA]',Ibeam)
        call write_real(2,'E-Beam [MeV]',Ebeam)
        call write_real(2,'E-Back [MeV]',Eback)
        string='Initial production rate [s^-1]'
        call write_real(2,string,prate(Zix, Nix, is))
        string='Decay rate [s^-1]'
        call write_real(2,string,lambda(Zix, Nix, is))
        string='Initial production yield ['//rstr//'/mAh]'
        call write_real(2,string,yield(Zix, Nix, is, 1))
        string='Total activity at EOI ['//rstr//']'
        call write_real(2,string,activity(Zix, Nix, is, Ntime))
        string=''
        write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tirrad(k), k = 1, 5)
        call write_char(2,'Irradiation time',string)
        write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tcool(k), k = 1, 5)
        call write_char(2,'Cooling time',string)
        if (Thalf(Zix, Nix, is) > 1.e17) then
          string='stable'
        else
          write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') &
 & (Td(Zix, Nix, is, k), k = 1, 5)
        endif
        call write_char(2,'Half-life',string)
        if (Tmax(Zix, Nix, is) > 1.e17) then
          string='infinity'
        else
          write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")')  &
 & (Tp(Zix, Nix, is, k), k = 1, 5)
        endif
        call write_char(2,'Maximum production at',string)
        call write_datablock(quantity,Ncol,numtime,col,un)
        do it = 1, numtime
          write(1, '(5es15.6)') Tgrid(it), activity(Zix, Nix, is, it), Niso(Zix, Nix, is, it), &
 &          yield(Zix, Nix, is, it), Nisorel(Zix, Nix, is, it)
        enddo
        close (unit = 1)
      enddo
    enddo
  enddo
  return
end subroutine prodout
! Copyright A.J. Koning 2021
