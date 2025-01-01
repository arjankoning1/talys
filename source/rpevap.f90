subroutine rpevap
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Evaporation of residual products
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
  use A0_talys_mod
!   dbl             ! double precision kind
! All global variables
!   idchannel       ! identifier for exclusive channel
!   idnum           ! counter for exclusive channel
!   maxN            ! maximal number of neutrons away from initial compound nucleus
!   maxpfns         ! maximum energy of prompt fission neutrons spectrum
!   maxZ            ! maximal number of protons away from initial compound nucleus
!   multiplicity    ! particle multiplicity
!   nin             ! counter for incident energy
!   Nlast           ! last discrete level
!   nubar           ! average nu
!   nuc             ! symbol of nucleus
!   numen           ! maximum number of outgoing energies
!   numen2          ! maximum number of outgoing energies
!   numia           ! maximum number of alphas in channel description
!   numid           ! maximum number of deuterons in channel description
!   numih           ! maximum number of helions in channel description
!   numin           ! maximum number of neutrons in channel description
!   numip           ! maximum number of protons in channel description
!   numit           ! maximum number of tritons in channel description
!   numneu          ! number of neutrons
!   numpar          ! number of particles
!   nupre           ! pre - neutron emission nu
!   Pdisnu          ! prompt fission neutrons distribution
!   pfns            ! prompt fission neutrons spectrum
!   Rfiseps         ! ratio for limit for fission cross section per nucleus
!   tau             ! lifetime of state in seconds
!   xschannel       ! channel cross section
!   xsfistot        ! total fission cross section
!   xsfpApost       ! post - neutron emission corrected cross section
!   xsfpex          ! excitation energy spectrum per fission fragment
!   xsfptotpost     ! post - neutron emission fission product cross section
!   xsfpZApost      ! post - neutron emission corrected isotopic cross section
!   xsinitpop       ! initial population cross section
!   xspop           ! population cross section
!   xspopex         ! population cross section summed over spin and parity
!   xspopnuc        ! population cross section per nucleus
!   xssumout        ! cross section summed over mechanisms
!   yield           ! yield of produced isotope in MBq / (mA.h)
!   yieldApost      ! post - neutron emission corrected fission yield
!   yieldApre       ! pre - neutron emission fission yield
!   yieldtotpost    ! post - neutron emission fission product yield
!   yieldZApost     ! post - neutron emission corrected isotopic yield
!   yieldZApre      ! pre - neutron emission isotopic yield
!   Ztarget         ! charge number of target nucleus
! Variables for output
!   flagmain        ! flag for main output
! Variables for compound reactions
!   flagcomp        ! flag for compound angular distribution calculation
! Variables for preequilibrium
!   epreeq          ! on - set incident energy for preequilibrium calculation
! Variables for OMP
!   flagompall      ! flag for new optical model calculation for all residual
! Variables for gamma rays
!   flagracap       ! flag for radiative capture model
! Variables for input energies
!   enincmax        ! maximum incident energy
! Variables for main input
!   Arp             ! A of residual product
!   Atarget0        ! mass number of target nucleus
!   Zrp             ! Z of residual product
!   Ztarget0        ! charge number of target nucleus
! Variables for basic reaction
!   flagrpruns      ! flag to designate that run is for residual product
! Variables for multiple emission
!   xspopnuc0       ! population cross section per nucleus
! Variables for nuclides
!   parinclude      ! logical to include outgoing particle
!
! *** Declaration of local data
!
  implicit none
  integer   :: ia                                    ! mass number from abundance table
  integer   :: in                                    ! counter for neutrons
  integer   :: iz                                    ! charge number of residual nucleus
!
! ********************** Loop over fission fragments *******************
!
! do a full TALYS calculation for each residual product
!
  flagrpruns = .true.
!   do type=0,6
!   do nen=0,numen2
!   pfns(type,nen)=0.
!   maxpfns(type,nen)=0.
!   enddo
!   enddo
!   fiseps=Rfiseps*xsfistot
  write(*, '(/" ########## Start of loop over residual products"/)')
  do ia = 1, Atarget0
    do iz = 1, Ztarget0
      in = ia - iz
      if (in < 1 .or. in > numneu) cycle
      if (xspopnuc0(iz, ia) < popeps) cycle
      Arp = ia
      Zrp = iz
      call evaptalys
!
! Add fission product cross sections
!
!   do 110 Zix=0,maxZ
!   do 120 Nix=0,maxN
!   xsfpZApost(iz,in)=xsfpZApost(iz,in)+xspopnuc(Zix,Nix)
!   xsfpApost(ia)=xsfpApost(ia)+xspopnuc(Zix,Nix)
!   do 130 nex=0,Nlast(Zix,Nix,0)
!   if (nex.eq.0.or.tau(Zix,Nix,nex).ne.0.)
!    +            xsfpex(iz,in,nex)=xsfpex(iz,in,nex)+
!    +            xspopex(Zix,Nix,nex)
! 130         continue
! 120       continue
! 110     continue
!
! Add prompt fission particle and gamma production and spectra
!
!   if (xsinitpop.gt.0.) then
!   do 202 npar=1,numin
!   do 204 type=0,6
!   xsexcpart(type,npar)=0.
! 204         continue
!   do 206 iaa=0,numia
!   do 206 ih=0,numih
!   do 206 it=0,numit
!   do 206 id=0,numid
!   do 206 ip=0,numip
!   do 206 inn=0,numin
!   if (inn+ip+id+it+ih+iaa.ne.npar) goto 206
!   ident=100000*inn+10000*ip+1000*id+100*it+10*ih+iaa
!   do 208 idc=0,idnum
!   if (idchannel(idc).eq.ident) then
!   xsc=xschannel(idc)
!   if (inn.gt.0) xsexcpart(1,inn)=xsexcpart(1,inn)+xsc
!   if (ip.gt.0) xsexcpart(2,ip)=xsexcpart(2,ip)+xsc
!   if (id.gt.0) xsexcpart(3,id)=xsexcpart(3,id)+xsc
!   if (it.gt.0) xsexcpart(4,it)=xsexcpart(4,it)+xsc
!   if (ih.gt.0) xsexcpart(5,ih)=xsexcpart(5,ih)+xsc
!   if (iaa.gt.0) xsexcpart(6,iaa)=xsexcpart(6,iaa)+xsc
!   endif
! 208           continue
! 206         continue
! 202       continue
!   yA=yieldApre(ia)
!   yZA=yieldZApre(iz,in)
!   do 210 type=0,6
!   nubar(type)=nubar(type)+yZA*multiplicity(type)
!   if (yA.gt.0.) then
!   nupre(type,ia)=nupre(type,ia)+yZA/yA*multiplicity(type)
!   do 220 npar=1,numin
!   Pdisnu(type,npar)=Pdisnu(type,npar)+
!    +              xsexcpart(type,npar)/xsinitpop
! 220           continue
!   do 230 nen=0,numen2
!   pfns(type,nen)=pfns(type,nen)+yZA*xssumout(type,nen)
! 230           continue
!   endif
! 210       continue
!   endif
    enddo
  enddo
  write(*, '(/" ########## End of loop over residual products"/)')
!
! Average energy and relation to Maxwellian
!
!   do 290 type=0,6
!   sumpfns=0.
!   Esumpfns=0.
!   do 292 nen=1,numen2
!   dE=espec(type,nen)-espec(type,nen-1)
!   sumpfns=sumpfns+pfns(type,nen)*dE
!   Esumpfns=Esumpfns+espec(type,nen)*pfns(type,nen)*dE
!   maxpfns(type,nen)=0.
! 292   continue
!   if (sumpfns.gt.0.) then
!   Eavpfns(type)=Esumpfns/sumpfns
!   else
!   Eavpfns(type)=0.
!   endif
!   summax=0.
!   do 294 nen=1,numen2
!   dE=espec(type,nen)-espec(type,nen-1)
!   E=espec(type,nen)
!   Eav=Eavpfns(type)
!   maxwell=sqrt(E)*exp(-E/Eav)
!   summax=summax+maxwell*dE
! 294   continue
!   do 296 nen=1,numen2
!   E=espec(type,nen)
!   Eav=Eavpfns(type)
!   maxwell=sqrt(E)*exp(-E/Eav)
!   if (maxwell.gt.0..and.sumpfns.gt.0.)
!    +      maxpfns(type,nen)=pfns(type,nen)/maxwell*summax/sumpfns
! 296   continue
! 290 continue
!   sumpost=0.
!   do 310 ia=1,Atarget0
!   sumpost=sumpost+xsfpApost(ia)
! 310 continue
!   sumpost=0.5*sumpost
!   if (sumpost.gt.0.) then
!   xsfptotpost=0.
!   yieldtotpost=0.
!   do 320 iz=1,Ztarget0
!   do 330 ia=iz+1,Atarget0
!   if (xsfpZApost(iz,ia).eq.0.) goto 330
!   in=ia-iz
!   if (in.gt.numneu) goto 330
!   yieldZApost(iz,in)=xsfpZApost(iz,in)/sumpost
!   yieldApost(ia)=yieldApost(ia)+yieldZApost(iz,in)
!   yieldnpost(in)=yieldnpost(in)+yieldZApost(iz,in)
!   xsfptotpost=xsfptotpost+xsfpZApost(iz,in)
!   yieldtotpost=yieldtotpost+yieldZApost(iz,in)
! 330     continue
! 320   continue
!   endif
!
! Reset variables to those of original target.
!
  flagrpruns = .false.
  call talysinput
  flagmain = .false.
  call talysinitial
  flagmain = .true.
  if ( .not. flagompall) call basicxs(0, 0)
  if (parinclude(0)) call gamma(0, 0)
  if (enincmax >= epreeq .or. flagracap) then
    call preeqinit
    call excitoninit
  endif
  if (flagracap) call racapinit
  if (flagcomp) call compoundinit
!
! Output
!
!   call massdisout
!   call nubarout
!   call nudisout
!   if (flagspec) call pfnsout
  return
end subroutine rpevap
! Copyright A.J. Koning 2021
