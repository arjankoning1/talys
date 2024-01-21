      MODULE om_retrieve

!To print change "ccalone" to "alone" in the DEFINE statement below
!MS$DEFINE ccalone
!
!======================================================================
!                       MODULE om_retrieve
!======================================================================
!
! GOAL: Code to retrieve optical model potentials from the RIPL optical
!       model potential (omp) library and to format them for input into
!       the SCAT2000, ECIS06 and OPTMAN computer codes.
!
!----------------------------------------------------------------------
!
! USE:
!
! To call the module:
!
!    > use om_retrieve
!
! Used I/O files:
! ki = 171 ko = 15 standard output: 6
! temporalily used files: 13,34,35,36,37,38,44,45
!
! Note: File 6 remain open !
!
c------------------------ PUBLIC INTERFACE -----------------------------
C PUBLIC :: retrieve
C PUBLIC  Number_Energies, Energies, Emin, Emax
C PUBLIC  Ztarget, Atarget, RIPL_Index, Calc_Type, Def_Index, Iflag
!
!****************************************************************************************************
! WARNING; This module does not produce correct OMP tables for soft-rotor model (imodel=3) potentials                                    :
!****************************************************************************************************
!
!  RETURN FLAG VALUES:
!
!  Iflag = 0  normal return: Both OMP and coupling scheme tables produced
!  Iflag > 0  ERROR: OMP table not printed
!  Iflag < 0  WARNING: OMP table produced but not the coupling scheme table
!
!  WARNINGS:
!  Iflag = -1  No coupling scheme available, but ECIS can be used
!  Iflag = -2  Inputted target Z outside range of RIPL potential
!  Iflag = -3  Inputted target A outside range of RIPL potential
!
!  ERRORS:
!  Iflag = +1  No coupling scheme available, OPTMAN is needed
!  Iflag = +2  Structure model not implemented yet in RIPL library
!  Iflag = +3  Coupled-channel model not implemented in SCAT2000
!  Iflag = +4  Only OPTMAN code can be used for soft rotor potentials
!  Iflag = +5  Producing a dispersive input from a non-dispersive potential, change modtyp to 2.
!  Iflag = +6  Gaussian form factor is not allowed in ECIS
!  Iflag = +7  Non-local potentials not allowed in ECIS
!  Iflag = +8  Wrong type of calculations requested (modtyp>6)
!  Iflag = +9  End of RIPL file: Usually it means that the requested potential does not exist
!  Iflag = +10 Potential dependence not implemented in OPTMAN interfase
!  Iflag = +11 OPTMAN cant not calculate CC vibrational model
!  Iflag = +12 Dispersive parameters Nv<>Ns in OPTMAN
!  Iflag = +13 VIB and DWBA options valid for even-even targets
!  Iflag = +14 Real and Imag radii diff. in OPTMAN
!  Iflag = +15 Real and Imag diffuness diff. in OPTMAN
!---------------------------------------------------------------
!
! COMMENTS:
!
! Originally coded within RIPL-2 project by
! P.G.Young, Group T-2, Los Alamos National Laboratory
! Mail Stop B243, email address pgy@lanl.gov
!
! R. Capote Noy, IAEA Nuclear Data Section
! email address R.CapoteNoy@iaea.org or rcapotenoy@yahoo.com
!
! P. Talou, Group T-2, Los Alamos National Laboratory
! email address talou@lanl.gov
!---------------------------------------------------------------
!
! RELEASE(S):
!! Date         Release #          Author(s)           Comments
! ----         ---------          ---------           --------
! 23/08/2005      0.1         A.Delgado, R.Capote  Original version
! 07/11/2005      0.2             R.Capote         OPTMAN added
! 31/07/2006      0.3             R.Capote         Format extended
!                                                  Dispersive ECIS input added
! 28/08/2008      0.4         R.Capote, P.Talou    RIPL-3 version
!
! 30/09/2008      0.5         R.Capote, P.Talou    RIPL-3 version (validated)
!
! 14/10/2015      0.6         R.Capote             Multi-band CC added
!
! 10/12/2019      0.7         R.Capote             Printing and finetuning with AK for TALYS
!
! 29/05/2020      0.8         R.Capote             NDIM6 dimension increased for multi-band potentials
!                                                  (done before in EMPIRE). Cleanup of I/O.
!                                                  Adapted to be used with TALYS for imodel=0,1,2,4
!                                                  (all structure models but the soft rotator imodel=3)
!---------------------------------------------------------------
!
! AUTHOR(S) INFORMATIONS:
!
c     Version date: August 2008, RCN & PT
c
c     Some details and protections added based on P.Talou's validation.
c     Vibrational model corrected to allow calculation at energies below
c     excited levels' energy (i.e. with closed channels). It is noted that
c     a true CC calculation with more than two levels coupled within the
c     vibrational model is not coded (RCN).
c
c------------------------------------------------------------------------
c     Version date: November 2007, RCN
c
c     Volume integrals are calculated if modtyp <= 2.
c     ECIS06 input implemented
c
c------------------------------------------------------------------------
c     Version date: September 23 2004, RCN
c
c     1. Format expanded to allow for:
c        -Energy dependent radius up to E**2    (r(i,j,12))
c        -Dispersive-like energy dependent radius ( r(i,j,13)>0 )
c        -Soukhovitski et al energy dependence of the real potential
c         See J.Phys.G. Nucl.Part.Phys. 30(2004) 905-920)
c        -Buck and Perey energy dependence of the local real potential
c         (As coded and employed by B. Morillon and P. Romain,
C         See Phys.Rev.C70(2004) 014601)
c
c     2. Analytical dispersive integrals are included
c         See Quesada JM et al, Comp. Phys. Comm. 153(2003) 97
c                               Phys. Rev. C67(2003) 067601
c
c     3. General numerical solution of the dispersive integral is
c        kept to deal with special Ws(E) and Wd(E) form factors.
c         See Capote et al, J.Phys.G. Nucl.Part.Phys. 27(2001) B15-B19
c
c        Numerical integration is used for Nagadi et al dispersive OMP
c         See Nagadi et al, Phys. Rev. C68 (2003) 044610
c
c     4. Relativistic kinematics is used to calculate the kinematical
c        conversion factor xkine and reduced mass amu (for irel>0)
c        Kinematical conversion factor for discrete levels is assumed
c        non-relativistic.
c
c     5. The isimple=1 option to send energy block to SCAT code was
c        eliminated (it could be reproduced by sending energies one
c        by one. The speed gain is negligible today)
c
c------------------------------------------------------------------------
c
c     P.G.Young, Group T-2, Los Alamos National Laboratory
c     Mail Stop B243, email address pgy@lanl.gov
c
c     R. Capote Noy, IAEA Nuclear Data Section
c     email address R.CapoteNoy@iaea.org or rcapotenoy@yahoo.com
c
c     Code compiled using :
c         f77 om-retrieve.f -C -O3 -o omretrieve
c     That is, -r8 should not be turned on in order for the dispersive
c     corrections to be computed accurately.
c
c     INPUT FILES
c       omp-parameter-u.dat = RIPL3 Optical Model Potential (OMP) file
c       gs-mass-sp.dat = ground-state mass, spin-parity file
c       ominput.inp = input instructions, defined below
c
c     OUTPUT FILES
c       sc2.inp   = input file for the scat2000 code (modtyp=1)
c       ecis.inp  = standard input file for the ECIS code (modtyp=2,5,6)
c       ecistc.inp= alternate input file for the ECIS code.
c                   (modtyp=2,3). Useful for generating transmission
c                   transmission coefficients. Obsolete.
c       ecisvib   = input file for the ECIS code with vibrational
c                   model activated (modtyp=2).
c       ecisdw.inp= input file for the ECIS code with DWBA
c                   model activated (modtyp=3).
c       optman.inp= standard input file for the OPTMAN code (modtyp=4)
c       massinfo.out = descriptive remarks about the ground-state mass
c                      data used in the gs-mass-sp.dat file.
c       omp-table.dat = a table of numerical values of each potential made
c                   at each incident energy when modtyp=2 (ECIS input).
c       ccomp-lev.dat = a table of the coupling scheme when imodel>0
c
c     INPUT DETAILS (ominput.inp file)
c
c     read(5,*) ne
c
c       ne = number of incident energies
c          = 0 to use a built-in array of incident energies
c
c     if(ne.gt.0) read(5,*) (en(n),n=1,ne)
c
c       en(n) = incident energies in MeV in the laboratory system
c
c     read(5,*,end=990) iztar,iatar,irefget,modtyp
c
c       abs(iztar) = Z of target nucleus
c                    Set iztar negative to provide integer projectile
c                    mass in input decks.
c       abs(iatar) = A of target nucleus
c                    Set iatar negative to provide integer target
c                    mass in input decks.
c
c       irefget = reference number of optical model potential to be
c                 retrieved from the RIPL library
c         Set irefget = negative izaproj (projectile) to retrieve all
c                       spherical potentials in the RIPL library for
c                       this izaproj and the inputted izatar (target).
c
c       iabs(modtyp)= 1 to generate SCAT2 input file (sc2.inp)
c       iabs(modtyp)= 2 to generate ECIS input file (ecis.inp)
c                       For dispersive potentials it is assuming that the Im(geom)=Real(geom).
c                       Therefore for dispersive potentials it is better to use modtyp = 5
c       iabs(modtyp)= 3 to generate ECIS DWBA input files
c                       (ecisdw.inp), using structure information from
c                       an external file (deform.dat), which may be the
c                       om-deformations.dat*file or a user-provided
c                       deform.dat file.
c       iabs(modtyp)= 4 to generate OPTMAN input file (optman.inp)
c       iabs(modtyp)= 5 to generate ECIS input file for dispersive
c                       potentials with externally calculated corrections
c                       (recommended option for dispersive potentials)
c       iabs(modtyp)= 6 to generate ECIS input file for dispersive
c                       potentials with ECIS internally calculated dispersive
c                       corrections. Avoid using it for charged particle potentials.
c
c     if(iabs(modtyp).eq.3)read(5,*)kdef
c
c       Set modtyp negative to force imaginary volume and surface
c       potentials to be zero or positive.  If iabs(modtyp)=1, then
c       all SCAT2000 inputs are single energy, i.e., the compact
c       energy representation of SCAT2000 is bypassed.
c
c
c         kdef=1 use JENDL-3.2 betas from om-deformations.dat file
c             =2 use ENSDF(Q)
c             =3 use ENSDF(BE2) and ENSDF(BE3)
c             =4 use Raman and Spear for BE2,BE3
c             =5 use deform.dat file provided by user
c
c       See om-deformations.readme file for more information on the
c       kdef options and for the format of deform.dat for using kdef=5.
c       For kdef=1-4, the om-deformations.dat file must be copied to
c       deform.dat.  Note that where multiple entries occur for a single
c       state in om-deformations.dat, the code makes an ECIS input
c       for each entry.
c
      implicit real*8 (A-H,O-Z)

      PRIVATE

c------------------------ PUBLIC INTERFACE -----------------------------

      PUBLIC :: retrieve

      PUBLIC  Number_Energies, Energies, Emin, Emax
      PUBLIC  Ztarget, Atarget, RIPL_Index, Calc_Type, Def_Index, Iflag

c--------------------------- DECLARATIONS -----------------------------

      integer Number_Energies, Iflag
      integer Ztarget, Atarget, RIPL_Index, Calc_Type, Def_Index
      real*8  Energies (500), Emin, Emax
C
C     RCN, 08/2004, to handle new extension to the OMP RIPL format
C     PARAMETER(NDIM1 = 10, NDIM2 = 13, NDIM3 = 24, NDIM4 = 30,
C    &          NDIM5 = 10, NDIM6 = 10, NDIM7 = 120)
C     RCN, 08/2008, to handle new extension to the OMP RIPL-3 format
C     PARAMETER(NDIM1 = 10, NDIM2 = 13, NDIM3 = 25, NDIM4 = 40,
C    &          NDIM5 = 10, NDIM6 = 10, NDIM7 = 120, NDIM8 = 50)
C     RCN, 09/2014, to increase NDIM6 for rigid-rotor potentials
C                   (multi-level coupling)
C     PARAMETER(NDIM1 = 10, NDIM2 = 13, NDIM3 = 25, NDIM4 = 40,
C    &          NDIM5 = 10, NDIM6 = 30, NDIM7 = 120, NDIM8 = 50)
      parameter (ndim1=10, ndim2=13, ndim3=25, ndim4=40, ndim5=10,
     + ndim6=30,ndim7=120,ndim8=50)

      character*1 author(80),refer(80),summary(320)
      integer iref
      integer izmin,izmax,iamin,iamax,imodel
      integer jrange(6)
      real*8 epot(6,ndim1),rco(6,ndim1,ndim2),aco(6,ndim1,ndim2),
     +  pot(6,ndim1,ndim3)
      integer ncoll(ndim4),nvib(ndim4),nisotop,iz(ndim4),ia(ndim4),
     +  lmaxOM(ndim4)
      real*8 bandk(ndim4),def(ndim4,ndim5)
      integer idef(ndim4),izproj,iaproj
      real*8 exv(ndim7,ndim4)
      integer iparv(ndim7,ndim4),irel,nph(ndim7,ndim4)
      real*8 defv(ndim7,ndim4),thetm(ndim7,ndim4),ex(ndim6,ndim4),
     + spin(ndim6,ndim4),spinv(ndim7,ndim4),ecoul(ndim1),rcoul(ndim1),
     + rcoul0(ndim1),beta(ndim1),rcoul1(ndim1),rcoul2(ndim1),
     + acoul(ndim1),rcoul3(ndim1),defr(ndim7,ndim4)
      integer ipar(ndim6,ndim4),jcoul

      integer ne, izatar, irefget, modtyp, iztar, iatar, izaproj, irndp,
     +  irndt, nonegw

      integer ko,ieof
      integer IDRs

      character*1 opname(14),tarpar,projpar
      character*16 relz
      character*11 modelz
      character*27 dispz

c     physical constants
      real*8 amu, amu0c2, hbarc, cso

      integer ipsc,lmaxec
      real*8 tarmas,projmas,projspi,ztar,atar,zproj,tarspi
      real*8 coul,eta,efermi,xratio,el,encoul

      real*8 v,rv,av,w,rw,aw,vd,rvd,avd,wd,rwd,awd,
     + vso,rvso,avso,wso,rwso,awso,rc,ac

      character*1 tgtp

      character*2 nuc(124)
      character*8 parname(7)

      character*50 becis1,becis2,ecis1,ecis2
      integer ntar,ncol
      integer koptom,kecis,jchk

      integer kdef,ndisx,ndis
      real*8 tgts,edis(ndim8),bdis(ndim8),sdis(ndim8)
      integer ipdis(ndim8),langmo(ndim8)
c
      integer ki,idr
      real*8 en(500)
c
      integer nesc,isimple,nsc2
      real*8 ensc(300)

      real*8 SR_hw(ndim4),SR_amb0(ndim4),SR_amg0(ndim4),
     + SR_gam0(ndim4),SR_bet0(ndim4),SR_bet4(ndim4),
     + SR_bb42(ndim4),SR_gamg(ndim4),SR_delg(ndim4),
     + SR_bet3(ndim4),SR_et0(ndim4),SR_amu0(ndim4),
     + SR_hw0(ndim4),SR_bb32(ndim4),SR_gamde(ndim4),
     + SR_dpar(ndim4),SR_gshape(ndim4)
      integer SR_ntu(ndim6,ndim4),SR_nnb(ndim6,ndim4),
     + SR_nng(ndim6,ndim4),SR_nno(ndim6,ndim4)

      data ki,idr/171,0/
      data (nuc(i),i=1,124) /
     +  'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     +  'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     +  'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     +  'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     +  'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     +  'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     +  'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +  'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     +  'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     +  'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     +  'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
     +  'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','B9','C0',
     +  'C1','C2','C3','C4'/
      data (parname(i),i=1,7) /'neutron ','proton  ','deuteron',
     +  'triton  ','he-3    ','alpha   ','gamma   '/
c------------------
      CONTAINS
c------------------


c-----------------------------------------------------------------------------
      SUBROUTINE retrieve
c-----------------------------------------------------------------------------
      logical optman,scat,ecis

      character*3 ifin
      data isuit,ifin/0,'FIN'/
        character*66 Diagnostic(19)
        data Diagnostic/
     +'WARNING:Inputted target A outside range of RIPL potential',       ! Iflag = -3
     +'WARNING:Inputted target Z outside range of RIPL potential',   ! Iflag = -2
     +'WARNING:No coupling scheme available, but ECIS can be used',           ! Iflag = -1
     +'RIPL OMP retrieval succeeded',                                                            ! Iflag =  0
     +'ERROR:No coupling scheme available, OPTMAN is needed',                       ! Iflag = +1
     +'ERROR:Structure model not implemented yet in RIPL library',            ! Iflag = +2
     +'ERROR:Coupled-channel model not implemented in SCAT2000',                   ! Iflag = +3
     +'ERROR:OPTMAN code should be used for soft rotor potential',           ! Iflag = +4
     +'ERROR:Non-dispersive potential, set modtyp to 2',                               ! Iflag = +5
     +'ERROR:Gaussian form factor not allowed in ECIS',                               ! Iflag = +6
     +'ERROR:Non-local potentials not allowed in ECIS',                               ! Iflag = +7
     +'ERROR:Wrong type of calculations requested (modtyp>6)',                   ! Iflag = +8
     +'ERROR:EOF: Requested potential (key) not found',                               ! Iflag = +9
     +'ERROR:Potential dependence not implemented in OPTMAN interfase',  ! Iflag = +10
     +'ERROR:OPTMAN cant not calculate CC vibrational model',                       ! Iflag = +11
     +'ERROR:Dispersive parameters Nv<>Ns in OPTMAN',                                       ! Iflag = +12
     +'ERROR:VIB and DWBA options valid for even-even targets',                   ! Iflag = +13
     +'ERROR:Diff. real and imag radii can not be used in OPTMAN',       ! Iflag = +14
     +'ERROR:Diff. real and imag diffunesses can not be used in OPTMAN'/ ! Iflag = +15

      data nee/21/
      data optman/.false./,scat/.false./,ecis/.false./

      open(unit=ki,file='om-parameter-u.dat')
      open(unit=35,file='omp-table.dat')
!MS$IF DEFINED (alone)
        open(unit=37,file='vint-table.dat')
!MS$ENDIF
      Iflag = 0 ! normal return
      IDRs  = 0
c----------------------------------------------------------------
c     Transfer input information
c
      ne = Number_Energies
      if (ne.eq.0) nee = 0
      do n=1,ne
       en(n)= Energies(n)
      enddo
      iztar = Ztarget
      iatar = Atarget
      irefget = RIPL_Index
      modtyp = Calc_Type
      nonegw=0
      if(modtyp.lt.0) nonegw=1
      modtyp=iabs(modtyp)
        kdef = 0
      if (modtyp.eq.3) kdef = Def_Index
c----------------------------------------------------------------
      nsc2=0
      nread = 0
      idr = 0

      irndp=0
      if(iztar.lt.0)irndp=1
      iztar=iabs(iztar)
      irndt=0
      if(iatar.lt.0)irndt=1
      iatar=iabs(iatar)
      izatar=1000*iztar+iatar
c
c     output unit
c
      ko=15
      ecis = .true. ! default
      if(modtyp.eq.1) then
        scat = .true.
        ecis = .false.
      endif

      if(modtyp.eq.4) then
        optman = .true.
        ecis = .false.
      endif
      if(modtyp.eq.1) open(unit=ko,file='sc2.inp')
      if(modtyp.eq.2) open(unit=ko,file='ecis.inp')
      if(modtyp.eq.3) open(unit=ko,file='ecisdw.inp')
      if(modtyp.eq.4) open(unit=ko,file='optman.inp')
      if(modtyp.ge.5) then
        open(unit=ko,file='ecis.inp')
C*****************************
C       Skipping CMS calculation
C       open(unit=16,file='ecistc.inp')
C*****************************
      endif
C*****************************
C     Skipping CMS calculation
C     if(modtyp.eq.2 .or. modtyp.eq.3) open(unit=16,file='ecistc.inp')
C*****************************
C     if(modtyp.eq.4) open(unit=16,file='optmantc.inp')
c
c     Get optical model parameters
c
  98  rewind ki
 100  call omin30
      if(ieof.eq.1) then
        close(ki)
        close(ko)
!MS$IF DEFINED (alone)
          close(44)
!MS$ENDIF
        write(6,102) Iflag,Diagnostic(4+Iflag)
 102    format(1x,90(1h*)/1x,'OMP retrieval flag=',i3,2x,A66/1x,90(1h*))
        return
      endif
!MS$IF DEFINED (alone)
      if( idr.ne.0 .and. modtyp.le.2) then
        open(unit=36,file='DOM-table.dat')
        open(unit=34,file='DOMtoW.dat')
      endif
      if(izproj.GT.0.and. modtyp.le.2)
     &    open(unit=38,file='coul-vint-table.dat')
!MS$ENDIF
c
c     Set energy grid up to emax if nee=0
c
      if(nee.eq.0) call egrid(ne,emin,emax,en)
      call optname
c
c     Check for desired reference and Z,A range
 108  if(irefget.lt.0.and.imodel.lt.1) go to 110
      if(irefget.eq.iref) go to 112
      go to 100
 110  izapr=iabs(irefget)
      izpr=izapr/1000
      iapr=mod(izapr,1000)
      if(izpr.ne.izproj.or.iapr.ne.iaproj) go to 100
 112  if(izmin.eq.0.and.izmax.eq.0) go to 114
      if(iztar.lt.izmin.or.iztar.gt.izmax)  then
        write(6,
     +'(" WARNING: Inputted target Z outside range of RIPL potential")')
        Iflag = -2
        endif
 114  if(iamin.eq.0.and.iamax.eq.0) go to 120
      if(iatar.lt.iamin.or.iatar.gt.iamax) then
        write(6,
     +'(" WARNING: Inputted target A outside range of RIPL potential")')
        Iflag = -3
        endif
c
c     Check if input file can be made for this potential
c
 120  if(imodel.le.4)go to 122
      write(6,
     + '(" ERROR: imodel=",i2," in RIPL library not implemented yet.")')
      Iflag = 2
      go to 1020
c
 122  if(modtyp.ne.1)go to 126
      if(imodel.eq.0)go to 126
      write(6,
     +'(" ERROR: Coupled-channel model not implemented for SCAT2000")')
      Iflag = 3
      go to 1020

 126  if(imodel.eq.2.and.modtyp.eq.3)write(6,'(" You are making an ECIS0
     +3 DWBA input from a vibrational potential.")')
c
      if(imodel.eq.1.and.modtyp.eq.3)write(6,'(" You are making an ECIS0
     +3 DWBA input from a rotational potential.")')

C     if(imodel.eq.4) write(6,'(" You are making an input from a rigid-s
C    +oft rotor potential.")')

      if(imodel.eq.3.and.modtyp.ne.4) then
       write(6,
     + '(" ERROR: You are making an input from a soft-rotor potential.")
     +')
       write(6,
     + '(" ERROR: Only OPTMAN code can be used, set modtyp = -4")')
       Iflag = 4
       go to 1020
      endif

      if(idr.eq.0 .and. ecis  .and. modtyp.gt.3) then
       write(6,
     + '(" ERROR: You are trying to produce a dispersive input",
     +   " from a non-dispersive potential.")')
       write(6,
     + '(" ERROR: Set modtyp to -2")')

       Iflag = 5
       go to 1020
      endif

      if(idr.le.-2.and. ecis  .and. modtyp.gt.3) then
       write(6,
     + '(" Dispersive integrals have to be calculated numerically",
     +   " for these potentials !")')
       write(6, '(" modtyp = 4/5/6 not available, resetted to -2")')
       modtyp = 2
       open(unit=36,file='DOM-table.dat')
       open(unit=34,file='DOMtoW.dat')
       if(izproj.GT.0) open(unit=38,file='coul-vint-table.dat')
      endif
c
c     Initialize program
 130  call setup
      kecis=0
c     Call output routines
      if(modtyp.ne.1)go to 150
c     Send one energy at a time to scatip2
      isimple=2
      nesc=1
      do 146 n=1,ne
      ensc(1)=en(n)
      call scatip2
 146  continue
      go to 200
c
 150  continue
c
c     Check for Gaussian form factor for WD
      krange=iabs(jrange(4))
      if(krange.eq.0) then
        do j=1,krange
          if(rco(4,j,1).lt.0.) then
            write(6,
     +       '(" ERROR: No ECIS inputs made when Gaussian form factor ",
     +         "used for Ws",/,"iref=",i5)')iref
            Iflag = 6
            goto 1020
          endif
        enddo
      endif
c     Check for nonlocal potentials
      if(beta(1).gt.0) then
        write(6,
     +   '(" ERROR: No ECIS inputs made for non-local potentials, ",
     +     "use SCAT2000",/,"iref=",i5)')iref
        Iflag = 7
        goto 1020
      endif

      if(modtyp.ne.2)go to 160
c
c     Set up for table print
      call tableset
c
      if(imodel.eq.2)go to 159
c
c     rotational model being used
      kecis=1
      ko=15
      call ecisip(kecis)
C*****************************
C     Skipping CMS calculation
C     kecis=2
C     ko=16
C     call ecisip(kecis)
C*****************************
      go to 200
c
c     vibrational model being used
 159  kecis=3
      call ecisdwip(kecis)
      go to 200
 160  if(modtyp.ne.3) goto 180
c
c     DWBA model being used
      kecis=4
      call ecisdwip(kecis)
      goto 200
 180  if(modtyp.ne.4) goto 190
c
c     soft rotor model being used
      kecis=1  ! if kecis = 2, then incident energies are converted to CMS
      ko=15
c     Set up for table print
c     call tableset
      call optmaninp
      goto 200
 190  if(modtyp.ne.5) goto 195
c
c     Set up for table print
      call tableset
c
c     Dispersive Coupled Channel rotational model being used
      kecis=1
      ko=15
      call ecisip(kecis)
C*****************************
C     Skipping CMS calculation
C     kecis=2
C     ko=16
C     call ecisip(kecis)
C*****************************
      goto 200
195   if(modtyp.ne.6) THEN
        write(6,*) 'ERROR: WRONG TYPE OF CALCULATION'
        Iflag= 8
        goto 1020
      endif
c
c     Set up for table print
      call tableset
c
c     Dispersive Coupled Channel rotational model being used
      kecis=1
      ko=15
      call ecisip(kecis)
C*****************************
C     Skipping CMS calculation
C     kecis=2
C     ko=16
C     call ecisip(kecis)
C*****************************
 200  if(irefget.le.0) go to 100
c
c     End computation
c
      if(scat) write(ko,'(i5)')isuit
      if(ecis) write(ko,'(a3)') ifin
 1020 close(ki)
      close(ko)
C*****************************
C     Skipping CMS calculation
C     if(ecis) backspace (16,ERR=1010)
C     if(ecis) read(16,*,END=1010)
C     if(ecis) write(16,'(a3)') ifin
C     if(ecis) close(16)
C*****************************
!MS$IF DEFINED (alone)
      if(idr.ne.0 .and. modtyp.le.2 ) then
        close(36)
        close(34)
      endif
      close(37)
      if(izproj.GT.0 .and. modtyp.ne.4 ) close(38)
      close(44)
!MS$ENDIF
      close(35)
      write(6,102) Iflag,Diagnostic(4+Iflag)
      return
C1020 stop 'ERROR'
C1010 close(16,status='delete')

      END SUBROUTINE retrieve


c-----------------------------------------------------------------------------
      SUBROUTINE omin30
c-----------------------------------------------------------------------------
c
c     routine to retrieve optical model parameters from RIPL library
c
  1   format(80a1)
c
c     Zero arrays
      call setr(0.d0,epot,5*ndim1)
      call setr(0.d0,rco,5*ndim1*ndim2)
      call setr(0.d0,aco,5*ndim1*ndim2)
      call setr(0.d0,pot,5*ndim1*ndim3)
      call setr(0.d0,ecoul,ndim1)
      call setr(0.d0,rcoul,ndim1)
      call setr(0.d0,acoul,ndim1)
      call setr(0.d0,beta,ndim1)
      call setr(0.d0,rcoul0,ndim1)
      call setr(0.d0,rcoul1,ndim1)
      call setr(0.d0,rcoul2,ndim1)
      call setr(0.d0,rcoul3,ndim1)
c
      ieof=0
      read(ki,*,end=999) iref
C     write(*,*) ' Reading RIPL OMP # ',iref
      read(ki,1) (author(m),m=1,80)
      read(ki,1) (refer(m),m=1,80)
      read(ki,1) (summary(m),m=1,320)
      read(ki,*) emin,emax
      read(ki,*) izmin,izmax
      read(ki,*) iamin,iamax
      read(ki,*) imodel,izproj,iaproj,irel,idr
c     if(idr.gt.0) write(*,'(2x,i5,4x,i3,3h<Z<,i3,5x,15a1)')
c    >        iref,izmin,izmax,(author(m),m=1,15)
c
c     write(*,'(1x,8I4)') iref,imodel,izproj,iaproj,irel,idr
c
      do 100 m=1,6
        read(ki,*) jrange(m)
        if(jrange(m).eq.0) go to 100
        krange=iabs(jrange(m))
        do 98 j=1,krange
          read(ki,*) epot(m,j)

          read(ki,*) (rco(m,j,n),n=1,ndim2)
          read(ki,*) (aco(m,j,n),n=1,ndim2)
          read(ki,*) (pot(m,j,n),n=1,ndim3)
c------------------ RCN


c         if((pot(m,j,20).ne.0.) .and.      (pot(m,j,24).ge.1) .and.
c    +   (pot(m,j,14) + pot(m,j,15) + pot(m,j,16)).ne.0. ) then
c            write (*,*) iref,izmin,izmax,' Cviso=',sngl(pot(m,j,20))
c            write (*,*)
c         endif
c------------------ RCN

  98    continue
 100  continue
      read(ki,*) jcoul
      if(jcoul.le.0) go to 110
      do 108 j=1,jcoul
 108  read(ki,*) ecoul(j),rcoul0(j),rcoul(j),rcoul1(j),rcoul2(j),
     >           beta(j),acoul(j),rcoul3(j)
 110  if(imodel.ne.1) go to 130
C     Reading rigid rotor parameters
      read(ki,*) nisotop

      do 120 n=1,nisotop
        read(ki,*) iz(n),ia(n),ncoll(n),lmaxOM(n),idef(n),bandk(n),
     +    (def(n,k),k=2,idef(n),2)

        nlvmax = 0
        do 124 k=1,ncoll(n)
          read(ki,*) ex(k,n),spin(k,n),ipar(k,n)
          if(ex(k,n).gt.100.d0 .and. nlvmax.eq.0) nlvmax = k-1
 124    continue

        if(modtyp.ne.4 .and. nlvmax.gt.0) ncoll(n) = nlvmax
 120  continue
      go to 200
 130  if(imodel.ne.2) go to 150
C     Reading vibrational rotor parameters
      read(ki,*) nisotop
      do 140 n=1,nisotop
        read(ki,*) iz(n),ia(n),nvib(n)
        do 138 k=1,nvib(n)
          read(ki,*) exv(k,n),spinv(k,n),iparv(k,n),nph(k,n),defv(k,n),
     +    thetm(k,n)
 138    continue
 140  continue
      go to 200

 150  if(imodel.ne.3)go to 170

C     Reading soft rotor parameters
      read(ki,*) nisotop
      do 165 n=1,nisotop
        read(ki,*) iz(n),ia(n),ncoll(n) !  Isotope atomic and mass number read
        read(ki,*)  ! Record 3 from OPTMAN (Hamiltonian parameters)
     +  SR_hw(n),SR_amb0(n),SR_amg0(n),SR_gam0(n),SR_bet0(n),SR_bet4(n)
        read(ki,*)  ! Record 3 from OPTMAN (Hamiltonian parameters)
     +  SR_bb42(n),SR_gamg(n),SR_delg(n),SR_bet3(n),SR_et0(n),SR_amu0(n)
        read(ki,*)  ! Record 3 from OPTMAN (Hamiltonian parameters)
     +  SR_hw0(n),SR_bb32(n),SR_gamde(n),SR_dpar(n),SR_gshape(n)
        do 152 k=1,ncoll(n)
          read(ki,*) exv(k,n),spinv(k,n),iparv(k,n),
     +               SR_ntu(k,n),SR_nnb(k,n),SR_nng(k,n),SR_nno(k,n)

 152    continue
 165  continue
      go to 200

 170  if(imodel.ne.4)go to 200

C     Reading rigid-soft rotor parameters
      read(ki,*) nisotop
      do 175 n=1,nisotop
        read(ki,*) iz(n),ia(n),ncoll(n),lmaxOM(n),idef(n),bandk(n),
     +    (def(n,k),k=2,idef(n),2)
        nlvmax = 0
        do 172 k=1,ncoll(n)
          read(ki,*) ex(k,n),spin(k,n),ipar(k,n),
     +               SR_ntu(k,n),SR_nnb(k,n),SR_nng(k,n),SR_nno(k,n),
     +               defv(k,n),defr(k,n)
C         write(*,*) iz(n),ia(n),ex(k,n)
          if(ex(k,n).gt.100.d0 .and. nlvmax.eq.0) nlvmax = k-1
 172    continue
C       write(*,*) iz(n),ia(n),ncoll(n),nlvmax
C       if(modtyp.ne.4 .and. nlvmax.gt.0) ncoll(n) = nlvmax
        if(                  nlvmax.gt.0) ncoll(n) = nlvmax
C       write(*,*) iz(n),ia(n),ncoll(n),nlvmax
 175    continue
C
 200  continue
      read(ki,1,end=999)idum
      return
 999  ieof=1

      Iflag = 9
      write(6,*)
     +  ' WARNING: End of RIPL optical potential library reached.'
      return
      end subroutine omin30

c-----------------------------------------------------------------------------
      subroutine setr (a,b,n)
c-----------------------------------------------------------------------------
c     ******************************************************************
c     set all elements of array b(n) to real number a.
c     ******************************************************************
      real*8 b(n),a
      do 100 k=1,n
  100 b(k)=a
      return
      end subroutine setr

c-----------------------------------------------------------------------------
      subroutine setup
c-----------------------------------------------------------------------------
c
c     Initialize program
      integer ifirst1,ifirst2,ifirst3,ifirst4,ifirst5,ifirst6
c
      becis1='FFFFFFFFFFFFFFFFFFFFFFFFFFFTFFFFFFFFFFFFFFFFFFFFFF'
      becis2='FFFFFFFFFFFFFTFFTTTFFTTFTFFFFFFFFFFFFFFFFFFFFFFFFF'
c
      ndisx=ndim8
      ifirst1=0
      ifirst2=0
      ifirst3=0
      ifirst4=0
      ifirst5=0
      ifirst6=0
c
      izaproj=1000*izproj + iaproj
      if(izaproj.eq.1   ) ipsc=1
      if(izaproj.eq.1001) ipsc=2
      if(izaproj.eq.1002) ipsc=3
      if(izaproj.eq.1003) ipsc=4
      if(izaproj.eq.2003) ipsc=5
      if(izaproj.eq.2004) ipsc=6
      zproj=real(izproj)
      atar=float(iatar)
      ztar=float(iztar)
      eta=1.-(2.*ztar/atar)
      encoul=0.4*ztar/atar**(1./3.)
c
c     Get Fermi energy, masses and kinematic factor (amu)
      call masses
      if(irndt.eq.1) tarmas=atar
      if(irndp.eq.1) projmas=iaproj
c
c     Internal Flags
      koptom=2
c     Set koptom=1 to compute excited state potentials at incident
c         energy only. (Shortcut that is not recommended.)
c     Set koptom=2 to compute excited state potentials at incident
c         energy minus the excitation energy (corrected to lab).
c         This is the recommended option.
c
      return
      end subroutine setup

c*******************
c     SCAT2000 input
c*******************
c-----------------------------------------------------------------------------
      subroutine scatip2
c-----------------------------------------------------------------------------
c     routine to write optical model parameters in parameter lib
c
      parameter(npcofx=7)
      dimension rsc(6),resc(6),asc(6),aesc(6),potsc(6,ndim1,npcofx),
     + npzen(6),epotsc(6,ndim1)

      real*8 b(6,ndim1,15),Visov,WVisov,WSisov
c
c     To use dispersive optical model package
c
c             variables
      real*8 As,Bs,Cs,AAv,Bv,AAvso,Bvso,EEE,Ep,Ea,Ef
      real*8 alpha_PB,beta_PB,gamma_PB,Vnonl,AlphaV
      real*8 DWS,DWV,DWVso,DerDWV,DerDWS,dtmp
      real*8 WDE,WVE
      integer n,iq,nnv,nns

c     Common blocks for numerical integration
      common /energy/EEE,Ef,Ep,Ea
      common /Wenerg/WDE,WVE
      common /pdatas/As,Bs,Cs,nns,iq
      common /pdatav/AAv,Bv,nnv

      data isuit/1/
C     data ipr,ida,iba,ipu/1,0,0,1/

      data ipr,ida,iba,ipu/1,-1,0,1/ ! TO HAVE LEGENDRE POLYNOMIALS

c
      call setr(0.d0,rsc,6)
      call setr(0.d0,resc,6)
      call setr(0.d0,asc,6)
      call setr(0.d0,aesc,6)
      call setr(0.d0,potsc,6*npcofx*ndim1)
c
      el=ensc(1)
      nnv = 0
      nns = 0
      VDcoul = 0.d0

c     Get coulomb radius, if needed
      rcoulsc=0.d0
      if(jcoul.lt.1)go to 100
      jc=1
c     if(isimple.ne.2)go to 90
      do 80 j=1,jcoul
      if(el.gt.ecoul(j))jc=j+1
  80  continue
      jc=min0(jc,jcoul)
  90  rcoulsc=rcoul(jc) + rcoul0(jc)*atar**(-1.d0/3.d0) +
     +       rcoul1(jc)*atar**(-2.d0/3.d0) +
     +       rcoul2(jc)*atar**(-5.d0/3.d0) +
c            RCN addition to consider new Morillon-Romain potential
     +       rcoul3(jc)*atar
      betasc=beta(jc)
 100  encoul2=0.
      if(rcoulsc.gt.0.) encoul2=1.73d0*ztar/(rcoulsc*atar**(1.d0/3.d0))
      ipot=0

      do 150 i=1,6
c
c     For Lane consistent potentials change the incident energy for proton potentials
c     by the specified Coulomb shift (pot(1,1,25)
c
      VCshift = 0.d0
      IF(izproj.eq.1 .and. pot(i,1,25).gt.0.d0)
     >   VCshift = pot(i,1,25)*ztar/atar**(1.d0/3.d0)

      jab=iabs(jrange(i))
      ii=i
      if(i.eq.2)ii=3
      if(i.eq.3)ii=2
      npzen(ii)=jrange(i)
      if(jrange(i).eq.0)go to 150
      do 104 j=1,jab
 104  epotsc(ii,j)=epot(i,j)

      vc = 0.d0
      VVcoul = 0.d0
      VScoul = 0.d0
      DerDWV = 0.d0
      DerDWS = 0.d0
      AlphaV = 0.d0
c
c     Find energy range that energy el falls in.
 110  jp=1
      do 112 j=1,jab
      if(el.gt.epot(i,j))jp=j+1
 112  continue
      j=min0(jp,jab)
      js=1
      epotsc(ii,js)=epot(i,j)
      if(npzen(ii).gt.0)npzen(ii)=+1
      if(npzen(ii).lt.0)npzen(ii)=-1

      Ef= efermi

      if(pot(i,j,18).ne.0.) Ef=pot(i,j,18) + pot(i,j,19)*atar

      elf = el - Ef - VCshift
c
c     Calculate radius and diffuseness parameters
      if(rco(i,j,13).eq.0.) then
        rsc(ii)=abs(rco(i,j,1)) + rco(i,j,3)*eta
     *       + rco(i,j,4)/atar + rco(i,j,5)/sqrt(atar)
     *       + rco(i,j,6)*atar**(2./3.) + rco(i,j,7)*atar
     *       + rco(i,j,8)*atar**2  + rco(i,j,9)*atar**3
     *       + rco(i,j,10)*atar**(1./3.)
     *       + rco(i,j,11)*atar**(-1./3.)
C--------------------------------------------------------------------
C     RCN, 08/2004, to handle new extension to the OMP RIPL-2 format
     *       + rco(i,j,2)*el + rco(i,j,12)*el*el
      else
C     RCN, 09/2004, to handle new extension to the OMP RIPL-2 format
        nn = int(rco(i,j,7))
        rsc(ii)= ( abs(rco(i,j,1)) + rco(i,j,2)*atar ) *
     *           ( 1.d0 - ( rco(i,j,3) + rco(i,j,4)*atar ) * elf**nn/
     *           ( elf**nn + ( rco(i,j,5) + rco(i,j,6)*atar )**nn ) )
      endif
      resc(ii)=0.

      asc(ii)=abs(aco(i,j,1)) + aco(i,j,2)*el + aco(i,j,3)*eta
     *        + aco(i,j,4)/atar + aco(i,j,5)/sqrt(atar)
     *        + aco(i,j,6)*atar**(2./3.) + aco(i,j,7)*atar
     *        + aco(i,j,8)*atar**2 + aco(i,j,9)*atar**3
     *        + aco(i,j,10)*atar**(1./3.) + aco(i,j,11)*atar**(-1./3.)
      aesc(ii)=0.
C--------------------------------------------------------------------
c
      if(pot(i,j,24).eq.0.) go to 120
c
c     Special Koning-type potential formulas
c
      if (pot(i,j,24).eq.1. .or. pot(i,j,24).eq.3.) then
c
c       Koning-type formulas
c
        if(i.eq.1) call bcoget(b,j,Visov,WVisov,WSisov)
c
        elf = el - Ef - VCshift
        vc=0.d0
        if(i.eq.1 .and. b(1,j,5).ne.0.d0) then
          vc = b(1,j,1)*encoul2*( b(1,j,2) - 2.*b(1,j,3)*elf +
     +    3.*b(1,j,4)*elf**2 + b(i,j,14)*b(i,j,13)*exp(-b(i,j,14)*elf) )
            VDcoul = b(i,j,5)*vc
        endif

        nn = int(pot(i,j,13))
c       Retrieving average energy of the particle states Ep
        Ep=Ef
        if( (i.eq.2) .or. (i.eq.4) ) Ep=pot(i,j,20)
        if(Ep.eq.0.) Ep=Ef
        elf = el - Ep - VCshift

        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))

        potsc(ii,js,1)=
     +     b(i,j,1)*( b(i,j,15) - b(i,j,2)*elf + b(i,j,3)*elf**2 -
     +     b(i,j,4)*elf**3 + b(i,j,13)*exp(-b(i,j,14)*elf) ) +
     +     b(i,j,5)*vc + b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +     b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +     (elf**nn + b(i,j,10)**nn)) + b(i,j,11)*exp(-b(i,j,12)*elf)

      endif

      if (pot(i,j,24).eq.2.) then
c
c       Morillon-Romain formulas
c
        if(i.eq.1) call bcoget(b,j,Visov,WVisov,WSisov)
c
        elf = el - Ef - VCshift
        nn = int(pot(i,j,13))
c
c       Vhf(E) calculated from nonlocal approximation
c          as suggested by Perey and Buck
c
        alpha_PB = b(i,j,1)
        beta_PB  = b(i,j,2)
        gamma_PB = b(i,j,3)
        EEE = el
        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))
        vc=0.d0
        Vnonl = 0.d0
        if(i.eq.1 .or. i.eq.5) then
          Vnonl = -Vhf(EEE,alpha_PB,beta_PB,gamma_PB)
          vc = 0.d0
          if(i.eq.1 .and. b(1,j,5).ne.0.d0) then
C           MR do not use derivarive of the potential
C
C           Numerical derivative of the Vhf
C           Vnonlm = -Vhf(EEE-0.05,alpha_PB,beta_PB,gamma_PB)
C           Vnonlp = -Vhf(EEE+0.05,alpha_PB,beta_PB,gamma_PB)
C           Coulomb correction for Hartree-Fock potential
C           vc = encoul2*(Vnonlm-Vnonlp)*10.d0
C
C           MR are using constant Coulomb correction
            vc = encoul2
            VDcoul = b(i,j,5)*vc
          endif
        endif

        potsc(ii,js,1)=
     +    Vnonl + b(i,j,5)*vc +
     +    b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +    b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +    (elf**nn + b(i,j,10)**nn)) +
     +    b(i,j,11)*exp(-b(i,j,12)*elf)

      endif
c
c     Nonlocality consideration
c
c     Retrieving energy above which nonlocality in the volume absorptive
c               potential is considered (Ea)
c
      Ea=pot(i,j,21)
      if(Ea.eq.0.) Ea=1000.1d0

      if(i.eq.2 .and. Ea.lt.1000.d0) THEN
        AlphaV=pot(i,j,22)
        if(AlphaV.eq.0.d0) AlphaV=1.65d0
        if(el.gt.(Ef+Ea) ) potsc(ii,js,1) = potsc(ii,js,1) +
     +       AlphaV*(sqrt(el)+(Ef+Ea)**1.5/(2.*el)-1.5*sqrt(Ef+Ea))
      endif
c
      go to 150
c
 120  if(pot(i,j,23).eq.0)go to 130
c
c     Special Varner-type potential formulas
      potsc(ii,js,1)= (pot(i,j,1) + pot(i,j,2)*eta)/
     *    (1.+ exp((pot(i,j,3) - el + pot(i,j,4)*encoul2)/pot(i,j,5)))
      if(pot(i,j,6).eq.0.)go to 150
      potsc(ii,js,1)= potsc(ii,js,1)
     *    + pot(i,j,6)*exp((pot(i,j,7)*el - pot(i,j,8))/pot(i,j,6))
      go to 150
c
 130  if(pot(i,j,22).eq.0)go to 140
c
c     Special Smith-type potential formulas
      pi=acos(-1.d0)
      potsc(ii,js,1)=pot(i,j,1) + pot(i,j,2)*eta
     *    + pot(i,j,6)*exp(pot(i,j,7)*el + pot(i,j,8)*el*el)
     *    + pot(i,j,9)*el*exp(pot(i,j,10)*el**pot(i,j,11))
      if(pot(i,j,5).ne.0.)potsc(ii,js,1)=potsc(ii,js,1)
     *    + pot(i,j,3)*cos(2.*pi*(atar - pot(i,j,4))/pot(i,j,5))
      go to 150
c
c     Standard potential formulas
 140  potsc(ii,js,1)=pot(i,j,1) + pot(i,j,7)*eta + pot(i,j,8)*encoul
     *    + pot(i,j,9)*atar + pot(i,j,10)*atar**(1./3.)
     *    + pot(i,j,11)*atar**(-2./3.) + pot(i,j,12)*encoul2
     *    + (pot(i,j,2) + pot(i,j,13)*eta + pot(i,j,14)*atar)*el
     *    + pot(i,j,3)*el*el + pot(i,j,4)*el*el*el + pot(i,j,6)*sqrt(el)
     *    + (pot(i,j,5) + pot(i,j,15)*eta + pot(i,j,16)*el)*log(el)
     *    + pot(i,j,17)*encoul/el**2
c
 150  continue
c
      if(rco(4,1,1).lt.0.) rsc(4)=-rsc(4)
      if(nonegw.eq.0) go to 152
      if(potsc(3,js,1).lt.0.0) potsc(3,js,1)=0.
      if(potsc(4,js,1).lt.0.0) potsc(4,js,1)=0.
c
c     To calculate dispersion relation contribution
c
 152  if(abs(idr).ge.2) then
c
c       Exact calculation of the dispersive contribution
c
        EEE = el

        i=2
c       Only one energy range
        j=1
c       Real volume contribution from Dispersive relation
        DWV=0.d0
c
        if(jrange(i).gt.0 .and. pot(2,1,24).ne.0) then
          AAv=b(i,j,6)
          Bv =b(i,j,7)
          n = nint( pot(i,j,13) )
          nnv = n
          if(n.eq.0 .or. mod(n,2).eq.1) stop
     +      'ERROR: Zero or odd exponent in Wv(E) for dispersive OMP'

c         Retrieving average energy of the particle states Ep
          Ep=pot(i,j,20)
          if(Ep.eq.0.) Ep=Ef

c         analytical DOM integral
          DWV=DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE,n,DerDWV)
          if (pot(1,1,24).ne.1) DerDWV = 0.d0
c         Coulomb correction for real volume potential
          DerDWV = -b(1,1,5)*encoul2*DerDWV
C         numerical DOM derivative (not needed for a time being)
C         DWVp = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE+0.1d0,n,dtmp)
C         DWVm = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE-0.1d0,n,dtmp)
C         DerDWV = -b(1,1,5)*encoul2*(DWVp-DWVm)*5.d0
c         if(idr.le.-2) then
c           numerical DOM integral (not needed for a time being)
c           WVE=WVf(AAv,Bv,Ep,Ef,EEE,n)
c           DWV=2*DOM_int(Delta_WV,WVf,Ef,Ef+5.*Bv,150000.d0,EEE,0.d0)
c         endif
c
c         Nonlocality correction to the DOM integral
c           (only used if Ea is non-zero)
          Ea=pot(i,j,21)
          if(Ea.eq.0.) Ea=1000.1d0

          Dwplus = 0.d0
          Dwmin = 0.d0
          T12der = 0.d0
          if(Ea.lt.1000.) THEN
             AlphaV=pot(i,j,22)
             if(AlphaV.eq.0.d0) AlphaV=1.65d0
             Dwplus = AlphaV*DOM_INT_T2(Ef,Ea,EEE)
             dtmp1 = Wvf(AAv,Bv,Ep,Ef,Ef+Ea,n)
             Dwmin = dtmp1*DOM_INT_T1(Ef,Ea,EEE)
             DWV = DWV + Dwplus + Dwmin
c            Coulomb correction for nonlocal dispersive contribution
c                to real volume potential
             if(b(1,1,5).ne.0.d0 .and. pot(1,1,24).eq.1) then
               if(eee.ne.0.05d0) then
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.05d0)
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.05d0)
                 T2der = AlphaV*(T2p-T2m)*10.d0
                 T1p = DOM_INT_T1(Ef,Ea,EEE+0.05d0)
                 T1m = DOM_INT_T1(Ef,Ea,EEE-0.05d0)
                 T1der = dtmp1*(T1p-T1m)*10.d0
                 T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
              else
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.1d0)
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.1d0)
                 T2der = AlphaV*(T2p-T2m)*5.d0
                 T1p = DOM_INT_T1(Ef,Ea,EEE+0.1d0)
                 T1m = DOM_INT_T1(Ef,Ea,EEE-0.1d0)
                 T1der = dtmp1*(T1p-T1m)*5.d0
                 T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
               endif
             endif

          endif
          VVcoul = DerDWV + T12der
        endif

        i=4
c       Only one energy range
        j=1
c       Real surface contribution from Dispersive relation
        DWS=0.d0

        if(jrange(i).gt.0 .and. pot(4,1,24).ne.0) then
          n = nint( pot(i,j,13) )
          nns = n
          As=b(i,j,8)
          Bs=b(i,j,10)
          Cs=b(i,j,9)

          if(n.eq.0 .or. mod(n,2).eq.1) stop
     +      'ERROR: Zero or odd exponent in Wd(E) for dispersive OMP'
            iq=1
          if(b(4,j,12).gt.0.) iq=nint(b(4,j,12))

c         Retrieving average energy of the particle states Ep
          Ep=pot(i,j,20)
          if(Ep.eq.0.) Ep=Ef

          if(idr.ge.2) then
c           analytical DOM integral
            DWS = DOM_INT_Ws(Ef,Ep,As,Bs,Cs,EEE,n,DerDWS)
            if (pot(1,1,24).ne.1) DerDWS = 0.d0

            VScoul = -b(1,1,5)*encoul2*DerDWS
          endif

          if(idr.le.-2) then
c           numerical DOM integral
            WDE=WDf(As,Bs,Cs,Ep,EEE,n,iq)
            DWS = 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE,WDE)
c           Coulomb correction for real surface potential
            if(b(1,1,5).ne.0.d0) then
              WDE=WDf(As,Bs,Cs,Ep,EEE+0.1d0,n,iq)
              DWSp =
     >         2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE+0.1d0,WDE)
              WDE=WDf(As,Bs,Cs,Ep,EEE-0.1d0,n,iq)
              DWSm =
     >         2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE-0.1d0,WDE)
c             Numerical derivative

              DerDWS = (DWSp-DWSm)*5.d0

              if (pot(1,1,24).ne.1) DerDWS = 0.d0

              VScoul = -b(1,1,5)*encoul2*DerDWS


            endif
          endif
        endif

        i=6
c       Only one energy range
        j=1
c       Real spin orbit contribution from Dispersive relation
        DWVso=0.d0
c
        if(jrange(i).gt.0.and.pot(6,1,24).ne.0.and.abs(idr).eq.3) then
          AAvso=b(i,j,6)
          Bvso =b(i,j,7)
          n = nint( pot(i,j,13) )

          if(n.eq.0 .or. mod(n,2).eq.1) stop
     +      'ERROR: Zero or odd exponent in Wso(E) for dispersive OMP'

c         analytical DOM integral
          DWVso=DOM_INT_Wv(Ef,Ef,AAvso,Bvso,EEE,n,dtmp)

        endif

c       Adding real volume dispersive and Coulomb contribution to the real potential
c       Geometry parameters are the same as for the volume potential(imag and real)
        potsc(1,1,1) = potsc(1,1,1) + DWV + VVcoul
c       Including real surface and Coulomb dispersive contribution
c       Geometry parameters are the same as for the imaginary surface potential
        potsc(2,1,1) = DWS  + VScoul
        rsc(2)=rsc(4)
        asc(2)=asc(4)
        npzen(2)=1
c       Adding real spin orbit dispersive contribution to the real spin orbit potential
c       Geometry parameters are the same as for the imaginary spin orbit potential(imag and real)
        potsc(5,1,1) = potsc(5,1,1) + DWVso

      endif
c
c     Write SCAT2 input
c
      if(nsc2.ge.1)go to 156
 154  itmp=irel
      if(irel.eq.2) itmp=1
c
c     Dispersion relations are considered explicitly in this code
c     so SCAT2 does not need to make any DR related calculation
c         => we are setting idr = 0 always
      write(ko,'(6i5)') ipr,ida,iba,ipu, 0 ,itmp
      nsc2=nsc2+1
 156  write(ko,'(i5)') nesc
      ensc(1)=-ensc(1)
      write(ko,'(7f10.5)') (ensc(j),j=1,nesc)
      ensc(1)=-ensc(1)
      write(ko,'(2i5)') iztar,iatar
      write(ko,'(2i5)') ipsc,ipot
c
c     Factor coming from Dirac equation reduction (relativistic approximation).
c     Using this factor force us to employ relativistic kinematics as well.
      gamma=1.d0
      if(irel.eq.2) then
C
C       Follwing Madland private comm. 28/09/2008
C       it should be applied to all potentials when used
C       including spin-orbit
C
c       Target system mass in MeV
        EMtar=tarmas*amu0c2
c       Total system mass in MeV
        EMtot=(tarmas+projmas)*amu0c2
c       Total kinetic energy in cm
        Tcm=sqrt(2*EMtar*el + EMtot**2) - EMtot
c       Relativistic correction to the potential (non relativistic target!!).
        gamma=1.d0+Tcm/(Tcm+2*projmas*amu0c2)
      endif

      do 200 i=1,6
      write(ko,'(4f9.5,i5)') rsc(i),resc(i),asc(i),aesc(i),npzen(i)
      do 160 j=1,iabs(npzen(i))
        write(ko,5)epotsc(i,j),(gamma*potsc(i,j,k),k=1,npcofx)
  5   format(f9.3,1p,4(1x,e12.5),/,9x,3(1x,e12.5))
 160  continue
 200  continue
C     write(ko,'(3f10.5)') rcoulsc,efermi,betasc
      write(ko,'(3f10.5)') rcoulsc,ef    ,betasc  ! Bug corrected on 20 Dec 2011
      write(ko,'(i5)') isuit
      return
      end subroutine scatip2

c*************************
c     OPTMAN input
c*************************
c-----------------------------------------------------------------------------
      subroutine optmaninp
c-----------------------------------------------------------------------------
c
c     ROUTINE TO PREPARE AND WRITE OPTMAN INPUT
      real*8 nang(150), AlphaV !,dtmp

      real*8 CAVRss,CARRss,CAARss,CARDss,CAACss
      integer k

      if((imodel.eq.1 .or. imodel.eq.3 .or. imodel.eq.4) .and.
     +     idr.ne.0  .and. pot(1,1,24).ne.1.d0 ) then
        write(ko,*) 'ERROR: RIPL to OPTMAN interface is not prepared yet
     + to deal with this potential dependence (Contact R.CapoteNoy@iaea.
     +org)'
        Iflag = 10
        return
      endif

      if (imodel.eq.2) then
        write(ko,*)
     +   'ERROR:OPTMAN cant not calculate CC vibrational model'
        Iflag = 11
        return
      endif

C     if (idr.eq.0 .and. izproj.ne.0 .and. imodel.ne.3) then
C       write(ko,*)
C    +   'OPTMAN cant not calculate non-dispersive CC proton OMPs'
C       stop 6
C     endif

      if((imodel.eq.1 .or. imodel.eq.3 .or. imodel.eq.4) .and.
     +     idr.ne.0  .and. nint(pot(2,1,13)).ne.nint(pot(4,1,13)) ) then
        Iflag = 12
        write(ko,*) 'ERROR: RIPL to OPTMAN interface can not calculate d
     +ispersive effects using Nv not equal Ns (Contact R.CapoteNoy@iaea.
     +org)'
        return
      endif

C==============================================================================
C     OPTMAN input parameters
C     OPTMAN 1st line
      mejob = 1 ! OMP calculation, no fitting
      mepot = 2 ! potential expanded by derivatives
      meham = 5 ! Davidov-Chaban(3), Davidov-Filippov(4), Soft rotor(5)
      mesho = 2
      mehao = 2
      mesha = 4

C mesol > 3 - solution using iterations with exact coupled channels
C solution with the number of coupled states equal MESOL as zero
C approximation, MESOL must be less than or equal to 20;
      mesol = 5 ! automatic selection
c     mesol = 2 ! automatic selection

      if(imodel.eq.4) then ! rigid + soft rotor (4)
        mepot = 1
        meham = 1
C mesol > 3 - solution using iterations with exact coupled channels
C solution with the number of coupled states equal MESOL as zero
C approximation, MESOL must be less than or equal to 20;
        mesol =20
        mesha = 1
        mesho = 0
        mehao = 0
C       exact calculation for odd and odd-odd targets
        if ( mod(iatar-iztar,2).ne.0 .or. mod(iztar,2).ne.0 ) mesol = 2
      endif

      if(imodel.le.1) then ! rigid rotor (1) or spherical (0)
        mepot = 1
        meham = 1
        mesha = 1
        mesho = 0
        mehao = 0
C       exact calculation for odd and odd-odd targets
        if ( mod(iatar-iztar,2).ne.0 .or. mod(iztar,2).ne.0 ) mesol = 2
      endif

      mepri = 0 ! short output
      meapp = 0 ! 1 (No energy losses)
      mevol = 0
      merel = 0
      if(irel.eq.1) merel = 2  ! relativistic kinematic only
      if(irel.eq.2) merel = 1  ! relativistic kinematic + potential dependence

      mecul = 0 ! Energy dependent i.e. DVc = - Ccoul*Z/A**(1/3) * dV/dE
C     mecul = 1 ! Energy independent Coulomb correction DVc = Ccoul*Z/A**(1/3)
C     mecul = 2 ! E = Ep - Ccoul*Z/A**(1/3) for real and imag potentials
C     mecul = 3 ! E = Ep - Ccoul*Z/A**(1/3) for real potential

      if(pot(1,1,25).ne.0.d0) mecul =3
      if(pot(2,1,25).ne.0.d0 .or. pot(4,1,25).ne.0.d0) mecul =2

      merzz = 1 ! Constant charge radius(0) Energy dependent(1)
      merrr = 1 ! Real potential radius is constant(0) Energy dependent(1)
      medis = 0 ! No dispersion
      if (abs(idr).eq.3) medis = 1
      if (abs(idr).eq.2) medis = 2
C     If merip = 1 then potentials defined at every energy (loop over energies)
      merip = 0
      if(idr.eq.0 .and. imodel.lt.4) merip = 1
c------------------------------------------------------------------------------
      jchk=1
      if(kecis.ne.1)jchk=0

      ncol=1
      if(imodel.ne.0) then
        call couch11
        ncol = ncoll(ntar)
      endif
C
      el=en(1)
      i=1
      call optmod(i,kecis)
          if(Iflag.gt.0) return
C
C     Checking if OPTMAN can deal with them.
C
      if((rvso.ne.0.d0 .and. rwso.ne.0.d0) .and.
     +     rvso.ne.rwso ) then
        write(6,*) 'ERROR: OPTMAN can not deal with different radii for
     +real and imaginary SO potentials'
           Iflag = 14
       return
      endif

      if((avso.ne.0.d0 .and. awso.ne.0.d0) .and.
     +     avso.ne.awso ) then
       write(6,*) 'ERROR: OPTMAN can not deal with different diffuseness
     +  for real and imaginary SO potentials'
           Iflag = 15
       return
      endif

      call optname

C-----Writing OPTMAN input

C-----CARD 1 : Title

      IF(imodel.eq.0) then
       write(ko,'(f10.5,"-",f10.5," MeV ",a8," on ",i3,a2,
     + ": rigid rotor, OMP ",14a1," RIPL REF#=",i5)')
     + en(1),en(ne),parname(ipsc),iatar,nuc(iztar),opname,iref
      ELSEIF(imodel.eq.1) then
       write(ko,'(f10.5,"-",f10.5," MeV ",a8," on ",i3,a2,
     + ": rigid rotor, OMP ",14a1," RIPL REF#=",i5)')
     + en(1),en(ne),parname(ipsc),iatar,nuc(iztar),opname,iref
      ELSEIF(imodel.eq.3) then
       write(ko,'(f10.5,"-",f10.5," MeV ",a8," on ",i3,a2,
     + ": soft rotor, OMP ",14a1," RIPL REF#=",i5)')
     + en(1),en(ne),parname(ipsc),iatar,nuc(iztar),opname,iref
      ELSEIF(imodel.eq.4) then
        write(ko,'(f10.5,"-",f10.5," MeV ",a8," on ",i3,a2,
     + ": rigid+soft rotor, OMP ",14a1," RIPL REF#=",i5)')
     + en(1),en(ne),parname(ipsc),iatar,nuc(iztar),opname,iref
      ENDIF

      write(ko,'(20i2.2)')
     +mejob,mepot,meham,mepri,mesol,mesha,mesho,mehao,
     +meapp,mevol,merel,mecul,merzz,merrr,medis,merip,0,0,0,0
C
C     Soft rotor hamiltonian
C
      if(meham.gt.1 .and. ntar.gt.0) then
        write(ko,'(6E12.5)') ! Record 3 from OPTMAN (Hamiltonian parameters)
     +    SR_hw(ntar),SR_amb0(ntar),SR_amg0(ntar),
     +    SR_gam0(ntar),SR_bet0(ntar),SR_bet4(ntar)
        write(ko,'(6E12.5)') ! Record 3 from OPTMAN (Hamiltonian parameters)
     +    SR_bb42(ntar),SR_gamg(ntar),SR_delg(ntar),
     +    SR_bet3(ntar),SR_et0(ntar),SR_amu0(ntar)
        write(ko,'(6E12.5)') ! Record 3 from OPTMAN (Hamiltonian parameters)
     +    SR_hw0(ntar),SR_bb32(ntar),SR_gamde(ntar),
     +    SR_dpar(ntar),SR_gshape(ntar)
      endif

      npd = 8     ! maximum multipole in use
      las = 8
      if(imodel.eq.1 .or. imodel.eq.4) npd = idef(ntar)

C     if(imodel.eq.1 .or. imodel.eq.4) then

C       npd = idef(ntar)

C       las = npd
C      endif

      if(imodel.eq.3) then
        npd = 4   ! maximum multipole in use
        las = 4
      endif

      if(imodel.eq.0) then
        npd = 0
        las = 0
      endif
C     mtet = 91
      mtet = 10
C     mtet = 0

C     KODMA 0 needed as coupled states are not ordered

C     write(ko,'(9I3)') ncol,ne,npd,las,mtet,90,200,180,0
      write(ko,'(9I3)') ncol,ne,npd,las,mtet,90,200,180,1  ! KODMA = 1

      write(ko,'(6e12.5)') (en(i),i=1,min(50,ne))
      write(ko,'(36I2.2)') (izproj,i=1,min(50,ne))
      if(mtet.gt.1) then
        nang(1) = 0.d0
        if(izproj.gt.0) nang(1)=0.5d0
        do i=2,mtet
          nang(i) = dble(180.d0/(mtet-1))*(i-1)
        enddo
        write(ko,'(6e12.5)') (nang(i),i=1,mtet)
      endif

      if(imodel.ge.1) then
        if(meham.eq.1 .and. imodel.ne.4) then ! rigid rotor

          do k=1,ncoll(ntar)
              if(ex(k,ntar).le.100.d0) then
C             Normal excited states (for n,n or p,p scattering)
              write(ko,'(F11.7,1x,4I2)')
     +          ex(k,ntar),nint(2*spin(k,ntar)),ipar(k,ntar),
     +                            nint(2*spin(1,ntar)), 0
            else
C             Isobar analogue states (for p,n scattering)
              write(ko,'(F11.7,1x,4I2)')
     +          ex(k,ntar)-100.d0,nint(2*spin(k,ntar)),ipar(k,ntar),
     +                            nint(2*spin(1,ntar)), 1
            endif
          enddo
        else                    ! soft rotor
          if(imodel.ne.4) then ! soft rotor

            do k=1,ncoll(ntar)
              write(ko,'(F11.7,1x,6I2)') exv(k,ntar),
     +         nint(2*spinv(k,ntar)),iparv(k,ntar),
     +            SR_ntu(k,ntar), SR_nnb(k,ntar),
     +            SR_nng(k,ntar), SR_nno(k,ntar)
            enddo
          else                 ! rigid (or rigid-soft) rotor
C     READ(20,3)(EL(I),JO(I),NPO(I),KO(I),NCA(I),
C    *     NUMB(I),BETB(I),I=1,NUR)
C   3 FORMAT(E12.7,5I2,E12.7)
            do k=1,ncoll(ntar)
C             write(ko,'(e12.5,5I2,2E12.5)') ex(k,ntar),
              write(ko,'(F11.7,1x,5I2,2E12.5)') ex(k,ntar),
     +         nint(2*spin(k,ntar)),ipar(k,ntar),
     +            SR_ntu(k,ntar), SR_nnb(k,ntar),
     +            SR_nng(k,ntar), defv(k,ntar), defr(k,ntar) ! SR_nno(k,ntar)
                ! SR_nno(k,ntar) is always zero in the rigid rotor model
            enddo
          endif
        endif
      else  ! spherical potential
        write(ko,'(F11.7,1x,4I2)')  0.d0,0,+1,0
      endif
C     This line is new for OPTMAN v.12 (from November 2011)
C     Flags for inclusion of resonances
      write(ko,'(I3)') 0   ! no resonances are considered

C     Ef=dble(int(100000*efermi))/100000
      Ef=efermi
      if(pot(1,1,18).ne.0.) Ef=pot(1,1,18) + pot(1,1,19)*atar
C    >  Ef=dble(int(100000*(pot(1,1,18) + pot(1,1,19)*atar)))/100000

      write(ko,'(6F12.5)') projmas,projspi,tarmas,ztar,Ef,Ef

      if(izproj.eq.0) then
        rc = 1.d0
        ac = 1.d0
      endif

      CAVRss  = 0.d0
      CARRss  = 0.d0
      CARDss  = 0.d0
      CAARss  = 0.d0
      CAACss  = 0.d0

      iaref = 0
      IF(pot(1,1,23).gt.0 .and. NINT(pot(1,1,24)).eq.1)
     >  iaref = NINT(pot(1,1,23))
      IF(pot(2,1,23).gt.0 .and. NINT(pot(2,1,24)).eq.1)
     >  iaref = NINT(pot(2,1,23))
      IF(pot(4,1,23).gt.0 .and. NINT(pot(4,1,24)).eq.1)
     >  iaref = NINT(pot(4,1,23))

C     ALFNEW,VRD,CAVR,CARR,CAAR,CARD,
C     CAAC,ATI
C     dtmp = ATIS(IIIS)-ATIS(1)
C     VRLA=VRG+CAVR*dtmp
C     RR=(RRG+CARR*dtmp)*ASQ
C     AR0=ARG+CAAR*dtmp
C     RD=(RDG+CARD*dtmp)*ASQ
C     AC0=ACG+CAAC*dtmp

      IF(iaref.gt.0) THEN
        CAVRss = pot(1,1,22) ! Real vol HF potential
        CARRss = rco(1,1,7)  ! Real vol radius
        CAARss = aco(1,1,7)  ! Real vol diffuseness
        CARDss = rco(4,1,7)  ! Real sur radius
        CAACss = aco(2,1,7)  ! Imag vol diffuseness
      ENDIF
C
C     Dispersive deformed potentials are treated exactly.
C     All others potentials are calculated at a fixed energy,
C     therefore loop over incident energies is needed.
C
      if(imodel.ge.3 .or. (idr.ne.0 .and. imodel.eq.1)) then
        write(ko,'(6F12.5)')
C               Vr0                         Vr1
     +    pot(1,1,14)+pot(1,1,15)*iatar, -pot(1,1,3)-pot(1,1,4)*iatar,
C               Vr2                         Vr3
     +    pot(1,1,5) +pot(1,1,6) *iatar, -pot(1,1,7),
C          VRLA           ALAVR
     +    pot(1,1,16),  pot(1,1,17)

        write(ko,'(6F12.5)')  0.d0, 0.d0, 0.d0,
C                        W0
     +    pot(4,1,1) +pot(4,1,7)*iatar +pot(4,1,9)*iatar**(-1.d0/3.d0),
C            Bs          Cs
     +    pot(4,1,6), pot(4,1,2)

        write(ko,'(6F12.5)')  0.d0, 0.d0, 0.d0,
C             Av=b(2,j,6)                  Bv=b(2,j,7)
     +    pot(2,1,1)+pot(2,1,2)*iatar, pot(2,1,3)+pot(2,1,4)*iatar, 0.d0

        write(ko,'(6F12.5)')
C               Vso                        Lso
     +     pot(5,1,10)+pot(5,1,11)*iatar, pot(5,1,12), 0.d0, 0.d0,
C                 Wso        Bso
     +     pot(6,1,1), pot(6,1,3)

        if(iaref.gt.0) THEN
          write(ko,'(6F12.5)')   rv - rco(1,1,7)*atar, 0.d0, 0.d0,
     +              pot(2,1,13), av - aco(1,1,7)*atar, 0.d0
          write(ko,'(6F12.5)')
     +        rwd - rco(4,1,7)*atar,  awd - aco(4,1,7)*atar, 0.d0,
     +        rw  - rco(2,1,7)*atar,   aw - aco(2,1,7)*atar, 0.d0
          write(ko,'(6F12.5)')1.d0, 1.d0, 0.d0, rvso, avso, 0.d0
        else
C                                                   nv
          write(ko,'(6F12.5)')  rv, 0.d0, 0.d0, pot(2,1,13), av, 0.d0
          write(ko,'(6F12.5)') rwd,  awd, 0.d0,   rw,   aw, 0.d0
          write(ko,'(6F12.5)')1.d0, 1.d0, 0.d0, rvso, avso, 0.d0
        endif

        Ccoul = 0.d0
        if(izproj.gt.0 .and. rc.gt.0.d0 .and. mecul.eq.0)
     +        Ccoul =  pot(1,1,9)*1.73/rc
        if(izproj.gt.0 .and. mecul.ge.2) Ccoul =  pot(1,1,25)
        write(ko,'(6F12.5)')  rc, 0.d0, 0.d0, ac,Ccoul, 1.d0

        write(ko,'(6F12.5)')
C                 Cviso          Cwiso                      Ea
     +    abs(pot(1,1,20)),abs(pot(4,1,8)), 0.d0, pot(2,1,21)

        AlphaV = pot(2,1,22)
        if(pot(2,1,21).ne.0.d0 .and. AlphaV.eq.0.d0) AlphaV=1.65d0

C       This line is new for OPTMAN v.10 (from March 2008)
C       write(ko,'(6g12.5)') AlphaV, 0.d0
C       Modified for OPTMAN v.12
        write(ko,'(6F12.5)') AlphaV,0.d0,CAVRss,CARRss,CAARss,CARDss

C       This line is new for OPTMAN v.12 (from November 2011)
        if(iaref.eq.0) then
          write(ko,'(6F12.5)') CAACss, dble(iatar)
        else
          write(ko,'(6F12.5)') CAACss, dble(iaref)
        endif
        if(meham.eq.1 .or. imodel.eq.4 ) ! rigid rotor
     +      write(ko,'(6g12.5)') (def(ntar,i),i=2,las,2)
        return

      endif

      write(ko,'(6F12.5)')   v, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      write(ko,'(6F12.5)')  wd, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      write(ko,'(6F12.5)')   w, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      write(ko,'(6F12.5)') vso, 0.d0,  wso, 0.d0, 0.d0, 0.d0
      write(ko,'(6F12.5)')  rv, 0.d0, 0.d0, 0.d0,   av, 0.d0
      write(ko,'(6F12.5)') rwd,  awd, 0.d0,   rw,   aw, 0.d0
      write(ko,'(6F12.5)')1.d0, 1.d0, 0.d0, rvso, avso, 0.d0
      write(ko,'(6F12.5)')  rc, 0.d0, 0.d0,   ac, 0.d0, 1.d0
C                         Cviso  Cwiso       Ea
      write(ko,'(6F12.5)') 0.d0, 0.d0, 0d0, 0.d0

C     These lines are new for OPTMAN v.10 (from March 2008)
C     write(ko,'(6g12.5)') AlphaV, 0.d0
C      Modified for OPTMAN v.12
      write(ko,'(6F12.5)') AlphaV,0.d0,CAVRss,CARRss,CAARss,CARDss
C     This line is new for OPTMAN v.12 (from November 2011)
      write(ko,'(6F12.5)') CAACss, dble(iaref)

      if(meham.eq.1 .and. imodel.gt.0) ! rigid rotor
     +      write(ko,'(6g12.5)') (def(ntar,i),i=2,las,2)
      if(merip.eq.0) return

      do i = 2, ne
        el=en(i)
        if(kecis.eq.2)el=en(i)/xkine(el)
        call optmod(i,kecis)
            if(Iflag.gt.0) return

        if(izproj.eq.0) then
          rc = 1.d0
          ac = 1.d0
        endif

        write(ko,*) 'POTENTIAL PARAMETERS at Eproj=',el
        write(ko,'(6F12.5)')   v, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
        write(ko,'(6F12.5)')  wd, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
        write(ko,'(6F12.5)')   w, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
        write(ko,'(6F12.5)') vso, 0.d0,  wso, 0.d0, 0.d0, 0.d0
        write(ko,'(6F12.5)')  rv, 0.d0, 0.d0, 0.d0,   av, 0.d0
        write(ko,'(6F12.5)') rwd,  awd, 0.d0,   rw,   aw, 0.d0
        write(ko,'(6F12.5)')1.d0, 1.d0, 0.d0, rvso, avso, 0.d0
        write(ko,'(6F12.5)')  rc, 0.d0, 0.d0,   ac, 0.d0, 1.d0
C                           Cviso  Cwiso       Ea
        write(ko,'(6F12.5)') 0.d0, 0.d0, 0d0, 0.d0
C       These lines are new for OPTMAN > v.10 (from March 2008)
C       write(ko,'(6e12.5)') AlphaV, 0.d0
C       Modified for OPTMAN v.12
        write(ko,'(6F12.5)') AlphaV,0.d0,CAVRss,CARRss,CAARss,CARDss
C       This line is new for OPTMAN v.12 (from November 2011)
        write(ko,'(6F12.5)') CAACss, dble(iaref)

       enddo

      return
      end subroutine optmaninp

c*************************
c     ECIS96-ECIS06 input
c*************************
c-----------------------------------------------------------------------------
      subroutine ecisip(kecis)
c-----------------------------------------------------------------------------
c
c     ROUTINE TO PREPARE AND WRITE ECIS INPUT
c
c     To use dispersive optical model package
c             variables
c
      real*8 As,Bs,Cs,AAv,Bv,EEE,Ep,Ea,Ef,AAvso,Bvso,alphav
      real*8 VHFnum,DWV,VSOnum,DWVso,ALhf,ALso,DWS,ALpha,RStep
      integer nv,iq,ns,nl,nn2, kecis
      real*8 wsref,wvref,wsoref,vsoref,vvref
      real*8 fv,fs,fso,tv,ts,tso

      common /energy/EEE,Ef,Ep,Ea
      common /pdatas/As,Bs,Cs,ns,iq
      common /pdatav/AAv,Bv,nv
      common /pdatao/AAvso,Bvso,nl
      common /dispcor/VHFnum,DWV,VSOnum,DWVso,ALhf,ALso,DWS,ALpha,RStep

C     astep=10.
      astep=20.

      ncol=1
      iterm=20
      npp=1

C     rmatch= 25.d0
      rmatch = 1.35d0*iatar**(1./3.) +  10.d0*0.65d0
      if(rmatch.lt.25.d0) rmatch = 25.d0
      if (izaproj.ne.1) rmatch = 40.d0

      ecis1=becis1
      ecis2=becis2

      ecis2(1:1)='T'
      ecis2(9:9)='T'

      ecis2(15:15)='T'
      ecis2(20:20)='T'
      if(imodel.eq.0)ecis2(40:40)='T'
      if(modtyp.eq.3) ecis1(30:30)='T' ! DWBA
      if (izaproj.ne.1) ecis2(10:10)='T'
      if (irel.ge.1) THEN
        ecis1(8:8)='T'
        ecis2(45:45)='T'  ! reduced energy is used
      endif

      if (modtyp.ge.5) then
C       Preparing dispersive input for ECIS06
        ecis1(10:10)='T'! ENERGY DEPENDENT POTENTIALS BY DISPERSION RELATIONS (GS)
        ecis1(20:20)='T'! ENERGY DEPENDENT POTENTIALS BY DISPERSION RELATIONS (EXC.LEV.)
C       ecis1(37:37)='T'! NEXT CALCULATION CHANGING ONLY ENERGY AND SOME OMP
        ecis2(9:9)='T'  ! TRANSMISSION COEFFICIENTS FOR THE GS WRITTEN ON FILE 56
      endif

C-----Smatrix output
      ecis2(6:6) = 'T'

C     23- LO(23) NO USE OF PADE APPROXIMANT RESULTS AND SHIFT TO USUAL  ECIS-101
C                COUPLED EQUATIONS WHEN CONVERGENCE IS NOT OBTAINED.    ECIS-102
      ecis1(23:23)='T'

C     if(imodel.eq.1) then
      if(imodel.ne.3) then  ! for all models except soft-rotor
        call couch11
        ncol=1
        if(ntar.gt.0)ncol=ncoll(ntar)
        npp=ncol
        if(koptom.eq.1)npp=1
      endif

C     21- LO(21) USUAL COUPLED EQUATIONS-(INVERSE: ITERATIONS).         ECIS-095
C                NOT ALLOWED WITH DIRAC EQUATION. WHEN IT IS USED WITH  ECIS-096
C                DEFORMED SPIN-ORBIT, THE DERIVATIVE TERMS ARE NOT TAKENECIS-097
C                INTO ACCOUNT AND THE COMPUTATION IS INCORRECT.         ECIS-098
      convg=1.0d-8
      ecis1(21:21)='F'
      if(el.lt.10.d0) then
        ecis1(21:21)='T'
        convg=1.0d-10
      else
        if(ntar.gt.0) then
            if(spin(1,ntar).ge.0.d0) then
            ecis1(21:21)='T'
            convg=1.0d-10
          endif
        endif
      endif

      ecis2(13:13)='T'  ! Tlj always calculated

      if(kecis.ne.1) then
        ecis2(13:13)='T'
        ecis2(14:14)='F'
        ecis2(15:15)='F'
        astep=180.
      endif
C
C     Loop over incident energies
C
      do n=1,ne
        el=en(n)
        if(kecis.eq.2)el=en(n)/xkine(el)
c
c       Check for gaussian ff for WD
        call gaussff(igauflg)
        if(igauflg.eq.1) then
          Iflag = 6
          write(6,*) ' ERROR: Gaussian form factor requested in ECIS'
          return
        endif
c
        if(imodel.eq.0) iterm=1  ! one iteration forced for spherical potentials
        if(imodel.ne.1.or.ntar.eq.0)go to 54
C       if(spin(1,ntar).ge.2.5.and.el.lt.0.1)convg=1.0d-10
  54    ldwmax=2.4*1.25*iatar**(1./3.)*0.22*sqrt(projmas*en(n))
        njmax=max(2*ldwmax,20)
        jchk=1
        if(kecis.ne.1)jchk=0
        call optmod(1,kecis)
            if(Iflag.gt.0) return
        call optname
        write(ko,'(f7.3," MeV ",a8," on ",i3,a2," - Optical model: ",
     >   14a1," REF#=",i5)') el,parname(ipsc),iatar,nuc(iztar),
     >   opname,iref
        write(ko,'(a50)') ecis1
        write(ko,'(a50)') ecis2
C       Only one potential for a full dispersive calculation
        nppaa = npp
        if (modtyp.eq.6) nppaa = 1
        if(modtyp.ne.5) then
          write(ko,'(4i5)') ncol,njmax,iterm,nppaa
        else
          write(ko,'(4i5)') ncol,njmax,iterm,1
        endif

C       write(ko,'(10x,f10.5,10x,1p,3(2x,e8.1))')
C    +    rmatch,convg,convg,convg

        write(ko,'(2f10.5,10x,1p,3(2x,e8.1))')
     +    RStep,rmatch,convg,convg,convg

        if (ecis2(15:15).eq.'T') write(ko,'( )')
        if(imodel.eq.1) then
          call couch22
        else
          tgts=0.
          tgtp='+'
          write(ko,'(f5.2,2i2,a1,5f10.5)') tgts,0,1,tgtp,el,
     +      projspi,projmas,tarmas,ztar*zproj
          write(ko,'( )')

        endif
C
C       LOOP OVER EXCITED LEVELS
C
        do j=1,nppaa
          if(nppaa.ne.1) then
            el=en(n)
            if(kecis.eq.2)el=el/xkine(el)
            xratio = tarmas / (tarmas+projmas)
            if(ntar.gt.0) el=el-ex(j,ntar)/xratio
            if(el.lt.0. .and. modtyp.ne.5) el=0.0001
            jchk=0
            if(j.gt.1) then
                          call optmod(j,kecis)
                  if(Iflag.gt.0) return
            endif
          endif

          if(modtyp.ge.5) then
            v   = VHFnum
            vd  = 0.d0
            vso = VSOnum
          endif

          if(.not.(modtyp.eq.5 .and. j.gt.1)) then
            write(ko,'(3f10.5)') v,rv,av
            write(ko,'(3f10.5)') w,rw,aw
            write(ko,'(3f10.5)') vd,rvd,avd
            write(ko,'(3f10.5)') wd,rwd,awd
            write(ko,'(3f10.5)') vso,rvso,avso
            write(ko,'(3f10.5)') wso,rwso,awso

C     COULOMB POTENTIAL
C      1-10   VAL(20)  REDUCED COULOMB RADIUS IN FERMIS.                ECIS-724
C     11-20   VAL(21)  DIFFUSENESS OF A WOODS-SAXON CHARGE DISTRIBUTION.ECIS-725

            write(ko,'(3f10.5)') rc,ac,0.
            write(ko,'(3f10.5)') 0.,0.,0.
          endif

          if (modtyp.eq.5 .and. j.eq.1) then
C           Storing reference depths for the first level
            wvref  = w
            wsref  = wd
            wsoref = wso
            vvref  = v
            vsoref = vso
          endif

          if (modtyp.ne.5 .and. j.eq.nppaa)
     +                write(ko,'(3f10.5)') 0.d0,astep,180.d0


          if (modtyp.eq.5) then
C           ECIS 06
C           Dispersive corrections read externally
C
            if(j.eq.1) then
              write(ko,'(3f10.5)') 0.,astep,180.
              nn2 = 2
              if(ea.ge.1000.d0) nn2=0 ! there is non non-locality
              write(ko,'(14I5)')  kecis,  -nn2,  nv,  ns,  -nl,  1
            endif

            fvv = 0.d0
            fvs = 0.d0
            fv  = 0.d0
            fs  = 0.d0
            fso = 0.d0
            tv  = 0.d0
            ts  = 0.d0
            tso = 0.d0

            if(vvref.ne.0.d0) fvv = (v-vvref)/vvref
            if(wvref.ne.0.d0) then
              fv = DWV/wvref
              tv = (w -wvref)/wvref
            endif

            if(wsref.ne.0.d0) then
              fs = DWS/wsref
              ts = (wd -wsref)/wsref
            endif
            if(wsoref.ne.0.d0) then
              fso = DWVso/wsoref
              tso = (wso -wsoref )/wsoref
            endif
            if(vsoref.ne.0.d0) fvs = (vso-vsoref)/vsoref

            write (ko,'(2(G10.4,F10.4,G10.4))')
     >        tv, fv  , fvv,
     >        ts, fs  , 0.d0
            write (ko,'(2(G10.4,F10.4,G10.4))')
     >        tso, fso, fvs
          endif

          if (modtyp.eq.6) then
C           ECIS 06
C           Full dispersive corrections internally calculated
            nn2 = 2
            alphav = ALpha

            if(ea.ge.1000.d0) then
              nn2=0
              ea = 0.d0
              alphav = 0.d0
            endif

            if(abs(idr).ge.2) nl = -nl
            write(ko,'(14I5)')  kecis ,  -nn2,  nv,  ns,  -nl,  0
            write(ko,'(7F10.5)')   el ,    ef,  ep,  ea, ALso, 0.d0,  Bv
            write(ko,'(7F10.5)')alphav, 0.d0, Bs,  Cs, 0.d0, Bvso, ALhf
          endif

        enddo ! LOOP OVER EXCITED STATES
      enddo ! LOOP OVER ENERGIES
      return
      end subroutine ecisip

c-----------------------------------------------------------------------------
      subroutine ecisdwip(kecis)
c-----------------------------------------------------------------------------
c
c     Routine to prepare and write ECIS inputs for vibrational and DWBA cases
c
      character*4 lmodel
      character*1 tgtpx
      real*8 xratio
      integer kkk,kecis

  14  ncol=2
      iterm=20
C     npp=2
C     astep=10.
      astep=20.
      if(koptom.eq.1)npp=1
C     rmatch= 25.d0
      rmatch = 1.35d0*iatar**(1./3.) +  10.d0*0.65d0
      if(rmatch.lt.25.d0) rmatch = 25.d0
      if (izaproj.ne.1) rmatch = 40.d0


      ecis1=becis1
      ecis2=becis2


C     RCN 08/2004  Output of the cross sections requested
      ecis2(9:9)='T'
      ecis1(12:12)='T'
      ecis2(20:20)='T'
      if(imodel.eq.0) ecis2(40:40)='T'
      if (izaproj.ne.1) ecis2(10:10)='T'
      if (irel.ge.1) THEN
        ecis1(8:8)='T'
        ecis2(45:45)='T'  ! reduced energy is used
      endif

      if(modtyp.eq.2) call vib1

      kexit=0
      if(modtyp.eq.3) iterm=1
      if(modtyp.eq.3) call dwba1(kexit)
      if(kexit.eq.1) return
      lmodel='VIBR'
      if(modtyp.eq.3)lmodel='DWBA'

      ecis2(14:14)='F'
      ecis2(15:15)='F'
      ecis2(30:30)='T'

C-----Smatrix output
      ecis2(6:6) = 'T'

      xratio = tarmas / (tarmas+projmas)
      tgts=sdis(1)
      tgtp='+'
      if(tgts.eq.0.0.and.ipdis(1).ge.0)go to 60
          Iflag = 13
      write(6,'(" ERROR: Nucleus izatar=",i6,": VIB and DWBA options onl
     +y set up for even-even targets.")') izatar
      return
  60  do kkk=1,ne
        eninc=en(kkk)
        ecm=xkine(eninc)*eninc
c
c       Check for gaussian ff for WD
        el=eninc
        call gaussff(igauflg)
        if(igauflg.eq.1) then
            Iflag = 6
            write(6,*) ' ERROR: Gaussian form factor requested in ECIS'
          return
        endif

        ldwmax=2.4*1.25*iatar**(1./3.)*0.22*sqrt(projmas*eninc)
        njmax=max(2*ldwmax,20)
        do j=2,ndis
          npp = 2  ! only one excited state is coupled at the time
c         if (ecm.le.edis(j)+0.0001) cycle
          if (ecm.le.edis(j)+0.0001) npp = 1
          tgtpx='+'
          if(ipdis(j).lt.0)tgtpx='-'
          write(ko,'(f6.2," MeV ",a8," on ",$)') eninc,parname(ipsc)
          write(ko,'(i3,a2,": ",a4," for ",$)')iatar,nuc(iztar),lmodel
          write(ko,'(f4.1,a1," level at ",f5.3," MeV")')
     +     sdis(j),tgtpx,edis(j)
          write(ko,'(a50)') ecis1
          write(ko,'(a50)') ecis2
          write(ko,'(4i5)') ncol,njmax,iterm,npp

          convg=1.0d-8
          if(el.lt.10.d0) convg=1.0d-10

          write(ko,'(10x,f10.5,10x,1p,3(2x,e8.1))') rmatch,convg,convg,
     +      convg
          if (ecis2(15:15).eq.'T') write(ko,'( )')
          write(ko,'(f5.2,2i2,a1,5f10.5)') tgts,0,1,tgtp,eninc,
     +      projspi,projmas,tarmas,ztar*zproj
          write(ko,'(2i5)') 0,0
          jj=2
          if(koptom.eq.1)jj=1
          if(npp.eq.1) jj=1
c         if(typdis(j).eq.'D')bdis(j)=bdis(j)/(rv*atar**(1./3.))
          write(ko,'(f5.2,2i2,a1,5f10.5)') sdis(j),0,jj,tgtpx,edis(j),
     +      projspi,projmas,tarmas,ztar*zproj
          write(ko,'(2i5)') 1,1
          xjdis=abs(sdis(1)-sdis(j)) + 0.001
          jdis=int(xjdis)
          if(modtyp.eq.3)jdis=langmo(j)
          write(ko,'(i5,5x,f10.5)') jdis,bdis(j)

          jchk=1
          do k=1,npp
            if (k.gt.1)go to 74
            el=eninc
C           jchk=0

            jchk=1
            go to 76
   74       jchk=0

            xratio = tarmas / (tarmas+projmas)
            el=eninc-edis(j)/xratio
            if(el.lt.0. .and. modtyp.ne.5) el=0.0001
   76       call optmod(k,kecis)
                 if(Iflag.gt.0) return
            write(ko,'(3f10.5)') v,rv,av
            write(ko,'(3f10.5)') w,rw,aw
            write(ko,'(3f10.5)') vd,rvd,avd
            write(ko,'(3f10.5)') wd,rwd,awd
            write(ko,'(3f10.5)') vso,rvso,avso
            write(ko,'(3f10.5)') wso,rwso,awso
C           COULOMB POTENTIAL
C      1-10   VAL(20)  REDUCED COULOMB RADIUS IN FERMIS.                ECIS-724
C     11-20   VAL(21)  DIFFUSENESS OF A WOODS-SAXON CHARGE DISTRIBUTION.ECIS-725
            write(ko,'(3f10.5)') rc,ac,0.
            write(ko,'(3f10.5)') 0.,0.,0.
          enddo
          write(ko,'(3f10.5)') 0.,astep,180.
        enddo
      enddo
      return
      end subroutine ecisdwip

c-----------------------------------------------------------------------------
      subroutine vib1
c-----------------------------------------------------------------------------
c
c     Routine to setup for vibrational model calculation
c
      ecis2(20:20)='T'
c
      ntar=0
      do 44 nis=1,nisotop
      if(iz(nis).ne.iztar.or.ia(nis).ne.iatar)go to 44
      ntar=nis
!     ntar>0 (Collective Levels exist )
      open(unit=45,file='cclev-omp.dat')
      write(45,'("# Reference Number =",i5," Incident Particle: ",a8,
     +  " imodel=",i1," Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),imodel,iztar,iatar,Efermi
C     write(45,'("Type potential: ",a11,2x,a27,2x,a16)')
C    +              modelz,dispz,relz
C     write(45,'("First Author: ",14a1)') (opname(i),i=1,14)
C     write(45,'("Z-Range=",i2,"-",i3,"  A-Range=",i3,"-",i3,
C    + "  E-Range=",f8.3,"-",f8.3," MeV         Ef =",f7.3)')
C    +  izmin,izmax,iamin,iamax,emin,emax,Efermi

      SELECT CASE (imodel)

        CASE (1) ! Rigid rotor parameters (imodel=1)
          write(35,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,'(A26)') '# LEVELS: rigid rotor     '
          do k=1,ncoll(ntar)
            write(45,99045) EX(k,ntar), SPIn(k,ntar), IPAr(k,ntar)
          enddo

        CASE (2) ! Vibrational rotor parameters (imodel=2)
          write(35,10) ncoll(ntar)
          write(45,10) ncoll(ntar)
          write(45,'(A26)') '# LEVELS: vibrational     '
          do k=1,nvib(ntar)
            write(45,99045) EXV(k,ntar), SPInv(k,ntar),IPArv(k,ntar),
     &                      NPH(k,ntar), DEFv(k,ntar), THEtm(k,ntar)
          enddo

        CASE (3) ! Soft rotor parameters (imodel=3)
          write(35,10) NCOll(ntar)
          write(45,10) NCOll(ntar)
          write(45,'(A26)') '# LEVELS: soft rotor      '
          ! Record 3 from OPTMAN (Hamiltonian parameters)
          write(45,99047)  SR_hw(ntar),SR_amb0(ntar),SR_amg0(ntar),
     +                     SR_gam0(ntar),SR_bet0(ntar),SR_bet4(ntar)
          ! Record 4 from OPTMAN (Hamiltonian parameters)
          write(45,99047)  SR_bb42(ntar),SR_gamg(ntar),SR_delg(ntar),
     +                     SR_bet3(ntar),SR_et0(ntar),SR_amu0(ntar)
          ! Record 5 from OPTMAN (Hamiltonian parameters)
          write(45,99047)  SR_hw0(ntar),SR_bb32(ntar),SR_gamde(ntar),
     +                     SR_dpar(ntar),SR_gshape(ntar)
          do k=1,NCOll(ntar)
           write(45,99049) EXV(k,ntar),SPInv(k,ntar),IPArv(k,ntar),
     +                     SR_ntu(k,ntar),SR_nnb(k,ntar),
     +                     SR_nng(k,ntar),SR_nno(k,ntar)
          enddo

        CASE (4) ! Rigid-soft rotor parameters (imodel=4)
          write(35,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,'(A26)') '# LEVELS: rigid-soft rotor'
          do k=1,NCOll(ntar)
            write(45,99049)  ex(k,ntar),spin(k,ntar),ipar(k,ntar),
     +               SR_ntu(k,ntar),SR_nnb(k,ntar),SR_nng(k,ntar),
     +               SR_nno(k,ntar),defv(k,ntar),defr(k,ntar)
          enddo

        CASE DEFAULT
        write(35,'(A21)') '# SPHERICAL POTENTIAL'

      END SELECT

      close(45)
      if(idr.ne.0) THEN
        write(35,
     +    '("# WARNING: V = VHF + DWV for dispersive potentials")')
      else
        write(35,'("# ")')
      endif
      write(35,'("# En(MeV)   V     RV     AV  ",
     + "    DWV     W      RW     AW  ",
     + "    VD     RVD    AVD     WD     RWD    AWD ",
     + "    VSO    RVSO   AVSO    WSO    RWSO   AWSO   VHFpn    WSpn",
     + "     DVSpn     WVpn    DVVpn   V(Lane)  W(Lane)    VHF" )')

      go to 48
  44  continue

!     ntar<=0 Collective Levels do not exist
      if(ntar.le.0) THEN
       write(35,'("# ncoll=",i2)') 0
       write(35,'(1H#/,"# En(MeV)   V     RV     AV  ",
     + "    DWV     W      RW     AW  ",
     + "    VD     RVD    AVD     WD     RWD    AWD ",
     + "    VSO    RVSO   AVSO    WSO    RWSO   AWSO   VHFpn    WSpn",
     + "     DVSpn     WVpn    DVVpn   V(Lane)  W(Lane)    VHF" )')
      endif

      ndis=0

      if(Iflag.ne.-1 .and. imodel.ne.0) then
        write(6,'(" WARNING: RIPL reference number: ",i6)') iref
        write(6,'(" WARNING: Target nucleus izatar=",i6,
     +    " does not have coupled-channel data.")') izatar
        Iflag=-1
      endif

      return
c
  48  ndis=nvib(ntar)

      do 54 k=1,ndis
      edis(k)=exv(k,ntar)
      sdis(k)=spinv(k,ntar)
      ipdis(k)=iparv(k,ntar)
  54  bdis(k)=defv(k,ntar)
      return
 10   FORMAT(2H# ,'ncoll=',i2,' lmax=',i2,' ndef=',i2,
     +       ' K=',f4.1,' betas:',8(f7.3,1x))
99045 FORMAT (f12.8,f7.1,2I4,1p,2(1x,e11.4))
99047 FORMAT (6(e11.5,1x))
99049 FORMAT (f12.8,f7.1,I4,4I2,2(F8.5,1x))
      end subroutine vib1

c-----------------------------------------------------------------------------
      subroutine dwba1(kexit)
c-----------------------------------------------------------------------------
c
c     Routine to setup for vibrational model calculation
c
      character*2 sym,idummy
      character*10 authz,audef1,audef2
c
      if(kdef.eq.1) audef1='JENDL-3.2 '
      if(kdef.eq.1) audef2='JENDL-3.2 '
      if(kdef.eq.2) audef1='ENSDF(Q)  '
      if(kdef.eq.2) audef2='ENSDF(Q)  '
      if(kdef.eq.3) audef1='ENSDF(BE2)'
      if(kdef.eq.3) audef2='ENSDF(BE3)'
      if(kdef.eq.4) audef1='Raman2    '
c     if(kdef.eq.4) audef1='Raman     '
      if(kdef.eq.4) audef2='Kibedi    '
c     if(kdef.eq.4) audef2='Spear     '
      if(kdef.eq.5) audef1='User      '
      if(kdef.eq.5) audef2='User      '
c
      ecis2(20:20)='T'
C     if(modtyp.eq.3)ecis2(42:42)='T'
      if(modtyp.eq.3)ecis1(30:30)='T'
c
      open(unit=17,file='deform.dat')
c
      klev=1
      edis(klev)=0.
      sdis(klev)=0.
      ipdis(klev)=0
      langmo(klev)=0
      bdis(klev)=0.
c
c     Position deform.dat file
      do 25 i=1,4
  25  read(17,'(a2)') idummy
c
  30  read(17,'(2i4,1x,a2,1x,f10.6,1x,f4.1,i3,i2,1x,f10.6,2x,a10)',
     + end=900) kz,ka,sym,edisz,sdisz,ipdisz,langmoz,bdisz,authz
      kza=1000*kz+ka
      if(kza.gt.izatar)go to 900
      if(kza.ne.izatar)go to 30
      if(kdef.eq.5)go to 40
      if(authz.ne.audef1.and.authz.ne.audef2) go to 30
c
  40  klev=klev+1
      edis(klev)=edisz
      sdis(klev)=sdisz
      ipdis(klev)=ipdisz
      langmo(klev)=langmoz
      bdis(klev)=bdisz
      go to 30
c
 900  ntar=0
      ndis=klev
      if(ndis.eq.1) go to 990
      ndis=min0(ndisx,ndis)
      go to 999
c
 990  write(6,'("DWBA deformations not found for izatar=",i7)') izatar
      kexit=1
 999  close (unit=17)
      return
      end subroutine dwba1

c-----------------------------------------------------------------------------
      subroutine gaussff(igauflg)
c-----------------------------------------------------------------------------
c
c     Routine to see if Gaussian form factor required for WD(el)
c
      igauflg=0
      jab=iabs(jrange(4))
      if(jrange(4).lt.1)go to 40
c
      jp=1
      do 30 j=1,jab
      if(el.gt.epot(4,j)) jp=j+1
  30  continue
      j=min0(jp,jab)
      if(rco(4,j,1).lt.0.)igauflg=1
c
  40  return
      end subroutine gaussff

c-----------------------------------------------------------------------------
      subroutine optmod(k,kecis)
c-----------------------------------------------------------------------------
c
c     Routine to generate input for ECIS from RIPL library
c
      real*8 rlib(6),alib(6),vlib(6),b(6,ndim1,15)
      real*8 OPTrlib(6),OPTalib(6),OPTvlib(6)
c
c     To use dispersive optical model package
c
c             variables
      real*8 As,Bs,Cs,AAv,Bv,AAvso,Bvso,EEE,Ep,Ea,Ef,VCshift
      real*8 alpha_PB,beta_PB,gamma_PB,Vnonl,VVcoul,VScoul,VDcoul
      real*8 DWS,DWV,DWVnonl,DWVso,DerDWV,DerDWS,dtmp,dtmp1
!MS$IF DEFINED (alone)
      real*8 dtmp2
      real*8 rvolint,wvolint,rsurint,wsurint,DWVintNONL,DWVint
      real*8 vcoldisp, vcolint, scolint, rvolintsph, wsoint, vsoint
!MS$ENDIF
      real*8 VHFnum,DWVcor,VSOnum,DSOcor,ALhf,ALso,DWScor,AHF
      real*8 WDE,WVE,Visov,WVisov,WSisov,AlphaV,ALpha,RSTep

      real*8 betass(ndim5)
      integer n,iq,nns,iqm,nnv,nnl,i
      integer iaref
      integer k,kecis

      common /energy/EEE,Ef,Ep,Ea
      common /Wenerg/WDE,WVE
      common /pdatas/As,Bs,Cs,nns,iq
      common /pdatav/AAv,Bv,nnv
      common /pdatao/AAvso,Bvso,nnl
      common /dispcor/VHFnum,DWVcor,VSOnum,DSOcor,ALhf,ALso,
     +                DWScor,ALpha,RSTep

      nnv = 0
      nns = 0
      nnl = 0
      VHFnum = 0.d0
      VSOnum = 0.d0
      DWVcor = 0.d0
      DWScor = 0.d0
      DSOcor = 0.d0
      VDcoul = 0.d0
      AHF  = 0.d0
      ALhf = 0.d0
      ALso = 0.d0
      Alpha = 0.d0

      OPTalib = 0.d0
      OPTrlib = 0.d0
      OPTvlib = 0.d0


      RStep = 10.d0
c
c     Generate optical model parameters for ECIS
c
      pi=acos(-1.d0)
      rc=0.d0
      ac=1.d0

      IDRs = 0 ! Imaginary and HF geometry in dispersive potentials are the same
      if(jcoul.lt.1)go to 194
      jc=1
      do 190 j=1,jcoul
      if(el.gt.ecoul(j)) jc=j+1
 190  continue
      jc=min0(jc,jcoul)
      rc=rcoul0(jc)*atar**(-1./3.) + rcoul(jc) +
     +   rcoul1(jc)*atar**(-2./3.) + rcoul2(jc)*atar**(-5./3.) +
c        RCN addition to consider new Morillon-Romain potential
     +   rcoul3(jc)
c    +   rcoul3(jc)*atar

      if(beta(jc).gt.0.d0) then
        Iflag = 6
        write(6,*)
     +'ERROR: ECIS can not deal with non-local potentials, use SCAT2000'
            return
      endif
      ac=acoul(jc)
 194  encoul2=0.
      if(rc.gt.0.) encoul2=1.73*ztar/(rc*atar**(1./3.))
c
      do i=1,6
c
c     For Lane consistent potentials change the incident energy for proton potentials
c     by the specified Coulomb shift (pot(1,1,25)
c
      VCshift = 0.d0
      IF(izproj.eq.1 .and. pot(i,1,25).gt.0.d0)
     >   VCshift = pot(i,1,25)*ztar/atar**(1.d0/3.d0)
      vc = 0.d0
      DerDWV = 0.d0
      DerDWS = 0.d0
      DWVnonl = 0.d0
      VVcoul = 0.d0
      VScoul = 0.d0
      AlphaV = 0.d0
      rlib(i)=0.d0
      alib(i)=0.d0
      vlib(i)=0.d0
      jab=iabs(jrange(i))
      if(jrange(i).lt.1)go to 300
      jp=1
      do 204 j=1,jab
      if(el.gt.epot(i,j)) jp=j+1
 204  continue
      j=min0(jp,jab)
      Ef=efermi
      if(pot(i,j,18).ne.0.) Ef=pot(i,j,18) + pot(i,j,19)*atar

      elf = el - Ef - VCshift

      iaref = 0
C     Reference nucleus defined
      if (NINT(pot(i,j,24)).eq.1 .and. pot(i,j,23).gt.0.)
     *  iaref = NINT(pot(i,j,23))
c
c     Calculate radius and diffuseness parameters
      if(rco(i,j,13).eq.0.) then
        rlib(i)=abs(rco(i,j,1)) + rco(i,j,3)*eta
     *       + rco(i,j,4)/atar + rco(i,j,5)/sqrt(atar)
     *       + rco(i,j,6)*atar**(2./3.) + rco(i,j,7)*atar
     *       + rco(i,j,8)*atar**2  + rco(i,j,9)*atar**3
     *       + rco(i,j,10)*atar**(1./3.)
     *       + rco(i,j,11)*atar**(-1./3.)
C--------------------------------------------------------------------
C     RCN, 08/2004, to handle new extension to the OMP RIPL-2 format
C    *       + rco(i,j,2)*(el-VCshift)
C    *       + rco(i,j,12)*(el-VCshift)*(el-VCshift)
     *       + rco(i,j,2)*el + rco(i,j,12)*el*el
      else
C     RCN, 09/2004, to handle new extension to the OMP RIPL-2 format
        nn = int(rco(i,j,7))
        rlib(i)= ( abs(rco(i,j,1)) + rco(i,j,2)*atar ) *
     *           ( 1.d0 - ( rco(i,j,3) + rco(i,j,4)*atar ) * elf**nn/
     *           ( elf**nn + ( rco(i,j,5) + rco(i,j,6)*atar )**nn ) )
      endif
      alib(i)=abs(aco(i,j,1)) + aco(i,j,2)*el + aco(i,j,3)*eta
     *        + aco(i,j,4)/atar + aco(i,j,5)/sqrt(atar)
     *        + aco(i,j,6)*atar**(2./3.) + aco(i,j,7)*atar
     *        + aco(i,j,8)*atar**2 + aco(i,j,9)*atar**3
     *        + aco(i,j,10)*atar**(1./3.) + aco(i,j,11)*atar**(-1./3.)

C     Dependence calculated inside OPTMAN
C     if(iaref.gt.0 .and. rco(i,j,13).eq.0.) THEN
      if(iaref.eq.0) THEN
        OPTrlib(i) =  rlib(i)
        OPTalib(i) =  alib(i)
        if (i.eq.2 .and. idr.ge.2) THEN
          IF( ABS(rlib(1)-rlib(2)).gt.0.0001 ) IDRs = 1
          IF( ABS(alib(1)-alib(2)).gt.0.0001 ) IDRs = 1
        endif
      else
        OPTrlib(i) =  rco(i,j,1) + rco(i,j,7)*(atar-iaref)
        OPTalib(i) =  aco(i,j,1) + aco(i,j,7)*(atar-iaref)
        if (i.eq.2 .and. idr.ge.2) then
          IF( ABS(OPTrlib(1)-OPTrlib(2)).gt.0.0001 ) IDRs = 1
          IF( ABS(OPTalib(1)-OPTalib(2)).gt.0.0001 ) IDRs = 1
        endif
      endif

      if(alib(i).gt.0.d0) RStep = min( 0.25d0*alib(i), RStep)
c
      if (pot(i,j,24).eq.0.) go to 210
c
c     Special Koning-type potential formulas
c
      if (pot(i,j,24).eq.1. .or. pot(i,j,24).eq.3.) then
c
c       Koning-type formulas
c
        elf = el - Ef - VCshift
        if(i.eq.1) call bcoget(b,j,Visov,WVisov,WSisov)

        if(i.eq.1 .and. b(1,j,5).ne.0.d0) then
          vc = b(1,j,1)*encoul2*( b(1,j,2) - 2.*b(1,j,3)*elf +
     +    3.*b(1,j,4)*elf**2 + b(i,j,14)*b(i,j,13)*exp(-b(i,j,14)*elf) )
          VDcoul = b(i,j,5)*vc
        endif

        nn = int(pot(i,j,13))
c       Retrieving average energy of the particle states Ep
        Ep=Ef
        if( (i.eq.2) .or. (i.eq.4) ) Ep=pot(i,j,20)
        if(Ep.eq.0.) Ep=Ef
        elf = el - Ep - VCshift

        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))

        if(i.eq.1) then
          if(b(i,j,1).ne.0.d0) then
            ALhf = b(i,j,14)
            AHF  = b(i,j,1)*b(i,j,13)
          else
            ALhf = b(i,j,12)
            AHF  = b(i,j,11)
          endif
        endif
        if(i.eq.5) ALso = b(i,j,12)

        vlib(i)=
     +    b(i,j,1)*( b(i,j,15) - b(i,j,2)*elf + b(i,j,3)*elf**2 -
     +    b(i,j,4)*elf**3 + b(i,j,13)*exp(-b(i,j,14)*elf) ) +
     +    b(i,j,5)*vc + b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +    b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +    (elf**nn + b(i,j,10)**nn)) +
     +    b(i,j,11)*exp(-b(i,j,12)*elf)
      endif

      if (pot(i,j,24).eq.2.) then
c
c       Morillon-Romain formulas
c
        elf = el - Ef - VCshift
        if(i.eq.1) call bcoget(b,j,Visov,WVisov,WSisov)
c
c       Vhf(E) calculated from nonlocal approximation
c          as suggested by Perey and Buck
        alpha_PB = b(i,j,1)
        beta_PB  = b(i,j,2)
        gamma_PB = b(i,j,3)
        EEE = el  - VCshift
        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))

        Vnonl = 0.d0
        if(i.eq.1 .or. i.eq.5) then
          Vnonl = -Vhf(EEE,alpha_PB,beta_PB,gamma_PB)
          vc = 0.d0
          if(i.eq.1 .and. b(1,j,5).ne.0.d0) then
C           MR do not use derivarive of the potential
C
C           Numerical derivative of the Vhf
C           Vnonlm = -Vhf(EEE-0.05,alpha_PB,beta_PB,gamma_PB)
C           Vnonlp = -Vhf(EEE+0.05,alpha_PB,beta_PB,gamma_PB)
C           Coulomb correction for Hartree-Fock potential
C           vc = encoul2*(Vnonlm-Vnonlp)*10.d0
C
C           MR are using constant Coulomb correction
            vc = encoul2
            VDcoul = b(i,j,5)*vc
          endif
        endif
        nn = int(pot(i,j,13))

        vlib(i)=
     +    Vnonl + b(i,j,5)*vc +
     +    b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +    b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +    (elf**nn + b(i,j,10)**nn)) +
     +    b(i,j,11)*exp(-b(i,j,12)*elf)
      endif
c
c     Nonlocality consideration
c
c     Retrieving energy above which nonlocality in the volume absorptive
c               potential is considered (Ea)
c
      Ea=pot(i,j,21)
      if(Ea.eq.0.) Ea=1000.1d0
      if(i.eq.2 .and. Ea.lt.1000.d0 ) then
        AlphaV=pot(i,j,22)
        if(AlphaV.eq.0.d0) AlphaV=1.65d0
        if(el.gt.(Ef+Ea))  vlib(i) = vlib(i) +
     +   AlphaV*(sqrt(el)+(Ef+Ea)**1.5d0/(2.d0*el)-1.5d0*sqrt(Ef+Ea))
      endif

      go to 300
c
 210  if (pot(i,j,23).eq.0.) go to 220
c
c     Special Varner-type potential formulas
c
      elf = el - VCshift
c
      vlib(i)= (pot(i,j,1) + pot(i,j,2)*eta)/
     +    (1.+ exp((pot(i,j,3) - elf + pot(i,j,4)*encoul2)/pot(i,j,5)))
      if(pot(i,j,6).eq.0.)go to 300
      vlib(i) = vlib(i)
     +    + pot(i,j,6)*exp((pot(i,j,7)*elf - pot(i,j,8))/pot(i,j,6))
      go to 300
c
 220  if (pot(i,j,22).eq.0.) go to 230
c
c     Special Smith-type potential formulas
c
      elf = el - VCshift
c
      vlib(i)=pot(i,j,1) + pot(i,j,2)*eta
     +    + pot(i,j,6)*exp(pot(i,j,7)*elf + pot(i,j,8)*elf*elf)
     +    + pot(i,j,9)*elf*exp(pot(i,j,10)*elf**pot(i,j,11))
      if(pot(i,j,5).ne.0.)vlib(i)=vlib(i)
     +    + pot(i,j,3)*cos(2.*pi*(atar - pot(i,j,4))/pot(i,j,5))
      go to 300
c
c     Standard potential formulas
c
 230  elf = el - VCshift
      vlib(i)=pot(i,j,1) + pot(i,j,7)*eta + pot(i,j,8)*encoul
     +   + pot(i,j,9)*atar + pot(i,j,10)*atar**(1./3.)
     +   + pot(i,j,11)*atar**(-2./3.) + pot(i,j,12)*encoul2
      if(elf.gt.0.) vlib(i) = vlib(i) +
     +      (pot(i,j,2) + pot(i,j,13)*eta + pot(i,j,14)*atar)*elf
     +   + pot(i,j,3)*elf**2 + pot(i,j,4)*elf**3 + pot(i,j,6)*sqrt(elf)
     +   + pot(i,j,17)*encoul/elf**2
     +   + (pot(i,j,5) + pot(i,j,15)*eta + pot(i,j,16)*elf)*log(elf)


 300  if(i.eq.1) VHFnum = vlib(1)
      if(i.eq.2) ALPha  = AlphaV
      if(i.eq.5) VSOnum = vlib(5)

      ENDDO
c
c     To calculate dispersion relation contribution
c
 152  if(abs(idr).ge.2) then
c
c       Exact calculation of the dispersive contribution
c
        i=2
c       Only one energy range
        j=1
c
c       Only shifted if the real potential is shifted
c
        VCshift = 0.d0
        IF(izproj.eq.1 .and. pot(1,1,25).gt.0.d0)
     >   VCshift = pot(1,1,25)*ztar/atar**(1.d0/3.d0)
c
        EEE = el - VCshift
c       Real volume contribution from Dispersive relation
        DWV=0.d0
c
        if(jrange(i).gt.0 .and. pot(2,1,24).ne.0) then
          AAv=b(i,j,6)
          Bv =b(i,j,7)
          n = nint( pot(i,j,13) )
          nnv = n

          if(n.eq.0 .or. mod(n,2).eq.1) stop
     +      'ERROR: Zero or odd exponent in Wv(E) for dispersive OMP'
          Ep=pot(i,j,20)
          if(Ep.eq.0.) Ep=Ef

          Ea=pot(i,j,21)
          if(Ea.eq.0.) Ea=1000.1d0

c         if(modtyp.eq.6) goto 154
c         Analytical DOM integral
          DWV=DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE,n,DerDWV)

C         Numerical DOM integral (used in RIPL-2 released interface)
C         WVE=WVf(AAv,Bv,Ep,Ef,EEE,n)
C         ftmp1=2*DOM_int(Delta_WV,WVf,Ef,Ef+5.*Bv,150000.d0,EEE,0.d0)

C         Coulomb correction for real volume potential
          if (pot(1,1,25).ne.0) DerDWV = 0.d0
          DerDWV = -b(1,1,5)*encoul2*DerDWV
C         numerical DOM derivative (not needed for a time being)
C         DWVp = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE+0.1d0,n,dtmp)

C         DWVm = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE-0.1d0,n,dtmp)
C         DerDWV = -b(1,1,5)*encoul2*(DWVp-DWVm)*5.d0
c         if(idr.le.-2) then
c           numerical DOM integral (not needed for a time being)
c           WVE=WVf(AAv,Bv,Ep,Ef,EEE,n)
c           DWV=2*DOM_int(Delta_WV,WVf,Ef,Ef+5.*Bv,150000.d0,EEE,0.d0)
c         endif
c
          T12der = 0.d0
          if(Ea.lt.1000.d0) THEN
             AlphaV=pot(i,j,22)
             if(AlphaV.eq.0.d0) AlphaV=1.65d0

             Dwplus = AlphaV*DOM_INT_T2(Ef,Ea,EEE)
             dtmp1 = Wvf(AAv,Bv,Ep,Ef,Ef+Ea,n)
             Dwmin = dtmp1*DOM_INT_T1(Ef,Ea,EEE)

             DWV = DWV + Dwplus + Dwmin
             DWVnonl = Dwplus + Dwmin
c            Coulomb correction for nonlocal dispersive contribution
c                to real volume potential
             if(b(1,1,5).ne.0.d0 .and. pot(1,1,25).eq.0) then
               if(eee.ne.0.05d0) then
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.05d0)
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.05d0)
                     T2der = AlphaV*(T2p-T2m)*10.d0
                 T1p = DOM_INT_T1(Ef,Ea,EEE+0.05d0)
                 T1m = DOM_INT_T1(Ef,Ea,EEE-0.05d0)
                 T1der = dtmp1*(T1p-T1m)*10.d0
                 T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
               else
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.1d0)
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.1d0)
                     T2der = AlphaV*(T2p-T2m)*5.d0
                 T1p = DOM_INT_T1(Ef,Ea,EEE+0.1d0)
                 T1m = DOM_INT_T1(Ef,Ea,EEE-0.1d0)
                 T1der = dtmp1*(T1p-T1m)*5.d0
                 T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
               endif
             endif
          endif
          VVcoul = DerDWV + T12der

        endif

154     i=4
c       Only one energy range
        j=1
c
c
c       Only shifted if the real potential is shifted
c
        VCshift = 0.d0
        IF(izproj.eq.1 .and. pot(1,1,25).gt.0.d0)
     >   VCshift = pot(1,1,25)*ztar/atar**(1.d0/3.d0)
c
        EEE = el - VCshift
c       Real surface contribution from Dispersive relation
        DWS=0.d0
        if(jrange(i).gt.0 .and. pot(4,1,24).ne.0) then
          As=b(i,j,8)
          Bs=b(i,j,10)
          Cs=b(i,j,9)
          n = nint( pot(i,j,13) )
          nns = n
          if(n.eq.0 .or. mod(n,2).eq.1) stop
     +      'ERROR: Zero or odd exponent in Wd(E) for dispersive OMP'
          iq=1
          if(b(4,j,12).gt.0.) iq=nint(b(4,j,12))

c         Retrieving average energy of the particle states Ep
          Ep=pot(i,j,20)
          if(Ep.eq.0.) Ep=Ef

C         if(modtyp.eq.6) goto 156

          if(idr.ge.2) then
c           analytical DOM integral
            DWS = DOM_INT_Ws(Ef,Ep,As,Bs,Cs,EEE,n,DerDWS)
c           Coulomb correction for real surface potential
            if (pot(1,1,25).ne.0) DerDWS = 0.d0
            VScoul = -b(1,1,5)*encoul2*DerDWS
          endif

          if(idr.le.-2) then
c           numerical DOM integral
            nns=n
            WDE=WDf(As,Bs,Cs,Ep,EEE,n,iq)
            DWS = 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE,WDE)
c           Coulomb correction for real surface potential
            if(b(1,1,5).ne.0.d0 .and. pot(1,1,24).eq.1) then
              WDE=WDf(As,Bs,Cs,Ep,EEE+0.1d0,n,iq)
              DWSp =
     >         2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE+0.1d0,WDE)
              WDE=WDf(As,Bs,Cs,Ep,EEE-0.1d0,n,iq)
              DWSm =
     >         2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE-0.1d0,WDE)
c             Numerical derivative
              DerDWS = (DWSp-DWSm)*5.d0
              if (pot(1,1,25).ne.0) DerDWS = 0.d0
              VScoul = -b(1,1,5)*encoul2*DerDWS
            endif
          endif

        endif

156     i=6
c       Only one energy range
        j=1
c       Real spin orbit contribution from Dispersive relation
        DWVso=0.d0
c
        if(jrange(i).gt.0.and.pot(6,1,24).ne.0.and.abs(idr).eq.3) then
          AAvso = b(i,j,6)
          Bvso  = b(i,j,7)
          nnl = nint( pot(i,j,13) )
          if(nnl.eq.0 .or. mod(nnl,2).eq.1) stop
     +      'ERROR: Zero or odd exponent in Wso(E) for dispersive OMP'

c         analytical DOM integral
          DWVso=DOM_INT_Wv(Ef,Ef,AAvso,Bvso,EEE,nnl,dtmp)

        endif

c       Adding real volume dispersive contribution to the real potential
c       Geometry parameters are the same as for the volume potential(imag and real).

        if(idr.ne.0) then
C         if (IDRs.gt.0) DWVcor = DWV + VVcoul
          DWVcor = DWV + VVcoul
          DWScor = DWS + VScoul
          DSOcor = DWVso
        endif
C       Dependence calculated inside OPTMAN
        iaref = 0
C       Reference nucleus defined
        if (NINT(pot(1,1,24)).eq.1 .and. pot(1,1,23).gt.0.)
     *    iaref = NINT(pot(1,1,23))

        if(iaref.gt.0) then
          OPTvlib(1) =  vlib(1) + pot(1,1,22)*(atar-iaref) + DWVcor
        else
          OPTvlib(1) =  vlib(1) + DWVcor
        endif
        vlib(1)= vlib(1) + DWVcor

      OPTvlib(2) = vlib(2)
      OPTvlib(4) = vlib(4)
c       Including real surface and Coulomb dispersive contribution
c       Geometry parameters are the same as for the imaginary surface potential.
        vlib(3)= DWScor
        alib(3)= alib(4)
        rlib(3)= rlib(4)
        OPTvlib(3)= DWScor
        OPTalib(3)= OPTalib(4)
        OPTrlib(3)= OPTrlib(4)
c       Adding real spin orbit dispersive contribution to the real spin orbit potential
c       Geometry parameters are the same as for the imaginary spin orbit potential(imag and real).
        vlib(5) = vlib(5) + DSOcor
        OPTvlib(5) = vlib(5)
      OPTvlib(6) = vlib(6)

      endif

      RStep = max( RStep, 0.1d0)
c
c     Factor coming from Dirac equation reduction (relativistic approximation).
c     Using this factor force us to employ relativistic kinematics as well.
c
      gamma=1.d0
      if(irel.eq.2) then
C
C       Follwing Madland private comm. 28/09/2008
C       it should be applied to all potentials when used
C       including spin-orbit
C
c       Target system mass in MeV
        EMtar=tarmas*amu0c2
c       Total system mass in MeV
        EMtot=(tarmas+projmas)*amu0c2
c       Total kinetic energy in cm
        Tcm=sqrt(2*EMtar*el + EMtot**2) - EMtot
c       Relativistic correction to the potential (non relativistic target!!).
        gamma=1.d0+Tcm/(Tcm+2*projmas*amu0c2)
      endif
c
      v=vlib(1)*gamma
      OPTvlib(1) = OPTvlib(1)*gamma
      rv=rlib(1)
      av=alib(1)

      w =vlib(2)*gamma
      OPTvlib(2) = OPTvlib(2)*gamma
      if(w.lt.0.0.and.nonegw.eq.1) OPTvlib(2) = 0.d0
      if(w.lt.0.0.and.nonegw.eq.1) w=0.
      rw=rlib(2)
      aw=alib(2)

      vd=vlib(3)*gamma
      OPTvlib(3) = OPTvlib(3)*gamma
      rvd=rlib(3)
      avd=alib(3)

      wd=vlib(4)*gamma
      OPTvlib(4) = OPTvlib(4)*gamma
      if(w.lt.0.0.and.nonegw.eq.1) OPTvlib(4) = 0.d0
      if(wd.lt.0.0.and.nonegw.eq.1) wd=0.
      rwd=rlib(4)
      awd=alib(4)

      vso=vlib(5)*gamma
      OPTvlib(5) = OPTvlib(5)*gamma
      rvso=rlib(5)
      avso=alib(5)

      wso=vlib(6)*gamma
      OPTvlib(6) = OPTvlib(6)*gamma
      rwso=rlib(6)
      awso=alib(6)


      IF(k.eq.1 .and. kecis.eq.1) then
C       IF(IDRs>0)
C    >    write(*,*)
C    >    'Dispersive potential with Real geom.<> Imag. geom. !!'
C       write(*,*) 'GS Dispers. OMP parameters'
C       write(*,*) 'E=',sngl(El),' Ef=',sngl(Ef),' Ep=',sngl(Ep)
C     write(*,*) 'Real HF volume  VHF,rHF,aHF:',
C    >    sngl(VHFnum*gamma),sngl(rlib(1)),sngl(alib(1))
C     write(*,*) 'Disper. volume  DWV,rWV,aWV:',
C    >    sngl(DWV*gamma),sngl(rlib(2)),sngl(alib(2))
C       write(*,*) '         DWV(sym), DW+, DW-:',
C    >    sngl((DWV-(Dwplus + Dwmin))*gamma),
C    >    sngl(Dwplus*gamma),sngl(Dwmin*gamma)
C       write(*,*) 'Imag volume  wv,rwv,awv:',
C    >    sngl(max(0.d0,vlib(2)*gamma)),sngl(rlib(2)),sngl(alib(2))
c       Real surface geometry parameters are the same as for the imaginary surface potential.
C     write(*,*) 'Real surface vd,rvd,avd:',
C    >    sngl(vlib(3)*gamma),sngl(rlib(3)),sngl(alib(3))
C     write(*,*) 'Imag surface wd,rwd,awd:',
C    >    sngl(max(0.d0,vlib(4)*gamma)),sngl(rlib(4)),sngl(alib(4))
C     write(*,*) 'Real SO pot vso,rvs,avs:',
C    >    sngl(vlib(5)*gamma),sngl(rlib(5)),sngl(alib(5))
C     write(*,*) 'Imag SO pot wso,rws,aws:',
C    >    sngl(vlib(6)*gamma),sngl(rlib(6)),sngl(alib(6))
C       write(*,*)
      ENDIF

C     if (modtyp.ge.4) return

      if(jchk.eq.1) then
       iqm = 1
       do j=1,ndim5
         betass(j)=0.d0
       enddo
       if(imodel.gt.0 .and. ntar.gt.0) then
         iqm=idef(ntar)
         do j=2,iqm,2
           betass(j)=dble(def(ntar,j))
         enddo
       endif

!MS$IF DEFINED (alone)
C
C      DEFORMED VOLUME INTEGRALS
C
C      write(*,'(1x,A3,F7.3,3x,A5,5(f5.3,1x))')
C    >    ' E=',el,' DEF:',(betass(j),j=2,iqm,2)
       call VOLIN(1,dble(VHFnum*gamma) , dble(rv), !dble(rv *atar**(1./3.)),
     >             dble(av),rvolintsph)
       ! rvolintsph = rvolintsph  / atar / iaproj
       rvolintsph = rvolintsph     / iaproj

       rvolint = VolumeInt2D(0,BETAss,iqm,dble(Vhfnum*gamma),
     >                     dble(rv *atar**(1./3.)),dble(av),1)
       rvolint = rvolint / atar     / iaproj

       ! write(6,*) ' VOL INTEGRAL =',rvolint,' SPH=',rvolintsph

       DWVint = VolumeInt2D(0,BETAss,iqm,dble(DWV*gamma),
     >                     dble(rw *atar**(1./3.)),dble(aw),1)
       DWVint = DWVint  / atar / iaproj
       rvolint = rvolint + DWVint

       vcoldisp = 0.d0
       if(VVcoul.gt.0) then
         vc = VVcoul*gamma
         vcoldisp = VolumeInt2D(0,BETAss,iqm,dble(vc),
     >                     dble(rw *atar**(1./3.)),dble(aw),1)
         vcoldisp = vcoldisp / atar / iaproj
       endif
       rvolint = rvolint + vcoldisp
C
C      Spin-orbit volume integrals
C
       vsoint = VolumeInt2D(0,BETAss,iqm,dble(vso*gamma),
     >                     dble(rvso *atar**(1./3.)),dble(avso),3)
       vsoint = vsoint / atar / iaproj * cso
       wsoint = VolumeInt2D(0,BETAss,iqm,dble(wso*gamma),
     >                     dble(rwso *atar**(1./3.)),dble(awso),3)
       wsoint = wsoint / atar / iaproj * cso
       DWVintNONL = VolumeInt2D(0,BETAss,iqm,dble(DWVnonl*gamma),
     >                     dble(rw *atar**(1./3.)),dble(aw),1)
       DWVintNONL = DWVintNONL  / atar / iaproj
       wvolint = VolumeInt2D(0,BETAss,iqm,dble(w),
     >                     dble(rw *atar**(1./3.)),dble(aw),1)
       wvolint = wvolint  / atar / iaproj
       rsurint = VolumeInt2D(0,BETAss,iqm,dble(vd),
     >                     dble(rvd *atar**(1./3.)),dble(avd),2)
       rsurint = rsurint / atar     / iaproj
       wsurint = VolumeInt2D(0,BETAss,iqm,dble(wd),
     >                     dble(rwd *atar**(1./3.)),dble(awd),2)
       wsurint = wsurint / atar     / iaproj
       vcolint = 0.d0
       scolint = 0.d0
       if((VDcoul+VScoul).gt.0) then
         vc = VDcoul*gamma
         vcolint = VolumeInt2D(0,BETAss,iqm,dble(vc),
     >                     dble(rv *atar**(1./3.)),dble(av),1)
         vcolint = vcolint / atar / iaproj
         vcolint = vcolint + vcoldisp
         vc = VScoul*gamma
         scolint = VolumeInt2D(0,BETAss,iqm,dble(vc),
     >                     dble(rvd *atar**(1./3.)),dble(avd),2)
         scolint = scolint / atar / iaproj
       endif
!MS$ENDIF

C
C       Printing potential values
C
C      if(modtyp.gt.2) return

        if(AS.ne.0.d0 .and. AHF.ne.0.d0 .and. AAV.ne.0.d0) then
C         write(35,'(f7.3,f9.3,f6.3,f6.3,f7.3,
C    +             6(f7.3,f6.3,f6.3),23(1x,f7.3,1x))')
C    +    el,v,rv,av,DWVcor,w,rw,aw,vd,rvd,avd,wd,rwd,awd,
C    +    vso,rvso,avso,wso,rwso,awso,
C    +    2*sqrt(eta/atar)*Visov *VHFnum/AHF,
C    +    2*sqrt(eta/atar)*WSisov*wd/AS,
C    +    2*sqrt(eta/atar)*WSisov*DWScor/AS,
C    +    2*sqrt(eta/atar)*WVisov*w/AAV,
C    +    2*sqrt(eta/atar)*WVisov*DWVcor/AAV,
C    +    Visov *VHFnum/AHF, WSisov*wd/AS, el
C         write(35,'(f7.3,f9.3,f6.3,f6.3,f7.3,
C    +             6(f7.3,f6.3,f6.3),23(1x,f7.3,1x))') el,
          write(35,'(f7.3,f8.3,f7.4,f7.4,f8.4,
     +             5(f8.4,f7.4,f7.4),
     +             26(1x,f7.3,1x))') el,
     +    OPTvlib(1),OPTrlib(1),OPTalib(1),DWVcor,
     +    OPTvlib(2),OPTrlib(2),OPTalib(2),
     +    OPTvlib(3),OPTrlib(3),OPTalib(3),
     +    OPTvlib(4),OPTrlib(4),OPTalib(4),
     +    OPTvlib(5),OPTrlib(5),OPTalib(5),
     +    OPTvlib(6),OPTrlib(6),OPTalib(6),
     +    2*sqrt(eta/atar)*Visov *VHFnum/AHF,
     +    2*sqrt(eta/atar)*WSisov*wd/AS,
     +    2*sqrt(eta/atar)*WSisov*DWScor/AS,
     +    2*sqrt(eta/atar)*WVisov*w/AAV,
     +    2*sqrt(eta/atar)*WVisov*DWVcor/AAV,
     +    Visov *VHFnum/AHF, WSisov*wd/AS, VHFnum

        else
C         write(35,'(f7.3,f9.3,f6.3,f6.3,f7.3,
C    +             6(f7.3,f6.3,f6.3),18(1x,f6.3,1x))')
          write(35,'(f7.3,f8.3,f7.4,f7.4,f8.4,
     +             5(f8.4,f7.4,f7.4),
     +             26(1x,f7.3,1x))') el,
     +    v,rv,av,DWVcor,w,rw,aw,vd,rvd,avd,wd,rwd,awd,
     +    vso,rvso,avso,wso,rwso,awso
        endif
!MS$IF DEFINED (alone)
        if( idr.ne.0) then
          write(36,'(f7.3,15(1x,f8.3,1x))')
     +    el,DWV,DWV-Dwplus-Dwmin,Dwplus,Dwmin,Dwvso,
     +    OPTvlib(1),Vhfnum,DWVcor,VDCoul,DerDWV,T12der,VScoul,DWS,
     +    OPTvlib(3)

          dtmp = w
          if(dtmp.eq.0.d0) dtmp=1.d0
          dtmp1 = wd
          if(dtmp1.eq.0.d0) dtmp1=1.d0
          dtmp2 = wso
          if(dtmp2.eq.0.d0) dtmp2=1.d0
          write(34,'(f7.3,1x,12(G13.7,1x))')
     +      el,DWVcor/dtmp,(DWV-Dwplus-Dwmin)/dtmp,Dwplus/dtmp,
     +      Dwmin/dtmp,OPTvlib(2),DWS/dtmp1,OPTvlib(4),Dwvso/dtmp2,
     +      OPTvlib(6)
        endif
        write(37,'(f7.3,2x,3(f8.1,2x),8(f8.3,2x))')
     +    el,rvolint+rsurint,rvolint,rsurint,vsoint,
     +    DWVint,DWVintNONL,DWVint-DWVintNONL,
     +    wvolint+wsurint,wvolint,wsurint,wsoint
        if(izproj.gt.0) then
          write(38,'(f7.3,2x,3(f8.1,2x),1x,4(f8.3,3x))')
     +    el,rvolint+rsurint,rvolint,rsurint,
     +    vcolint+scolint,vcolint,scolint,vcoldisp
        endif
!MS$ENDIF
      endif

      return
      end subroutine optmod

c-----------------------------------------------------------------------------
      subroutine bcoget(b,j,Visov,WVisov,WSisov)
c-----------------------------------------------------------------------------
c
c     Routine to compute b coefficients for Koning global potential
c
c     Modified by R.Capote to extend RIPL format
c
c     Version of Aug. 20, 2008
c
      real*8 b(6,ndim1,15),Visov,WVisov,WSisov
c
       Visov  = 0.d0
       WVisov = 0.d0
       WSisov = 0.d0

       call setr(0.d0,b,90*ndim1)
C      Default Coulomb correction = 1.73 Z/(rc*A**(1/3)) = Ecoul2
       b(1,j,5)  =  pot(1,j,9)

C      Original Koning dependence
       b(1,j,1)  =  pot(1,j,1) + pot(1,j,2)*atar + pot(1,j,8)*eta
C      Visov = 2*pot(1,j,8)*sqrt(eta/atar)
       Visov =   pot(1,j,8)

       if(pot(1,j,24).lt.3) then

c        Soukhovitski dependence
         if((pot(1,j,20).ne.0.) .and.
     +    (pot(1,j,14) + pot(1,j,15)*atar + pot(1,j,16)).ne.0. ) then
           b(1,j,1)  =  pot(1,j,1) + pot(1,j,2)*atar + pot(1,j,8)*eta +
     +                pot(1,j,20)*eta/
     +               (pot(1,j,14) + pot(1,j,15)*atar + pot(1,j,16))
c          Visov = 2 * ( pot(1,j,8) + pot(1,j,20) ) * sqrt(eta/atar)
C          Visov =       pot(1,j,8) +
C    +                   pot(1,j,20)/
C    +                  (pot(1,j,14) + pot(1,j,15)*atar + pot(1,j,16))

           Visov =       pot(1,j,20)
         endif

         b(1,j,2)  =  pot(1,j,3) + pot(1,j,4)*atar
         b(1,j,3)  =  pot(1,j,5) + pot(1,j,6)*atar
         b(1,j,4)  =  pot(1,j,7)
         b(1,j,5)  =  pot(1,j,9)
         b(1,j,11) =  pot(1,j,10) + pot(1,j,11)*atar
         b(1,j,12) =  pot(1,j,12)

c        b coefficients from 13 to 15 added for Soukhovitski potential
c        V^DISP_R
C        b(1,j,13)  =  pot(1,j,16)
         b(1,j,13)  =  pot(1,j,16) + pot(1,j,2)*atar

c        Lambda_R
         b(1,j,14)  =  pot(1,j,17)

c        V^0_R + V^A_R*(A-232)
         b(1,j,15)  =  pot(1,j,14) + pot(1,j,15)*atar

c        To preserve compatibility with RIPL-2 Koning database
c        b(i,j,15) must be equal to 1., b(1,j,13)  =  0 !!! for Koning OMP
         if((pot(1,j,20).eq.0.) .and.
     >      (pot(1,j,14) + pot(1,j,15)*atar + pot(1,j,16)).eq.0.) then
           b(1,j,13) = pot(1,j,16)
           b(1,j,15) = 1.d0
         endif
       endif

       if(pot(1,j,24).eq.3) then

c        Li & Cai OMP

         b(1,j,15)= 1.d0
         b(1,j,5) = 0.d0   ! No Coulomb correction

         if(b(1,j,1).eq.0.d0) stop 'ERROR in the RIPL pot. compilation'
c                         beta0         beta1
         b(1,j,2)  = -(pot(1,j,3) + pot(1,j,4)*eta)/b(1,j,1)
       endif

c      Wv( Av )
C      b(2,j,6)  =  pot(2,j,1) + pot(2,j,2)*atar
c      added eta dependence for Av parameter (RCN, 12/2006) if pot(2,j,8)<>0
       b(2,j,6)  =  pot(2,j,1) + pot(2,j,2)*atar + pot(2,j,8)*eta

c      Wv( Bv )
c      added eta dependence for Bv parameter (RCN, 12/2006, corrected 08/08) if pot(2,j,9)<>0
       b(2,j,7)  =  pot(2,j,3) + pot(2,j,4)*atar + pot(2,j,9)*eta

C    factor 2*sqrt(eta/atar) explicitly quoted in the definition of the isovector potential
C      WVisov = 2*pot(2,j,8)*sqrt(eta/atar)
       WVisov =   pot(2,j,8)

c      Wd
c      b(4,j,8)  =  pot(4,j,1) + pot(4,j,8)*eta
c      added A dependence for As parameter (RCN, 09/2004) if pot(4,j,7)<>0
C      b(4,j,8)  =  pot(4,j,1) + pot(4,j,8)*eta + pot(4,j,7)*atar
c      added A**(-1/3) dependence for As parameter (RCN, 11/2005) if pot(4,j,9)<>0
c      Wd( As )
       b(4,j,8)  =  pot(4,j,1) + pot(4,j,8)*eta + pot(4,j,7)*atar +
     >                pot(4,j,9)*atar**(-1.d0/3.d0)

C    factor 2*sqrt(eta/atar) explicitly quoted in the definition of the potential
C      WSisov = 2*pot(4,j,8)*sqrt(eta/atar)
       WSisov =   pot(4,j,8)

c      Wd( Cs )
       if(pot(4,j,3).ne.0.) then
c        added A dependence for Cs parameter (RCN, 05/2008) if pot(4,j,10)<>0
         b(4,j,9)  =  pot(4,j,2) + pot(4,j,10)*atar +
     +                pot(4,j,3)/(1.+ exp((atar-pot(4,j,4))/pot(4,j,5)))
       else
         b(4,j,9)  =  pot(4,j,2) + pot(4,j,10)*atar
       endif
c      Wd( Bs )
c      added A dependence for Bs parameter (RCN, 05/2008) if pot(4,j,10)<>0
       b(4,j,10) =  pot(4,j,6) + pot(4,j,11)*atar
C      Wd( q )
       b(4,j,12) =  pot(4,j,12)

c      Vso
       b(5,j,11) =  pot(5,j,10) + pot(5,j,11)*atar
       b(5,j,12) =  pot(5,j,12)

       if(pot(1,j,24).eq.2) then
c
c       Vso(E) calculated from nonlocal approximation
c            as suggested by Perey and Buck
         b(5,j,1)  =  pot(5,j,1) + pot(5,j,2)*atar + pot(5,j,8)*eta
         b(5,j,2)  =  pot(5,j,3) + pot(5,j,4)*atar
         b(5,j,3)  =  pot(5,j,5) + pot(5,j,6)*atar
       endif

c      Wso
       b(6,j,6)  =  pot(6,j,1)
       b(6,j,7)  =  pot(6,j,3)

      return
      end subroutine bcoget

c-----------------------------------------------------------------------------
      subroutine optname
c-----------------------------------------------------------------------------
c
c     Routine to load opname with first author name (or first
c     14 characters).
c
      character*1 opnamex
      dimension opnamex(14)
c
      do 20 i=1,14
  20  opnamex(i)=author(i)
      do 22 i=2,14
      j=i
      if(opnamex(i).eq.' '.and.opnamex(i-1).ne.'.')go to 24
      if(opnamex(i).eq.',')go to 24
  22  continue
  24  do 26 ii=1,j
      i=j-ii+1
      if(opnamex(i).ne.' '.and.opnamex(i).ne.',')go to 28
  26  continue
  28  jj=i
      do 30 i=1,14
  30  opname(i)=' '
      do 32 i=1,jj
  32  opname(i)=opnamex(i)
c
c     Set names for table print
      if(imodel.eq.0)modelz='spherical  '
      if(imodel.eq.1)modelz='rotational '
      if(imodel.eq.2)modelz='vibrational'
      if(imodel.eq.3)modelz='soft rotor '

      if(imodel.eq.4)modelz='rigid+soft '

      if(idr.eq.0) dispz='non-dispersive'
      if(idr.lt.0) dispz='dispersive   (Num.Integr.) '
      if(idr.gt.0) dispz='dispersive   (Anal.Integr.)'
      if(idr.eq.-3)dispz='dispersive+SO(Num.Integr.) '
      if(idr.eq. 3)dispz='dispersive+SO(Anal.Integr.)'
      if(irel.eq.0)relz='non-relativistic'
      if(irel.eq.1)relz='relativ.Schoedr.'
      if(irel.eq.2)relz='relativ + gamma '
      return
      end subroutine optname

c-----------------------------------------------------------------------------
      subroutine couch11
c-----------------------------------------------------------------------------
c
c     First setup for coupled channels parameters
c
c
      ecis1(1:1)='T'
      ecis1(12:12)='T'
c     ecis1(21:21)='T'
      ntar=0
      do 44 nis=1,nisotop
      if(iz(nis).ne.iztar.or.ia(nis).ne.iatar)go to 44
      ntar=nis
!     ntar>0 (Collective Levels exist )
      open(unit=45,file='cclev-omp.dat')
      write(45,'("# Reference Number =",i5," Incident Particle: ",a8,
     +  " imodel=",i1," Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),imodel,iztar,iatar,Efermi
C     write(45,'("Type potential: ",a11,2x,a27,2x,a16)')
C    +              modelz,dispz,relz
C     write(45,'("First Author: ",14a1)') (opname(i),i=1,14)
C     write(45,'("Z-Range=",i2,"-",i3,"  A-Range=",i3,"-",i3,
C    + "  E-Range=",f8.3,"-",f8.3," MeV         Ef =",f7.3)')
C    +  izmin,izmax,iamin,iamax,emin,emax,Efermi

      SELECT CASE (imodel)
        CASE (1) ! Rigid rotor parameters (imodel=1)
          write(35,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,'(A26)') '# LEVELS: rigid rotor     '
          do k=1,ncoll(ntar)
            write(45,99045) EX(k,ntar), SPIn(k,ntar), IPAr(k,ntar)
          enddo

        CASE (2) ! Vibrational rotor parameters (imodel=2)
        write(35,10) ncoll(ntar)
        write(45,10) ncoll(ntar)
          write(45,'(A26)') '# LEVELS: vibrational     '
          do k=1,nvib(ntar)
          write(45,99045) EXV(k,ntar), SPInv(k,ntar),IPArv(k,ntar),
     &                      NPH(k,ntar), DEFv(k,ntar), THEtm(k,ntar)
          enddo

        CASE (3) ! Soft rotor parameters (imodel=3)
        write(35,10) NCOll(ntar)
        write(45,10) NCOll(ntar)
          write(45,'(A26)') '# LEVELS: soft rotor      '
          ! Record 3 from OPTMAN (Hamiltonian parameters)
          write(45,99047)  SR_hw(ntar),SR_amb0(ntar),SR_amg0(ntar),
     +                     SR_gam0(ntar),SR_bet0(ntar),SR_bet4(ntar)
          ! Record 4 from OPTMAN (Hamiltonian parameters)
          write(45,99047)  SR_bb42(ntar),SR_gamg(ntar),SR_delg(ntar),
     +                     SR_bet3(ntar),SR_et0(ntar),SR_amu0(ntar)
          ! Record 5 from OPTMAN (Hamiltonian parameters)
          write(45,99047)  SR_hw0(ntar),SR_bb32(ntar),SR_gamde(ntar),
     +                     SR_dpar(ntar),SR_gshape(ntar)
          do k=1,NCOll(ntar)
           write(45,99049) EXV(k,ntar),SPInv(k,ntar),IPArv(k,ntar),
     +                     SR_ntu(k,ntar),SR_nnb(k,ntar),
     +                     SR_nng(k,ntar),SR_nno(k,ntar)
          enddo

        CASE (4) ! Rigid-soft rotor parameters (imodel=4)
          write(35,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,10) NCOll(ntar),lmaxOM(ntar),
     +      idef(ntar), bandk(ntar),(def(ntar,k),k=2,idef(ntar),2)
          write(45,'(A26)') '# LEVELS: rigid-soft rotor'
          do k=1,NCOll(ntar)
            write(45,99049)  ex(k,ntar),spin(k,ntar),ipar(k,ntar),
     +               SR_ntu(k,ntar),SR_nnb(k,ntar),SR_nng(k,ntar),
     +               SR_nno(k,ntar),defv(k,ntar),defr(k,ntar)
          enddo

        CASE DEFAULT
        write(35,'(A21)') '# SPHERICAL POTENTIAL'

      END SELECT
      close(45)
      if(idr.ne.0) THEN
        write(35,
     +    '("# WARNING: V = VHF + DWV for dispersive potentials")')
      else
        write(35,'("# ")')
      endif
      write(35,'("# En(MeV)   V     RV     AV  ",
     + "    DWV     W      RW     AW  ",
     + "    VD     RVD    AVD     WD     RWD    AWD ",
     + "    VSO    RVSO   AVSO    WSO    RWSO   AWSO   VHFpn    WSpn",
     + "     DVSpn     WVpn    DVVpn   V(Lane)  W(Lane)    VHF" )')
      go to 48
  44  continue
!     ntar<=0 Collective Levels do not exist
      if(ntar.le.0) THEN
        write(35,'("# ncoll=",i2)') 0
        write(35,'(1H#/,"# En(MeV)   V     RV     AV  ",
     + "    DWV     W      RW     AW  ",
     + "    VD     RVD    AVD     WD     RWD    AWD ",
     + "    VSO    RVSO   AVSO    WSO    RWSO   AWSO   VHFpn    WSpn",
     + "     DVSpn     WVpn    DVVpn   V(Lane)  W(Lane)    VHF" )')
      endif
      if(imodel.gt.0 .and. imodel.lt.4 .and. Iflag.ne.-1) then

        write(6,'(" WARNING: RIPL reference number: ",i6)')iref

        write(6,'(" WARNING: Target nucleus izatar=",i6,
     +  " does not have coupled-channel data.")') izatar
          Iflag = -1

        endif
      if(imodel.eq.4 .and. Iflag.ne.1) then
        write(6,'(" ERROR: RIPL reference number: ",i6)')iref

            write(6,'(" ERROR: Target nucleus izatar=",i6,
     +  " does not have soft rotor hamiltonian data.")') izatar
          Iflag = +1
      endif
  48  return
 10   FORMAT(2H# ,'ncoll=',i2,' lmax=',i2,' ndef=',i2,
     +       ' K=',f4.1,' betas:',8(f7.3,1x))
99045 FORMAT (f12.8,f7.1,2I4,1p,2(1x,e11.4))
99047 FORMAT (6(e11.5,1x))
99049 FORMAT (f12.8,f7.1,I4,4I2,2(F8.5,1x))

      end subroutine couch11

c-----------------------------------------------------------------------------
      subroutine couch22
c-----------------------------------------------------------------------------
c
c     Second setup for coupled channels parameters
c
      if(ntar.ne.0)go to 54
      write(ko,'(" PUT ROTATIONAL BAND STATES IN HERE")')
      go to 62
  54  do 56 k=1,ncol
      elx=ex(k,ntar)
      if(k.eq.1) elx=el
      tarspi=spin(k,ntar)
      tarpar='+'
      if(ipar(k,ntar).eq.-1)tarpar='-'
      kom=k
      if(koptom.eq.1)kom=1
      if(k.eq.1) then
        write(ko,'(f5.2,2i2,a1,5f10.5)') tarspi,0,kom,tarpar,elx,
     +  projspi,projmas,tarmas,ztar*zproj
      else
        if (modtyp.eq.6 .or. modtyp.eq.5) kom = 1
        write(ko,'(f5.2,2i2,a1,5f10.5)') tarspi,0,kom,tarpar,elx
      endif
  56  continue
      iqm=idef(ntar)
      iqmax=lmaxOM(ntar)
      aspin=bandk(ntar)
      ik=1
      write(ko,'(2i5,f10.5,i5)') iqm,iqmax,aspin,ik
c     iq=iqm/2
      write(ko,'(7f10.5)') (def(ntar,j),j=2,iqm,2)
  62  return
      end subroutine couch22

c-----------------------------------------------------------------------------
      subroutine masses
c-----------------------------------------------------------------------------
c
c     Routine to retrieve masses and compute separation energies
c     and Fermi energy
c
      real*8 zatar,zaproj
      character*1 dp
      real*8 ck2
      real*8 ampipm,ampi0,ampi
c-------------------------------------------------------------------------------
c  Physical  constants (following ECIS03 definitions and ENDF manual 2001)
c-------------------------------------------------------------------------------
C CONSTANTS COMPUTED FROM THE FUNDAMENTAL CONSTANTS, ATOMIC MASS, HBAR*C
C AND ALPHA, AS GIVEN IN THE EUROPEAN PHYSICAL JOURNAL, PAGE 73, VOLUME
C 15 (2000) REFERRING FOR THESE VALUES TO THE 1998 CODATA SET WHICH MAY
C BE FOUND AT http://physics.nist.gov/constants
C     amu0c2=931.494013 +/- 0.000037 MeV
      amu0c2=931.494013D0
C     hbarc=197.3269601 +/- 0.0000078 (*1.0E-9 eV*cm)
      hbarc=197.3269601D0
C
      ampipm = 1.395688d+02
      ampi0  = 1.349645d+02
      ampi   = (2.d0*ampipm + ampi0)/3.d0
      cso    = (hbarc/ampi)**2
c
      ck2    = (2.d0*amu0c2)/(hbarc**2)
c
c     Compute lab to cm factor
      zatar=float(izatar)
      call energy2(zatar,tarmas,tarspi,tarpar,tarex)
      zaproj=float(izaproj)
      call energy2(zaproj,projmas,projspi,projpar,projex)
c     Compute Fermi energy
      call energy2(zatar-zaproj,dm,ds,dp,tarexm1)
      call energy2(zatar+zaproj,dm,ds,dp,tarexp1)
      efermi=-0.5*(tarexm1-tarexp1+2.*projex)
      etmp = nint(efermi*10000)
      efermi = etmp/10000  ! Fermi energy redued to four numbers after the dot
      return
      end subroutine masses


c-----------------------------------------------------------------------------
      subroutine energy2(za,exactmas,spin,parity,excess)
c-----------------------------------------------------------------------------
c
c     Routine energy2 looks up values of ground-state mass excess
c     (in MeV), spin, and parity. Missing data produces a flag.
c     Note that negative za omits Duflo for missing data.
c
      common/xener/nmass,inpgrd,lza(9200),ener(9200),spnpar(9200)
      character*4 bcd
      character*1 parity
      integer ia
      dimension bcd(120)
      data k13/13/
c
  1   format(28h0***** ground-state data for,i6,19h not in table *****)
  5   format(5x,'+++++ ground state of ',f6.0,' is incompletely describe
     +d,     spin,parity = ',f6.2,a1,2x,'+++++')
  6   format(5x,'+++++',28x,'assignments changed to,   spin,parity =
     +',f6.2,a1,2x,'+++++')
  8   format(///' M A S S   D A T A   I N P U T   T O   S U B R O U T I
     +N E   E N E R G Y',/(20a4))
c
c     FIRST CALL CAUSES DATA TO BE READ IN FROM UNIT 13
      if(inpgrd.eq.12345) go to 10
      open(unit=k13,status='unknown',file='gs-mass-sp.dat')
!MS$IF DEFINED (alone)
      open(unit=44,status='unknown',file='massinfo.out')
!MS$ENDIF
      read(k13,'(20a4)')bcd
!MS$IF DEFINED (alone)
      write(44,8)bcd
!MS$ENDIF
      read(k13,'(i7)')nmass
      read(k13,'(5(i7,f11.6,f7.1))')(lza(n),ener(n),spnpar(n),n=1,nmass)
      rewind k13
      inpgrd = 12345
      close (unit=k13)
c
c     SET FLAG FOR DUFLO OPTION
 10   imiss=1
      if(za) 12,15,20
 12   za=abs(za)
      imiss=0
      go to 20
c
c     Z=0,A=0 IS CONSIDERED A PHOTON
   15 excess=0.
      spin=0.
      parity='-'
      exactmas=0.
      return
c
c     FIND REQUESTED NUCLEUS IN TABLE
C  20 iza=jfix (za)
   20 iza=ifix (real(za))
      ia=mod(iza,1000)
      amass=dble(ia)
      if(iza.lt.lza(1).or.iza.gt.lza(nmass))go to 40
      in=1
      if(iza.eq.lza(1))go to 50
      il=1
      iu=nmass
 122  ii=(il+iu)/2
      if(lza(ii)-iza)124,130,125
 124  il=ii
      go to 126
 125  iu=ii
 126  if(iu-il-1)135,140,122
 130  in=ii
      go to 50
 135  in=il
      go to 150
 140  in=iu
 150  if(lza(in).ne.iza)go to 40
      go to 50
c
c     REQUESTED ISOTOPE IS NOT IN TABLES
 40   continue
!MS$IF DEFINED (alone)
      write(44,1) iza
!MS$ENDIF
      excess= -999999.
      if(imiss.eq.0) stop 7776
      npro=iza/1000
      nneu=ia-npro
      call duflo(nneu,npro,excess)
      exactmas=amass+excess/amu0c2
      go to 200
c
c     PREPARE RETRIEVED DATA FOR EXIT
 50   continue
      excess=ener(in)
      exactmas=amass+excess/amu0c2
      if(spnpar(in).ge.9900.)spin=spnpar(in)-9900.
      if(spnpar(in).ge.9900.)parity='u'
      if(spnpar(in).ge.9900.)go to 200
      if(spnpar(in).ge.100.)parity='+'
      if(spnpar(in).ge.100.)spin=spnpar(in)-100.
      if(spnpar(in).lt.0.)parity='-'
      if(spnpar(in).lt.0.)spin=spnpar(in)+100.
 200  kspin=spin+0.01
      if((parity.ne.'u').and.(kspin.ne.99))return
!MS$IF DEFINED (alone)
      write(44,5) za,spin,parity
!MS$ENDIF
      if(parity.eq.'u')parity='+'
      if(kspin.eq.99)spin=.25*(1.-(-1.)**ia)
!MS$IF DEFINED (alone)
      write(44,6) spin,parity
!MS$ENDIF
      return
      end subroutine energy2

c-----------------------------------------------------------------------------
      subroutine duflo(nn,nz,excess)
c-----------------------------------------------------------------------------
      implicit real*8(A-H,O-Z)
c
c                --------------------------------------------
c                |Nuclear mass formula of Duflo-Zuker (1992)|
c                --------------------------------------------
c
      dimension xmag(6),zmag(5),a(21)
      data zmag/14.,28.,50.,82.,114./
      data xmag/14.,28.,50.,82.,126.,184./
      data a/16.178,18.422,120.146,202.305,12.454,0.73598,5.204,1.0645,
     +1.4206,0.0548,0.1786,.6181,.0988,.0265,-.1537,.3113,-.6650,-.0553,
     +-.0401,.1774,.4523/
c
      x=nn
      z=nz
      t=abs(x-z)*.5
      v=x+z
      s=v**(2./3.)
      u=v**(1./3.)
c     a5=a(5)
c     if(z.gt.x) a5=0.
      E0=a(1)*v-a(2)*s-a(3)*t*t/v+a(4)*t*t/u/v-a(5)*t/u-a(6)*z*z/u
      esh=0.
      esh1=0.
      do 10 k=2,5
        f1=zmag(k-1)
        f2=zmag(k)
        dfz=f2-f1
        if(z.ge.f1.and.z.lt.f2) then
          roz=(z-f1)/dfz
          pz=roz*(roz-1)*dfz
          do 20 l=2,6
            f3=xmag(l-1)
            f4=xmag(l)
            dfn=f4-f3
            if(x.ge.f3.and.x.lt.f4) then
              ron=(x-f3)/dfn
              pn=ron*(ron-1)*dfn
              esh=(pn+pz)*a(8)+a(10)*pn*pz
              xx=2.*ron-1.
              zz=2.*roz-1.
              txxx=pn*xx
              tzzz=pz*zz
              txxz=pn*zz
              tzzx=pz*xx
              kl=l-k
              if(kl.eq.0) esh1=a(k+10)*(txxx+tzzz)+a(k+15)*(txxz+tzzx)
              if(kl.eq.1)
     +          esh1=a(k+11)*txxx-a(k+16)*txxz+a(k+10)*tzzz-a(k+15)*tzzx
              if(kl.eq.2)
     +          esh1=a(k+12)*txxx+a(k+17)*txxz+a(k+10)*tzzz+a(k+15)*tzzx
                edef=a(9)*(pn+pz)+a(11)*pn*pz
                if(esh.lt.edef) esh=edef
            end if
   20     continue
        end if
   10 continue
      ebin=e0+esh+esh1
      nn2=nn/2
      nz2=nz/2
      nn2=2*nn2
      nz2=2*nz2
      if(nn2.ne.nn) ebin=ebin-a(7)/u
      if(nz2.ne.nz) ebin=ebin-a(7)/u
      rneumas=1.008665
      rpromas=1.007825
      amass=nz*rpromas+nn*rneumas-ebin/amu0c2
      excess=(amass-nz-nn)*amu0c2
      return
      end subroutine duflo

c-----------------------------------------------------------------------------
      subroutine egrid(ne,eminx,emaxx,en)
c-----------------------------------------------------------------------------
c     Routine to generate auto energy grid.
c     Currently setup to produce grid for A. Koning
c
      implicit real*8(A-H,O-Z)
      dimension en(500),enx(500)
c
      enx(1)=0.d0
      enx(2)=0.001d0
      enx(3)=0.002d0
c
c     energy range 0-0.010 MeV, delta E = 0.002 MeV
      do 10 i=4,7
  10  enx(i)=enx(i-1)+0.002d0
c
c    energy range 0.01-0.10 MeV, delta E = 0.005 MeV
      do 20 i=8,25
  20  enx(i)=enx(i-1)+0.005d0
c
c    energy range 0.1-1.0 MeV, delta E = 0.05 MeV
      do 30 i=26,43
  30  enx(i)=enx(i-1)+0.05d0
c
c    energy range 1.0-10.0 MeV, delta E = 0.2 MeV
      do 40 i=44,88
  40  enx(i)=enx(i-1)+0.2d0
c
c    energy range 10-200 MeV, delta E = 0.5 MeV
      do 50 i=89,468
      nex=i
  50  enx(i)=enx(i-1)+0.5d0
c
c     Now put energy grid between potential limits.
      ne=1
      do 80 n=2,nex
      if(enx(n).lt.eminx.or.enx(n).gt.emaxx)go to 80
      ne=ne+1
      en(ne)=enx(n)
  80  continue
      return
      end subroutine egrid

c-----------------------------------------------------------------------------
      subroutine tableset
c-----------------------------------------------------------------------------
c
c     Routine to setup for table print
c
      write(35,'(1H#/"# Reference Number =",i5," Incident Particle: ",
     + a8," imodel=",i1)') iref,parname(ipsc),imodel
      write(35,'("# Target: Z=",i3," A=",i4," Number of energies=",i4)')
     +     iztar,iatar,ne
      write(35,'("# Type potential: ",a11,2x,a27,2x,a16)')
     +              modelz,dispz,relz
      write(35,'("# First Author: ",14a1)') (opname(i),i=1,14)
      write(35,'("# Z-Range=",i2,"-",i3,"  A-Range=",i3,"-",i3,
     + "  E-Range=",f8.3,"-",f8.3," MeV         Ef =",f7.3)')
     +  izmin,izmax,iamin,iamax,emin,emax,Efermi
C     if(imodel.eq.0) then
C       write(35,'("# ncoll=",i2)') 0
C       write(35,'(1H#/,"# En(MeV)      V    RV    AV ",
C    +          "   DWV      W    RW    AW ",
C    +          "    VD   RVD   AVD     WD   RWD   AWD ",
C    +          "   VSO   RVSO  AVSO   WSO   RWSO  AWSO   VHFpn WSpn",
C    +          " DVSpn    WVpn    DVVpn   V(Lane)  W(Lane)    E" )')
C     endif
!MS$IF DEFINED (alone)
      if(idr.ne.0 .and. modtyp.le.2) then
        write(34,'("# Reference Number =",i5," Incident Particle: ",a8,
     +  " imodel=",i1," Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),imodel,iztar,iatar,Efermi

        write(36,'("# Reference Number =",i5," Incident Particle: ",a8,
     +  " imodel=",i1," Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),imodel,iztar,iatar,Efermi

        write(34,'("# En(MeV) DWV/w        DWVsym/w       DW>/w      ",
     +             "   DW</w            w         DWS/wd           wd",
     +             "         DWvso/wso       wso")')

        write(36,'("# En(MeV)   DWV    DWV(symm)     DW>      DW<     ",
     +             " DWvso      Vv        Vhf      Vdisp    VHFder    ",
     +             " DWVder    T12der   DVSder  +   DWS   =   Vs")')

      endif
      if(modtyp.le.2) then
        write(37,'("# Reference Number =",i5," Incident Particle: ",a8,
     +  " imodel=",i1," Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +   iref,parname(ipsc),imodel,iztar,iatar,Efermi
      write(37,'("# En(MeV)     RealPot   RVolPot   RSurPot    VSOPot",
     +   "   DVolPot    DVolNonL  DVolSym    ImPot    ImVolPot",
     +   "   ImSurPot   WSOPot")')
      endif
      if(izproj.GT.0 .and. modtyp.le.2) then
        write(38,'("# Reference Number =",i5," Incident Particle: ",a8,
     +  " imodel=",i1," Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),imodel,iztar,iatar,Efermi
        write(38,'("# En(MeV)   RealPot   RVolPot    RSurPot ",
     +  "  Int(Coul)  VInt(Coul) SInt(Coul) VIntDisp(Coul)")')
      endif
!MS$ENDIF
      return
      end subroutine tableset

c-----------------------------------------------------------------------------
      REAL*8 FUNCTION Vhf(einp,alpha_PB,beta_PB,gamma_PB)
c-----------------------------------------------------------------------------
c
c     According to Morillon B, Romain P, PRC70(2004)014601
c
c     Originally coded in c++ by Morillon B. and Romain P.
c
c     Coded in FORTRAN and tested by RCN, August 2004.
c
      implicit real*8(A-H,O-Z)
      real*8 einp,alpha_PB,beta_PB,gamma_PB
      real*8 Vtmp,Etmp,miu_sur_hbar2, coef1, coef2
      integer niter

c     getting amu
      xtmp=xkine(einp)
      miu_sur_hbar2 = amu / hbarc**2
      coef1 = -0.5d0 * beta_PB**2 * miu_sur_hbar2
      coef2 =  4.0d0 * (gamma_PB * miu_sur_hbar2)**2

      niter = 0
C     Vhf = -45.d0

      Vhf = alpha_PB
10    niter = niter + 1
      Vtmp = Vhf
      Etmp = einp - Vtmp
      Vhf = alpha_PB * dexp(coef1 * Etmp + coef2 * Etmp**2)
      if( abs(Vhf - Vtmp) .GT. 0.0001 .AND.  niter.LT.10000) goto 10
      return
      end FUNCTION Vhf

c-----------------------------------------------------------------------------
      REAL*8 FUNCTION xkine(ei)
c-----------------------------------------------------------------------------
c***********************************************************************
c     From lab to CM (the input quantity is el = Elab)
c***********************************************************************
c     RCN 08/2004, xkine calculated by relativistic kinematics when needed

      mtot = (tarmas+projmas)
      xratio = tarmas / mtot

      if(irel .eq. 0) then
c
c-----------------------------------------------------------------------
c  Classical    kinematics (energy independent amu and xkine)
c-----------------------------------------------------------------------
c
          amu = projmas*tarmas / mtot * amu0c2
          xkine = tarmas / mtot
c         e1  = el*xkine
c         w2  = ck2*amu
c         ak2 = w2*e1
      else
c
c-----------------------------------------------------------------------
c  Relativistic kinematics
c-----------------------------------------------------------------------
c
c         e1  = amu0c2*mtot*
c    * (        DSQRT(1.d0 +
c    *       2.d0*el/(amu0c2*tarmas*((1.d0+projmas/tarmas)**2))) - 1.d0)
          p2  = (ei*(ei + 2.d0*amu0c2*projmas)) /
     *          ((1.d0+projmas/tarmas)**2 + 2.d0*ei/(amu0c2*tarmas))
c         ak2 = p2 / (hbarc*hbarc)
          etoti = DSQRT((amu0c2*projmas)**2 + p2)
          etott = DSQRT((amu0c2*tarmas)**2  + p2)
          amu   = etoti*etott / (etoti + etott)
c         amu   = amu / amu0c2
          xkine = etott / (etoti + etott)
      endif
      return
      end FUNCTION xkine

C
C==========================================================================
C     AUTHOR: Dr. Roberto Capote Noy
c
C     e-mail: rcapotenoy@yahoo.com; r.capotenoy@iaea.org
C
C     DISPERSIVE OPTICAL MODEL POTENTIAL PACKAGE
C
c     Analytical dispersive integrals are included
c     see Quesada JM et al, Computer Physics Communications 153(2003) 97
c                           Phys. Rev. C67(2003) 067601
C
      REAL*8 FUNCTION DOM_INT_Wv(Ef,Ep,Av,Bv,Einc,n,DerivIntWv)
C
C      Analytical dispersive integral and its derivative for
C      Wv(E)=Av*(E-Ep)**n/( (E-Ep)**n + Bv**n )  for E>Ep
C      Wv(E)=Wv(2*Ef-E)                          for E<2Ef-Ep
C      Wv(E)=0                                     OTHERWISE
C
      IMPLICIT NONE
      REAL*8 Ef,Ep,Av,Bv,E,pi,Einc
      REAL*8 E0,Ex,Eplus,Emin,Rs,ResEmin,ResEplus
      REAL*8 DerEmin, DerEplus, Rds, DerivIntWv
      COMPLEX*16 Pj,I,Zj,Ztmp
      COMPLEX*16 Fs,Ds
      INTEGER N,j,IS
      DATA I/(0.d0,1.d0)/
      pi=4.d0*atan(1.d0)

      IS = 1
      E = Einc
      IF(Einc.LE.Ef) THEN
        E=2.d0*Ef-Einc
C       Odd function
        IS = -1
      ENDIF
      E0 = Ep - Ef
      Ex = E  - Ef
      Eplus = Ex + E0
      Emin  = Ex - E0
      DOM_INT_Wv = 0.d0
      DerivIntWv = 0.d0

      ResEmin  =  Emin**n / (Emin**n + Bv**n)

      DerEmin  =  Emin**(n-1) *
     >           ( Emin**n + Bv**n*(1.d0 + n*log(dabs(Emin)) ) )
     >           / (Emin**n + Bv**n)**2

      ResEplus = -Eplus**n / (Eplus**n + Bv**n)

      DerEplus = -Eplus**(n-1) *
     >           ( Eplus**n + Bv**n*(1.d0+n*log(Eplus)) )
     >           / (Eplus**n + Bv**n)**2

C----------------------------------
C     Complex arithmetic follows
      Fs = (0.d0,0.d0)
      Ds = (0.d0,0.d0)
      do j=1,n
       Ztmp = I*(2*j-1)/dble(n)*pi
       Pj = Bv*exp(Ztmp)
       Zj = Pj * (2*Pj +Eplus -Emin) * Ex
       Zj = Zj / ( (Pj+E0) * (Pj+Eplus) * (Pj-Emin) )
       Fs = Fs + Zj*log(-Pj)
       Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*log(-Pj)
     >           /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
      enddo

      IF(ABS(IMAG(Fs)).gt.1.e-4) STOP
     > 'ERROR: (F) Too big imag part in DWv'
      Rs  = REAL(Fs)
      IF(ABS(IMAG(Ds)).gt.1.e-4) STOP
     > 'ERROR: (D) Too big imag part in DWv'
      Rds = REAL(Ds)
C----------------------------------

      DOM_INT_Wv = -Av/pi*IS*
     &  ( Rs/n  + (ResEplus*log(Eplus) + ResEmin*log(dabs(Emin))) )

      DerivIntWv = -Av/pi*IS*( Rds/n + (DerEplus + DerEmin) )

      RETURN
      END FUNCTION DOM_INT_Wv

c-----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION
     +          DOM_INT_Ws(Ef,Ep,As,Bs,Cs,Einc,m,DerivIntWs)
c-----------------------------------------------------------------------------
C
C      Analytical dispersive integral and its derivative for
C      Ws(E)=As*(E-Ep)**m/( (E-Ep)**m + Bs**m ) * exp(-Cs*(E-Ep)) for E>Ep
C      Ws(E)=Ws(2*Ef-E)                                           for E<2Ef-Ep
C      Ws(E)=0                                                    OTHERWISE
C
      IMPLICIT NONE
      REAL*8 Ef,Ep,As,Bs,Cs,E,Einc
      COMPLEX*16 I,Pj,Zj,Ztmp
      COMPLEX*16 Fs,Ds
      REAL*8 E0,Ex,Eplus,Emin,pi
      REAL*8 Rs,ResEmin,ResEplus
      REAL*8 DerivIntWs,DerEmin,DerEplus,Rds
      INTEGER m,j,IS

C     To be uncommented for release, problem with MS FORTRAN compilers (link error)
C     COMPLEX*16 zfi
C     REAL*8 Ein

      DATA I/(0.d0,1.d0)/
      pi=4.d0*atan(1.d0)

      IS = 1
      E = Einc
      IF(Einc.LE.Ef) THEN
        E=2.d0*Ef-Einc
C       Odd function
        IS = -1
      ENDIF
      E0 = Ep - Ef
      Ex = E  - Ef
      Eplus = Ex + E0
      Emin  = Ex - E0
      DOM_INT_Ws = 0.d0
      DerivIntWs = 0.d0

      ResEmin  =  Emin**m / (Emin**m + Bs**m)

      DerEmin  = -Emin**(m-1) *
     >          ( Emin**m + Bs**m + ( -Cs*Emin**(m+1) +
     >            Bs**m *(-Cs*Emin+m) ) * exp(-Cs*Emin)*EIn(Cs*Emin) )
     >           / (Emin**m + Bs**m)**2

      ResEplus = -Eplus**m / (Eplus**m + Bs**m)

      DerEplus =  Eplus**(m-1) *
     >          ( Eplus**m + Bs**m + ( Cs*Eplus**(m+1) +
     >            Bs**m *(Cs*Eplus+m) ) * exp(Cs*Eplus)*EIn(-Cs*Eplus) )
     >           / (Eplus**m + Bs**m)**2

C----------------------------------
C     Complex arithmetic follows
      Fs = (0.d0,0.d0)
      Ds = (0.d0,0.d0)
      do j=1,m
       Ztmp = I*(2*j-1)/dble(m)*pi
       Pj = Bs*exp(Ztmp)
       Zj = Pj * (2*Pj +Eplus -Emin) * Ex
       Zj = Zj / (Pj+E0) / (Pj+Eplus) / (Pj-Emin)
       Fs = Fs + Zj* zfi(-Pj*Cs)
       Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*zfi(-Pj*Cs)
     >           /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
      enddo

      IF(ABS(IMAG(Fs)).gt.1.e-4) STOP
     >  'ERROR: (F) Too big imag part in DWs'
      Rs  = REAL(Fs)
      IF(ABS(IMAG(Ds)).gt.1.e-4) STOP
     >  'ERROR: (D) Too big imag part in DWs'
      Rds = REAL(Ds)
C----------------------------------

      DOM_INT_Ws = As/pi*IS*(Rs/m
     &                  - ResEplus*exp(Cs*Eplus)*EIn(-Cs*Eplus)
     &                  - ResEmin*exp(-Cs*Emin)*EIn(Cs*Emin) )

      DerivIntWs = As/pi*IS*( Rds/m + DerEplus + DerEmin)

      RETURN
      END FUNCTION DOM_INT_Ws

C
C-----FUNCTION TO EVALUATE exp(Z)*E1(Z)
C
c-----------------------------------------------------------------------------
      COMPLEX*16 function zfi(za)
c-----------------------------------------------------------------------------
C
C Complex exponential integral function multiplied by exponential
C
C AUTHOR: J. Raynal
C
      IMPLICIT NONE
      real*8 aj
      complex*16 za,y
      integer m,i
      zfi=0.d0
      if (za.eq.0.) return
      if (dabs(dreal(za)+18.5d0).ge.25.d0) go to 3
      if (dsqrt(625.d0-(dreal(za)+18.5d0)**2)/1.665d0.lt.dabs(dimag(za))
     1) go to 3
      zfi=-.57721566490153d0-cdlog(za)
      y=1.d0
      do 1 m=1,2000
      aj=m
      y=-y*za/aj
      if (cdabs(y).lt.1.d-15*cdabs(zfi)) go to 2
    1 zfi=zfi-y/aj
    2 zfi=cdexp(za)*zfi
      return
    3 do 4 i=1,20
      aj=21-i
      zfi=aj/(za+zfi)
    4 zfi=aj/(1.d0+zfi)
      zfi=1.d0/(zfi+za)
      return
      end function zfi

C
C-----FUNCTION TO EVALUATE Ei(X)
C
c-----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION EIn(X)
c-----------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 FAC, H, X
      INTEGER N
      EIn = 0.57721566490153d0+LOG(ABS(X))
      FAC = 1.0
      DO N = 1,100
      H = FLOAT(N)
      FAC = FAC*H
      EIn = EIn + X**N/(H*FAC)
      ENDDO
      RETURN
      END FUNCTION EIn

c-----------------------------------------------------------------------------
      real*8 function WV(A,B,Ep,Ef,E,n)
c-----------------------------------------------------------------------------
      IMPLICIT NONE
      real*8 A,B,Ep,Ef,E,ee
      integer n

      WV=0.d0
      if(E.LE.Ef) E=2.d0*Ef-E
      if(E.LT.Ep) return

      ee=(E-Ep)**n
      WV=A*ee/(ee+B**n)

      return
      end function WV

      real*8 function WDD(A,B,C,Ep,Ef,E,m)
      IMPLICIT NONE
      real*8 A,B,C,Ep,Ef,E,ee,arg
      integer m

      WDD=0.d0
      if(E.LE.Ef) E=2.d0*Ef-E
      if(E.LT.Ep) return

      arg=C*(E-Ep)
      IF(arg.GT.15) return
      ee=(E-Ep)**m
      WDD=A*ee/(ee+B**m)*EXP(-arg)
      return
      end function WDD



c-----------------------------------------------------------------------------
      REAL*8 FUNCTION DOM_int_T1_VDK(Av,Ef,Ea,E)
c-----------------------------------------------------------------------------
C
C     Integral over E' corresponding to nonlocal additions T1(E'<<0)
C
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 Eax, Ex, Ea
      Pi=4.d0*ATAN(1.d0)
      Ex=E-Ef
      Ea2=Ea**2
      Ea3=Ea2*Ea
      Ea4=Ea2*Ea2
      Eax=Ex+Ea

C     Following formula (16) of the Vanderkam paper
      DDEN1 = 2*Pi*Ea2*(Eax**2+Ea2)
      DNUM1 = 0.5d0*Pi*Ex*(Ea3-Ea2*Eax) - Ea2*Ex*(Ex+2*Ea)*dlog(Ea)
     >   + Ea2*Eax**2*dlog(dabs(Eax/Ea)) + Ea2*Eax**2*dlog(dabs(Eax))
     >   - Ea4*dlog(Ea)

      DOM_int_T1_VDK = -Av*DNUM1/DDEN1

      RETURN
      END FUNCTION DOM_int_T1_VDK

c-----------------------------------------------------------------------------
      REAL*8 FUNCTION DOM_int_T1(Ef,Ea,E)
c-----------------------------------------------------------------------------
C
C     Integral over E' corresponding to nonlocal additions T1(E'<<0)
C
C     Following formula (20) CPC
C
      IMPLICIT NONE

      real*8 E,Ea,Ef,Ex,Ea2,Eax,Pi,T11,T12,T13
      Pi=4.d0*ATAN(1.d0)

      Ex=E-Ef
      Ea2=Ea**2
      Eax=dabs(Ex+Ea)

      T11 = 0.5d0*log(Ea)/Ex
      T12 =  ( (2*Ea+Ex)*log(Ea)+0.5d0*pi*Ex )
     >      /(2.*(Eax**2 + Ea2))
      T13 = -Eax**2*log(Eax)/(Ex*(Eax**2+Ea2))

      DOM_int_T1 = Ex/Pi*(T11+T12+T13)


      RETURN
      END FUNCTION DOM_int_T1
C
c-----------------------------------------------------------------------------
      REAL*8 FUNCTION DOM_int_T2_JM(Ef,Ea,E)
c-----------------------------------------------------------------------------
C
C     Integral over E' corresponding to nonlocal additions T2(E'>>0)
C
      IMPLICIT NONE
      real*8 E,Ea,Ef,EL,Pi,DOM_int

      Pi=4.d0*ATAN(1.d0)
      EL=Ef+Ea
      DOM_int= 1.d0 / Pi * (
     >      sqrt(abs(Ef)) * atan( (2*sqrt(EL*abs(Ef)))/(EL-abs(Ef)) )
     > +    EL**1.5d0/(2*Ef)*log(Ea/EL) )

      IF(E.GT.EL) THEN

      DOM_int = DOM_int + 1.d0/Pi* (
     >  sqrt(E) * log( (sqrt(E)+sqrt(EL)) / (sqrt(E)-sqrt(EL)) ) +
     >  1.5d0*sqrt(EL)*log((E-EL)/Ea) + EL**1.5d0/(2*E)*log(EL/(E-EL)) )

      ELSEIF(E.EQ.EL) THEN

      DOM_int = DOM_int + 1.d0/Pi*1.5d0*sqrt(EL)
     > *log((2**(4.d0/3.d0)*EL)/Ea)

      ELSEIF(E.GT.0.d0 .AND. E.LE.EL) THEN

      DOM_int = DOM_int + 1.d0/Pi * (
     > sqrt(e) * log( (sqrt(E)+sqrt(EL)) / (sqrt(EL)-sqrt(E)) ) +
     > 1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )

      ELSEIF(E.EQ.0.d0) THEN

C     CPC formula
C     DOM_int = DOM_int + 1.d0/Pi*( 0.5*EL**(1./3.)
C    > + log(EL/Ea) + 0.5d0*sqrt(EL) )

      DOM_int = DOM_int + 1.d0/Pi*1.5d0*sqrt(EL)
     > *log(EL/Ea)+0.5d0*sqrt(EL)

      ELSE

      DOM_int = DOM_int + 1.d0/Pi * (
     > -sqrt(abs(E))*atan( 2*(sqrt(EL-abs(E))) / (EL-abs(E)) ) +
     > 1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )

      ENDIF

      DOM_int_T2_JM = DOM_int

      RETURN
      END FUNCTION DOM_int_T2_JM

      DOUBLE PRECISION FUNCTION DOM_int_T2(Ef,Ea,E)
C
C     Integral over E' corresponding to nonlocal additions T2(E'>>0)
C
C     Following Eq.(24) PRC94(2016)064605
C     (CPC equation (20) contains typos !!)
C
      IMPLICIT NONE
      DOUBLE PRECISION E,Ea,Ef,EL,Pi

      Pi=4.d0*ATAN(1.d0)
      EL=Ef+Ea
      DOM_int_T2= 1.d0 / Pi * (
     >      sqrt(abs(Ef)) * atan( (2*sqrt(EL*abs(Ef)))/(EL-abs(Ef)) )
     > +    EL**1.5d0/(2*Ef)*log(Ea/EL) )

      IF(E.GT.EL) THEN

      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi* (
     >  sqrt(E) * log( (sqrt(E)+sqrt(EL)) / (sqrt(E)-sqrt(EL)) ) +
     >  1.5d0*sqrt(EL)*log((E-EL)/Ea) + EL**1.5d0/(2*E)*log(EL/(E-EL)) )

      ELSEIF(E.EQ.EL) THEN

      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi*1.5d0*sqrt(EL)
     > *log((2**(4.d0/3.d0)*EL)/Ea)

C     ELSEIF(E.GT.0.d0 .AND. E.LE.EL) THEN
      ELSEIF(E.GT.0.d0 .AND. E.LT.EL) THEN

      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi * (
     > sqrt(E) * log( (sqrt(E)+sqrt(EL)) / (sqrt(EL)-sqrt(E)) ) +
     > 1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )

      ELSEIF(E.EQ.0.d0) THEN

C     CPC (wrong) formula
C     DOM_int = DOM_int + 1.d0/Pi*( 0.5*EL**(1./3.)
C    > + log(EL/Ea) + 0.5d0*sqrt(EL) )

C     PRC94(2016)064605
      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi*1.5d0*sqrt(EL)
     > *log(EL/Ea)+0.5d0*sqrt(EL)

      ELSE ! E<0

C     PRC94(2016)064605
      DOM_int_T2 = DOM_int_T2 + 1.d0/Pi * (
     > -sqrt(abs(E))*atan( 2*sqrt(EL*abs(E)) / (EL-abs(E)) ) +
     > 1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )

      ENDIF
      RETURN
      END FUNCTION DOM_int_T2

c-----------------------------------------------------------------------------
      REAL*8 FUNCTION DOM_int_T2_CPC(Ef,Ea,E)
c-----------------------------------------------------------------------------
C
C     Integral over E' corresponding to nonlocal additions T2(E'>>0)
C
      IMPLICIT REAL*8(A-H,O-Z)

      Pi=4.d0*ATAN(1.d0)
      Eltt=Ef+Ea

      R1=1.5*DSQRT(Eltt)*dLOG(abs((Eltt-E)/Ea))

      IF(E.eq.0.d0) THEN
        R2=0.5*Eltt**1.5d0*(1.d0/El-dlog(abs(Eltt/Ea))/Ef)
      ELSE
        R2=0.5*Eltt**1.5d0/(E*Ef)*
     >      ( Ef*dlog(dabs(Eltt/(Eltt-E))) -E*dlog(dabs(Eltt/Ea)) )
      ENDIF

      R3=2*DSQRT(dABS(Ef))*( 0.5d0*Pi - atan( DSQRT(Eltt/dabs(Ef)) ) )

      IF(E.GE.0.d0) THEN
        R4=DSQRT(E)*
     >      dlog(dabs( (dsqrt(Eltt)+dsqrt(E))/(dsqrt(Eltt)-dsqrt(E)) ))
      ELSE
        R4=-2.d0*DSQRT(dabs(E))*
     >     ( 0.5d0*Pi - atan( DSQRT(dabs(Eltt/E)) ) )
      ENDIF

      DOM_int_T2_CPC =  1.d0/Pi*(R1+R2+R3+R4)

      RETURN
      END FUNCTION DOM_int_T2_CPC

c*******************************************
c-----------------------------------------------------------------------------
      REAL*8 FUNCTION DOM_int(DELTAF,F,Ef,Eint,Ecut,E,WE_cte)
c-----------------------------------------------------------------------------
C
C     DOM integral (20 points Gauss-Legendre)
C
C     Divided in two intervals for higher accuracy
C     The first interval corresponds to peak of the integrand
C
      DOUBLE PRECISION Eint,Ef,Ecut,WE_cte,E
      DOUBLE PRECISION F,WG,XG,WWW,XXX,DELTAF
      DOUBLE PRECISION ABSC1,CENTR1,HLGTH1,RESG1
      DOUBLE PRECISION ABSC2,CENTR2,HLGTH2,RESG2
      INTEGER J
      EXTERNAL F
      DIMENSION XG(10),WG(10)
C
C     THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C     BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C     CORRESPONDING WEIGHTS ARE GIVEN.
C
C     XG - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
C     WG - WEIGHTS OF THE 20-POINT GAUSS RULE
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      DATA WG  (  1) / 0.0176140071 3915211831 1861962351 853 D0 /
      DATA WG  (  2) / 0.0406014298 0038694133 1039952274 932 D0 /
      DATA WG  (  3) / 0.0626720483 3410906356 9506535187 042 D0 /
      DATA WG  (  4) / 0.0832767415 7670474872 4758143222 046 D0 /
      DATA WG  (  5) / 0.1019301198 1724043503 6750135480 350 D0 /
      DATA WG  (  6) / 0.1181945319 6151841731 2377377711 382 D0 /
      DATA WG  (  7) / 0.1316886384 4917662689 8494499748 163 D0 /
      DATA WG  (  8) / 0.1420961093 1838205132 9298325067 165 D0 /
      DATA WG  (  9) / 0.1491729864 7260374678 7828737001 969 D0 /
      DATA WG  ( 10) / 0.1527533871 3072585069 8084331955 098 D0 /
C
      DATA XG( 1) / 0.9931285991 8509492478 6122388471 320 D0 /
      DATA XG( 2) / 0.9639719272 7791379126 7666131197 277 D0 /
      DATA XG( 3) / 0.9122344282 5132590586 7752441203 298 D0 /
      DATA XG( 4) / 0.8391169718 2221882339 4529061701 521 D0 /
      DATA XG( 5) / 0.7463319064 6015079261 4305070355 642 D0 /
      DATA XG( 6) / 0.6360536807 2651502545 2836696226 286 D0 /
      DATA XG( 7) / 0.5108670019 5082709800 4364050955 251 D0 /
      DATA XG( 8) / 0.3737060887 1541956067 2548177024 927 D0 /
      DATA XG( 9) / 0.2277858511 4164507808 0496195368 575 D0 /
      DATA XG(10) / 0.0765265211 3349733375 4640409398 838 D0 /
C
      CENTR1 = 0.5D+00*(Ef+Eint)
      HLGTH1 = 0.5D+00*(Eint-Ef)
      CENTR2 = 0.5D+00*(Ecut+Eint)
      HLGTH2 = 0.5D+00*(Ecut-Eint)
C
C     COMPUTE THE 20-POINT GAUSS-KRONROD APPROXIMATION
C     TO THE INTEGRAL in TWO INTERVALS (Ef - Eint, Eint - Ecut)
C
      RESG1 = 0.0D+00
      RESG2 = 0.0D+00
      DO J=1,10
      XXX=XG(J)
      WWW=WG(J)
        ABSC1 = HLGTH1*XXX
      RESG1=RESG1 + WWW*(DELTAF(F,CENTR1-ABSC1)+DELTAF(F,CENTR1+ABSC1))
        ABSC2 = HLGTH2*XXX
      RESG2=RESG2 + WWW*(DELTAF(F,CENTR2-ABSC2)+DELTAF(F,CENTR2+ABSC2))
      ENDDO
c
      CORR=0.5d0*WE_cte/(E-Ef)*dlog((Ecut-(E-Ef))/(Ecut+(E-Ef)))
      DOM_int = ( RESG1*HLGTH1+RESG2*HLGTH2 + CORR) *
     >           (E-Ef)/(acos(-1.d0))
c
      RETURN
      END FUNCTION DOM_int

c-----------------------------------------------------------------------------
      subroutine VolIn(itype,depth,rad,dif,volint)
c-----------------------------------------------------------------------------
c**********************************************************
c  Calculate volume integrals (from SCAT2000, O.Bersillon)*
c---------------------------------------------------------*
c  ITYPE  = 1 Woods-Saxon potential                       *
c           2 Woods-Saxon derivative potential            *
c  DEPTH  = depth       of the potential                  *
c  RAD    = radius      of the potential                  *
c  DIF    = diffuseness of the potential                  *
c  VOLINT = volume integral                               *
c**********************************************************
c
      implicit none
c
      double precision depth,rad,dif,volint
      integer itype
      double precision arg
      double precision zero,one,three,four, pi
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             three /3.0d+00/
      data             four  /4.0d+00/
      data             pi    /3.1415926d0/
c=======================================================================
c
      volint = 0.d0
      if(rad.eq.zero .or. depth.eq.zero) return
c
      arg = (pi*dif/rad)**2
c
      if(itype .eq. 1) then
        volint = four*pi*(rad**3)*depth*(one + arg) / three
      else if(itype .eq.  2) then
        volint = four*four*pi*dif*(rad**2)*depth*(one + arg/three)
      endif
      return
      end subroutine VolIn

      REAL*8 FUNCTION VolumeInt2d
     >   (LAMBDA,BETAFF,MAX_LAMBDA,vdd,rdd,add,KEY) !INPUT

      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION BETAFF(MAX_LAMBDA)
      DIMENSION X(100),W(100)

      common/pgaus/x,w,n1
      common/parameters/vdepth,radius,diffuss  ! common for potential parameters
      common/armon/LANDA,LANDAM
c
c Radius is the theta dependent deformed radius (carried through a common statement),
c which takes different values for the varius shapes (HF, volume,surface & spin-orbit)
c XR is the increment applied to RD for the upper limit of radial integration
      DATA PI/3.1415926d0/,XR/30.d0/

C     EXTERNAL VolumeR,SurfaceR,SpinOrbitR

      n1 = 60             ! 60 points for integration both in r and theta
      call gauleg(x,w,n1) ! Calculating Gauss-Legendre abscissae and weight

      LANDA  = LAMBDA

      XX=0.5d0*pi
      VolumeInt2d = 0.d0
      IF(VDD.eq.0.d0) RETURN

      SSr = 0.d0
C
C     Integration over theta
C
      DO J=1,(n1+1)/2
        DX=XX*X(J)
        THETAP = XX + DX
        THETAM = XX - DX
        XP = cos(THETAP)
        XM = cos(THETAM)
        SUMP=0.d0
        SUMM=0.d0
        DO L=1,MAX_LAMBDA
          if(betaFF(L).eq.0.d0) cycle
          SUMP=SUMP+SQRT((2.*L+1.)/4./PI)*PLGNDR(L,0,XP)*betaFF(L)
          SUMM=SUMM+SQRT((2.*L+1.)/4./PI)*PLGNDR(L,0,XM)*betaFF(L)
        ENDDO
C
        YL0P=SQRT((2*LANDA+1)/(4*PI))*PLGNDR(LANDA,0,XP)
        YL0M=SQRT((2*LANDA+1)/(4*PI))*PLGNDR(LANDA,0,XM)
        FTETAP = 2.d0*PI*YL0P*SIN(THETAP)
        FTETAM = 2.d0*PI*YL0M*SIN(THETAM)
C
C       Integration over radius(theta)
C
C       Real parts
C
        vdepth  = vdd
        diffuss = add
C
        radius  = rdd*(1.d0 + SUMP)
        if(key.eq.1) XVp  = qgauss(VolumeR   ,0.d0,radius+XR) * FTETAP
        if(key.eq.2) XVp  = qgauss(SurfaceR  ,0.d0,radius+XR) * FTETAP
        if(key.eq.3) XVp  = qgauss(SpinOrbitR,0.d0,radius+XR) * FTETAP
C
        radius  = rdd*(1.d0 + SUMM)
        if(key.eq.1) XVm  = qgauss(VolumeR   ,0.d0,radius+XR) * FTETAM
        if(key.eq.2) XVm  = qgauss(SurfaceR  ,0.d0,radius+XR) * FTETAM
        if(key.eq.3) XVm  = qgauss(SpinOrbitR,0.d0,radius+XR) * FTETAM

        SSR = SSR + W(J)*(XVp + XVm)
      ENDDO
      VolumeInt2d = Sqrt(4*pi)*XX*SSR
      return
      end FUNCTION VolumeInt2d
C
C-----------------------------------------------------------------------
C
c-----------------------------------------------------------------------------
      REAL*8 FUNCTION qgauss(VVHF,A,B)
c-----------------------------------------------------------------------------
C
C     Gauss-Legendre quadrature integration method
C
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION X(100),W(100)

      common/pgaus/x,w,n1
      EXTERNAL VVHF

      XM=0.5*(B+A)
      XR=0.5*(B-A)
      SS=0.d0
      DO J=1,(n1+1)/2
         DX=XR*X(J)
         SS=SS+W(J)*(VVHF(XM+DX)+VVHF(XM-DX))
      ENDDO
      qgauss = XR*SS
      RETURN
      END   FUNCTION qgauss
C
C-----------------------------------------------------------------------
C
c-----------------------------------------------------------------------------
      SUBROUTINE GAULEG(X,W,N)
c-----------------------------------------------------------------------------
C
C     ABCISSAS AND WEIGHTS ARE CALCULATED
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),W(N)
      PARAMETER (EPS=3. D-14)
      M=(N+1)/2
      DO I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
 1      P1=1.D0
        P2=0.D0
        DO J=1,N
          P3=P2
          P2=P1
          P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
        ENDDO
        PP=N*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS) GO TO 1
        X((N+1)/2+1-I) = Z
        W((N+1)/2+1-I) = 2.D0/((1.D0-Z*Z)*PP*PP)
      ENDDO
      RETURN
      END   SUBROUTINE GAULEG
C
C --------------------------------------------------------------------------
C
c-----------------------------------------------------------------------------
      REAL*8 FUNCTION PLGNDR(L,M,X)
c-----------------------------------------------------------------------------
      IMPLICIT real*8 (A-H,O-Z)
C
C Computes the associated Legendre Polinomials  P(lm). Here m and l are
C integers satisfying 0<=m<=l, while x lies in the range -1<=x<=1.
C
      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1) then
        write(6,*) 'ERROR: bad arguments in PLGNDR'
        stop
      endif
      PMM=1.d0
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.d0-X)*(1.d0+x))
        RFACT=1.d0
        DO I=1,M
          PMM=PMM*RFACT*SOMX2
          RFACT=RFACT+2.d0
        ENDDO
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
          ENDDO
          PLGNDR=PLL
        ENDIF
      ENDIF
      RETURN
      END   FUNCTION PLGNDR

C ---------------------------------------------------------
C
c-----------------------------------------------------------------------------
      REAL*8 FUNCTION FACT(N)
c-----------------------------------------------------------------------------
C
C     Returns  the value N! as a floating point number
C
      IMPLICIT real*8 (A-H,O-Z)

      DIMENSION A(6)
      DATA NTOP,A/0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/
      IF (N.LT.0) THEN
        write(6,*) 'ERROR: NEGATIVE FACTORIAL CALL'
        stop
      ELSE IF (N.LE.NTOP) THEN
         FACT=A(N+1)
      ELSE IF (N.LE.6) THEN
         DO J=NTOP+1,N
           A(J+1)=J*A(J)
         ENDDO
         NTOP=N
         FACT=A(N+1)
      ELSE
        write(6,*) 'ERROR: BAD LAMBDA AND/OR MU'
        stop
      ENDIF
      RETURN
      END   FUNCTION FACT

c-----------------------------------------------------------------------------
      REAL*8 FUNCTION DELTA_WV(WVf,y)
c-----------------------------------------------------------------------------

      REAL*8  E,Ef,A,B,Ep,y,WVf,WDE,WVE,Ea
      INTEGER n
      COMMON /energy/E,Ef,Ep,Ea
      COMMON /Wenerg/WDE,WVE
      COMMON /pdatav/A,B,n

      DELTA_WV = (WVf(A,B,Ep,Ef,y,n) - WVE)
     >           /((y-Ef)**2-(E-Ef)**2)
      RETURN
      END FUNCTION DELTA_WV


c-----------------------------------------------------------------------------
      REAL*8 FUNCTION WVf(A,B,Ep,Ef,E,n)
c-----------------------------------------------------------------------------

      REAL*8  A,B,Ep,E,Ef
      INTEGER n

      WVf = 0.d0
      IF (E.LE.Ef) E = 2.d0*Ef-E
      IF (E.LE.Ep) RETURN

      ee = (E-Ep)**n
      WVf = A*ee/(ee+B**n)

      RETURN
      END FUNCTION WVf

c-----------------------------------------------------------------------------
      REAL*8 FUNCTION DELTA_WD(WDf,y)
c-----------------------------------------------------------------------------

      REAL*8 E,Ef,A,B,C,Ep,y,WDf,WDE,WVE,Ea
      INTEGER m,iq
      COMMON /energy/E,Ef,Ep,Ea
      COMMON /Wenerg/WDE,WVE
      COMMON /pdatas/A,B,C,m,iq

      DELTA_WD = (WDf(A,B,C,Ep,y,m,iq) - WDE)
     >            /((y-Ef)**2-(E-Ef)**2)

      RETURN
      END FUNCTION DELTA_WD

c-----------------------------------------------------------------------------
      REAL*8 FUNCTION WDf(A,B,C,Ep,E,m,iq)
c-----------------------------------------------------------------------------
      REAL*8 A,B,C,Ep,E,ee,arg
      INTEGER m,iq

      WDf = 0.d0
      IF (E.Lt.Ep) RETURN
      arg = C*(E-Ep)**iq
      IF (arg.GT.15) RETURN
      ee = (E-Ep)**m
      WDf = A*ee/(ee+B**m)*EXP(-arg)

      END FUNCTION WDf

c==   START integrands in r for the various terms
c-----------------------------------------------------------------------------
      real*8 function VolumeR(r)
c-----------------------------------------------------------------------------
      IMPLICIT real*8 (A-H,O-Z)
      common/parameters/vdepth,radius,diffuss
      common/armon/LANDA,LANDAM
      VolumeR=vdepth/(1+exp((r-radius)/diffuss))*r**(landa+2)
      return
      end   function VolumeR

c-----------------------------------------------------------------------------
      real*8 function SurfaceR(r)
c-----------------------------------------------------------------------------
      IMPLICIT real*8 (A-H,O-Z)
      common/parameters/vdepth,radius,diffuss
      common/armon/LANDA,LANDAM
      SurfaceR = 4./(1+exp((r-radius)/diffuss))**2*
     > vdepth*exp((r-radius)/diffuss)*r**(landa+2)
      return
      end   function SurfaceR

c-----------------------------------------------------------------------------
      real*8 function SpinOrbitR(r)
c-----------------------------------------------------------------------------
      IMPLICIT real*8 (A-H,O-Z)
      common/parameters/vdepth,radius,diffuss
      common/armon/LANDA,LANDAM
C     SpinOrbitR =  -vdepth/(diffuss*radius)*exp((r-radius)/diffuss)
C    >             /  (1+exp((r-radius)/diffuss))**2*r**(landa+2)
      SpinOrbitR =  -vdepth/diffuss*exp((r-radius)/diffuss)
     >             /  (1+exp((r-radius)/diffuss))**2*r**(landa+1)
      return
      end   function SpinOrbitR



      END MODULE om_retrieve
