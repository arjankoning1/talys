      subroutine geftalys(ZCN,ACN,NE,excit,fsig,iprint,NEVT)
c
c +---------------------------------------------------------------------
c | Author: Vasily Simutkin and Michail Onegin
c | Date  : December 18, 2013
c | Task  : FORTRAN version of GEF: Fission yields and neutrons
c +---------------------------------------------------------------------
c
        include "gef.cmb"
         real*4 JFraglight,JFragheavy,Irigidspher,IfragEff
         real*4 Irigid,Ieff,Jrms
         real*4 JFRAGpre(200,150,0:300)
         real*4 Lymass
         real*8 WN(20:320),sq2pi,WNT(20:320),WN2P(30:320)
         real*8 WN0(20:320),WN1(20:320),WN2(20:320),WN3(20:320),
     ;          WN11(20:320),WN22(20:320)
      real*8 rani,ranj,rijk,ranb,rans,r
      integer AUCD,EMODE,ZT,AT
c     common /outgef/ysum(350),ysump(350),yAZ(200,150),yAZp(200,150),
c    ; anpre_sum(350),anpost_sum(350),ann_sum(50),anMean_sum
      common /rand/rani,ranj,rijk,ranb,rans,nrand
         integer PaWidthS2,PaWidtS2
         logical BSUB
         external Erf
         common /comvar/ZUCD
         common /cgamma/Egamma(1000),ArrayEn(100),ArrayTn(100)
         dimension Beta(0:5,2,150),Edefo(0:3,2,150),Zmean(0:3,2,350),
     ;  Zshift(0:3,2,350),Temp(0:3,2,350),TempFF(0:3,2,350),
     ;  Eshell(0:3,2,350),PEOZ(0:5,2,350),PEON(0:5,2,350),
     ;  EPART(0:5,2,350),SpinRMSNZ(0:5,2,200,150)
         dimension APRE(350),AMPRE(0:5,350)
         dimension ANZPRE(200,150),ANPRE(200),ANMPRE(0:5,200),
     ;   ZISOPRE(350,150),ZISOPOST(350,150)
        dimension APOST(350),AMPOST(0:5,350),ZPOST(200),ZMPOST(0:5,200)
        dimension ANPOST(200),ANMPOST(0:5,200),ANZPOST(200,150)
        dimension AN2dpost(350,0:20),AN2dpre(350,0:20),
     ;  ANApre(350),ANApost(350),aNN(0:50)
         dimension yAZ_temp(350,150),yAZp_temp(350,150)
         dimension excit(1000), fsig(1000)
        dimension Ecollsaddle(0:5),Edefo2d(350,800),Ecoll2d(350,500),
     ;  Eintr2d(350,1000),Eexc2d(350,1000),AEkin(200,300)
         Data sq2pi/0.3989422804D0/,pi/3.14159265/
         Data PZMeanS1/52.8/,PZMeanS2/55/,PZMeanS3/65/
         do i=0,5
           do j=1,2
             do k=1,150
               Beta(i,j,k)=0.
               do im=1,200
                 SpinRMSNZ(i,j,im,k)=0.
               enddo
             enddo
           enddo
         enddo
         do i=0,3
           do j=1,2
             do k=1,150
               Edefo(i,j,k)=0.
             enddo
             do k=1,350
               Zshift(i,j,k)=0.
             enddo
           enddo
         enddo

c         BSUB=.FALSE.
         BSUB=.TRUE.
         MODEC=2
c         NEVT=50000
c   MODEC=1 - analitical calculations
c   MODEC=2 - Monte Carlo calculations

           nrand=0
      r=5.d0**19
      call setijk(r)
      call setran
         do i=1,350
          YSUM(i)=0.0
          YSUMp(i)=0.0
          anpre_sum(i)=0.0
          anpost_sum(i)=0.0
        end do

         do i=0,50
          ann_sum(i)=0.0
         end do

         anMean_sum=0.0

c   27.09.2013
        do i=1,350
          do j=1,150
            yAZ(i,j)=0.0
            yAZp(i,j)=0.0
          end do
         end do
c   end 27.09.2013

         do i=1,350
          apre(i)=0.0
          apost(i)=0.0
           do im=0,5
            ampre(im,i)=0.0
            ampost(im,i)=0.0
            do j=1,2
              PEOZ(im,j,i)=0.0
              PEON(im,j,i)=0.0
              EPART(im,j,i)=0.0
            end do
           end do
           do j=1,800
            Edefo2d(i,j)=0.0
           end do
           do j=1,150
            ZISOPRE(i,j)=0.0
            ZISOPOST(i,j)=0.0
           end do
           do j=1,500
            Ecoll2d(i,j)=0.0
           end do
           do j=1,1000
            Eintr2d(i,j)=0.0
            Eexc2d(i,j)=0.0
            Egamma(j)=0.0
           end do

         end do

         do i=1,200
          anpre(i)=0.0
          anpost(i)=0.0
          ZPOST(i)=0.0
          APOST(i)=0.0
           do j=1,150
            anzpre(i,j)=0.0
            anzpost(i,j)=0.0
             do k=0,300
              JFRAGpre(i,j,k)=0.0
             end do
           end do
           do j=1,300
            AEkin(i,j)=0.0
           end do
           do im=0,5
            anmpre(im,i)=0.0
            anmpost(im,i)=0.0
            ZMPOST(im,i)=0.0
           end do
         end do

         do i=1,350
          do n=0,20
           AN2dpost(i,n)=0.0
           an2dpre(i,n)=0.0
          end do
         end do

         do i=0,50
          ann(i)=0.0
         end do
         if (ZCN.gt.136.) return
         if (ACN.gt.350.) return
         if (ACN-ZCN.gt.203.) return

         Emaxvalid = 100
c   /' Maximum allowed excitation energy '/
         PDZMeanS1 = -0.1
         PZCurvS1 = 0.75
         PDZMeanS2 = -0.45
         PZCurvS2 = 0.5
         PAWidtS2 = 7
         PAWidthS2 = 14
c   /' A width of Mode 2 (box) '/
         PDZMeanS3 = 0
c   /' Shift of mean Z of Mode 3 '/
         PZCurvS3 = 0.22
c   /' Curvature in Z of Mode 3 '/
         PShellS1 = -2.8
c   /' Shell effect for Mode 1 '/
         PZS1olapShell = -0.5
         PZS1olappos = 44.85
c   /' Shell in light fragment for enhanced S1 around Pu '/
         PZS1olapcurv = 0.1
c   /' Shell in light fragment for enhanced S1 around Pu '/
         PShellS2 = -4.0
c    /' Shell effect for Mode 2 '/
         PShellS3 = -5.8
c   /' Shell effect for Mode 3 '/
         PZS3olappos = 37
c   /' Pos. of S3 shell in light fragment (in N!) '/
         PZS3olapcurv = 0.005
c   /' for width of S3 shell in light fragment '/
         DeltaS0 = 0
c       /' Shell effect for SL, for individual systems '/
         TlowS1 = 0.30
         TlowS2 = 0.31
c    /' Slope parameter for tunneling '/
         TlowS3 = 0.32
c   /' Slope parameter for tunneling '/
         TlowSL = 0.31
c    /' Slope parameter for tunneling '/
         TlowS11 = 0.36
c   /' Slope parameter for tunneling '/
         TlowS1offset = 6.7
c   /' Offset energy for variation of TlowS1 '/
         TlowS1slope = 0
c     /' not used '/
         Pattpol = 0.055
c   /' Attenuation of 132Sn shell '/
         dEDefoS1 = -3
c      /' Deformation energy expense for Mode 1 '/
         dEDefoS2 = 0
c    /' Deformation energy expense for Mode 2 '/
         dEDefoS3 = 0
c      /' Deformation energy expense for Mode 3 '/
         betaL0 = 26.6
c   /' Offset for deformation of light fragment '/
         betaL1 = 0.80
c   /' Slope for deformation of light fragment '/
         betaH0 = 48.0
c     /' Offset for deformation of heavy fragment '/
         betaH1 = 0.70
c    /' Slope for deformation of heavy fragment '/
         kappa = 0
c    /' N/Z dedendence of A-asym. potential '/
         TCOLLFRAC = 0.13
c     /' Tcoll per energy gain from saddle to scission '/
         ECOLLFRAC = 0.0
c     /' Ecoll per energy gain from saddle to scission '/
         TFCOLL = 0.1
c       /' Tcoll per energy above barrier '/
         TCOLLMIN = 0.12
c     /' Minimum coll. temperature '/
         ESHIFTSASCI = -30
c   /' Shift of saddle-scission energy '/
         EDISSFRAC = 0.3
         SIGDEFO = 0.165
         EexcSIGrel = 0.5
c   /' Relative sigma of coll. and intr. energy '/
         DNECK = 3.
c   /' Tip distance at scission / fm '/
         FTRUNC50 = 1
c       /' Truncation near Z = 50 '/
         ZTRUNC50 = 50
c      /' Z value for truncation '/
         FTRUNC28 = 0.56
c      /' Truncation near Z = 28 '/
         ZTRUNC28 = 30.5
c     /' Z value for truncation '/
         ZMAXS2 = 60
c        /' Maximum Z of S2 channel in light fragment '/
         NTRANSFEREO = 6
c    /' Steps for E sorting for even-odd effect '/
         NTRANSFERE = 12
c     /' Steps for E sorting for energy division '/
         Csort = 0.1
c         /' Smoothing of energy sorting '/
         PZEOsymm = 2.25
c   /' Even-odd effect in Z at symmetry '/
         PNEOSymm = 0.5
c   /' Even-odd effect in N at symmetry '/
         REOTHRESH = 0.04
c  /' Threshold for asymmetry-driven even-odd effect'/
         REOSIGMA = 0.35
         REOMAX = 0.40
c    /' Maximum even-odd effect '/
         POLARadd = 0.3
c   /' Offset for enhanced polarization '/
         POLARfac = 1
c   /' Enhancement of polarization of ligu. drop '/
         TPOLRED = 0.01
c    /' Reduction of temperature for sigma(Z) '/
         HOMPOL = 2.0
c   /' hbar omega of polarization oscillation '/
         ZPOL1 = 0
c         /' Extra charge polarization of S1 '/
         ETHRESHSUPPS1 = 0.5
         ESIGSUPPS1 = 0.3
c    /' Sigma of suppression threshold '/
         Pnx = 0
c   /' Enhanced inverse neutron x section '/
         Tscale = 0.85
         Econd = 2
         Emode = 0
c    /' 0: E over BFB; 1: E over gs; 2: Eneutron '/
         Jscaling = 1
c    /' General scaling of fragment angular momenta '/
         Spinodd = 0.4
c    /' RMS Spin enhancement for odd Z '/

c      start of calculations

c         ZCN=92.
c         ACN=236.
c         ENCM=14.
c         ENCM=0.0253E-6
         IEiso=0
         EMODE=0

         RNCN=ACN-ZCN

c   /' Shell effects for the symmetric fission channel '/

        DeltaS0 = 0
       If (NInt(ZCN).eq. 95.And.NInt(ACN).eq.242) Then
          DeltaS0 = -0.1
       end if
       If (NInt(ZCN).eq.95.And.NInt(ACN).eq.244) Then
          DeltaS0 = -0.1
       end if
       If (NInt(ZCN).eq.96.And.NInt(ACN).eq.244) Then
          DeltaS0 = 0.1
       end if
       If (NInt(ZCN).eq.93.And.NInt(ACN).eq.238) Then
          DeltaS0 = 0.2
       end if
       If (NInt(ZCN).eq.94.And.NInt(ACN).eq.240) Then
          DeltaS0 = -0.1
       end if
       If (NInt(ZCN).eq.94.And.NInt(ACN).eq.241) Then
          DeltaS0 = -0.35
       end if
       If (NInt(ZCN).eq.94.And.NInt(ACN).eq.243) Then
          DeltaS0 = -0.3
       end if
       If (NInt(ZCN).eq.90.And.NInt(ACN).eq.228) Then
          DeltaS0 = 0.65
       end if
       If (NInt(ZCN).eq.90.And.NInt(ACN).eq.230) Then
          DeltaS0 = 0.65
       end if
       If (NInt(ZCN).eq.90.And.NInt(ACN).eq.233) Then
          DeltaS0 = -0.2
       end if
       If (NInt(ZCN).eq.92.And.NInt(ACN).eq.233) Then
          DeltaS0 = 0.65
       end if
       If (NInt(ZCN).eq.92.And.NInt(ACN).eq.234) Then
          DeltaS0 = 0.5
       end if
       If (NInt(ZCN).eq.92.And.NInt(ACN).eq.235) Then
          DeltaS0 = 0.35
       end if
       If (NInt(ZCN).eq.92.And.NInt(ACN).eq.236) Then
          DeltaS0 = 0.3
       end if


c    /' Mean deformation as a function of mass '/
c    /' Mode 0: liquid drop '/
      beta1prev = 0.3
      beta2prev = 0.3
      beta1opt = beta1prev
      beta2opt = beta2prev
      Do I = 10 , nInt(ZCN) - 10
       Z1 = I
       Z2 = ZCN - Z1
       A1 = Z1 / ZCN * ACN
       A2 = ACN - A1

      call Beta_Equi(A1,A2,Z1,Z2,dneck,beta1prev,beta2prev,
     ;  beta1opt,beta2opt)

      If (Bsub) Then

      call  BetaEquiA(A1,A2,Z1,Z2,dneck,beta1prev,
     ;  beta2prev,beta1opt,beta2opt)

      End If

c       write(*,*) ' Z1=',I,' beta1=',beta1opt

      Beta(0,1,I) = beta1opt
c   /' "light" fragment '/
      Beta(0,2,I) = beta1opt
c   /' "heavy" fragment '/
      beta1prev = beta1opt
      beta2prev = beta2opt
      Edefo1 = LyMass(Z1,A1,beta1opt) - LyMass(Z1,A1,0.)
      Edefo(0,1,I) = Edefo1
c   /' "light" fragment '/
      Edefo(0,2,I) = Edefo1
c   /' "heavy" fragment '/
      End do

c    /' Mode 1: deformed shells (light) and spherical (heavy) '/
       Do I = 10 ,  nInt(ZCN) - 10
      Z1 = I
      Z2 = ZCN - Z1
      A1 = (Z1 - ZPOL1) / ZCN * ACN
c   /' polarization roughly considered '/
      A2 = ACN - A1
c    ' Betaoptlight(A1,A2,Z1,Z2,dneck,0,rbetald)
c      /' numean of Cf requires shells in the light fragment: '/
      rbetashell = betalight(Real(I))
c    ' rbeta = 0.5 * (rbetald + rbetashell)
      rbeta = rbetashell
      Beta(1,1,I) = rbeta
c  /' "light" fragment '/
      Edefo1 = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.)
      Edefo(1,1,I) = Edefo1
c   /' "light" fragment '/
       End do

      Do I = 10 , nInt(ZCN) - 10
      rbeta = 0
      Beta(1,2,I) = rbeta
c  /' "heavy" fragment '/
      Edefo(1,2,I) = 0
c    /' "heavy" fragment '/
      End do

c    /' Mode 2: deformed shells (light and heavy) '/
      Do I = 10 , nInt(ZCN) - 10
      Z1 = I
      Z2 = ZCN - Z1
      A1 = (Z1 - 0.5E0) / ZCN * ACN
c   /' polarization roughly considered '/
      A2 = ACN - A1
c    ' Betaoptlight(A1,A2,Z1,Z2,dneck,betaheavy(Z2),rbetald)
      rbetashell = betalight(Real(I))
c    ' rbeta = 0.5 * (rbetald + rbetashell)
      rbeta = rbetashell
      Beta(2,1,I) = rbeta
      Edefo1 = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.)
      Edefo(2,1,I) = Edefo1
      End do

      Do I = 10 ,nInt(ZCN) - 10
      rbeta = betaheavy(Real(I))
      Beta(2,2,I) = rbeta
      Z1 = I
      A1 = (Z1 + 0.5E0) / ZCN * ACN
c   /' polarization rougly considered '/
      Edefo1 = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.)
      Edefo(2,2,I) = Edefo1
      End do

c    /' Mode 3 '/
      Do I = 10 , nInt(ZCN) - 10
      Z1 = I
      Z2 = ZCN - Z1
      A1 = (Z1 - 0.5E0) / ZCN * ACN
c  /' polarization roughly considered '/
      A2 = ACN - A1
      rbeta = betalight(Real(I))
c  /'  rbeta = 0 '/
c  /'  Call Betaoptlight(A1,A2,Z1,Z2,dneck,betaheavy(Z2),rbeta) '/
      Beta(3,1,I) = rbeta
      Edefo1 = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.)
      Edefo(3,1,I) = Edefo1
      End do

      Do I = 10 , nInt(ZCN) - 10
      rbeta = betaheavy(Real(I))
c   /' beta = 0.5 '/
      Beta(3,2,I) = rbeta
      Z1 = I
      A1 = (Z1 + 0.5E0) / ZCN * ACN
c  /' polarization rougly considered '/
      Edefo1 = LyMass(Z1,A1,rbeta) - LyMass(Z1,A1,0.)
      Edefo(3,2,I) = Edefo1
      End do

c    /' Mode 4: (Channel ST1 in both fragments) '/
      Do I = 10 , nInt(ZCN) - 10
      Z1 = I
      Z2 = ZCN - Z1
      rbeta = Beta(1,2,I)
      if( rbeta .lt. 0) Then
       rbeta = 0
      end if
      Beta(4,1,I) = rbeta
      Beta(4,2,I) = rbeta
      End do

c    /' Mode 5: (Channel ST2 in both fragments) '/
      Do I = 10 , nInt(ZCN) - 10
      Z1 = I
      Z2 = ZCN - Z1
      rbeta = Beta(2,2,I)
      if (rbeta .lt. 0) Then
       rbeta = 0
      end if
      Beta(5,1,nInt(Z1)) = rbeta
      Beta(5,2,nInt(Z1)) = rbeta
      end do

c    /' Mean Z as a function of mass '/

c    /' Mode 0 '/
      Do I = 10 , nInt(ACN) - 10
      ZUCD = I / ACN * ZCN
      beta1 = Beta(0,1,nInt(ZUCD))
      beta2 = Beta(0,2,nInt(ZCN - ZUCD))
c       write(*,*) i,' beta=',beta1,beta2
      Z1 = Zequi(ZCN,Real(I), ACN - I, beta1, beta2, dneck,0)
c       write(*,*) ' Z1=',Z1
      Zmean(0,1,I) = Z1
      Zmean(0,2,nInt(ACN - I)) = ZCN - Z1
      Zshift(0,1,I) = Z1 - ZUCD
      Zshift(0,2,nInt(ACN - I)) = ZUCD - Z1
c      If(nInt(ACN-I).eq.119) Then
c       write(*,*) ' ZUCD,Z1=',ZUCD,Z1
c      end if

      end do

c      write(*,*) ' Zshift=', ZShift(0,2,119)

c    /' Mode 1 '/
      Do I = 10 , nInt(ACN) - 10
      ZUCD = I / ACN * ZCN
      Z = ZUCD + ZPOL1
c /' Charge polarisation is considered in a crude way '/
      beta1 = Beta(1,1,nInt(Z))
c  /' "light" fragment '/
      Z = ZUCD - ZPOL1
      beta2 = Beta(1,2,nInt(ZCN - Z))
c  /' "heavy" fragment '/
      Z1 = Zequi(ZCN,Real(I), ACN - I, beta1, beta2, dneck,1)
c       write(*,*) i,beta1,beta2,Z1
      Z1 = Z1 + ZPOL1
c  /' Charge polarization by shell '/

      If (ZCN - Z1 .lt. 50) Then
        Z1 = ZCN - 50
c  /' Z of mean heavy fragment not below 50 '/
      EndIf

      Zmean(1,1,I) = Z1
      Zmean(1,2,nInt(ACN - I)) = ZCN - Z1
      Zshift(1,1,I) = Z1 - ZUCD
      Zshift(1,2,nInt(ACN - I)) = ZUCD - Z1
      end do


c    /' Mode 2 '/
      Do I = 10 , nInt(ACN) - 10
      ZUCD = I / ACN * ZCN
      Z = ZUCD
c  /' Charge polarisation is here neglected '/
      beta1 = Beta(2,1,nInt(Z))
      beta2 = Beta(2,2,nInt(ZCN - Z))
      Z1 = Zequi(ZCN,Real(I), ACN-I, beta1, beta2, dneck,2)
      Zmean(2,1,I) = Z1
      Zshift(2,1,I) = Z1 - ZUCD
      Zmean(2,2,nInt(ACN - I)) = ZCN - Z1
      Zshift(2,2,nInt(ACN - I)) = ZUCD - Z1
      end do

c      write(*,*) ' Zshift=',(Zshift(2,2,I),i=1,350)

c    /' Mode 3 '/
      Do I = 10 , nInt(ACN) - 10
      ZUCD = I / ACN * ZCN
      Z = ZUCD
c  /' Charge polarisation is here neglected '/
      beta1 = Beta(3,1,nInt(Z))
      beta2 = Beta(3,2,nInt(ZCN - Z))
      Z1 = Zequi(ZCN,Real(I), ACN - I, beta1, beta2, dneck,3)
      Zmean(3,1,I) = Z1
      Zshift(3,1,I) = Z1 - ZUCD
      Zmean(3,2,nInt(ACN - I)) = ZCN - Z1
      Zshift(3,2,nInt(ACN - I)) = ZUCD - Z1
      end do


c    /' Central Z values of fission modes '/

c    /' Fit to positions of fission channels (Boeckstiegel et al., 2008) '/
c    /' PDZMeanS1 and PDZMeanS2 allow for slight adjustments '/
       ZCMode1 = (52.9E0 - 51.5E0) / (1.56E0 - 1.50E0) *
     ;            (ZCN**1.3E0 / ACN - 1.50E0) + 51.5E0 + PDZMeanS1
       ZCMode2 = (55.8E0 - 54.5E0) / (1.56E0 - 1.50E0) *
     ;            (ZCN**1.3E0 / ACN - 1.50E0) + 54.5E0 + PDZMeanS2
       ZCMode3 = ZCMode2 + 4.5E0 + PDZMeanS3
       ZCMode0 = ZCN * 0.5E0
c      /' Central Z value of SL mode '/

       PZMeanS1 = ZCMode1
       PZMeanS2 = ZCMode2
       PZMeanS3 = ZCMode3

c    /' General relations between Z and A of fission channels '/
      RZpol = 0
      Do I = 1 , 3
      RA = (ZCMode0 - RZPol) * ACN / ZCN
      RZpol = Zshift(0,2,nInt(RA))

c       write(*,*) ' RA,RZPol',RA,RZPol

      end do

      ACMode0 = (ZCMode0 - RZPol) * ACN / ZCN

c       write(*,*) ' ZC=',ZCMode0,RZPol,ACN,ZCN
c       write(*,*) ' A0=',ACMode0

c  /' mean position in mass '/
      NCMode0 = ACMode0 - ZCMode0

      RZpol = 0
      do I = 1 , 3
      RA = (ZCMode1 - RZPol) * ACN / ZCN
      RZpol = Zshift(1,2,nInt(RA))

c       write(*,*) ' RA,RZPol',RA,RZPol

      end do
c       write(*,*) ' ZC=',ZCMode1,RZPol

      ACMode1 = (ZCMode1 - RZPol) * ACN / ZCN
      NCMode1 = ACMode1 - ZCMode1

      RZpol = 0
      do I = 1 , 3
      RA = (ZCMode2 - RZPol) * ACN / ZCN
      RZpol = Zshift(2,2,nInt(RA))
      end do
c       write(*,*) ' ZC=',ZCMode2,RZPol

      ACMode2 = (ZCMode2 - RZPol) * ACN / ZCN
      NCMode2 = ACMode2 - ZCMode2

      RZpol = 0
      do I = 1 , 3
      RA = (ZCMode3 - RZPol) * ACN / ZCN
      RZpol = Zshift(3,2,nInt(RA))
      end do

      ACMode3 = (ZCMode3 - RZPol) * ACN / ZCN
      NCMode3 = ACMode3 - ZCMode3


c    /' Potential curvatures of fission modes '/

       RI = (RNCN - ZCN)/ACN
c   /' measure of neutron excess '/
       RZCurvS0 = 8.E0 / ZCN**2 * AMasscurv(ZCN, ACN, RI)
       RACurvS0 = 8.E0 / ACN**2 * AMasscurv(ZCN, ACN, RI)

c        temp=Masscurv(ZCN, ACN, RI)

c       write(*,*) ' Mc=',temp

         SN=UMass(ZCN,ACN-1.)+Aypair(ZCN,ACN-1.)-
     ;   (UMass(ZCN,ACN)+Aypair(ZCN,ACN))
c         write(*,*) ' UMass=',UMass(ZCN,ACN-1.)


c    ---the version without the loop over the energy range---
c         Eex=ENCM+SN-BFTFB(ZCN,ACN,1)
c         write(*,*) ' Eex=',Eex

c    ---the new version with the energy loop---
c      write(*,*) 'test1'

       do i=1,350
         ysum(i)=0
         ysump(i)=0
       end do

       do 2000 nex=1,ne

         do i=1,350
          apre(i)=0.0
          apost(i)=0.0
           do im=0,5
            ampre(im,i)=0.0
            ampost(im,i)=0.0
           end do
           do j=1,800
            Edefo2d(i,j)=0.0
           end do
           do j=1,500
            Ecoll2d(i,j)=0.0
           end do
           do j=1,1000
            Eintr2d(i,j)=0.0
            Eexc2d(i,j)=0.0
           end do

         end do

         do i=1,200
          anpre(i)=0.0
          anpost(i)=0.0
           do j=1,150
            anzpre(i,j)=0.0
            anzpost(i,j)=0.0
             do k=0,300
              JFRAGpre(i,j,k)=0.0
             end do
           end do
           do j=1,300
            AEkin(i,j)=0.0
           end do
           do im=0,5
            anmpre(im,i)=0.0
            anmpost(im,i)=0.0
           end do
         end do

         do i=1,350
          do n=0,20
           AN2dpost(i,n)=0.0
           an2dpre(i,n)=0.0
          end do
         end do

         do i=0,50
          ann(i)=0.0
         end do

c   27.09.2013
        do i=1,350
          do j=1,150
            yAZ_temp(i,j)=0.0
            yAZp_temp(i,j)=0.0
          end do
         end do
c   end 27.09.2013

       if (iprint.ne.0) then
       write(*,*) 'Detailed GEF output'
       write(*,*)
       write(*,*) 'Before: Excit(', nex, ')=', Excit(nex)
       endif

        PEexc=Excit(nex)
        REexcGS=Excit(nex)

       Eex=Excit(nex)-BFTFB(ZCN,ACN,1)
       if (iprint.ne.0) then
       write(*,*) 'After: Eex =', Eex
       endif

         EexcS0prov=Eex

c    /' Additional influence of N=82 assumed '/
       DeltaNPol = 82.E0 - 50.E0/ZCN * RNCN
       RShellS1eff = PShellS1 * (1.E0 - PAttPol * Abs(DeltaNPol))

c        write(*,*) ' RShellS1Eff=',RShellS1eff

c    /' In Pu, the Z=50 shell meets Z=44 in the light fragment. '/
c    /' A deformed shell at Z=44 is assumed to explain the enhancement
c       of the S1 channel around Pu '/
       S1enhance = PZS1olapShell +
     ;         (ZCN - 50 - PZS1olapPos)**2 * PZS1olapcurv
c        write(*,*) ' S1enhance=',s1enhance
       If(S1enhance.gt.0.) Then
        S1enhance=0.
       endif
       RShellS1eff = RShellS1eff + S1enhance

c        write(*,*) ' RShellS1Eff=',RShellS1eff
       RShellS3eff = PShellS3 * (1.E0 - PZS3olapcurv
     ;   * (ZCN - 60.5E0 - PZS3olappos)**2)
       RShellS3eff = Min(RShellS3eff,0.)

       EldS1 = RACurvS0 * (ACN/ZCN*(PZMeanS1 - ZCMODE0) )**2
       BS1 = EldS1 + RShellS1eff
       EexcS1prov = EExcS0prov - BS1

c       write(*,*)  RACurvS0,PZMeanS1,ZCMODE0

       EldS2 = RACurvS0 * (ACN/ZCN*(PZMeanS2 - ZCMODE0) )**2
       BS2 = EldS2 + PShellS2
       EexcS2prov = EExcS0prov - BS2

       EldS3 = RACurvS0 * (ACN/ZCN*(PZMeanS3 - ZCMODE0) )**2
       BS3 = EldS3 + RShellS3eff
       EexcS3prov = EExcS0prov - BS3

c    /' Mode 11 (overlap of channel 1 in light and heavy fragment '/
c    /' Potential depth with respect to liquid-drop potential: BS11 '/
       BS11 = 2.E0 * (RShellS1eff
     ;        + PZCurvS1 * (ZCMode1 - ZCMode0)**2 )
       If (BS11.gt.RShellS1eff.And.BS11.lt.RShellS1eff + 7.0E0) Then
         BS11 = Min(BS11,RShellS1eff)
       End If
c    /' If the ridge between the two valleys is less than 7 MeV, it is assumed
c       that it acts like one broad valley. Note that the calculated
c       ridge is sharp and rather thin! '/

c        write(*,*) ' BS1=',BS1

c    /' Lowering of effective barrier by lower ZPM due to larger width in
c       partial overlap region (shells in light and heavy fragment) '/
       DES11ZPM = - 0.2E0 * (2.E0 * Abs(ZCMode1 - ZCMode0) )
       DES11ZPM = Min(0.,DES11ZPM)
       BS11 = BS11 + DES11ZPM

       EexcS11prov = EExcS0prov - BS11

c    /' Mode 22 (overlap of channel 2 in light and heavy fragment '/
c    /' Potential depth with respect to liquid-drop potential: BS22 '/

       BS22 = 2.E0 * (EldS2 + PShellS2)
     ;   + 2.E0 * PZCurvS2 * (ZCMode2 - ZCMode0)**2
c  /' Parabola '/

       EexcS22prov = EExcS0prov - BS22

      EMinBarr = Min(0.,BS1)
      EMinBarr = Min(EMinBarr,BS2)
      EMinBarr = Min(EMinBarr,BS3)
      EMinBarr = Min(EMinBarr,BS11)
      EMinBarr = Min(EMinBarr,BS22)


      EexcS0 = EexcS0prov + EMinBarr - DeltaS0
      EexcS1 = EexcS1prov + EMinBarr
      EexcS2 = EexcS2prov + EMinBarr
      EexcS3 = EexcS3prov + EMinBarr
      EexcS11 = EexcS11prov + EMinBarr
      EexcS22 = EexcS22prov + EMinBarr

c       write(*,*) ' EexcS0=',EexcS0
c       write(*,*) ' EexcS1=',EexcS1
c       write(*,*) ' EexcS2=',EexcS2
c       write(*,*) ' EexcS3=',EexcS3
c       write(*,*) ' EexcS11=',EexcS11
c       write(*,*) ' EexcS22=',EexcS22


c    /' Energy above the lowest fission saddle '/
      EexcBarr = Max(EExcS0,EExcS1)
      EexcBarr = Max(EexcBarr,EExcS2)
      EexcBarr = Max(EexcBarr,EExcS3)
      EexcBarr = Max(EexcBarr,EexcS11)
      EexcBarr = Max(EexcBarr,EexcS22)

c    /' Collective temperature used for calculating the widths
c       in mass asymmetry and charge polarization '/

      If(EExcS0.lt.0) Then
       Etunn = -EExcS0
      Else
       Etunn = 0
      end if
       REexceff = Max(0.1,EExcS0)
c  '  TCollMode0 = TFCOLL * REexceff +   /' empirical, replaced by TRusanov '/
      TCollMode0 = TCOLLFRAC * (DeSaddleScission(ZCN**2
     ;  / ACN**0.33333E0) - Etunn)
      TCollMode0 = Max(TCollMode0,0.)

c ' Print "TColl ";DeSaddleScission(PZCN^2/PACN^0.3333),Etunn,TCollMode0

c    /' Temperature description fitting to the empirical systematics of Rusanov et al. '/
c    /' Here from Ye. N. Gruzintsev et al., Z. Phys. A 323 (1986) 307 '/
c    /' Empirical description of the nuclear temperature according to the '/
c    /' Fermi-gas description. Should be valid at higher excitation energies '/
      TRusanov = TRusanovf(REexceff,ACN)
c  '  Print "Temperatures, (GEF, Total, Rusanov): ", TCollMode0, TFCOLL * REexceff, TRusanov
      TCollMode0 = Max(TCollMode0,TRusanov)
c    /' Transition vom const. temp. to Fermi gas occurs around 20 MeV by MAX function '/
c  '   Print "Effective: ",TCollMode0

      TPolMode0 = TPolRed * TCollMode0
      TAsymMode0 = Sqrt(TCollMode0**2 + TCOLLMIN**2)

      Epotscission = (DeSaddleScission(ZCN**2 / ACN**0.33333E0) - Etunn)

c    /' Suppression of S1 fission channel due to reduced pairing in 132Sn '/
c    /' At very low excitation energy on the fission path, the binding energy at the
c       S1 fission channel does not profit as much from pairing as SL and S2,
c       because pairing is reduced in magic nuclei. This leads to a reduction of
c       the yield in S1 in the case that the fully paired ground-state configuration
c       is populated on the fission path with a considerable probability. '/
       EeffS2 = Max(EexcS2,0.) + EDISSFRAC * Epotscission - 2.3E0
       EeffS2 = Max(0.,EeffS2)
c       /' -2.3 MeV, because fission channels are assumed to be chosen before scission '/
       If(EeffS2.lt. ETHRESHSUPPS1 + 2.E0 * ESIGSUPPS1) Then
       EexcS1 = EexcS1 -
     ;    0.5E0 * 4.E0 * 12.E0 / Sqrt(132.E0)
     ;  * Gaussintegral(ETHRESHSUPPS1 - EeffS2,ESIGSUPPS1)
       End If

c   /' The relative yield of S1 becomes smaller for spontaneous fission of heavier nuclei. '/
c   /' The physical reason is not clear. '/

      TlowS1used = TlowS1 - TlowS1slope
     ;   * (EDISSFRAC * Epotscission - 3.E0-TlowS1offset)

      TCollMode1 = TFCOLL * EexcS1 +
     ; TCOLLFRAC * (DeSaddleScission(ZCN**2 / ACN**0.33333E0) - Etunn)
      TCollMode1 = Max(TCollmode1,0.)
      TPolMode1 = TPolRed * TCollMode1
      TAsymMode1 = Sqrt(TCollMode1**2 + TCOLLMIN**2)

      TCollMode2 = TFCOLL * EexcS2 +
     ; TCOLLFRAC * (DeSaddleScission(ZCN**2 / ACN**0.33333E0) - Etunn)
      TCollMode2 = Max(TCollmode2,0.)
      TPolMode2 = TPolRed * TCollMode2
      TAsymMode2 = Sqrt(TCollMode2**2 + TCOLLMIN**2)

      TCollMode3 = TFCOLL * EexcS3 +
     ; TCOLLFRAC * (DeSaddleScission(ZCN**2 / ACN**0.33333E0) - Etunn)
      TCollMode3 = Max(TCollmode3,0.)
      TPolMode3 = TPolRed * TCollMode3
      TAsymMode3 = Sqrt(TCollMode3**2 + TCOLLMIN**2)

c    /' Stiffness in polarization '/

      RZ = ZCN * 0.5E0
      RA = ACN * 0.5E0
      beta1 = Beta(0,1,nInt(RZ+0.5))
      beta2 = Beta(0,2,nInt(RZ+0.5))
      RPolCurvS0 = ( Lymass( RZ - 1.E0, RA, beta1 ) +
     ;        Lymass( RZ + 1.0E0, RA, beta2 ) +
     ;        Lymass( RZ + 1.0E0, RA, beta1 ) +
     ;        Lymass( RZ - 1.0E0, RA, beta2 ) +
     ;        ecoul( RZ - 1.0E0, RA, beta1,
     ;               RZ + 1.0E0, RA, beta2, dneck) +
     ;        ecoul( RZ + 1.0E0, RA, beta1,
     ;               RZ - 1.0E0, RA, beta2, dneck) -
     ;    2.0E0*ecoul( RZ, RA, beta1, RZ, RA, beta2, dneck) -
     ;    2.0E0*Lymass( RZ, RA, beta1 ) -
     ;    2.0E0*Lymass( RZ, RA, beta2) ) * 0.5E0

      PPolCurvS0 = RPolCurvS0

      RPolCurvS1 = RPolCurvS0
      RPolCurvS2 = RPolCurvS0
      RPolCUrvS3 = RPolCurvS0

c    /' Mean values and standard deviations of fission modes '/

c    Dim As Single REintr_S1, R_E_intr_S2, R_E_intr_S3   ' intrinsic exc. energies at barrier
c    Dim As Single R_Att_S1, R_Att_S2, R_Att_S3            ' attenuation of shell
c    Dim As Single E_backshift
      Ebackshift = -3.

      SIGZMode0 = Sqrt(0.5E0 * TAsymMode0/RZCurvS0)
      expo=HOMPOL/(2.E0 * TPolMode0)
      If (TCollMode0 .gt. 1.E-2.and.abs(expo).lt.80.) Then
      SigPolMode0 = Sqrt(0.25E0 * HOMPOL / RPolCurvS0 /
     ;                fTanh(expo))
      Else
      SigPolMode0 = Sqrt(0.25E0 * HOMPOL / RPolCurvS0)
c        /' including influence of zero-point motion '/
      Endif

      REintrS1 = Max(EExcS1+Aypair(ZCN,ACN)+EBackshift,0.)
      RAttS1 = exp(-REintrS1/18.5E0)
      SIGZMode1 = Sqrt(0.5E0 * TAsymMode1/(PZCurvS1*RAttS1))
      expo=HOMPOL/(2.E0 * TPolMode1)
      If (TCollMode1 .gt. 1.E-2.and.abs(expo).lt.80.) Then
      SigPolMode1 = Sqrt(0.25E0 * HOMPOL / RPolCurvS1 /
     ;                fTanh(expo))
      Else
      SigPolMode1 = Sqrt(0.25E0 * HOMPOL / RPolCurvS1)
      Endif

      REintrS2 = Max(EExcS2+Aypair(ZCN,ACN)+EBackshift,0.)
      RAttS2 = exp(-REintrS2/18.5E0)
      SIGZMode2 = Sqrt(0.5E0 * TAsymMode2/(PZCurvS2*RAttS2))
      expo=HOMPOL/(2.E0 * TPolMode2)
      If (TCollMode2 .gt. 1.E-2.and.abs(expo).lt.80.) Then
      SigPolMode2 = Sqrt(0.25E0 * HOMPOL / RPolCurvS2 /
     ;                fTanh(expo))
      Else
      SigPolMode2 = Sqrt(0.25E0 * HOMPOL / RPolCurvS2)
      Endif

      REintrS3 = Max(EexcS3+Aypair(ZCN,ACN)+EBackshift,0.)
      RAttS3 = exp(-REintrS3/18.5E0)
      SIGZMode3 = Sqrt(0.5E0 * TAsymMode3/(PZCurvS3*RAttS3))
      expo=HOMPOL/(2.E0 * TPolMode3)
      If (TCollMode3 .gt. 1.E-2.and.abs(expo).lt.80.) Then
      SigPolMode3 = Sqrt(0.25E0 * HOMPOL / RPolCurvS3 /
     ;                fTanh(expo))
      Else
      SigPolMode3 = Sqrt(0.25E0 * HOMPOL / RPolCurvS3)
      End if


c    /' Energy-dependent shift of fission channels '/
c    Scope
c      Dim As Single DZ_S1,DZ_S2,DZ_S3

c      DZ_S1 =  ZC_Mode_1 * _
c               (P_Z_Curv_S1*R_Att_S1 / (R_Z_Curv_S0 + P_Z_Curv_S1*R_Att_S1) _
c             - (P_Z_Curv_S1 / (R_Z_Curv_S0 + P_Z_Curv_S1) ) )
c      DZ_S2 =  ZC_Mode_2 * _
c               (P_Z_Curv_S2*R_Att_S2 / (R_Z_Curv_S0 + P_Z_Curv_S2*R_Att_S2) _
c             - (P_Z_Curv_S2 / (R_Z_Curv_S0 + P_Z_Curv_S2) ) )
c      DZ_S3 =  ZC_Mode_3 * _
c               (P_Z_Curv_S3*R_Att_S3 / (R_Z_Curv_S0 + P_Z_Curv_S3*R_Att_S3) _
c             - (P_Z_Curv_S3 / (R_Z_Curv_S0 + P_Z_Curv_S3) ) )


      DZ_S1 = 0.
      DZ_S2 = 0.
      DZ_S3 = 0.


      PZMeanS0 = ZCMode0
      PZMeanS1 = ZCMode1 + DZ_S1
c  /' Copy to global parameter '/
      PZMeanS2 = ZCMode2 + DZ_S2
c  /'             "            '/
      PZMeanS3 = ZCMode3 + DZ_S3
c    End Scope


c    /' General relations between Z and A of fission channels '/
c    /' 2nd iteration '/

      RZpol = 0
      do I = 1 , 3
      RA = (ZCMode0 - RZPol) * ACN / ZCN
      RZpol = Zshift(0,2,nInt(RA))
      end do

      ACMode0 = (ZCMode0 - RZPol) * ACN / ZCN
c   /' mean position in mass '/
      NCMode0 = ACMode0 - ZCMode0

      RZpol = 0
      do I = 1 , 3
      RA = (ZCMode1 - RZPol) * ACN / ZCN
      RZpol = Zshift(1,2,nInt(RA))
      end do
      ACMode1 = (ZCMode1 - RZPol) * ACN / ZCN
      NCMode1 = ACMode1 - ZCMode1

      RZpol = 0
      do I = 1 , 3
      RA = (ZCMode2 - RZPol) * ACN / ZCN
      RZpol = Zshift(2,2,nInt(RA))
      end do
      ACMode2 = (ZCMode2 - RZPol) * ACN / ZCN
      NCMode2 = ACMode2 - ZCMode2

      RZpol = 0
      do I = 1 , 3
      RA = (ZCMode3 - RZPol) * ACN / ZCN
      RZpol = Zshift(3,2,Int(RA))
      end do
      ACMode3 = (ZCMode3 - RZPol) * ACN / ZCN
      NCMode3 = ACMode3 - ZCMode3


c   /' Yields of the fission modes '/

      YieldMode0 = Getyield(EexcS0,EexcS0,TlowSL,TEgidy(ACN,0.E0,0.8E0))
      YieldMode0 = Max(YieldMode0,0.)

      YieldMode1 =
     ; Getyield(EexcS1,EexcS0,TlowS1used,TEgidy(ACN,RShellS1eff
     ; + dEDefoS1,0.8E0))
c                  /'  - Getyield(EexcS0 - EldS1,Tlow,Thigh); '/
      YieldMode1 = Max(YieldMode1,0.)

      YieldMode2 = Getyield(EexcS2,EexcS0,TlowS2,TEgidy(ACN,PShellS2
     ; + dEDefoS2,0.8E0))
c                  /'  - Getyield(EexcS0 - EldS2,Tlow,Thigh); '/
      YieldMode2 = Max(YieldMode2,0.)

      YieldMode3 = Getyield(EexcS3,EexcS0,TlowS3,TEgidy(ACN,PShellS3
     ; + dEDefoS3,0.8E0))
c                  /'  - Getyield(EexcS0 - EldS2,Tlow,Thigh); '/
      YieldMode3 = Max(YieldMode3,0.)

      If(BS11 .gt. 0) Then
       YieldMode11 = 0
      Else
       YieldMode11 = Getyield(EexcS11,EexcS0, TlowS11,
     ;     TEgidy(ACN,RShellS1eff + 2.E0 * dEDefoS1,0.8E0))
      End If

      If(BS22.gt.0) Then
       YieldMode22 = 0
      Else
       YieldMode22 = Getyield(EexcS22,EexcS0, TlowS2,
     ;     TEgidy(ACN,PShellS2,0.8E0))
      End If

      YieldNorm = YieldMode0 + YieldMode1 + YieldMode2 + YieldMode3
     ;              + YieldMode11 + YieldMode22
      YieldMode0 = YieldMode0 / YieldNorm
      YieldMode1 = YieldMode1 / YieldNorm
      YieldMode2 = YieldMode2 / YieldNorm
      YieldMode3 = YieldMode3 / YieldNorm
      YieldMode11 = YieldMode11 / YieldNorm
      YieldMode22 = YieldMode22 / YieldNorm

      if (iprint.ne.0) then
      write(*,*) ' Yield0=',YieldMode0
      write(*,*) ' Yield1=',YieldMode1
      write(*,*) ' Yield2=',YieldMode2
      write(*,*) ' Yield3=',YieldMode3
      write(*,*) ' Yield11=',YieldMode11
      write(*,*) ' Yield22=',YieldMode22
      end if

c    /' Mass widhts of the fission channels '/

      SigAMode0 = SigZMode0 * ACN / ZCN
c  /' width in mass '/
      SigAMode1 = SigZMode1 * ACN / ZCN
      SigAMode2 = SigZMode2 * ACN / ZCN
      SigAMode3 = SigZMode3 * ACN / ZCN
      SigAMode11 = SigZMode1 * sqrt(2.E0) * ACN / ZCN
      SigAMode22 = SigZMode2 * sqrt(2.E0) * ACN / ZCN

c       write(*,*) ' Sig0=', SigPolMode0
c       write(*,*) ' Sig1=', SigPolMode1
c       write(*,*) ' Sig2=', SigPolMode2
c       write(*,*) ' Sig3=', SigPolMode3

c    /' Shell effects of different fission channels '/
c    /' This is the "real" microscopic shell effect, not the effective shell-correction energy '/
c    /' EShell acts on the level density and determines the T parameter '/

      do I = 1 , nInt(ACN) - 1
       do J = 0 , 3
        EShell(J,1,I) = 0
c  /' Shells in "light" fragment assumed to be zero '/
       end do
      DU0 = 0.
      EShell(0,2,I) = 0
c /' Shell = 0 in symmetric mode '/
      DU1 = RShellS1eff + dEDefoS1
c /' + RACurvS1 * (ACMode1 - Float(I,6))**2; '/
      DU1 = MIN(DU1,0.E0)
c /' Technical limit '/
      EShell(1,2,I) = DU1

      DU2 = PShellS2 + dEDefoS2
c /' + RACurvS2 * (ACMode2 - Float(I,6))**2; '/
      DU2 = Min(DU2,0.E0)
c /' Technical limit '/
      EShell(2,2,I) = DU2

      DU3 = RShellS3eff + dEDefoS3
c /' + RACurvS3 * (ACMode3 - Float(I,6))**2; '/
      DU3 = Min(DU3,0.E0)
c /' Technical limit '/
      EShell(3,2,I) = DU3

      end do

c    /' Intrinsic temperatures of fragments at scission '/

c    /' Mean values '/
      TintrMode0 = TEgidy(ACMode0,0.0,0.8)
      TintrMode1heavy = TEgidy(ACMode1,RShellS1eff + dEDefoS1,0.8)
      TintrMode1light = TEgidy(ACN - ACMode1,0.0,0.8)
      TintrMode2heavy = TEgidy(ACMode2,PShellS2 + dEDefoS2,0.8)
      TintrMode2light = TEgidy(ACN - ACMode2,0.0,0.8)
      TintrMode3heavy = TEgidy(ACMode3,RShellS3eff + dEDefoS3,0.8)
      TintrMode3light = TEgidy(ACN - ACMode3,0.0,0.8)

c        write(*,*) ' TintrMode1h=',TintrMode1heavy

c    /' Mass-dependent values of individual fragments '/
c    /' Mode 0 '/
      do IA = 1 , nInt(ACN) - 1
      T = TEgidy(Float(IA),EShell(0,1,IA),0.8)
      Temp(0,1,IA) = T
c /' "light" fragment '/
      T = TEgidy(Float(IA),EShell(0,2,IA),0.8)
      Temp(0,2,IA) = T

      T = TEgidy(Float(IA),0.,1.)
      TempFF(0,1,IA) = T
      TempFF(0,2,IA) = T
      end do

c    /' Mode 1 '/
      do IA = 1 , nInt(ACN) - 1
      T = TEgidy(Float(IA),EShell(1,1,IA),0.8)
      Temp(1,1,IA) = T

c       write(*,*) ' IA=',IA,EShell(1,1,IA)

c /' "light" fragment '/
      T = TEgidy(Float(IA),EShell(1,2,IA),0.8)
      Temp(1,2,IA) = T
c /' "heavy" fragment '/

      T = TEgidy(Float(IA),0.,1.)
      TempFF(1,1,IA) = T
      TempFF(1,2,IA) = T
      end do

c    /' Mode 2 '/
      do IA = 1 , nInt(ACN) - 1
      T = TEgidy(Float(IA),EShell(2,1,IA),0.8)
      Temp(2,1,IA) = T
c /' "light" fragment '/
      T = TEgidy(Float(IA),EShell(2,2,IA),0.8)
      Temp(2,2,IA) = T
c /' "heavy" fragment '/

c   /' The next section is introduced, because energy sorting is not strong enough,
c      when shells are only introduced in the heavy fragment.
c      Ad hoc assumption: For Mode 2 there are shells in both fragments of about
c      equal size. Technically, we neglect the shells in both fragments.
c      This has about the same effect for the energy sorting. '/

       T = TEgidy(Float(IA),0.,0.8)
c  ' FF at scssion
      Temp(2,1,IA) = T
c /' "light" fragment '/
      T = TEgidy(Float(IA),0.,0.8)
c  ' FF at scission
      Temp(2,2,IA) = T
c /' "heavy" fragment '/

      T = TEgidy(Float(IA),0.,1.)
c   ' shell effect neglected
      TempFF(2,1,IA) = T
c ' FFs in their ground state
      TempFF(2,2,IA) = T
c  ' FFs in their ground state
      end do

c    /' Mode 3 '/
      do IA = 1 , nInt(ACN) -1
      T = TEgidy(Float(IA),0.,0.8)
      Temp(3,1,IA) = T
      T = TEgidy(Float(IA),0.,0.8)
      Temp(3,2,IA) = T

      T = TEgidy(Float(IA),0.,1.)
      TempFF(3,1,IA) = T
      TempFF(3,2,IA) = T
      end do

      Do IMode = 0 , 5
      Ecollsaddle(IMode) = 0
      If (IMode .eq. 0) Then
       Etot = EexcS0
      end if
      If (IMode .eq. 1) Then
       Etot = EexcS1
      end if
      If (IMode .eq. 2) Then
       Etot = EexcS2
      end if
      If (IMode .eq. 3) Then
       Etot = EexcS3
      end if
      If (IMode .eq. 4) Then
       Etot = EexcS11
      end if
      If (IMode .eq. 5) Then
       Etot = EexcS22
      end if
      If ((nInt(ZCN)-(nInt(ZCN)/2)*2
     ;  + nInt(RNCN)-(nInt(RNCN)/2)*2) .eq. 0) Then
c  /' Even-even CN '/
        If( Etot .gt. 0. .And. Etot .lt. 2.E0 * 12.E0/SQRT(ACN)) Then
          Ecollsaddle(IMode) = Etot
          Etot = 0
c         /' Excitation below the pairing gap in even-even CN goes into collective excitations '/
        End If
      End If
      Etot = Max(Etot,0.0)
      Etot = Etot + EDISSFRAC * Epotscission
c      /' All excitation energy at saddle and part of the potential-energy gain to scission
c         go into intrinsic excitation energy at scission '/

c       write(*,*) ' Etot=',Etot

       do IA1 = 40 , nInt(ACN) - 40

        IA2 = nInt(ACN) - IA1
        If (IMode .le. 3) Then
          T1 = Temp(IMode,1,IA1)
          T2 = Temp(IMode,2,IA2)
        End If
        If (IMode .eq. 4) Then
          T1 = Temp(1,2,IA1)
          T2 = Temp(1,2,IA2)
        End If
        If (IMode .eq. 5) Then
          T1 = Temp(2,2,IA1)
          T2 = Temp(2,2,IA2)
        End If
        DT = ABS(T2 - T1)

c            write(*,*) ' T1,T2,DT=',T1,T2,DT

c          /' Even-odd effect '/
          expo=Etot/PZEOsymm
          IF (nInt(ZCN) -(nInt(ZCN)/2)*2.eq.0.and.abs(expo).lt.80.) Then
           Rincr1P = Exp(-expo)
        Else
           Rincr1P = 0.0
        End If
          expo=Etot/PNEOsymm
        If (nInt(RNCN)- (nInt(RNCN)/2)*2.eq.0.and.abs(expo).lt.80.) Then
           Rincr1N = Exp(-expo)
        Else
           Rincr1N = 0
        End If
        If (IMode .le. 3) Then
          PEOZ(IMode,1,IA1) = Rincr1P
          PEOZ(IMode,2,IA2) = Rincr1P
          PEON(IMode,1,IA1) = Rincr1N
          PEON(IMode,2,IA2) = Rincr1N
        End If
        If (IMode .eq. 4) Then
          PEOZ(4,2,IA1) = Rincr1P
          PEOZ(4,2,IA2) = Rincr1P
          PEON(4,2,IA1) = Rincr1N
          PEON(4,2,IA2) = Rincr1N
        End If
        If (IMode .eq. 5) Then
          PEOZ(5,2,IA1) = Rincr1P
          PEOZ(5,2,IA2) = Rincr1P
          PEON(5,2,IA1) = Rincr1N
          PEON(5,2,IA2) = Rincr1N
        End If

        Rincr2 = Gaussintegral(DT/Etot-REOThresh,
     ;            REOSigma*(DT+0.0001))
c                  /' even-odd effect due to asymmetry '/

c        write(*,*)  DT/Etot-REOThresh,REOSigma*(DT+0.0001)

        Rincr2P = (REOMAX - Rincr1P) * Rincr2

c        write(*,*) ' PEOZ=',PEOZ(IMode,1,IA1),Rincr2,Rincr2P

        PEOZ(IMode,1,IA1) =
     ;      PEOZ(IMode,1,IA1) + Rincr2P
        IF (nInt(ZCN) -(nInt(ZCN)/2)*2 .eq. 0) Then
           PEOZ(IMode,2,IA2) =
     ;         PEOZ(IMode,2,IA2) + Rincr2P
        Else
           PEOZ(IMode,2,IA2) =
     ;         PEOZ(IMode,2,IA2) - Rincr2P
        End if

        Rincr2N = (REOMAX - Rincr1N) * Rincr2
        PEON(IMode,1,IA1) =
     ;      PEON(IMode,1,IA1) + Rincr2N
        IF (nInt(RNCN)-(nInt(RNCN)/2)*2 .eq. 0) Then
           PEON(IMode,2,IA2) =
     ;         PEON(IMode,2,IA2) + Rincr2N
        Else
           PEON(IMode,2,IA2) =
     ;         PEON(IMode,2,IA2) - Rincr2N
        End if

c          /' Energy sorting '/
c     /' E1 = Etot * Gaussintegral(T2-T1,0.03); '/
        If( Abs(T1-T2) .lt. 1.E-6) Then
          E1 = 0.5E0 * Etot
        Else
          E1ES = Csort * T1 * T2 / ( Abs(T1 - T2) )
          E1ES = Min(E1ES,0.5E0*Etot)
c           /' Asymptotic value after "complete" energy sorting '/
          E1FG = Etot * A1 / ACN
c  /' in Fermi-gas regime '/
          If (Etot .lt. 13) Then
            E1 = E1ES
          end if
c ' complete energy sorting
          If (Etot .ge. 13 .and. Etot .le. 20) Then
c ' transition region
            E1 = E1ES + (Etot-13)/7*(E1FG-E1ES)
          End If
          If (Etot .gt. 20) Then
            E1 = E1FG
          end if
c  ' Fermi-gas regime
        End If
        E2 = Etot - E1
        EPART(IMode,1,IA1) = E1
c  /' Mean E* in light fragment '/
        EPART(IMode,2,IA2) = E2
c  /' Mean E* in heavy fragment '/

c        write(*,*) ' PEOZ=',PEOZ(IMode,1,IA1)

        end do
       end do


      EINTRSCISSION = Etot
c /' (For Mode 2) Global parameter '/

c   /'** RMS angular momentum of fission fragments **'/
c   /' Following Naik et al., EPJ A 31 (2007) 195 and  '/
c   /' S. G. Kadmensky, Phys. At. Nucl. 71 (2008) 1193 '/


c    /' CN spin '/
      ZT = ZCN
      AT = ACN
       If (Emode .eq. 2) Then
        AT = AT -1.
       end if
      IMAT = IMATENDF(ZT,AT)
      IF (IEiso .eq. 0) Then
c ' fissioning nucleus in ground state
       SpinCN =  RNucTab(IMAT,5)
      Else
c ' fissioning nucleus in isomeric state
c      IF NISOMAT(IMAT) < IEiso Then
c         Print "The isomer is not in the table of nuclear properties."
c         Print "Z, A, #iso ",ZT,AT,IEiso
c         Print "Please restart GEF."
c         End
c      End If
c       SpinCN = NucTab(IMAT + IEiso).RSPI
      End If

c      Write(*,*) 'ZT, AT, IMAT',ZT, AT, IMAT
c      Write(*,*) 'SPINCN',SpinCN
c '   Sleep


      Do IZ1 = 20 , nInt(ZCN) - 20
       AUCD = nInt(IZ1 * ACN / ZCN)
      Do IA1 = AUCD - 15 , AUCD + 15
c        /' Rigid momentum of inertia for spherical nucleus '/
        Irigidspher = 1.16E0**2 * IA1**1.6667E0 / 103.8415E0
c                /' unit: hbar^2/MeV '/
        Do IMode = 0 , 5

c          /' First (normally light) fission fragment: '/

          beta1 = Beta(IMode,1,IZ1)
          alph = beta1 / sqrt(4.E0 * pi / 5.E0)
          Irigid = Irigidspher * (1.E0 + 0.5E0*alph + 9.E0/7.E0*alph**2)
c                /' From Hasse & Myers, Geometrical Relationships ... '/
          Eexc = EPART(IMode,1,IA1)
          If (Eexc .lt. 0) Then
            Eexc = 0.
          end if
          T = UTemp(IZ1,IA1,Eexc)
          T = sqrt(T**2 + 0.8**2)
c ' For ZPM
          Ieff = Irigid * (1.E0 - 0.8E0 * exp(-0.693E0 * Eexc / 5.E0))
          Jrms = sqrt(2.E0 * Ieff * T)

c          /' Influence of CN spin '/
          Jrms = sqrt(Jrms**2 + 1./3. * SpinCN**2)

c          /' Incoming neutron (spin + orbital) '/
          If (Emode .eq. 2) Then
c            ' 2/3 * 1.16 * sqrt(2 * 939.65) / 197.33 = 0.1699
            Jrms = sqrt(Jrms**2 + 1./3 * 0.5**2 +
     ;         1./3. * (0.1699 * AT**0.333333 * sqrt(PEexc))**2)
          End If

          If (IZ1-(IZ1/2)*2 .gt. 0.5) Then
            Jrms = Jrms + Spinodd * (IA1/140.)**0.66667
c /' empirical '/
          End if

c           /' Additional angular momentum of unpaired proton. '/
c           /' See also Tomar et al., Pramana 68 (2007) 111 '/

          Jrms = Jrms * Jscaling
c      Write(*,*) IZ1,IMode,beta1,T,Eexc,SpinCN
c      Write(*,*) ' ',Irigidspher,Irigid,Ieff,Jrms

          SpinRMSNZ(IMode,1,IA1-IZ1,IZ1) = Jrms

c       write(*,*) 'spin',IMode,IZ1,IA1,Jrms


c      Write(*,*)  IA1,T,Eexc,Irigidspher,Irigid,Ieff,Jrms

c          /' Second (normally heavy) fission fragment: '/

          beta2 = Beta(IMode,2,IZ1)
          alph = beta2 / sqrt(4.E0 * pi / 5.E0)
          Irigid = Irigidspher * (1.E0 + 0.5E0*alph + 9.E0/7.E0*alph**2)
c                  /' From Hasse & Myers, Geometrical Relationships ... '/
          Eexc = EPART(IMode,2,IA1)
          If (Eexc .lt. 0) Then
           Eexc = 0
          end if

          T = UTemp(IZ1,IA1,Eexc)
          T = sqrt(T**2 + 0.8**2)
c      ' For ZPM
          Ieff = Irigid * (1.E0 - 0.8E0 * exp(-0.693E0 * Eexc / 5.E0))
          Jrms = sqrt(2.E0 * Ieff * T)

c          /' Influence of CN spin '/
          Jrms = sqrt(Jrms**2 + 1./3. * SpinCN**2)

c          /' Incoming neutron (spin + orbital) '/
          If(Emode .eq. 2) Then
c            ' 2/3 * 1.16 * sqrt(2 * 939.65) / 197.33 = 0.1699
            Jrms = sqrt(Jrms**2 + 1./3. * 0.5**2 +
     ;        1./3. * (0.1699 * AT**0.333333 * sqrt(PEexc))**2)
          End If

          If (IZ1 -(IZ1/2)*2 .gt. 0.5) Then
            Jrms = Jrms + Spinodd * (A1/140)**0.66667
          End if

c /' empirical '/
c           /' Additional angular momentum of unpaired proton. '/
c           /' See also Tomar et al., Pramana 68 (2007) 111 '/

          Jrms = Jrms * Jscaling

          SpinRMSNZ(IMode,2,IA1-IZ1,IZ1) = Jrms

        End do
        End do
        End do


       if(MODEC.eq.2) goto 1000

c       write(*,*) ' PAWidth2=',PAWidthS2

      do I=30,200
       WN0(i)=sq2pi*DEXP(-(I-ACMODE0)**2/(2.D0*SigAMode0**2))/SigAMode0
       WN1(i)=sq2pi*DEXP(-(I-ACMode1)**2/(2.D0*SigAMode1**2))/SigAMode1
       WN2(i)=sq2pi*DEXP(-(I-ACMode2)**2/(2.D0*SigAMode2**2))/SigAMode2
       WN3(i)=sq2pi*DEXP(-(I-ACMode3)**2/(2.D0*SigAMode3**2))/SigAMode3
        if(YieldMode11.ne.0.0) then
       WN11(i)=sq2pi*DEXP(-(I-ACMode0)**2/
     ;      (2.D0*SigAMode11**2))/SigAMode11
        else
       WN11(I)=0.0D0
        end if
        if(YieldMode22.ne.0.0) then
       WN22(i)=sq2pi*DEXP(-(I-ACMode0)**2/
     ;      (2.D0*SigAMode22**2))/SigAMode22
        else
       WN22(I)=0.0D0
        end if
      end do

       do i=30,200
        k1=i-30
         if(k1.gt.PAWidtS2) k1=PAWidtS2
        k2=200-i
         if(k2.gt.PAWidtS2) k2=PAWidtS2
         dk=k1+k2+1.
         i1=i-k1
         i2=i+k2
        do ik=i1,i2
         WN2P(i)=WN2P(i)+WN2(ik)
        end do
         WN2P(i)=WN2P(i)/dk
       end do

       do i=30,200
        WN2(i)=WN2P(i)
       end do

       MA0=ACN
c        write(*,*) ' MA0=',MA0

       do I=30,200
        WN(I)=WN0(I)*YieldMode0
     ;       +WN1(I)*YieldMode1
     ;       +WN2(I)*YieldMode2
     ;       +WN3(I)*YieldMode3
     ;       +WN11(I)*YieldMode11
     ;       +WN22(I)*YieldMode22

       end do

       sum=0.0
       do I=30,200
        WNT(I)=WN(I)+WN(MA0-I)
       sum=sum+WNT(I)
       end do

       do I=30,200
        WNT(I)=WNT(I)*100.
        write(*,*) I,WNT(I)
       end do

       write(*,*) ' Sum=', Sum*100.

        goto 3000

1000   continue

      if (iprint.ne.0) then
       write(*,*)
       Write(*,*) NEVT,' events will be calculated.'
      end if

       Racc = 1.E0 / NEVT


c   ' Print " "
c   ' Print "Pre-routines finished at ";time;"."

      RSum0 = YieldMode0
      RSum1 = RSum0 + YieldMode1
      RSum2 = RSum1 + YieldMode2
      RSum3 = RSum2 + YieldMode3
      RSum4 = RSum3 + YieldMode11
      RSum5 = RSum4 + YieldMode22

       Do ILoop = 1 , NEVT
c   /' Repeated loop '/

c      /' Chosing fission mode'/


      RChoice = Rndm(-1.)

c       write(*,*) RChoice,RSum0,RSum1,RSum2,RSum3

      IMode = 5
      IF (RChoice .lt. RSum0) Then
       IMode = 0
      end if
      If (RChoice .ge. RSum0 .And. RChoice .lt. RSum1) Then
       IMode = 1
      end if
      If (RChoice .ge. RSum1 .And. RChoice .lt. RSum2) Then
       IMode = 2
      end if
      If (RChoice .ge. RSum2 .And. RChoice .lt. RSum3) Then
       IMode = 3
      end if
      If (RChoice .ge. RSum3 .And. RChoice .lt. RSum4) Then
       IMode = 4
      end if
      If (RChoice .ge. RSum4 .And. RChoice .lt. RSum5) Then
       IMode = 5
      end if

c       write(*,*) ' IMode=',IMode

c ' IMode = 1  /'***********************************************'/

c    /' Chosing Z and A values '/

       If(IMode.eq. 0) then
        RAhelp = PGauss(ACMode0,SigAMode0)
c /' random choice of mass '/
        If (RAhelp .gt. 0.5 * ACN) Then
          RAheavy = RAhelp
        Else
          RAheavy = ACN - RAhelp
        End If
        RZpol = Zshift(0,2,nInt(RAheavy))
c /' local polarization '/
        RZ = RAheavy * ZCN / ACN + RZpol
c /' mean position in Z for given mass '/
        RZheavy = PGauss(RZ,SigPolMode0)
c /' random choice of Z '/
       end if

       If(IMode.eq. 1) then
c      Case 1
        RAheavy = PGauss(ACMode1,SigAMode1)
        RZpol = Zshift(1,2,nInt(RAheavy))
        RZ = RAheavy * ZCN / ACN + RZpol
        RZheavy = PGauss(RZ,SigPolMode1)
       end if

       If(IMode.eq. 2) then
c      Case 2
900     continue
        RAheavy = PGauss(ACMode2,SigAMode2)
        RAheavy = PBox(ACMode2,SigAMOde2,PAWidthS2)
        RAheavy = max(RAheavy,1.)
        RZpol = Zshift(2,2,nInt(RAheavy))
        RZ = RAheavy * ZCN / ACN + RZpol
        Rtest = RNDM(-1.)
        If (Rtest.gt.0.5*ERF((RZ-ZTRUNC50)/(FTRUNC50*SigZMode2))+0.5E0
     ;  .Or.
     ;   Rtest .gt.0.5*ERF((ZCN-RZ-ZTRUNC28)/(FTRUNC28*SigZMode2))+0.5)
     ;  Then
         Goto 900
        else
c        /' truncation below Z = 35 and below Z = 50 due to properties of deformed shells '/
         RZheavy = PGauss(RZ,SigPolMode2)
         end if
        end if

       If(IMode.eq. 3) then
c      Case 3
        RAheavy = PGauss(ACMode3,SigAMode3)
        RZpol = Zshift(3,2,nInt(RAheavy))
        RZ = RAheavy * ZCN / ACN + RZpol
        RZheavy = PGauss(RZ,SigPolMode3)
       end if

       If(IMode.eq. 4) then
c      Case 4
        RAhelp = PGauss(ACMode0,SigAMode11)
        If (RAhelp .gt. 0.5 * ACN) Then
          RAheavy = RAhelp
        Else
          RAheavy = ACN - RAhelp
        End If
        RZpol = 0
        RZ = RAheavy * ZCN / ACN
        RZheavy = PGAUSS(RZ,SigPolMode0)
       end if

       If(IMode.eq. 5) then
c      Case 5
        RAhelp = PGauss(ACMode0,SigAMode22)
        If (RAhelp .gt. 0.5 * ACN) Then
          RAheavy = RAhelp
        Else
          RAheavy = ACN - RAhelp
        End If
        RZpol = 0
        RZ = RAheavy * ZCN / ACN
        RZheavy = PGauss(RZ,SigPolMode0)
       end if
       RZheavy=max(RZheavy,1.)
       RZheavy=min(RZheavy,ZCN-1.)
       RAheavy=max(RAheavy,1.)
       RAheavy=min(RAheavy,ACN-1.)

       RZlight = ZCN - RZheavy
       RAlight = ACN - RAheavy

c        write(*,*) ' RAlight,RAheavy=',RAlight,RAheavy
c        write(*,*) ' RZlight,RZheavy=',RZlight,RZheavy

      RNheavy = RAheavy - RZheavy
       nRA=min(nInt(RAheavy),350)
       INheavy = Ivenodd(RNheavy,PEON(IMode,2,nRA))
       IZheavy = Ivenodd(RZheavy,PEOZ(Imode,2,nRA))

c        write(*,*) ' RZheavy ',RZheavy, PEOZ(Imode,2,nInt(RAheavy))

c       INheavy=RNheavy
c       IZheavy=RZheavy

      IAheavy = INheavy + IZheavy

      INlight = RNCN - INheavy
      IZlight = ZCN - IZheavy
      IAlight = ACN - IAheavy

      INlight=max(INlight,1)
      IZlight=max(IZlight,1)
      IAlight=max(IAlight,1)

      INheavy=min(INheavy,200)
      IAheavy=min(IAheavy,350)
      IZheavy=min(IZheavy,150)

c    /' Correct mass distribution, pre-neutron '/

      APRE(IAheavy) = APRE(IAheavy) + Racc
      APRE(IAlight) = APRE(IAlight) + Racc
      AMPRE(IMode,IAheavy) = AMPRE(IMode,IAheavy) + Racc
      AMPRE(IMode,IAlight) = AMPRE(IMode,IAlight) + Racc

c    /' Nuclide distribution, pre-neutron '/

      ANZPRE(INheavy,IZheavy) = ANZPRE(INheavy,IZheavy) + 1.
      ANZPRE(INlight,IZlight) = ANZPRE(INlight,IZlight) + 1.
      ANPRE(INheavy) = ANPRE(INheavy) + Racc
      ANPRE(INlight) = ANPRE(INlight) + Racc
      ANMPRE(IMode,INheavy) = ANMPRE(IMode,INheavy) + Racc
      ANMPRE(IMode,INlight) = ANMPRE(IMode,INlight) + Racc

      ZISOPRE(IAHeavy,IZHeavy) = ZISOPRE(IAHeavy,IZHeavy) + Racc
      ZISOPRE(IALight,IZLight) = ZISOPRE(IALight,IZLight) + Racc

c    27.09.2013
c     filling of the array of P(A,Z) distribition
      yAZ_temp(IAheavy,IZheavy)=yAZ_temp(IAheavy,IZheavy)+1
      yAZ_temp(IAlight,IZlight)=yAZ_temp(IAlight,IZlight)+1

c    /' Excitation energy of fragments '/

      If (IMode .le. 3) Then
       Eexcheavymean = Edefo(IMode,2,IZheavy)
c /' Only deformation energy '/
       Eexclightmean = Edefo(IMode,1,IZlight)
c /' Only deformation energy '/
       ESIGDEFOheavy =
     ;  ( LyMass(Float(IZheavy),Float(IAheavy),beta(IMode,2,IZheavy)
     ; + SIGDEFO)-
     ;    LyMass(Float(IZheavy),Float(IAheavy),beta(IMode,2,IZheavy)))
       ESIGDEFOlight =
     ;  ( LyMass(Float(IZlight),Float(IAlight),beta(IMode,1,IZlight)
     ; + SIGDEFO)-
     ;   LyMass(Float(IZlight),Float(IAlight),beta(IMode,1,IZlight) ))

c ' If beta(IMOde,1,IZlight) > 0.68 Then
c '   ESIGDEFOlight = 0.4*ESIGDEFOlight
c ' End If

       End If

       If (IMode .eq. 4) Then
        Eexcheavymean = Edefo(1,2,IZheavy)
        Eexclightmean = Edefo(1,2,IZlight)
c /' Shell effect stored for "heavy" fragment '/
      ESIGDEFOheavy =
     ;   ( LyMass(Float(IZheavy),Float(IAheavy),beta(1,2,IZheavy)
     ; + SIGDEFO) -
     ;   LyMass(Float(IZheavy),Float(IAheavy),beta(1,2,IZheavy) ))
      ESIGDEFOlight =
     ;   ( LyMass(Float(IZlight),Float(IAlight),beta(2,2,IZlight)
     ; + SIGDEFO) -
     ;   LyMass(Float(IZlight),Float(IAlight),beta(2,2,IZlight) ))
       End If
       If (IMode .eq. 5) Then
      Eexcheavymean = Edefo(2,2,IZheavy)
      Eexclightmean = Edefo(2,2,IZlight)
c /' Shell effect stored for "heavy" fragment '/

      ESIGDEFOheavy =
     ;  ( LyMass(Float(IZheavy),Float(IAheavy),beta(2,2,IZheavy)
     ; + SIGDEFO) -
     ;   LyMass(Float(IZheavy),Float(IAheavy),beta(2,2,IZheavy) ))
      ESIGDEFOlight =
     ; ( LyMass(Float(IZlight),Float(IAlight),beta(2,2,IZlight)
     ; + SIGDEFO) -
     ;   LyMass(Float(IZlight),Float(IAlight),beta(2,2,IZlight) ))
      End If
      If (Eexcheavymean .lt. 0) Then
        Eexcheavymean = 0
      End if
      If (Eexclightmean .lt. 0) Then
        Eexclightmean = 0
      End if

1010  continue
      Eexcheavy = PGAUSS(Eexcheavymean,ESIGDEFOheavy)
      if (Eexcheavy .lt. 0) goto 1010
1020  continue
      Eexclight = PGAUSS(Eexclightmean,ESIGDEFOlight)
      if (Eexclight .lt. 0) goto 1020

c    /' Assumption: width in TKE is the '/
c    /' quadratic sum of widths in defo and coll. '/
c    /' Remark of caution: The width of the TKE for fixed mass contains often
c       several fission modes. The width in Lang et al. for fixed Z contains several A,
c       which contributes already with about 3% to the width. Therefore, the
c       width in TXE (or TKE) for fixed A and fixed mode may be much smaller than
c       4 or 5 percent! '/

c        write(*,*) ' A,Z=',IAheavy,IAlight
c        write(*,*) ' Ubound=', 800,10*Eexcheavy+0.5,
c    ;  10*Eexclight+0.5

       If (10*Eexcheavy+0.5 .le. 800) Then
        Edefo2d(IAheavy,nInt(10*Eexcheavy+0.5))=
     ;    Edefo2d(IAheavy,nInt(10*Eexcheavy+0.5))+Racc
       End If
       If (10*Eexclight+0.5 .le. 800) Then
         Edefo2d(IAlight,nInt(10*Eexclight+0.5))=
     ;    Edefo2d(IAlight,nInt(10*Eexclight+0.5))+Racc
       End If

c    /' Temperatures of fragments '/

      If (IMode .le. 3) Then
       Tlight = Temp(IMode,1,IAlight)
       Theavy = Temp(IMode,2,IAheavy)
      End If
      If (IMode .eq. 4) Then
       Tlight = Temp(1,2,IAlight)
       Theavy = Temp(1,2,IAheavy)
      End If
      If (IMode .eq. 5) Then
       Tlight = Temp(2,2,IAlight)
       Theavy = Temp(2,2,IAheavy)
      End If


c    /' Intrinsic excitation energies of fragments '/

      Eintrlightmean = EPART(IMode,1,IAlight)
      If (Eintrlightmean .lt. 0) Then
       Eintrlightmean = 0.
      End if
      Eintrheavymean = EPART(IMode,2,IAheavy)
      If (Eintrheavymean .lt.0.) Then
        Eintrheavymean = 0.
      End if


1030   continue
       Eintrlight=PGauss(Eintrlightmean,EexcSIGrel*Eintrlightmean+0.5)
       if (Eintrlight .lt. 0) goto 1030
1040   continue
       Eintrheavy=PGauss(Eintrheavymean,EexcSIGrel*Eintrheavymean+0.5)
       if (Eintrheavy .lt. 0) goto 1040

      Eintrlight = Eintrlight - Aypair(Float(IZlight),Float(IAlight))
      Eintrheavy = Eintrheavy - Aypair(Float(IZheavy),Float(IAheavy))
c          /' Staggering of BE by pairing '/
c          /' Assumption: pairing only felt in the lowest nuclear levels, '/
c          /' at the end of the evaporation cascade '/
c          /' (This should be a good assumption for higher excitation energies. '/
c          /' Some deviations occur due to the even-odd effect at low exc. energies. '/

      If (Eintrheavy .lt. 100.) Then
       Eintr2d(IAheavy,nInt(10*Eintrheavy+0.5))=
     ;    Eintr2d(IAheavy,nInt(10*Eintrheavy+0.5))+Racc
      End If
      If (Eintrlight .lt. 100.) Then
       Eintr2d(IAlight,nInt(10*Eintrlight+0.5))=
     ;    Eintr2d(IAlight,nInt(10*Eintrlight+0.5))+Racc
      End If

      Eexcheavy = Eexcheavy + Eintrheavy
      Eexclight = Eexclight + Eintrlight
c    /' Now: deformation + intrinsic excitation energy '/

      Ecollmean=ECOLLFRAC*(DeSaddleScission(ZCN**2/ACN**0.33333E0)
     ;  -Etunn)
      If (Ecollmean .lt. 0) Then
       Ecollmean = 0
      End if

c    /' Experimental data of prompt neutron yields show an enhancement of the '/
c    /' neutron yield for odd-Z CN, corresponding to an enhanced E* by about 1.6 MeV '/
c    /' The enhancement is equally shared to the light and the heavy fragment. '/
c    /' The neutron number of the CN has no influence. '/
c    /' The origin of this effect is not clear. '/
c    /' By technical reasons, this additional energy is introduced here into the '/
c    /' collective energy at scission, because this energy is divided equally between both fragments. '/
c    /'  KHS, 31. Jan. 2012 '/
      If ( (nInt(ZCN)-(nInt(ZCN)/2)*2) .eq. 1) Then
        Ecollmean = Ecollmean + 1.6
      end if

       Ecoll = -1.
1050    continue
       Ecoll = PGauss(Ecollmean,EexcSIGrel* Ecollmean)
       if (Ecoll .lt. 0.) goto 1050
      Ecoll = Ecoll + Ecollsaddle(IMode)
      Ecollheavy = 0.5E0 * Ecoll
      Ecolllight = 0.5E0 * Ecoll


      Ecoll2d(IAheavy,nInt(10*Ecollheavy + 0.5)) =
     ;   Ecoll2d(IAheavy,nInt(10*Ecollheavy + 0.5)) + Racc
      Ecoll2d(IAlight,nInt(10*Ecolllight + 0.5)) =
     ;   Ecoll2d(IAlight,nInt(10*Ecolllight + 0.5)) + Racc

      Eexcheavy = Eexcheavy + Ecollheavy
      Eexclight = Eexclight + Ecolllight
c    /' Now: also collective excitation energy added '/

      If (Eexcheavy .lt. 100) Then
       Eexc2d(IAheavy,nInt(10*Eexcheavy+0.5))=
     ;   Eexc2d(IAheavy,nInt(Eexcheavy+0.5))+Racc
      End If
      If (Eexclight .lt. 100) Then
       Eexc2d(IAlight,nInt(10*Eexclight+0.5))=
     ;   Eexc2d(IAlight,nInt(Eexclight+0.5))+Racc
      End If


c    /'** Angular momentum of fragments **'/



       JFraglight=0.
      If (IMode .le. 3) Then

c        write(*,*) ' Spin=',SpinRMSNZ(IMode,1,INLight,IZlight)

       JFraglight=PLinGauss(SpinRMSNZ(IMode,1,INLight,IZlight)/
     ; sqrt(2.))-0.5
       JFragheavy=PLinGauss(SpinRMSNZ(IMode,2,INHeavy,IZheavy)/
     ; sqrt(2.)) - 0.5
      End If
      If (IMode .eq. 4) Then
       JFraglight=PLinGauss(SpinRMSNZ(1,2,INLight,IZlight)/
     ; sqrt(2.)) - 0.5
       JFragheavy = PLinGauss(SpinRMSNZ(1,2,INHeavy,IZheavy)/
     ; sqrt(2.)) - 0.5
      End If
      If (IMode .eq. 5) Then
       JFraglight = PLinGauss(SpinRMSNZ(2,2,INLight,IZlight)/
     ; sqrt(2.)) - 0.5
       JFragheavy = PLinGauss(SpinRMSNZ(2,2,INHeavy,IZheavy)/
     ; sqrt(2.)) - 0.5
      End If

      If (JFraglight .lt. 0) Then
       JFraglight = 0.
      Endif
      JFRAGpre(INlight,IZlight,nInt(JFraglight)) =
     ;  JFRAGpre(INlight,IZlight,nInt(JFraglight)) + 1.

      If (JFragheavy .lt. 0) Then
       JFragheavy = 0
      Endif

      JFRAGpre(INheavy,IZheavy,nInt(JFragheavy)) =
     ; JFRAGpre(INheavy,IZheavy,nInt(JFragheavy)) + 1.


      Irigidspher = 1.16E0**2 * IAlight**1.6667E0 / 103.8415E0

c       write(*,*) ' Irigsph=',Irigidspher,JFraglight

      IfragEff = 0.45 * Irigidspher
      Erotlight =  JFraglight*(JFraglight+1)/(2*IfragEff)

      Irigidspher = 1.16E0**2 * IAheavy**1.6667E0 / 103.8415E0
      IfragEff = 0.45 * Irigidspher
      Erotheavy =  JFragheavy*(JFragheavy+1)/(2*IfragEff)

c  /' Kinetic energies of fragments '/

c      ' TXE includes E*CN
      Qvalue = UMass(ZCN,ACN) + Aypair(ZCN,ACN) -
     ;     (UMass(Float(IZheavy),Float(IAheavy))+
     ; Aypair(Float(IZheavy),Float(IAheavy))
     ;   + UMass(Float(IZlight),Float(IAlight))+
     ; Aypair(Float(IZlight),Float(IAlight)))
      TXE = Eexclight + Eexcheavy + Erotlight + Erotheavy

c       write(*,*) ' Ex=',Eexclight,Eexcheavy,Erotlight,Erotheavy
c       write(*,*) ' Qvalue,TXE=',Qvalue,TXE

      Etotal = Qvalue + REexcGS
      TKE = Etotal - TXE
      If (TKE .lt. 0) Then

        Write(*,*) ' <E> Event with excessive excitation energy found.'
        TXEcorr = TXE + TKE - 1
        TKE = 1.
        Eexclight = Eexclight * TXEcorr/TXE
        Eexcheavy = Eexcheavy * TXEcorr/TXE
        TXE = TXEcorr
      End If

      Ekinlight = TKE * IAheavy / ACN
      Ekinheavy = TKE * IAlight / ACN

      If(IAheavy.lt.200.and.IAheavy.gt.1) then
      If(Ekinheavy+1.lt.300.and.
     ; Ekinheavy.gt.1) Then
        AEkin(IAheavy,nInt(Ekinheavy+0.5))=
     ;  AEkin(IAheavy,nInt(Ekinheavy+0.5)) + 1.
      End If
      End if

      If(IAlight.lt.200.and.IAlight.gt.1) then
      If(Ekinlight+1.lt.300.and.
     ; Ekinlight.gt.1) Then
        AEkin(IAlight,nInt(Ekinlight+0.5)) =
     ; AEkin(IAlight,nInt(Ekinlight+0.5)) + 1.
      End If
      End if

      If (Bsub) Then

c   /' Test output '/
c        Write(*,*) ' IZ_light,I_A_light,I_Z_heavy,I_A_heavy,TKE'
c        Write(*,*) IZlight,' ',IAlight,' ',IZheavy,' ',IAheavy,' ',TKE

      End If

c ' Print " "
c    /'** Neutron evaporation **'/

      Ngtot = 0
      Nglight = 0
      Ngheavy = 0
      Egtot10 = 0

c     ' Pre-scission kinetic energy
      Epre = Epotscission - Ecoll - Eintrlight - Eintrheavy -
     +  Erotlight - Erotheavy
      If(Epre.lt.0) Then
        Epre = 0
      end if
c       write(*,*) ' Epre=',Epre,' Imode=',Imode
c       write(*,*) ' IAheavy, IAlight=',IAheavy, IAlight

      If(IMode .le. 3) Then
       TheavyFF = TempFF(IMode,2,IAheavy)
       TlightFF = TempFF(IMode,1,IAlight)
      End If
      If(IMode .eq. 4) Then
c   ' S11
       TheavyFF = TempFF(1,2,IAheavy)
       TlightFF = TempFF(1,2,IAlight)
      End If
      If(IMode .eq. 5) Then
c  ' S22
       TheavyFF = TempFF(2,2,IAheavy)
       TlightFF = TempFF(2,2,IAlight)
      End If

c       write(*,*) ' TheavyFF, TlightFF=',TheavyFF, TlightFF

       do I = 1 , 100
        ArrayEn(I) = 0
        ArrayTn(I) = 0
       end do

      call Eva(2,IZheavy,IAheavy,Eexcheavy,TheavyFF,RZHeavyPost,
     ; RAHeavyPost,EFinalHeavy,EKIN)

c      write(*,*) ' E_Final_Heavy',EFinalHeavy,' EKIN=',EKIN
c      write(*,*) ' ZheavePost, Aheavypost=',RZHeavyPost,RAHeavyPost

      RNHeavyPost = RAHeavyPost - RZHeavyPost
      IZHeavyPost = nInt(RZHeavyPost)
      INHeavyPost = nInt(RNHeavyPost)
      IAHeavyPost = IZHeavyPost + INHeavyPost

      INheavyPost=min(INheavyPost,200)
      IAheavyPost=min(IAheavyPost,350)
      IZheavyPost=min(IZheavyPost,150)

      APOST(IAheavyPost)=APOST(IAheavyPost)+Racc
      AMPOST(IMODE,IAheavyPost)=AMPOST(IMODE,IAheavyPost)+Racc

      ZPOST(IZheavyPost)= ZPOST(IZheavyPost)+Racc
      ZMPOST(IMode,IZheavypost)=ZMPOST(IMode,IZheavypost)+Racc

      ANPOST(INHeavyPost)=ANPOST(INHeavyPost)+Racc
      ANMPOST(IMode,INHeavyPost)=ANMPOST(IMode,INHeavyPost)+Racc

      ANZPOST(INheavyPost,IZheavypost)=
     ;     ANZPOST(INheavyPost,IZheavypost)+Racc
      iRN = IAheavy - IAHeavyPost
      If (iRN .gt. 0 .and. iRN .le. 20) Then
       AN2dpost(IAHeavyPost,iRN)=AN2dpost(IAHeavyPost,iRN)+1.
       AN2dpre(IAHeavy,iRN)=AN2dpre(IAHeavy,iRN)+1.
      End If

       do I = 1 , 100
        ArrayEn(I) = 0
        ArrayTn(I) = 0
       end do

      call Eva(1,IZlight,IAlight,Eexclight,TlightFF,RZlightPost,
     ; RAlightPost,EFinallight,EKIN)

c      write(*,*) ' E_Final_light',EFinallight,' EKIN=',EKIN
c      write(*,*) ' ZlightPost, Alightpost=',RZlightPost,RAlightPost

      RNlightPost = RAlightPost - RZlightPost
      IZlightPost = nInt(RZlightPost)
      INlightPost = nInt(RNlightPost)
      IAlightPost = IZlightPost + INlightPost

      IZlightPost=max(IZlightPost,1)
      INlightPost=max(INlightPost,1)
      IAlightPost=max(IAlightPost,1)

      APOST(IAlightPost)=APOST(IAlightPost)+Racc
      AMPOST(IMODE,IAlightPost)=AMPOST(IMODE,IAlightPost)+Racc

      ZPOST(IZlightPost)= ZPOST(IZlightPost)+Racc
      ZMPOST(IMode,IZlightpost)=ZMPOST(IMode,IZlightpost)+Racc

      ANPOST(INlightPost)=ANPOST(INlightPost)+Racc
      ANMPOST(IMode,INlightPost)=ANMPOST(IMode,INlightPost)+Racc

      ANZPOST(INlightPost,IZlightpost)=
     ;     ANZPOST(INlightPost,IZlightpost)+Racc

      iRNtot = iRN + IALight - IALightPost
      If (iRNtot .gt. 0 .and. iRNtot .le. 50) Then
      aNN(iRNtot)=aNN(iRNtot)+1
      End If


      iRN = IAlight - IAlightPost
      If (iRN .gt. 0 .and. iRN .le. 20) Then
       AN2dpost(IAlightPost,iRN)=AN2dpost(IAlightPost,iRN)+1.
       AN2dpre(IAlight,iRN)=AN2dpre(IAlight,iRN)+1.
      End If

c    26.10.2013
c     filling of the array of Post(A,Z) distribition
      yAZp_temp(IAheavyPost,IZheavyPost)=
     ;  yAZp_temp(IAheavyPost,IZheavyPost)+1
      yAZp_temp(IAlightPost,IZlightPost)=
     ;  yAZp_temp(IAlightPost,IZlightPost)+1


      end do

      yield_norm=0
      yieldp_norm=0
      if (iprint.ne.0) then
        write(*,*)
        write(*,*) '         N   Yield(pre,post-neutron)'
        write(*,*)
       do i=20,120
        write(*,9000) i,anpre(i),anpost(i)
       end do

        write(*,*)
        write(*,*) '         A   Yield(pre,post-neutron)'
        write(*,*)
       do i=50,180
        write(*,9000) i,apre(i),apost(i)
       end do
       endif

       do i=1,350
        yield_norm=yield_norm+apre(i)
        yieldp_norm=yieldp_norm+apost(i)
       end do

      if (iprint.ne.0) then
       write(*,*)'Yield normalized to ', yield_norm,yieldp_norm
      end if

       yield_norm=0
       yieldp_norm=0
       do i=1,200
         do j=1,150
            yield_norm=yield_norm+yAZ_temp(i,j)
            yieldp_norm=yieldp_norm+yAZp_temp(i,j)
         end do
       end do

c       do i=1,350
c          apre(i)=apre(i)*2.0/yield_norm
c       end do

       do i=1,200
         do j=1,150
            yAZ_temp(i,j)=yAZ_temp(i,j)*2.0/yield_norm
            yAZp_temp(i,j)=yAZp_temp(i,j)*2.0/yield_norm
         end do
       end do


C   /' Calculation of mean neutron yield per mass '/
       Do I = 1 , 350
       Zaehler = 0
       aNenner = 0
       Do J = 0 , 20
         Zaehler = Zaehler + J * aN2Dpre(I,J)
         aNenner = aNenner + aN2Dpre(I,J)
       end do
        if(aNenner.ne.0.0) then
       ANApre(I) = Zaehler / aNenner
        else
       ANApre(I)=0.0
        end if
       end do
       Do I = 1 , 350
       Zaehler = 0
       aNenner = 0
       Do J = 0 , 20
         Zaehler = Zaehler + J * aN2Dpost(I,J)
         aNenner = aNenner + aN2Dpost(I,J)
       end do
        if(aNenner.ne.0.0) then
       ANApost(I) = Zaehler / aNenner
        else
       ANApost(I)=0.0
        end if
       end do

      if (iprint.ne.0) then
       write(*,*)
       write(*,*) '  Neutron evaporation multiplicity'
       write(*,*)

       do i=50,180
        write(*,9000) i,anapre(i),anapost(i)
       end do

      end if

       aNsum = 0
       aNmean = 0
      do I = 0 , 50
       aNsum = aNsum + aNN(I)
      end do
      do I = 0 , 50
      if(aNsum.gt.0.) then
       aNN(I) = aNN(I) / aNsum
      else
       aNN(I)=0.0
      end if
      end do
      do I = 0 , 50
       aNmean = aNmean + I * aNN(I)
      end do

      if (iprint.ne.0) then
       write(*,*)
       write(*,*) '--- Neutron-multiplicity distribution ---'
       write(*,*)
      do I = 0 , 50
      If (aNN(I) .gt. 0) Then
        write(*,*) I, aNN(I)
      End If
      end do
       write(*,*)
      write(*,*) 'Mean value  ',aNmean

      end if


c      27.09.2013
c      adding data into the array of the total P(A,Z) yield
       do i=1,200
         do j=1,150
            yAZ(i,j)=yAZ(i,j)+yAZ_temp(i,j)*fsig(nex)
            yAZp(i,j)=yAZp(i,j)+yAZp_temp(i,j)*fsig(nex)
         end do
       end do

c      --- here we add data into the array with the total yield
       do i=1,350
c    write(*,*)'fsig(',nex,')=',fsig(nex)
         ysum(i)=ysum(i)+apre(i)*fsig(nex)
         ysump(i)=ysump(i)+apost(i)*fsig(nex)

         anpre_sum(i)=anpre_sum(i)+anapre(i)*fsig(nex)
         anpost_sum(i)=anpost_sum(i)+anapost(i)*fsig(nex)
       end do

       do i=1,50
        ann_sum(i)=ann_sum(i)+ann(i)*fsig(nex)
       end do
       anMean_sum=anMean_sum+anMean*fsig(nex)

3000   continue

2000   continue

9000   format(i5,F10.5,5x,F10.5)

         return
         end


           Function AyPair(Z,A)
           Epair=-12.0/sqrt(A)
           iz=nInt(Z)
           ia=nInt(A)
           iz2=iz/2
           ia2=ia/2
           if(iz.ne.iz2*2) then
            mz=0
           else
            mz=1
           endif
           if(ia.ne.ia2*2) then
            ma=0
           else
            ma=1
           endif
            m=mz+ma
            AyPair=Epair*m
           return
           end

      Subroutine BetaEquiA(A1,A2,Z1,Z2,d,beta1prev,beta2prev,
     ;    beta1opt,beta2opt)
c    /' Analytical approximation to the numerical calculations in BetaEqui '/

       data x/0./,x0/1.26/, Zrel/1./
       data ZsqrAmem/1./
       data P0/0.1339/,P1/1.452/,P2/-3.09276/,P3/4.26162/,P4/-2.545356/

       ZsqrAin = (Z1+Z2)**2 / (A1+A2)**0.333333
       If (ZsqrAin .ne. ZsqrAmem) Then
         x0 = 1.26
         x = ZsqrAin * 0.001 - x0
         P0 = 0.1339
         If (x .lt. x0) Then
           P1 = 1.452
c '          P2 = -2.13 * P1
c '          P3 = 2.935 * P1
c '          P4 = -1.753 * P1
         Else
           P1 = 1.452 + 1.8387 * x**2
c '          P2 = -(3.13325 - 0.6697 * x + 5.375 * x^2)
c '          P3 = 4.305  - 1.3566 * x + 9.839 * x^2
c '          P4 = -(2.503 + 0.318 * x + 3.027 * x^2)
         End If
         P2 = -2.13 * P1
         P3 = 2.935 * P1
         P4 = -1.753 * P1
       End If

       Zrel = Z1 / (Z1 + Z2)

       beta1opt = P0 + P1 * Zrel + P2 * Zrel**2
     ;          + P3 * Zrel**3 + P4 * Zrel**4
       beta2opt = beta1opt
c  ' This is wrong, but beta2opt is not used after call of BetaEqui (!!!)

       ZsqrAmem = ZsqrAin

      return
      End

       Subroutine Eva(Ilh,ZCN,ACN,EINIT,T,ZRES,ARES,EFINAL,EKIN)
c       real*4 Ni,Nf
        integer ZCN,ACN
         common /cgamma/Egamma(1000),ArrayEn(100),ArrayTn(100)

c            /' Z_CN,A_CN,E_init       Parameters of initial nucleus '/
c            /' Z_RES,A_RES,E_FINAL    after evaporation '/
c            /' T temperature coefficient of level density '/
c            /' E_kin kinetic energy of neutron '/
       data EMIN/0.0/,Tscale/0.85/
c   ' Final energy for evaporation chain '/
c            Dim As Single SN,SNeff       /' Neutron separation energy '/

          Tnacc = 0.0
      Ngtot = 0
      Nglight = 0
      Ngheavy = 0
      Egtot10 = 0

          Ifold = 1

c           write(*,*) ' Z,A,E=',ZCN,ACN,EINIT

          Ai = ACN
          Zi = ZCN
          Ei = EINIT

          SN = UMass(Zi,Ai-1.E0) + Aypair(Zi,Ai-1.E0)
     ;       - (UMass(Zi,Ai) + Aypair(Zi,Ai))
c  '                    /' Lypair: even-odd staggering of masses '/
          SNeff = SN

c           write(*,*) ' SN=',SN

          Zf = Zi
          Af = Ai
          Ef = Ei

50      continue

c        /' Treat gamma competition '/
            Tm = Utemp(nInt(Zi),nInt(Ai),Ei)
c ' emitting nucleus
            Td = Utemp(nInt(Zi),nInt(Ai),Ei-SNeff)

            Gammag = 0.624 * Ai**1.6 / ( Tm**5 )
c  ' in meV (Ignatyuk, Bologna)
            Gammag = Gammag * 1.E-9
c ' in MeV

            Tmean = (Tm + Td)/2

c '          Gamma_n = (Ai-1)^0.66667 * 0.0302 * Td^2 / exp(SN/Tmean)   ' in MeV (Moretto)
            Gamman = (Ai-1.)**0.66667 * 0.13 * Td**2 / exp(SN/Tmean)
c ' in MeV (Mor. PRC 54,3062)

            Tn = 0.658 / Gamman
c ' in units of 10^-21 s (hbar=0.658zs*MeV)
c '  Tn = Tn + 100
            Tnacc = Tnacc + Tn

            If( (Ei-Sneff) .lt. Abs(Aypair(Zi,Ai-1))) Then
c ' rest energy below mean g.s. of odd-odd nuclei
              If((nInt(Zi)-(nInt(Zi)/2)*2).eq.0 .or.
     ;           (nInt(Ai-Zi-1)-(nInt(Ai-Zi-1)/2)*2).eq.0) Then
c ' even Z or even N
                Gamman = 0.1 * Gamman
              End If
              If((nInt(Zi)-(nInt(Zi)/2)*2).eq.0 .and.
     ;           (nInt(Ai-Zi-1)-(nInt(Ai-Zi-1)/2)*2).eq.0) Then
c              If Zi Mod 2 < 0.5 and (Ai-Zi-1) Mod 2 < 0.5 Then
c ' even Z and even N
c              ' For low level density below pairing gap in even-even nuclei
                Gamman = 0.1 * Gamman
              End If
            End If
c              /' Reduces the even-odd effect in neutron number
c                 due to low level density below the pairing gap '/

            Pgamma = Gammag / (Gammag + Gamman)

            Pgamma = 0.3 * Pgamma
c ' correction to Ignatyuk-Moretto estimation

            If (RNDM(-1.) .lt. Pgamma) Then
c ' gamma will be emitted

                Eg = PEgammahigh(Zi,Ai,Ei)

c                ' Accumulate E1 gammas
                N = nInt(Eg*10)
                If (N .gt. 0) Then
                  Ngtot = Ngtot + 1
                  If (Ilh .eq. 1) Then
                   Nglight = Nglight + 1
                  endif
                  If (Ilh .eq. 2) Then
                    Ngheavy = Ngheavy + 1
                  endif
                  Egtot10 = Egtot10 + EG*10.
                  If (N .le. 1000) Then
                    Egamma(N) = Egamma(N) + 1.
                  End If
                End If
                Ei = Ei - Eg
            End If

            IF (Ei-Sneff .le. EMIN) THEN
              Zf = Zi
              Af = Ai
              Ef = Ei
              Goto 200
            End If

            ITry = 0
100          continue
            ITry = ITry + 1
            If (ITry .lt. 100) Then
                Td = UTemp(nInt(Zi),nInt(Ai-1),Ei-SNeff)
                TCT = TEgidy(Ai-1,0.,Tscale)
c          '     E_kin = PMaxwellv(Td)  /' Maxwell with 1/v '/
c          '     E_kin = PEN0(Td,Ai-1)    /' Maxwell with opt. model Koning et al. '/
                Ekin = PEN1(Td,TCT,Ei-SNeff,Ai-1)
c                           /' modified Maxwell with opt. model Koning et al. '/
c          '     E_kin = PEN2(Td,Ai-1)    /' Maxwell with inv. x sect. Blatt & W. '/
c          '     E_kin = PENint(Td)     /' Maxwell with inv. x sect. Dostrovsky '/
c          '     E_kin = PMaxwell(Td)   /' Maxwell, const. x section '/
            Else
c              /' E_kin too high after several attemps '/
c              /' no neutron emitted '/
              Af = Ai
              Zf = Zi
              Ef = Ei

              GOTO 200
            EndIf

            If (Ekin .gt. Ei-SNeff) Then
c              /' E_kin from PMaxwell is not available '/
c              /' Try again '/
              Goto 100
            End If

            Af = Ai - 1
            Zf = Zi
            Ef = Ei - Ekin - SNeff

c           /' ANAL(EN,E_kin);
c              ANAL(ENM(I_MODE),E_kin); '/

            SN = (UMass(Zf,Af-1.E0) + Aypair(Zf,Af-1.E0))
     ;             - (UMass(Zf,Af) + Aypair(Zf,Af))

            SNeff = SN
            Ai = Af
            Zi = Zf
            Ei = Ef

            ArrayEn(Ifold) = Ekin
            ArrayTn(Ifold) = Tnacc
            Ifold = Ifold + 1


          if (Ei-SNeff .gt. EMIN) then
            goto 50
          end if

200    continue
          ARES = Af
          ZRES = Zf
          EFINAL = MAX(Ef,0.)
        return
       End


            Function ALDMas(Z,A,beta)
      common /datain/beldm(203,136),ushel(203,136),RNucTab(3885,8)
         real*4 Lymass
            AN = A - Z
            AN=max(AN,1.)
            BEtab = BELDM(nInt(AN),nInt(Z)) + 2 * 12 / sqrt(A)
     ;                  - 0.00001433*Z**2.39
c           ' The values in BEtab are the negative binding energies!
c           ' Pairing in Thomas Fermi masses is zero for Z,N even !
         If (BEtab.eq.0) Then
           BEtab = LyMass(Z,A,0.)
c         Print "Warning: Binding energy of Z=";Z;", A=";A;" not in mass table,";
c                        " replaced by LYMASS"
         End If
          ALDMAS = BEtab + FEDEFOLys(Z,A,beta)
            return
            end

            Function Ushell(Z,A)
      common /datain/beldm(203,136),ushel(203,136),RNucTab(3885,8)
            AN = A - Z
            AN = max(AN,1.)
            Res = UShel(nInt(AN),nInt(Z))
            If(Res.gt.0) Then
              Res = 0.3 * Res
            endif
c     ' KHS (12. Feb. 2012)
c     '      ' The positive shell effects for deformed nuclei seem to be too positive
c            ' This gives too many high-energetic prompt neutrons.
          USHELL = Res
          return
          End

            Function UMass(Z,A)
            BE = ALdmas(Z,A,0.) + USHELL(Z,A)
            UMASS = BE
            return
            End

        Function Zequi(ZCN,A1,A2,beta1,beta2,d,Imode)
c    /' Determines the minimum potential of the scission-point configuration
c       represented by two deformed nuclei divided by a tip distance d.
c       A1, A2, beta1, beta2, d are fixed, Z1 is searched for and returned on output.  '/

c       /' ZCN: Z of fissioning nucleus '/
c       /' A1: A of first fission fragment '/
c       /' A2: A of second fission fragment '/
c       /' beta1: deformation of first fission fragment '/
c       /' beta2: deformation of second fission fragment '/
c       /' d: tip distance '/

         real*4 Lymass
         common /comvar/ZUCD
        data POLARadd/0.3/,POLARfac/1./

c             Dim As Single RZequi
c             Dim As Single ACN
c             Dim As Single Z1UCD,Z2UCD
c             Dim As Single re1,re2,re3,eps1,eps2,DZPol /' help variables '/

c         write(*,*) ZCN,A1,A2,beta1,beta2,d,Imode

c         ACN = A1 + A2
          Z1UCD = A1 / (A1 + A2) * ZCN
          Z2UCD = ZCN - Z1UCD
          re1 = Lymass( Z1UCD-1.E0, A1, beta1 ) +
     ;          Lymass( Z2UCD+1.E0, A2, beta2 ) +
     ;          ECoul( Z1UCD-1.E0, A1, beta1,Z2UCD+1.E0, A2, beta2, d)
          re2 = Lymass( Z1UCD, A1, beta1) +
     ;          Lymass( Z2UCD, A2, beta2) +
     ;          ECoul( Z1UCD, A1, beta1,Z2UCD, A2, beta2, d )
          re3 = Lymass( Z1UCD+1.E0, A1, beta1 ) +
     ;          Lymass( Z2UCD-1.E0, A2, beta2 ) +
     ;           ECoul( Z1UCD+1.E0, A1, beta1,Z2UCD-1.E0, A2, beta2, d)
          eps2 = ( re1 - 2.E0*re2 + re3 ) / 2.E0
          eps1 = ( re3 - re1 ) / 2.E0
          DZPol = -eps1 / ( 2.E0 * eps2 )

c           write(*,*) ' DZPol=',DZPol

          If (Imode .gt. 0) Then
c            /' Purely empirical enhancement of charge polarization '/
            DZPOL = DZPOL * POLARfac + POLARadd
          End If

          RZequi = ZUCD + DZPOL
          Zequi = RZequi
        return
        End

        Function betalight(Z)
c      /' Deformation of light fission fragment for S1 and S2 '/
c      /' Systematic correlation Z vs. beta for deformed shells '/
c      /' Z of fission fragment '/
        data betaL0/26.6/, betaL1/0.80/
        beta = (Z - betaL0) * betaL1/20.E0
        betalight = beta
        return
        End


        Function betaheavy(Z)
c      /' Deformation of heavy fission fragment for S2 '/
c      /' Systematic correlation Z vs. beta for deformed shells '/
c      /' Z of fission fragment '/
        data betaH0 /48.0/,betaH1/0.70/
        beta = (Z - betaH0) * betaH1/20.E0
        betaheavy = beta
        return
        End

        Function AMasscurv(Z,A,RI)
        Data kappa/0/
c     /'  Fit to  Data of Fig. 7 of                                             '/
c     /'  "Shell effects in the symmetric-modal fission of pre-actinide nuclei" '/
c     /'  S. I. Mulgin, K.-H. Schmidt, A. Grewe, S. V. Zhdanov                  '/
c     /'  Nucl. Phys. A 640 (1998) 375                                          '/

      ZsquareoverA = Z**2/A
      ZsqrA = ZsquareoverA * (1.E0 - kappa * RI**2) /
     ;  (1.E0 - kappa * ((226.E0 - 2.E0*91.E0)/226.E0)**2)

      Result1 = F1(ZsqrA)
      Result2 = F2(ZsqrA)

c      write(*,*) result1,result2

      Result = Min(Result1,Result2)
      AMasscurv = Result
      return
      End

      Function DeSaddleScission(ZsquareoverAthird)
c    /' Energy release between saddle and scission '/
c    /' M. Asghar, R. W. Hasse, J. Physique C 6 (1984) 455 '/
      data ESHIFTSASCI/-30./
      Result = (31.E0 - 11.E0) / (1550.E0 - 1300.E0) *
     ;        (ZsquareoverAthird - 1300.E0 + ESHIFTSASCI) + 11.E0
      Result = max(Result,0.)
      DeSaddleScission = Result
      return
      End

        Function TRusanovf(E, A)
c     /' Fermi-gas level density, parameterisation of Rusanov et al. '/
       If( E .gt. 0) Then
         TRusanovf = sqrt(E / (0.094E0 * A) )
       Else
         TRusanovf = 0.
       End If
       return
       End

      Function IVENODD(RORIGIN,REVENODD)

c     /' Procedure to calculate IOUT from RIN in a way that         '/
c     /' on the average a flat distribution in RIN results in a      '/
c     /' fluctuating distribution in IOUT with an even-odd effect as '/
c     /' given by REVENODD                                          '/
c
c     /' ------------------------------------------------------------ '/
c     /' EXAMPLES :                                                   '/
c     /' ------------------------------------------------------------ '/
c     /'    If REVENODD = 0 :                                       '/
c     /'           CEIL(RIN)  ----                                   '/
c     /'                                                              '/
c     /'              RIN ->                                         '/
c     /'            (somewhere in between CEIL(RIN) and FLOOR(RIN)) '/
c     /'                                                              '/
c     /'           FLOOR(RIN) ----       --> IOUT                   '/
c     /' ------------------------------------------------------------ '/
c     /'    If REVENODD > 0 :                                       '/
c     /'      The interval for the above treatment is                 '/
c     /'         larger for FLOOR(RIN) = even and                    '/
c     /'         smaller for FLOOR(RIN) = odd                        '/
c     /'    For REVENODD < 0 : just opposite treatment              '/
c     /' ------------------------------------------------------------ '/
c
c     /' ------------------------------------------------------------ '/
c     /' On input:   RORIGIN    nuclear charge (real number)         '/
c     /'             REVENODD  requested even-odd effect            '/
c     /' Intermediate quantity: RIN = RORIGIN + 0.5                 '/
c     /' On output:  IOUT       nuclear charge (integer)             '/
c     /' ------------------------------------------------------------ '/
        integer rfloor

        REVENODD = MIN(REVENODD,1.)
        RIN = RORIGIN + 0.5E0
        RFLOOR = INT(RIN)
         If (Abs(REVENODD) .lt. 1.E-3) Then
          IOUT = RFLOOR
         Else
        RREST = RIN - RFLOOR
        RMIDDLE = RFLOOR + 0.5E0
        IF (RFLOOR - (RFLOOR/2)*2 .eq. 0) THEN
c /' even before modif. '/
         RHELP = RMIDDLE + (RREST - 0.5E0)
     ;      * 1.E0 / Max(0.01,(1.E0 + REVENODD))
         RHELP = Min(RHELP,RMIDDLE+1)
         RHELP = Max(RHELP,RMIDDLE-1)
        ELSE
c  /' odd before modification '/
          RHELP = RMIDDLE + (RREST - 0.5E0)
     ;      * 1.E0 / Max(0.01,(1.E0 - REVENODD))
         RHELP = Min(RHELP,RMIDDLE+1)
         RHELP = Max(RHELP,RMIDDLE-1)
        EndIf
         IOUT = INT(RHELP)
       EndIf
        Ivenodd = IOUT
       return
       End

       Function Erf(x)
c   /' Sergei Winitzki, 2008: relative accuracy < 1.4E-4 '/
c   '  Dim As Single a = 0.147
c   '  Const As Single b = 1.27324  ' 4/pi,
c   '  Const As Single pi = 3.14159
c   '  Dim As Single Result
c   '  Result = Sqr(1.E0 - exp(-x^2 * (b + a * x^2) / (1.E0 - a * x^2) ) )
c   '  If x < 0 Then Result = - Result
c   '  WiErf = Result
       external Erfc
       Erf = 1.E0 - Erfc(x)
       return
       End

       Function Erfc(x)
c   /' Complementary error function from numerical recipes '/
       z = Abs(x)
       t = 1.E0/(1.E0 + 0.5E0 * z)
       r = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196
     ; +t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398
     ; +t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
       If (x .lt. 0) Then
        Erfc = 2.E0 - r
       Else
        Erfc = r
       Endif
       return
       End

       Function fTanh(x)
       Result = (exp(x) - exp(-x)) / (exp(x) + exp(-x))
       fTanh = Result
       return
       End

      Function PBox(AMean,Sigma,Bottom)
      integer Bottom
c   ' Rectangular distribution folded with a Gaussian distribution
      R = PGauss(AMean,Sigma)
      R = R + (Rndm(-1.)-0.5)*Bottom
      PBox = R
      return
      End

c     FUNCTION  PGAUS(AV,SD)
c     A=-6.
c     DO 1 I=1,12
c1    A=A+RNDM(-1.)
c     PGAUS=A*SD+AV
c     RETURN
c     END

      Function PGauss(AMean,Sigma)
c  ' Box-Mueller method
      data Iset/1/,GSet/0.5/
      If (ISet .eq. 0) then
10       continue
        V1 = 2.E0 * Rndm(-1.) - 1.E0
        V2 = 2.E0 * Rndm(-1.) - 1.E0
        R = V1**2 + V2**2
       If (R .ge. 1.E0) Then
        Goto 10
       end if
        Fac = Sqrt(-2.E0 * aLog(R)/R)
        GSet = V1 * Fac
        GasDev = V2 * Fac
        ISet = 1
       Else
        GasDev = GSet
        ISet = 0
      End if
       Result = Sigma * GasDev
       PGauss = AMean + Result
      return
      End

      Function PLinGauss(RSigma)
c  /' Random-number generator for linear * Gaussian function '/
c  /' Distribution of nuclear angular momenta '/
      Brms = RSigma / sqrt(2.)
c    ' Because
c    ' the sum of two PGauss functions increases the width.
      RRes = Abs(PGauss(0.,Brms)) + Abs(PGauss(0.,Brms))
      if (RSigma.gt.0.) RRes=RRes+Brms/4.E0*(1.E0 - exp(-RRes/RSigma))
c    ' correction of shape (approximative)
      PLinGauss = RRes
      return
      End

      Function Ualev(Z,A)
c    '  U_alev = 0.073 * A + 0.095 * A^0.666667  'Ignatyuk (1970's)
       Ualev = 0.078 * A + 0.115 * A**0.6666667
c ' Ignatyuk (Bologna 2000)
c    '  U_alev = 0.089 * A    ' only volume term
      return
      End

       Function UTemp(IZ, IA, E)
c       ' Temperature (modified Gilbert-Cameron composite level density)
c       ' KHS (10. 2. 2012)
       data fgamma /0.055/,Tscale/0.85/,Econd/2./
c       ' Used global parameters: Tscale
c    '  alev = U_alev(Z,A) * 1.1   ' Factor adjusted to high-energy prompt neutrons in U235(nth,f)
       alev = Ualev(Float(IZ),Float(IA)) * 0.86
c  ' " with the correction for non-constant T (FG range)
c    '  alev = U_alev(Z,A)

       RShell = UShell(Float(IZ),Float(IA))
       TCT = TEgidy(Float(IA),RShell,Tscale)
       expo=-fgamma*E
       expo=min(expo,80.)
       Eeff0 = E - Econd + Aypair(Float(IZ),Float(IA))
     ; + Rshell*(1-exp(expo))

       TFG=0.
       If (Eeff0 .gt. 0.5) Then
         Eeff1 = Eeff0 + 0.1
         expo = 2.E0 * sqrt(alev * Eeff0)
         if (expo.lt.80.) then
           Rho0 = 1.E0/Eeff0**1.25 * exp(expo)
           Rho1 = 1.E0/Eeff1**1.25 * exp(2.E0 * sqrt(alev * Eeff1))
c '         Rho0 = 1.E0/Eeff0 * exp(2.E0 * sqrt(alev * Eeff0))
c '         Rho1 = 1.E0/Eeff1 * exp(2.E0 * sqrt(alev * Eeff1))
           TFG = 0.1E0 / (alog(Rho1) - alog(Rho0))
         endif
       Else
         TFG = 0
       End If
       Res = TCT
       If (TFG .gt. Res) Then
         Res = TFG
       end if

c ' If Res > 1.4 Then Res = 1.4

       UTemp = Res
       return
       End

      Function PEN1(Td,TCT,Ed,Ad)
c   /' Random-number generator for a surface Maxwell distribution '/
c   /' y = x * exp(-x/T)  '/
c   /' With correction for non-constant T in FG regime and '/
c   /' with energy-dependent inverse n cross section '/
c   /' From optical model of Koning and Delaroche '/
c   /' (analytical approximation by KHS, 2012) '/
       Data Ecrit/6./
c  ' Ecrit above basis of FG level-density, average value
c   Repeat_PEN1:
 100   continue
 200   continue
c   /' Pure Maxwell '/
      E = -Td * (aLog(Rndm(-1.)) + aLog(Rndm(-1.)))
c  ' pure Maxwell

c   /' Correction for non-constant T in FG regime '/
      If (E .lt. Ed) Then
        If (Ed .gt. Ecrit) Then
c   ' reduction for non-constant T in FG regime
        Efinal = Ed-E
         If (Efinal .gt. Ecrit) Then
          Fred = exp(2./Td*sqrt(Ed*Efinal)+(E-2*Ed)/Td)
         Else
       Fred = exp(2*sqrt(Ed*Ecrit)/Td-(Ecrit-Efinal)/TCT +(E-2*Ed)/Td)
         End If
       If (Rndm(-1.) .gt. Fred) Then
         Goto 200
       end if
       End If
      End If
c   /' Optical model: '/
      x = E / Ad**0.6666667
      if (x .lt. 4.5) then
       y =155.*x**(-0.14)+15.*x**(-0.14)*sin((3*x**0.14+0.06)/0.115-4.3)
      else
       Y = 290 * exp(-alog10(x)/alog10(6.))
      end if
      R = 400 * Rndm(-1.)
      If (R .gt. y) Then
       Goto 100
      endif

      PEN1 = E
      return
      end

      Function PEgammahigh(Zi, Ai, Ei)
c   ' Random function, returns gamma energy in MeV
c   ' From PRL 49 (1982) 434
c   ' For energies above Sn: competition with neutrons included
      dimension sigma(5000)

       N = min(nInt(Ei*10),5000)

       do i=1,5000
         sigma(i)=0.
       enddo
        If (N .le. 0) Then
         N = 1
        end if

c     ReDim As Single sigma(N)  ' sigma is not normalized
c                               ' (Normalization is done by Monte-Carlo procedure.)

       E0 = E0GDR(Zi,Ai)
       G0 = WidthGDR(E0)
       Tm = Utemp(nInt(Zi),nInt(Ai),Ei)

c     ' Establish distribution
       sigMax = 0.
c       For Eg = 0.1 To Ei Step 0.1
        Eg=0.1
       do i=1,N
        M = nInt(Eg*10.)
        sigma(M) = Eg**3 / Tm**2 * exp(-Eg/Tm) *
     ;      (G0 * Eg) / ( (Eg**2 - E0**2)**2 + G0**2 * E0**2 )
       if (sigma(M) .gt. sigMax) then
        sigMax = sigma(M)
       end if
        Eg=0.1*i
       end do

c    ' Dice gamma energy from distribution
100    continue
        xran = max( rndm(-1.) * Ei * 10, 1.)
        M = min (nint(xran),N)
c   ' in units of 100 keV
        yran = rndm(-1.) * sigMax
      if (yran .gt. sigma(M)) then
       goto 100
      endif

      PEgammahigh = xran/10.
c ' convert to MeV
      return
      End

      Function E0GDR(Z,A)
c     ' Calculates the centroid energy of the GDR for spherical nucleus
c     ' according to the FRDM (ADNDT 59 (1995) 185 and PLB 670 (2008) 200)
c      Data epsilon /0.0768/
c     Static As Single J = 32.7
c     Static As Single Q = 29.2
c     Static As Single R0 = 1.16
c     Static As Single mstar = 874
c     Static As Single hbar = 197.3
c     Dim As Single Aonethird,u,N,E0

c     ' according to [9] in Phys. Lett. B 690 (2010) 473:
       E0GDR = 18.0/A**0.333333 + 25.0 / A**0.1666667

c     ' according to the FRDM (ADNDT 59 (1995) 185 and PLB 670 (2008) 200):
c  '    Aonethird = A^0.333333
c  '    N = A - Z
c  '    u = (1-epsilon)/Aonethird * 3*J/Q
c  '    E0_GDR = hbar /(R0*Aonethird)*sqr(8*J*A^2/ (mstar*4*N*Z) ) * _
c  '      (1 + u - epsilon * (1+epsilon+3*u)/(1+epsilon+u))^(-1/2)
      return
      End

      Function WidthGDR(E0)
c     ' Spreading width of the GDR (Nucl. Phys. A 531 (1991) 27)
      WidthGDR = 1.99 * (E0/10.)**1.6
      return
      End

      Function IMATENDF(IZ,IA)
      common /datain/beldm(203,136),ushel(203,136),RNucTab(3885,8)
      Do I = 1,3885
      IMAT = I
      If (IZ .eq. nInt(RNucTab(IMAT,2)) .And.
     ; IA .eq. nInt(RNucTab(IMAT,3))) Then
       Exit
      End if
      End do
      IMATENDF = IMAT
      Return
      End
      function rndm(r)
      implicit double precision (a-h,o-z)
      common /rand/   rani,ranj,rijk,ranb,rans,nrand
      real*4 r,rndm
      real*8 ran2
      rndm=sngl(ran2())
C      write(6,*) rndm
      return
      end
*-------------------------------------
c     subroutine rdmout(k)
c     return
c     end
*-------------------------------------
c     subroutine rdmin(k)
c     return
c     end
*-------------------------------------

*
c     common /rand/   rani,ranj,rijk,ranb,rans,nrand
*
*
*
c     nrand=0
c     r=5.d0**19
c     call setijk(r)
c     call setran
*
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  r a n d o m  *-*-*-*-*-*-*-*-*-*-*-*-*-*-
*
*     dieser generator wird auch in mcnp benutzt
*
      subroutine random
      implicit double precision (a-h,o-z)
      common /rand/   rani,ranj,rijk,ranb,rans,nrand
      parameter (p=2d0**24,q=2d0**(-24),fb=4867484d0,fs=10256733d0)
*
      entry setijk(r)
      rani=aint(r*q)
      ranj=r-rani*p
      rijk=rani*p+ranj
      return
*
      entry advijk
      a=fs*ranj
      b=fb*ranj+fs*rani+aint(a*q)
      ranj=a-aint(a*q)*p
      rani=b-aint(b*q)*p
      rijk=rani*p+ranj
      return
*
      entry setran
      ranb=rani
      rans=ranj
      return
      end
*
*-----------------------------------------------------------------------
*
      double precision function ran2()
      implicit double precision (a-h,o-z)
      common /rand/   rani,ranj,rijk,ranb,rans,nrand
      parameter (p=2d0**24,q=2d0**(-24),r=2d0**(-48),
     *           gb=1136868d0,gs=6328637d0)
*
      a=gs*rans
      b=gb*rans+gs*ranb+aint(a*q)
      rans=a-aint(a*q)*p
      ranb=b-aint(b*q)*p
      ran2 =(ranb*p+rans)*r
      nrand=nrand+1
      return
      end
*--------------------------------------------------------
c     n=n+1
c     read(lunc)
c    *      rani1,ranj1,rijk1,ranb1,rans1,nrand1,
c     if (cont.eq.1) then
c       rani=rani1
c       ranj=ranj1
c       rijk=rijk1
c       ranb=ranb1
c       rans=rans1
c       nrand=nrand1
c     endif
