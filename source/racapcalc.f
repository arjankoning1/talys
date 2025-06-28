      subroutine racapcalc(knp,ein,kz1,ka1,beta2tar,
     &  spinexptar,iparexptar,spinthtar,iparthtar,
     &  a1exp,a1th,kz2,ka2,spinexpproj,a2,
     &  spinexpcn,iparexpcn,spinthcn,iparthcn,afexp,afth,Sn,
     &  levmax,iexplvnum,explveng,explvsp,iexplvpar,exf,
     &  nmax,potjlmini,potjlmfin,jspmax,nldmax,rhonldposj,rhonldnegj,
     &  vncapini,rncapini,ancapini,vncapfnl,rncapfnl,ancapfnl,
     &  xsall,esig,xspopex,xspop,iopt,
     &  pi,e2,amu,hc,nlevtot,spfacst,ispect)
c
      implicit double precision(a-h,o-z)
      real ein
      integer kz1,ka1,kz2,ka2
      integer knp
      real spinexptar,spinthtar
      real spinexpproj
      real spinexpcn,spinthcn
      integer iparexptar,iparthtar
      integer iparexpcn,iparthcn
      double precision a1,a1exp,a1th,a2,af,afth,afexp
      real Sn
      integer levmax
      integer iexplvnum,iexplvpar(levmax)
      real explveng(levmax),explvsp(levmax)
      real exf(levmax)
      integer nmax,nlevtot
      real   potjlmini(nmax),potjlmfin(nmax)
      integer jspmax,nldmax
      double precision rhonldposj(nldmax,0:jspmax)
      double precision rhonldnegj(nldmax,0:jspmax)
      real   vncapini,rncapini,ancapini
      real   vncapfnl,rncapfnl,ancapfnl
      integer iopt
      real pi,e2,hc
      double precision amu
      real beta2tar
      real xspopex(levmax),xspop(levmax,0:jspmax,2)
      real ecm,rad1,fnorm
      real xsall(3),esig
      real spfacst(levmax)
      real rr(nmax)
      double precision rhodens(nmax)
      dimension potf(0:nmax),wff(0:nmax),poti(nmax),wfi(nmax)
      dimension vi(nmax),vf(nmax),vinit(nmax)
      dimension soinit(nmax),sof(nmax)
      dimension rho(2,3),rho1(2,3),rhof(2,3),nz(2)
      dimension crossx(3),sfactor(3)
      dimension spinff(levmax),ipaff(levmax)
      dimension rhobin(levmax,0:jspmax,2)
      dimension n0n(33),ln(33),nn(33)
      dimension lnf(33)
c
      data n0n/
     &0,0,0,0,0,1,0,0,1,1,0,0,1,1,0,2,0,1,0,1,2,2,0,1,0,1,2,2,3,0,1,0,1/
      data ln/
     &0,1,1,2,2,0,3,3,1,1,4,4,2,2,5,0,5,3,6,3,1,1,6,4,7,4,2,2,0,7,5,8,5/
      data nn/
     &  2,  6,  8, 14, 18, 20, 28, 34, 38, 40, 50, 58, 64, 68, 80, 82,
     & 92,100,114,120,124,126,138,148,164,172,178,182,184,198,210,228,
     &238/
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     set some control value
c
      ipot=5
      xsall(1)=0.d0
      xsall(2)=0.d0
      xsall(3)=0.d0
      nlevtot=0
      facv1=0.
      ebound=0.
c
c     ipot is NN interaction
c     ipot=1 - Gaussian potential
c     ipot=2 - M3Y (Paris)
c     ipot=3 - M3Y (OPEP)
c     ipot=4 - M3Y (Reid-standard, energy and density independent) DEFAULT
c     ipot=5 - DDM3Y (Reid, energy and density dependent)
c
      iopnorm=0
c
c     iopnorm is the normalized on volume interaction of selected 'iopt',
c     here, iopnorm=0 means no normalization
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     set the radial step and limit for potential calculation,
c     radial wavefunction, and other functions related to radial etc...
c     n is the number of discretization points.
c     h is the step in fm.
c     rlim=20fm is the total range.
c     here, n, h, rlim should be exactly the same as the discretization
c     points of JLMB potential
c
      n=nmax
      rlim=20.
      h=rlim/dble(n)
      do i=1,n
        rr(i)=rlim*dble(i)/dble(n)
      end do
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     A+a=gamma+B
c     ka1 and kz1 is the target nuclei(A),ka2 and kz2 is the projectile nuclei(a)
c     kaf and kzf is the final nuclei(B=A+a)
c     a1,a2=masses (real)
c     kz1,kz2=charges number (integer)
c     ka1,ka2=mass number (integer)
c
      a1=a1exp
      af=afexp
      kaf=ka1+ka2
      kzf=kz1+kz2
      if (knp.eq.1) nnf=kaf-kzf
      if (knp.eq.2) nnf=kzf
c     knp=1 is for neutron capture, nnf should be the total neutron number in the final nuclei;
c     knp=2 is for proton capture, nnf should be the total proton number in the final nuclei.
c     nnf is used for the selection rule in the following code.
c
c     set of experimental and theoretical real masses: experimental mass frist.
c     if experimental mass not available, theoretical mass is used.
c
      if (a1.eq.0.0) a1=a1th
      if (af.eq.0.0) af=afth
      if (a1.eq.0.0.or.af.eq.0.0) return
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     ecm=cm incident energies r
c
      ecm=ein*a1/(a2+a1)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     set the coulomb radii and charged, reduced mass, and odd-even functions(lz,lz1,lz2).
c
c     hc=hbar*c  (in Mev*fm)
c     fsc=inverse of the fine-structure constant
c     amn=amu*a2, mean nucleon mass (in Mev c**(-2)), a2 is real neutron or proton mass
c     hm=hc**2/2*amn
c
      rc=1.2*(dble(ka1)**(1./3.)+dble(ka2)**(1./3.))
      n1=a1+0.5d0
      n2=a2+0.5d0
      lz=2-mod(n1+n2,2)
      lz1=2-mod(n1,2)
      lz2=2-mod(n2,2)
c
      fsc=hc/e2
      amn=amu*a2
      hm=hc*hc/(2.0*amn)
      ama=a1+a2
      rmu=a1*a2/ama
      rm=rmu/hm
      ze=kz1*kz2*e2
      eta0=kz1*kz2*sqrt(rmu*amn/2.)/fsc
c
c     set the spin and parity of target,projectile and G.S. of final nuclei.
c
      spin1=90.0
      ipa1=9
      spin2=0.5
      spinf=90.0
      ipaf=9
c
      spin1=spinthtar
      if (spinexptar.lt.90.) spin1=spinexptar
      ipa1=iparthtar
      if (iparexptar.lt.9) ipa1=iparexptar
      spin2=spinexpproj
c
      i1=nint(spin1*(mod(ka1,2)+1))
      i2=nint(spin2*(mod(ka2,2)+1))
      i1i2=(lz1*i1+1)*(lz2*i2+1)
c
      spinf=spinthcn
      if (spinexpcn.lt.90.) spinf=spinexpcn
      ipaf=iparthcn
      if (iparexpcn.lt.9) ipaf=iparexpcn
c
c     check the Q-value, only Q-value > 0 will be considered for calculation
c
      if (Sn.le.0.0) return
c
c Initialization of the spin-orbit part of the potential
c
      do 385 j=1,n
      soinit(j)=0.
  385 sof(j)=0.
c
c     if (iopt.eq.1.or.iopt.eq.3) goto 13
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     determination of the matter density distributions
c     and calculation of the nuclear potential (only for M3Y)
c
      do 290 i=1,2
      do 290 j=1,3
        rho(i,j)=0.
        rho1(i,j)=0.
  290   rhof(i,j)=0.
c     use of R=1.322 A^1/3 and a=0.66 for matter density
      do 325 iif=1,2
        aa=dble(ka1)+iif-1
        do 325 inz=1,2
          rad=(1.322-7.6d-4*aa+4.d-6*aa*aa-8.d-9*aa*aa*aa)*aa**(1./3.)
        dif=0.66
        if (iif.eq.1) then
          rho1(inz,2)=rad
          rho1(inz,3)=dif
        endif
        if (iif.eq.2) then
          rhof(inz,2)=rad
          rhof(inz,3)=dif
        endif
  325 continue
c
c normalization of density distribution
c
      do 340 iif=1,2
        nz(1)=(ka1-kz1)+iif-1
        nz(2)=kz1
        do 350 inz=1,2
          do 360 iij=1,3
          if (iif.eq.1) rho(inz,iij)=rho1(inz,iij)
  360     if (iif.eq.2) rho(inz,iij)=rhof(inz,iij)
          x=0.
          do 370 ir=1,n
            rhodens(ir)=1./(1.+dexp((x-rho(inz,2))/rho(inz,3)))
  370     x=x+h
          call intu1(rhodens,n,h,volj,rx,bull,3)
          rho(inz,1)=nz(inz)/volj
          if (iif.eq.1) rho1(inz,1)=rho(inz,1)
          if (iif.eq.2) rhof(inz,1)=rho(inz,1)
  350   continue
  340 continue
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     computing the initial nuclear potential (vinit)
c
c     normalization of folding potential to the volume integral of iopnorm
c     volj corresponds to the volume integral of the chosen potential
c     volj3 corresponds to the volume integral of the normalizing potential
c
      do 380 inz=1,2
      do 380 ii=1,3
  380 rho(inz,ii)=rho1(inz,ii)
c
c eave=100 keV average energy to calculate the folding potential
c assumed to be energy-independent, as well as its volume integral !!
c
      eave=0.100
      if (iopt.eq.4.or.iopnorm.eq.4) then
        call foldpot(eave,rho,ka1,kz1,vi,h,volj,ipot)
        if (iopnorm.eq.4) volj3=volj
      endif
      e0=ecm
      if (iopt.eq.4) then
        do 395 j=1,n
  395   vinit(j)=vi(j)
      endif
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  iopt=1 - WS-type potential
c  iopt=3 - JLMB potential
c
c  13 continue
c
      if (iopt.eq.1) then
        aa=dble(ka1)
        rad1=rncapini*aa**(1./3.)
        do j=1,n
           vi(j)=dble(-vncapini/(1.0+
     &       exp((rr(j)-rad1)/ancapini)))
           vinit(j)=vi(j)
        enddo
        call intu1(vi,n,h,volj,rx,bull,3)
        volj=volj/dble(ka1)
      endif
c
      if (iopt.eq.3) then
        do j=1,n
          vi(j)=potjlmini(j)
          vinit(j)=potjlmini(j)
        enddo
        call intu1(vi,n,h,volj,rx,bull,3)
        volj=volj/dble(ka1)
      endif
c
      if (iopnorm.ne.0.and.iopnorm.ne.4.and.iopnorm.ne.5.
     &  and.iopt.ne.iopnorm) then
        if (iopt.eq.1) then
          aa=dble(ka1)
          rad1=rncapini*aa**(1./3.)
          do j=1,n
            vi(j)=dble(-vncapini/(1.0
     &        +exp((rr(j)-rad1)/ancapini)))
            vinit(j)=vi(j)
          enddo
        endif
        if (iopt.eq.3) then
          do j=1,n
            vi(j)=potjlmini(j)
          enddo
        endif
        call intu1(vi,n,h,volj3,rx,bull,3)
        volj3=volj3/dble(ka1)
      endif
c
      fnorm=1.
      if (iopnorm.ne.0.and.iopnorm.ne.iopt) fnorm=volj3/volj
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     computing the final potential (vf)
c
c     M3Y potential
      if (iopt.eq.4) then
        e0=0.0
        do 400 inz=1,2
        do 400 ii=1,3
  400   rho(inz,ii)=rhof(inz,ii)
        call foldpot(e0,rho,kaf,kzf,vf,h,volj,ipot)
      endif
c     Woods-Saxon Potential
      if (iopt.eq.1) then
        aa=dble(kaf)
        radf=rncapfnl*aa**(1./3.)
        do j=1,n
          vf(j)=dble(-vncapfnl/(1.0+exp((rr(j)-radf)/ancapfnl)))
        enddo
      endif
c     JLMB potential
      if (iopt.eq.3) then
        do j=1,n
          vf(j)=potjlmfin(j)
        enddo
      endif
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     here begin to set the experimental levels and levels from theoretical level density,
c     also set the quantum numbers (ln, nn, n0n from shell model data), for calculating the Schodinger Equation.
c
c     Determination of the experimental or theoretical excitation spectrum, with
c     excitation energy, spin, parity and spectroscopic factor of final nucleus states.
c
      do 700 ken=1,levmax
      do 700 jspin=0,jspmax
      do 700 jparity=1,2
  700 rhobin(ken,jspin,jparity)=0.d0
c
c   spinff(i), ipaff(i), and exf(i) are array, for all exited levels of final states (spin,parity,exited energy).
c   Here, spinf and ipaf are the number for only G.S. of final states (spin,parity).
c   spinff(1), ipaff(1), and exf(1) are the G.S spin, parity and level position (exf(1)=0.00).
c
      if (spinf.lt.90.) spinff(1)=spinf
      if (ipaf.lt.9) ipaff(1)=ipaf
      nlevf=1
      exf(1)=0.0
      do i=1,iexplvnum-1
        if (explveng(i).lt.Sn) then
          exf(i+1)=explveng(i)
          spinff(i+1)=explvsp(i)
          ipaff(i+1)=iexplvpar(i)
          nlevf=nlevf+1
          if (nlevf.ge.levmax) exit
        endif
      enddo
      if (ispect.eq.2) nlevf=1
      emaxex=exf(nlevf)
      nlevf=min(nlevf,levmax)
      iexplvnum=nlevf
c
      do 710 jlev=1,nlevf
      jspin=idnint(spinff(jlev)-dble(mod(kaf,2))/2.)
      jspin=min(jspin,jspmax)
      jspin=max(jspin,0)
      jparity=min(idnint(ipaff(jlev)/2.d0+1.5d0),2)
  710 rhobin(jlev,jspin,jparity)=1.d0
      nlevfmax=nlevf
c
c only experimental levels are considered, no theoretical NLD levels are used.
c
      if (ispect.eq.1) goto 205
      if (nlevfmax.ge.levmax) goto 205
c
c construction of combinatorial level density in bin of Qvalue/nbin from ethmin to ethmax.
c
      bin=0.25
      kethmin=idnint(emaxex/bin+0.499)+1
      ethmin=dble(kethmin-1)*bin
      kethmax=int(Sn/bin)
      kethmax=min(kethmax,nldmax)
      nlevfmax=min(nlevf+kethmax-kethmin+1,levmax)
c
c Theoretical NLD from TALYS above the experimental level "nlevf"
c Energy range splitted in "nbin" bins of dE="bin"
c
c exclude ground state from theoretical NLD:
c
      do j=0,jspmax
        if (rhonldposj(1,j).ge.1.) rhonldposj(1,j)=rhonldposj(1,j)-1.0
        if (rhonldnegj(1,j).ge.1.) rhonldnegj(1,j)=rhonldnegj(1,j)-1.0
      enddo
c
c set the theoretical levels from NLD in given bin:
c
      do 730 keth=kethmin,kethmax
        jlev=nlevf+(keth-kethmin+1)
        if (jlev.gt.levmax) then
          jlev=levmax
          goto 730
        endif
        exf(jlev)=ethmin+bin*(keth-kethmin+0.5)
        do j=0,jspmax
          rhobin(jlev,j,1)=rhonldposj(keth,j)*bin
          rhobin(jlev,j,2)=rhonldnegj(keth,j)*bin
c         if (j.lt.jspmax-1) then
c           rhobin(jlev,j,1)=rhonldposj(keth,j)*bin-
c    +      rhonldposj(keth,j+1)*bin
c           rhobin(jlev,j,2)=rhonldnegj(keth,j)*bin-
c    +      rhonldnegj(keth,j+1)*bin
c         else
c           rhobin(jlev,j,1)=0.0
c           rhobin(jlev,j,2)=0.0
c         endif
          if (rhobin(jlev,j,1).lt.0.0) rhobin(jlev,j,1)=0.0
          if (rhobin(jlev,j,2).lt.0.0) rhobin(jlev,j,2)=0.0
        enddo
  730 continue
  205 continue
      nlevtot=nlevfmax
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     Loop on all the available excited states in final nucleus for racap calculation
c
c     lam=multipole order
c     i1,i2=spins of target and projetile (exact if integer, twice if half-integer)
c     rc=coulomb radius (in fm)
c     lf,li=orbital momentum of the final and initial state (integer)
c     ii=channel spin (arising from the coupling of i1 and i2)
c     ji,jf=total spin
c     n0=number of nodes in the final wave function
c
      numf=0
      esig=0.d0
      sfactor(1)=0.d0
      sfactor(2)=0.d0
      sfactor(3)=0.d0
      do jlev=1,nlevfmax
        xspopex(jlev)=0.
        do jspin=0,jspmax
        do jparity=1,2
          xspop(jlev,jspin,jparity)=0.
        enddo
        enddo
      enddo
c
      do 210 jlev=1,nlevfmax
      esigex=0.
      ef=exf(jlev)
      crossx(1)=0.d0
      crossx(2)=0.d0
      crossx(3)=0.d0
      if (ef.gt.Sn) then
        write(*,*) 'E*=',ef,' >  Q-value=',Sn
        return
      endif
c
      do 221 jspin=0,jspmax
      do 221 jparity=1,2
      esigjp=0.
      spinf=dble(jspin+dble(mod(kaf,2))/2.)
      ipaf=int((dble(jparity)-1.5)*2.)
      if (jlev.gt.1) then
        weight=rhobin(jlev,jspin,jparity)*(0.1d0+0.33d0*dexp(-0.8*ef))
      else
        weight=rhobin(jlev,jspin,jparity)*spfacst(jlev)
      endif
      if (weight.le.1.d-20) goto 221
c
c     remove the fake contributions that does not match to the known experimental levels
c
      if(ispect.ne.2) then
        if (jlev.le.iexplvnum) then
          if (spinff(jlev).ne.spinf.or.ipaff(jlev).ne.ipaf) goto 221
        end if
      end if
      if(ispect.eq.2) then
        if (jlev.eq.1) then
          if (spinff(1).ne.spinf.or.ipaff(1).ne.ipaf) goto 221
        end if
      end if
c
      if (abs(ipaf).ne.1)
     +  write(*,*) ' Problem in Parity of (',kzf,',',kaf,') G.S.=',ipaf
      jf=int(spinf*(mod(kaf,2)+1))
      jf=jf*lz
c      write(8964,*) jspin,jparity,weight
c
c  Loop on different transitions
c  lam=1 for E1
c  lam=2 for E2
c  lam=3 for M1
c
      lammax=3
      ipasslam=0
      nlam1=3
      nlam2=15
      do 220 jlam=1,lammax
      lam=jlam
      if (lam.eq.1) fac0=8.*pi*(lam+1)/dble(nlam1)**2./(lam*100*i1i2)
      if (lam.eq.2) fac0=8.*pi*(lam+1)/dble(nlam2)**2./(lam*100*i1i2)
      if (lam.eq.3) fac0=0.011*8.*pi*(1+1)/dble(nlam1)**2./(1*100*i1i2)
c
c     selection rules on spin distribution
c
      icont=1
      lfmin=int(abs(spinf-spin1)-0.5)
      palf=(-1)**lfmin
      ipalf=int(palf)
      if (ipa1*ipaf.ne.ipalf) lfmin=lfmin+1
      if (lfmin.lt.0) lfmin=lfmin+2
      if (lam.eq.1) limin=lfmin-lam
      if (lam.eq.2) limin=lfmin-lam
      if (lam.eq.3) limin=lfmin
      if (limin.lt.0) limin=limin+2
      if (lfmin.ne.int(abs(spinf-spin1)-0.5)) icont=2
      if (i1.eq.0.or.jf.eq.0) icont=1
      li=limin
      lf=lfmin
c
c  cut-off at li=licut. No contribution considered for partial waves with li>licut !
c
      if (lam.eq.1) licut=6-2*lam
      if (lam.eq.2) licut=6-2*lam
      if (lam.eq.3) licut=4
      if (licut.gt.100) then
        write(*,*) ' dimension in dephase and coufra must be extended'
        stop
      endif
      if (li.gt.licut) goto 220
      if(lf.lt.0) goto 220
      if(li.lt.0) goto 220
c
c     determination of the number of nodes of the radial wave function
c
      n0=0
      do 460 j=1,33
        lnf(j)=ln(j)
        if (lnf(j).eq.lfmin.and.nn(j).le.nnf) n0=n0n(j)+1
        if (lnf(j).eq.lfmin.and.nn(j).gt.nnf) then
          n0=n0n(j)
          goto 470
        endif
  460 continue
  470 continue
      n0init=n0
c
c     loop on the number of possible contribution icont
c
c     Ii     =  If = I1+I2 = spii = spin1+spin2
c     Ji=Ii+li, Jf=If+lf,
c     Ji     =  Jf       +    Lam  (EL or ML)
c     Ji=spini, Jf=spinf
c     par1*par2*parli=parJi=parJf*par(EL or ML)
c     par1*par2*parlf=parJf
c
      ipass=0
      jjmax=2*li+1
      do 500 ipos=1,icont
      ipassf=-1
      spii=spin1+spin2
      if (ipos.eq.2) spii=spin1-spin2
      if (icont.eq.1.and.(spii.lt.dabs(spinf-lf).or.spii.gt.spinf+lf))
     &  spii=spin1-spin2
      ii=int(spii*(mod(kaf,2)+1))
      ii=ii*lz
      ifv=ii
c
      do 501 jj=1,jjmax
      spini=spii+dble(li-jj+1)
      if (spini.lt.0) goto 501
      if (lam.eq.1) then
        if (spini.gt.(spinf+lam).or.spini.lt.dabs(spinf-lam)) goto 501
      endif
      if (lam.eq.2) then
        if (spini.gt.(spinf+lam).or.spini.lt.dabs(spinf-lam)) goto 501
      endif
      if (lam.eq.3) then
        if (spini.gt.(spinf+1).or.spini.lt.dabs(spinf-1)) goto 501
      endif
      ipassf=ipassf+1
      ji=int(spini*(mod(kaf,2)+1))
      ji=ji*lz
c
c     add the constrain of the coupling of li with Ii to Ji
c     add the constrain of the coupling of lf with If to Jf
c
      if (ji.gt.(ifv+2*li).or.ji.lt.abs(ifv-2*li).or.n0.lt.0) goto 501
      if (jf.lt.abs(2*lf-ifv).or.jf.gt.2*lf+ifv.or.n0.lt.0) goto 501
      if (ji.lt.abs(2*li-ii).or.ji.gt.2*li+ii.or.n0.lt.0) goto 501
c
      if (lam.eq.1) then
        if(lam*2.lt.abs(ji-jf).or.lam*2.gt.ji+jf) goto 501
      endif
c
      if (lam.eq.2) then
        if(lam*2.lt.abs(ji-jf).or.lam*2.gt.ji+jf) goto 501
      endif
c
      if (lam.eq.3) then
        if(2.lt.abs(ji-jf).or.2.gt.ji+jf) goto 501
      endif
c
      if (lam.eq.1) then
        call clebs(2*li,2*1,2*lf,0,0,0,c1)
        call sixj(jf,ji,2*1,2*li,2*lf,ii,s1)
        facree1=(kz1*(a2/ama)**1.+kz2*(-a1/ama)**1.)**2.
     &   *(2.*1.+1.)*(jf+1)*(ji+1)*(2*li+1)*(c1*s1)**2.
        fac1=fac0*facree1*weight
      endif
c
      if (lam.eq.2) then
        call clebs(2*li,2*2,2*lf,0,0,0,c2)
        call sixj(jf,ji,2*2,2*li,2*lf,ii,s2)
        facree2=(kz1*(a2/ama)**2.+kz2*(-a1/ama)**2.)**2.
     &  *(2.*2.+1.)*(jf+1)*(ji+1)*(2*li+1)*(c2*s2)**2.
        fac1re=fac0*facree2*weight
        if (lz1*i1.ge.2.and.li.eq.lf) then
          call sixj(lz1*i1,lz2*i2,ifv,ii,2*2,lz1*i1,estar1)
          call sixj(ji,jf,2*2,ifv,ii,2*lf,estar2)
          call clebs(lz1*i1,4,lz1*i1,lz1*i1,0,lz1*i1,ectar)
          call calmagelc (kz1,ka1,dmag1,qelc1,beta2tar)
          facine2=(ji+1)*(lz1*i1+1)*(ii+1)*(ifv+1)*5./16./pi
     &      *qelc1*qelc1*estar1**2.*estar2**2./ectar**2.
        elseif (lz1*i1.lt.2.or.li.ne.lf) then
          facine2=0.
        endif
        fac1in=fac0*facine2*weight
      endif
c
      if (lam.eq.3) then
        if (lf.eq.li) then
          call sixj(ji,jf,2,2*li,2*lf,ii,srem)
          facrem1=(kz1/a1*(a2/ama)**1.-kz2/a2*(-a1/ama)**1.)**2.
     &      *(jf+1)*3.*li*(li+1)*(2*li+1)*(ji+1)*srem**2.
          if(lz1*i1.eq.0) then
            factarm1=0.
          else
            call sixj(lz1*i1,lz2*i2,ifv,ii,2,lz1*i1,starm)
            call clebs(lz1*i1,2,lz1*i1,lz1*i1,0,lz1*i1,ctarm)
            call sixj(ji,jf,2,ifv,ii,2*lf,sinallm)
            call calmagelc (kz1,ka1,dmag1,qelc1,beta2tar)
            factarm1=(ji+1)*(lz1*i1+1)*(ii+1)*(ifv+1)*3./4./3.14159
     &        *dmag1*dmag1*starm**2.*sinallm**2./ctarm**2.
          endif
          if(lz2*i2.eq.0) then
            facprojm1=0.
          else
            call sixj(lz2*i2,lz1*i1,ifv,ii,2,lz2*i2,sprojm)
            call clebs(lz2*i2,2,lz2*i2,lz2*i2,0,lz2*i2,cprojm)
            call sixj(ji,jf,2,ifv,ii,2*lf,sinallm)
            call calmagelc (kz2,ka2,dmag2,qelc1,beta2tar)
            facprojm1=(ji+1)*(lz2*i2+1)*(ii+1)*(ifv+1)*3./4./3.14159
     &        *dmag2*dmag2*sprojm**2.*sinallm**2./cprojm**2.
          endif
          fac1=fac0*(facrem1+factarm1+facprojm1)*weight
        else
          fac1=0.
        endif
      endif
c      write(8965,*) jspin,jparity,weight,fac1
c
c     Determination of the final state characteristics.
c     Normalization of the final potential on the excitation energy
c
      if (ipass.ne.0.or.ipassf.ne.0) goto 570
      if (ipasslam.ne.0) goto 570
c
      zf=jf*(jf+2)-ifv*(ifv+2)-4*lf*(lf+1)
      s2=0.
      e=0.
      umax=0.
      ee0=0.
      facv=1.
      xfac=0.05
      fvmax=1.3
      if (knp.eq.2.and.Sn.le.0.9) fvmax=1.9
      jtestmax=int((fvmax-1.)/xfac)
      itest=0
      itestmax=20
      jtest=0
      itrue=1
  430 continue
      itest=itest+1
      ee1=ee0
      potf(0)=0.
      do 420 j=1,n
        x=j*h
        vn=facv*(vf(j)+sof(j)*zf/8.)
        vc=ze/x
        if(x.le.rc) vc=ze*(3-(x/rc)**2.)/2./rc
        vce=lf*(lf+1)/x/x
        potf(j)=(vn+vc)*rm+vce
        if(lam.ne.0.or.j.le.2) goto 420
        if(j.gt.2.and.potf(j-1).le.max(potf(j-2),potf(j)))
     &    goto 420
        umax=potf(j-1)
        imax=j-1
  420 continue
      if(lam.eq.0) then
        do j=1,imax
          potf(j)=potf(j)-umax
        enddo
        do j=imax+1,nmax
          potf(j)=0
        enddo
      endif
c
c     computing the final energy (ebound) and wave function (wff)
c
      eps=1.d-4
      call num1l(n,h,e,s2,potf,wff,n0,eps)
c     write(8976,*) n,h,e,s2,n0
c     write(8977,*),itest,facv,n0
      ebound=(e+umax)/rm
      eexcit=Sn+ebound
      ee0=eexcit
c
      if (dabs(ee0-ef).lt.100.*eps.and.n0.ne.-1) goto 490
      if (itest.gt.itestmax) then
c        write(8972,*) "No solution found at E*=',ef"
c  No solution found at E*=',ef
        goto 500
      endif
      if (n0.eq.-1) then
        jtest=jtest+1
        n0=n0init
        ee0=Sn
        facv=facv+xfac*2.
        if (jtest.gt.jtestmax) then
          goto 500
        endif
        goto 430
      endif
      facv0=facv
      if (itest.eq.1) ee1=ee0
      if (itrue.eq.1) then
        if (ee0.gt.ef.and.ee1.gt.ef) facv=facv+xfac
        if (ee0.lt.ef.and.ee1.lt.ef) facv=facv-xfac
        if ((ee0-ef)*(ee1-ef).lt.0..and.ee1.ne.ee0) itrue=0
      endif
      if (itrue.ne.1) then
        if (ee0.eq.ee1) then
          xfac=xfac/2.
          if (ee0.lt.ef) facv=facv-xfac
          if (ee0.gt.ef) facv=facv+xfac
        else
          facv=facv0+(facv1-facv0)/(ee1-ee0)*(ef-ee0)
        endif
      endif
      facv1=facv0
      goto 430
c      write(8973,*) facv1 
c
  490 continue
      if (facv.gt.fvmax.or.facv.lt.1./fvmax) then
c        write(8974,*) "facv.gt.fvmax.or.facv.lt.1./fvmax" 
        goto 500
      endif
c
      ipass=ipass+1
      ipasslam=1
  570 continue
c
c    Calculation at the incident energies Einc
c
      e0=ecm
      eta=eta0/sqrt(e0)
      qk=sqrt(e0*rm)
      if (lam.eq.1) then
        fac2=fac1*((e0-ebound)/hc)**3./(qk*qk*fsc)*sqrt(rmu*amn/e0/2.)
      endif
      if (lam.eq.2) then
        fac2in=fac1in*((e0-ebound)/hc)**5.
     &       /(qk*qk*fsc)*sqrt(rmu*amn/e0/2.)
        fac2re=fac1re*((e0-ebound)/hc)**5.
     &       /(qk*qk*fsc)*sqrt(rmu*amn/e0/2.)
      endif
      if (lam.eq.3) then
        fac2=fac1*((e0-ebound)/hc)**3./(qk*qk*fsc)*sqrt(rmu*amn/e0/2.)
      endif
      fac3=e0*exp(2.*pi*eta)
c      write(8966,*) jspin,jparity,weight,fac1,fac2,fac3
c
c     Determination of the initial state wave function
c
      zi=ji*(ji+2)-ifv*(ifv+2)-4*li*(li+1)
      do 520 j=1,n
        x=j*h
        vn=fnorm*(vinit(j)+soinit(j)*zi/8.)
        vc=ze/x
        if (x.le.rc) vc=ze*(3-(x/rc)**2.)/2./rc
        vce=li*(li+1)/x/x
        poti(j)=(vn+vc)*rm+vce
  520 continue
c
c     computing the initial wave function (wfi) and phase shift (dep)
c
      call dephase(n,h,poti,wfi,eps,dep,li,eta,qk,rfin,ifail)
      if (ifail.eq.1) goto 501
      elmat=0.
      elmatin=0.
      elmatre=0.
      wff(0)=0.
      do 530 j=1,n
        if (lam.eq.1) then
          z=wff(j)*wfi(j)*(j*h)**1.
          elmat=elmat+z
        endif
        if (lam.eq.2) then
          zin=wff(j)*wfi(j)
          zre=wff(j)*wfi(j)*(j*h)**2.
          elmatin=elmatin+zin
          elmatre=elmatre+zre
        endif
        if (lam.eq.3) then
          z=wff(j)*wfi(j)
          elmat=elmat+z
        endif
c        write(8967,*) jspin,jparity,j,wff(j),wfi(j)
  530 continue
      if (lam.eq.1) then
        elmat=(elmat*h)**2.*fac2
        crossx(1)=crossx(1)+elmat
        xsall(1)=xsall(1)+elmat
        sfactor(1)=sfactor(1)+elmat*fac3
      endif
      if (lam.eq.2) then
        elmatin=(elmatin*h)**2.*fac2in
        elmatre=(elmatre*h)**2.*fac2re
        crossx(2)=crossx(2)+elmatin+elmatre
        xsall(2)=xsall(2)+elmatin+elmatre
        sfactor(2)=sfactor(2)+(elmatin+elmatre)*fac3
      endif
      if (lam.eq.3) then
        elmat=(elmat*h)**2.*fac2
        crossx(3)=crossx(3)+elmat
        xsall(3)=xsall(3)+elmat
        sfactor(3)=sfactor(3)+elmat*fac3
      endif
      esig=esig+elmat+elmatin+elmatre
      esigjp=esigjp+elmat+elmatin+elmatre
      esigex=esigex+elmat+elmatin+elmatre
  501 continue
      if (ipassf.eq.-1) goto 500
  500 continue
c
  220 continue
      xspop(jlev,jspin,jparity)=esigjp
  221 continue
c
      xspopex(jlev)=esigex
      numf=numf+1
  210 continue
c
c     end of loop all the available excited states in final nucleus
c
      return
      end
c
      subroutine sixj(j1,j2,j3,l1,l2,l3,q)
c
c     sixj
      implicit double precision (a-h,o-z)
      dimension m(7),m1(4),m2(4),m3(4),ft(0:100)
      save
      data lmem,(ft(i),i=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/
      entry sixji(j1,j2,j3,l1,l2,l3,q)
      i1=j1
      i2=j2
      i3=j3
      k1=l1
      k2=l2
      k3=l3
      is=0
      q=0
      m(1)=i1+i2+i3
      m(2)=i1+k2+k3
      m(3)=k1+i2+k3
      m(4)=k1+k2+i3
      do 17 i=1,4
      if(mod(m(i),2).eq.1) goto 8
   17 continue
      l=max(i1+i2+k1+k2,i1+i3+k1+k3,i2+i3+k2+k3)
      l=l/2
      if(l.le.lmem) goto 6
      do 10 i=lmem,l
   10 ft(i+1)=ft(i)*(i+1)
      lmem=l
    6 if(i1.lt.abs(i2-i3).or.i1.gt.i2+i3) return
      if(i1.lt.abs(k2-k3).or.i1.gt.k2+k3) return
      if(k1.lt.abs(i2-k3).or.k1.gt.i2+k3) return
      if(k1.lt.abs(k2-i3).or.k1.gt.k2+i3) return
      if(i1) 8,2,1
    2 if(i2.lt.0) goto 8
    9 if(i3.lt.0) goto 8
   14 if(k1.lt.0) goto 8
   19 if(k2.lt.0) goto 8
   23 if(k3.lt.0) goto 8
   27 q=sqrt(1.0d0/(i2+1)/(k2+1))
      is=(i2+k2+k1)/2+is
      if(mod(is,2).eq.1) q=-q
      return
    1 if(i1.gt.1) goto 3
      if(i2.lt.0) return
   12 if(i3.lt.0) return
   16 if(k1.lt.0) return
   21 if(k2.lt.0) return
   25 if(k3.lt.0) return
   28 if(i2.lt.i3) goto 4
      ic=i2
      i2=i3
      i3=ic
      ic=k2
      k2=k3
      k3=ic
    4 if(k2.gt.k3) goto 5
      i11=i1+k1+i2-k2
      i11=i11/2
      i12=i11-i2+k2
      q=sqrt(i11*i12*1.0d0/i3/(i3+1)/k3/(k3+1))
      is =i11+k2+is
      if(mod(is,2).eq.1) q=-q
      return
    5 i11=k3-k1+i2
      i11=i11/2+1
      i12=i11+k1+1
      q=sqrt(i11*i12*1.0d0/i3/(i3+1)/k2/(k2+1))
      is =i12-1+is
      if(mod(is ,2).eq.1) q=-q
      return
    3 if(i2.ge.i1) goto 7
      if(i2.lt.0) goto 8
      ic=i2
      i2=i1
      i1=ic
      ic=k1
      k1=k2
      k2=ic
      if(i1.eq.0) goto 9
      if(i1.eq.1) goto 12
    7 if(i3.ge.i1) goto 13
      if(i3.lt.0) goto 8
      ic=i3
      i3=i1
      i1=ic
      ic=k3
      k3=k1
      k1=ic
      if(i1.eq.0) goto 14
      if(i1.eq.1) goto 16
   13 if(k1.ge.i1) goto 18
      if(k1.lt.0) goto 8
      ic=k1
      k1=i1
      i1=ic
      ic=k2
      k2=i2
      i2=ic
      if(i1.eq.0) goto 19
      if(i1.eq.1) goto 21
   18 if(k2.ge.i1) goto 22
      if(k2.lt.0) goto 8
      ic=k2
      k2=i1
      i1=ic
      ic=k1
      k1=i2
      i2=ic
      if (i1.eq.0) goto 23
      if(i1.eq.1) goto 25
   22 if(k3.ge.i1) goto 26
      if(k3.lt.0) goto 8
      ic=k3
      k3=i1
      i1=ic
      ic=k1
      k1=i3
      i3=ic
      if(i1.eq.0) goto 27
      if(i1.eq.1) goto 28
   26 m1(4)=i3
      m1(1)=i3
      m1(3)=k3
      m1(2)=k3
      m2(2)=i1
      m2(1)=i1
      m2(4)=k1
      m2(3)=k1
      m3(3)=i2
      m3(1)=i2
      m3(4)=k2
      m3(2)=k2
      m(1)=i1+i2+i3
      m(2)=i1+k2+k3
      m(3)=k1+i2+k3
      m(4)=k1+k2+i3
      q1=1
      do 11 i=1,4
      m(i)=m(i)/2
   11 q1=ft(m(i)-m1(i))*ft(m(i)-m2(i))*ft(m(i)-m3(i))*q1/ft(m(i)+1)
      q1=sqrt(q1)
      m1(1)=i1+k1
      m1(2)=i2+k2
      m1(3)=i3+k3
      ic=m1(1)+m1(2)
      m(5)=ic/2
      ic=m1(2)+m1(3)
      m(6)=ic/2
      ic=m1(1)+m1(3)
      m(7)=ic/2
      maxz=min(m(5),m(6),m(7))
      minz=max(m(1),m(2),m(3),m(4))
      x=0
      do 15 i=minz,maxz
      q2=1
      do 20 j=1,7
      ij=i-m(j)
      if(j.gt.4) ij=-ij
   20 q2=q2*ft(ij)
      q2=ft(i+1)/q2
   15 x=-x+q2
      q=x*q1
      is=maxz+is
      if(mod(is,2).eq.1) q=-q
      return
    8 print 1010,j1,j2,j3,l1,l2,l3
 1010 format(10h erreur 6j,2(3x,3i3))
      return
      end
c
      subroutine clebs (l1,l2,l3,m1,m2,m3,q)
c
c     clebs
      implicit double precision (a-h,o-z)
      dimension ft(0:100)
      save
      datalmem,(ft(i),i=0,10)/9,2*1.0d0,2.0d0,6.0d0,24.0d0,
     1 120.0d0,720.0d0,5040.0d0,40320.0d0,362880.0d0,3628800.0d0/
      is=(l1+l2+m1-m2)/2
      k3=-m3
      q1=l3+1
      q=0
      i1=l1
      i2=l2
      i3=l3
      k1=m1
      k2=m2
      if(k1+k2+k3.ne.0)return
      l=i1+i2+i3
      if(mod(l,2).eq.1)goto8
      l=l/2
      if(l.le.lmem)goto6
      do10i=lmem,l
   10 ft(i+1)=(i+1)*ft(i)
      lmem=l
    6 j1=abs(k1)
      j2=abs(k2)
      j3=abs(k3)
      if(i1.lt.abs(i2-i3).or.i1.gt.i2+i3)return
      if(j1+j2.eq.0)goto11
      j1=i1-j1
      j2=i2-j2
      j3=i3-j3
      if(j1)8,2,1
    2 if(i1.ne.0)goto1
      if(j2.lt.0)goto8
   13 if(j3.lt.0)goto8
    4 q=sqrt(q1/(i2+1))
      is=is+(i2-k2)/2
      if(mod(is,2).eq.1)q=-q
      return
    1 if(j2.gt.j1)goto3
      if(j2.lt.0)goto8
      is=is+l
      j1=j2
      k1=k2
      k2=m1
      i1=i2
      i2=l1
      if(i1.eq.0)goto13
    3 if(j3.gt.j1)goto5
      if(j3.lt.0)goto8
      is=is+l
      j1=k3
      k3=k1
      k1=j1
      i3=i1
      i1=l3
      j1=j3
      if(i1.eq.0)goto4
    5 if(k1.ge.0)goto9
      k1=-k1
      k2=-k2
      k3=-k3
      is=is+l
    9 continue
      q1=q1*ft(l-i3)/ft(l-i1)/ft(l-i2)/ft(l+1)
      i2=(i2+k2)/2
      i3=(i3+k3)/2
      k2=i2-k2
      k3=i3-k3
      j1=j1/2
      i1=j1+k1
      j2=i3-k2
      j3=max(j2,0)
      is=is+i1+k2
      x=0
      do7i=j3,j1
    7 x=-x+ft(i1+i)*ft(i2+i3-i)/ft(j1-i)/ft(i3-i)/ft(i-j2)/ft(i)
      q=x*       sqrt(q1*ft(j1)*ft(k2)*ft(k3)*ft(i3)/ft(i1)/ft(i2))
      if(mod(is,2).eq.1)q=-q
      return
    8 q=0
      print1010,l1,l2,l3,m1,m2,m3
 1010 format(10h erreur 3j,2(3x,3i3))
      return
   11 if(mod(l,2).eq.1)return
      i1=l-i1
      i2=l-l2
      i3=l-l3
      q=sqrt(ft(i1)*ft(i2)*ft(i3)/ft(l+1)*q1)
      i1=i1/2
      i2=i2/2
      i3=i3/2
      l =l/2
      q=q*ft(l )/ft(i1)/ft(i2)/ft(i3)
      if(mod(l +is,2).eq.1)q=-q
      return
      end
c
      subroutine calmagelc (kzmag,kamag,dmagxy,qelcxytar,b2tar)
c     subroutine for calculating the magnetic dipole and electric quadrupole
c
      implicit double precision (a-h,o-z)
      real b2tar
      dimension lnmag(33),snmag(33),nncom(33)
      dimension lpmag(33),spmag(33),npcom(33)
      data lnmag/
     &0,1,1,2,2,0,3,3,1,1,4,4,2,2,5,0,5,3,6,3,1,1,6,4,7,4,2,2,0,7,5,8,5/
      data snmag/
     &0.5,1.5,0.5,2.5,1.5,0.5,3.5,2.5,1.5,0.5,4.5,3.5,2.5,1.5,5.5,0.5,
     &4.5,3.5,6.5,2.5,1.5,0.5,5.5,4.5,7.5,3.5,2.5,1.5,0.5,6.5,5.5,8.5,
     &4.5/
      data nncom/
     &  2,  6,  8, 14, 18, 20, 28, 34, 38, 40, 50, 58, 64, 68, 80, 82,
     & 92,100,114,120,124,126,138,148,164,172,178,182,184,198,210,228,
     &238/
      data lpmag/
     &0,1,1,2,2,0,3,3,1,4,1,4,5,2,2,0,5,6,3,3,1,1,6,7,4,4,2,2,7,0,8,5,5/
      data spmag/
     &0.5,1.5,0.5,2.5,1.5,0.5,3.5,2.5,1.5,4.5,0.5,3.5,5.5,2.5,1.5,0.5,
     &4.5,6.5,3.5,2.5,1.5,0.5,5.5,7.5,4.5,3.5,2.5,1.5,6.5,0.5,8.5,5.5,
     &4.5/
      data npcom/
     &  2,  6,  8, 14, 18, 20, 28, 34, 38, 48, 50, 58, 70, 76, 80, 82,
     & 92,106,114,120,124,126,138,154,164,172,178,182,196,198,216,228,
     &238/
c
c     begin magnetic dipole
c
      gp=5.5856
      gn=-3.8263
      knmag=kamag-kzmag
c
c     even Z and even N
c
cAK   if (mod(kzmag,2).eq.0.and.mod(knmag,2).eq.0) then
      dmagxy=0.0
      dmagxyn=0.0
      dmagxyp=0.0
cAK   endif
c
c     odd Z and even N
c
      if (mod(kzmag,2).ne.0.and.mod(knmag,2).eq.0) then
        do i=1,33
         if (kzmag.lt.npcom(i)) then
          if (spmag(i)-lpmag(i).eq.0.5) then
          dmagxy=(spmag(i)-0.5)+0.5*gp
          exit
          endif
          if (spmag(i)-lpmag(i).eq.-0.5) then
          dmagxy=spmag(i)/(spmag(i)+1.0)*(spmag(i)+1.5-0.5*gp)
          exit
          endif
         endif
        enddo
      endif
c
c     even Z and odd N
c
      if (mod(kzmag,2).eq.0.and.mod(knmag,2).ne.0) then
       do i=1,33
        if (knmag.lt.nncom(i)) then
         if (snmag(i)-lnmag(i).eq.0.5) then
         dmagxy=0.5*gn
         exit
         endif
         if (snmag(i)-lnmag(i).eq.-0.5) then
         dmagxy=snmag(i)/(snmag(i)+1.0)*(0.0-0.5*gn)
         exit
         endif
        endif
       enddo
      endif
c
c      odd Z and odd N
c
      if (mod(kzmag,2).ne.0.and.mod(knmag,2).ne.0) then
       do i=1,33
        if (knmag.lt.nncom(i)) then
         if (snmag(i)-lnmag(i).eq.0.5) then
         dmagxyn=0.5*gn
         exit
         endif
         if (snmag(i)-lnmag(i).eq.-0.5) then
         dmagxyn=snmag(i)/(snmag(i)+1.0)*(0.0-0.5*gn)
         exit
         endif
        endif
       enddo
c
       do i=1,33
        if (kzmag.lt.npcom(i)) then
         if (spmag(i)-lpmag(i).eq.0.5) then
         dmagxyp=(spmag(i)-0.5)+0.5*gp
         exit
         endif
         if (spmag(i)-lpmag(i).eq.-0.5) then
         dmagxyp=spmag(i)/(spmag(i)+1.0)*(spmag(i)+1.5-0.5*gp)
         exit
         endif
        endif
       enddo
      dmagxy=dmagxyp+dmagxyn
      endif
c
c     finish  magnetic dipole
c
c     begin electric quadrupole
c
      aa1=dble(kamag)
      qelcxytar=dble(b2tar)*aa1**(5./3.)/0.9174368
c
c     finish electric quadrupole
c
c
c     end of subroutine for calculating the magnetic dipole and electric quadrupole
c
c==========================================================================
      return
      end
c
      subroutine num1l(n,h,e,s2,u,s,no,eps)
c
c     num1l
c     version corrigee le 21 nov 72
c     integration de l"equation de schroedinger par la methode de numero
c     pour e negatif
c     recherche de l"energie propre par la methode de raphson-newton
      implicit double precision (a-h,o-z)
      dimension u(0:200),s(0:200)
      save
      data rap1,rap2/0,0/
      h12=h*h/12
c     controle des conditions asymptotiques
      if(e.gt.0) e=0
      dei=0
      epss=.1d-10
      if(u(n-1).gt.epss) goto 10
      dei=u(n-1)-epss
      do 8 k=1,n
    8 u(k)=u(k)-dei
   10 u(n)=u(n-1)
c     calcul du nombre d"etats lies par integration a energie nulle
      s(0)=1.d-10
      s(1)=1.d-10
      b0=0
      aa=h12*u(1)
      if (s2) 16,18,16
   16 b0=-s(1)*aa
   18 b1=s(1)*(1-aa)
      do 38 k=2,n
      b2=12*s(k-1)-10*b1-b0
      if (abs(b2).lt.1.d+10) goto 22
      b2=b2*1.d-20
      b1=b1*1.d-20
   22 aa=h12*u(k)
      s(k)=b2/(1-aa)
      b0=b1
   38 b1=b2
      do 42 k=5,n
      n0=k
      if(u(k).lt.0) goto 44
   42 continue
   44 nel=0
      do 52 k=n0,n
      if (s(k-1)*s(k)) 46,50,52
   46 nel=nel+2
      goto 52
   50 nel=nel+1
   52 continue
      nel=nel/2
      if(nel.gt.no) goto 64
      if(nel.eq.no) goto 60
   62 no=-1
      return
   60 rap1=s(n-1)/s(n)
      rap2=exp(h*sqrt(u(n-1)-e))
      if(rap1.lt.rap2) goto 62
c     calcul de emin et emax entre lesquelles se trouve l"energie propre
   64 umin=u(1)
      do 70 k=2,n
      if(u(k).lt.umin) umin=u(k)
   70 continue
      emin=umin
      emax=0
c     debut de la recherche de l"energie propre dans l"intervalle maximu
      te=emax-emin
c     rejet de l"energie d"essai e proposee si elle est a l"exterieur de
c     bornes (emin,emax)
      if((e.lt.emin).or.(e.gt.emax)) e=emin+te/2
      e1=emin
      e2=emax
      j=2
      i=1
      goto 102
c     reduction des bornes emin et emax
   90 emin=e1
      emax=e2
      te=emax-emin
      j=2
   98 i=1
  100 e=emin+te*i/j
  102 de=0
  104 e=e+de
      if(e.gt.0) goto 204
      s(n)=1.d-10
      n1=n-1
      expo=min(h*sqrt((u(n-1)+u(n))/2-e),80.d0)
      rap2=exp(expo)
      s(n1)=s(n)*rap2
      aa=h12*(u(n1)-e)
      b0=s(n)*(1-aa)
      b1=s(n1)*(1-aa)
      n1=n-2
      do 138 kaux=1,n1
      k=n1-kaux+1
      b2=12*s(k+1)-10*b1-b0
      aa=h12*(u(k)-e)
      s(k)=b2/(1-aa)
      b0=b1
      b1=b2
      if(u(k).lt.e) goto 140
  138 continue
  140 n1=k
c     normalisation de la fonction d"onde a s(n1)
      do 146 kaux=n1,n
      k=n-kaux+n1
  146 s(k)=s(k)/s(n1)
c     debut de l"integration vers l"exterieur jusqu"a n1
      s(1)=1.d-10
      b0=0
      aa=h12*(u(1)-e)
      if(s2) 156,158,156
  156 b0=-s(1)*aa
  158 b1=s(1)*(1-aa)
      do 170 k=2,n1
      b2=12*s(k-1)-10*b1-b0
      aa=h12*(u(k)-e)
      s(k)=b2/(1-aa)
      b0=b1
  170 b1=b2
c     normalisation de la fonction a s(n1)
      do 174 k=1,n1
  174 s(k)=s(k)/s(n1)
c     calcul de la correction d"energie
      som=0
      do 180 k=1,n
  180 som=som+s(k)*s(k)
      de=((-s(n1-1)+2-s(n1+1))/(h*h)+u(n1)-e)/som
      if(abs(de).gt.eps) goto 104
c     calcul du nombre de noeuds de l"etat propre trouve
      do 182 k=5,n
      if(u(k).lt.e) goto 184
  182 continue
  184 n0=k
      nel=0
      do 192 k=n0,n1
      if(s(k-1)*s(k)) 186,190,192
  186 nel=nel+2
      goto 192
  190 nel=nel+1
  192 continue
      nel=nel/2
c     l"etat propre trouve est-il le bon
      if(nel-no) 198,214,202
  198 if(e.gt.e1) e1=e
      goto 204
  202 if(e.lt.e2) e2=e
  204 i=i+2
      if (i.le.j)  goto 100
      j=2*j
      if(abs(e1-emin).gt.eps.or.abs(emax-e2).gt.eps) goto 90
      goto 98
c     normalisation de la fonction propre
  214 som=1/sqrt(som*h)
      do 218 k=1,n
  218 s(k)=s(k)*som
      e=e+dei
      return
c     debut formats
c2000 format(/,35x,56hl"etat demande n"est pas lie. retour de num1l avec
c    1 no=-1,/)
c     fin formats
      end
c
      subroutine dephase(maxv,h,w,y,eps,delta,l,eta,qk,rfin,ifail)
c
c     dephase
      implicit double precision (a-h,o-z)
      dimension fg(2),dfc(101),gc(101),dgc(101),fc(101),w(200),y(200)
      ifail=0
      dd=0
      pas=h*eps
      qq=qk*qk
      ll=l+1
      l1=l*ll
      r=2.0*eta*qk
      rmax=maxv*h
      do20 k=2,maxv
      rmax=rmax-h
      if(abs(w(maxv-k+1)-(l1/rmax+r)/rmax).lt.eps)goto20
      if(k.gt.3)goto21
      print2001
      ifail=1
      return
   20 continue
      rmax=h*(maxv-2)
   21 h2 = h**2
      h212 = h2/12.0
      aa=h212*(qq-w(1))
      y(1)=h**ll
      b0=0.0
      b1=y(1)*(1.0+aa)
      do22 k=2,maxv
      aa=h212*(qq-w(k))
      b2=12.0*y(k-1)-10.0*b1-b0
      y(k)=b2/(1.0+aa)
      b0=b1
      b1=b2
   22 continue
      r=rmax
      n1=int((rmax+1.0d-5)/h)+1
      call coufra(qk*r,eta,l,l,fc,dfc,gc,dgc)
      fg(1)=fc(ll)
      fg(2)=gc(ll)
      do23 n=n1,maxv
      r=r+h
      call coufra(qk*r,eta,l,l,fc,dfc,gc,dgc)
      c=fg(2)*y(n)-gc(ll)*y(n-1)
      d=fc(ll)*y(n-1)-fg(1)*y(n)
      z=fg(1)*gc(ll)-fg(2)*fc(ll)
      del=abs(z)/sqrt(c*c+d*d)
      if(abs(1-dd/del).lt.pas) goto 24
      dd = del
      fg(1)=fc(ll)
      fg(2)=gc(ll)
   23 continue
c      print2000,l,qk
      ifail=1
   24 rfin=r
      d=d/c
      if(y(min(n,200))*(fc(ll)+gc(ll)*d).lt.0.)del=-del
      delta=atan(d)
      do25 k=1,maxv
   25 y(k)=y(k)*del
      return
c2000 format(15x,29h no conv. in dephase() for l=,i3,3h k=,f10.6)
 2001 format(15x,' the potential has no coulomb asymptotic form')
      end
c
      subroutine coufra(rho,eta,minl,maxl,fc,fcp,gc,gcp)
c
c     coufra
      implicit double precision (a-h,o-z)
      double precision k,k1,k2,k3,k4,m1,m2,m3,m4
c     dimension fc(maxl),fcp(maxl),gc(maxl),gcp(maxl)
      dimension fc(101),fcp(101),gc(101),gcp(101)
      save
      data accur,step/1.d-7,100.0d0/
      pace = step
      acc = accur
      r = rho
      ktr = 1
      lmax = maxl
      lmin1 = minl+1
      xll1 = minl*lmin1
      eta2 = eta**2
      turn = eta+sqrt(eta2+xll1)
      if(r.lt.turn.and.abs(eta).ge.1.d-6) ktr = -1
      ktrp = ktr
      goto 2
    1 r = turn
      tf = f
      tfp = fp
      lmax = minl
      ktrp = 1
    2 etar = eta*r
      rho2=r*r
      pl = lmax+1
      pmx = pl+0.5d0
c     fraction continue pour fp(maxl)/f(maxl) ; xl=f ; xlprime=fp ********
      fp = eta/pl+pl/r
      dk = etar+etar
      del = 0
      d = 0
      f = 1
      k = (pl*pl-pl+etar)*(pl+pl-1)
      if(pl*pl+pl+etar.ne.0.) goto 3
      r = r*1.0000001d0
      goto 2
    3 h = (pl*pl+eta2)*(1-pl*pl)*rho2
      k = k+dk+pl*pl*6
      d = 1/(d*h+k)
      del = del*(d*k-1)
      if(pl.lt.pmx) del = -r*(pl*pl+eta2)*(pl+1)*d/pl
      pl = pl+1
      fp = fp+del
      if(d.lt.0) f = -f
      if(pl.gt.20000.0d0) goto 11
      if(abs(del/fp).ge.acc) goto 3
      fp = f*fp
      if(lmax.eq.minl) goto 5
      fc(lmax+1) = f
      fcp(lmax+1) = fp
c     recurrence arriere pour f et fp ; gc,gcp utilises pour stockage ***
      l = lmax
      do 4 lp=lmin1,lmax
      pl = l
      gc(l+1) = eta/pl+pl/r
      gcp(l+1) = sqrt(eta2+pl*pl)/pl
      fc(l) =(gc(l+1)*fc(l+1)+fcp(l+1))/gcp(l+1)
      fcp(l) = gc(l+1)*fc(l)-gcp(l+1)*fc(l+1)
    4 l = l-1
      f = fc(lmin1)
      fp = fcp(lmin1)
    5 if(ktrp.eq.-1) goto 1
c     meme calcul pour r = turn si rho.lt.turn
c     p + i.q calcule en minl , equation (32)
      p = 0.0
      q = r-eta
      pl = 0.0
      ar = -(eta2+xll1)
      ai = eta
      br = q+q
      bi = 2.0
      wi = eta+eta
      dr = br/(br*br+bi*bi)
      di = -bi/(br*br+bi*bi)
      dp = -(ar*di+ai*dr)
      dq = ar*dr-ai*di
    6 p = p+dp
      q = q+dq
      pl = pl+2.0
      ar = ar+pl
      ai = ai+wi
      bi = bi+2.0
      d = ar*dr-ai*di+br
      di = ai*dr+ar*di+bi
      t = 1.0/(d*d+di*di)
      dr = t*d
      di = -t*di
      h = br*dr-bi*di-1.0
      k = bi*dr+br*di
      t = dp*h-dq*k
      dq = dp*k+dq*h
      dp = t
      if(pl.gt.46000.0d0) goto 11
      if(abs(dp)+abs(dq).ge.(abs(p)+abs(q))*acc) goto 6
      p = p/r
      q = q/r
c     calcul de fp,g,gp, et normalisation de f en l = minl **************
      g = (fp-p*f)/q
      gp = p*g-q*f
      w = 1.0/sqrt(abs(fp*g-f*gp))
      g = w*g
      gp = w*gp
      if(ktr.eq.1) goto 8
      f = tf
      fp = tfp
      lmax = maxl
c     calcul de g(minl) et gp(minl) par integration runge-kutta a partir
c     voir fox et mayers(1968) pg 202
      if(rho.lt.0.2d0*turn) pace = 999.0d0
      r3=1.0d0/3.0d0
      h = (rho-turn)/(pace+1.0)
      h2 = h/2.0
      i2 = int(pace+0.001d0)
      etah = eta*h
      h2ll = h2*xll1
      s = (etah+h2ll/r)/r-h2
    7 rh2 = r+h2
      t = (etah+h2ll/rh2)/rh2-h2
      k1 = h2*gp
      m1 = s*g
      k2 = h2*(gp+m1)
      m2 = t*(g+k1)
      k3 = h*(gp+m2)
      m3 = t*(g+k2)
      m3 = m3+m3
      k4 = h2*(gp+m3)
      rh = r+h
      s = (etah+h2ll/rh)/rh-h2
      m4 = s*(g+k3)
      g = g+(k1+k2+k2+k3+k4)*r3
      gp = gp+(m1+m2+m2+m3+m4)*r3
      r = rh
      i2 = i2-1
      if(abs(gp).gt.1.d300) goto 11
      if(i2.ge.0) goto 7
      w = 1.0/(fp*g-f*gp)
c     recurrence avant a partir de gc(minl) et gcp(minl)
c     renormalisation de fc et fcp pour chaque valeur de l **************
    8 gc(lmin1) = g
      gcp(lmin1) = gp
      if(lmax.eq.minl) goto 10
      do 9 l=lmin1,lmax
      t = gc(l+1)
      gc(l+1) = (gc(l)*gc(l+1)-gcp(l))/gcp(l+1)
      gcp(l+1) = gc(l)*gcp(l+1)-gc(l+1)*t
      fc(l+1) = w*fc(l+1)
    9 fcp(l+1) = w*fcp(l+1)
      fc(lmin1) = w*fc(lmin1)
      fcp(lmin1) = w*fcp(lmin1)
      return
   10 fc(lmin1) = w*f
      fcp(lmin1) = w*fp
      return
   11 w = 0.0
      g = 0.0
      gp = 0.0
      goto 8
      end
c
      subroutine intu1(u,nr,du,volj,rms,rin,key)
c**************************************************************
c
c     moments m[n] calculation for u(r) (n=-1,1,2) in  [ 0. , (nr-1)*du
c
c                        (nr-1)*du
c           m[n] = 4*pi* integral [u(r)*r**(n+2)] dr
c                         0.                            (3 < nr < 5001)
c
c     by neuton-kortes method (k=3)
c
c      m[-1]=rin ,  m[1]=volj  ,  m[2]=rms
c
c                 n        key
c                -1        -1
c               -1,1        0
c                 1         1
c              -1,1,2       2
c                1,2        3
c
c**************************************************************
      implicit double precision(a-h,o-z)
      dimension u(nr)
      fact=0.375*du*12.5663706144d0
      volj=0.
      rms=0.
      rin=0.
      sum0=0.
      sum1=0.
      sum2=0.
      a0=0.
      a1=0.
      a2=0.
      r=0.
      nr2=nr-1
      nr1=(nr2/3)*3
      do 400 kr=3,nr1,3
      r1=r+du
      r2=r1+du
      r3=r2+du
      if(key.eq.1.or.key.eq.3) goto 2
      b0=u(kr-1)*r1+u(kr)*r2
      c0=u(kr+1)*r3
      sum0=sum0+a0+3.*b0+c0
      a0=c0
      if(key.eq.-1) goto 1
   2  rr1=r1*r1
      rr2=r2*r2
      rr3=r3*r3
      b1=u(kr-1)*rr1+u(kr)*rr2
      c1=u(kr+1)*rr3
      sum1=sum1+a1+3.*b1+c1
      a1=c1
      if(key.lt.2) goto 1
      rrr1=rr1*rr1
      rrr2=rr2*rr2
      rrr3=rr3*rr3
      b2=u(kr-1)*rrr1+u(kr)*rrr2
      c2=u(kr+1)*rrr3
      sum2=sum2+a2+3.*b2+c2
      a2=c2
   1  r=r3
  400 continue
      if(nr2.le.nr1) goto 3
      nr1=nr1+1
      c5=du*4.d0/3.d0
      do 4 kr=nr1,nr2
      if(key.ne.1.and.key.ne.3) sum0=sum0+c5*(u(kr)*(kr-1)+u(kr+1)*kr)
      if(key.ne.-1) sum1=sum1+c5*du*(u(kr)*(kr-1)**2+u(kr+1)*kr**2)
      if(key.gt.1) sum2=sum2+c5*du**3*(u(kr)*(kr-1)**4+u(kr+1)*kr**4)
    4 continue
    3 rin=fact*sum0
      volj=fact*sum1
      if(sum1.ne.0.) rms=sum2/sum1
      return
      end
c
      subroutine foldpot(e,rho,ka,kz,vfold,hu,volj,ipot)
      implicit double precision(a-h,o-z)
      parameter(nmx=1000,ndimax=5001)
ctest common/sp/dq,drint(2),nr(2),nq,nu,du,fro(ndimax,2),ro(ndimax,2),
ctest*am(4),u(ndimax),roe(ndimax,2),beta
      dimension drint(2),nr(2),fro(ndimax,2),ro(ndimax,2),
     *am(4),u(ndimax),roe(ndimax,2)
      common /const/ pi,sqrpi,pi4,pi32
      common /opt/ iopt(16)
      common /iso/ iso
      dimension rho(2,3)
      dimension vfold(nmx)
      dimension bcof(ndimax),bint(9),b(1),ue(ndimax)
      dimension a1(6),a2(6),b1(6),b2(6)
      equivalence (am(1),amp),(am(2),amt),(drint(1),d1),(drint(2),d2),
     *  (nr(1),n1),(nr(2),n2)
      data a1/1577.1071d0,34751.04d0,19839.156d0,25129.599d0,25129.599d0
     *,-15348.25d0/,b1/0.16,16.0,16.,16.,16.,16./,a2/-1239.9272d0,
     * -12754.865d0,-9857.0604d0,-10726.653d0,-10726.653d0,5908.71d0/,
     * b2/0.0625,6.2500,6.25,6.25,6.25,6.25/
      pi=3.1415926536d0
      sqrpi=dsqrt(pi)
      pi4=4.*3.1415926536d0
      pi32=sqrpi*pi
      zero=0.d0
      lone=1
      k3=3
      iopt(1)=ipot
      iopt(2)=1
      iopt(3)=0
      iopt(4)=0
      iopt(5)=3
      iopt(6)=0
      iopt(7)=0
      iopt(8)=0
      iopt(9)=0
      iopt(10)=0
      iopt(11)=0
      iopt(12)=2
      iopt(13)=1
      iopt(14)=0
      iopt(15)=0
      iopt(16)=0
c
      key=iopt(1)
      iopt2=iopt(2)
      alpha=0.
      beta=0.
      ce=1.
      if(key.eq.5) then
c
c     parameters of Kobos et al. (1984) NPA 425,205
c     ce=0.444
c     alpha=4.10
c     beta=10.67
c     parameters of Oberhummer et al. (Private communication, 1996)
c     fitted on G-matrix effective interaction of Jeukenne et al. 77
c
      ce=0.44073
      alpha=4.3529
      beta=10.639
      endif
      am(1)=1.
      am(2)=dble(ka)
      am(3)=0.
      am(4)=dble(kz)
      d1=0.01
      r1=10.
      d2=d1
      r2=r1
      dq=0.005
      rq=5.
      du=hu
      nu=nmx
      nu=((nu+1)/2)*2+1
      do 34 i=1,nu
      u(i)=0.d0
      ue(i)=0.d0
  34  continue
      akoff=0.5/(pi*pi)
      nr(1)=r1/d1
      nr(2)=r2/d2
      nr(1)=((nr(1)+1)/2)*2+1
      nr(2)=((nr(2)+1)/2)*2+1
      nq=rq/dq
      nq=((nq+1)/2)*2 +1
      if(iopt(12).ge.1) goto 30
      is1=1
      is2=1
      goto 32
  30  if(iopt(12).eq.2) goto 31
      is1=2
      is2=2
      goto 35
  31  is1=1
      is2=2
  35  anmzp=dabs(amp-2.*am(3))
      anmzt=dabs(amt-2.*am(4))
      if(anmzp.lt.0.00001d0.or.anmzt.lt.0.00001d0) is2=1
  32  do 33 iso=is1,is2
      if(iso.eq.2) key=6
      de=0.
      goto (11,12,13,14,15,16),key
   11 write(6,981)
      goto 20
   12 if(iopt(13).ne.0) de=-276.*(1.-0.005*e/amp)
      write(6,982)de
      goto 20
   13 if(iopt(13).ne.0) de=-81.
      write(6,983)de
      goto 20
   14 if(iopt(13).ne.0) de=-276.*(1.-0.005*e/amp)
c
      goto 20
   15 de=-276.*(1.-0.005*e/amp)
c
      goto 20
   16 if(iopt(13).ne.0) de=227*(1.-0.004*e/amp)
c
   20 continue
c
      no=2
      do 22 i=1,ndimax
      ro(i,no)=0.
   22 roe(i,no)=0.
      tiso=1.
      if(iso.gt.1) tiso=-1.
      x=drint(no)
      do 21 i=2,nr(no)
      xx=x*x
      rhon=rho(1,1)/(1.+dexp((x-rho(1,2))/rho(1,3)))
      rhop=rho(2,1)/(1.+dexp((x-rho(2,2))/rho(2,3)))
      ro(i,no)=rhop*tiso+rhon
      if(iopt(1).eq.5) roe(i,no)=pi4*xx*ro(i,no)*dexp(-ro(i,no)*beta)
      ro(i,no)=ro(i,no)*xx*pi4
      x=x+drint(no)
   21 continue
c
      nintV=max0(nr(1),nr(2))
      if(amp.eq.1.) nintV=nr(2)
      if(iopt(6).eq.1) write(6,993)
      iflg=0
      if(key.ne.5.and.iopt(4).eq.2.and.iopt(5).eq.2) iflg=1
      if(key.ne.5.and.amp.eq.1..and.iopt(5).eq.2) iflg=1
      q=0.
      kf=1
      do 200 k=1,nq
      qq=q*q
      sum1=0.
      sum2=0.
      s1=0.
      s2=0.
      if(d1.eq.d2.and.iflg.eq.0)
     *call ftrans(bcof,lone,nintV,zero,d1,q,bint,b,ndimax)
      if(amp.eq.1) goto 120
      if(d1.ne.d2.and.iflg.eq.0)
     *call ftrans(bcof,lone,n1,zero,d1,q,bint,b,ndimax)
      if(key.ne.5) goto 101
      do 100 j=2,n1
  100 s1=s1+bcof(j)*roe(j,1)
  101 if(iopt(4).ge.6) goto 110
      if(iopt(4).eq.2.and.iopt(10).eq.0) goto 110
      do 102 j=2,n1
  102 sum1=sum1+bcof(j)*ro(j,1)
      goto 126
  110 sum1=fro(k,1)
      goto 126
  120 s1=1.
      sum1=1.
  126 if(d1.ne.d2.and.iflg.eq.0)
     *call ftrans(bcof,lone,n2,zero,d2,q,bint,b,ndimax)
      if(key.ne.5) goto 131
      do 130 j=2,n2
  130 s2=s2+bcof(j)*roe(j,2)
  131 if(iopt(5).ge.6) goto 133
      if(iopt(5).eq.2.and.iopt(11).eq.0) goto 133
      do 132 j=2,n2
  132 sum2=sum2+bcof(j)*ro(j,2)
      goto 140
  133 sum2=fro(k,2)
  140 if(key.gt.1) goto 172
      c=a1(key)*dexp(-qq*b1(key))+a2(key)*dexp(-qq*b2(key))+de
      goto 180
  172 c=a1(key)/(qq+b1(key))+a2(key)/(qq+b2(key))+de
  180 if(iopt(7).eq.1) goto 191
      fro(k,1)=(sum1*sum2+s1*s2*alpha)*c*qq
  191 if(iopt(6).eq.0) goto 200
      if(k.ne.kf) goto 200
      write(6,920) q,c,sum1,sum2,s1,s2
      kf=kf+5
  200 q=q+dq
      if(iopt(7).eq.1) goto 1
  920 format(f7.3,5d14.6)
      r=0.
      do 300 i=1,nu
      sum=0.
      call ftrans(bcof,lone,nq,zero,dq,r,bint,b,ndimax)
      do 400 k=1,nq
  400 sum=sum+bcof(k)*fro(k,1)
      u(i)=akoff*sum*ce+u(i)
  300 r=r+du
  33  continue
c
      r=du
      do 520 i=2,nu
      ud=u(i)
      u(i)=ud+ue(i)
      if (i-1.le.nmx) vfold(i-1)=u(i)
  520 r=r+du
      call intu1(u,nu,du,volj,rms,vmi33,k3)
      volj=volj/(amp*amt)
      rms=dsqrt(rms)
      write(6,911) volj,rms
      if(iopt2.eq.0) goto 750
c
 2097 format(/,'   r    ','     v(r)')
      r=du
      do 2098 i=2,nu
 2099 format(f8.3,f12.4)
 2098 r=r+du
  750 continue
  555 continue
    1 continue
  800 format(i4,f8.4,8a8)
  950 format(5d16.8)
 1990 format(e11.4)
 1991 format(8e11.4)
 2991 format(5d12.4)
 1999 format(f8.2,2x,f6.3,8a8)
  900 format(///'  folding model potential',/,'   pc version of ii-94')
  901 format(16i1,8a8)
  803 format(/'      dr = ',f7.4,'      nr = ',i4,
     1 '    dq = ',f7.4,'    nq = ',i4,
     2 '    du = ',f7.4,'   nu=',i4)
  804 format(/'      dr1 = ',f7.4,'     nr1 = ',i4,'     dr2 = ',f7.4,
     1 '     nr2 = ',i4,'     dq = ',f7.4,'     nq = ',i4,
     2 '     du = ',f7.4,'   nu=',i4)
  805 format(1h0,/,' option  1234567890123456',/9x,16i1,8a8)
  801 format(4d15.8)
  802 format(/'      dr = ',f7.4,10x,8a8,/(5d15.7))
  902 format(9f8.4)
  907 format(i4,8f8.4)
c  97 format(//'  folding potential',//,
c    *4x,'r',9x,'u',10x,'ud',10x,'ue',/)
c  98 format(f8.3,10d12.4)
  911 format(/'      volume integral = ',f7.1,' mev*fm**3      <rms>  =
     1',f7.3,' fm')
  981 format(/'      n - n  forces :  gauss    v(r) = - 553.18 * exp(-r*
     1r/0.64) + 1781.4 * exp(-r*r/0.25) ')
  993 format(/'      fourie transforms'/5x,'q',5x,'v(k)',11x,
     * 'ro1(k)',8x,'ro2(k)',8x,'ro1"(k)',7x,'ro2"(k)')
  936 format(13f6.3)
  982 format(/'      n - n  forces :  m3y   v(r) = 11061.6 * exp(-4r)/4r
     1 - 2537.5 * exp(-2.5r)/2.5r ',f7.1,'*del(r)')
  983 format(/'      n - n  forces :  m3y      v(r) = 6315 * exp(-4r)/4r
     1 - 1961 * exp(-2.5r)/2.5r ',f7.1,'*del(r)')
  984 format(/'      n - n  forces :  m3y      v(r) = 7999 * exp(-4r)/4r
     1 - 2134 * exp(-2.5r)/2.5r ',f7.1,'*del(r)')
  985 format(/
     *'      n - n  forces :'/'       ddm3y  v(r) = [ 7999 * exp(-4r)/4r
     1 - 2134 * exp(-2.5r)/2.5r ',f7.1,'*del(r)] * ce[1+a*exp(-beta*ro(r
     2))]')
  986 format(/'  n - n  forces :  m3y-isovec.  v(r) =-4886 * exp(-4r)/4r
     1 + 1176 * exp(-2.5r)/2.5r ',f7.1,'*del(r)')
  988 format('      ce = ',f7.4,'      alpha = ',f7.4,'      beta = ',
     *  f7.4)
  989 format(/'        exact snke potential')
  903 format(/'      e = ',f9.2,'    ap,zp = ',2f5.2,'  at,zt= ',2f5.1)
  990 format(/'      ro',i1,' = ',f8.5,' / [ 1 + exp( ( r - ',f7.4,
     1 ' )/',f7.4,') ]')
      return
      end
