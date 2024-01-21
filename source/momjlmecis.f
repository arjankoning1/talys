      subroutine momjlmecis(zi,ni,zpi,einci,ilv1,ilw1,
     &    itjlmr,itjlmi,rhomomn,rhomomp,vpot,rcjlm)
c     version as of 29/03/07.
c
c     ******************************************************************
c
c     get working parameters from TALYS(ecisinput) and sets the default
c     values.
c     Folds the density held in file namegd with the modified JLM
c     interaction of PRC 63 024607 2001,
C     and outputs the potentials in ECIS format
c
c     initial values for many of the parameters are established in
c     the data statements below and in the momblockdata common block.
c
c     ******************************************************************
c
      implicit none
      include 'mom.cmb'
      integer ni,zi
      real*8 zpi,einci,ilv1,ilw1,itjlmr,itjlmi,rc
      integer j
      real rhomomn(numjlm,6),rhomomp(numjlm,6),vpot(numjlm,6),rcjlm
      real*8 einc0,einc1,rden0
      integer i
      data einc0,einc1/0.d0,0.d0/
      data rden0/0.d0/
c
      if(zpi.eq.0.d0) then
c     neutron
         emp=1.0087d0
         zp=0.0d0
      else
c     proton
         emp=1.0073d0
         zp=1.0d0
      endif
      z      = zi*1.d0
      a      = (ni+zi)*1.d0
      einc0  = einci
      irel   = 1
      al0r   = 1.0d0
      altr   = ilv1
      al0i   = 1.d0
      alti   = ilw1
      tjlmr  = itjlmr
      tjlmi  = itjlmi
      icorv  = 1
      icorw  = 1
      einc=einc0+einc1
      kmw    = 1
      rmax=npts*h
      em=a
      rden=rden0
      if(rden.eq.0.d0) rden=20.d0
      lamax=ilamaxi
c     calculate the coulomb radius based on Elton's formula
      rc=1.123+2.352*(a**(-.666666))-2.07*(a**(-1.333333))
c     do preliminary calculations
      call prelim
      call prepgl
c     Initialize the densities
      do i=1,numjlm
        do j=1,5
          rhoni(i,j)=rhomomn(i,j)
          rhopi(i,j)=rhomomp(i,j)
        enddo
        rhop(i)=rhopi(i,1)
        rhon(i)=rhoni(i,1)
      enddo
c     calculate coulomb field for coulomb corrections if necessary
      if (zp.ne.0.d0) call coul
c     do the folding for the central potential
      call jlmcen
c     do the calculation of the S.O. potential form factor
      call simpls
c
c    Put results in arrays to be exported
c
      do i=1,numjlm
        vpot(i,1)=vcr(i)
        vpot(i,2)=vci(i)
        vpot(i,3)=0.
        vpot(i,4)=0.
        vpot(i,5)=vlsr(i)
        vpot(i,6)=vlsi(i)
      enddo
      rcjlm=rc
      return
      end
      subroutine prelim
c     version as of 26/3/2007.
c
c     ******************************************************************
c     do preliminary calculations for kinematics & coulomb
c     ******************************************************************
c
      implicit none

      include 'mom.cmb'
      real*8 aamu,aksq,grel,xemp,xem,xemt
      real*8 aux,wcm,temp,x,wm,wmp
      real*8 amen
      data amen/931.47d0/
c
c     constants used below are ---
c
c     .047846 = 2 * ( 1 aamu ) / (hbar) ** 2
c     .034447 = (1 aamu) * e ** 2 / (hbar) ** 2
c
      aamu=emp*em/(emp+em)
      emu=.047846d0*aamu
c
c     choose between relativistic and nonrelativistic kinematics.
c
      if(irel.ne.0) goto 5
c
c     nonrelativistic kinematics.
c
      ecm=einc*em/(emp+em)
      aksq = emu * ecm
      ak = sqrt (aksq)
      grel=1.d0
      goto 22
c
c     relativistic kinematics.  schrodinger equation is written with
c     reduced mass and ak.  ecm does not appear.  the choice for the
c     relativistic correction concerns the choice made for the
c     reduced mass and ak.  the factor grel is the ratio of the
c     relativistic reduced mass to the nonrelativistic reduced mass.
c
c     model 1 ---
c
c     ak is chosen as ingemarsson, physica scripta vol. 9, p. 156
c     (1974), expression (17).  the reduced mass is formed from
c     the relativistic masses of target and projectile in the c.m.
c     system.
c
    5 xemp=amen*emp
      xem =amen*em
      xemt=xemp+xem
      aux=2.d0*xem*einc
      wcm=sqrt(xemt*xemt+aux)
      temp=einc*(einc+2.d0*xemp)
      ak=sqrt(temp)*xem/(197.32d0*wcm)
      aksq=ak*ak
      x=aux/(xemt*xemt)
      if(x.gt.1.e-3) goto 8
      ecm=(.0625d0*x*x*x-.125d0*x*x+.5d0*x)*xemt
      goto 10
    8 ecm=wcm-xemt
c
   10 wm=xem*(einc+xemt)/wcm
      wmp=wcm-wm
      grel=wmp*wm/((wmp+wm)*amen*aamu)
c
      aamu=grel*aamu
      emu=grel*emu
c
   22 zzp=.034447d0*aamu*z*zp
c     eta = zzp/ak
      return
      end
      subroutine prepgl
c     ********************************************8
c     calculation of the Gauss-Legendre coeff.
c     calculation of  Ylm on the  GL mesh
c     for spherical npgl=1  -----> trivial
c     *********************************************
      implicit none
      include 'mom.cmb'
      real*8 sqr,sqr1
      real*8 cosi,cosi2,cosi4,cosi6,cosi8
      real*8 xpl5,xpl6,xpl7,xpl8,xpl9,xpl10
      integer ita,ii,jj
      npgl=1
      if(npgl.eq.1)then
         xg(1)=0.0d0
         wg(1)=1.0d0
      endif
      if(npgl.eq.16)then
         XG(1)=0.04830766568773800d0
         XG(2)=0.14447196158279600d0
         XG(3)=0.23928736225213700d0
         XG(4)=0.33186860228212700d0
         XG(5)=0.42135127613063500d0
         XG(6)=0.50689990893222900d0
         XG(7)=0.58771575724076200d0
         XG(8)=0.66304426693021500d0
         XG(9)=0.73218211874028900d0
         XG(10)=0.79448379596794200d0
         XG(11)=0.84936761373256900d0
         XG(12)=0.89632115576605200d0
         XG(13)=0.93490607593773900d0
         XG(14)=0.96476225558750600d0
         XG(15)=0.98561151154526800d0
         XG(16)=0.99726386184948100d0
         WG(1)=0.09654008851472700d0
         WG(2)=0.09563872007927400d0
         WG(3)=0.09384439908080400d0
         WG(4)=0.09117387869576300d0
         WG(5)=0.08765209300440300d0
         WG(6)=0.08331192422694600d0
         WG(7)=0.07819389578707000d0
         WG(8)=0.07234579410884800d0
         WG(9)=0.06582222277636100d0
         WG(10)=0.05868409347853500d0
         WG(11)=0.05099805926237600d0
         WG(12)=0.04283589802222600d0
         WG(13)=0.03427386291302100d0
         WG(14)=0.02539206530926200d0
         WG(15)=0.01627439473090500d0
         WG(16)=0.00701861000947000d0
      endif
      SQR=DSQRT(4.25d0/PI)
      sqr1=dsqrt(21.d0/(4.d0*pi))
      DO 300 ita=1,npgl
      COSI=XG(ita)
      COSI2=COSI*COSI
      COSI4=COSI2*COSI2
      COSI6=COSI2*COSI4
      COSI8=COSI2*COSI6
      PLL(ita,1)=1.0d0
      PLL(ita,2)=0.63078310d0*(1.5d0*COSI2-0.50d0)
      PLL(ita,3)=0.846284367d0*(4.375d0*COSI4-3.75d0*COSI2+0.375d0)
      PLL(ita,4)=1.017107d0*(14.4375d0*COSI6
     &  -19.6875d0*COSI4+6.5625d0*COSI2
     &   -0.3125d0)
      PLL(ita,5)=SQR*(50.2734375d0*COSI8-93.84375d0*COSI6
     &   +54.140625d0*COSI4
     &   -9.84375d0*COSI2+0.2734375d0)
      xpl5=cosi*(63.*cosi4-70.*cosi2+15.)/8.d0
      xpl6=(231.*cosi6-315.*cosi4+105.*cosi2-5.)/16.d0
      xpl7= (13d0*cosi*xpl6-6d0*xpl5)/7.d0
      xpl8= (15d0*cosi*xpl7-7d0*xpl6)/8.d0
      xpl9= (17d0*cosi*xpl8-8d0*xpl7)/9.d0
      xpl10=(19d0*cosi*xpl9-9d0*xpl8)/10.d0
      PLL(ita,6)=sqr1*xpl10
  300 CONTINUE
c     normalisation des Ylm
      do jj=2,6
        do ii=1,npgl
           pll(ii,jj)=pll(ii,jj)*dsqrt(4.0*pi)
        enddo
      enddo
      return
      end
      subroutine coul
c     version as of 3/19/88.
c
c     ******************************************************************
c
c     computes coulomb potential by integration.
c     potential is computed from density rhop.
c
c
c     ******************************************************************
c
      implicit none
      include 'mom.cmb'
      real*8 q(numjlm),r,aux,xnorm
      integer i,j
c
c     first find the charge (unnormalized) contained within a
c     given radius.  outward integration is used when density rhop
c     is employed.
c
c
      do i=1,npts
         q(i)=0.0
         vcoul(i)=0.0
      enddo

      q(1)=rhop(1)*h*h*h/3.d0
      r=2.d0*h
      do   i=2,npts
         aux=r*r*rhop(i)+(r-h)*(r-h)*rhop(i-1)
         q(i)=q(i-1)+aux*h/2.d0
         r=r+h
      enddo
c
c     now find the unnormalized coulomb potential by inward inte-
c     gration.
c
      vcoul(npts)=q(npts)/rmax
      r=rmax-h
      do  j=2,npts
         i=npts-j+1
         aux=q(i)/(r*r)+q(i+1)/((r+h)*(r+h))
         vcoul(i)=vcoul(i+1)+aux*h/2.d0
         r=r-h
      enddo
c
c     normalize to correct coulomb potential at outer boundary.
c
      xnorm=(2.d0*zzp)/(rmax*emu*vcoul(npts))
c
      do i=1,npts
         vcoul(i)=xnorm*vcoul(i)
      enddo
      return
      end
      subroutine jlmcen
c     version as of 3/19/88.
c
c     ******************************************************************
c
c     calculates real and imaginary parts of central potential
c     using jlm prescription with improved lda
c     ref: phys. rev. c16(1977)80
c     uses the improved lda evaluated  at the projectile position (r')
c     as in PRC 58, 1118 (1998) and PRC 63, 024607 (2001)
c
c     if kmw is nonzero, imaginary potentials are multiplied by
c     a k-mass factor, which is a correction to the original jlm work.
c     see ref 41,42 of PRC 58, 1118 (1998)
c
c     ******************************************************************
c
      implicit none

      include 'mom.cmb'
      real*8 temp(32)
      real*8 str(numjlm,32,4),st(4)
      real*8 dc(4,4),fc(4,4)
      real*8 fro(numjlm),fio(numjlm),gro(numjlm),gio(numjlm)
      real*8 capd,r,cen,vzero,ren,redms,akmass,wzero,rho,alf
      real*8 aaf,eef,xfr1,xfr2,eps1,eps2,eps,e,vcorr,f,emass,sdot
      real*8 tau3
      integer i,j,ita,n,ii,jj,ipol
      integer ix,ity,nblk8
c
c
c     input arrays of coefficients in block data
      tau3=2.d0*(.5d0-zp)
      call xform
      do  i=1,numjlm
         vci(i)=0.d0
         vcr(i)=0.d0
      enddo
c     typo in PRC 58 1118 (1998) v1.01
      capd=126.25d0
      do  i=1,4
         do  j=1,4
            dc(i,j)=dhc(i,j)
            fc(i,j)=fhc(i,j)
         enddo
      enddo
c     loop on the angular coordinate theta (ita)
c     for this spherical version npgl=1
c     this loop is here only for historical reasons since the
c     present spherical code is a version of our deformed code
c     that is restricted to the spherical symmetry (l=0)
      do 132 ita=1,npgl
      r=h
c     loop on the radial coordinate r
      do 32 n=1,npts
      cen=0.d0
      vzero=0.d0
      ren=0.d0
      redms=1.d0
      akmass=1.d0
      wzero=0.d0
      rho=vang(n,ita)
      alf=valf(n,ita)
c
c      calculate the Fermi energy eps using the modified expression
c     of PRC 58,1118 (1998). bug in v1.0, corrected in v1.01
      aaf=2.d0
      eef=9.d0
      xfr1=1./(1.+exp((ecm-eef)/aaf))
      xfr2=1.-xfr1
      eps1=-22.d0-rho*(298.52d0-3760.23d0*rho+12435.82d0*rho**2)
      eps2=rho*(-510.8d0+3222.d0*rho-6250.d0*rho**2)
      eps=eps1*xfr1+eps2*xfr2
      e=ecm-vcorr(r)
      f=e
      if(icorv.eq.0)e=ecm
      if(icorw.eq.0)f=ecm
      do  ii=1,3
         do  jj=1,3
            vzero=vzero+ac(ii,jj)*(rho**ii)*(e**(jj-1))
            redms=redms-(jj-1)*ac(ii,jj)*(rho**ii)*(e**(jj-2))
            ren=ren+bc(ii,jj)*(rho**ii)*(e**(jj-1))
            akmass=akmass-cc(ii,jj)*(rho**ii)*(e**(jj-1))
         enddo
      enddo
      emass=redms/akmass
      do i=1,4
         do j=1,4
            wzero=wzero
     &            +dc(i,j)*(rho**i)*(f**(j-1))/(1.d0+capd/((f-eps)**2))
            cen=cen+fc(i,j)*(rho**i)*(f**(j-1))/(1.d0+1.d0/(f-eps))
         enddo
      enddo
      fr(n)=al0r*xlv*vzero
      fi(n)=al0i*xlw*wzero
      gr(n)=altr*xlv*akmass*ren
      gi(n)=alti*xlw*cen/emass
      if(kmw.eq.0) goto 30
      fi(n)=akmass*fi(n)
      gi(n)=akmass*gi(n)
 30   continue
 32   r=r+h
      do n=1,npts
         rho=vang(n,ita)
         alf=valf(n,ita)

         str(n,ita,1)=fr(n)
         str(n,ita,3)=fi(n)
         str(n,ita,2)=alf*gr(n)
         str(n,ita,4)=alf*gi(n)
      enddo
 132  continue
c     pot diagonal l=0
      do ix=1,npts
         do ity=1,4
            do ita=1,npgl
               temp(ita)=str(ix,ita,ity)
            enddo
            st(ity)=sdot(npgl,temp,1,wg,1)
         enddo
         vcentr(ix,1)=st(1)+tau3*st(2)
         vcenti(ix,1)=st(3)+tau3*st(4)
      enddo
      do ipol=1,1
         do i=1,npts
               fr(i)=vcentr(i,1)
               gr(i)=vcentr(i,1)
               fi(i)=vcenti(i,1)
               gi(i)=vcenti(i,1)
         enddo
         nblk8 = npts/8
         call convolbesl0(fr,gr,fro,gro,nblk8,tjlmr)
         call convolbesl0(fi,gi,fio,gio,nblk8,tjlmi)
         do i=1,npts
               vcentr(i,1)=fro(i)
               vcenti(i,1)=fio(i)
               vcr(i)=vcentr(i,1)
               vci(i)=vcenti(i,1)
         enddo
      enddo
      return
      end
      subroutine xform
c     26/03/07
c     **************************************************************
c     calculate the density in (r,theta) coordinate from the (r,l) ones
c     trivial for l=0
c     really useful in deformed densities
c     included in this program only for code base compatibility with
c     full-fledged code
c     **************************************************************
      implicit none
      include 'mom.cmb'
      real*8 alf,rho,anorm
      integer lamax1,ix,ita,ipol,i,j
      lamax1 = lamax
      DO ita=1,32
         DO IX=1,numjlm
            valf(ix,ita)=0.0d0
            VANG(IX,ita)=0.0d0
         enddo
      enddo
      anorm = 1.0d0
      do ita=1,npgl
        do ix=1,npts
          do ipol=1,lamax1
             rho=rhopi(ix,ipol)+rhoni(ix,ipol)
             alf=rhoni(ix,ipol)-rhopi(ix,ipol)
             vang(ix,ita)=vang(ix,ita)+rho*pll(ita,ipol)*anorm
             valf(ix,ita)=valf(ix,ita)+alf*pll(ita,ipol)*anorm
          enddo
        enddo
      enddo
      do i=1,npts
         do j=1,npgl
            valf(i,j)=valf(i,j)/vang(i,j)
         enddo
      enddo
      return
      end
      subroutine convolbesl0(frc,grc,ofr,ogr,nblk8,t)
c     ****************************************************************
c
c     does the convolution using the bessel function expansion of exp.
c     restricted to l=0 : trivial
c     ****************************************************************
      implicit none
      include 'mom.cmb'
      real*8 frc(numjlm), grc(numjlm), ofr(numjlm), ogr(numjlm),t
      integer nblk8
      real*8 mu,mu1,ma,newtfac
      integer iz(numjlm),ifac(12),ibad(12),npo,il,i,jj,indic,j,k
      integer ir,ii,npx,incx
      real*8 cho(11), tho(11),big(numjlm),bag(8)
      real*8 zz(numjlm),bes(numjlm,numjlm),dup(12),cxp(numjlm)
      real*8 axp(numjlm),bxp(numjlm),r(numjlm),rprim(numjlm),r2(numjlm)
      real*8 rprim2(numjlm),theta(numjlm),cprim(numjlm)
      real*8 dx,quatpi,pimu,sqpimu,fnorm,rj,ra,rfac,rbad,term
      real*8 expo,ta,tb,ca,cb,signe,soms,somv
      real*8 sdot,ssum
c   coeff. integ. open newton-cotes 8 pts see Abramowitz & Stegun, Dover
      data bag/0.,460.,-954.,2196.,-2459.,2196.,-954.,460./
      npo=8*nblk8
      il=0
      quatpi=4.0d0*pi
      dx= 0.1
      mu=t**2
      mu1=1.d0/mu
      pimu=mu*pi
      sqpimu=dsqrt(pimu)
      fnorm=quatpi/(sqpimu**3)
      newtfac=8.d0*dx/945.d0
c     preparation of the arrays: r, r', etc...
      do i=1,numjlm
         iz(i)=i
c        indic=mod(i,8)+1
         indic=i-8*(i/8)+1
         big(min(i+1,numjlm))=bag(indic)
         r(i)=(real(i)-.1d0)*dx
         rprim(i)=(real(i)-0.0d0)*h
         r2(i)=mu1*r(i)**2
         rprim2(i)=mu1*rprim(i)**2
         axp(i)=r2(i)
         theta(i)=1.
      enddo
      do i=120,numjlm
         theta(i)=dexp(-120.d0)
         axp(i)=r2(i)-120.d0
         rprim2(i)=rprim2(i)-120.d0
      enddo
      do i=1,numjlm
         cxp(i)=dexp(-axp(i))
         cprim(i)=dexp(-rprim2(i))
      enddo
c     coeff of the bessel function expansion, trivial for l=0
      do i=1,8
         cho(i)=0.0d0
         tho(i)=0.0d0
      enddo
      if (il.eq.0) then
           cho(1)=1.0d0
           tho(1)=-1.0d0
      else
         write(6,*) 'Only for spherical nuclei'
         stop
      endif

      ibad(1)=2*il+3
      ifac(1)=1
      jj=1
      do i=2,5
         ifac(i)=ifac(i-1)*i
         ibad(i)=ibad(i-1)*(2*il+2*i+1)
      enddo
      do i=1,il
         jj=jj*(2*iz(i)+1)
      enddo
      rj=real(jj)
c     calculate the bessel fct.
      axp(1)=0.
      bxp(1)=0.
      do j=1,npts
         do i=2,npo
           zz(i)=2*r(i)*rprim(j)*mu1
           if(zz(i).le.0.1d0) then
              ra=1.d0
              signe=-1.d0
              do ir=1,3
                 rfac=real(ifac(ir))
                 rbad=real(ibad(ir))
                 term=(-0.5d0*zz(i)**2)**ir/(rfac*rbad)
                 if(dabs(term).lt.1.d-70)term=0.d0
                 ra=ra+term*signe
                 signe=-signe
              enddo
              bes(i,j)=ra*(zz(i)**il)/rj
           else
              expo=0.
              if(zz(i).ge.170.d0)then
                 expo=-150.
                 zz(i)=zz(i)-150.d0
              endif
              axp(i)=dexp(zz(i))
              bxp(i)=1./axp(i)
              if(expo.eq.-150.)zz(i)=zz(i)+150.d0
              do ii=1,il+1
                 dup(ii)=zz(i)**(-iz(ii))
              enddo
              ta=sdot(il+1,cho(1),1,dup(1),1)
              tb=sdot(il+1,tho(1),1,dup(1),1)
              ca=ta*axp(i)*dexp(-expo)
              cb=tb*bxp(i)*dexp(expo)
              bes(i,j)=0.5d0*(ca+cb)
              if(dabs(bes(i,j)).lt.1.d-70) bes(i,j)=0.d0
           endif
         enddo
      enddo
c     do the folding
      npx=npo-1
      do j=1,npts
         ma=theta(j)*mu
         do k=2,npo
            axp(k)=bes(k,j)*r2(k)*big(k)*cxp(k)*theta(k)*frc(k)
            bxp(k)=bes(k,j)*r2(k)*big(k)*cxp(k)*theta(k)*grc(k)
         enddo
         incx=1
         soms=ssum(npx,axp(2),incx)
         somv=ssum(npx,bxp(2),incx)
         ofr(j)=soms*newtfac*fnorm*cprim(j)*ma
         ogr(j)=somv*newtfac*fnorm*cprim(j)*ma
      enddo
      return
      end
c     ******************************************************************
      subroutine simpls
c     26/10/95 Ebauge
c     corrigee 10/11/95
c     version 26/03/07
c     ******************************************************************
c     calculate the  SO form factor  Scheerbaum-style
c      Nucl. Phys. A257,77,(1986)
c     ******************************************************************
c
      implicit none
      include 'mom.cmb'
      real*8 vlt(numjlm),rr,dvdr,vsot
      integer i
      do  i=1,numjlm
         vlsr(i)=0.d0
         vlsi(i)=0.d0
      enddo
      do  i=1,npts
         if(zp.eq.0) vsot=1.d0/3.d0*(2.d0*rhop(i)+rhon(i))
         if(zp.eq.1) vsot=1.d0/3.d0*(2.d0*rhon(i)+rhop(i))
         vlt(i)=vsot
      enddo
      rr=npts*h
      dvdr=(vlt(npts)-vlt(npts-1))/h
      vlsr(npts)=dvdr/rr
      vlsi(npts)=vlsr(npts)
      rr=h
      dvdr=(vlt(2)-vlt(1))/h
      vlsr(1)=dvdr/rr
      vlsi(1)=vlsr(1)
      do i=2,npts-1
         rr=i*h
         dvdr=(vlt(i+1)-vlt(i-1))/(2*h)
         vlsr(i)=dvdr/rr
         vlsi(i)=vlsr(i)
      enddo
      return
      end
      function vcorr(r)
c     version as of 26/3/2007.
c
c     ******************************************************************
c     calc. the coulomb corr. to the energy of the incident particle
c     ******************************************************************
c
      implicit none
      include 'mom.cmb'
      real*8 vcorr,xr,fracp,fracm,r
      integer id
      vcorr=0.d0
      if(zzp.eq.0.d0) return
      if(r.lt.rmax) goto 10
      vcorr=2.d0*zzp/(emu*r)
      return
   10 xr=10.d0*r
      id=xr
      if(id.gt.0) goto 14
      id=1
      fracp=0.d0
      fracm=1.d0
      goto 16
   14 fracp=xr-id
      fracm=1.d0-fracp
   16 vcorr=fracm*vcoul(id)+fracp*vcoul(id+1)
      return
      end
      DOUBLE PRECISION FUNCTION SDOT(N,SX,INCX,SY,INCY)
      implicit none
C**** CALCULATION OF DOT PRODUCT *******
C**** DOUBLE PRECISION *****************************************
      real*8  SX(1),SY(1),STEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
      STEMP=0.0D0
      SDOT=0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GOTO 20
      IX=1
      IY=1
      IF(INCX.LT.0)IX=(-N+1)*INCX+1
      IF(INCY.LT.0)IY=(-N+1)*INCY+1
      DO 10 I=1,N
      STEMP=STEMP+SX(IX)*SY(IY)
      IX=IX+INCX
      IY=IY+INCY
 10   CONTINUE
       SDOT=STEMP
       RETURN
c 20  M=MOD(N,5)
  20  M=N-5*(N/5)
      IF(M.EQ.0)GOTO 40
      DO 30 I=1,M
      STEMP=STEMP+SX(I)*SY(I)
  30  CONTINUE
      IF(N.LT.5)GOTO 60
  40  MP1=M+1
      DO 50 I=MP1,N,5
      STEMP=STEMP+SX(I)*SY(I)+SX(I+1)*SY(I+1)+SX(I+2)*SY(I+2)
     1                       +SX(I+3)*SY(I+3)+SX(I+4)*SY(I+4)
  50  CONTINUE
  60  SDOT=STEMP
      RETURN
      END
      DOUBLE PRECISION FUNCTION SSUM (N,SX,INCX)
      implicit none
      integer N,INCX,I
      REAL*8                 SX(N),STEMP
      integer NINCX,M,MP1
      sSUM=0.0D0
      STEMP=0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
      NINCX=N*INCX
      DO 10 I=1,NINCX,INCX
      STEMP=STEMP+SX(I)
  10  CONTINUE
       SSUM=STEMP
      RETURN
c 20  M=MOD(N,6)
  20  M=N-6*(N/6)
      IF(M.EQ.0)GOTO 40
      DO 30 I=1,M
      STEMP=STEMP+SX(I)
   30 CONTINUE
      IF(N.LT.6)GOTO 60
  40  MP1=M+1
      DO 50 I=MP1,N,6
      STEMP=STEMP+SX(I)+SX(I+1)+SX(I+2)+SX(I+3)
     1                 +SX(I+4)+SX(I+5)
  50  CONTINUE
  60   SSUM=STEMP
      RETURN
      END
      block data
c     version as of 12/09/01.
c     holds the default values for input parameters
c     and the jlm polynomial parameterisation as revised in PRC 58, 1118
c     (1998)
      implicit none
      include 'mom.cmb'
      data h/.1d0/
      data pi/3.141592653589793d0/
      data npts/182/
      data zp/0.d0/
      data z,a/0.d0,0.d0/
      data irel/0/
      data al0r,altr/1.d0,1.d0/
      data al0i,alti/1.d0,1.d0/
      data tjlmr,tjlmi,efermi,icorv,icorw
     &                    /1.25d0,1.35d0,0.d0,1,1/
      data kmw/0/
      data xlv,xlw/1.d0,1.d0/
      data ilamaxi/3/
      data ac/-974.d0,    7097.d0, -19530.d0,
     &        11.26d0,   -125.7d0,    418.d0,
     &       -.0425d0,    .5853d0,  -2.054d0/
      data bc/360.1d0,   -2691.d0,   7733.d0,
     &       -5.224d0,     51.3d0,  -171.7d0,
     &       .02051d0,    -.247d0,   .8846d0/
      data cc/4.557d0,   -2.051d0,  -65.09d0,
     &     -.005291d0,   -.4906d0,   3.095d0,
     &        .6108d-5, .001812d0,  -.0119d0/
c    cut  120 mev parconv-9 with D=625.0 MeV**2
      data dhc/-.65986d3,   .11437d5, -.74505d5,   .17609d6,
     &          .10768d2,  -.29076d3,  .22068d4,  -.54579d4,
     &         -.78863d-1,  .24430d1, -.19926d2,   .51127d2,
     &          .18755d-3, -.62028d-2, .51754d-1, -.13386d0/
c    cut  115 mev 2par fixes par derivee en 2 pts 1parconv-4.15
      data fhc/.45959d3,   -.76929d4,   .55250d5,   -.14373d6,
     &        -.64399d1,    .14639d3,  -.11121d4,    .30382d4,
     &         .40403d-1,  -.10244d1,   .79667d1,   -.22202d2,
     &        -.90086d-4,   .23367d-2, -.18008d-1,   .50258d-1/
      end
