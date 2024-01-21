      subroutine gnrl(galpha,gbeta,gam,mu,nu,lamda,s,df,id)
c
c +---------------------------------------------------------------------
c | Author: Gilles Noguere
c | Date  : December 23, 2012
c | Task  : Fluctuation integrals for unresolved resonances from NJOY
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer mu,nu,lamda,id,j,k,l
      real    galpha,gbeta,gam,s,df,xj,qp(10,4),qw(10,4)
c
c **********************************************************************
c
      data qw/1.1120413E-01,2.3546798E-01,2.8440987E-01,2.2419127E-01,
     + .10967668,.030493789,.0042930874,2.5827047E-04,
     + 4.9031965E-06,1.4079206E-08,.033773418,.079932171,.12835937,
     + .17652616,.21347043,.21154965,.13365186,.022630659,
     + 1.6313638E-05,2.745383E-31,3.3376214E-04,.018506108,.12309946,
     + .29918923,.33431475,.17766657,.042695894,4.0760575E-03,
     + 1.1766115E-04,5.0989546E-07,1.7623788E-03,.021517749,
     + .080979849,.18797998,.30156335,.29616091,.10775649,
     + 2.5171914E-03,8.9630388E-10,0.0/
      data qp/3.0013465E-03,7.8592886E-02,4.3282415E-01,1.3345267,
     + 3.0481846,5.8263198,9.9452656,1.5782128E+01,23.996824,
     + 36.216208,1.3219203E-02,7.2349624E-02,.19089473,.39528842,
     + .74083443,1.3498293,2.5297983,5.2384894,13.821772,
     + 75.647525,1.0004488E-03,.026197629,.14427472,.44484223,
     + 1.0160615,1.9421066,3.3150885,5.2607092,7.9989414,
     + 12.072069,.013219203,.072349624,.19089473,.39528842,
     + .74083443,1.3498293,2.5297983,5.2384894,13.821772,
     + 75.647525/
c
c ****************** calculation of s  *********************************
c
      s=0.
      if (galpha.le.0.) return
      if (gam.le.0.) return
      if (gbeta.lt.0.) return
      if (gbeta.gt.0..and.df.lt.0.) return
c
      if (gbeta.eq.0..and.df.eq.0.) then
        if (id.eq.1) then
          do j=1,10
            s=s+qw(j,mu)*qp(j,mu)*qp(j,mu)/(galpha*qp(j,mu)+gam)
          enddo
        else if (id.eq.2) then
          do j=1,10
            s=s+qw(j,mu)*qp(j,mu)/(galpha*qp(j,mu)+gam)
          enddo
        endif
        return
      endif
c
      if (gbeta.eq.0..and.df.gt.0.) then
        if (id.eq.1) then
          do j=1,10
            xj=qp(j,mu)
            do k=1,10
              s=s+qw(j,mu)*qw(k,lamda)*xj*xj/(galpha*xj+gam
     +           +df*qp(k,lamda))
            enddo
          enddo
         else if (id.eq.2) then
            do j=1,10
               xj=qp(j,mu)
               do k=1,10
                  s=s+qw(j,mu)*qw(k,lamda)*xj/(galpha*xj+gam
     +              +df*qp(k,lamda))
               enddo
            enddo
         endif
         return
      endif
c
      if (gbeta.gt.0..and.df.eq.0.) then
         if (id.eq.1) then
            do j=1,10
               xj=qp(j,mu)
               do k=1,10
                  s=s+qw(j,mu)*qw(k,nu)*xj*xj/(galpha*xj
     +              +gbeta*qp(k,nu)+gam)
               enddo
            enddo
         else if (id.eq.2) then
            do j=1,10
               xj=qp(j,mu)
               do k=1,10
                  s=s+qw(j,mu)*qw(k,nu)*xj/(galpha*xj
     &              +gbeta*qp(k,nu)+gam)
               enddo
            enddo
         else if (id.eq.3) then
            do j=1,10
               do k=1,10
                  s=s+qw(j,mu)*qw(k,nu)*qp(j,mu)*qp(k,nu)
     +              /(galpha*qp(j,mu)+gbeta*qp(k,nu)+gam)
               enddo
            enddo
         endif
         return
      endif
c
      if (gbeta.gt.0..and.df.gt.0.) then
         if (id.eq.1) then
            do j=1,10
               xj=qp(j,mu)
               do k=1,10
                  do l=1,10
                     s=s+qw(j,mu)*qw(k,nu)*qw(l,lamda)*xj*xj
     +                 /(galpha*xj+gbeta*qp(k,nu)+gam+df*qp(l,lamda))
                  enddo
               enddo
            enddo
         else if (id.eq.2) then
            do j=1,10
               do k=1,10
                  do l=1,10
                     s=s+qw(j,mu)*qw(k,nu)*qw(l,lamda)*qp(j,mu)/
     +                 (galpha*qp(j,mu)+gbeta*qp(k,nu)+gam
     +                 +df*qp(l,lamda))
                  enddo
               enddo
            enddo
         else if (id.eq.3) then
            do j=1,10
               do k=1,10
                  do l=1,10
                     s=s+qw(j,mu)*qw(k,nu)*qw(l,lamda)*qp(j,mu)
     +                 *qp(k,nu)/(galpha*qp(j,mu)+gbeta*qp(k,nu)
     +                 +gam+df*qp(l,lamda))
                  enddo
               enddo
            enddo
         endif
         return
      endif
      return
      end
