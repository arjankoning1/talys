c=======================================================================
c
c     RECENT: All parameters for LRF=7 are here
c
c=======================================================================
C
C     PARAMETERS
C
c-----------------------------------------------------------------------
      PARAMETER (MaxNgroup = 100)    ! (L,J) sequences = cumulative
      PARAMETER (MaxNpp    = 10)     ! 10 channels = Not cumulative
      PARAMETER (MaxNchan  = 10)     ! 10 channels = Not cumulative
      PARAMETER (MaxNres   = 100000) ! # of resonances = cumulative
      PARAMETER (MaxNtriag = 55)     ! 10      x 11         /2 = 55
c     MaxNtriag                      = MaxNchan*(MaxNchan+1)/2
c-----------------------------------------------------------------------
c
c     COMMON
c
c     Labelled common follows
c     (blank common is ALL in recent.h)
c
c-----------------------------------------------------------------------
c-----Options from RECENT - use IMEDIT for output listing
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
C-----MRMLV.F
      COMMON/MRMLVCOM/Awr7,Su
c-----MRMLW.F
c-----2017/4/13 - Deleted Lrf
      COMMON/MRMLWCOM/Npp,Ngroup,Nres,Nro7,Naps7,Kg,Minr,Maxr
c-----------------------------------------------------------------------
c
c     ALL arrays follow - this is the fixed field equivalent of
c     the dummy array A used by SAMRML - see the above parameters
c     for the definitions of ALL of the below fixed field arrays.
c
c-----------------------------------------------------------------------
      common/rmldataf/Ema(MaxNpp),Emb(MaxNpp),Spina(MaxNpp),
     1 Spinb(MaxNpp),Qqq(MaxNpp),Pa(MaxNpp),Pb(MaxNpp),
     2 Spin(MaxNgroup),Parity(MaxNgroup),Eres(MaxNres),
     3 Gamgam(MaxNres),Goj(MaxNgroup),Crss(MaxNpp)
c
      common/rmldatai/Nchan(MaxNgroup),Kza(MaxNpp),Kzb(MaxNpp),
     1 Ishift(MaxNpp),Lpent(MaxNpp),Mt7(MaxNpp),Nresg(MaxNgroup),
     2 Nent(MaxNgroup),Next(MaxNgroup)
c
      common/rmldat2f/Chspin7(MaxNchan,MaxNgroup),
     1                Bndry  (MaxNchan,MaxNgroup),
     1                Rdeff  (MaxNchan,MaxNgroup),
     1                Rdtru  (MaxNchan,MaxNgroup),
     1                Zke    (MaxNchan,MaxNgroup),
     1                Zkfe   (MaxNchan,MaxNgroup),
     1                Zkte   (MaxNchan,MaxNgroup),
     1                Zeta7  (MaxNchan,MaxNgroup),
     1                Echan  (MaxNchan,MaxNgroup),
     1                Gamma  (MaxNchan,MaxNres),
     1                Betapr (MaxNchan,MaxNres),
     1                Gbetpr (3       ,MaxNres)
c
      common/rmldat2i/Ipp   (MaxNchan,MaxNgroup),
     1                Lspin (MaxNchan,MaxNgroup)
c
      common/rmldat3f/Beta7(MaxNtriag ,MaxNres),
     1 Alphar(MaxNres),Alphai(MaxNres),Difen(MaxNres),
     2 Xden(MaxNres),Sinsqr(MaxNchan),Sin2ph(MaxNchan),
     3 Sinphi(MaxNchan),Cosphi(MaxNchan),Cscs(2,MaxNtriag),
     4 Rootp(MaxNchan),Elinvr(MaxNchan),Elinvi(MaxNchan),
     5                  Xxxxr(MaxNtriag),Xxxxi(MaxNtriag),
     6 Xqr(MaxNchan,MaxNchan),Xqi(MaxNchan,MaxNchan),
     7 Yinv(2,MaxNtriag),Rmat(2,MaxNtriag),Ymat(2,MaxNtriag)
