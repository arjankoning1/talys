C=======================================================================
C
C     GROUPIE COMMON
C
C     2019/5/26 - Added Entire Self-Shielded Table
C                 XCALL = 25 X 6 X MAXGROUP = 3000000 8 byte words
C                 Earlier only 1 group results at a time were in memory.
C                 Saving ALL greatly simplifies the logic
C
C=======================================================================
C
C     PARAMETERS
C
C-----------------------------------------------------------------------
c-----2017/3/07 - INCREASED TO 3,000,000 FROM 600,000.
c-----2023/3/07 - Decreased page size from 3,000,000 to 120,000
c     PARAMETER (MAXPOINT = 3000000) ! Data points in memory
      PARAMETER (MAXPOINT =  120000) ! Data points in memory
c----- 3,000,000 X 8 bytes X 10 arrays = 240,000,000 bytes
c-----   120,000 X 8 bytes X 10 arrays =   9,600,000 bytes
c-----REPORT FOR EACH MAT
      PARAMETER (MAXMAT  =  1000)
c-----MAX. NUMBER OF BANDS
      PARAMETER (MAXBAND = 2)
C-----MAX. NUMBER OF GROUPS
c-----2015/12/22 - Increased to 20,000 from 1,000
      PARAMETER (MAXGROUP = 20000)
c----- 20,000 X 8 bytes X 72 arrays = 11,520,000 bytes
c-----MULTI-GROUP OUTPUT TO ENDF/B FORMAT
c-----2015/12/22 - Increased to 20,000 from 3,000
c-----2019/06/20 - Increased to 100,000 for MF=2 output
      PARAMETER (MAXOUT   = 100000)
c----- 20,000 X 8 bytes X 2 arrays = 320,000 bytes
c-----2020/8/6 - Max # of TART groups (for URR calculation)
      PARAMETER (MAXTART = 1000) ! Actually 616, but close enough
c----- 10 + 1000 + 18,000 + 2000 + 12,000 = 33,100 x 8 bytes
C-----------------------------------------------------------------------
C
C     STORAGE
C
C-----------------------------------------------------------------------
c
C     PAGED STORAGE FOR
C     1 = SPECTRUM
C     2 = TOTAL
C     3 = ELASTIC
C     4 = CAPTURE
C     5 = FISSION
C     6 = OTHER - NOT IN PAGING SYSTEM, ONLY IN MEMORY ARRAYS
C
C-----------------------------------------------------------------------
      COMMON/REPORT/ERBTAB(6,MAXMAT),ERLIB(25,MAXBAND),NBNTAB(MAXMAT),
     1 IZATAB(MAXMAT),LZA
      COMMON/GROUPR/EGROUP(MAXGROUP+1),TOTAV(MAXGROUP+1),
     1 EAV(MAXGROUP+1),AVN(MAXGROUP+1),NGR,NGRP1,NGRP2,IGR
C-----MULTI-BAND STORAGE
      COMMON/BANDID/WTBAND  (MAXBAND,MAXGROUP),
     1              XCBAND(6,MAXBAND,MAXGROUP)
      COMMON/INTNRM/XCINT(25,6),XCNORM(25)
      COMMON/MATERR/ERMAT(25,MAXBAND),ERNOW(25,MAXBAND)
      COMMON/SIMPLE/XCFI(25,6),AVEXP(25),AVNORM(25),SIGMAB(25)
      COMMON/MISSIT/FST(25),DFST(25),DDFST(25),DEAVST(25)
      COMMON/ELPAS1/SHIELD(25),YLOW(5),YHIGH(5),YLOWP1(5),YLOWP2(5),
     1 YLOWP3(5),YHIGHP1(5),YHIGHP2(5),YHIGHP3(5)
      CHARACTER*4 POINT,REACT2,REACT3
      COMMON/ELPASC/POINT(25),REACT2(2,6),REACT3(2,6)
c-----2020/8/6 - Add for MF=2/152 & 153 output including boundaries
      COMMON/URRCOM/XCURRLIM(2,5),EAVURR(MAXTART),XCUALL(3,6,MAXTART),
     1 WTUBAND(MAXBAND,MAXTART),XCUBAND(6,MAXBAND,MAXTART)
C-----ENDF Formatted Output.
      COMMON/GROUPOUT/XOUT(MAXOUT),YOUT(MAXOUT)
C-----BLANK COMMON
      COMMON XPAGE(MAXPOINT,5),YPAGE(MAXPOINT,5),XCALL(25,6,MAXGROUP)
