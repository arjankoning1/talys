subroutine pseudo_resonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Compound nucleus Resonances
!
! Author    : Arjan Koning
!
! 2024-05-08: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! All global variables
!   numlev2        ! maximum number of levels
! Variables for levels
!   edis           ! energy of level
!   jdis           ! spin of level
!   parlev         ! parity of level
!
! *** Declaration of local data
!
  implicit none
  integer           :: i                 ! counter
  integer           :: k                 ! counter
  integer           :: nres
  integer           :: Pres(0:numlev2)     ! 
  real(sgl)         :: fgr(0:numresgrid)    ! 
  real(sgl)         :: Eres(0:numlev2)     ! 
! real(sgl)         :: Jres(0:numlev2)     ! 
! real(sgl)         :: sdis(0:numlev2)     ! 
  real(sgl)         :: Ec
  real(sgl)         :: Er
  real(sgl)         :: E
  real(sgl)         :: degrid
  real(sgl)         :: Emax
  real(sgl)         :: peak
  real(sgl)         :: reswidth
  real(sgl)         :: weight
! real(sgl)         :: J0
! real(sgl)         :: sc
  real(sgl)         :: maxdis
  real(sgl)         :: fac1
  real(sgl)         :: hgam
  real(sgl)         :: lorentz
!
! ******************** Default nuclear levels **************************
!
  peak = 1.
  reswidth = 0.1
! J0 = targetspin + parspin(k0)
  Eres = 0.
! Jres = 0.
! Pres = 0.
  k = 0
! sc = 1.
  do i = 1, nlevmax2(0, 0)
    Ec = edis(0,0,i)
    Er = Ec - Q(0)
    if (Er > 0.) then
      k = k + 1
      Eres(k) = Er
!     Jres(k) = jdis(0,0,i)
!     Pres(k) = parlev(0,0,i)
!     sdis(k) = spindis(sc,Jres(k))
    endif
  enddo
  Nres = k
  Emax= 10.
  degrid = Emax / numresgrid
  Eresgrid = 0.
  resgrid = 0.
  fgr = 0.
  fac1 = 1. / pi
  hgam = 0.5 * reswidth
  maxdis = 0.
  do i = 1, numresgrid
    Eresgrid(i) = i*degrid
    do k = 1, Nres
      E = abs(Eresgrid(i) - Eres(k))
      lorentz = fac1 * hgam / ( E*E + hgam*hgam)
!     fgr(i) = fgr(i) + sdis(k) * lorentz
!     weight = 1. / (sqrt(real(k)))
      weight = 1. / real(k)
      fgr(i) = fgr(i) + lorentz * weight
    enddo
!   maxdis = max(maxdis, fgr(i))
    maxdis = 1.
  enddo
  do i = 1, numresgrid
    resgrid(i) = 1. +  peak / maxdis * fgr(i)
  enddo
  return
end subroutine pseudo_resonance
! Copyright A.J. Koning 2024
