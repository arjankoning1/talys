subroutine cascade(Zcomp, Ncomp, nex)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Gamma-ray cascade
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
! Definition of single and double precision variables
!   dbl             ! double precision kind
! All global variables
!   numNchan        ! maximum number of outgoing neutron in individual channel description
!   numZchan        ! maximum number of outgoing protons in individual channel description
! Variables for discrete levels
!   flagelectron    ! flag for application of electron conversion coefficient
! Variables for output
!   flaggamdis      ! flag for output of discrete gamma - ray intensities
! Variables for multiple emission
!   mcontrib        ! contribution to emission spectrum
!   xsgamdis        ! discrete gamma - ray cross section
!   xsgamdistot     ! total discrete gamma - ray cross section
!   xspartial       ! emitted cross section flux per energy bin
! Variables for incident channel
!   popdecay        ! decay from population
!   xspop           ! population cross section
!   xspopex         ! population cross section summed over spin and parity
! Variables for levels
!   branchlevel     ! level to which branching takes place
!   branchratio     ! gamma - ray branching ratio to level
!   conv            ! conversion coefficient
!   jdis            ! spin of level
!   nbranch         ! number of branching levels
!   parlev          ! parity of level
!
! *** Declaration of local data
!
  implicit none
  integer   :: i       ! counter
  integer   :: J       ! spin of level
  integer   :: Jres    ! spin of level
  integer   :: k       ! designator for particle
  integer   :: Ncomp   ! neutron number index for compound nucleus
  integer   :: nex     ! excitation energy bin of compound nucleus
  integer   :: parity  ! parity
  integer   :: Pres    ! parity
  integer   :: Zcomp   ! proton number index for compound nucleus
  real(dbl) :: intens  ! total gamma intensity
  real(dbl) :: xsgamma ! help variable
  real(dbl) :: xsJP    ! population cross section for spin and parity
!
! ******************* Gamma-ray cascade ********************************
!
  J = int(jdis(Zcomp, Ncomp, nex))
  parity = parlev(Zcomp, Ncomp, nex)
  xsJP = xspop(Zcomp, Ncomp, nex, J, parity)
  do i = 1, nbranch(Zcomp, Ncomp, nex)
    k = branchlevel(Zcomp, Ncomp, nex, i)
    Jres = int(jdis(Zcomp, Ncomp, k))
    Pres = parlev(Zcomp, Ncomp, k)
    intens = xsJP * branchratio(Zcomp, Ncomp, nex, i)
    xspop(Zcomp, Ncomp, k, Jres, Pres) = xspop(Zcomp, Ncomp, k, Jres, Pres) + intens
    popdecay(0, k, Jres, Pres) = popdecay(0, k, Jres, Pres) + intens
    xspopex(Zcomp, Ncomp, k) = xspopex(Zcomp, Ncomp, k) + intens
    xspopex(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) - intens
    xspartial(0, nex) = xspartial(0, nex) + intens
    if (Zcomp <= numZchan .and. Ncomp <= numNchan) mcontrib(0, nex, k) = intens
!
! ************ Storage of discrete gamma line intensities **************
!
    if (flaggamdis) then
      if (flagelectron) then
        xsgamma = intens / (1. + conv(Zcomp, Ncomp, nex, i))
      else
        xsgamma = intens
      endif
      xsgamdis(Zcomp, Ncomp, nex, k) = xsgamma
      xsgamdistot(Zcomp, Ncomp) = xsgamdistot(Zcomp, Ncomp) + xsgamma
    endif
  enddo
  return
end subroutine cascade
! Copyright A.J. Koning 2021
