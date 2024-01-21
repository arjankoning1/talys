subroutine buenhance(Z00, N00, BFsum, Ebu)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Inelastic breakup enhancement brought by breakup neutrons
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
!   sgl          ! single precision kind
! All global variables
!   numenout     ! maximum number of outgoing energies
! Constants
!   nuc          ! symbol of nucleus
! Variables for main input
!   Atarget      ! mass number of target nucleus
!   Ninit       ! neutron number of initial compound nucleus
!   Zinit       ! charge number of initial compound nucleus
!   Ztarget      ! charge number of target nucleus
! Variables for preequilibrium
!   ebubin            ! outgoing breakup nucleon energy bin for integration
!   ENHratio          !  breakup nucleons enhancing reaction cross
!   xsBF              ! nucleon inelastic breakup cross section
!   xsBUnuc           ! nucleon breakup cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: N                     ! neutron number of residual nucleus
  integer   :: N00                   ! neutron number for residual nucleus
  integer   :: nen                   ! energy counter
  integer   :: Nix                   ! neutron number index for residual nucleus
  integer   :: type                  ! particle type
  integer   :: Z                     ! charge number of target nucleus
  integer   :: Z00                   ! charge number for residual nucleus
  integer   :: Zix                   ! charge number index for residual nucleus
  real(sgl) :: BCin                  ! effective incident/outgoing Coulomb barrier
  real(sgl) :: BCout(2)              ! effective incident/outgoing Coulomb barrier
  real(sgl) :: Bdeut                 ! deuteon binding energy
  real(sgl) :: BFenhance(2)          ! inelastic breakup enhancement brought by breakup  nucleons to any (Z,N) residual n
  real(sgl) :: BFsum                 ! total inelastic breakup enhancement to (Z,N) residual  nucleus population
  real(sgl) :: conv03(2)             ! help variable for breakup enhancement  calculations
  real(sgl) :: E0n03                 ! break-up energy
  real(sgl) :: Ebreak(2)             ! Centroid energy of the breakup nucleon energy  distributions, Kalbach 2003, in Lab
  real(sgl) :: Ebu                   !
  real(sgl) :: Emaxn                 ! maximum breakup nucleons energy in Laboratory System
  real(sgl) :: En                    ! outgoing energy
  real(sgl) :: enhance03(2)          ! help variable for breakup enhancement  calculations
  real(sgl) :: Gauss                 ! function for Gaussian
  real(sgl) :: sumGauss(2)           ! help variable for ckhecking the  calculation procedures
  real(sgl) :: sumtest(2)            ! check of breakup spectrum
  real(sgl) :: sumtest0(2, numenout) ! help variable for ckhecking the  calculation procedures
  real(sgl) :: term03                ! help variable for breakup enhancement  calculations
  real(sgl) :: w03                   ! width
  real(sgl) :: width                 ! Full width at half maximum of the breakup nucleon  energy distribution, Kalbach 20
!
! ************************** Avrigeanu model ***************************
!
! DEUTERON Break-up model by M. and V. Avrigeanu, PRC95,024607(2017)
! Inelastic breakup enhancement, M.Avrigeanu et al., PRC94,014606(2016)
!
!   distributions, Kalbach 2003, in Laboratory System
!   energy distribution, Kalbach 2003
!   ratios, Eq. (2) from PRC94,014606(2016)
!   n + Atarget  sig(n,Z,A,Eout)/sig_Total(n,Enout);
!   p + Atarget  sig(p,Z,A,Eout)/sig_Reaction(p,Epout)
!   nucleons to any (Z,N) residual nucleus from d+Atarge
!   interaction
!   nucleus population
!
  BFsum = 0.
!
  Ebreak = 0.
  conv03 = 0.
  enhance03 = 0.
  sumtest = 0.
  sumGauss = 0.
  BFenhance = 0.
!
  sumtest0 = 0.
!
  Z = Z00
  N = N00
  Bdeut = 2.225
  BCin = Ztarget / 9.5
  BCout(1) = 0.
  BCout(2) = BCin
  Emaxn = Ebu * (Atarget + 1.) / (Atarget + 2.) - Bdeut * (Atarget + 1.) / Atarget
!
!   Breakup threshold
!
  if(Emaxn <= 0.d+00) then
    Emaxn = 0.
    write(8, * )" deuteron energy lower than Bd, no breakup"
    go to 499
  endif
!
!----------------------     breakup nucleon type  ----------------------
!
  do type = 1, 2
!
    Ebreak(type) = 0.5 * (Atarget + 1.) / (Atarget + 2.) * Ebu + &
        0.5 * (Atarget + 1.) / Atarget * ( - Bdeut - BCin + 2. * BCout(type))
    width = 1.15 + 0.12 * Ebu - Atarget / 140.
!
    if(width <= 0.d+00) then
      width = 0.
      go to 499
    endif
!
    if(Ebreak(type) <= 0.d+00) then
      Ebreak(type) = 0.01
    endif
!
    E0n03 = Ebreak(type)
    w03 = width
!
    conv03(type)     = 0.
    sumGauss(type)   = 0.
    BFenhance(type)  = 0.
!
!----------------------     breakup nucleon energy   -------------------
!
!   energy bin for integration in Eq.(2), PRC94,014606.
!
     Zix = Zinit - Z
     Nix = Ninit - N
     do nen = 1, numenout
       En = ebubin * nen
       if(En <= Emaxn) then
         sumtest0(type, nen) = GAUSS(En, E0n03, w03)
         term03 = ENHratio(type, Zix, Nix, nen) * GAUSS(En, E0n03, w03) * ebubin
!
!   BREAKUP THRESHOLD:  En>Emax BU
!
      else
        sumtest0(type, nen) = 0.
        term03 = 0.
      endif
!
      sumGauss(type) = sumGauss(type) + sumtest0(type, nen) * ebubin
      conv03(type) = conv03(type) + term03
     enddo
!
!--------------------   end   breakup nucleon energy   -----------------
!
    enhance03(type) = xsBF(type) * conv03(type)
    if(sumGauss(type) == 0.d+00) cycle
    BFenhance(type) = enhance03(type) / sumGauss(type)
    BFsum = BFsum + BFenhance(TYPE)
!
  enddo
!
!----------------------     end   breakup nucleon type  ----------------
!
499      continue
!
    do type = 1, 2
     if (sumGauss(type) == 0.d+00) cycle
     do nen = 1, numenout
       sumtest0(type, nen) = xsBUnuc(type) * sumtest0(type, nen) / sumGauss(type)
       sumtest(type) = sumtest(type) + sumtest0(type, nen) * ebubin
     enddo
    write(8, 999)Ebu, type, BFenhance(type)
    enddo
     write(8, *)'                 BFsum=', BFsum, '  for  Z=', Z, ' A=', Z+N, nuc(Z)
!
999     Format(2x, f8.3, 1x, I3, 3x, f12.5, 3x, f8.3, 1x, f8.3, 1x, f8.3, 1x, f8.3, &
         4x, f8.3, 3x, f8.3, 4x, f8.3, 4x, f8.3, 4x, f8.3, 4x, f8.3)
!
  return
end subroutine buenhance
! Copyright A.J. Koning 2021
