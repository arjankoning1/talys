subroutine dtheory(Zix, Nix, E)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Theoretical calculation of average neutron spacings
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
!   sgl        ! single precision kind
!   dbl        ! double precision kind
! All global variables
!   numJ       ! maximum J - value
!   numl       ! number of l values
!   numN       ! maximum number of neutrons from initial compound nucleus
! Variables for main
!   Ltarget    ! excited level of target
! Variables for level density
!   ldmodel    ! level density model
! Variables for incident channel
!   lmaxinc    ! maximal l - value for transm. coeff. for incident channel
! Variables to normalize compound nucleus cross section
!   pardif     ! difference between target and compound nucleus parity
! Variables for resonance  parameters
!   Dl         ! mean resonance spacing per l value
!   Dlj        ! mean resonance spacing per J, l value
! Variables for levels
!   jdis       ! spin of level
!   parlev     ! parity of level
! Variables for masses
!   S          ! separation energy
!
! *** Declaration of local data
!
  implicit none
  integer   :: J                    ! spin of level
  integer   :: J2                   ! 2 * J
  integer   :: J2b                  ! 2 * start of J summation
  integer   :: J2e                  ! 2 * end of J summation
  integer   :: jj2                  ! 2 * j
  integer   :: jj2beg               ! 2 * start of j summation
  integer   :: jj2end               ! 2 * end of j summation
  integer   :: l                    ! multipolarity
  integer   :: L0                   ! excited level of target
  integer   :: l2                   ! 2 * l
  integer   :: l2beg                ! 2 * start of l summation
  integer   :: l2end                ! 2 * end of l summation
  integer   :: lmaxdth              ! maximal l-value for transmission coefficients for  incident channel
  integer   :: Nix                  ! neutron number index for residual nucleus
  integer   :: Nres                 ! maximal neutron number index for residual nucleus
  integer   :: parity               ! parity
  integer   :: parspin2i            ! 2 * particle spin for incident channel
  integer   :: tpar                 ! target parity
  integer   :: tspin2               ! 2 * target spin
  integer   :: Zix                  ! charge number index for residual nucleus
  real(sgl) :: E                    ! incident energy
  real(sgl) :: tspin                ! target spin
  real(dbl) :: density              ! level density
  real(dbl) :: rho(0:numl, 0:numJ)  ! integrated level density
  real(dbl) :: rhosum               ! help variable
!
! ********************** Level density parameters **********************
!
! levels   : subroutine for discrete levels
!
  do l = 0, numl
    Dl(l) = 0.
    do J = 0, numJ
      Dlj(l, J) = 0.
      rho(l, J) = 0.
    enddo
  enddo
  Nres = min(numN, Nix + 1)
  call levels(Zix, Nres)
  if (Zix == 0 .and. Nix == 0) then
    L0 = Ltarget
  else
    L0 = 0
  endif
  tspin = jdis(Zix, Nres, L0)
  tspin2 = int(2. * jdis(Zix, Nres, L0))
  tpar = parlev(Zix, Nres, L0)
  parspin2i = 1
  lmaxdth = max(lmaxinc, 5)
  J2b = mod(int(2. * (tspin + 0.5)), 2)
  J2e = int(2 * (lmaxdth + 0.5 + tspin))
  J2e = min(J2e, 2 * numJ)
  do parity = - 1, 1, 2
    pardif = abs(tpar - parity) / 2
    do J2 = J2b, J2e, 2
      J = J2 / 2
      jj2beg = abs(J2 - tspin2)
      jj2end = J2 + tspin2
      do jj2 = jj2beg, jj2end, 2
        l2beg = abs(jj2 - parspin2i)
        l2end = jj2 + parspin2i
        l2end = min(l2end, 2 * lmaxdth)
        do l2 = l2beg, l2end, 2
          l = l2 / 2
          if (mod(l, 2) /= pardif) cycle
          rho(l, J) = density(Zix, Nix, max(0., real(S(Zix, Nix, 1)) + E), 0.5 * J2, parity, 0, ldmodel(Zix, Nix))
          Dlj(l, J) = real(1.e6 / rho(l, J))
        enddo
      enddo
    enddo
  enddo
  do l = 0, numl
    rhosum = 0.
    do J = 0, numJ
      rhosum = rhosum + rho(l, J)
    enddo
    if (rhosum >= 1.e-10) Dl(l) = real(1.e6 / rhosum)
  enddo
  return
end subroutine dtheory
! Copyright A.J. Koning 2021
