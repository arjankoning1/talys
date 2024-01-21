subroutine bonetti
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Determination of effective absorption optical potential
!
! Author    : Marieke Duijvestijn and Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl         ! single precision kind
! All global variables
!   numen       ! maximum number of outgoing energies
! Variables for input energies
!   enincmax    ! maximum incident energy
! Variables for main input
!   Atarget     ! mass number of target nucleus
!   k0          ! index of incident particle
! Variables for nuclides
!   Nindex      ! neutron number index for residual nucleus
!   Zindex      ! charge number index for residual nucleus
!   ZZ          ! charge number of residual nucleus
! Constants
!   onethird    ! 1 / 3
!   parZ        ! charge number of particle
! Variables for ECIS
!   prodZ       ! product of charges of projectile and target nucleus
! Variables for optical model
!   aw          ! imaginary volume diffuseness
!   awd         ! imaginary surface diffuseness
!   rw          ! imaginary volume radius
!   rwd         ! imaginary surface radius
!   w           ! imaginary volume depth
!   wd          ! imaginary surface depth
! Variables for optical model
!   jlmexist    ! flag for existence of tabulated radial matter density
! Variables for JLM
!   normjlm     ! JLM potential normalization factors
!   potjlm      ! JLM potential depth values
!   rhojlmn     ! density for neutrons
!   rhojlmp     ! density for protons
! Variables for masses
!   S           ! separation energy
! Variables for exciton model initialization
!   wvol        ! absorption part of the optical potential averaged over the volume
!
! *** Declaration of local data
!
  implicit none
  integer   :: i       ! counter
  integer   :: k       ! designator for particle
  integer   :: nen     ! energy counter
  integer   :: nenend  ! help variable
  integer   :: Nix     ! neutron number index for residual nucleus
  integer   :: nrbins  ! number of continuum excitation energy bins
  integer   :: Z       ! charge number of target nucleus
  integer   :: Zix     ! charge number index for residual nucleus
  real(sgl) :: dr      ! integration bin width
  real(sgl) :: e       ! energy
  real(sgl) :: emax    ! maximal emission energy within bin decay
  real(sgl) :: expo    ! help variable
  real(sgl) :: fwssurf ! derivative of the imaginary surface optical potential  form factor
  real(sgl) :: fwsvol  ! imaginary volume optical potential form factor
  real(sgl) :: radd    ! imaginary surface nuclear radius
  real(sgl) :: radv    ! imaginary volume nuclear radius
  real(sgl) :: Rjlm    ! JLM normalization factor
  real(sgl) :: rr      ! running variable in integration over the radius
  real(sgl) :: sum1    ! help variable
  real(sgl) :: sum2    ! help variable
  real(sgl) :: term1   ! help variable
  real(sgl) :: term2   ! help variable
  real(sgl) :: wjlm    ! help variable
!
! ********** Calculation of volume absorptive optical potential ********
!
! mom        : subroutine for microscopic optical model (Eric Bauge)
! optical    : subroutine for determination of optical potential
!
  wvol = 0.
  Zix = Zindex(0, 0, k0)
  Nix = Nindex(0, 0, k0)
  emax = enincmax + S(0, 0, k0)
  nenend = 10 * min(numen, int(emax + 1.))
  Rjlm = 1.
  do k = 1, 2
!
! A. Phenomenological potential (always done, also for the case of normalizing the JLM OMP)
!
    nrbins = 50
    dr = 20. / nrbins
    do nen = - 80, nenend
      if (jlmexist(Zix, Nix, k) .and. nen > 1) cycle
      e = 0.1 * real(nen)
      call optical(Zix, Nix, k, e)
      radv = rw * Atarget **onethird
      radd = rwd * Atarget **onethird
      sum1 = 0.
      sum2 = 0.
      do i = 1, nrbins
        rr = (i - 0.5) * dr
        expo = (rr - radv) / aw
        if (expo <= 80.) then
          fwsvol = 1. / (1. + exp(expo))
        else
          fwsvol = 0.
        endif
        expo = (rr - radd) / awd
        if (expo <= 80.) then
          fwssurf = - exp(expo) / (awd * (1. + exp(expo) **2))
        else
          fwssurf = 0.
        endif
        term2 = fwsvol * (rr **2) * dr
        term1 = term2 * (w * fwsvol - 4. * awd * wd * fwssurf)
        sum1 = sum1 + term1
        sum2 = sum2 + term2
      enddo
      if (sum2 /= 0.) wvol(k, nen) = sum1 / sum2
    enddo
!
! B. JLM potential
!
    if (jlmexist(Zix, Nix, k)) then
      nrbins = 122
      dr = 12.2 / nrbins
      Z = ZZ(Zix, Nix, k)
      prodZ = Z * parZ(k)
      do nen = nenend, - 80, - 1
!
! For negative energies, the phenomenological OMP is adopted, normalized at the first positive energy point.
!
        if (nen < 1) then
          wvol(k, nen) = Rjlm * wvol(k, nen)
          cycle
        endif
        e = 0.1 * real(nen)
        call mom(Zix, Nix, dble(prodZ), dble(e))
        sum1 = 0.
        sum2 = 0.
        do i = 1, nrbins
          rr = (i - 0.5) * dr
          if (k == 1) then
            term2 = rhojlmn(Zix, Nix, i, 1) * (rr **2) * dr
          else
            term2 = rhojlmp(Zix, Nix, i, 1) * (rr **2) * dr
          endif
          term1 = - term2 * normjlm(Zix, Nix, 2) * potjlm(Zix, Nix, i, 2)
          sum1 = sum1 + term1
          sum2 = sum2 + term2
        enddo
        if (nen == 1) then
          Rjlm = 1.
          if (sum2 /= 0.) then
            wjlm = sum1 / sum2
            Rjlm = wjlm / wvol(k, nen)
          endif
        endif
        if (sum2 /= 0.) wvol(k, nen) = sum1 / sum2
      enddo
    endif
    do nen = - 200, - 81
      wvol(k, nen) = wvol(k, - 80)
    enddo
  enddo
  return
end subroutine bonetti
! Copyright A.J. Koning 2021
