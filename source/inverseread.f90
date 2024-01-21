subroutine inverseread(Zcomp, Ncomp)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read ECIS results for outgoing particles and energy grid
!
! Author    : Arjan Koning, Stephane Hilaire, Eric Bauge and Pascal Romain
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
!   dbl           ! double precision kind
! All global variables
!   numl          ! number of l values
! Variables for numerics
!   transeps      ! absolute limit for transmission coefficient
! Variables for direct reactions
!   flagrot       ! flag for use of rotational optical model per outgoing p
! Variables for energy grid
!   ebegin        ! first energy point of energy grid
!   ecisstatus    ! status of ECIS file
!   eendmax       ! last energy point of energy grid
!   translimit    ! limit for transmission coefficient
! Variables for energies
!   eend          ! last energy point of energy grid
! Variables for inverse channel data
!   csfile        ! file with inverse reaction cross sections
!   Tjl           ! transmission coefficient per particle, energy, spin and l - value
!   Tl            ! transmission coefficients per particle, energy and l - value
!   transfile     ! file with transmission coefficients
!   xselas        ! total elastic cross section (shape + compound)
!   xsopt         ! optical model reaction cross section
!   xsreac        ! reaction cross section
!   xstot         ! total cross section (neutrons only)
!  Variables for gamma-ray strength functions
!   lmax          ! maximal l - value for transmission coefficients
! Variables for nuclides
!   Nindex        ! neutron number index for residual nucleus
!   parskip       ! logical to skip outgoing particle
!   Zindex        ! charge number index for residual nucleus
! Constants
!   parspin       ! spin of particle
! Variables for levels
!   jdis          ! spin of level
! Variables for deformation parameters
!   colltype      ! type of collectivity (D, V or R)
!
! *** Declaration of local data
!
  implicit none
  logical   :: lexist      ! logical to determine existence
  integer   :: i           ! level
  integer   :: ispin       ! spin index
  integer   :: k           ! designator for particle
  integer   :: l           ! multipolarity
  integer   :: lev         ! level number
  integer   :: Ncomp       ! neutron number index for compound nucleus
  integer   :: nen         ! energy counter
  integer   :: Nix         ! neutron number index for residual nucleus
  integer   :: nJ          ! number of total J values for transmission coefficients
  integer   :: nS          ! number of states
  integer   :: type        ! particle type
  integer   :: Zcomp       ! proton number index for compound nucleus
  integer   :: Zix         ! charge number index for residual nucleus
  real(sgl) :: factor      ! multiplication factor
  real(sgl) :: groundspin2 ! 2 * spin of ground state
  real(sgl) :: jres        ! j-value
  real(sgl) :: rj          ! help variable
  real(dbl) :: Tcoef       ! transmission coefficients as a function of spin and  l-value
  real(dbl) :: teps        ! help variable
  real(dbl) :: xs          ! help variable
!
! ************ Read total, reaction and elastic cross section **********
!
  Tjl = 0.
  Tl = 0.
  xselas = 0.
  xsopt = 0.
  xsreac = 0.
  xstot = 0.
  inquire (file = csfile, exist = lexist)
  if ( .not. lexist) then
    write(*, '(" TALYS-error: The first calculation of a run", " should always be done with ecissave y and", &
 &    " eciscalc y. Non-existent file: ", a13)') csfile
    stop
  endif
  open (unit = 3, file = csfile, status = 'unknown')
  do type = 1, 6
    if (parskip(type)) cycle
    do nen = ebegin(type), eendmax(type)
      read(3, '()')
      if (type == 1) then
        read(3, * ) xs
        xstot(type, nen) = max(real(xs), 0.)
      endif
      read(3, * ) xs
      xsreac(type, nen) = max(real(xs), 0.)
      xsopt(type, nen) = xsreac(type, nen)
      if (type == 1) then
        read(3, * ) xs
        xselas(type, nen) = max(real(xs), 0.)
      endif
    enddo
  enddo
  close (unit = 3, status = ecisstatus)
!
! ******************* Read transmission coefficients *******************
!
! The transmission coefficient Tjl has four indices: particle type, energy, spin and  l-value.
! For spin-1/2 particles, we use the indices -1 and 1 for the two spin values.
! For spin-1 particles, we use -1, 0 and 1 and for spin-0 particles we use 0 only.
!
! For rotational nuclei, the rotational transmission coefficients are transformed into into spherical equivalents.
!
  open (unit = 7, file = transfile, status = 'unknown')
  do type = 1, 6
    if (parskip(type)) cycle
    Zix = Zindex(Zcomp, Ncomp, type)
    Nix = Nindex(Zcomp, Ncomp, type)
    groundspin2 = int(2. * jdis(Zix, Nix, 0))
    do nen = ebegin(type), eendmax(type)
      read(7, '(55x, i5)') nJ
      do i = 1, nJ
        read(7, '(f10.1, 5x, i5)') rj, nS
        do k = 1, nS
          read(7, '(i3, i6, f9.1, e20.8)') lev, l, jres, Tcoef
          if (l > numl) cycle
          if (lev == 1) then
            if (colltype(Zix, Nix) /= 'S' .and. flagrot(type)) then
              factor = (2. * rj + 1.) / (2. * jres + 1.) / (groundspin2 + 1.)
            else
              factor = 1.
            endif
            if (parspin(type) == 0.5) then
              ispin = int(2. * (jres - real(l)))
            else
              ispin = int(jres - real(l))
            endif
            Tjl(type, nen, ispin, l) = Tjl(type, nen, ispin, l) + factor * max(real(Tcoef), 0.)
          endif
        enddo
      enddo
    enddo
  enddo
!
! ************** Processing of transmission coefficients ***************
!
! Transmission coefficients averaged over spin and determination of maximal l-value.
! ECIS stops its output of transmission coefficients somewhat too early.
! For the highest l values the transmission coefficient for (l+spin) is not written in the output.
! Since these are small numbers we put them equal to the value for (l-spin).
!
  do type = 1, 6
    if (parskip(type)) cycle
!
! 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
!
    if (type /= 3 .and. type /= 6) then
Loop1:      do nen = ebegin(type), eend(type)
        do l = 0, numl
          if (Tjl(type, nen, - 1, l) /= 0 .and. Tjl(type, nen, 1, l) == 0) &
            Tjl(type, nen, 1, l) = Tjl(type, nen, - 1, l)
          if (Tjl(type, nen, - 1, l) == 0 .and. Tjl(type, nen, 1, l) /= 0 &
            .and. l > 0) Tjl(type, nen, - 1, l) = Tjl(type, nen, 1, l)
          Tl(type, nen, l) = ((l + 1) * Tjl(type, nen, 1, l) + l * Tjl(type, nen, - 1, l)) / (2 * l + 1)
          teps = Tl(type, nen, 0) * translimit / (2 * l + 1)
          teps = max(teps, transeps)
          if (Tjl(type, nen, - 1, l) < teps .and. Tjl(type, nen, 1, l) < teps) then
            lmax(type, nen) = l - 1
            cycle Loop1
          endif
          lmax(type, nen) = l
        enddo
      enddo Loop1
    endif
!
! 2. Spin 1 particles: Deuterons
!
    if (type == 3) then
Loop2:      do nen = ebegin(type), eend(type)
        do l = 0, numl
          if (Tjl(type, nen, - 1, l) /= 0 .and. Tjl(type, nen, 0, l) == 0) Tjl(type, nen, 0, l) = Tjl(type, nen, - 1, l)
          if (Tjl(type, nen, - 1, l) /= 0 .and. Tjl(type, nen, 1, l) == 0) Tjl(type, nen, 1, l) = Tjl(type, nen, - 1, l)
          if (Tjl(type, nen, - 1, l) == 0 .and. Tjl(type, nen, 1, l) /= 0 .and. l > 0) Tjl(type, nen, - 1, l) = Tjl(type, nen, 1, l)
          Tl(type, nen, l) = ((2 * l + 3) * Tjl(type, nen, 1, l) + (2 * l + 1) * Tjl(type, nen, 0, l) + &
 &          (2 * l - 1) * Tjl(type, nen, - 1, l)) / (3 * (2 * l + 1))
          teps = Tl(type, nen, 0) * translimit / (2 * l + 1)
          teps = max(teps, transeps)
          if (Tjl(type, nen, - 1, l) < teps .and. Tjl(type, nen, 0, l) < teps .and. Tjl(type, nen, 1, l) < teps) then
            lmax(type, nen) = l - 1
            cycle Loop2
          endif
          lmax(type, nen) = l
        enddo
      enddo Loop2
    endif
!
! 3. Spin 0 particles: Alpha-particles
!
    if (type == 6) then
Loop3:      do nen = ebegin(type), eend(type)
        do l = 0, numl
          Tl(type, nen, l) = Tjl(type, nen, 0, l)
          teps = Tl(type, nen, 0) * translimit / (2 * l + 1)
          teps = max(teps, transeps)
          if (Tl(type, nen, l) < teps) then
            lmax(type, nen) = l - 1
            cycle Loop3
          endif
          lmax(type, nen) = l
        enddo
      enddo Loop3
    endif
  enddo
  close (unit = 7, status = ecisstatus)
  open (unit = 10, file = 'ecis.invin', status = 'unknown')
  close (unit = 10, status = ecisstatus)
  return
end subroutine inverseread
! Copyright A.J. Koning 2021
