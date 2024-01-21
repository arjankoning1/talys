subroutine totalrecoil
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Recoil results
!
! Author    : Arjan Koning and Stephane Hilaire
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! Variables for numerics
!   maxenrec       ! number of recoil energies
!   maxN           ! maximal number of neutrons away from initial compound nucleus
!   maxZ           ! maximal number of protons away from initial compound nucleus
!   nanglecont     ! number of angles for continuum
!   nanglerec      ! number of recoil angles
!   xseps          ! limit for cross sections
! Variables for basic reaction
!   flaglabddx     ! flag for calculation of DDX in LAB system
! Variables for incident channel
!   xspopnuc       ! population cross section per nucleus
! Variables for nuclides
!   parskip        ! logical to skip outgoing particle
! Constants
!   twopi          ! 2 * pi
! Variables for level density
!   Nlast          ! last discrete level
! Variables for recoil initialization
!   dcosangcont    ! width of cosine bin
! Variables for recoil
!   areareclab     ! Total surface of LAB ddx bins
!   ddxejlab       ! array containing the ddx spectrum of light part
!   ddxrec         ! array containing the lab double differential xs
!   dEejlab        ! width of ejectile lab bin
!   Erecmax        ! minimal energy limit of recoil bin
!   Erecmin        ! minimal energy limit of recoil bin
!   iejlab         ! number of ejectile lab bins
!   recoilint      ! total recoil integrated over spectrum
!   specrecoil     ! recoil spectrum
!   xsejlab        ! LAB ejectile spectrum
!   xsejlabint     ! LAB energy - integrated spectrum
!
! *** Declaration of local data
!
  implicit none
  integer   :: iang     ! running variable for angle
  integer   :: ixl      ! help variable
  integer   :: iyl      ! help variable
  integer   :: Ncomp    ! neutron number index for compound nucleus
  integer   :: nen      ! energy counter
  integer   :: nexrec   ! energy index for recoils
  integer   :: type     ! particle type
  integer   :: Zcomp    ! proton number index for compound nucleus
  real(sgl) :: dErecoil ! width of recoil bin
  real(sgl) :: sum      ! help variable
  real(sgl) :: sumen    ! help variable
!
! ******************** Angle-integrated spectra ************************
!
  if (flaglabddx) then
    do type = 0, 6
      if (parskip(type)) cycle
      sumen = 0.
      do nen = 1, iejlab(type)
        sum = 0.
        do iang = 0, nanglecont
          sum = sum + ddxejlab(type, nen, iang) * dcosangcont(iang)
        enddo
        xsejlab(type, nen) = sum * twopi
        sumen = sumen + sum * twopi * dEejlab(type, nen)
      enddo
      xsejlabint(type) = sumen
    enddo
  endif
!
! ************************* Recoils in LAB system **********************
!
  do Zcomp = 0, maxZ
    do Ncomp = 0, maxN
      if (xspopnuc(Zcomp, Ncomp) < xseps) cycle
      sumen = 0.
      do ixl = 0, maxenrec
        dErecoil = abs(Erecmax(Zcomp, Ncomp, ixl) - Erecmin(Zcomp, Ncomp, ixl))
        if (dErecoil <= 1.0e-14) dErecoil = 1.
        sum = 0.
        do nexrec = 0, max(Nlast(Zcomp, Ncomp, 0), 1)
          do iyl = 0, nanglerec
            sum = sum + ddxrec(Zcomp, Ncomp, nexrec, ixl, iyl) * areareclab(Zcomp, Ncomp, ixl, iyl)
          enddo
        enddo
        specrecoil(Zcomp, Ncomp, ixl) = sum * twopi / dErecoil
        sumen = sumen + sum * twopi
      enddo
      recoilint(Zcomp, Ncomp) = sumen
    enddo
  enddo
  return
end subroutine totalrecoil
! Copyright A.J. Koning 2021
