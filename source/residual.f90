subroutine residual
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Residual production cross sections
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
! All global variables
!   numNastro      ! maximal number of neutrons away from initial CN for astroph. calcs
!   numZastro      ! maximal number of protons away from initial CN for astroph. calcs
! Variables for numerics
!   maxN           ! maximal number of neutrons away from initial compound nucleus
!   maxZ           ! maximal number of protons away from initial compound nucleus
!   xseps          ! limit for cross sections
! Variables for input energies
!   nin            ! counter for incident energy
! Variables for basic reaction
!   flagastro      ! flag for calculation of astrophysics reaction rate
! Variables for OMP
!   flagomponly    ! flag to execute ONLY an optical model calculation
! Variables for energies
!   Qres           ! Q - value for residual nucleus
! Variables for incident channel
!   maxA           ! maximal number of nucleons away from initial compound nucleus
!   xsbranch       ! branching ratio for isomeric cross section
!   xsmassprod     ! residual production cross section per mass unit
!   xsresprod      ! total residual production ( = reaction) cross section
!   xspopex        ! population cross section summed over spin and parity
!   xspopnuc       ! population cross section per nucleus
! Variables for levels
!   tau            ! lifetime of state in seconds
! Variables for level density
!   Nlast          ! last discrete level
! Variables for astro
!   xsastro        ! cross section for astrophysical calculatio
!   xsastroex      ! cross section for astrophysical calculati
!
! *** Declaration of local data
!
  implicit none
  integer :: Acomp  ! mass number index for compound nucleus
  integer :: Ncomp  ! neutron number index for compound nucleus
  integer :: nex    ! excitation energy bin of compound nucleus
  integer :: Zcomp  ! proton number index for compound nucleus
!
! ************************ Cross sections ******************************
!
  if (flagomponly) return
  xsresprod = 0.
  do Acomp = 0, maxA
    xsmassprod(Acomp) = 0.
    do Zcomp = 0, maxZ
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxN) cycle
      if (xspopnuc(Zcomp, Ncomp) /= 0.) then
        xsresprod = xsresprod + xspopnuc(Zcomp, Ncomp)
        xsmassprod(Acomp) = xsmassprod(Acomp) + xspopnuc(Zcomp, Ncomp)
        do nex = 0, Nlast(Zcomp, Ncomp, 0)
          if (nex == 0 .or. tau(Zcomp, Ncomp, nex) /= 0.) &
 &          xsbranch(Zcomp, Ncomp, nex) = xspopex(Zcomp, Ncomp, nex) / xspopnuc(Zcomp, Ncomp)
        enddo
      endif
!
! For non-threshold reactions (positive Q-value) we always assign a minimum value to the exclusive cross section.
! (The transmission coefficients for these reactions might have been zero (from ECIS), but non-threshold reactions
! theoretically have a non-zero cross section.)
!
       if (Qres(Zcomp, Ncomp, 0) > 0..and. xspopnuc(Zcomp, Ncomp) <= xseps) xspopnuc(Zcomp, Ncomp) = xseps
       if (flagastro .and. Zcomp <= numZastro .and. Ncomp <= numNastro) then
         xsastro(Zcomp, Ncomp, nin) = xspopnuc(Zcomp, Ncomp)
         do nex = 0, Nlast(Zcomp, Ncomp, 0)
           if (nex == 0 .or. tau(Zcomp, Ncomp, nex) /= 0.) &
             xsastroex(Zcomp, Ncomp, nin, nex) = xspopex(Zcomp, Ncomp, nex)
         enddo
       endif
    enddo
  enddo
  return
end subroutine residual
! Copyright A.J. Koning 2021
