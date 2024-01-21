subroutine fissionparout(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output for fission parameters
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
! Variables for fission
!   axtype           ! type of axiality of barrier
!   betafiscor       ! adjustable factor for fission path width
!   betafiscoradjust ! adjustable factor for fission path width
!   fbarrier         ! height of fission barrier
!   fismodelx        ! fission model
!   flagclass2       ! flag for class2 states in fission
!   fwidth           ! width of fission barrier
!   vfiscor          ! adjustable factor for fission path height
!   vfiscoradjust    ! adjustable factor for fission path height
!   widthc2          ! width of class2 states
! Variables for level density
!   Rclass2mom    ! norm. constant for moment of inertia for class 2 states
!   Rtransmom     ! norm. constant for moment of inertia for transition states
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   NN            ! neutron number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Constants
!   cparity       ! parity (character)
!   nuc           ! symbol of nucleus
! Variables for fission parameters
!   efisc2hb      ! energy of class2 states
!   efisc2rot     ! energy of class2 rotational transition states
!   efistrhb      ! energy of head band transition states
!   efistrrot     ! energy of rotational transition states
!   fecont        ! start of continuum energy
!   jfisc2hb      ! spin of class2 states
!   jfisc2rot     ! spin of class2 rotational transition states
!   jfistrhb      ! spin of head band transition states
!   jfistrrot     ! spin of rotational transition states
!   minertc2      ! moment of inertia for class2 states
!   minertia      ! moment of inertia of fission barrier deformation
!   nclass2       ! number of sets of class2 states
!   nfisbar       ! number of fission barrier parameters
!   nfisc2hb      ! number of class2 states for barrier
!   nfisc2rot     ! number of class2 rotational transition states for barrier
!   nfistrhb      ! number of head band transition states for barrier
!   nfistrrot     ! number of rotational transition states for barrier
!   pfisc2hb      ! parity of class2 states
!   pfisc2rot     ! parity of class2 rotational transition states
!   pfistrhb      ! parity of head band transition states
!   pfistrrot     ! parity of rotational transition states
!
! *** Declaration of local data
!
  implicit none
  integer :: A   ! mass number of target nucleus
  integer :: i   ! counter
  integer :: j   ! counter
  integer :: N   ! neutron number of residual nucleus
  integer :: Nix ! neutron number index for residual nucleus
  integer :: Z   ! charge number of target nucleus
  integer :: Zix ! charge number index for residual nucleus
!
! ****************** Output of fission barrier parameters **************
!
! 1. Main fission parameters
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  write(*, '(/" Fission information for Z=", i3, " N=", i3, " (", i3, a2, ") "/)')  Z, N, A, nuc(Z)
  write(*, '(" Number of fission barriers           :", i3)') nfisbar(Zix, Nix)
  if (flagclass2) write(*, '(" Number of sets of class2 states      :", i3)') nclass2(Zix, Nix)
  if (fismodelx(Zix, Nix) == 5) then
    write(*, '(" Correction factor betafiscor:", f8.3)') betafiscor(Zix, Nix)
    write(*, '(" Correction factor vfiscor   :", f8.3)') vfiscor(Zix, Nix)
    write(*, '(" Adjustable factor betafiscoradjust:", f8.3)') betafiscoradjust(Zix, Nix)
    write(*, '(" Adjustable factor vfiscoradjust   :", f8.3)') vfiscoradjust(Zix, Nix)
  endif
  do i = 1, nfisbar(Zix, Nix)
    write(*, '(/" Parameters for fission barrier", i3/)') i
      write(*, '(" Type of axiality                     :", i3)') axtype(Zix, Nix, i)
    write(*, '(" Height of fission barrier ", i1, "          :", f8.3)') i, fbarrier(Zix, Nix, i)
    write(*, '(" Width of fission barrier ", i1, "           :", f8.3)') i, fwidth(Zix, Nix, i)
    write(*, '(" Rtransmom                            :", f8.3)') Rtransmom(Zix, Nix, i)
    write(*, '(" Moment of inertia                    :", f8.3)') minertia(Zix, Nix, i)
    write(*, '(" Number of head band transition states:", i3)') nfistrhb(Zix, Nix, i)
    write(*, '(" Start of continuum energy            :", f8.3)') fecont(Zix, Nix, i)
!
! 2. Head band transition states
!
    write(*, '(/" Head band transition states"/)')
    write(*, '("  no.    E    spin    parity"/)')
    do j = 1, nfistrhb(Zix, Nix, i)
      write(*, '(1x, i4, f8.3, f6.1, 3x, a1)') j, efistrhb(Zix, Nix, i, j), jfistrhb(Zix, Nix, i, j), &
 &      cparity(pfistrhb(Zix, Nix, i, j))
    enddo
!
! 3. Rotational bands transition states
!
    write(*, '(/" Rotational bands"/)')
    write(*, '("  no.    E    spin    parity"/)')
    do j = 1, nfistrrot(Zix, Nix, i)
      write(*, '(1x, i4, f8.3, f6.1, 3x, a1)') j, efistrrot(Zix, Nix, i, j), jfistrrot(Zix, Nix, i, j), &
 &      cparity(pfistrrot(Zix, Nix, i, j))
    enddo
  enddo
!
! 4. Class2 states
!
  if (flagclass2) then
    do i = 1, nclass2(Zix, Nix)
      write(*, '(/" Parameters for set", i3, " of class2 states"/)') i
      write(*, '(" Rclass2mom                         :", f8.3)') Rclass2mom(Zix, Nix, i)
      write(*, '(" Moment of inertia                  :", f8.3)') minertc2(Zix, Nix, i)
      write(*, '(" Number of class2 states            :", i3)') nfisc2hb(Zix, Nix, i)
      write(*, '(" Width of class2 states (MeV)       :", f8.3)') widthc2(Zix, Nix, i)
      write(*, '(/" Class 2 states"/)')
      write(*, '("  no.    E    spin    parity"/)')
      do j = 1, nfisc2hb(Zix, Nix, i)
        write(*, '(1x, i4, f8.3, f6.1, 3x, a1)') j, efisc2hb(Zix, Nix, i, j), &
 &        jfisc2hb(Zix, Nix, i, j), cparity(pfisc2hb(Zix, Nix, i, j))
      enddo
!
! 5. Rotational bands
!
      write(*, '(/" Rotational bands"/)')
      write(*, '("  no.    E    spin    parity"/)')
      do j = 1, nfisc2rot(Zix, Nix, i)
        write(*, '(1x, i4, f8.3, f6.1, 3x, a1)') j, efisc2rot(Zix, Nix, i, j), &
 &        jfisc2rot(Zix, Nix, i, j), cparity(pfisc2rot(Zix, Nix, i, j))
      enddo
    enddo
  endif
  return
end subroutine fissionparout
! Copyright A.J. Koning 2021
