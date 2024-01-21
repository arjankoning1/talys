subroutine rotclass2(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Build rotational bands on class2 states
!
! Author    : Stephane Hilaire and Pascal Romain
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
! All global variables
!   numrot        ! number of rotational states
! Variables for fission
!   fbaradjust    ! adjustable factor for fission parameters
!   fbarrier      ! height of fission barrier
!   widthc2       ! width of class2 states
! Variables for fission parameters
!   efisc2hb      ! energy of class2 states
!   efisc2rot     ! energy of class2 rotational transition states
!   Emaxclass2    ! maximum energy for class2 states
!   fecont        ! start of continuum energy
!   jfisc2hb      ! spin of class2 states
!   jfisc2rot     ! spin of class2 rotational transition states
!   minertc2      ! moment of inertia for class2 states
!   nclass2       ! number of sets of class2 states
!   nfisbar       ! number of fission barrier parameters
!   nfisc2hb      ! number of class2 states for barrier
!   nfisc2rot     ! number of class2 rotational transition states for barrier
!   pfisc2hb      ! parity of class2 states
!   pfisc2rot     ! parity of class2 rotational transition states
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                 ! counter
  integer   :: its               ! counter
  integer   :: itstot            ! help variable
  integer   :: jts               ! counter
  integer   :: nbi               ! counter
  integer   :: Nix               ! neutron number index for residual nucleus
  integer   :: pa1               ! help variable
  integer   :: pa2               ! help variable
  integer   :: Zix               ! charge number index for residual nucleus
  real(sgl) :: Eband             ! help variable
  real(sgl) :: Ecut              ! help variable
  real(sgl) :: en1               ! start energy of local adjustment
  real(sgl) :: en2               ! end energy of local adjustment
  real(sgl) :: Erk10             ! energy shift between k=0- and k=1-
  real(sgl) :: Erot              ! rotational energy
  real(sgl) :: jstart            ! help variable
  real(sgl) :: jstep             ! help variable
  real(sgl) :: rj                ! help variable
  real(sgl) :: rj1               ! help variable
  real(sgl) :: rj2               ! help variable
!
! ******************* Rotational bands building ************************
!
! Calculation of maximum energy of class2 states
!
  if (nfisbar(Zix, Nix) == 2) then
    Emaxclass2(Zix, Nix, 1) = fbaradjust(Zix, Nix, 1) * fbarrier(Zix, Nix, 1) + fecont(Zix, Nix, 1) + 0.5 * widthc2(Zix, Nix, 1)
  endif
  if (nfisbar(Zix, Nix) == 3) then
    Emaxclass2(Zix, Nix, 1) = fbaradjust(Zix, Nix, 1) * fbarrier(Zix, Nix, 1) + fecont(Zix, Nix, 1) + 0.5 * widthc2(Zix, Nix, 1)
    Emaxclass2(Zix, Nix, 2) = fbaradjust(Zix, Nix, 2) * fbarrier(Zix, Nix, 2) + fecont(Zix, Nix, 2) + 0.5 * widthc2(Zix, Nix, 2)
  endif
!
! Rotational bands construction
!
Loop1:  do nbi = 1, nclass2(Zix, Nix)
    nfisc2rot(Zix, Nix, nbi) = 0
    Ecut = Emaxclass2(Zix, Nix, nbi)
    Loop2: do i = 1, nfisc2hb(Zix, Nix, nbi)
      jstart = jfisc2hb(Zix, Nix, nbi, i)
      jstep = 1.
      Erk10 = 0.
      if (jfisc2hb(Zix, Nix, nbi, i) == 0) jstep = 2.
      if ((jfisc2hb(Zix, Nix, nbi, i) == 0) .and. (pfisc2hb(Zix, Nix, nbi, i) ==  - 1)) then
        jstart = 1.
        Erk10 = 1. / minertc2(Zix, Nix, nbi)
      endif
      rj = jstart - jstep
      do
        rj = rj + jstep
        Erot = (rj * (rj + 1.) - jstart * (jstart + 1.)) / (2. * minertc2(Zix, Nix, nbi))
        Eband = efisc2hb(Zix, Nix, nbi, i) + Erot + Erk10
        if (Eband > Ecut) cycle Loop2
        nfisc2rot(Zix, Nix, nbi) = nfisc2rot(Zix, Nix, nbi) + 1
        itstot = nfisc2rot(Zix, Nix, nbi)
        efisc2rot(Zix, Nix, nbi, itstot) = Eband
        jfisc2rot(Zix, Nix, nbi, itstot) = rj
        pfisc2rot(Zix, Nix, nbi, itstot) = pfisc2hb(Zix, Nix, nbi, i)
        if (itstot == numrot) cycle Loop1
      enddo
    enddo Loop2
  enddo Loop1
!
! ****** Reorganize class2 states by increasing excitation energy ******
!
  do nbi = 1, nclass2(Zix, Nix)
    do its = 1, nfisc2rot(Zix, Nix, nbi) - 1
      do jts = its + 1, nfisc2rot(Zix, Nix, nbi)
        en1 = efisc2rot(Zix, Nix, nbi, its)
        rj1 = jfisc2rot(Zix, Nix, nbi, its)
        pa1 = pfisc2rot(Zix, Nix, nbi, its)
        en2 = efisc2rot(Zix, Nix, nbi, jts)
        rj2 = jfisc2rot(Zix, Nix, nbi, jts)
        pa2 = pfisc2rot(Zix, Nix, nbi, jts)
        if (en1 > en2) then
          efisc2rot(Zix, Nix, nbi, its) = en2
          jfisc2rot(Zix, Nix, nbi, its) = rj2
          pfisc2rot(Zix, Nix, nbi, its) = pa2
          efisc2rot(Zix, Nix, nbi, jts) = en1
          jfisc2rot(Zix, Nix, nbi, jts) = rj1
          pfisc2rot(Zix, Nix, nbi, jts) = pa1
        endif
      enddo
    enddo
  enddo
  if (nfisbar(Zix, Nix) == 2) Emaxclass2(Zix, Nix, 1) = Emaxclass2(Zix, Nix, 1) - 0.5 * widthc2(Zix, Nix, 1)
  if (nfisbar(Zix, Nix) == 3) then
    Emaxclass2(Zix, Nix, 1) = Emaxclass2(Zix, Nix, 1) - 0.5 * widthc2(Zix, Nix, 1)
    Emaxclass2(Zix, Nix, 2) = Emaxclass2(Zix, Nix, 2) - 0.5 * widthc2(Zix, Nix, 2)
  endif
  return
end subroutine rotclass2
! Copyright A.J. Koning 2021
