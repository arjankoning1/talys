subroutine rotband(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Build rotational bands on transition states
!
! Author    : Stephane Hilaire
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
!   numrot       ! number of rotational states
! Variables for fission parameters
!   efistrhb     ! energy of head band transition states
!   efistrrot    ! energy of rotational transition states
!   fecont       ! start of continuum energy
!   jfistrhb     ! spin of head band transition states
!   jfistrrot    ! spin of rotational transition states
!   minertia     ! moment of inertia of fission barrier deformation
!   nfisbar      ! number of fission barrier parameters
!   nfistrhb     ! number of head band transition states for barrier
!   nfistrrot    ! number of rotational transition states for barrier
!   pfistrhb     ! parity of head band transition states
!   pfistrrot    ! parity of rotational transition states
!
! *** Declaration of local data
!
  implicit none
  integer   :: i      ! level
  integer   :: its    ! counter
  integer   :: itstot ! help variable
  integer   :: jts    ! counter
  integer   :: nbi    ! counter
  integer   :: Nix    ! neutron number index for residual nucleus
  integer   :: pa1    ! help variable
  integer   :: pa2    ! help variable
  integer   :: Zix    ! charge number index for residual nucleus
  real(sgl) :: Eband  ! help variable
  real(sgl) :: en1    ! start energy of local adjustment
  real(sgl) :: en2    ! end energy of local adjustment
  real(sgl) :: Erk10  ! energy shift between k=0- and k=1-
  real(sgl) :: Erot   ! rotational energy
  real(sgl) :: jstart ! help variable
  real(sgl) :: jstep  ! help variable
  real(sgl) :: rj     ! help variable
  real(sgl) :: rj1    ! help variable
  real(sgl) :: rj2    ! help variable
!
! ******************* Rotational bands building ************************
!
!
Loop1:  do nbi = 1, nfisbar(Zix, Nix)
    nfistrrot(Zix, Nix, nbi) = 0
Loop2: do i = 1, nfistrhb(Zix, Nix, nbi)
      jstart = jfistrhb(Zix, Nix, nbi, i)
      jstep = 1.
      Erk10 = 0.
      if (jfistrhb(Zix, Nix, nbi, i) == 0) jstep = 2.
      if ((jfistrhb(Zix, Nix, nbi, i) == 0) .and. (pfistrhb(Zix, Nix, nbi, i) ==  - 1)) then
        jstart = 1.
        Erk10 = 1. / minertia(Zix, Nix, nbi)
      endif
      rj = jstart - jstep
      do
        rj = rj + jstep
        Erot = (rj * (rj + 1.) - jstart * (jstart + 1.)) / (2. * minertia(Zix, Nix, nbi))
        Eband = efistrhb(Zix, Nix, nbi, i) + Erot + Erk10
        if (Eband > fecont(Zix, Nix, nbi)) cycle Loop2
        nfistrrot(Zix, Nix, nbi) = nfistrrot(Zix, Nix, nbi) + 1
        itstot = nfistrrot(Zix, Nix, nbi)
        efistrrot(Zix, Nix, nbi, itstot) = Eband
        jfistrrot(Zix, Nix, nbi, itstot) = rj
        pfistrrot(Zix, Nix, nbi, itstot) = pfistrhb(Zix, Nix, nbi, i)
        if (itstot > numrot) cycle Loop1
      enddo
    enddo Loop2
  enddo Loop1
!
! *** Reorganize transition states by increasing excitation energy *****
!
  do nbi = 1, nfisbar(Zix, Nix)
    do its = 1, nfistrrot(Zix, Nix, nbi) - 1
      do jts = its + 1, nfistrrot(Zix, Nix, nbi)
        en1 = efistrrot(Zix, Nix, nbi, its)
        rj1 = jfistrrot(Zix, Nix, nbi, its)
        pa1 = pfistrrot(Zix, Nix, nbi, its)
        en2 = efistrrot(Zix, Nix, nbi, jts)
        rj2 = jfistrrot(Zix, Nix, nbi, jts)
        pa2 = pfistrrot(Zix, Nix, nbi, jts)
        if (en1 > en2) then
          efistrrot(Zix, Nix, nbi, its) = en2
          jfistrrot(Zix, Nix, nbi, its) = rj2
          pfistrrot(Zix, Nix, nbi, its) = pa2
          efistrrot(Zix, Nix, nbi, jts) = en1
          jfistrrot(Zix, Nix, nbi, jts) = rj1
          pfistrrot(Zix, Nix, nbi, jts) = pa1
        endif
      enddo
    enddo
  enddo
  return
end subroutine rotband
! Copyright A.J. Koning 2021
