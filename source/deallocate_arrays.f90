subroutine deallocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_talys_mod
  implicit none
  integer :: Nix
  integer :: Zix
!
  if (allocated(Erescue)) deallocate(Erescue)
  if (allocated(frescue)) deallocate(frescue)
  if (allocated(ddxrec)) deallocate(ddxrec)
  if (allocated(phtable1)) deallocate(phtable1)
  if (allocated(phtable2)) deallocate(phtable2)
! if (allocated(fqrpa)) deallocate(fqrpa)
! if (allocated(fqrpaJP)) deallocate(fqrpaJP)
  do Nix = 0, numN
    do Zix = 0, numZ
      if (allocated(qrpa(Zix,Nix)%e)) deallocate(qrpa(Zix,Nix)%e)
      if (allocated(qrpa(Zix,Nix)%f)) deallocate(qrpa(Zix,Nix)%f)
      if (allocated(qrpa(Zix,Nix)%fJP)) deallocate(qrpa(Zix,Nix)%fJP)
    enddo
  enddo
!
  return
end subroutine deallocate_arrays
