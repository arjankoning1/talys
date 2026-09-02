subroutine deallocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_talys_mod
  implicit none
!
  if (allocated(Erescue)) deallocate(Erescue)
  if (allocated(frescue)) deallocate(frescue)
  if (allocated(ddxrec)) deallocate(ddxrec)
  if (allocated(phtable1)) deallocate(phtable1)
  if (allocated(phtable2)) deallocate(phtable2)
  if (allocated(fqrpa)) deallocate(fqrpa)
  if (allocated(fqrpaJP)) deallocate(fqrpaJP)
!
  return
end subroutine deallocate_arrays
