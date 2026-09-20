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
  if (allocated(ldtable)) deallocate(ldtable)
  if (allocated(ldtableT)) deallocate(ldtableT)
  if (allocated(ldtableN)) deallocate(ldtableN)
  if (allocated(ldtottable)) deallocate(ldtottable)
  if (allocated(ldtottableP)) deallocate(ldtottableP)
  do Nix = 0, numN
    do Zix = 0, numZ
      if (allocated(qrpa(Zix,Nix)%e)) deallocate(qrpa(Zix,Nix)%e)
      if (allocated(qrpa(Zix,Nix)%f)) deallocate(qrpa(Zix,Nix)%f)
      if (allocated(qrpa(Zix,Nix)%fJP)) deallocate(qrpa(Zix,Nix)%fJP)
    enddo
  enddo
  if (allocated(phexist1)) deallocate(phexist1)
  if (allocated(phexist2)) deallocate(phexist2)
  if (allocated(xspopph)) deallocate(xspopph)
  if (allocated(xspopph2)) deallocate(xspopph2)
  if (allocated(wemission)) deallocate(wemission)
  if (allocated(wemission2)) deallocate(wemission2)
  if (allocated(feedexcl)) deallocate(feedexcl)
  if (allocated(phdensjp)) deallocate(phdensjp)
  if (allocated(fxsgamdischan)) deallocate(fxsgamdischan)
  if (allocated(fxsgamchannel)) deallocate(fxsgamchannel)
  if (allocated(fisfeedJP)) deallocate(fisfeedJP)
  if (allocated(xspop)) deallocate(xspop)
  if (allocated(rhogrid)) deallocate(rhogrid)
!
  return
end subroutine deallocate_arrays
