subroutine allocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_talys_mod
  implicit none
!
  if (flagrescue) then
    allocate(Erescue(nummt,-1:numisom,numen6))
    allocate(frescue(nummt,-1:numisom,numen6))
    Erescue = 0.
    frescue = 0.
  endif
  if (flagrecoil) then
    allocate(ddxrec(0:maxZ,0:maxN,0:numex,0:maxenrec,0:numangrec))
    ddxrec = 0.
  endif
  if (phmodel == 2) then
    allocate(phtable1(0:1, 0:1, 0:numexc, 0:numexc, 0:numdens))
    allocate(phtable2(0:1, 0:1, 0:numexc, 0:numexc, 0:numexc, 0:numexc, 0:numdens))
    phtable1 = 0.
    phtable2 = 0.
  endif
! allocate(fqrpa(0:maxZ+2, 0:maxN+2, 0:numgamqrpa, numTqrpa, 0:1, numgam))
! fqrpa = 0.
! allocate(fqrpaJP(0:min(maxZ+2,numZph), 0:min(maxN+2,numNph), 0:numgamqrpa, &
!&  numTqrpa, 0:1, 0:9, 0:1))
! fqrpaJP = 0.
!
  return
end subroutine allocate_arrays
