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
! ddxrec is also filled for the residual reached after charged-particle
! emission.  Zindex/Nindex can then be two indices beyond the parent
! compound-nucleus range (alpha emission), so retain those boundary bins.
    allocate(ddxrec(0:maxZ+2,0:maxN+2,0:numex,0:maxenrec,0:numangrec))
    ddxrec = 0.
  endif
  if (phmodel == 2) then
    allocate(phtable1(0:1, 0:1, 0:numexc, 0:numexc, 0:numdens))
    allocate(phtable2(0:1, 0:1, 0:numexc, 0:numexc, 0:numexc, 0:numexc, 0:numdens))
    phtable1 = 0.
    phtable2 = 0.
  endif
  if (flagchannels) then
    allocate(feedexcl(0:min(maxZ,numZchan),0:min(maxN,numNchan),0:numpar,0:numex+1,0:numex+1))
    feedexcl = 0.
  endif
  if (flagracap) then
    allocate(phdensjp(0:numZ,0:numN,0:numdens,0:numJph,-1:1))
  endif
  return
end subroutine allocate_arrays
