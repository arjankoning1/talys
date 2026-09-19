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
    allocate(ddxrec(0:maxZ+2,0:maxN+2,0:numex,0:maxenrec,0:nanglerec))
    ddxrec = 0.
  endif
  if (phmodel == 2) then
    if (flag2comp) then
      allocate(phexist2(0:numZ,0:numN,0:numexc,0:numexc,0:numexc,0:numexc))
      allocate(phtable2(0:1, 0:1, 0:numexc, 0:numexc, 0:numexc, 0:numexc, 0:numdens))
      phexist2 = .false.
      phtable2 = 0.
    else
      allocate(phexist1(0:numZ,0:numN,0:numexc,0:numexc))
      allocate(phtable1(0:1, 0:1, 0:numexc, 0:numexc, 0:numdens))
      phexist1 = .false.
      phtable1 = 0.
    endif
  endif
  if (flagfission) then
    allocate(fisfeedJP(0:maxZ,0:maxN,0:numex+1,0:numJ,-1:1))
    fisfeedJP = 0.
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
