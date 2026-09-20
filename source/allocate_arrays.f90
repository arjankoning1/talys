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
! Emission rates are required for both particle-hole state density models.
  if (flag2comp) then
    allocate(wemission2(0:numpar,0:numparx,0:numparx,0:numen))
    wemission2 = 0.
  else
    allocate(wemission(0:numpar,0:numparx,0:numen))
    wemission = 0.
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
  if (flagracap) then
    allocate(phdensjp(0:numZ,0:numN,0:numdens,0:numJph,-1:1))
  endif
  if (flagfission) then
    allocate(rhofis(1:numbinfis,0:numJ,-1:1,1:numbar))
    rhofis = 0.d0
!
! Head-band transition states
!
    if (flaghbstate) then
      allocate(efistrhb(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numlev))
      allocate(jfistrhb(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numlev))
      allocate(pfistrhb(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numlev))

      efistrhb = 0.
      jfistrhb = 0.
      pfistrhb = 1
    endif
!
! Class-2 transition states
!
    if (flagclass2) then
      allocate(efisc2hb(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numlev))
      allocate(jfisc2hb(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numlev))
      allocate(pfisc2hb(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numlev))

      efisc2hb = 0.
      jfisc2hb = 0.
      pfisc2hb = 1
    endif
!
    allocate(efistrrot(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numrot))
    allocate(jfistrrot(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numrot))
    allocate(pfistrrot(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numrot))

    efistrrot = 0.
    jfistrrot = 0.
    pfistrrot = 1

    if (flagclass2) then
      allocate(efisc2rot(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numrot))
      allocate(jfisc2rot(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numrot))
      allocate(pfisc2rot(0:min(numZ,maxZ+2),0:min(numN,maxN+2),1:numbar,0:numrot))

      efisc2rot = 0.
      jfisc2rot = 0.
      pfisc2rot = 1
    endif
  endif
  return
end subroutine allocate_arrays
