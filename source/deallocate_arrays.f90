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
  if (allocated(ENHratio)) deallocate(ENHratio)
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
  if (allocated(rhofis)) deallocate(rhofis)
  if (allocated(efistrrot)) deallocate(efistrrot)
  if (allocated(jfistrrot)) deallocate(jfistrrot)
  if (allocated(pfistrrot)) deallocate(pfistrrot)

  if (allocated(efisc2rot)) deallocate(efisc2rot)
  if (allocated(jfisc2rot)) deallocate(jfisc2rot)
  if (allocated(pfisc2rot)) deallocate(pfisc2rot)
  if (allocated(efistrhb)) deallocate(efistrhb)
  if (allocated(jfistrhb)) deallocate(jfistrhb)
  if (allocated(pfistrhb)) deallocate(pfistrhb)

  if (allocated(efisc2hb)) deallocate(efisc2hb)
  if (allocated(jfisc2hb)) deallocate(jfisc2hb)
  if (allocated(pfisc2hb)) deallocate(pfisc2hb)
!
! Medical isotope production
!
  if (allocated(Nenrp)) deallocate(Nenrp)
  if (allocated(prate)) deallocate(prate)
  if (allocated(Erp)) deallocate(Erp)
  if (allocated(xsrp)) deallocate(xsrp)

  if (allocated(Tmaxactivity)) deallocate(Tmaxactivity)
  if (allocated(Tp)) deallocate(Tp)
  if (allocated(Tgrid)) deallocate(Tgrid)

  if (allocated(Niso)) deallocate(Niso)
  if (allocated(activity)) deallocate(activity)
  if (allocated(yield)) deallocate(yield)
  if (allocated(Nisorel)) deallocate(Nisorel)
  if (allocated(Nisotot)) deallocate(Nisotot)
  if (allocated(Tmax)) deallocate(Tmax)
  if (allocated(popdecay)) deallocate(popdecay)
  if (allocated(preeqpopex)) deallocate(preeqpopex)

  if (allocated(xspopex)) deallocate(xspopex)
  if (allocated(xspopexP)) deallocate(xspopexP)

  if (allocated(maxJ)) deallocate(maxJ)
  if (allocated(deltaEx)) deallocate(deltaEx)
  if (allocated(Ex)) deallocate(Ex)

  if (allocated(Dmulti)) deallocate(Dmulti)
  if (allocated(fisfeedex)) deallocate(fisfeedex)
  if (allocated(mcontrib)) deallocate(mcontrib)
  if (allocated(mpecontrib)) deallocate(mpecontrib)
  if (allocated(popexcl)) deallocate(popexcl)
  if (allocated(xsbinspec)) deallocate(xsbinspec)
  if (allocated(xsmpe)) deallocate(xsmpe)
  if (allocated(xspartial)) deallocate(xspartial)
  if (allocated(Emsd)) deallocate(Emsd)
  if (allocated(nangleint)) deallocate(nangleint)

  if (allocated(xsdw)) deallocate(xsdw)
  if (allocated(xsdwin)) deallocate(xsdwin)

  if (allocated(xscont1)) deallocate(xscont1)
  if (allocated(xscontad1)) deallocate(xscontad1)

  if (allocated(msdstep1)) deallocate(msdstep1)
  if (allocated(msdstepad1)) deallocate(msdstepad1)

  if (allocated(msdstep)) deallocate(msdstep)
  if (allocated(msdstepad)) deallocate(msdstepad)

  if (allocated(xscont)) deallocate(xscont)
  if (allocated(xscontad)) deallocate(xscontad)
  if (allocated(msdstep0)) deallocate(msdstep0)
  if (allocated(msdstepad0)) deallocate(msdstepad0)

  if (allocated(msdstepint)) deallocate(msdstepint)
  if (allocated(msdstepintad)) deallocate(msdstepintad)
  if (allocated(msdsum)) deallocate(msdsum)
  if (allocated(msdtot)) deallocate(msdtot)
  if (allocated(msdtotad)) deallocate(msdtotad)
  if (allocated(msdtotintad)) deallocate(msdtotintad)
!
  return
end subroutine deallocate_arrays
