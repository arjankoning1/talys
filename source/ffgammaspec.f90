subroutine ffgammaspec
!
! Add discrete gamma transitions from the current fission-fragment
! evaporation calculation to the prompt-fission gamma spectrum.
  use A0_talys_mod
  implicit none
  integer   :: Zix
  integer   :: Nix
  integer   :: i1
  integer   :: i2
  integer   :: nen
  real(sgl) :: Egam
  real(sgl) :: xsline
  real(sgl) :: wlow
  real(sgl) :: whigh
  real(sgl) :: dE
  do Zix = 0, maxZ
    do Nix = 0, maxN
      do i1 = 1, Nlast(Zix, Nix, 0)
        do i2 = 0, i1 - 1
          xsline = xsgamdis(Zix, Nix, i1, i2)
          if (xsline <= 0.) cycle
          Egam = edis(Zix, Nix, i1) - edis(Zix, Nix, i2)
          if (Egam < Epfns(1)) cycle
          if (Egam > Epfns(NEpfns)) cycle
!
! Locate the two surrounding PFGS grid points.
!
          call locate(Epfns, 0, NEpfns - 1, Egam, nen)
          nen = nen + 1
          nen = max(1, min(nen, NEpfns - 1))
          dE = Epfns(nen + 1) - Epfns(nen)
          if (dE > 0.) then
            whigh = (Egam - Epfns(nen)) / dE
            whigh = max(0., min(1., whigh))
          else
            whigh = 0.
          endif
          wlow = 1. - whigh
!
! xsgamdis is an integrated line cross section [mb].
! pfns is a differential spectrum [mb/MeV].
!
          pfns(0, nen) = pfns(0, nen) + &
            wlow * xsline / dEpfns(nen)
          pfns(0, nen + 1) = pfns(0, nen + 1) + &
            whigh * xsline / dEpfns(nen + 1)
!
! At present TALYS does not transform gamma spectra from the fragment
! frame to the laboratory frame, so use the same treatment here.
!
          pfnscm(0, nen) = pfnscm(0, nen) + &
            wlow * xsline / dEpfns(nen)
          pfnscm(0, nen + 1) = pfnscm(0, nen + 1) + &
            whigh * xsline / dEpfns(nen + 1)
        enddo
      enddo
    enddo
  enddo
  return
end subroutine ffgammaspec
