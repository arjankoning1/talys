subroutine normalization
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Normalize cross sections to experimental or evaluated data
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl              ! single precision kind
! All global variables
!   numen            ! maximum number of outgoing energies
!   numisom          ! number of isomers
!   numlev           ! maximum number of discrete levels
!   nummt            ! number of MT numbers
! Variables for energies
!   Ninclow        ! number of incident energies below Elow
! Variables for normalization
!   xseladjust       ! elastic cross section adjustment
!   xsnonadjust      ! nonelastic cross section adjustment
!   xstotadjust      ! total cross section adjustment
! Variables for output
!   flagcompo        ! flag for output of cross section components
!   flagspec         ! flag for output of spectra
! Variables for fission
!   flagfission      ! flag for fission
! Variables for input energies
!   nin              ! counter for incident energy
! Variables for total cross sections
!   xsexclcont       ! exclusive single channel cross section for continuum
!   xsexclusive      ! exclusive single channel cross section
!   xsfistot         ! total fission cross section
!   xsfistot0        ! total fission cross section
! Variables for energy grid
!   Crescue          ! adjustment factor for this incident energy
!   Einc             ! incident energy in MeV
!   Erescue          ! energy grid for adjustment factors
!   frescue          ! adjustment factor
!   Nrescue          ! number of energies for adjustment factors
! Variables for energies
!   idchannel        ! identifier for exclusive channel
! Variables for exclusive channels
!   idnum            ! counter for exclusive channel
!   xschaniso        ! channel cross section per isomer
!   xschannel        ! channel cross section
!   xschannelsp      ! channel cross section spectra
!   xsfischannel     ! fission channel cross section
!   xsgamchannel     ! gamma channel cross section
!   xsgamdischan     ! discrete gamma channel cross section
! Variables for multiple emission
!   xsngn            ! total (projectile, gamma - ejectile) cross section
!   xspopcomp        ! compound population cross section per nucleus
!   xspoppreeq       ! preequilibrium population cross section per nucleus
! Variables for binary reactions
!   xscompdisc       ! compound cross section for discrete state
!   xscompdisctot    ! compound cross section summed over discrete states
!   xsconttot        ! total cross section for continuum
!   xsdircont        ! direct cross section for continuum
!   xsdirect         ! total direct cross section
!   xsdisc           ! total cross section for discrete state
!   xsdisctot        ! total cross section summed over discrete states
!   xselastot        ! total elastic cross section (shape + compound)
!   xsnonel          ! non - elastic cross
!   xspopdir         ! direct population cross section per nucleus
! Variables for incident channel
!   channelsum       ! sum over exclusive channel cross sections
!   multiplicity     ! particle multiplicity
!   xsbinary         ! cross section from initial compound to residual nucleus
!   xscompcont       ! compound cross section for continuum
!   xsdirdisc        ! direct cross section for discrete state direct cross section
!   xsdirdisctot     ! direct cross section summed over discrete states
!   xsparticle       ! total particle production cross section
!   xspopex          ! population cross section summed over spin and parity
!   xspopnuc         ! population cross section per nucleus
!   xsreacinc        ! reaction cross section for incident channel
!   xstotinc         ! total cross section (neutrons only) for incident channel
! Constants
!   parN             ! neutron number of particle
!   parZ             ! charge number of particle
!
! *** Declaration of local data
!
  implicit none
  integer   :: i1            ! value
  integer   :: i2            ! value
  integer   :: idc           ! help variable
  integer   :: iiso          ! isomeric state
  integer   :: is            ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: is2           ! help variable
  integer   :: iyield        ! particle yield
  integer   :: mt            ! MT number
  integer   :: mt0           ! MT number
  integer   :: mtc           ! MT number for continuum
  integer   :: MTchan(nummt) ! channel index for MT number
  integer   :: mtd           ! MT number for discrete states
  integer   :: mtf           ! MT number for fission
  integer   :: nen           ! energy counter
  integer   :: nex           ! excitation energy bin of compound nucleus
  integer   :: Nix           ! neutron number index for residual nucleus
  integer   :: type          ! particle type
  integer   :: Zix           ! charge number index for residual nucleus
  real(sgl) :: Efac          ! help variable
  real(sgl) :: R             ! radius
  real(sgl) :: ratio         ! if 0<ratio<1 then x is between xl1 and xl2
  real(sgl) :: ratiogs       ! adjustment factor for ground state
  real(sgl) :: ratioiso      ! adjustment factor for isomer
  real(sgl) :: xsadd         ! total difference in cross section after normalization
  real(sgl) :: xsdif         ! difference between in-flux and out-flux per bin
  real(sgl) :: xsdifelas     ! difference in elastic cross section
  real(sgl) :: xsdifgs       ! difference in ground state cross section
  real(sgl) :: xsdifiso      ! difference in isomeric cross section
  real(sgl) :: xsdiftot      ! difference in total cross section
  real(sgl) :: xsfrac        ! fractional cross section
!
! ************************ Indices for MT numbers **********************
!
! Indices for MT numbers
!
  do mt = 1, nummt
    MTchan(mt) = - 1000
  enddo
  MTchan(4) = 100000
  MTchan(11) = 201000
  MTchan(16) = 200000
  MTchan(17) = 300000
  MTchan(22) = 100001
  MTchan(23) = 100003
  MTchan(24) = 200001
  MTchan(25) = 300001
  MTchan(28) = 110000
  MTchan(29) = 100002
  MTchan(30) = 200002
  MTchan(32) = 101000
  MTchan(33) = 100100
  MTchan(34) = 100010
  MTchan(35) = 101002
  MTchan(36) = 100102
  MTchan(37) = 400000
  MTchan(41) = 210000
  MTchan(42) = 310000
  MTchan(44) = 120000
  MTchan(45) = 110001
  MTchan(102) = 000000
  MTchan(103) = 010000
  MTchan(104) = 001000
  MTchan(105) = 000100
  MTchan(106) = 000010
  MTchan(107) = 000001
  MTchan(108) = 000002
  MTchan(109) = 000003
  MTchan(111) = 020000
  MTchan(112) = 010001
  MTchan(113) = 000102
  MTchan(114) = 001002
  MTchan(115) = 011000
  MTchan(116) = 010100
  MTchan(117) = 001001
  MTchan(152) = 500000
  MTchan(153) = 600000
  MTchan(154) = 200100
  MTchan(155) = 000101
  MTchan(156) = 410000
  MTchan(157) = 301000
  MTchan(158) = 101001
  MTchan(159) = 210001
  MTchan(160) = 700000
  MTchan(161) = 800000
  MTchan(162) = 510000
  MTchan(163) = 610000
  MTchan(164) = 710000
  MTchan(165) = 400001
  MTchan(166) = 500001
  MTchan(167) = 600001
  MTchan(168) = 700001
  MTchan(169) = 401000
  MTchan(170) = 501000
  MTchan(171) = 601000
  MTchan(172) = 300100
  MTchan(173) = 400100
  MTchan(174) = 500100
  MTchan(175) = 600100
  MTchan(176) = 200010
  MTchan(177) = 300010
  MTchan(178) = 400010
  MTchan(179) = 320000
  MTchan(180) = 300002
  MTchan(181) = 310001
  MTchan(182) = 001100
  MTchan(183) = 111000
  MTchan(184) = 110100
  MTchan(185) = 101100
  MTchan(186) = 110010
  MTchan(187) = 101010
  MTchan(188) = 100110
  MTchan(189) = 100101
  MTchan(190) = 220000
  MTchan(191) = 010010
  MTchan(192) = 001010
  MTchan(193) = 000011
  MTchan(194) = 420000
  MTchan(195) = 400002
  MTchan(196) = 410001
  MTchan(197) = 030000
  MTchan(198) = 130000
  MTchan(199) = 320001
  MTchan(200) = 520000
!
! ************************ Adjustment factors **************************
!
! Set incident energy dependent adjustment factors (purely for fitting purposes).
!
  do mt = 1, nummt
    do is = - 1, numisom
      Crescue(mt, is) = 1.
      if (Nrescue(mt, is) == 0) cycle
      if (nin == Ninclow + 1 .or. Einc <= Erescue(mt, is, 1)) then
        Crescue(mt, is) = frescue(mt, is, 1)
        goto 140
      endif
      if (Einc >= Erescue(mt, is, Nrescue(mt, is))) then
        Crescue(mt, is) = frescue(mt, is, Nrescue(mt, is))
        goto 140
      endif
      do nen = 1, Nrescue(mt, is) - 1
        if (Einc > Erescue(mt, is, nen) .and. Einc <= Erescue(mt, is, nen + 1)) then
          Efac = (Einc - Erescue(mt, is, nen)) / (Erescue(mt, is, nen + 1) - Erescue(mt, is, nen))
          Crescue(mt, is) = frescue(mt, is, nen) + Efac * (frescue(mt, is, nen + 1) - frescue(mt, is, nen))
          if (frescue(mt, is, nen + 1) > 1.e10) Crescue(mt, is) = frescue(mt, is, nen)
          if (frescue(mt, is, nen) > 1.e10) Crescue(mt, is) = frescue(mt, is, nen + 1)
          goto 140
        endif
      enddo
  140 if (Crescue(mt, is) == 0.) Crescue(mt, is) = 1.
    enddo
  enddo
!
! *************************** Normalization ****************************
!
! Application of normalization from rescue files.
! For each partial cross section the applied difference is stored, so that later it can be redistributed.
!
  if (Crescue(1, - 1) /= 1.) xsdiftot = xstotinc * (1. / Crescue(1, - 1) - 1.)
  if (Crescue(2, - 1) /= 1.) xsdifelas = xselastot * (1. / Crescue(2, - 1) - 1.)
  xsadd = 0.
  do mt = 4, nummt
    do is = - 1, numisom
      if (Crescue(mt, is) /= 1.) goto 217
    enddo
    cycle
  217   iiso = 0
    xsdifiso = 0.
    ratioiso = 1.
    do is = numisom, - 1, - 1
      if (is >= 0) then
        do is2 = - 1, numisom
          if (Crescue(mt, is2) /= 1.) goto 227
        enddo
        cycle
      endif
  227     xsdif = 0.
      ratio = 1. / Crescue(mt, is)
!
! Fission
!
      if (flagfission .and. mt == 18 .and. xsfistot > 0) then
        xsdif = xsfistot * (ratio - 1.)
        xsadd = xsadd + xsdif
        xsfistot = xsfistot * ratio
        xsfistot0 = xsfistot0 * ratio
        do mtf = 19, 38
          if (mtf > 21 .and. mtf < 38) cycle
          do idc = 0, idnum
            if ((mtf == 19 .and. idchannel(idc) == 000000) .or. (mtf == 20 .and. idchannel(idc) == 100000) .or. &
              (mtf == 21 .and. idchannel(idc) == 200000) .or. (mtf == 22 .and. idchannel(idc) == 300000)) then
              R = ratio
              if (Crescue(mtf, is) /= 1.) R = 1. / Crescue(mtf, is)
              xsfischannel(idc) = xsfischannel(idc) * R
            endif
          enddo
        enddo
        cycle
      endif
!
! Partial channels
!
      do idc = 0, idnum
        if (idchannel(idc) /= MTchan(mt)) cycle
        if (xschannel(idc) <= 0.) cycle
!
! Find associated residual nucleus
!
        Zix = 0
        Nix = 0
        do type = 1, 6
          iyield = mod(idchannel(idc), 10 **(7 - type)) / (10 **(6 - type))
          Zix = Zix + iyield * parZ(type)
          Nix = Nix + iyield * parN(type)
        enddo
        if (xspopnuc(Zix, Nix) > 0.) then
          xsfrac = xschannel(idc) / xspopnuc(Zix, Nix)
        else
          xsfrac = 1.
        endif
!
! Normalize ground state
!
        if (is == 0 .and. Crescue(mt, is) /= 1.) then
          ratiogs = 1. / Crescue(mt, is)
          xsdifgs = xschaniso(idc, 0) * (ratiogs - 1.)
          xschaniso(idc, 0) = xschaniso(idc, 0) * ratiogs
          R = 1. + (ratioiso - 1.) * xsfrac
          xspopex(Zix, Nix, 0) = xspopex(Zix, Nix, 0) * R
        endif
!
! Normalize isomer
!
        if (is == 1) then
          do i1 = 1, numlev
            if (xschaniso(idc, i1) /= 0.) then
              if (Crescue(mt, is) /= 1.) then
                ratioiso = 1. / Crescue(mt, is)
                xsdifiso = xschaniso(idc, i1) * (ratioiso - 1.)
                xschaniso(idc, i1) = xschaniso(idc, i1) * ratioiso
                R = 1. + (ratioiso - 1.) * xsfrac
                xspopex(Zix, Nix, i1) = xspopex(Zix, Nix, i1) * R
              endif
              iiso = i1
              exit
            endif
          enddo
        endif
!
! Correct ground state, isomer or total cross section
!
        if (is ==  - 1) then
          R = 1. + (ratioiso - 1.) * xsfrac
          if (Crescue(mt, 1) == 1.) then
            if (Crescue(mt, 0) /= 1.) then
              xschaniso(idc, iiso) = max(xschannel(idc) - xschaniso(idc, 0), 0.)
              xspopex(Zix, Nix, iiso) = max(xspopnuc(Zix, Nix) - xspopex(Zix, Nix, 0), 0.d0)
            else
              if (Crescue(mt, - 1) /= 1.) then
                xsdif = xschannel(idc) * (ratio - 1.)
                xschannel(idc) = xschannel(idc) * ratio
                xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) * R
                do i1 = 0, numlev
                  xschaniso(idc, i1) = xschaniso(idc, i1) * ratio
                  xspopex(Zix, Nix, i1) = xspopex(Zix, Nix, i1) * R
                enddo
                if (flagcompo) then
                  xspopdir(Zix, Nix) = xspopdir(Zix, Nix) * R
                  xspoppreeq(Zix, Nix) = xspoppreeq(Zix, Nix) * R
                  xspopcomp(Zix, Nix) = xspopcomp(Zix, Nix) * R
                endif
              endif
            endif
          else
            if (Crescue(mt, 0) /= 1.) then
              xsdif = xsdifiso + xsdifgs
              ratio = 1. + xsdif / xschannel(idc)
              xschannel(idc) = xschannel(idc) * ratio
              xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) * R
              if (flagcompo) then
                xspopdir(Zix, Nix) = xspopdir(Zix, Nix) * R
                xspoppreeq(Zix, Nix) = xspoppreeq(Zix, Nix) * R
                xspopcomp(Zix, Nix) = xspopcomp(Zix, Nix) * R
              endif
            else
              if (Crescue(mt, - 1) /= 1.) then
                xsdif = xschannel(idc) * (ratio - 1.)
                xschannel(idc) = xschannel(idc) * ratio
                xspopnuc(Zix, Nix) = xspopnuc(Zix, Nix) * R
                if (flagcompo) then
                  xspopdir(Zix, Nix) = xspopdir(Zix, Nix) * R
                  xspoppreeq(Zix, Nix) = xspoppreeq(Zix, Nix) * R
                  xspopcomp(Zix, Nix) = xspopcomp(Zix, Nix) * R
                endif
              endif
              xschaniso(idc, 0) = max(xschannel(idc) - xschaniso(idc, iiso), 0.)
              xspopex(Zix, Nix, 0) = xspopex(Zix, Nix, 0) * R
            endif
          endif
        endif
        if (is >= 0) cycle
        xsadd = xsadd + xsdif
!
! Normalization of related cross sections and/or spectra, by the rescue factors.
!
        xsgamchannel(idc) = xsgamchannel(idc) * ratio
        do i1 = 0, numlev
          do i2 = 0, numlev
            xsgamdischan(idc, i1, i2) = xsgamdischan(idc, i1, i2) * ratio
          enddo
        enddo
        if (flagspec) then
          do type = 0, 6
            do nen = 0, numen
              xschannelsp(idc, type, nen) = xschannelsp(idc, type, nen) * ratio
            enddo
          enddo
        endif
        do type = 1, 6
          iyield = mod(idchannel(idc), 10 **(7 - type)) / (10 **(6 - type))
          xsparticle(type) = xsparticle(type) + iyield * xsdif
          multiplicity(type) = xsparticle(type) / xsreacinc
        enddo
!
! Normalization of discrete level and continuum cross sections
!
! Due to possible memory limitation, we allow only individual discrete level adjustment for inelastic scattering.
!
        if (mt == 4 .or. (mt >= 103 .and. mt <= 107)) then
          if (mt == 4) then
            type = 1
            mt0 = 50
            mtc = 91
          else
            type = mt - 101
          endif
          do nex = 0, numlev
            R = ratio
            if (type == 1) then
              mtd = mt0 + nex
              if (Crescue(mtd, is) /= 1.) R = 1. / Crescue(mtd, is)
            endif
            xsdisc(type, nex) = xsdisc(type, nex) * R
            xscompdisc(type, nex) = xscompdisc(type, nex) * R
            xsdirdisc(type, nex) = xsdirdisc(type, nex) * R
          enddo
          R = ratio
          if (type == 1 .and. Crescue(mtc, is) /= 1.) R = 1. / Crescue(mtc, is)
          xsexclcont(type) = xsexclcont(type) * R
          xsngn(type) = xsngn(type) * ratio
          xsexclusive(type) = xsexclusive(type) * ratio
          xsdisctot(type) = xsdisctot(type) * ratio
          xsdirdisctot(type) = xsdirdisctot(type) * ratio
          xscompdisctot(type) = xscompdisctot(type) * ratio
          xscompcont(type) = xscompcont(type) * ratio
          xsdircont(type) = xsdircont(type) * ratio
          xsconttot(type) = xsconttot(type) * ratio
          xsdirect(type) = xsdirect(type) * ratio
          xsbinary(type) = xsbinary(type) * ratio
        endif
      enddo
    enddo
  enddo
!
! Put difference in the elastic (or total) cross section
!
  if (Crescue(2, - 1) /= 1.) then
    xseladjust(nin) = xsdifelas
    xselastot = xselastot + xseladjust(nin)
    xstotadjust(nin) = xsdifelas + xsadd
    xstotinc = xstotinc + xstotadjust(nin)
  else
    xstotadjust(nin) = xsdiftot
    xstotinc = xstotinc + xstotadjust(nin)
    xseladjust(nin) = xsdiftot - xsadd
    xselastot = xselastot + xseladjust(nin)
  endif
  xsnonadjust(nin) = xsadd
  xsnonel = xsnonel + xsnonadjust(nin)
  channelsum = channelsum + xsadd
  return
end subroutine normalization
! Copyright A.J. Koning 2021
