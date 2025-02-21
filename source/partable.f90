subroutine partable(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write model parameters per nucleus to separate file
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
! Variables for level density
!   D0               ! s - wave resonance spacing in eV
! Variables for preequilibrium
!   g                ! single - particle level density parameter
!   gadjust          ! adjustable factor for single - particle particle-hole states
!   gn               ! single - particle neutron level density parameter
!   gnadjust         ! adjustable factor for single - particle proton par
!   gp               ! single - particle proton level density parameter
!   gpadjust         ! adjustable factor for single - particle neutron pa
!   phmodel          ! particle - hole state density model
! Variables for fission
!   bdamp            ! fission partial damping parameter
!   bdampadjust      ! correction for fission partial damping parameter
!   betafiscor       ! adjustable factor for fission path width
!   betafiscoradjust ! adjustable factor for fission path width
!   fbaradjust       ! adjustable factor for fission parameters
!   fbarrier         ! height of fission barrier
!   fismodelx        ! fission model
!   flagfission      ! flag for fission
!   fwidth           ! width of fission barrier
!   fwidthadjust     ! adjustable factor for fission parameters
!   vfiscor          ! adjustable factor for fission path height
!   vfiscoradjust    ! adjustable factor for fission path height
!   widthc2          ! width of class2 states
! Variables for gamma rays
!   egr              ! energy of GR
!   egradjust        ! adjustable factor for energy of GR
!   etable           ! constant to adjust tabulated strength functions
!   etableadjust     ! correction to adjust tabulated strength functions
!   ftable           ! constant to adjust tabulated strength functions
!   ftableadjust     ! correction to adjust tabulated strength functions
!   gamgam           ! total radiative width in eV
!   gamgamadjust     ! adjustable factor for radiative parameter
!   gammax           ! number of l - values for gamma multipolarity
!   ggr              ! width of GR
!   ggradjust        ! adjustable factor for width of GR
!   sgr              ! strength of GR
!   sgradjust        ! adjustable factor for strength of GR
!   strength         ! E1 strength function model
!   strengthM1       ! model for M1 gamma - ray strength function
!   upbend           ! upbend of M1
!   wtable           ! constant to adjust tabulated strength functions
!   wtableadjust     ! correction to adjust tabulated strength functions
! Variables for level density
!   aadjust          ! adjustable factor for level density parameter
!   alev             ! level density parameter
!   ctable           ! constant to adjust tabulated level densities
!   ctableadjust     ! correction to adjust tabulated level densities
!   deltaW           ! shell correction in nuclear mass
!   E0               ! particle constant of temperature formula
!   E0adjust         ! adjustable factor for E0
!   Exmatch          ! matching point for Ex
!   Exmatchadjust    ! adjustable factor for matching energy
!   flagcolall       ! flag for collective enhancement of level density
!   gammald          ! gamma - constant for asymptotic level density para
!   Krotconstant     ! normalization constant for rotational enhancemen
!   ldmodel          ! level density model
!   Nlow             ! lowest discrete level for temperature matching
!   Ntop             ! highest discrete level for temperature matching
!   pair             ! pairing energy
!   Pshift           ! adjustable pairing shift
!   Pshiftadjust     ! adjustable correction to pairing shift
!   ptable           ! constant to adjust tabulated level densities
!   ptableadjust     ! correction to adjust tabulated level densities
!   Rclass2mom       ! norm. constant for moment of inertia for class 2
!   Rtransmom        ! norm. constant for moment of inertia for transit
!   s2adjust         ! adjustable constant (Z, A, barrier - dependent) for
!   T                ! temperature
!   Tadjust          ! adjustable factor for temperature
!   Ufermi           ! energy of Fermi distribution for damping of rotational effects
!   cfermi           ! width of Fermi distribution for damping of rotational effects
!
! Variables for nuclides
!   AA               ! mass number of residual nucleus
!   ZZ               ! charge number of residual nucleus
! Constants
!   nuc              ! symbol of nucleus
!  Variables for gamma-ray strength functions
!   ngr              ! number of GR
! Variables for fission parameters
!   nfisbar          ! number of fission barrier parameters
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring
  character(len=6)  :: finalnuclide
  integer :: A                 ! mass number of target nucleus
  integer :: ibar              ! fission barrier
  integer :: l                 ! multipolarity
  integer :: Nix               ! neutron number index for residual nucleus
  integer :: Z                 ! charge number of target nucleus
  integer :: Zix               ! charge number index for residual nucleus
!
! ****************************** Z and A of nucleus ********************
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  massstring = '   '
  write(massstring,'(i3)') A
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  write(51, '("## residual:")')
  write(51, '("##   Z: ", i0)') Z
  write(51, '("##   A: ", i0)') A
  write(51, '("##   nuclide: ", a)') finalnuclide
! ********************** Level density parameters **********************
!
  write(51, '("## parameters: ")')
  write(51, '("##   level density")')
  write(51, '("a              ", 2i4, f10.5)') Z, A, alev(Zix, Nix)
  write(51, '("aadjust        ", 2i4, f10.5)') Z, A, aadjust(Zix, Nix)
  write(51, '("gammald        ", 2i4, f10.5)') Z, A, gammald(Zix, Nix)
  write(51, '("pair           ", 2i4, f10.5)') Z, A, pair(Zix, Nix)
  do ibar = 0, nfisbar(Zix, Nix)
    write(51, '("Pshift         ", 2i4, f10.5, i4)') Z, A, Pshift(Zix, Nix, ibar), ibar
    write(51, '("Pshiftadjust   ", 2i4, f10.5, i4)') Z, A, pshiftadjust(Zix, Nix, ibar), ibar
    write(51, '("deltaW         ", 2i4, f10.5, i4)') Z, A, deltaW(Zix, Nix, ibar), ibar
    if (ldmodel(Zix, Nix) == 1) then
      write(51, '("T              ", 2i4, f10.5, i4)') Z, A, T(Zix, Nix, ibar), ibar
      write(51, '("E0             ", 2i4, f10.5, i4)') Z, A, E0(Zix, Nix, ibar), ibar
      write(51, '("Exmatch        ", 2i4, f10.5, i4)') Z, A, Exmatch(Zix, Nix, ibar), ibar
      write(51, '("Tadjust        ", 2i4, f10.5, i4)') Z, A, Tadjust(Zix, Nix, ibar), ibar
      write(51, '("E0adjust       ", 2i4, f10.5, i4)') Z, A, E0adjust(Zix, Nix, ibar), ibar
      write(51, '("Exmatchadjust  ", 2i4, f10.5, i4)') Z, A, Exmatchadjust(Zix, Nix, ibar), ibar
    endif
    write(51, '("Ntop           ", 2i4, 2i4)') Z, A, Ntop(Zix, Nix, ibar), ibar
    write(51, '("Nlow           ", 2i4, 2i4)') Z, A, Nlow(Zix, Nix, ibar), ibar
    write(51, '("s2adjust       ", 2i4, f10.5, i4)') Z, A, s2adjust(Zix, Nix, ibar), ibar
    write(51, '("ctable         ", 2i4, f10.5, i4)') Z, A, ctable(Zix, Nix, ibar), ibar
    write(51, '("ptable         ", 2i4, f10.5, i4)') Z, A, ptable(Zix, Nix, ibar), ibar
    write(51, '("ctableadjust   ", 2i4, f10.5, i4)') Z, A, ctableadjust(Zix, Nix, ibar), ibar
    write(51, '("ptableadjust   ", 2i4, f10.5, i4)') Z, A, ptableadjust(Zix, Nix, ibar), ibar
    write(51, '("Ufermi         ", 2i4, f10.5, i4)') Z, A, Ufermi(Zix, Nix, ibar), ibar
    write(51, '("cfermi         ", 2i4, f10.5, i4)') Z, A, cfermi(Zix, Nix, ibar), ibar
    if (flagcolall) write(51, '("Krotconstant   ", 2i4, f10.5, i4)') Z, A, Krotconstant(Zix, Nix, ibar), ibar
  enddo
  if (D0(Zix, Nix) /= 0.) write(51, '("D0             ", 2i4, es12.5)') Z, A, D0(Zix, Nix) * 0.001
  if (phmodel == 1) then
    write(51, '("g              ", 2i4, f10.5)') Z, A, g(Zix, Nix)
    write(51, '("gp             ", 2i4, f10.5)') Z, A, gp(Zix, Nix)
    write(51, '("gn             ", 2i4, f10.5)') Z, A, gn(Zix, Nix)
    write(51, '("gnadjust       ", 2i4, f10.5)') Z, A, gnadjust(Zix, Nix)
    write(51, '("gpadjust       ", 2i4, f10.5)') Z, A, gpadjust(Zix, Nix)
    write(51, '("gadjust        ", 2i4, f10.5)') Z, A, gadjust(Zix, Nix)
  endif
  write(51, '("Risomer        ", 2i4, f10.5)') Z, A, Risomer(Zix, Nix)
!
! ************************ Gamma-ray parameters ************************
!
  write(51, '("## parameters: ")')
  write(51, '("##   photon strength function")')
  write(51, '("gamgam         ", 2i4, f10.5)') Z, A, gamgam(Zix, Nix)
  write(51, '("gamgamadjust   ", 2i4, f10.5)') Z, A, gamgamadjust(Zix, Nix)
  do l = 1, gammax
    if (strength <= 2 .or. strength == 5) then
      write(51, '("sgr            ", 2i4, f8.3, " E", i1)') Z, A, sgr(Zix, Nix, 1, l, 1), l
      write(51, '("egr            ", 2i4, f8.3, " E", i1)') Z, A, egr(Zix, Nix, 1, l, 1), l
      write(51, '("ggr            ", 2i4, f8.3, " E", i1)') Z, A, ggr(Zix, Nix, 1, l, 1), l
      write(51, '("sgradjust      ", 2i4, f8.3, " E", i1)') Z, A, sgradjust(Zix, Nix, 1, l, 1), l
      write(51, '("egradjust      ", 2i4, f8.3, " E", i1)') Z, A, egradjust(Zix, Nix, 1, l, 1), l
      write(51, '("ggradjust      ", 2i4, f8.3, " E", i1)') Z, A, ggradjust(Zix, Nix, 1, l, 1), l
    else
      write(51, '("etable         ", 2i4, f10.5, " E", i1)') Z, A, etable(Zix, Nix, 1, l), l
      write(51, '("ftable         ", 2i4, f10.5, " E", i1)') Z, A, ftable(Zix, Nix, 1, l), l
      write(51, '("wtable         ", 2i4, f10.5, " E", i1)') Z, A, wtable(Zix, Nix, 1, l), l
      write(51, '("etableadjust   ", 2i4, f10.5, " E", i1)') Z, A, etableadjust(Zix, Nix, 1, l), l
      write(51, '("ftableadjust   ", 2i4, f10.5, " E", i1)') Z, A, ftableadjust(Zix, Nix, 1, l), l
      write(51, '("wtableadjust   ", 2i4, f10.5, " E", i1)') Z, A, wtableadjust(Zix, Nix, 1, l), l
    endif
    if (strengthM1 == 8 .or. strengthM1 == 10) then
      write(51, '("etable         ", 2i4, f10.5, " M", i1)') Z, A, etable(Zix, Nix, 0, l), l
      write(51, '("ftable         ", 2i4, f10.5, " M", i1)') Z, A, ftable(Zix, Nix, 0, l), l
      write(51, '("wtable         ", 2i4, f10.5, " M", i1)') Z, A, wtable(Zix, Nix, 0, l), l
      write(51, '("etableadjust   ", 2i4, f10.5, " M", i1)') Z, A, etableadjust(Zix, Nix, 0, l), l
      write(51, '("ftableadjust   ", 2i4, f10.5, " M", i1)') Z, A, ftableadjust(Zix, Nix, 0, l), l
      write(51, '("wtableadjust   ", 2i4, f10.5, " M", i1)') Z, A, wtableadjust(Zix, Nix, 0, l), l
    endif
    if (strengthM1 >= 3) then
      write(51, '("upbendc        ", 2i4, es12.5, " M", i1)') Z, A, upbend(Zix, Nix, 0, l, 1), l
      write(51, '("upbende        ", 2i4, es12.5, " M", i1)') Z, A, upbend(Zix, Nix, 0, l, 2), l
      write(51, '("upbendf        ", 2i4, es12.5, " M", i1)') Z, A, upbend(Zix, Nix, 0, l, 3), l
      write(51, '("upbendcadjust  ", 2i4, es12.5, " M", i1)') Z, A, upbendadjust(Zix, Nix, 0, l, 1), l
      write(51, '("upbendeadjust  ", 2i4, es12.5, " M", i1)') Z, A, upbendadjust(Zix, Nix, 0, l, 2), l
      write(51, '("upbendfadjust  ", 2i4, es12.5, " M", i1)') Z, A, upbendadjust(Zix, Nix, 0, l, 3), l
    endif
    if (ngr(Zix, Nix, 1, l) == 2) then
      write(51, '("sgr            ", 2i4, f8.3, " E", i1, " 2")') Z, A, sgr(Zix, Nix, 1, l, 2), l
      write(51, '("egr            ", 2i4, f8.3, " E", i1, " 2")') Z, A, egr(Zix, Nix, 1, l, 2), l
      write(51, '("ggr            ", 2i4, f8.3, " E", i1, " 2")') Z, A, ggr(Zix, Nix, 1, l, 2), l
      write(51, '("sgradjust      ", 2i4, f8.3, " E", i1, " 2")') Z, A, sgradjust(Zix, Nix, 1, l, 2), l
      write(51, '("egradjust      ", 2i4, f8.3, " E", i1, " 2")') Z, A, egradjust(Zix, Nix, 1, l, 2), l
      write(51, '("ggradjust      ", 2i4, f8.3, " E", i1, " 2")') Z, A, ggradjust(Zix, Nix, 1, l, 2), l
    endif
    write(51, '("sgr            ", 2i4, f8.3, " M", i1)') Z, A, sgr(Zix, Nix, 0, l, 1), l
    write(51, '("egr            ", 2i4, f8.3, " M", i1)') Z, A, egr(Zix, Nix, 0, l, 1), l
    write(51, '("ggr            ", 2i4, f8.3, " M", i1)') Z, A, ggr(Zix, Nix, 0, l, 1), l
    write(51, '("sgradjust      ", 2i4, f8.3, " M", i1)') Z, A, sgradjust(Zix, Nix, 0, l, 1), l
    write(51, '("egradjust      ", 2i4, f8.3, " M", i1)') Z, A, egradjust(Zix, Nix, 0, l, 1), l
    write(51, '("ggradjust      ", 2i4, f8.3, " M", i1)') Z, A, ggradjust(Zix, Nix, 0, l, 1), l
  enddo
!
! ************************** Fission parameters ************************
!
  if (flagfission) then
    write(51, '("##")')
    write(51, '("## parameters: ")')
    write(51, '("##   fission")')
    write(51, '("##")')
    do ibar = 1, nfisbar(Zix, Nix)
      if (fismodelx(Zix, Nix) == 5) then
        if (ibar == 1) then
          write(51, '("betafiscor     ", 2i4, f10.5)') Z, A, betafiscor(Zix, Nix)
          write(51, '("vfiscor        ", 2i4, f10.5)') Z, A, vfiscor(Zix, Nix)
          write(51, '("betafiscoradjust ", 2i4, f10.5)') Z, A, betafiscoradjust(Zix, Nix)
          write(51, '("vfiscoradjust  ", 2i4, f10.5)') Z, A, vfiscoradjust(Zix, Nix)
        endif
        write(51, '("bdamp          ", 2i4, f10.5, i3)') Z, A, bdamp(Zix, Nix, ibar), ibar
        write(51, '("bdampadjust    ", 2i4, f10.5, i3)') Z, A, bdampadjust(Zix, Nix, ibar), ibar
      else
        write(51, '("fisbar         ", 2i4, f10.5, i3)') Z, A, fbarrier(Zix, Nix, ibar), ibar
        write(51, '("fishw          ", 2i4, f10.5, i3)') Z, A, fwidth(Zix, Nix, ibar), ibar
        write(51, '("fisbaradjust   ", 2i4, f10.5, i3)') Z, A, fbaradjust(Zix, Nix, ibar), ibar
        write(51, '("fishwadjust    ", 2i4, f10.5, i3)') Z, A, fwidthadjust(Zix, Nix, ibar), ibar
      endif
      if (ibar < nfisbar(Zix, Nix)) write(51, '("class2width    ", 2i4, f10.5, i3)') Z, A, widthc2(Zix, Nix, ibar), ibar
      write(51, '("Rtransmom      ", 2i4, f10.5, i3)') Z, A, Rtransmom(Zix, Nix, ibar), ibar
      write(51, '("Rclass2mom     ", 2i4, f10.5, i3)') Z, A, Rclass2mom(Zix, Nix, ibar), ibar
    enddo
  endif
  return
end subroutine partable
! Copyright A.J. Koning 2021
