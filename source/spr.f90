subroutine spr
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : S, P and R' resonance parameters
!
! Author    : Arjan Koning
!
! 2025-12-13: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
! Variables for OMP
!   Rprime         ! potential scattering radius
!   Sstrength      ! s, p, d, etc - wave strength function
! Variables for basic reaction
!   flagendf       ! flag for information for ENDF - 6 file
!   flagendfdet    ! flag for detailed ENDF - 6 information per channel
! Variables for input energies
!   nin            ! counter for incident energy
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   Ztarget        ! charge number of target nucleus
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Variables for energies
!   Ninclow      ! number of incident energies below Elow
!   wavenum        ! wave number
! Variables for incident channel
!   Tlinc          ! transm. coeff. as a function of l for incident channel
!   xselasinc      ! total elastic cross section (neutrons only) for inc. channel
! Constants
!   fourpi         ! 4. * pi
!   onethird       ! 1 / 3
!   twopi          ! 2 * pi
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: topline    ! topline
  character(len=80) :: quantity   ! quantity
  real(sgl) :: Efac ! help variable
  real(sgl) :: r2k2 ! help variable
  real(sgl) :: Rpo  ! standard value for R
  integer   :: indent
  integer   :: id2
  integer   :: id4
!
! **************** S, P and R' resonance parameters ********************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  Rprime = 10.*sqrt(0.001*max(xselasinc, 0.)/fourpi)
  Efac = 1. / (sqrt(1.e6 * Einc) * twopi)
  Rpo = 1.35 * Atarget **onethird
  r2k2 = Rpo * Rpo * wavenum * wavenum
  Sstrength(0) = Tlinc(0) * Efac
  Sstrength(1) = Tlinc(1) * Efac * (1. + r2k2) / r2k2
  Sstrength(2) = Tlinc(2) * Efac * (9. + 3. * r2k2 + r2k2 * r2k2) / (r2k2 * r2k2)
  if ((flagoutomp .or. (flagendf .and. flagendfdet)) .and. (Einc <= 0.1 .or. nin == Ninclow + 1)) then
    open (unit = 1, file = 'spr.opt', status = 'replace')
    quantity='basic OMP quantities'
    topline=trim(targetnuclide)//' '//trim(quantity)
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_char(id2,'parameters','')
    call write_real(id4,'S0',Sstrength(0)*1.e4)
    if (S0(0,0) > 0.) then
      call write_real(id4,'experimental S0',S0(0,0))
      call write_real(id4,'experimental S0 unc.',dS0(0,0))
      call write_real(id4,'C/E S0',Sstrength(0)*1.e4/S0(0, 0))
    endif
    call write_real(id4,'S1',Sstrength(1)*1.e4)
    if (S1(0,0) > 0.) then
      call write_real(id4,'experimental S1',S1(0,0))
      call write_real(id4,'experimental S1 unc.',dS1(0,0))
      call write_real(id4,'C/E S1',Sstrength(1)*1.e4/S1(0, 0))
    endif
    call write_real(id4,'Rprime [fm]',Rprime)
    if (Rscat(0,0) > 0.) then
      call write_real(id4,'experimental Rprime [fm]',Rscat(0,0))
      call write_real(id4,'experimental Rprime unc. [fm]',dRscat(0,0))
      call write_real(id4,'C/E Rprime',Rprime/Rscat(0, 0))
    endif
!   write(1, '(2i4, 3f8.4)') Atarget, Ztarget, Sstrength(0)*1.e4, Sstrength(1) * 1.e4, Rprime
    close (unit = 1)
  endif
  return
end subroutine spr
! Copyright A.J. Koning 2021
