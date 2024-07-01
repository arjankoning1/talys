subroutine talysinput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : User input, defaults and setting variables
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2024-06-25: Added input_fit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ******************** Set defaults and read input *********************
!
  call readinput
  call input_path
  call input_main
  call input_best
  call input_basicreac
  call read_energies
  call initial_best
  call input_basicpar
  call initial_adjust
  call input_astro
  call input_numerics
  call input_mass
  call input_levels
  call input_densitymodel
  call input_densitypar
  call input_gammamodel
  call input_gammapar
  call input_ompmodel
  call input_omppar
  call input_compoundmodel
  call input_compoundpar
  call input_directmodel
  call input_directpar
  call input_fissionmodel
  call input_fissionpar
  call input_preeqmodel
  call input_preeqpar
  call input_medical
  call input_fit
  call input_output
  call input_outfiles
  call checkkeyword
  call checkvalue
  return
end subroutine talysinput
! Copyright A.J. Koning 2021
