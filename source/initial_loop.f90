subroutine initial_loop
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Set variables for possible loop over nuclides
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
! Variables for input energies
!   energyfile     ! file with energies for OMP calculation
! Variables for main input
!   Arp            ! A of residual product
!   Atarget        ! mass number of target nucleus
!   Atarget0       ! mass number of target nucleus
!   Ltarget        ! excited level of target
!   ptype0         ! type of incident particle
!   Starget        ! symbol of target nucleus
!   Zrp            ! Z of residual product
!   Ztarget        ! charge number of target nucleus
!   Ztarget0       ! charge number of target nucleus
! Variables for basic reaction
!   flagffruns     ! flag to designate subsequent evaporation of fission products
!   flagrpruns     ! flag to designate that run is for residual product
! Constants
!   nuc            ! symbol of nucleus
! Variables for mass distribution
!   Aff            ! mass number of fission fragment
!   Zff            ! charge number of fission fragment
!
! *** Declaration of local data
!
  implicit none
!
! ************ Read first set of variables from input lines ************
!
! 1. Initializations
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified and the corresponding values are read.
! Erroneous input is immediately checked.
! The keywords and number of values on each line are retrieved from the input.
!
! ************* Process first set of input variables *******************
!
! Special case for loop over fission fragments or residual products to be evaporated.
!
  if (flagffruns) then
    Ztarget = Zff
    Atarget = Aff
    energyfile = 'ff000000.ex'
  endif
  if (flagrpruns) then
    Ztarget = Zrp
    Atarget = Arp
    energyfile = 'rp000000.ex'
  endif
  if (flagffruns .or. flagrpruns) then
    Starget = nuc(Ztarget)
    Ltarget = 0
    ptype0 = '0'
    write(energyfile(3:5), '(i3.3)') Ztarget
    write(energyfile(6:8), '(i3.3)') Atarget
  else
    Ztarget0 = Ztarget
    Atarget0 = Atarget
  endif
  return
end subroutine initial_loop
! Copyright A.J. Koning 2021
