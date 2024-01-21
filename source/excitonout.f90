subroutine excitonout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of exciton model parameters
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
!   sgl           ! single precision kind
! Variables for preequilibrium
!   M2constant    ! constant for matrix element in exciton model
! Variables for main input
!   Ainit         ! mass number of initial compound nucleus
! Constants
!   hbar          ! Planck's constant / 2.pi in MeV.s
!   parname       ! name of particle
! Variables for exciton model initialization
!   Qfactor       ! Q - factor for neutron / proton distinction
! Variables for preequilibrium
!   Ecomp         ! total energy of composite system
!   p0            ! initial particle number
! Variables for preequilibrium initialization
!   maxpar        ! maximal particle number
! Variables for exciton model
!   depletion     ! depletion factor at each stage
!   M2            ! square of matrix element
!   tauexc        ! lifetime of exciton state
!   wemispart     ! emission rate per particle and exciton number
!   wemistot      ! total emission rate per exciton number
!
! *** Declaration of local data
!
  implicit none
  integer   :: h          ! help variable
  integer   :: n          ! exciton number
  integer   :: p          ! particle number
  integer   :: type       ! particle type
  real(sgl) :: lambdaplus ! transition rate for n --> n+2
!
! ************************ Exciton model *******************************
!
  write(*, '(/" ++++++++++ EXCITON MODEL ++++++++++")')
!
! 1. Output of matrix element
!
! matrix    : subroutine for matrix element for exciton model
!
  write(*, '(/" 1. Matrix element for E= ", f8.3/)') Ecomp
  write(*, '(" Constant for matrix element: ", f7.3/)') M2constant
  write(*, '("  p h       M2"/)')
  do p = p0, maxpar
    h = p - p0
    n = p + h
    call matrix(Ainit, n)
    write(*, '(1x, 2i2, 2x, es12.5)') p, h, M2
  enddo
!
! 2. Output of Q-factors
!
  write(*, '(/" 2. Q-factors"/)')
  write(*, '("  p h  ", 6(a8, 1x), /)') (parname(type), type = 1, 6)
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, 7f9.5)') p, h, (Qfactor(type, p), type = 1, 6)
  enddo
!
! 3. Output of emission rates or escape widths
!
  write(*, '(/" 3. Emission rates or escape widths"/)')
  write(*, '(" A. Emission rates ( /sec)"/)')
  write(*, '("  p h", 3x, 7(a8, 4x), "Total"/)') (parname(type), type = 0, 6)
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, 8es12.5)') p, h, (wemispart(type, p, h), type = 0, 6), wemistot(p, h)
  enddo
  write(*, '(/" B. Escape widths (MeV)"/)')
  write(*, '("  p h", 2x, 7(a8, 4x), "Total"/)') (parname(type), type = 0, 6)
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, 8es12.5)') p, h, (wemispart(type, p, h) * hbar, type = 0, 6), wemistot(p, h) * hbar
  enddo
!
! 4. Output of transition rates or damping widths and total widths
!
  write(*, '(/" 4. Internal transition rates or damping widths, total widths"/)')
  write(*, '(" A. Internal transition rates ( /sec)"/)')
  write(*, '("  p h    lambdaplus"/)')
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, es15.6)') p, h, lambdaplus(0, 0, p, h)
  enddo
  write(*, '(/" B. Damping widths (MeV)"/)')
  write(*, '("  p h     gammaplus"/)')
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, es15.6)') p, h, lambdaplus(0, 0, p, h)*hbar
  enddo
  write(*, '(/" C. Total widths (MeV)"/)')
  write(*, '("  p h     gammatot"/)')
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, es15.6)') p, h, (lambdaplus(0, 0, p, h) + wemistot(p, h)) * hbar
  enddo
!
! 5. Output of depletion factors
!
  write(*, '(/" 5. Depletion factors"/)')
  write(*, '("  p h  depletion"/)')
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, f10.5)') p, h, depletion(p, h)
  enddo
!
! 6. Output of lifetimes of exciton states
!
  write(*, '(/" 6. Lifetimes")')
  write(*, '(/"  p h   mean lifetime"/)')
  do p = p0, maxpar
    h = p - p0
    write(*, '(1x, 2i2, es15.6)') p, h, tauexc(p, h)
  enddo
  return
end subroutine excitonout
! Copyright A.J. Koning 2021
