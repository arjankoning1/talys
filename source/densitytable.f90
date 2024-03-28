subroutine densitytable(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Tabulated level densities
!
! Author    : Arjan Koning and Marieke Duijvestijn
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
!   dbl            ! double precision kind
! All global variables
!   numJ           ! maximum J - value
! Variables for fission
!   axtype         ! type of axiality of barrier
!   flagfission    ! flag for fission
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Variables for level density
!   flagparity     ! flag for non - equal parity distribution
!   ldmodel        ! level density model
! Variables for files
!   path           ! directory containing files to be read
! Constants
!   nuc            ! symbol of nucleus
!   twopi          ! 2 * pi
! Variables for level density
!   edens          ! energy grid for tabulated level densities
!   ldexist        ! flag for existence of level density table
!   ldtable        ! level density from table
!   ldtottable     ! total level density per parity from table
!   ldtottableP    ! total level density per parity from table
!   nendens        ! number of energies for level density grid
! Variables for fission parameters
!   nfisbar        ! number of fission barrier parameters
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist           ! logical to determine existence
  character(len=2)  :: denchar          ! string for level density file
  character(len=132):: denfile          ! level density parameter file
  integer           :: A                ! mass number of target nucleus
  integer           :: ia               ! mass number from abundance table
  integer           :: ibar             ! fission barrier
  integer           :: istat            ! error code
  integer           :: J                ! spin of level
  integer           :: JJ               ! help variable
  integer           :: Jmid             ! help variable
  integer           :: ldmod            ! level density model
  integer           :: nex              ! excitation energy bin of compound nucleus
  integer           :: Nix              ! neutron number index for residual nucleus
  integer           :: nloop            ! help variable
  integer           :: parity           ! parity
  integer           :: ploop            ! help variable
  integer           :: Z                ! charge number of target nucleus
  integer           :: Zix              ! charge number index for residual nucleus
  real(sgl)         :: ald              ! level density parameter
  real(sgl)         :: Eex              ! excitation energy
  real(sgl)         :: ignatyuk         ! function for energy dependent level density parameter a
  real(sgl)         :: Ktriax           ! level density enhancement factor for triaxial shapes
  real(sgl)         :: spincut          ! spin cutoff factor
  real(sgl)         :: term             ! help variable
  real(sgl)         :: dJ               ! help variable
  real(sgl)         :: Jold             ! help variable
  real(sgl)         :: factor           ! help variable
  real(dbl)         :: sumJ             ! help variable
  real(dbl)         :: Rsum             ! help variable
  real(dbl)         :: Rd               ! help variable
  real(dbl)         :: jt               ! help variable
  real(dbl)         :: ldsum            ! help variable
  real(dbl)         :: Rdis(0:numJ)     ! spin dependent level density
  real(dbl)         :: Rdisnew(0:numJ)  ! spin dependent level density
  real(dbl)         :: ld2j1(0:numJ)    ! spin dependent level density
  real(dbl)         :: ldtot            ! total level density
  real(dbl)         :: pardisloc        ! variable to account for parity distribution
!
! *********** Tabulated level densities from Goriely *******************
!
  Z = ZZ(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  denchar = trim(nuc(Z))
  nloop = 0
  if (flagfission) nloop = nfisbar(Zix, Nix)
  ldmod = ldmodel(Zix, Nix)
  if (ldmod >= 5) then
    ploop = - 1
    pardisloc = 1.
  else
    ploop = 1
    pardisloc = 0.5
  endif
  do ibar = 0, nloop
    ldexist(Zix, Nix, ibar) = .false.
    denfile = '                                                      '
!
! Ground state
!
    if (ibar == 0) then
      if (ldmod == 4) then
        denfile = trim(path)//'density/ground/goriely/'//trim(denchar)//'.tab'
      endif
      if (ldmod == 5) then
        denfile = trim(path)//'density/ground/hilaire/'//trim(denchar)//'.tab'
      endif
      if (ldmod == 6) then
        denfile = trim(path)//'density/ground/hilaireD1M/'// trim(denchar)//'.tab'
      endif
      if (ldmod == 7) then
        denfile = trim(path)//'density/ground/bskg3/'// trim(denchar)//'.tab'
      endif
      if (densfile(Zix, Nix)(1:1) /= ' ') denfile = densfile(Zix,Nix)
    endif
!
! First barrier
!
    if (ibar == 1) then
      if (ldmod == 4) then
        denfile = trim(path)//'density/fission/goriely/inner/' //trim(denchar)//'.ld'
      else
        denfile = trim(path)//'density/fission/hilaire/Max1/' //trim(denchar)//'.ld'
      endif
    endif
!
! Second barrier
!
    if (ibar == 2) then
      if (ldmod == 4) then
        denfile = trim(path)//'density/fission/goriely/outer/' //trim(denchar)//'.ld'
      else
        denfile = trim(path)//'density/fission/hilaire/Max2/' //trim(denchar)//'.ld'
      endif
    endif
!
! Third barrier
!
    if (ibar == 3) then
      if (ldmod == 5) then
        denfile = trim(path)//'density/fission/hilaire/Max3/' //trim(denchar)//'.ld'
      else
        cycle
      endif
    endif
!
! Check existence of file and read data from the tables.
!
    inquire (file = denfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = denfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(denfile, istat)
Loop1: do parity = 1, ploop, -2
        do
          read(2, '(/31x, i3//)', iostat = istat) ia
          if (istat == -1) exit Loop1
          if (istat /= 0) call read_error(denfile, istat)
          if (A == ia) then
            ldexist(Zix, Nix, ibar) = .true.
            do nex = 1, nendens(Zix, Nix)
              ld2j1=0.
              read(2, '(24x, e9.2, 9x, 30e9.2)', iostat = istat) ldtot, (ld2j1(J), J = 0, 29)
              if (istat /= 0) call read_error(denfile, istat, ival = A, xval = real(ld2j1(0)))
!
! Extend tabulated level densities up to numJ using linear extrapolation
!
              ldsum = ldtot
              do J = 30, numJ
                factor = 1. - (J - 29.) / (numJ - 29.)
                ld2j1(J) = ld2j1(29) * factor
                ldsum = ldsum + ld2j1(J)
              enddo
              if (ldsum > 0.) then
                do J = 0, numJ
                  ld2j1(J) = ld2j1(J) * ldtot / ldsum
                enddo
              endif
!
! Adjust spin distribution if requested
!
              if ((Rspincut /= 1. .or. s2adjust(Zix,Nix,ibar) /= 1.) .and. ldtot > 0.) then
                jt=Rspincut * s2adjust(Zix,Nix,ibar)
                Jmid=0
                Rsum=0.
                do J=0,numJ
                  Rdis(J)=ld2j1(J)/ldtot
                  Rsum=Rsum+Rdis(J)
                enddo
                if (Rsum > 0.) then
                  do J=0,numJ
                    Rdis(J)=Rdis(J)/Rsum
                  enddo
                  sumJ=0.
                  do J=0,numJ
                    sumJ=sumJ+Rdis(J)
                    if (sumJ >= 0.5) then
                      Jmid=max(J-1,0)
                      exit
                    endif
                  enddo
                  sumJ=0.
                  Rdisnew=0.
                  do J = 0, numJ - 1
                    dJ=real(J-Jmid)
                    Jold=Jmid+dJ/jt
                    JJ=int(Jold)
                    if (JJ < 0) cycle
                    if (JJ > 28) cycle
                    Rd=Rdis(JJ)+(Jold-JJ)*(Rdis(JJ+1)-Rdis(JJ))
                    Rdisnew(J)=Rd
                    sumJ=sumJ+Rdisnew(J)
                  enddo
                  if (sumJ > 0.) then
                    do J=0,numJ
                      Rdisnew(J)=Rdisnew(J)/sumJ
                    enddo
                    do J=0,numJ
                      ld2j1(J)=ldtot*Rdisnew(J)
                    enddo
                  endif
                endif
              endif
!
! Determination of the mass-asymmetric enhancement factor for fission barrier
!
! ignatyuk: function for energy dependent level density parameter a
!
              Ktriax = 1.
              if (ibar > 0) then
                if (axtype(Zix, Nix, ibar) == 2) Ktriax = 2.
                if (axtype(Zix, Nix, ibar) >= 3) then
                  Eex = edens(nex)
                  ald = ignatyuk(Zix, Nix, Eex, ibar)
                  term = sqrt(spincut(Zix, Nix, ald, Eex, ibar, 0))
                  if (axtype(Zix, Nix, ibar) == 3) Ktriax = 0.5 * sqrt(twopi) * term
                  if (axtype(Zix, Nix, ibar) == 4) Ktriax = sqrt(twopi) * term
                  if (axtype(Zix, Nix, ibar) == 5) Ktriax = 2. * sqrt(twopi) * term
                endif
              endif
              ldtottableP(Zix, Nix, nex, parity, ibar) = pardisloc * ldtot * Ktriax
              if (ploop == 1) ldtottableP(Zix, Nix, nex, - 1, ibar) = pardisloc * ldtot * Ktriax
              ldtottable(Zix, Nix, nex, ibar) = ldtottable(Zix, Nix, nex, ibar) + ldtot
              do J = 0, numJ
                ldtable(Zix, Nix, nex, J, parity, ibar) = pardisloc * ld2j1(J) * Ktriax
              enddo
              if (ploop == 1) then
                do J = 0, numJ
                  ldtable(Zix, Nix, nex, J, - 1, ibar) = pardisloc * ld2j1(J) * Ktriax
                enddo
              endif
            enddo
            read(2, '()')
            exit
          else
            do nex = 1, nendens(Zix, Nix) + 1
              read(2, '()')
            enddo
          endif
        enddo
      enddo Loop1
      close (unit = 2)
    endif
!
! Special case: make parity-independent level densities from parity-dependent tables (e.g. for testing the impact of
! parity-dependence).
!
    if (ldmod >= 5 .and. ldexist(Zix, Nix, ibar) .and. .not. flagparity) then
      do nex = 1, nendens(Zix, Nix)
        ldtottableP(Zix, Nix, nex, 1, ibar) = 0.5 * (ldtottableP(Zix, Nix, nex, - 1, ibar) + ldtottableP(Zix, Nix, nex, 1, ibar))
        do J = 0, numJ
          ldtable(Zix, Nix, nex, J, 1, ibar) = 0.5 * (ldtable(Zix, Nix, nex, J, - 1, ibar) + ldtable(Zix, Nix, nex, J, 1, ibar))
        enddo
      enddo
    endif
  enddo
  return
end subroutine densitytable
! Copyright A.J. Koning 2021
