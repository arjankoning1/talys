subroutine densitypar(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Level density parameters
!
! Author    : Arjan Koning and Stephane Hilaire
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
!   sgl               ! single precision kind
!   dbl               ! double precision kind
! Variables for preequilibrium
!   g                 ! single - particle level density parameter
!   gadjust           ! adjustable factor for single - particle particle-hole states
!   gn                ! single - particle neutron level density parameter
!   gnadjust          ! adjustable factor for single - particle proton par
!   gp                ! single - particle proton level density parameter
!   gpadjust          ! adjustable factor for single - particle neutron pa
!   Kph               ! constant for single - particle level density par.
! Variables for fission
!   axtype            ! type of axiality of barrier
!   fbarrier          ! height of fission barrier
!   flagfission       ! flag for fission
! Variables for discrete levels
!   nlev              ! number of levels for nucleus
! Variables for level density
!   aadjust           ! adjustable factor for level density parameter
!   alev              ! level density parameter
!   alimit            ! asymptotic level density parameter
!   alphald           ! alpha - constant for asymptotic level density para
!   betald            ! beta - constant for asymptotic level density param
!   ctable            ! constant to adjust tabulated level densities
!   ctableadjust      ! correction to adjust tabulated level densities
!   deltaW            ! shell correction in nuclear mass
!   flagasys          ! flag for all level density parameters a fr
!   flagcol           ! flag for collective enhancement of level density
!   flagcolldamp      ! flag for damping of coll. effects in eff. level density (without explici
!   gammald           ! gamma - constant for asymptotic level density para
!   gammashell1       ! gamma - constant for asymptotic level density para
!   gammashell2       ! gamma - constant for asymptotic level density para
!   ldmodel           ! level density model
!   Nlow              ! lowest discrete level for temperature matching
!   Ntop              ! highest discrete level for temperature matching
!   pair              ! pairing energy
!   pairconstant      ! constant for pairing energy systematics
!   Pshift            ! adjustable pairing shift
!   Pshiftadjust      ! adjustable correction to pairing shift
!   Pshiftconstant    ! global constant for pairing shift
!   ptable            ! constant to adjust tabulated level densities
!   ptableadjust      ! correction to adjust tabulated level densities
!   shellmodel        ! model for shell correction energies
! Variables for nuclides
!   AA                ! mass number of residual nucleus
!   NN                ! neutron number of residual nucleus
!   ZZ                ! charge number of residual nucleus
! Variables for files
!   path              ! directory containing files to be read
! Constants
!   amu               ! atomic mass unit in MeV
!   nuc               ! symbol of nucleus
!   onethird          ! 1 / 3
!   pi                ! pi
!   pi2               ! pi **2
!   twothird          ! 2 / 3
! Variables for levels
!   edis              ! energy of level
!   jdis              ! spin of level
! Variables for level density
!   aldcrit           ! critical level density parameter
!   Dcrit             ! critical determinant
!   delta             ! energy shift
!   delta0            ! systematical pairing energy
!   Econd             ! condensation energy
!   Ediscrete         ! energy of middle of discrete level region
!   ldparexist        ! flag for existence of tabulated level density
!   Nlast             ! last discrete level
!   Scrit             ! critical entropy
!   scutoffdisc       ! spin cutoff factor for discrete level region
!   Tcrit             ! critical temperature
!   Ucrit             ! critical U
! Variables for fission parameters
!   efistrrot         ! energy of rotational transition states
!   jfistrrot         ! spin of rotational transition states
!   nfisbar           ! number of fission barrier parameters
!   nfistrrot         ! number of rotational transition states for barr
! Variables for masses
!   nucmass           ! mass of nucleus
!   S                 ! separation energy
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: inpalev       ! logical to determine existence of input value for a
  logical           :: inpalimit     ! logical to determine existence of input value for alimit
  logical           :: inpdeltaW     ! logical to determine existence of input value for deltaW
  logical           :: inpgammald    ! logical to determine existence of input value for gammald
  logical           :: lexist        ! logical to determine existence
  character(len=5)  :: denchar       ! string for level density file
  character(len=22) :: denformat     ! format specifier
  character(len=132):: denfile       ! level density parameter file
  integer           :: A             ! mass number of target nucleus
  integer           :: i             ! counter
  integer           :: ia            ! mass number from abundance table
  integer           :: ibar          ! fission barrier
  integer           :: iloop         ! counter
  integer           :: imax          ! help variable
  integer           :: imin          ! help variable
  integer           :: istat         ! logical for file access
  integer           :: ldmod         ! level density model
  integer           :: N             ! neutron number of residual nucleus
  integer           :: Nix           ! neutron number index for residual nucleus
  integer           :: Nlow0         ! help variable
  integer           :: Ntop0         ! highest discrete level for temperature matching
  integer           :: oddN          ! help variable
  integer           :: oddZ          ! help variable
  integer           :: Z             ! charge number of target nucleus
  integer           :: Zix           ! charge number index for residual nucleus
  real(sgl)         :: ald           ! level density parameter
  real(sgl)         :: ald0          ! level density parameter
  real(sgl)         :: argum         ! help variable
  real(sgl)         :: denom         ! help variable
  real(sgl)         :: difprev       ! difference with previous result
  real(sgl)         :: expo          ! help variable
  real(sgl)         :: factor        ! multiplication factor
  real(sgl)         :: fU            ! help variable
  real(sgl)         :: pshift0       ! adjustable pairing shift
  real(sgl)         :: rj            ! help variable
  real(sgl)         :: scutoffsys    ! spin cutoff factor for discrete level from systematics
  real(sgl)         :: sd            ! spin cutoff factor for discrete level region
  real(sgl)         :: sigsum        ! help variable
  real(sgl)         :: Spair         ! help variable
  real(dbl)         :: mldm          ! liquid drop mass
  real(dbl)         :: mliquid1      ! function for liquid drop mass (Myers-Swiatecki)
  real(dbl)         :: mliquid2      ! function for liquid drop mass (Goriely)
!
! *************************** Initialization ***************************
!
! ldmodel 1: Gilbert and Cameron
! ldmodel 2: Back-shifted Fermi gas
! ldmodel 3: Superfluid model
! ldmodel 4: Statistical HFB model (Goriely)
! ldmodel 5: Combinatorial HFB model (Hilaire and Goriely)
! ldmodel 6: Combinatorial HFB model - Gogny force (Hilaire and Goriely)
! ldmodel 7: BSKG3 HFB model - Skyrme force (Ryssens and Goriely)
! ldmodel 8: QRPA model 
!
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  ldmod = ldmodel(Zix, Nix)
!
! *********** Read values from level density parameter file ************
!
! Level density parameters from the table can always be overruled by values given in the input file.
! With flagasys, all experimental level density parameters a from the table can be overruled by the systematics.
! We allow a maximum of Ntop=50
!
  denchar = trim(nuc(Z))//'.ld'
  if (ldmod == 1) denfile = trim(path)//'density/ground/ctm/'//denchar
  if (ldmod == 2) denfile = trim(path)//'density/ground/bfm/'//denchar
  if (ldmod == 3) denfile = trim(path)//'density/ground/gsm/'//denchar
  if (ldmod == 4) denfile = trim(path)//'density/ground/goriely/'//denchar
  if (ldmod == 5) denfile = trim(path)//'density/ground/hilaire/'//denchar
  if (ldmod == 6) denfile = trim(path)//'density/ground/hilaireD1M/'//denchar
  if (ldmod == 7) denfile = trim(path)//'density/ground/bskg3/'//denchar
  if (ldmod == 8) denfile = trim(path)//'density/ground/qrpa/'//denchar
  inquire (file = denfile, exist = lexist)
  if (lexist) then
    if (flagcol(Zix, Nix) .and. ldmod <= 3) then
      denformat='(4x,i4,32x,2i4,2f12.5)'
    else
      denformat='(4x,3i4,2f12.5)'
    endif
    open (unit = 2, file = denfile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(denfile, istat)
    do
      read(2, fmt = denformat, iostat = istat) ia, Nlow0, Ntop0, ald0, pshift0
      if (istat == -1) exit
      if (istat /= 0) call read_error(denfile, istat)
      if (A == ia) then
        ldparexist(Zix, Nix) = .true.
        if (Nlow(Zix, Nix, 0) ==  - 1) Nlow(Zix, Nix, 0) = Nlow0
        if (Ntop(Zix, Nix, 0) ==  - 1) Ntop(Zix, Nix, 0) = min(Ntop0, 50)
        if ( .not. flagasys) then
          if (ldmod <= 3) then
            if (alev(Zix, Nix) == 0.) alev(Zix, Nix) = aadjust(Zix, Nix) * ald0
            do ibar = 0, nfisbar(Zix, Nix)
              if (Pshift(Zix, Nix, ibar) == 1.e-20) Pshift(Zix, Nix, ibar) = pshift0 + Pshiftadjust(Zix, Nix, ibar)
            enddo
          else
            if (ctable(Zix,Nix,0) == 1.e-20) ctable(Zix, Nix, 0) = ald0
            if (ptable(Zix,Nix,0) == 1.e-20) ptable(Zix, Nix, 0) = pshift0
          endif
          ctable(Zix, Nix, 0) = ctable(Zix, Nix, 0) + ctableadjust(Zix, Nix, 0)
          ptable(Zix, Nix, 0) = ptable(Zix, Nix, 0) + ptableadjust(Zix, Nix, 0)
        endif
      endif
    enddo
    close (unit = 2)
  endif
!
! Matching levels
!
  do ibar = 0, nfisbar(Zix, Nix)
    if (ibar == 0) then
      if (Ntop(Zix, Nix, ibar) == -1) then
        Nlast(Zix, Nix, ibar) = nlev(Zix, Nix)
      else
        Nlast(Zix, Nix, ibar) = min(Ntop(Zix, Nix, ibar), nlev(Zix, Nix))
      endif
    else
      Nlast(Zix, Nix, ibar) = max(nfistrrot(Zix, Nix, ibar), 1)
    endif
    if (Ntop(Zix, Nix, ibar) ==  -1) Ntop(Zix, Nix, ibar) = Nlast(Zix, Nix, ibar)
    if (Nlow(Zix, Nix, ibar) ==  -1) Nlow(Zix, Nix, ibar) = 2
    if (Ntop(Zix, Nix, ibar) <= 2) Nlow(Zix, Nix, ibar) = 0
  enddo
!
! Determine spin cut-off parameter for discrete level region
!
! First assign the systematics value, then overrule in case of enough discrete level info.
!
  scutoffsys = (0.83 * (A **0.26)) **2
  do ibar = 0, nfisbar(Zix, Nix)
    scutoffdisc(Zix, Nix, ibar) = scutoffsys
    if (ldparexist(Zix, Nix)) then
      imax = Ntop(Zix, Nix, ibar)
      if (ibar == 0) then
        imin = Nlow(Zix, Nix, 0)
        Ediscrete(Zix, Nix, 0) = 0.5 * (edis(Zix, Nix, imin) + edis(Zix, Nix, imax))
      else
        imin = 1
        Ediscrete(Zix, Nix, ibar) = 0.5 * (efistrrot(Zix, Nix, ibar, imin) + efistrrot(Zix, Nix, ibar, imax))
      endif
      sigsum = 0.
      denom = 0.
      do i = imin, imax
        if (ibar == 0) then
          rj = jdis(Zix, Nix, i)
        else
          rj = jfistrrot(Zix, Nix, ibar, i)
        endif
        sigsum = sigsum + rj * (rj + 1) * (2 * rj + 1)
        denom = denom + 2 * rj + 1
      enddo
      sd = 0.
      if (denom /= 0.) sd = sigsum / (3. * denom)
      if (sd > scutoffsys / 3..and.sd < scutoffsys * 3.) scutoffdisc(Zix, Nix, ibar) = sd
    endif
  enddo
!
! Check input of various level density parameters
!
! mliquid1   : function for liquid drop mass (Myers-Swiatecki)
! mliquid2   : function for liquid drop mass (Goriely)
!
! shellmodel 1: Myers-Swiatecki
! shellmodel 2: Goriely
!
  if (alev(Zix, Nix) == 0.) then
    inpalev = .false.
  else
    inpalev = .true.
    if (ldmod == 3 .and. alimit(Zix, Nix) == 0.) alimit(Zix, Nix) = alev(Zix, Nix)
  endif
  inpdeltaW = .true.
  if (deltaW(Zix, Nix, 0) == 0.) then
    inpdeltaW = .false.
    if (shellmodel == 1) then
      mldm = mliquid1(Z, A)
    else
      mldm = mliquid2(Z, A)
    endif
    deltaW(Zix, Nix, 0) = real((nucmass(Zix, Nix) - mldm) * amu)
  endif
  inpalimit = .true.
  if (alimit(Zix, Nix) == 0.) then
    inpalimit = .false.
    alimit(Zix, Nix) = alphald(Zix, Nix) * A + betald(Zix, Nix) * (A **twothird)
  endif
  inpgammald = .true.
  if (gammald(Zix, Nix) ==  - 1.) then
    inpgammald = .false.
    gammald(Zix, Nix) = gammashell1(Zix, Nix) / (A **onethird) + gammashell2
  endif
!
! The Ignatyuk formula implies that alev, deltaW, gammald and alimit can not all be given as input.
! In that case we re-determine alev.
!
  if (inpalev .and. inpdeltaW .and. inpalimit .and. inpgammald) then
    inpalev = .false.
    alev(Zix, Nix) = 0.
  endif
!
! Pairing corrections
!
  oddZ = mod(Z, 2)
  oddN = mod(N, 2)
  delta0(Zix, Nix) = pairconstant / sqrt(real(A))
!
! Defaults
!
  if (pair(Zix, Nix) == 1.e-20) then
    if (ldmod == 3) then
      pair(Zix, Nix) = (oddZ + oddN) * delta0(Zix, Nix)
    else
      if (ldmod == 2) then
        pair(Zix, Nix) = (1. - oddZ - oddN) * delta0(Zix, Nix)
      else
        pair(Zix, Nix) = (2. - oddZ - oddN) * delta0(Zix, Nix)
      endif
    endif
  endif
  do ibar = 0, nfisbar(Zix, Nix)
    if (Pshift(Zix, Nix, ibar) == 1.e-20) Pshift(Zix, Nix, ibar) = Pshiftconstant(Zix, Nix) + Pshiftadjust(Zix, Nix, ibar)
  enddo
!
! ************************** Fission ***********************************
!
! Determine deltaW on the fission barrier.
! Note that this overrules input parameters for deltaW.
!
  if (flagfission) then
    do ibar = 1, nfisbar(Zix, Nix)
      if (deltaW(Zix, Nix, ibar) == 0.) then
        if (flagcolldamp) then
          deltaW(Zix, Nix, ibar) = abs(deltaW(Zix, Nix, 0)) * twothird
        else
          if (ibar == 1) then
            if (axtype(Zix, Nix, 1) == 1) then
              deltaW(Zix, Nix, ibar) = 1.5
            else
              deltaW(Zix, Nix, ibar) = 2.5
            endif
          else
            deltaW(Zix, Nix, ibar) = 0.6
          endif
        endif
      endif
    enddo
    if (nfisbar(Zix, Nix) == 1 .and. fbarrier(Zix, Nix, 1) == 0.) deltaW(Zix, Nix, 1) = deltaW(Zix, Nix, 2)
  endif
!
! Generalized superfluid model. The critical functions are calculated here.
!
  if (ldmod == 3) then
    Tcrit(Zix, Nix) = 0.567 * delta0(Zix, Nix)
    ald = alimit(Zix, Nix)
    difprev = 0.
    do ibar = 0, nfisbar(Zix, Nix)
      iloop = 0
      do
        factor = (1. - exp( - gammald(Zix, Nix) * ald * Tcrit(Zix, Nix) **2)) / (ald * (Tcrit(Zix, Nix) **2))
        aldcrit(Zix, Nix, ibar) = alimit(Zix, Nix) * (1. + deltaW(Zix, Nix, ibar) * factor)
        if (abs(aldcrit(Zix, Nix, ibar) - ald) > 0.001 .and. abs(aldcrit(Zix, Nix, ibar) - ald) /= difprev .and. iloop <= 1000) then
          difprev = abs(aldcrit(Zix, Nix, ibar) - ald)
          ald = aldcrit(Zix, Nix, ibar)
          iloop = iloop + 1
          if (ald <= 1.) exit
        else
          exit
        endif
      enddo
      if (aldcrit(Zix, Nix, ibar) < alimit(Zix, Nix) / 3.) then
        expo = min( - gammald(Zix, Nix) * S(Zix, Nix, 1), 80.)
        fU = 1. - exp(expo)
        factor = 1. + fU * deltaW(Zix, Nix, ibar) / S(Zix, Nix, 1)
        aldcrit(Zix, Nix, ibar) = max(alimit(Zix, Nix) * factor, 1.)
      endif
      Econd(Zix, Nix, ibar) = 1.5 / pi2 * aldcrit(Zix, Nix, ibar) * delta0(Zix, Nix) **2
      Ucrit(Zix, Nix, ibar) = aldcrit(Zix, Nix, ibar) * Tcrit(Zix, Nix) **2 + Econd(Zix, Nix, ibar)
      Scrit(Zix, Nix, ibar) = 2. * aldcrit(Zix, Nix, ibar) * Tcrit(Zix, Nix)
      Dcrit(Zix, Nix, ibar) = 144. / pi * (aldcrit(Zix, Nix, ibar) **3) * (Tcrit(Zix, Nix) **5)
      delta(Zix, Nix, ibar) = Econd(Zix, Nix, ibar) - pair(Zix, Nix) - Pshift(Zix, Nix, ibar)
    enddo
  else
!
! Constant temperature and back-shifted Fermi gas model
!
    do ibar = 0, nfisbar(Zix, Nix)
      delta(Zix, Nix, ibar) = pair(Zix, Nix) + Pshift(Zix, Nix, ibar)
    enddo
  endif
!
! 1. If no experimental level density parameter is available, i.e. as determined from the neutron resonance spacing, use the
!    Ignatyuk formula to derive the level density parameter at the separation energy.
!
  Spair = S(Zix, Nix, 1) - delta(Zix, Nix, 0)
  Spair = max(Spair, 1.)
  if ( .not. inpalev) then
    fU = 1. - exp( - gammald(Zix, Nix) * Spair)
    factor = 1. + fU * deltaW(Zix, Nix, 0) / Spair
    alev(Zix, Nix) = aadjust(Zix, Nix) * alimit(Zix, Nix) * factor
    alev(Zix, Nix) = max(alev(Zix, Nix), 1.)
  else
!
! 2. If an experimental level density parameter is available, then we impose the extra boundary boundary condition that it should be
!    equal to the energy dependent level density parameter at the neutron separation energy.
!    There are various possibilities. If alimit is not given as input, we re-adjust it.
!
    if (ldmod /= 3) then
      if ( .not. inpalimit) then
        fU = 1. - exp( - gammald(Zix, Nix) * Spair)
        factor = 1. + fU * deltaW(Zix, Nix, 0) / Spair
        alimit(Zix, Nix) = alev(Zix, Nix) / factor
      else
!
! If both alev and alimit are explicitly given as input, we re-adjust deltaW, provided it is not given as input.
!
        if ( .not. inpdeltaW) then
          fU = 1. - exp( - gammald(Zix, Nix) * Spair)
          factor = alev(Zix, Nix) / alimit(Zix, Nix) - 1.
          deltaW(Zix, Nix, 0) = Spair * factor / fU
        else
!
! Determine gammald if alev, alimit and deltaW are given by input.
!
          argum = 1. - Spair / deltaW(Zix, Nix, 0) * (alev(Zix, Nix) / alimit(Zix, Nix) - 1.)
          if (argum > 0..and.argum < 1.) then
            gammald(Zix, Nix) = - 1. / Spair * log(argum)
          else
!
! If gammald can not be solved or is unphysical (this may happen for certain parameter combinations)
! we re-adjust the shell correction.
!
            fU = 1. - exp( - gammald(Zix, Nix) * Spair)
            factor = alev(Zix, Nix) / alimit(Zix, Nix) - 1.
            deltaW(Zix, Nix, 0) = Spair * factor / fU
          endif
        endif
      endif
    endif
  endif
!
! ************** Single-particle level density parameter g *************
!
! One component
!
  if (g(Zix, Nix) == 0.) g(Zix, Nix) = A / Kph
  g(Zix, Nix) = gadjust(Zix, Nix) * g(Zix, Nix)
!
! Two component
!
  if (gp(Zix, Nix) == 0.) gp(Zix, Nix) = Z / Kph
  if (gn(Zix, Nix) == 0.) gn(Zix, Nix) = N / Kph
  gn(Zix, Nix) = gadjust(Zix, Nix) * gnadjust(Zix, Nix) * gn(Zix, Nix)
  gp(Zix, Nix) = gadjust(Zix, Nix) * gpadjust(Zix, Nix) * gp(Zix, Nix)
  return
end subroutine densitypar
! Copyright A.J. Koning 2021
