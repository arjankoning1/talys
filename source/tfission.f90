subroutine tfission(Zcomp, Ncomp, nex, J2, parity)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Fission transmission coefficients
!
! Author    : Stephane Hilaire and Pascal Romain
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
!   dbl           ! double precision kind
! All global variables
!   numhill       ! maximum number of Hill - Wheeler points
! Variables for level density
!   ldmodel       ! level density model
! Variables for fission
!   flagclass2    ! flag for class2 states in fission
!   flagfispartdamp    ! flag for fission partial damping
!   widthc2       ! width of class2 states
! Variables for numerics
!   transeps      ! absolute limit for transmission coefficient
! Constants
!   hbar           ! Planck's constant / 2.pi in MeV.s
!   twopi          ! 2 * pi
! Variables for excitation energy grid
!   deltaEx       ! excitation energy bin for population arrays
!   Exmax         ! maximum excitation energy for residual nucleus
! Variables for compound nucleus from target
!   dExinc        ! excitation energy bin for mother nucleus
!   Exinc         ! excitation energy of entrance bin
!   Fnorm         ! multiplication factor
! Variables for fission transmission coefficients
!   denfis        ! fission level density
!   gamfis        ! fission width
!   rhofisA       ! integrated level density corresponding to tfisA
!   taufis        ! fission lifetime
!   tfis          ! fission transmission coefficient for Hill - Wheeler magnitude
!   tfisA         ! transmission coefficient for Hill - Wheeler magnitude
!   tfisdown      ! fission transmission coefficients
!   tfisup        ! fission transmission coefficients
! Variables for fission parameters
!   efisc2rot     ! energy of class2 rotational transition states
!   Emaxclass2    ! maximum energy for class2 states
!   jfisc2rot     ! spin of class2 rotational transition states
!   nfisbar       ! number of fission barrier parameters
!   nfisc2rot     ! number of class2 rotational transition states for barrier
!   pfisc2rot     ! parity of class2 rotational transition states
!
! *** Declaration of local data
!
  implicit none
  integer           :: ic2      ! help variable
  integer           :: ihill    ! counter for Hill-Wheeler magnitude
  integer           :: iloop    ! loop counter
  integer           :: J        ! spin of level
  integer           :: J2       ! 2 * J
  integer           :: jc2      ! help variable
  integer           :: Ncomp    ! neutron number index for compound nucleus
  integer           :: nex      ! excitation energy bin of compound nucleus
  integer           :: parity   ! parity
  integer           :: pc2      ! help variable
  integer           :: Zcomp    ! proton number index for compound nucleus
  real(sgl)         :: boost    ! energy boost
  real(sgl)         :: boostmax ! maximum energy boost
  real(sgl)         :: damper   ! energy damping function
  real(sgl)         :: damper1  ! energy damping function
  real(sgl)         :: damper2  ! energy damping function
  real(sgl)         :: diffnrj  ! energy difference
  real(sgl)         :: ec2      ! help variable
  real(sgl)         :: Ecut     ! help variable
  real(sgl)         :: Ecut1    ! cutoff energy
  real(sgl)         :: Ecut2    ! cutoff energy
  real(sgl)         :: Eex      ! excitation energy
  real(sgl)         :: expo     ! help variable
  real(sgl)         :: term1    ! help variable
  real(sgl)         :: term11   ! help variable
  real(sgl)         :: term12   ! help variable
  real(sgl)         :: term2    ! help variable
  real(sgl)         :: term21   ! help variable
  real(sgl)         :: term22   ! help variable
  real(sgl)         :: twkbtransint
  real(sgl)         :: wo2      ! help variable
  real(sgl)         :: wo2damp  ! energy damping function
  real(sgl)         :: wo2damp1 ! energy damping function
  real(sgl)         :: wo2damp2 ! energy damping function
  real(dbl)         :: addnrj   ! added energy
  real(dbl)         :: density  ! level density
  real(dbl)         :: Rnorm    ! help variable
  real(dbl)         :: rnfb1    ! help variable
  real(dbl)         :: rnfb2    ! help variable
  real(dbl)         :: rnfb3    ! help variable
  real(dbl)         :: sumt2    ! help variable
  real(dbl)         :: sumt3    ! help variable
  real(dbl)         :: ta12     ! help variable
  real(dbl)         :: ta13     ! help variable
  real(dbl)         :: ta23     ! help variable
  real(dbl)         :: ta32     ! help variable
  real(dbl)         :: tgam2    ! help variable
  real(dbl)         :: tgam3    ! help variable
  real(dbl)         :: tf       ! help variable
  real(dbl)         :: tf12     ! help variable
  real(dbl)         :: tfb1     ! help variable
  real(dbl)         :: tfb2     ! help variable
  real(dbl)         :: tfb3     ! help variable
  real(dbl)         :: tdir     ! help variable
  real(dbl)         :: tdir12   ! help variable
  real(dbl)         :: tdir13   ! help variable
  real(dbl)         :: tdir21   ! help variable
  real(dbl)         :: tdir23   ! help variable
  real(dbl)         :: tfii     ! help variable
  real(dbl)         :: tfiii    ! help variable
  real(dbl)         :: ti2      ! help variable
  real(dbl)         :: ti3      ! help variable
  real(dbl)         :: trans    ! help variable
  real(dbl)         :: trans2   ! help variable
  real(dbl)         :: trans3   ! help variable
  real(dbl)         :: tsum123  ! help variable
  external twkbtransint
!
! ********** Calculation of fission transmission coefficients **********
!
! The fission transmission coefficients decrease very rapidly with excitation energy.
! Therefore, we calculate them at the endpoints and at the middle of each excitation energy bin.
! With this information, we can do logarithmic integration in subroutine compound.
!
! J and parity are in loops outside this subroutine
!
  J = J2 / 2
  dExinc = deltaEx(Zcomp, Ncomp, nex)
  do iloop = 1, 3
    tf = 0.
    if (iloop == 1) then
      Eex = max(Exinc - 0.5 * dExinc, 0.)
      tfisdown(J, parity) = 0.
    endif
    if (iloop == 2) then
      Eex = Exinc
      tfis(J, parity) = 0.
      gamfis(J, parity) = 0.
      taufis(J, parity) = 0.
      denfis(J, parity) = 0.
      do ihill = 0, numhill
        tfisA(J, parity, ihill) = 0.
        rhofisA(J, parity, ihill) = 1.
      enddo
    endif
    if (iloop == 3) then
      Eex = min(Exinc + 0.5 * dExinc, Exmax(Zcomp, Ncomp))
      tfisup(J, parity) = 0.
    endif
!
! 1. One barrier
!
! t1barrier : subroutine for fission transmission coefficient for one barrier
!
    if (nfisbar(Zcomp, Ncomp) == 1) then
      call t1barrier(Zcomp, Ncomp, J2, parity, 1, tfb1, rnfb1, Eex, iloop)
      tf = tfb1
    endif
!
! 2. Two barriers
!
    if (nfisbar(Zcomp, Ncomp) == 2) then
      if (flagfispartdamp) call tdirbarrier(Zcomp, Ncomp, J2, parity, 1, 2, tdir, rnfb1, Eex)
      call t1barrier(Zcomp, Ncomp, J2, parity, 1, tfb1, rnfb1, Eex, iloop)
      if (tfb1 >= transeps) then
        call t1barrier(Zcomp, Ncomp, J2, parity, 2, tfb2, rnfb2, Eex, iloop)
        if (tfb2 >= transeps) then
          if (flagfispartdamp) then
            trans = Twkbtransint(Eex, 1, Zcomp, Ncomp)
            tf = tfb1 * tfb2 / (tfb1 + tfb2) * trans + tdir * (1. - trans)
          else
            tf = tfb1 * tfb2 / (tfb1 + tfb2)
          endif
!
! ****************** Special treatment for class2 states ***************
!
          Ecut = Emaxclass2(Zcomp, Ncomp, 1) + 0.5 * widthc2(Zcomp, Ncomp, 1)
          if (flagclass2 .and. (Eex <= Ecut)) then
            term1 = - Eex + 0.5 * (efisc2rot(Zcomp, Ncomp, 1, 1) + efisc2rot(Zcomp, Ncomp, 1, nfisc2rot(Zcomp, Ncomp, 1)))
            term2 = efisc2rot(Zcomp, Ncomp, 1, nfisc2rot(Zcomp, Ncomp, 1)) - efisc2rot(Zcomp, Ncomp, 1, 1)
            damper = 1.
            if (term2 > 0.) then
              expo = 24. * term1 / term2
              if (abs(expo) <= 80.) damper = 1. / (1. + exp(expo))
            endif
            wo2 = 0.5 * widthc2(Zcomp, Ncomp, 1)
            wo2damp = wo2 * damper
            tfii = 0.
            do ic2 = nfisc2rot(Zcomp, Ncomp, 1), 1, - 1
              ec2 = efisc2rot(Zcomp, Ncomp, 1, ic2)
              diffnrj = abs(Eex - ec2)
              addnrj = diffnrj / wo2damp
              boost = 1. / (1. + addnrj **2)
              if (boost >= 0.25) then
                jc2 = int(2. * jfisc2rot(Zcomp, Ncomp, 1, ic2))
                pc2 = pfisc2rot(Zcomp, Ncomp, 1, ic2)
                if ((jc2 == J2) .and. (pc2 == parity)) then
                  boostmax = 4. / (tfb1 + tfb2)
                  boost = boostmax * boost
                  tfii = tfii + tf * boost
                endif
              endif
            enddo
            if (tfii > 0.) tf = tfii
          endif
        endif
      endif
    endif
!
! 3. Three barriers
!
    if (nfisbar(Zcomp, Ncomp) == 3) then
      call t1barrier(Zcomp, Ncomp, J2, parity, 1, tfb1, rnfb1, Eex, iloop)
      if (tfb1 >= transeps) then
        call t1barrier(Zcomp, Ncomp, J2, parity, 2, tfb2, rnfb2, Eex, iloop)
        if (tfb2 >= transeps) then
          call t1barrier(Zcomp, Ncomp, J2, parity, 3, tfb3, rnfb3, Eex, iloop)
          if (tfb3 >= transeps) then
          if (flagfispartdamp) then
            call tdirbarrier(Zcomp, Ncomp, J2, parity, 1, 2, tdir12, rnfb1, Eex)
            call tdirbarrier(Zcomp, Ncomp, J2, parity, 2, 3, tdir23, rnfb1, Eex)
            call tdirbarrier(Zcomp, Ncomp, J2, parity, 1, 3, tdir13, rnfb1, Eex)
            trans2 = Twkbtransint(Eex, 1, Zcomp, Ncomp)
            trans3 = Twkbtransint(Eex, 2, Zcomp, Ncomp)
            tdir21 = (1 - trans2) * tdir12
            tdir12 = (1 - trans2) * tdir12
            tdir23 = (1 - trans3) * tdir23
            tdir13 = (1 - trans2) * (1 - trans3) * tdir13
            Ta12 = tfb1 * trans2
            Ta23 = tfb2 * trans3
            Ta13 = trans3 * tdir12
            Ta32 = trans2 * tfb2
            tgam2 = 0.
            tgam3 = 0.
            sumT2 = tfb1 + tdir23 + Ta23 + Tgam2
            sumT3 = tdir21 + tfb3 + Ta32 + Tgam3
            Ti2 = Ta12 * (tdir23 / sumT2 + Ta23 * tfb3 / (sumT2 * sumT3))
            Ti3 = Ta13 * (tfb3 / sumT3 + Ta32 * tdir23 / (sumT2 * sumT3))
            if ( abs((Ta23 * Ta32) / (sumT2 * sumT3) - 1.)  <=  1.e-8 )  then
              Rnorm =  sumT2 * sumT3 / ( (tfb1 + tdir23 + Tgam2) * sumT3 + Ta23 * (tdir21 + tfb3 + Tgam3) )
            else
              Rnorm = 1. / (1. - (Ta23 * Ta32) / (sumT2 * sumT3))
            endif
            tf = tdir13 + Rnorm * (Ti2 + Ti3)
          else
            tf12 = tfb1 * tfb2 / (tfb1 + tfb2)
            tsum123 = tf12 + tfb3
            tf = tf12 * tfb3 / tsum123
          endif
!
! *********** Special treatment for class2 and class3 states ***********
!
            Ecut1 = Emaxclass2(Zcomp, Ncomp, 1) + 0.5 * widthc2(Zcomp, Ncomp, 1)
            Ecut2 = Emaxclass2(Zcomp, Ncomp, 2) + 0.5 * widthc2(Zcomp, Ncomp, 2)
            Ecut = max(Ecut1, Ecut2)
            if (flagclass2 .and. (Eex <= Ecut)) then
              term11 = - Eex + 0.5 * (efisc2rot(Zcomp, Ncomp, 1, 1) + efisc2rot(Zcomp, Ncomp, 1, nfisc2rot(Zcomp, Ncomp, 1)))
              term21 = efisc2rot(Zcomp, Ncomp, 1, nfisc2rot(Zcomp, Ncomp, 1)) - efisc2rot(Zcomp, Ncomp, 1, 1)
              term12 = - Eex + 0.5 * (efisc2rot(Zcomp, Ncomp, 2, 1) + efisc2rot(Zcomp, Ncomp, 2, nfisc2rot(Zcomp, Ncomp, 2)))
              term22 = efisc2rot(Zcomp, Ncomp, 2, nfisc2rot(Zcomp, Ncomp, 2)) - efisc2rot(Zcomp, Ncomp, 2, 1)
              damper1 = 1.
              if (term21 > 0.) then
                expo = 24. * term11 / term21
                if (abs(expo) <= 80.) damper1 = 1. / (1. + exp(expo))
              endif
              damper2 = 1.
              if (term22 > 0.) then
                expo = 24. * term12 / term22
                if (abs(expo) <= 80.) damper2 = 1. / (1. + exp(expo))
              endif
              wo2 = 0.5 * widthc2(Zcomp, Ncomp, 1)
              wo2damp1 = wo2 * damper1
              tfii = 0.
              do ic2 = nfisc2rot(Zcomp, Ncomp, 1), 1, - 1
                ec2 = efisc2rot(Zcomp, Ncomp, 1, ic2)
                diffnrj = abs(Eex - ec2)
                addnrj = diffnrj / wo2damp1
                boost = 1. / (1. + addnrj **2)
                if (boost >= 0.25) then
                  jc2 = int(2. * jfisc2rot(Zcomp, Ncomp, 1, ic2))
                  pc2 = pfisc2rot(Zcomp, Ncomp, 1, ic2)
                  if ((jc2 == J2) .and. (pc2 == parity)) then
                    boostmax = 4. / (tfb1 + tfb2)
                    boost = boostmax * boost
                    tfii = tfii + tfb1 * tfb2 / (tfb1 + tfb2) * boost
                  endif
                endif
              enddo
              if (tfii > 0.) then
                 tf12 = tfii
                 tsum123 = tf12 + tfb3
                 tf = tf12 * tfb3 / tsum123
              endif
              wo2 = 0.5 * widthc2(Zcomp, Ncomp, 2)
              wo2damp2 = wo2 * damper2
              tfiii = 0.
              do ic2 = nfisc2rot(Zcomp, Ncomp, 2), 1, - 1
                ec2 = efisc2rot(Zcomp, Ncomp, 2, ic2)
                diffnrj = abs(Eex - ec2)
                addnrj = diffnrj / wo2damp2
                boost = 1. / (1. + addnrj **2)
                if (boost >= 0.25) then
                  jc2 = int(2. * jfisc2rot(Zcomp, Ncomp, 2, ic2))
                  pc2 = pfisc2rot(Zcomp, Ncomp, 2, ic2)
                  if ((jc2 == J2) .and. (pc2 == parity)) then
                    boostmax = 4. / tsum123
                    boost = boostmax * boost
                    tfiii = tfiii + tf * boost
                  endif
                endif
              enddo
              if (tfiii > 0.) tf = tfiii
            endif
          endif
        endif
      endif
    endif
    tf = tf * Fnorm(-1)
    if (iloop == 1) tfisdown(J, parity) = tf
    if (iloop == 2) tfis(J, parity) = tf
    if (iloop == 3) tfisup(J, parity) = tf
  enddo
!
! Partial fission widths and lifetimes
!
  denfis(J, parity) = density(Zcomp, Ncomp, Exinc, real(J), parity, 0, ldmodel(Zcomp, Ncomp))
  if (denfis(J, parity) > 0.) gamfis(J, parity) = tfis(J, parity) / (twopi * denfis(J, parity))
  if (gamfis(J, parity) > 0.) taufis(J, parity) = hbar / gamfis(J, parity)
  return
end subroutine tfission
! Copyright A.J. Koning 2021
