subroutine racap
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Radiative capture model
!
! Author    : Stephane Goriely, Xu Yi, and Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
!   dbl             ! double precision kind
! All global variables
!   numdensracap    ! number of energy points for tabulated level densities for direct cap.
!   numJ            ! maximum J - value
!   numjlm          ! maximum number of radial points
!   numJph          ! maximum spin for particle - hole states
!   numlev2         ! maximum number of levels
! Variables for input energies
!   nin             ! counter for incident energy
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ztarget         ! charge number of target nucleus
! Variables for masses
!   beta2           ! deformation parameter
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for excitation energy grid
!   deltaEx         ! excitation energy bin for population arrays
!   Ex              ! excitation energy
!   maxex           ! maximum excitation energy bin for residual nucleus
! Constants
!   amu             ! atomic mass unit in MeV
!   e2              ! square of elementary charge in MeV.fm
!   hbarc           ! hbar.c in MeV.fm
!   parA            ! mass number of particle
!   parmass         ! mass of particle in a.m.u.
!   parN            ! neutron number of particle
!   parspin         ! spin of particle
!   parZ            ! charge number of particle
!   pi              ! pi
! Variables for levels
!   edis            ! energy of level
!   jdis            ! spin of level
!   parlev          ! parity of level
! Variables for direct capture initialization
!   avncap2         ! real volume diffuseness for JLM
!   chglnegj        ! help variable
!   chglposj        ! help variable
!   ispect          ! model for discrete levels
!   jlmracap2       ! JLM potential for direct capture
!   nlevexpracap    ! number of experimental levels in the final nucleus
!   nlevracap       ! number of levels in the final nucleus
!   racopt          ! OMP for radiative capture
!   rvncap2         ! real volume radius for JLM
!   spectfac        ! spectroscopic factor
!   vncap2          ! real volume depth for JLM
!   xsracap         ! direct radiative capture cross section
!   xsracapEM       ! direct - semidirect radiative capture cross section as
! Variables for direct capture
!   xsracape        ! direct radiative capture cross section
!   xsracapecont    ! direct radiative capture continuum cross section
!   xsracapedisc    ! direct radiative capture discrete cross section
!   xsracappop      ! population cross section for radiative capture
!   xsracappopex    ! population cross section for radiative capture
! Variables for optical model
!   av              ! real volume diffuseness
!   rv              ! real volume radius
!   v               ! real volume depth
! Variables for JLM
!   normjlm         ! JLM potential normalization factors
!   potjlm          ! JLM potential depth values
!   rhojlmn         ! density for neutrons
!   rhojlmp         ! density for protons
! Variables for masses
!   expmass         ! flag for using experimental nuclear mass if available
!   gsparity        ! ground state parity
!   gsspin          ! ground state spin
!   S               ! separation energy
!   thmass          ! theoretical mass
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                         ! counter
  integer   :: iopt                      ! index
  integer   :: J                         ! spin of level
  integer   :: nex                       ! excitation energy bin of compound nucleus
  integer   :: Nix                       ! neutron number index for residual nucleus
  integer   :: par                       ! parameter
  integer   :: parity                    ! parity
  integer   :: prac(numlev2)             ! parity
  integer   :: Zix                       ! charge number index for residual nucleus
  real(sgl) :: avncap1                   ! diffuseness
  real(sgl) :: erac(numlev2)             ! ebergy
  real(sgl) :: exfin(numlev2)            ! energy
  real(sgl) :: jlmracap1(numjlm)         ! JLMB potential for the radial grid: [0.1, 20], step of 0.1
  real(sgl) :: jrac(numlev2)             ! spin
  real(sgl) :: qq                        ! help variable
  real(sgl) :: rvncap1                   ! radius
  real(sgl) :: spfacst(numlev2)          ! spin
  real(sgl) :: vncap1                    ! potential
  real(sgl) :: xsall(3)                  ! cross section
  real(sgl) :: xsp(numlev2, 0:numJph, 2) ! (n,p) cross section
  real(sgl) :: xspex(numlev2)            ! help variable
  real(dbl) :: pZ                        ! product of charges
!
! ********************** General initializations ***********************
!
! Calculation the initial wood saxon potential parameter
!
  xsall = 0.
  xsp = 0.
  xspex = 0.
  exfin = 0.
  iopt = racopt
  Zix = parZ(k0)
  Nix = parN(k0)
  call optical(Zix, Nix, k0, Einc)
  vncap1 = v
  rvncap1 = rv
  avncap1 = av
!
! iopt=1 - Woods-Saxon potential
! iopt=2 - No.2, no use
! iopt=3 - JLMB, Jeukenne et al. (1977) and Burge PRC, 1998
! iopt=4 - Folding potential M3Y (now not available)
!
! Calculation the initial state JLMB potential
!
  pZ = dble(real(Zix) * ZTarget)
  call mom(Zix, Nix, pZ, dble(Einc))
  do i = 1, numjlm
    jlmracap1(i) = normjlm(Zix, Nix, 1) * (potjlm(Zix, Nix, i, 1) + &
      normjlm(Zix, Nix, 3) * (rhojlmn(Zix, Nix, i, 1) + &
      rhojlmp(Zix, Nix, i, 1)) / (rhojlmn(Zix, Nix, i, 1) - &
      rhojlmp(Zix, Nix, i, 1)) * potjlm(Zix, Nix, i, 3))
  enddo
!
! Set the experimental level position, spin and parity.
!
  spfacst = 0.
  jrac = 0.
  prac = 0
  erac = 0.
  spfacst = 0.
  do i = 1, nlevexpracap
    jrac(i) = jdis(0, 0, i)
    prac(i) = parlev(0, 0, i)
    erac(i) = edis(0, 0, i)
    spfacst(i) = spectfac(0, 0, i - 1)
  enddo
  do i = nlevexpracap, maxex(0, 0) + 1
    spfacst(i) = spectfac(0, 0, i - 1)
  enddo
  qq = S(0,0,k0)
!
! end of initial setting and begin to racap calculation.
!
  call racapcalc(k0, Einc, Ztarget, Atarget, beta2(Zix, Nix, 0), jdis(Zix, Nix, 0), parlev(Zix, Nix, 0), gsspin(Zix, Nix), &
    gsparity(Zix, Nix), expmass(Zix, Nix), thmass(Zix, Nix), &
    Zix, parA(k0), parspin(k0), parmass(k0), jdis(0, 0, 0), parlev(0, 0, 0), &
    gsspin(0, 0), gsparity(0, 0), expmass(0, 0), thmass(0, 0), qq, numlev2, nlevexpracap, erac, jrac, prac, exfin, &
    numjlm, jlmracap1, jlmracap2, numJph, numdensracap, chglposj, chglnegj, &
    vncap1, rvncap1, avncap1, vncap2, rvncap2, avncap2, xsall, xsracape, &
    xspex, xsp, iopt, pi, e2, amu, hbarc, nlevracap(0, 0), spfacst, ispect)
!
! Racap output cross section,transfer the unit of the cross section to millibarns.
!
  xsracape = xsracape * 1000.0
  xsracap(nin) = xsracape
  xsracapEM(nin, 1, 1) = xsall(1) * 1000.0
  xsracapEM(nin, 1, 2) = xsall(2) * 1000.0
  xsracapEM(nin, 0, 1) = xsall(3) * 1000.0
  xsracapedisc = 0.
  xsracapecont = 0.
  do i = 1, nlevracap(0, 0)
    if (i <= nlevexpracap) then
!   spectfac(0,0,i-1)=spfacst(i)
      xsracapedisc = xsracapedisc + xspex(i) * 1000.0
      xsracappopex(i - 1) = xspex(i) * 1000.0
      do J = 0, numJph
      do parity = - 1, 1, 2
        if (parity ==  - 1) par = 1
        if (parity == 1) par = 2
        xsracappop(i - 1, J, parity) = xsp(i, J, par) * 1000.0
      enddo
      enddo
    else
      xsracapecont = xsracapecont + xspex(i) * 1000.0
      do nex = nlevexpracap, maxex(0, 0)
!   spectfac(0,0,nex)=spfacst(i)
        if (exfin(i) < Ex(0, 0, nex) + deltaEx(0, 0, nex) / 2..and. &
          exfin(i) >= Ex(0, 0, nex) - deltaEx(0, 0, nex) / 2.) then
          xsracappopex(nex) = xsracappopex(nex) + xspex(i) * 1000.0
          do J = 0, numJph
          do parity = - 1, 1, 2
            if (parity ==  - 1) par = 1
            if (parity == 1) par = 2
            xsracappop(nex, J, parity) = xsracappop(nex, J, parity) + xsp(i, J, par) * 1000.0
          enddo
          enddo
        endif
      enddo
    endif
  enddo
  return
end subroutine racap
! Copyright A.J. Koning 2021
