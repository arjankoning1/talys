subroutine racapinit
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization of radiative capture model
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
!   numdens         ! number of energy points for tabulated level densities
!   numdensracap    ! number of energy points for tabulated level densities for direct cap.
!   numex           ! maximum number of excitation energies
!   numjlm          ! maximum number of radial points
!   numJph          ! maximum spin for particle - hole states
! Variables for level density
!   alev            ! level density parameter
! Variables for preequilibrium
!   flaggshell      ! flag for energy dependence of single particle level den
!   gn              ! single - particle neutron level density parameter
!   gp              ! single - particle proton level density parameter
! Variables for gamma rays
!   ldmodelracap    ! level density model for direct radiative capture
!   spectfacexp     ! experimental spectroscopic factor
!   spectfacth      ! theoretical spectroscopic factor
! Variables for main input
!   Ainit           ! mass number of initial compound nucleus
!   k0              ! index of incident particle
! Variables for OMP
!   flagjlm         ! flag for using semi - microscopic JLM OMP
! Variables for energies
!   Etotal          ! total energy of compound system (target + projectile)
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   nuc             ! symbol of nucleus
!   parsym          ! symbol of particle
! Variables for files
!   path            ! directory containing files to be read
! Variables for levels
!   edis            ! energy of level
!   jdis            ! spin of level
!   parlev          ! parity of level
! Variables for level density
!   nendens         ! number of energies for level density grid
!   Nlast           ! last discrete level
! Variables for direct capture initialization
!   avncap2         ! real volume diffuseness for JLM
!   chglnegj        ! help variable
!   chglposj        ! help variable
!   edensphjp       ! energy grid of ph spin - and parity - dependent level d
!   ispect          ! model for discrete levels
!   jlmracap2       ! JLM potential for direct capture
!   nlevexpracap    ! number of experimental levels in the final nucleus
!   nlevracap       ! number of levels in the final nucleus
!   phdensjp        ! ph spin - and parity - dependent level density from tab
!   phdenstot       ! total ph level density from table
!   racopt          ! OMP for radiative capture
!   rvncap2         ! real volume radius for JLM
!   spectfac        ! spectroscopic factor
!   vncap2          ! real volume depth for JLM
! Variables for optical model
!   av              ! real volume diffuseness
!   rv              ! real volume radius
!   v               ! real volume depth
! Variables for JLM
!   normjlm         ! JLM potential normalization factors
!   potjlm          ! JLM potential depth values
!   rhojlmn         ! density for neutrons
!   rhojlmp         ! density for protons
!
! *** Declaration of local data
!
  implicit none
  integer            :: acp                                  ! mass number
  integer            :: holen                                ! neutron hole number
  integer            :: holep                                ! proton hole number
  integer            :: i                                    ! counter
  integer            :: i0                                   ! particle type
  integer            :: istat                                ! logical for file access
  integer            :: j                                    ! counter
  integer            :: jlev                                 ! spin of level
  integer            :: ka                                   ! spectroscopic index
  integer            :: parity                               ! parity
  integer            :: phnumracap                           ! number of particle-holes
  integer            :: ptcln                                ! neutron particle number
  integer            :: ptclp                                ! proton particle number
  integer            :: wprty                                ! parity
  integer            :: xprty                                ! parity
  real(sgl)          :: e                                    ! energy
  real(sgl)          :: eopt                                 ! incident energy
  real(sgl)          :: Ephracap                             ! particle-hole energy
  real(sgl)          :: ignatyuk                             ! function for energy dependent level density parameter a
  real(sgl)          :: phdens2                              ! function for two-component particle-hole state density
  real(sgl)          :: racapdamp                            ! shell damping
  real(sgl)          :: sf                                   ! parameter
  real(sgl)          :: sigparleldenn                        ! neutron ph density
  real(sgl)          :: sigparleldenp                        ! proton ph density
  real(sgl)          :: spinf                                ! spin
  real(sgl)          :: surfwellE                            ! energy well depth for surface damping
  logical            :: lexist                               ! logical to determine existence
  logical            :: surfwellgo                           ! flag for surface damping
  integer            :: nexphjp                              ! energy index
  integer            :: nex                                  ! excitation energy bin of compound nucleus
  real(sgl)          :: eb                                   ! help variable
  real(sgl)          :: ee                                   ! energy
  real(sgl)          :: Egridphjp(0:numdens)                 ! energy of particle-hole pair
  real(dbl)          :: density                              ! level density
  real(dbl)          :: ldbneg                               ! lower level density
  real(dbl)          :: ldbpos                               ! upper level density
  real(dbl)          :: ldeneg                               ! lower level density
  real(dbl)          :: ldepos                               ! upper level density
  real(dbl)          :: ldtabneg                             ! lower level density
  real(dbl)          :: ldtabpos                             ! upper level density
  real(dbl)          :: lldbneg                              ! lower level density
  real(dbl)          :: lldbpos                              ! upper level density
  real(dbl)          :: lldeneg                              ! lower level density
  real(dbl)          :: lldepos                              ! upper level density
  character(len=132) :: filespec                             ! spectrum file
  real(dbl)          :: ldmdposj(numdensracap, 0:numJph)     ! upper level density
  real(dbl)          :: ldmdnegj(numdensracap, 0:numJph)     ! lower level density
  real(dbl)          :: phmdposj(numdensracap, 0:numJph)     ! upper particle-hole density
  real(dbl)          :: phmdnegj(numdensracap, 0:numJph)     ! lower particle-hole density
  real(dbl)          :: phjpmdposj(numdensracap, 0:numJph)   ! upper particle-hole density
  real(dbl)          :: phjpmdnegj(numdensracap, 0:numJph)   ! lower particle-hole density
  real(dbl)          :: ldmdpos(numdensracap)                ! upper level density
  real(dbl)          :: ldmdneg(numdensracap)                ! lower level density
  real(dbl)          :: phmdpos(numdensracap)                ! upper particle-hole density
  real(dbl)          :: phjpmdpos(numdensracap)              ! upper particle-hole density
  real(dbl)          :: phjpmdneg(numdensracap)              ! lower particle-hole density
  real(sgl)          :: Eldpd(numdensracap)                  ! energy of particle-hole pair
  real(sgl)          :: Sldpd(0:numJph)                      ! spin of particle-hole pair
!
! ********************** General initializations ***********************
!
!   accept the level density and ph level density, with simple calculation and gridlization.
!
!   set the grid of spin, for odd A and even A:
!
  if (mod(Ainit, 2) == 1) then
    do j = 0, numJph
      Sldpd(j) = real(j) + 0.5
      Sldpd(j) = real(Sldpd(j))
    end do
  end if
  if (mod(Ainit, 2) == 0) then
    do j = 0, numJph
      Sldpd(j) = real(j)
      Sldpd(j) = real(Sldpd(j))
    end do
  end if
!
!   set the grid of exitation energy points: 0.125(point1), 0.375(point2)  ...
!
  do i = 1, numdensracap
    Eldpd(i) = 0.125 + (dble(i) - 1.0) * 0.25
    Eldpd(i) = real(Eldpd(i))
  enddo
!
!   choose the NLD model for racap calculation (ldmodelracap=1, 2, or 3)
!
!   edens and nendens need to be read from phdensitytablejp.f. Therefore, if ldmodelracap
!   need to be set as a keyword in the input file, phdensitytablejp.f should be changed,
!   and talys.cmb should be changed accordingly.

  nlevracap = 0
  if (ldmodelracap == 1) then
    phdenstot = 0.
    edensphjp = 0.
    phdensjp = 0.
    call phdensitytablejp(0, 0)
  endif
!
! ldmodelracap may have been changed in phdensitytablejp
!
  if (ldmodelracap == 1) then
    do nex = 1, nendens(0, 0)
      Egridphjp(nex) = real(edensphjp(0, 0, nex))
    enddo
    Egridphjp(0) = 0.0
    do j = 0, numJph
      do i = 1, numdensracap
        if (Eldpd(i) < Egridphjp(nendens(0, 0))) then
          call locate(Egridphjp, 0, nendens(0, 0), Eldpd(i), nexphjp)
          eb = edensphjp(0, 0, nexphjp)
          ee = edensphjp(0, 0, nexphjp + 1)
          ldbpos = phdensjp(0, 0, nexphjp, j, 1)
          ldbneg = phdensjp(0, 0, nexphjp, j, - 1)
          ldepos = phdensjp(0, 0, nexphjp + 1, j, 1)
          ldeneg = phdensjp(0, 0, nexphjp + 1, j, - 1)
        else
          eb = edensphjp(0, 0, nendens(0, 0) - 1)
          ee = edensphjp(0, 0, nendens(0, 0))
          ldbpos = phdensjp(0, 0, nendens(0, 0) - 1, j, 1)
          ldbneg = phdensjp(0, 0, nendens(0, 0) - 1, j, - 1)
          ldepos = phdensjp(0, 0, nendens(0, 0), j, 1)
          ldeneg = phdensjp(0, 0, nendens(0, 0), j, - 1)
        endif
        if (ldbpos > 1..and.ldepos > 1.) then
          lldbpos = log(ldbpos)
          lldepos = log(ldepos)
          ldtabpos = exp(lldbpos + (Eldpd(i) - eb) / (ee - eb) * (lldepos - lldbpos))
        else
          ldtabpos = ldbpos + (Eldpd(i) - eb) / (ee - eb) * (ldepos - ldbpos)
        endif
        if (ldbneg > 1..and.ldeneg > 1.) then
          lldbneg = log(ldbneg)
          lldeneg = log(ldeneg)
          ldtabneg = exp(lldbneg + (Eldpd(i) - eb) / (ee - eb) * (lldeneg - lldbneg))
        else
          ldtabneg = ldbneg + (Eldpd(i) - eb) / (ee - eb) * (ldeneg - ldbneg)
        endif
        phjpmdposj(i, j) = ldtabpos
        phjpmdnegj(i, j) = ldtabneg
      enddo
    enddo
    do i = 1, numdensracap
      phjpmdpos(i) = 0.0
      phjpmdneg(i) = 0.0
      do j = 0, numJph
        phjpmdpos(i) = phjpmdpos(i) + phjpmdposj(i, j)
        phjpmdneg(i) = phjpmdneg(i) + phjpmdnegj(i, j)
      enddo
    enddo
  endif
!
  if (ldmodelracap == 2) then
!
!   calculation of total ph-NLD with phmodel from 1 to 2
!
!   begin to prepare the input parameters for the function "phdens2" (phmodel=2, multi-pphh)
!
!   define the particle-hole numbers for proton and neutron in the function "phdens2"
!
    phnumracap = k0
    if (phnumracap == 1) then
      ptclp = 0
      holep = 0
      ptcln = 1
      holen = 1
    endif
    if (phnumracap == 2) then
      ptclp = 1
      holep = 1
      ptcln = 0
      holen = 0
    endif
    if (phnumracap == 3) then
      ptclp = 1
      holep = 1
      ptcln = 1
      holen = 1
    endif
    if (phnumracap == 4) then
      ptclp = 1
      holep = 1
      ptcln = 2
      holen = 2
    endif
    if (phnumracap == 5) then
      ptclp = 2
      holep = 2
      ptcln = 1
      holen = 1
    endif
    if (phnumracap == 6) then
      ptclp = 2
      holep = 2
      ptcln = 2
      holen = 2
    endif
!
!   define the other parameters in the function "phdens2"
!
    sigparleldenp = real(gp(0, 0))
    sigparleldenn = real(gn(0, 0))
    Ephracap = Etotal
    if (flaggshell) then
      racapdamp = ignatyuk(0, 0, Ephracap, 0) / alev(0, 0)
      sigparleldenp = sigparleldenp * racapdamp
      sigparleldenn = sigparleldenn * racapdamp
    endif
    surfwellE = 60.0
    surfwellgo = .false.
!
!   end of prepare the input parameters for the function "phdens2" (phmodel=2, multi-pphh)
!
    do i = 1, numdensracap
      phmdpos(i) = phdens2(0, 0, int(ptclp), int(holep), int(ptcln), int(holen), &
        sigparleldenp, sigparleldenn, Eldpd(i), surfwellE, surfwellgo) / 2.0
!test     phmdneg(i)=phmdpos(i)
    end do
!
!   calculation of spin-dependent ph-NLD with phmodel from 1 to 2
!
    do j = 0, numJph
      do i = 1, numdensracap
        phmdposj(i, j) = phmdpos(i) / (4.0 * dble(j + 1) - 2.0)
        phmdnegj(i, j) = phmdposj(i, j)
      end do
    end do
  end if
!
  if (ldmodelracap == 3) then
!
!   calculation of spin-dependent NLD with ldmodel=5
!
    xprty = - 1
    wprty = 1
    do j = 0, numJph
      do i = 1, numdensracap
        ldmdposj(i, j) = density(0, 0, Eldpd(i), Sldpd(j), int(wprty), int(0), 5)
        ldmdnegj(i, j) = density(0, 0, Eldpd(i), Sldpd(j), int(xprty), int(0), 5)
      enddo
    enddo
!
!   calculation of total NLD with ldmodel=5
!
    do i = 1, numdensracap
      ldmdpos(i) = 0.0
      ldmdneg(i) = 0.0
      do j = 0, numJph
        ldmdpos(i) = ldmdpos(i) + ldmdposj(i, j)
        ldmdneg(i) = ldmdneg(i) + ldmdnegj(i, j)
      enddo
    enddo
  endif
!
!   store values to chglpos[j] and chglneg[j]
!
!AK   do i=1,numdensracap
!AK     if (ldmodelracap.eq.2) then
!AK       chglpos(i)=dble(phmdpos(i))
!AK       chglneg(i)=dble(phmdneg(i))
!AK    elseif (ldmodelracap.eq.3) then
!AK       chglpos(i)=dble(ldmdpos(i))
!AK       chglneg(i)=dble(ldmdneg(i))
!AK     else
!AK       chglpos(i)=dble(phjpmdpos(i))
!AK       chglneg(i)=dble(phjpmdneg(i))
!AK     endif
!AK   enddo

  do j = 0, numJph
    do i = 1, numdensracap
      if (ldmodelracap == 2) then
        chglposj(i, j) = dble(phmdposj(i, j))
        chglnegj(i, j) = dble(phmdnegj(i, j))
      elseif (ldmodelracap == 3) then
        chglposj(i, j) = dble(ldmdposj(i, j))
        chglnegj(i, j) = dble(ldmdnegj(i, j))
      else
        chglposj(i, j) = dble(phjpmdposj(i, j))
        chglnegj(i, j) = dble(phjpmdnegj(i, j))
      endif
    enddo
  enddo
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   set the potential type (JLMB or Wood-Saxon) for racap.f:
!   JLMB and Wood-Saxon are available for neutron and proton;
!   Wood-Saxon is available for d,t,h, and a.
!
  if (k0 == 1 .or. k0 == 2) then
    if (flagjlm) then
       racopt = 3
     else
       racopt = 1
     endif
   else
     racopt = 1
   endif
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   calculation the final wood saxon potential parameter
!
  vncap2 = 0.
  rvncap2 = 0.
  avncap2 = 0.
  i0 = 0
  eopt = real(0.000001)
  if (racopt == 1) then
    call optical(i0, i0, k0, eopt)
    vncap2 = v
    rvncap2 = rv
    avncap2 = av
  endif
!
!   calculation the final state JLMB potential
!
  do i = 1, numjlm
    jlmracap2(i) = 0.
  enddo
  if (racopt == 3) then
    call mom(0, 0, dble(0.), dble(0.000001))
    do i = 1, numjlm
      jlmracap2(i) = normjlm(0, 0, 1) * (potjlm(0, 0, i, 1) + &
       normjlm(0, 0, 3) * (rhojlmn(0, 0, i, 1) + rhojlmp(0, 0, i, 1)) / &
       (rhojlmn(0, 0, i, 1) - rhojlmp(0, 0, i, 1)) * potjlm(0, 0, i, 3))
    enddo
  endif
!
!   set the levels scheme for direct capture
!
!   nlevexpracap=nlevmax2(0,0)
  nlevexpracap = Nlast(0, 0, 0) + 1
  ispect = 3
  if (ispect == 2) nlevexpracap = 1
!
!   ispect=1, only available experimental levels are used
!   ispect=2, only theoretical levels deduced from NLDs are used,
!   except the G.S.
!   ispect=3, firstly available experimental levels are used, and then
!   theoretical levels deduced from NLDs are used for the rest energy
!   range. (default)
!
!   set the initial value of spectroscopic factor firstly, then read the experimental spectroscopic factor, if any
!
!   The initial spectroscopic factor for experimental levels is 0.06, which is
!   the mean value of all available experimental data in the database(Nucl. Data Sheet).
!   Therefore, for all nuclei to be calculated, the minimum of global deviation between
!   the experimental compiled spectroscopic factor (read from database) and the estimated
!   one (the mean value 0.06 here) is guaranteed.
!
!   The initial spectroscopic factor for theoretical levels (deduced from NLDs) is 0.5.
!   This estimation minimizes the deviation between the cross section computed by
!   ispect=2 (only theoretical levels deduced from NLDs) and the one computed by
!   ispect=3 (available experimental levels and theoretical levels deduced from NLDs) for all nuclei.
!   Note that compared to ispect=3, only theoretical transition schemes deduced by the spin and
!   parity dependent NLDs are take into account for ispect=2, i.e., the allowed transitions to all
!   discrete experimental known final states are replaced in the energy range from zero to $E_{x}$.
!
  do nex = 0, numex
    if (nex <= nlevexpracap - 1) then
      spectfac(0, 0, nex) = spectfacexp(0, 0, nex)
    else
      spectfac(0, 0, nex) = spectfacth(0, 0)
    endif
  enddo
  if (ispect /= 2) then
    write(filespec, '("levels/spect", a1, "/", a, ".spect", a1)') parsym(k0), trim(nuc(ZZ(0, 0, 0))), parsym(k0)
  acp = AA(0, 0, 0)
  filespec = trim(path) //filespec
  inquire (file = filespec, exist = lexist)
  if ( .not. lexist) return
  open(unit = 10, file = filespec, status = 'old')
  do
    read(10, * , iostat = istat) i, ka, jlev
    if (istat == -1) exit
    if (ka == acp) then
Loop1:  do i = 1, jlev
        read(10, '(2x, i3, f8.4, f6.2, 2x, i2, f10.5)') j, e, spinf, parity, sf
        do nex = 0, Nlast(0, 0, 0)
          if (abs(edis(0, 0, nex) - e) < 1.d-1 .and. jdis(0, 0, nex) == spinf .and. parlev(0, 0, nex) == parity) then
            spectfac(0, 0, nex) = sf
            cycle Loop1
          endif
        enddo
      enddo Loop1
      exit
    else
      do i = 1, jlev + 1
        read(10, * , iostat = istat)
        if (istat == -1) exit
      enddo
    endif
  enddo
  close(10)
  endif
  return
end subroutine racapinit
! Copyright A.J. Koning 2021
