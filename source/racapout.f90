subroutine racapout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output for racap into racap.out and racap.tot files
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
! All global variables
!   numex           ! maximum number of excitation energies
! Variables for energies
!   Ninclow       ! number of incident energies below Elow
! Variables for input energies
!   eninc           ! incident energy in MeV
!   nin             ! counter for incident energy
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Ztarget         ! charge number of target nucleus
! Variables for gamma rays
!   ldmodelracap    ! level density model for direct radiative capture
! Variables for energy grid
!   Einc            ! incident energy in MeV
! Variables for excitation energy grid
!   Ex              ! excitation energy
!   maxex           ! maximum excitation energy bin for residual nucleus
! Variables for energies
!   Qres            ! Q - value for residual nucleus
! Variables for incident channel
!   xspopnuc        ! population cross section per nucleus
! Variables for nuclides
!   Q               ! Q - value
!   ZZ              ! charge number of residual nucleus
! Constants
!   cparity         ! parity (character)
!   nuc             ! symbol of nucleus
!   parA            ! mass number of particle
!   parN            ! neutron number of particle
!   parsym          ! symbol of particle
!   parZ            ! charge number of particle
! Variables for levels
!   edis            ! energy of level
!   jdis            ! spin of level
!   parlev          ! parity of level
! Variables for direct capture initialization
!   nlevexpracap    ! number of experimental levels in the final nucleus
!   nlevracap       ! number of levels in the final nucleus
!   racopt          ! OMP for radiative capture
!   spectfac        ! spectroscopic factor
!   xsracap         ! direct radiative capture cross section
!   xsracapEM       ! direct - semidirect radiative capture cross section as
! Variables for direct capture
!   xsracape        ! direct radiative capture cross section
!   xsracapecont    ! direct radiative capture continuum cross section
!   xsracapedisc    ! direct radiative capture discrete cross section
!   xsracappop      ! population cross section for radiative capture
!   xsracappopex    ! population cross section for radiative capture
! Variables for masses
!   nucmass         ! mass of nucleus
!   S               ! separation energy
!   specmass        ! specific mass for residual nucleus
!
! *** Declaration of local data
!
  implicit none
  character(len=50) :: racapfile    ! file with direct capture cross sections
  character(len=1)  :: cpar         ! symbol of parity
  character(len=3)  :: massstring
  character(len=6)  :: finalnuclide
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(6)    ! header
  character(len=15) :: un(6)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer           :: Acpracap     ! compound nucleus A
  integer           :: i            ! counter
  integer           :: Ncol
  integer           :: istat        ! error code
  integer           :: JJJ          ! spin
  integer           :: nracapout    ! energy index
  integer           :: partar       ! parity of target
  integer           :: Zcpracap     ! compound nucleus Z
  real(sgl)         :: ecm          ! energy in C.M. frame
  real(sgl)         :: spintar      ! spin of target
!
!************************************************************************
!
  nracapout = nin
  Acpracap = Atarget + parA(k0)
  Zcpracap = ZZ(0, 0, 0)
  ecm = real(Einc * specmass(parZ(k0), parN(k0), k0))
  spintar = jdis(parZ(k0), parN(k0), 0)
  partar = parlev(parZ(k0), parN(k0), 0)
  cpar = '+'
  if (partar == -1) cpar = '-'
  if (nin == Ninclow + 1) then
    open (unit = 2, file = 'racap.tot', status = 'unknown')
    write(2, '("Direct capture reaction on target ", a, " (Z =", i3, ") with projectile ", a1, " (Qvalue=", f7.3, " MeV)"/)') &
 &    trim(targetnuclide), Ztarget, parsym(k0), Q(0)
    write(2, '("Characteristics of the Direct Radiative Capture", " calculation:")')
    if (ldmodelracap == 1) write(2, '("  ldmodelracap=1: spin-, parity-dependent pphh NLDs used for racap calculation")')
    if (ldmodelracap == 2) write(2, '("  ldmodelracap=2: spin-, parity-independent pphh NLDs used for racap calculation")')
    if (ldmodelracap == 3) write(2, '("  ldmodelracap=3: spin-, parity-dependent total NLDs used for racap calculation ", &
 &    "(ldmodel=5)")')
    if (racopt == 1) write(2, '("  racopt=1: Woods-Saxon potential", " KD used for racap calculation")')
    if (racopt == 3) write(2, '("  racopt=3: JLMB potential", " used for racap calculation")')
    write(2, '("  Total number of transitions from", f5.1, a1, " GS to ", i3, a2, " experimental levels:", i3)') spintar, &
 &    cpar, Acpracap, nuc(Zcpracap), nlevexpracap
    write(2, '("  Total number of transitions from", f5.1, a1, " GS to ", i3, a2, " levels:", i3)') spintar, cpar, Acpracap, &
 &    nuc(Zcpracap), maxex(0, 0) - nlevexpracap + 1
    write(2, '(/"  Spectroscopic factors ", /)')
    do i = 0, numex
      if (i == 0 .or. edis(0, 0, i) > 0.) then
        write(2, '(1p, g12.4, 0p, f5.1, 1x, a2, 1p, g12.4)') edis(0, 0, i), jdis(0, 0, i), cparity(parlev(0, 0, i)), &
 &        spectfac(0, 0, i)
      endif
    enddo
    write(2, * )
  else
    open (unit = 2, file = 'racap.tot', status = 'unknown')
    do
      read(2, * , iostat = istat)
      if (istat /= 0) exit
    enddo
  endif
  backspace 2
  write(2, * )
  write(2, '("==========  Direct capture at Elab=", es10.3, "  ==========")') Einc
  write(2, '("   Ecm =", es10.3, " MeV: ", /, "   Direct   radiative capture xs =", es12.5, " mb", &
 &  "  (Discrete=", es10.3, " - Continuum=", es10.3, ")", / , "   HF+Preeq radiative capture xs =", es12.5, " mb", / , &
 &  "   Total    radiative capture xs =", es12.5, " mb")') ecm, xsracape, xsracapedisc, xsracapecont, &
    xspopnuc(0, 0) - xsracape, xspopnuc(0, 0)
  write(2, * )
  if (mod(Acpracap, 2) == 0) then
    write(2, 1222) (real(JJJ), JJJ = 0, 10), (real(JJJ), JJJ = 0, 10)
 1222   format('Nlvl', 2x, ' Ex  ', 2x, ' J ', '   Pi ', 2x, ' Sp ', 2x, &
 &  'J/p=tot', 3x, 11('Jp=', f4.1, '+', 1x), 11('Jp=', f4.1, '-', 1x))
  else
    write(2, 1223) (real(JJJ + 0.5), JJJ = 0, 10), (real(JJJ + 0.5), JJJ = 0, 10)
 1223   format('Nlvl', 2x, ' Ex  ', 2x, ' J ', '   Pi ', 2x, ' Sp ', 2x, &
 &  'J/p=tot', 3x, 11('Jp=', f4.1, '+', 1x), 11('Jp=', f4.1, '-', 1x))
  endif
  do i = 0, maxex(0, 0)
       if (i <= nlevexpracap - 1) then
        write(2, 1225) i, edis(0, 0, i), jdis(0, 0, i), parlev(0, 0, i), spectfac(0, 0, i), xsracappopex(i)
      else
        write(2, 1226) i, Ex(0, 0, i), spectfac(0, 0, i), xsracappopex(i), (xsracappop(i, JJJ, 1), JJJ = 0, 10), &
 &        (xsracappop(i, JJJ, - 1), JJJ = 0, 10)
      endif
  enddo
 1225 format(1x, i3, 1x, f6.3, 1x, f5.1, 1x, i3, 1x, f6.2, 1x, 1p, e10.3, 32e9.2)
 1226 format(1x, i3, 1x, f6.3, 1x, ' 99.9', 1x, '  0', 1x, f6.2, 1x, 1p, e10.3, 32e9.2)
  close(2)
!
! write output racap.out with summary of reaction cross section
!
  racapfile = 'racap.out'
  un = 'mb'
  if (nin == Ninclow + 1) then
    open(unit = 1, file = racapfile, status = 'unknown')
    quantity='cross section'
    reaction='('//parsym(k0)//',g)'
    massstring='   '
    write(massstring,'(i3)') Acpracap
    finalnuclide=trim(nuc(Zcpracap))//trim(adjustl(massstring))
    topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' direct capture '//trim(quantity)
    col(1)='E'
    un(1)='MeV'
    col(2)='xs'
    col(3)='xs(E1)'
    col(4)='xs(E2)'
    col(5)='xs(M1)'
    col(6)='xs(tot)'
    Ncol=6
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,Qres(0, 0, 0),0.D0,0,0)
    call write_residual(Zcpracap,Acpracap,finalnuclide)
    call write_real(2,'E-max [MeV]',S(0, 0, k0))
    call write_integer(2,'number of levels',nlevracap(0, 0))
    call write_datablock(quantity,Ncol,Ninc-Ninclow,col,un)
!   write(1, '("# ", a1, " + ", a, ": Direct Capture to ", i3, a2)') &
!&    parsym(k0), trim(targetnuclide), Acpracap, nuc(Zcpracap)
!   write(1, '("# Q-value    =", es12.5, " mass=", f11.6, " Emax=", f11.6)') Qres(0, 0, 0), nucmass(0, 0), S(0, 0, k0)
!   write(1, '("# # transitions from ", f5.1, a1, " GS to ", i3, a2, " levels:", i3)') spintar, cpar, Acpracap, &
!&    nuc(Zcpracap), nlevracap(0, 0)
!   write(1, '("# # energies =", i6)') nracapout
!   write(1, '("#    E         xs                    ", " xs(E1)      xs(E2)      xs(M1)     xs(tot)")')
    close(1)
  endif
  open(unit = 1, file = racapfile, status = 'unknown', position = 'append')
  write(1, '(6es15.6)') eninc(nin), xsracap(nin), xsracapEM(nin, 1, 1), xsracapEM(nin, 1, 2), &
 &  xsracapEM(nin, 0, 1), xspopnuc(0, 0)
  close(1)
  return
end subroutine racapout
! Copyright A.J. Koning 2021
