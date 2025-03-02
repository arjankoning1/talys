subroutine fissionout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of fission cross sections
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
! Variables for fission
!   filefission    ! flag for fission cross sections on separate file
! Variables for existence libraries
!   fisexist       ! flag for existence of fission cross section
! Variables for numerics
!   maxN           ! maximal number of neutrons away from initial compound nucleus
!   maxZ           ! maximal number of protons away from initial compound nucleus
! Variables for input energies
!   eninc          ! incident energy in MeV
!   nin            ! counter for incident energy
!   Ninc         ! number of incident energies
! Variables for main input
!   Atarget        ! mass number of target nucleus
!   k0             ! index of incident particle
! Variables for energies
!   Ninclow      ! number of incident energies below Elow
! Variables for incident channel
!   maxA           ! maximal number of nucleons away from initial compound nucleus
! Variables for energy grid
!   Einc           ! incident energy in MeV
! Variables for multiple emission
!   xsfeed         ! cross section from compound to residual nucleus
! Variables for nuclides
!   AA             ! mass number of residual nucleus
!   ZZ             ! charge number of residual nucleus
! Constants
!   nuc            ! symbol of nucleus
!   parsym         ! symbol of particle
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=18) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=15) :: col(2)     ! header
  character(len=15) :: un(2)     ! header
  character(len=80) :: quantity   ! quantity
  character(len=12) :: rpfile      ! file with residual production cross sections
  integer           :: A           ! mass number of target nucleus
  integer           :: Acomp       ! mass number index for compound nucleus
  integer           :: Ncol        ! number of columns
  integer           :: Ncomp       ! neutron number index for compound nucleus
  integer           :: nen         ! energy counter
  integer           :: Z           ! charge number of target nucleus
  integer           :: Zcomp       ! proton number index for compound nucleus
  real(sgl)         :: x1(numpfns) ! help variable
  real(sgl)         :: x2(numpfns) ! help variable
  real(sgl)         :: x3(numpfns) ! help variable
!
! *********************** Fission cross sections ***********************
!
  write(*, '(/" 4b. Fission cross section per fissioning", " nuclide"/)')
  do Acomp = 0, maxA
    do Zcomp = 0, maxZ
      Ncomp = Acomp - Zcomp
      if (Ncomp < 0 .or. Ncomp > maxN) cycle
      if (xsfeed(Zcomp, Ncomp, - 1) /= 0.) then
        Z = ZZ(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        write(*, '(1x, 2i4, " (", i3, a2, ")", es12.5)') Z, A, A, nuc(Z), xsfeed(Zcomp, Ncomp, - 1)
      endif
    enddo
  enddo
!
! Write results to separate file
!
  if (filefission) then
    do Acomp = 0, maxA
      do Zcomp = 0, maxZ
        Ncomp = Acomp - Zcomp
        if (Ncomp < 0 .or. Ncomp > maxN) cycle
        if (xsfeed(Zcomp, Ncomp, - 1) == 0..and. .not. fisexist(Zcomp, Ncomp)) cycle
        rpfile = 'rp000000.fis'
        Z = ZZ(Zcomp, Ncomp, 0)
        A = AA(Zcomp, Ncomp, 0)
        write(rpfile(3:8), '(2i3.3)') Z, A
        massstring='   '
        write(massstring,'(i3)') A
        finalnuclide=trim(nuc(Z))//adjustl(massstring)
        if ( .not. fisexist(Zcomp, Ncomp)) then
          fisexist(Zcomp, Ncomp) = .true.
          open (unit = 1, file = rpfile, status = 'replace')
          quantity='cross section'
          reaction='('//parsym(k0)//',f)'
          topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)//' per fissioning nuclide'
          col(1)='E'
          un(1)='MeV'
          col(2)='xs'
          un(2)='mb'
          Ncol=2
          call write_header(topline,source,user,date,oformat)
          call write_target
          call write_reaction(reaction,0.D0,0.D0,0,0)
          call write_residual(Z,A,finalnuclide)
          call write_quantity(quantity)
          call write_datablock(Ncol,Ninc,col,un)
          do nen = 1, Ninclow
            write(1, '(2es15.6)') eninc(nen), 0.
          enddo
          do nen = Ninclow + 1, nin - 1
            write(1, '(2es15.6)') eninc(nen), 0.
          enddo
        else
          open (unit = 1, file = rpfile, status = 'old', position = 'append')
        endif
        write(1, '(2es15.6)') Einc, xsfeed(Zcomp, Ncomp, -1)
        close (unit = 1)
        call write_outfile(rpfile,flagoutall)
      enddo
    enddo
  endif
!
! Phenomenological PFNS (Iwamoto model)
!
  if (pfnsmodel == 1) then
    call iwamoto(Zinit, Ainit, real(S(0,0,1)), Einc, Tmadjust, Fsadjust, Epfns, NEpfns, x1, x2, x3, Eavpfns(1))
    do nen = 1, NEpfns
      pfns(1,nen) = x1(nen)
      maxpfns(1,nen) = x2(nen)
      pfnscm(1,nen) = x3(nen)
    enddo
    call pfnsout
  endif
  return
end subroutine fissionout
! Copyright A.J. Koning 2021
