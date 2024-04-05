subroutine endfread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read ECIS results for incident particle on ENDF-6 energy
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
!   sgl             ! single precision kind
!   dbl             ! double precision kind
! Variables for basic reaction
!   flagendfecis    ! flag for new ECIS calculation for ENDF - 6 files
! Variables for best files
!   flagrescue      ! flag for final rescue: normalization to data
! Variables for direct reactions
!   flagsys         ! flag for reaction cross section from systematics
! Variables for input energies
!   eninc           ! incident energy in MeV
!   Ninc            ! number of incident energies
! Variables for main input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
! Variables for energies
!   Ninclow       ! number of incident energies below Elow
! Variables for energy grid
!   coullimit       ! energy limit for charged particle OMP calculation
!   ecisstatus      ! status of ECIS file
!   Crescue         ! adjustment factor for this incident energy
!   Erescue         ! energy grid for adjustment factors
!   frescue         ! adjustment factor
!   Nrescue         ! number of energies for adjustment factors
! Variables for inverse channel data
!   threshnorm      ! normalization factor at threshold
! Variables for nuclides
!   AA              ! mass number of residual nucleus
!   ZZ              ! charge number of residual nucleus
! Constants
!   parA            ! mass number of particle
!   parsym          ! symbol of particle
!   parZ            ! charge number of particle
! Variables for thermal cross sections
!   fxselastot      ! total elastic cross section (neutrons only) for i
!   fxsnonel        ! non - elastic cross section for incident channel
!   fxstotinc       ! total cross section (neutrons only) for incident
! Variables for ENDF data
!   e6              ! energies of ENDF - 6 energy grid in MeV
!   nen6            ! total number of energies
!   xscompel6       ! compound elastic cross section
!   xsnonel6        ! non - elastic cross section
!   xselas6         ! total elastic cross section (neutrons only) for ENDF - 6 file
!   xselassh6       ! shape elastic cross section (neutrons only) for ENDF - 6 file
!   xsnon6          ! non - elastic cross section for ENDF - 6 file
!   xsopt6          ! optical model reaction cross section for ENDF - 6 file
!   xsreac6         ! reaction cross section for ENDF - 6 file
!   xstot6          ! total cross section (neutrons only) for ENDF - 6 file
!
! *** Declaration of local data
!
  implicit none
  character(len=18) :: reaction   ! reaction
  character(len=15) :: col(4)    ! header
  character(len=15) :: un(4)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: topline   ! topline
  character(len=72) :: line          ! input line
  integer           :: A             ! mass number of target nucleus
  integer           :: infileendf    ! file with cross sections for ENDF file
  integer           :: istat         ! logical for file access
  integer           :: mt            ! MT number
  integer           :: Ncol          ! number of columns
  integer           :: nen           ! energy counter
  integer           :: nen2          ! energy counter
  integer           :: nend          ! help variable
  integer           :: Nxs           ! number of cross sections
  integer           :: Z             ! charge number of target nucleus
  real(sgl)         :: e             ! energy
  real(sgl)         :: Ea            ! start energy of local adjustment
  real(sgl)         :: Eb            ! end energy of local adjustment
  real(sgl)         :: Efac          ! help variable
  real(sgl)         :: enuc          ! incident energy in MeV per nucleon
  real(sgl)         :: tripathi      ! function for semi-empirical reaction cross section of
  real(sgl)         :: xsa           ! help variable
  real(sgl)         :: xsb           ! help variable
  real(sgl)         :: xsc           ! interpolated cross section
  real(sgl)         :: xsd           ! help variable
  real(sgl)         :: xsdife        ! difference in elastic cross section
  real(sgl)         :: xsdift        ! difference in total cross section
  real(dbl)         :: xs            ! help variable
!
! ************ Read total, reaction and elastic cross section **********
!
  if (flagendfecis) then
    open (unit = 3, file = 'ecis.endfcs', status = 'unknown')
    infileendf = 3
    open (unit = 23, file = 'endf.cs', status = 'unknown')
    do
      read(3, '(a72)', iostat = istat) line
      if (istat == -1) exit
      write(23, '(a72)') line
    enddo
    rewind 3
  else
    infileendf = 23
  endif
  do nen = 1, nen6
    e = real(e6(nen))
    if (k0 > 1 .and. e < coullimit(k0)) cycle
    read(infileendf, '(57x, i3)') Nxs
    if (Nxs > 1) then
      read(infileendf, *) xs
      xstot6(nen) = max(real(xs), 0.)
    endif
    read(infileendf, * ) xs
    xsreac6(nen) = max(real(xs), 0.)
    xsopt6(nen) = xsreac6(nen)
    if (Nxs == 3) then
      read(infileendf, * ) xs
      xselassh6(nen) = max(real(xs), 0.)
    endif
  enddo
  close (unit = 3, status = ecisstatus)
  close (unit = 23, status = ecisstatus)
  open (unit = 10, file = 'ecis.endfin', status = 'unknown')
  close (unit = 10, status = ecisstatus)
!
! ********** Compound elastic contribution and normalization ***********
!
! locate     : subroutine to find value in ordered table
! pol1       : subroutine for polynomial interpolation of first order file
!
    if (k0 >= 1) then
      do nen = 1, nen6
        e = real(e6(nen))
        if (e <= eninc(1)) then
          xsc = xscompel6(1)
          xsd = xsnonel6(1)
        else
          call locate(eninc, 1, Ninc, e, nend)
          Ea = eninc(nend)
          Eb = eninc(nend + 1)
          xsa = xscompel6(nend)
          xsb = xscompel6(nend + 1)
          call pol1(Ea, Eb, xsa, xsb, e, xsc)
          xsa = xsnonel6(nend)
          xsb = xsnonel6(nend + 1)
          call pol1(Ea, Eb, xsa, xsb, e, xsd)
        endif
        xselas6(nen) = xselassh6(nen) + xsc
        xsnon6(nen) = xsd
        xstot6(nen) = xselas6(nen) + xsnon6(nen)
!
! ************************ Adjustment factors **************************
!
! Set incident energy dependent adjustment factors (purely for fitting purposes).
!
        if (flagrescue) then
Loop1:    do mt = 1, 3
            if (Nrescue(mt, - 1) == 0) cycle
            Crescue(mt, - 1) = 1.
            if (e <= Erescue(mt, - 1, 1)) then
              Crescue(mt, - 1) = frescue(mt, - 1, 1)
              cycle
            endif
            if (e >= Erescue(mt, - 1, Nrescue(mt, - 1))) then
              Crescue(mt, - 1) = frescue(mt, - 1, Nrescue(mt, - 1))
              cycle
            endif
            do nen2 = 1, Nrescue(mt, - 1) - 1
              if (e > Erescue(mt, - 1, nen2) .and. e <= Erescue(mt, - 1, nen2 + 1)) then
                Efac = (e - Erescue(mt, - 1, nen2)) / (Erescue(mt, - 1, nen2 + 1) - Erescue(mt, - 1, nen2))
                Crescue(mt, - 1) = frescue(mt, - 1, nen2) + &
                  Efac * (frescue(mt, - 1, nen2 + 1) - frescue(mt, - 1, nen2))
                if (frescue(mt, - 1, nen2 + 1) > 1.e10) Crescue(mt, - 1) = frescue(mt, - 1, nen2)
                if (frescue(mt, - 1, nen2) > 1.e10) Crescue(mt, - 1) = frescue(mt, - 1, nen2 + 1)
                cycle Loop1
              endif
            enddo
          enddo Loop1
!
! Put difference in the elastic (or total) cross section
!
          if (Crescue(1, - 1) /= 1..and.Crescue(1, - 1) /= 0.) xsdift = xstot6(nen) * (1. / Crescue(1, - 1) - 1.)
          if (Crescue(2, - 1) /= 1..and.Crescue(2, - 1) /= 0.) xsdife = xselas6(nen) * (1. / Crescue(2, - 1) - 1.)
          if (Crescue(2, - 1) /= 1..and.Crescue(2, - 1) /= 0.) then
            xselas6(nen) = xselas6(nen) + xsdife
            xstot6(nen) = xstot6(nen) + xsdife
          else
            xselas6(nen) = xselas6(nen) + xsdift
            xstot6(nen) = xstot6(nen) + xsdift
          endif
        endif
      enddo
    else
      do nen = 1, nen6
        xsnon6(nen) = xsreac6(nen)
      enddo
    endif
!
! ************ Normalization with semi-empirical results ***************
!
! tripathi  : function for semi-empirical reaction cross section of Tripathi et al.
!
  if (flagsys(k0)) then
    Z = ZZ(0, 0, k0)
    A = AA(0, 0, k0)
    do nen = 1, nen6
      if (xsopt6(nen) == 0.) cycle
      e = real(e6(nen))
      enuc = e / parA(k0)
      xs = tripathi(parZ(k0), parA(k0), Z, A, enuc)
      if (xs == 0.) xs = xsopt6(nen) * threshnorm(k0)
      xsnon6(nen) = xs
      if (k0 == 1) xselas6(nen) = xselas6(nen) + xsopt6(nen) - xs
    enddo
  endif
!
! **************** Write total cross sections to file ******************
!
  quantity='cross section'
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='Non-elastic'
  col(3)='Elastic'
  col(4)='Total'
  Ncol=4
  open (unit = 1, file = 'endf.tot', status = 'replace')
  reaction='('//parsym(k0)//',all)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_datablock(quantity,Ncol,nen6+Ninclow,col,un)
  do nen = 1, Ninclow
    write(1, '(4es15.6)') eninc(nen), fxsnonel(nen), fxselastot(nen), fxstotinc(nen)
  enddo
  do nen = 1, nen6
    write(1, '(4es15.6)') e6(nen), xsnon6(nen), xselas6(nen), xstot6(nen)
  enddo
  close (unit = 1)
!
! Total cross sections only
!
  col(2)='xs'
  Ncol=2
  open (unit = 1, file = 'endftot.tot', status = 'replace')
  reaction='('//parsym(k0)//',tot)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.D0,0.D0,3,1)
  call write_datablock(quantity,Ncol,nen6+Ninclow,col,un)
  do nen = 1, Ninclow
    write(1, '(2es15.6)') eninc(nen), fxstotinc(nen)
  enddo
  do nen = 1, nen6
    write(1, '(2es15.6)') e6(nen), xstot6(nen)
  enddo
  close (unit = 1)
!
! Elastic cross sections only
!
  open (unit = 1, file = 'endfel.tot', status = 'replace')
  reaction='('//parsym(k0)//',el)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.D0,0.D0,3,2)
  call write_datablock(quantity,Ncol,nen6+Ninclow,col,un)
  do nen = 1, Ninclow
    write(1, '(2es15.6)') eninc(nen), fxselastot(nen)
  enddo
  do nen = 1, nen6
    write(1, '(2es15.6)') e6(nen), xselas6(nen)
  enddo
  close (unit = 1)
!
! Nonelastic cross sections only
!
  open (unit = 1, file = 'endfnon.tot', status = 'replace')
  reaction='('//parsym(k0)//',non)'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(reaction,0.D0,0.D0,3,3)
  call write_datablock(quantity,Ncol,nen6+Ninclow,col,un)
  do nen = 1, Ninclow
    write(1, '(2es15.6)') eninc(nen), fxsnonel(nen)
  enddo
  do nen = 1, nen6
    write(1, '(2es15.6)') e6(nen), xsnon6(nen)
  enddo
  close (unit = 1)
  return
end subroutine endfread
! Copyright A.J. Koning 2021
