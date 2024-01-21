    subroutine sacs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Statistical analysis of cross sections
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
! All global variables
!   numia         ! maximum number of alphas in channel description
!   numid         ! maximum number of deuterons in channel description
!   numih         ! maximum number of helions in channel description
!   numin         ! maximum number of neutrons in channel description
!   numip         ! maximum number of protons in channel description
!   numit         ! maximum number of tritons in channel description
!   numP          ! number of points in reconstructed data file
! Variables for numerics
!   maxchannel    ! maximal number of outgoing particles in individual channel description
! Variables for main input
!   Atarget       ! mass number of target nucleus
!   k0            ! index of incident particle
!   Ztarget       ! charge number of target nucleus
! Variables for existence libraries
!   chanopen      ! flag to open channel with first non - zero cross s
! Variables for energies
!   Ethrexcl      ! threshold incident energy for exclusive channel
!   idchannel     ! identifier for exclusive channel
! Variables for exclusive channels
!   idnum         ! counter for exclusive channel
! Constants
!   parsym        ! symbol of particle
!
! *** Declaration of local data
!
  implicit none
  character(len=16) :: xsfile         ! file with channel cross sections
  character(len=132) :: line
  character(len=132) :: key
  integer           :: ia             ! mass number from abundance table
  integer           :: id             ! counter for deuterons
  integer           :: idc            ! help variable
  integer           :: ident          ! exclusive channel identifier
  integer           :: ih             ! hole number
  integer           :: in             ! counter for neutrons
  integer           :: ip             ! particle number
  integer           :: keyix
  integer           :: istat          ! logical for file access
  integer           :: it             ! counter for tritons
  integer           :: nen            ! energy counter
  integer           :: npart          ! number of particles in outgoing channel
  integer           :: Nxs            ! number of incident energies
  real(sgl)         :: Ea             ! start energy of local adjustment
  real(sgl)         :: Eb             ! end energy of local adjustment
  real(sgl)         :: Eh1            ! help variable
  real(sgl)         :: Eh2            ! help variable
  real(sgl)         :: Emax           ! maximal emission energy for particle channel
  real(sgl)         :: Ewidth         ! full width at half maximum
  real(sgl)         :: Exs(0:numP)    ! energy of cross section file
  real(sgl)         :: xs(0:numP)     ! help variable
  real(sgl)         :: xsa            ! help variable
  real(sgl)         :: xsb            ! help variable
  real(sgl)         :: xshalf         ! cross section at half maximum
  real(sgl)         :: xsmax          ! maximum cross sections
!
! ********* Read exclusive cross sections ******************************
!
  Exs = 0.
  xs = 0.
  open (unit = 1, file = 'sacs.dat', status = 'replace')
  write(1, '("# ", a1, " + ", a, ": Statistical analysis of cross sections")') parsym(k0), trim(targetnuclide)
  write(1, '("# Z   A     channel    E-thresh.   Emax    ", " xsmax    width ")')
  do npart = 0, maxchannel
  do ia = 0, numia
    do ih = 0, numih
      do it = 0, numit
        do id = 0, numid
          do ip = 0, numip
            do in = 0, numin
              if (in + ip + id + it + ih + ia /= npart) cycle
              if ( .not. chanopen(in, ip, id, it, ih, ia)) cycle
              ident = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
              do idc = 0, idnum
                if (idchannel(idc) == ident) then
                  xsfile = 'xs000000.tot'
                  write(xsfile(3:8), '(6i1)') in, ip, id, it, ih, ia
                  open (unit = 3, file = xsfile, status = 'old', iostat = istat)
                  if (istat == 0) then
                    Exs(0) = 0.
                    xs(0) = 0.
                    Emax = 0.
                    xsmax = 0.
                    do 
                      read(3,'(a)',iostat = istat) line
                      if (istat == -1) exit
                      key='entries'
                      keyix=index(line,trim(key))
                      if (keyix > 0) then
                        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nxs
                        read(3,'(/)')
                        do nen = 1, Nxs
                          read(3, *, iostat=istat ) Exs(nen), xs(nen)
                          if (istat == -1) exit
                        enddo
                        exit
                      endif
                    enddo
                    close (unit = 3)
                    do nen = 1, Nxs
                      if (xs(nen) > xsmax) then
                        Emax = Exs(nen)
                        xsmax = xs(nen)
                      endif
                    enddo
                    xshalf = 0.5 * xsmax
                    Eh1 = 0.
                    Eh2 = 0.
                    Ewidth = 0.
                    do nen = 1, Nxs
                      if (xs(nen - 1) < xshalf .and. xs(nen) >= xshalf) then
                        Ea = Exs(nen - 1)
                        Eb = Exs(nen)
                        xsa = xs(nen - 1)
                        xsb = xs(nen)
                        call pol1(xsa, xsb, Ea, Eb, xshalf, Eh1)
                      endif
                      if (xs(nen - 1) > xshalf .and. xs(nen) <= xshalf) then
                        Ea = Exs(nen - 1)
                        Eb = Exs(nen)
                        xsa = xs(nen - 1)
                        xsb = xs(nen)
                        call pol1(xsa, xsb, Ea, Eb, xshalf, Eh2)
                        exit
                      endif
                    enddo
                    if (Eh1 > 0..and.Eh2 > 0.) Ewidth = Eh2 - Eh1
                    write(1, '(2i4, 1x, a15, 2f8.3, es12.5, f10.3)') Ztarget, Atarget, xsfile, Ethrexcl(idc, 0), Emax, xsmax, Ewidth
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    enddo
  enddo
  close (unit = 1)
  return
end subroutine sacs
! Copyright A.J. Koning 2021
