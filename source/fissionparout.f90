subroutine fissionparout(Zix, Nix)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output for fission parameters
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
!   axtype           ! type of axiality of barrier
!   betafiscor       ! adjustable factor for fission path width
!   betafiscoradjust ! adjustable factor for fission path width
!   fbarrier         ! height of fission barrier
!   fismodelx        ! fission model
!   flagclass2       ! flag for class2 states in fission
!   fwidth           ! width of fission barrier
!   vfiscor          ! adjustable factor for fission path height
!   vfiscoradjust    ! adjustable factor for fission path height
!   rmiufiscor       ! adjustable factor for inertia mass along fission path 
!   widthc2          ! width of class2 states
! Variables for level density
!   Rclass2mom    ! norm. constant for moment of inertia for class 2 states
!   Rtransmom     ! norm. constant for moment of inertia for transition states
! Variables for nuclides
!   AA            ! mass number of residual nucleus
!   NN            ! neutron number of residual nucleus
!   ZZ            ! charge number of residual nucleus
! Constants
!   cparity       ! parity (character)
!   nuc           ! symbol of nucleus
! Variables for fission parameters
!   efisc2hb      ! energy of class2 states
!   efisc2rot     ! energy of class2 rotational transition states
!   efistrhb      ! energy of head band transition states
!   efistrrot     ! energy of rotational transition states
!   fecont        ! start of continuum energy
!   jfisc2hb      ! spin of class2 states
!   jfisc2rot     ! spin of class2 rotational transition states
!   jfistrhb      ! spin of head band transition states
!   jfistrrot     ! spin of rotational transition states
!   minertc2      ! moment of inertia for class2 states
!   minertia      ! moment of inertia of fission barrier deformation
!   nclass2       ! number of sets of class2 states
!   nfisbar       ! number of fission barrier parameters
!   nfisc2hb      ! number of class2 states for barrier
!   nfisc2rot     ! number of class2 rotational transition states for barrier
!   nfistrhb      ! number of head band transition states for barrier
!   nfistrrot     ! number of rotational transition states for barrier
!   pfisc2hb      ! parity of class2 states
!   pfisc2rot     ! parity of class2 rotational transition states
!   pfistrhb      ! parity of head band transition states
!   pfistrrot     ! parity of rotational transition states
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=15) :: col(4)    ! header
  character(len=15) :: un(4)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=80) :: fisfile
  character(len=200) :: topline   ! topline
  integer :: A   ! mass number of target nucleus
  integer :: i   ! counter
  integer :: j   ! counter
  integer :: Ncol
  integer :: indent
  integer :: id2
  integer :: id4
  integer :: N   ! neutron number of residual nucleus
  integer :: Nix ! neutron number index for residual nucleus
  integer :: Z   ! charge number of target nucleus
  integer :: Zix ! charge number index for residual nucleus
!
! ****************** Output of fission barrier parameters **************
!
! 1. Main fission parameters
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  Z = ZZ(Zix, Nix, 0)
  N = NN(Zix, Nix, 0)
  A = AA(Zix, Nix, 0)
  write(*, '(/" Fission information for Z=", i3, " N=", i3, " (", i3, a2, ") "/)')  Z, N, A, nuc(Z)
  massstring='   '
  write(massstring,'(i3)') A
  finalnuclide=trim(nuc(Z))//adjustl(massstring)
  fisfile = 'fis000000.txt'
  write(fisfile(4:9), '(2i3.3)') Z, A
  open (unit = 1, file = fisfile, status = 'replace')
  topline=trim(finalnuclide)//' fission parameters'
  call write_header(indent,topline,source,user,date,oformat)
  call write_residual(indent,Z,A,finalnuclide)
  call write_char(id2,'parameters','')
  call write_integer(id4,'Number of fission barriers',nfisbar(Zix, Nix))
  if (flagclass2) call write_integer(id4,'Number of sets of class2 states',nclass2(Zix, Nix))
  if (fismodelx(Zix, Nix) >= 5) then
    call write_real(id4,'Correction factor betafiscor',betafiscor(Zix, Nix))
    call write_real(id4,'Correction factor vfiscor',vfiscor(Zix, Nix))
    call write_real(id4,'Correction factor rmiufiscor',rmiufiscor(Zix, Nix))
    call write_real(id4,'Adjustable factor betafiscoradjust',betafiscoradjust(Zix, Nix))
    call write_real(id4,'Adjustable factor vfiscoradjust',vfiscoradjust(Zix, Nix))
    call write_real(id4,'Adjustable factor rmiufiscoradjust',rmiufiscoradjust(Zix, Nix))
  endif
  col = ''
  un = ''
  col(1) = 'Number'
  col(2) = 'Energy'
  un(2) = 'MeV'
  col(3) = 'Spin'
  col(4) = 'Parity'
  Ncol = 4
  do i = 1, nfisbar(Zix, Nix)
    quantity='Head band transition states'
    call write_quantity(id2,quantity)
    call write_integer(id4,'Fission barrier',i)
    call write_integer(id4,'Type of axiality',axtype(Zix, Nix, i))
    call write_real(id4,'Height of fission barrier',fbarrier(Zix, Nix, i))
    call write_real(id4,'Width of fission barrier',fwidth(Zix, Nix, i))
    call write_real(id4,'Rtransmom',Rtransmom(Zix, Nix, i))
    call write_real(id4,'Moment of inertia',minertia(Zix, Nix, i))
    call write_integer(id4,'Number of head band transition states',nfistrhb(Zix, Nix, i))
    call write_real(id4,'Start of continuum energy',fecont(Zix, Nix, i))
!
! 2. Head band transition states
!
    if (nfistrhb(Zix, Nix, i) > 0) then
      call write_datablock(id2,Ncol,nfistrhb(Zix, Nix, i),col,un)
      do j = 1, nfistrhb(Zix, Nix, i)
        write(1, '(6x, i6, 3x, es15.6, 6x,f6.1, 3x, 6x, a1)') j, efistrhb(Zix, Nix, i, j), jfistrhb(Zix, Nix, i, j), &
 &        cparity(pfistrhb(Zix, Nix, i, j))
      enddo
    endif
!
! 3. Rotational bands transition states
!
    if (nfistrrot(Zix, Nix, i) > 0) then
      quantity='Rotational bands'
      call write_quantity(id2,quantity)
      call write_datablock(id2,Ncol,nfistrrot(Zix, Nix, i),col,un)
      do j = 1, nfistrrot(Zix, Nix, i)
        write(1, '(6x, i6, 3x, es15.6, 6x,f6.1, 3x, 6x, a1)') j, efistrrot(Zix, Nix, i, j), jfistrrot(Zix, Nix, i, j), &
 &        cparity(pfistrrot(Zix, Nix, i, j))
      enddo
    endif
  enddo
!
! 4. Class2 states
!
  if (flagclass2) then
    do i = 1, nclass2(Zix, Nix)
      call write_integer(id4,'Set  of class2 states',i)
      call write_real(id4,'Rclass2mom',Rclass2mom(Zix, Nix, i))
      call write_real(id4,'Moment of inertia',minertc2(Zix, Nix, i))
      call write_integer(id4,'Number of class2 states',nfisc2hb(Zix, Nix, i))
      call write_real(id4,'Width of class2 states (MeV)',widthc2(Zix, Nix, i))
      if (nfisc2hb(Zix, Nix, i) > 0) then
        quantity='Class 2 states'
        call write_quantity(id2,quantity)
        do j = 1, nfisc2hb(Zix, Nix, i)
          write(1, '(6x, i6, 3x, es15.6, 6x,f6.1, 3x, 6x, a1)') j, efisc2hb(Zix, Nix, i, j), &
 &          jfisc2hb(Zix, Nix, i, j), cparity(pfisc2hb(Zix, Nix, i, j))
        enddo
      endif
!
! 5. Rotational bands
!
      if (nfisc2rot(Zix, Nix, i) > 0) then
        quantity='Rotational bands'
        call write_quantity(id2,quantity)
        call write_datablock(id2,Ncol,nfisc2rot(Zix, Nix, i),col,un)
        do j = 1, nfisc2rot(Zix, Nix, i)
          write(1, '(6x, i6, 3x, es15.6, 6x,f6.1, 3x, 6x, a1)') j, efisc2rot(Zix, Nix, i, j), &
 &          jfisc2rot(Zix, Nix, i, j), cparity(pfisc2rot(Zix, Nix, i, j))
        enddo
      endif
    enddo
  endif
  close(1)
  call write_outfile(fisfile,flagoutall)
  return
end subroutine fissionparout
! Copyright A.J. Koning 2021
