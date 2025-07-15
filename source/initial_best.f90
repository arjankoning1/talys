subroutine initial_best
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for best files
!
! Author    : Arjan Koning
!
! 2025-07-10: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl               ! single precision kind
! All global variables
!   numisom     ! number of isomers
!   numlines    ! number of input lines
! Variables for best files
!   flagbest         ! flag to use best set of adjusted parameters
!   flagbestend      ! flag to put best set of parameters at end of input file
! Variables for files
!   path        ! directory containing files to be read
! Constants
!   iso            ! counter for isotope
!   nuc            ! symbol of nucleus
! Variables for reading input lines
!   inline         ! input line
!   nlines         ! number of input lines
!   nlines0        ! number of input lines
! Error handling
!   range_index_error    ! Test if index is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=10)  :: afile     ! TALYS file with mass number
  character(len=15)  :: bestchar  ! help variable
  character(len=132) :: key       ! keyword
  logical            :: lexist    ! logical to determine existence
  integer            :: i         ! counter
  integer            :: istat     ! logical for file access
  integer            :: inum      ! counter
  integer            :: lenbest   ! length of best file
  integer            :: nbest     ! number of lines in best file
!
! **************** Set best files *******************************
!
  if (.not.flagbest) return
  bestpath = ' '
  if (flagnatural) then
    if (iso /= 1) then
      nbest = nlines - nlines0
      do i = nbest + 1, nlines
        inline(i - nbest) = inline(i)
      enddo
      nlines = nlines0
    endif
  endif
!
! If requested by input: retrieve best set of adjusted input parameters
!
! convert : subroutine to convert input line from upper case to lowercase
!
  if (bestfile == ' ') then
    afile = '000.talys'
    write(afile(1:3), '(i3.3)') Atarget
    bestchar = ptype0//'-'//trim(nuc(Ztarget))//afile
    if (bestpath(1:1) == ' ') bestpath = 'best/                                   '
    do i = 1, 40
      if (bestpath(i:i) == ' ') then
        if (bestpath(i-1:i-1) /= '/') then
          bestpath(i:i) = '/'
          lenbest = i
        else
          lenbest = i - 1
        endif
        exit
      endif
    enddo
    if (Starget(2:2) == ' ') then
      if (Ltarget == 0) then
        write(bestpath(lenbest+1:lenbest+5), '(a1, i3.3, "/")') Starget(1:1), Atarget
        write(bestpath(lenbest+6:lenbest+20), '(a15)') bestchar
      else
        write(bestpath(lenbest+1:lenbest+6), '(a1, i3.3, "m/")') Starget(1:1), Atarget
        write(bestpath(lenbest+7:lenbest+21), '(a15)') bestchar
      endif
    else
      if (Ltarget == 0) then
        write(bestpath(lenbest+1:lenbest+6), '(a2, i3.3, "/")') Starget(1:2), Atarget
        write(bestpath(lenbest+7:lenbest+21), '(a15)') bestchar
      else
        write(bestpath(lenbest+1:lenbest+7), '(a2, i3.3, "m/")') Starget(1:2), Atarget
        write(bestpath(lenbest+8:lenbest+22), '(a15)') bestchar
      endif
    endif
    bestfile = trim(path) // trim(bestpath)
  endif
  inquire (file = bestfile, exist = lexist)
  if ( .not. lexist) then
    write(*, '(" TALYS-warning: best file does not exist: ", a)') trim(bestfile)
    return
  endif
  open (unit = 3, file = bestfile, status = 'old')
  inum = 0
  do
    read(3, '(a132)', iostat = istat) key
    if (istat /= 0) exit
    inum = inum + 1
    i = numlines - inum
    inline(i) = key
    call convert(i)
  enddo
  close (unit = 3)
  if (inum > 0) then
    if (flagbestend) then
      do i = 1, inum
        inline(nlines + i) = inline(numlines - i)
      enddo
    else
      do i = nlines, 1, -1
        inline(i + inum) = inline(i)
      enddo
      do i = 1, inum
        inline(i) = inline(numlines - i)
      enddo
    endif
    nlines = nlines + inum
  endif
  return
end subroutine initial_best
! Copyright A.J. Koning 2021
