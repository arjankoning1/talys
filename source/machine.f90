subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Machine dependent statements
!
! Author    : Arjan Koning
!
! 2023-07-28: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Variables for files
!   nulldev     ! null device
!   path        ! directory containing files to be read
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist  ! logical to determine existence
  character(len=132):: codedir ! code directory
  character(len=132):: TALYSDIR! code directory runtime defined
  character(len=132):: OS      ! OS:windows
  integer           :: i       ! counter
  integer           :: year    ! year
  integer           :: month   ! month
  integer           :: day     ! day
  integer           :: values(8) ! date and time values
!
! ********************* Set directory for structure data ***************
!
! The code directory can be changed here via a script, or manually.
!
  codedir = '/Users/koning/talys/'
  i = len_trim(codedir)
  if (codedir(i:i) /= '/') codedir = codedir//'/'
!
! Another option is to set an environment variable TALYSDIR, e.g. put
! export TALYSDIR=/Users/koning/talys/     
! in your ~/.profile file.
! If TALYSDIR is not set, getenv will simply return an empty string
! (option provided by Viktor Zerkin: CODEDIR is in uppercase to
! prevent the path_change script changing this variable).
!
  call getenv('TALYSDIR',TALYSDIR)
  if (TALYSDIR.ne.' ') then
    CODEDIR=TALYSDIR
    i = len_trim(CODEDIR)
    if (CODEDIR(i:i) /= '/') CODEDIR(i+1:i+1)='/'
  endif
!
! Structure database
!
  path = trim(codedir)//'structure/'
  i = len_trim(path)
  if (path(i:i) /= '/') then
    i = i + 1
    path(i:i) = '/'
  endif
!
! The null device is a "black hole" for output that is produced, but not of interest to the user.
! Some ECIS output files fall in this category.
! To ensure compatibility with Unix, Linux, Windows and other systems a null device string is used,
! of which the default setting is given here.
! The input file may also be used to alter this setting, through the nulldev keyword
! Windows option provided by Viktor Zerkin.
!
  nulldev = '/dev/null'
  call getenv('OS',OS)
  if (OS.eq.'Windows_NT') nulldev='NUL'
!
! Test to check accessibility of structure files
!
  inquire (file = trim(path)//'abundance/H.abun', exist = lexist)
  if (.not. lexist) then
    write(*,*) 'codedir:[',trim(codedir),']'
    write(*,*) 'TALYSDIR:[',trim(TALYSDIR),']'
    write(*,*) 'expected file:',trim(path)//'abundance/H.abun'
    write(*, '(" TALYS-error: Structure database not installed: change path in machine.f90")')
    call exit(77)
  endif
!
! Set date
!
  call date_and_time(VALUES=values)
  year=values(1)
  month=values(2)
  day=values(3)
  date='xxxx-xx-xx'
  write(date(1:4),'(i4.4)') year
  write(date(6:7),'(i2.2)') month
  write(date(9:10),'(i2.2)') day
  user = 'Arjan Koning'
  return
end subroutine machine
! Copyright A.J. Koning 2023
