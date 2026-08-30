subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Machine dependent statements
!
! Author    : Arjan Koning
!
! 2026-08-24: Original code
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
  logical            :: lexist    ! logical to determine existence
  character(len=1024):: code_dir  ! code directory
  character(len=1024):: talys_dir ! code directory runtime defined
  character(len=1024):: talys_user
  character(len=132) :: OS        ! OS:windows
  integer            :: i         ! counter
  integer            :: envstat
  integer            :: n
  integer            :: year      ! year
  integer            :: month     ! month
  integer            :: day       ! day
  integer            :: values(8) ! date and time values
!
! ********************* Set directory for structure data ***************
!
!
! The preferred option is to set an environment variable TALYS_DIR, e.g. put in your ~/.profile or ~/.zshrc file.
!
! export TALYS_DIR=/path/to/talys/     
!
! If TALYS_DIR is not set, get_environment_variable will simply return an empty string
!
  call get_environment_variable('TALYS_DIR', talys_dir, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    code_dir = trim(talys_dir)
  else
!
! If for some reason the above does not work, the code directory can be changed here manually.
!
    code_dir = '/path/to/talys/'
  endif
  i = len_trim(code_dir)
  if (i > 0) then
    if (code_dir(i:i) /= '/') code_dir = trim(code_dir)//'/'
  endif
!
! Structure database
!
  path = trim(code_dir)//'structure/'
!
! The null device is a "black hole" for output that is produced, but not of interest to the user.
! Some ECIS output files are written to it.
! To ensure compatibility with macOS, Linux, Windows and other systems a null device string is used,
! of which the default setting is given here.
! The input file may also be used to alter this setting, through the nulldev keyword
!
  nulldev = '/dev/null'
  OS = ''
  call get_environment_variable('OS',OS)
  if (OS.eq.'Windows_NT') nulldev='NUL'
!
! Test to check accessibility of structure files
!
  inquire (file = trim(path)//'abundance/H.abun', exist = lexist)
  if (.not. lexist) then
    write(*, '(a)') 'TALYS error: structure database not found.'
    write(*, '(2a)') 'Expected file: ', trim(path)//'abundance/H.abun'
    write(*, '(a)') 'Set the TALYS_DIR environment variable:'
    write(*, '(a)') '  export TALYS_DIR=/path/to/talys'
    write(*, '(a)') 'Alternatively, edit code_dir in source/machine.f90'
    write(*, '(a)') 'and rebuild TALYS.'
    error stop 77
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
!
! Set user name for output files
! The best option is to set an environment variable TALYS_USER, e.g. put in your ~/.profile or ~/.zshrc file.
!
! export TALYS_USER="Your Name" 
!
  call get_environment_variable('TALYS_USER', talys_user, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    user = trim(talys_user)
  else
    user = 'Unknown User'
  endif
  return
end subroutine machine
! Copyright A.J. Koning 2026
