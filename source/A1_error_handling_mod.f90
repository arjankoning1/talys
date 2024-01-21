module A1_error_handling_mod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Subroutines for error handling of file opening and reading and parameter ranges
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
contains

subroutine read_error(errfile, istat, eor, eof, error, ival, xval, cval)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Message for file reading error
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  character(len=*)           :: errfile ! file with read error
  character(len=*), optional :: eor     ! action for end of record
  character(len=*), optional :: eof     ! action for end of file
  character(len=*), optional :: error   ! action for error in value
  character(len=*), optional :: cval    ! string near read error
  integer                    :: istat   ! error code
  integer, optional          :: ival    ! integer near read error
  real, optional             :: xval    ! real near read error
!
! ***************************** Error message **************************
!
  write(*, '(" TALYS-error: Error in ",a)') trim(errfile)
  write(*, '("              IOSTAT = ",i6)') istat
  if (istat == 0) write(*, '("              Wrong value")')
  if (istat == -1) then
    write(*, '("              End of file")')
    if (present(eof)) then
      if (eof == 'continue') then
        write(*, '("              Continuing...")')
        return
      endif
    endif
  endif
  if (istat == -2) then
    write(*, '("              End of record")')
    if (present(eor)) then
      if (eor == 'continue') then
        write(*, '("              Continuing...")')
        return
      endif
    endif
  endif
  if (istat > 0) then
    write(*, '("              Error in value")')
    if (present(ival)) then
      write(*, '("              Error near ",i9)') ival
    endif
    if (present(xval)) then
      write(*, '("              Error near ",es12.5)') xval
    endif
    if (present(cval)) then
      write(*, '("              Error near ",a)') cval
    endif
    if (present(error)) then
      if (error == 'continue') then
        write(*, '("              Continuing...")')
        return
      endif
    endif
  endif
  stop
  return
end subroutine read_error
! Copyright A.J. Koning 2021

subroutine range_integer_error(varname, variable, vmin, vmax, default, unit, flagoutside, error, &
 & index1, name1, index2, name2, index3, name3, index4, name4, index5, name5)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Test if integer variable is out of range
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  logical, optional          :: flagoutside! flag for value outside range
  character(len=*)           :: varname    ! variable name
  character(len=*), optional :: unit       ! unit for range description
  character(len=*), optional :: error      ! action for error in value
  character(len=*), optional :: name1      ! index name
  character(len=*), optional :: name2      ! index name
  character(len=*), optional :: name3      ! index name
  character(len=*), optional :: name4      ! index name
  character(len=*), optional :: name5      ! index name
  character(len=132)         :: unitstring ! unit for range description
  integer                    :: variable   ! variable
  integer                    :: vmin       ! minimum value
  integer                    :: vmax       ! maximum value
  integer, optional          :: default    ! default value
  integer, optional          :: index1     ! index
  integer, optional          :: index2     ! index
  integer, optional          :: index3     ! index
  integer, optional          :: index4     ! index
  integer, optional          :: index5     ! index
!
! ***************************** Error message **************************
!
  if (present(flagoutside)) flagoutside = .false.
  if (present(default)) then
    if (variable == default) return
  endif
  if (present(unit)) then
    unitstring=unit
  else
    unitstring=' '
  endif
  if (variable < vmin .or. variable > vmax) then
    if (present(error)) then
      if (error == 'continue') then
        write(*, '(" TALYS-warning: Variable out of range")')
      else
        write(*, '(" TALYS-error: Variable out of range")')
      endif
    endif
    write(*, '("             ",a," = ",i9)') trim(varname),variable
    write(*, '("       Range: ",i9," <= ",a," <= ",i9," ",a)') vmin,trim(varname),vmax,trim(unitstring)
    if (present(index1) .and. present(name1)) write(*, '("             ",a," = ",i9)') trim(name1),index1
    if (present(index2) .and. present(name2)) write(*, '("             ",a," = ",i9)') trim(name2),index2
    if (present(index3) .and. present(name3)) write(*, '("             ",a," = ",i9)') trim(name3),index3
    if (present(index4) .and. present(name4)) write(*, '("             ",a," = ",i9)') trim(name4),index4
    if (present(index5) .and. present(name5)) write(*, '("             ",a," = ",i9)') trim(name5),index5
    if (present(error)) then
      flagoutside = .true.
      if (error == 'continue') then
        write(*, '("              Continuing...")')
        return
      endif
    endif
    stop
  endif
  return
end subroutine range_integer_error
! Copyright A.J. Koning 2021

subroutine range_index_error(varname, indexname, variable, vmin, vmax, flagoutside, error)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Test if index (Z, A, fission barrier, etc.) is out of range
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  logical, optional          :: flagoutside! flag for value outside range
  character(len=*)           :: varname    ! variable name
  character(len=*)           :: indexname  ! index name
  character(len=*), optional :: error      ! action for error in value
  integer                    :: variable   ! variable
  integer                    :: vmin       ! minimum value
  integer                    :: vmax       ! maximum value
!
! ***************************** Error message **************************
!
  if (present(flagoutside)) flagoutside = .false.
  if (variable < vmin .or. variable > vmax) then
    if (present(error)) then
      if (error == 'continue') then
        write(*, '(" TALYS-warning: Index out of range")')
      else
        write(*, '(" TALYS-error: Index out of range")')
      endif
    endif
    write(*, '("             ",a)') trim(varname)
    write(*, '("       Index: ",a,i9)') trim(indexname),variable
    write(*, '("       Range: ",i9," <= ",a," <= ",i9)') vmin,trim(indexname),vmax
    if (present(error)) then
      flagoutside = .true.
      if (error == 'continue') then
        write(*, '("              Continuing...")')
        return
      endif
    endif
    stop
  endif
  return
end subroutine range_index_error
! Copyright A.J. Koning 2021

subroutine range_real_error(varname, variable, vmin, vmax, default, unit, index1, name1, index2, name2, index3, name3, &
 & index4, name4, index5, name5)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Test if real variable is out of range
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  character(len=*)           :: varname    ! variable name
  character(len=*), optional :: unit       ! unit for range description
  character(len=*), optional :: name1      ! index name
  character(len=*), optional :: name2      ! index name
  character(len=*), optional :: name3      ! index name
  character(len=*), optional :: name4      ! index name
  character(len=*), optional :: name5      ! index name
  character(len=132)         :: unitstring ! unit for range description
  integer, optional          :: index1     ! index
  integer, optional          :: index2     ! index
  integer, optional          :: index3     ! index
  integer, optional          :: index4     ! index
  integer, optional          :: index5     ! index
  real, optional             :: default    ! default value
  real                       :: variable   ! variable
  real                       :: vmin       ! minimum value
  real                       :: vmax       ! maximum value
!
! ***************************** Error message **************************
!
  if (present(default)) then
    if (variable == default) return
  endif
  if (present(unit)) then
    unitstring=unit
  else
    unitstring=' '
  endif
  if (variable < vmin .or. variable > vmax) then
    write(*, '(" TALYS-error: Variable out of range")')
    write(*, '("             ",a," = ",es12.5)') trim(varname),variable
    write(*, '("       Range: ",es12.5," <= ",a," <= ",es12.5," ",a)') vmin,trim(varname),vmax,trim(unitstring)
    if (present(index1) .and. present(name1)) write(*, '("             ",a," = ",i9)') trim(name1),index1
    if (present(index2) .and. present(name2)) write(*, '("             ",a," = ",i9)') trim(name2),index2
    if (present(index3) .and. present(name3)) write(*, '("             ",a," = ",i9)') trim(name3),index3
    if (present(index4) .and. present(name4)) write(*, '("             ",a," = ",i9)') trim(name4),index4
    if (present(index5) .and. present(name5)) write(*, '("             ",a," = ",i9)') trim(name5),index5
    stop
  endif
  return
end subroutine range_real_error
! Copyright A.J. Koning 2021

end module A1_error_handling_mod
