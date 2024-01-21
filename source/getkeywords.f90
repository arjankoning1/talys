subroutine getkeywords(line, word)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Retrieve keywords and values from input line
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch         ! character
  character(len=1)   :: chprev     ! character
  character(len=132) :: word(40)   ! words on input line
  character(len=132) :: line       ! input line
  integer            :: i          ! counter
  integer            :: ibeg       ! index to mark begin of word
  integer            :: iend       ! index to mark end of word
  integer            :: nkey       ! word counter
!
! **************** Read keywords and values from input line ************
!
! From each input line we retrieve the keyword, the number of values, and the values themselves.
! These are all stored in strings (the word array).
!
  word = ' '
  chprev = ' '
  nkey = 0
  do i = 1, 132
    if (i > 1) chprev = line(i - 1:i - 1)
    ch = line(i:i)
    if (ch /= ' ' .and. chprev == ' ') then
      nkey = nkey + 1
      ibeg = i
    endif
    if (ch == ' ' .and. chprev /= ' ') then
      iend = i - 1
      word(nkey) = line(ibeg:iend)
    endif
  enddo
  return
end subroutine getkeywords
! Copyright A.J. Koning 2021
