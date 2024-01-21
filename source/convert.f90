subroutine convert(i)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Convert input line from upper case to lowercase
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
! Variables for reading input lines
!   inline           ! input line
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numkey=20       ! number of keywords
  character(len=132) :: keyword(numkey) ! keyword
  character(len=132) :: str             ! input line
  integer            :: i               ! counter
  integer            :: k               ! counter
  integer            :: lkey            ! length of keyword
  integer            :: m               ! counter
!
! ************** Convert uppercase to lowercase characters *************
!
! For easy handling of all the input parameters, the whole input, both keywords and values, is converted to lowercase characters,
! with the exception of filenames or other character strings.
!
  data (keyword(m), m = 1, numkey) / 'abundance', 'bestpath', 'class2file', 'deformfile', 'e1file', 'hbtransfile', &
    'm1file', 'energy', 'integral', 'levelfile', 'nulldev', 'ompenergyfile', 'optmod', 'optmodfilen', 'optmodfilep', &
    'radialfile', 'rescuefile', 'strucpath', 'tjadjust', 'yieldfile' /
  str = inline(i)
  do k = 1, 132
    if (inline(i)(k:k) >= 'A' .and. inline(i)(k:k) <= 'Z') inline(i)(k:k) = achar(iachar(inline(i)(k:k)) + 32)
  enddo
  do k = 0, 110
    do m = numkey, 1, -1
      lkey = len_trim(keyword(m))
      if (inline(i)(k+1:k+lkey) == trim(keyword(m))) then
        inline(i)(k + lkey + 1:132) = str(k + lkey + 1:132)
        return
      endif
    enddo
  enddo
end subroutine convert
! Copyright A.J. Koning 2021
