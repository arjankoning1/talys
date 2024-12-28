subroutine kd03(k, pr)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read global optical model parameters
!
! Author    : Arjan Koning
!
! 2024-07-18: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
!
  use A0_talys_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
! Variables for optical model
!   v1_0, etc.      ! 
!
! *** Declaration of local data
!
  implicit none
  character(len=1)  :: pot           ! help variable
  character(len=1)  :: pr            ! help variable
  character(len=3)  :: setstring     ! help variable
  character(len=30) :: par           ! help variable
  character(len=132):: ompfile       ! optical model parameter file
  integer           :: istat         ! logical for file access
  integer           :: k             ! designator for particle
  integer           :: set           ! 
  real(sgl)         :: x             !
!
! ************** Read KD03 optical model parameters from database ***********
!
  par='kd03'
  set=0
  if (pr == 'f' .or. pr == 'y') par='pruitt_federal'
  if (pr == 'd') par='pruitt_democratic'
  if (pr /= 'n') set=pruittset
  write(setstring,'(i3)') set
  pot = ''
  ompfile = trim(path)//'optical/global/'//trim(par)//'/'//trim(adjustl(setstring))//'/parameters.json'
  open (unit = 2, file = ompfile, status = 'old', iostat = istat)
   if (istat /= 0) call read_error(ompfile, istat)
  read(2,'()')
  read(2,'()')
  call jsonline(v1_0)
  call jsonline(v1_asymm)
  call jsonline(v1_A)
  call jsonline(x)
  if (k == 1) v2_0 = x
  call jsonline(x)
  if (k == 1) v2_A = x
  call jsonline(x)
  if (k == 1) v3_0 = x
  call jsonline(x)
  if (k == 1) v3_A = x
  call jsonline(x)
  if (k == 2) v2_0 = x
  call jsonline(x)
  if (k == 2) v2_A = x
  call jsonline(x)
  if (k == 2) v3_0 = x
  call jsonline(x)
  if (k == 2) v3_A = x
  call jsonline(v4_0)
  call jsonline(rv_0)
  call jsonline(rv_A)
  call jsonline(av_0)
  call jsonline(av_A)
  read(2,'()')
  read(2,'()')
  call jsonline(rc_0)
  call jsonline(rc_A)
  call jsonline(rc_A2)
  read(2,'()')
  read(2,'()')
  call jsonline(vso1_0)
  call jsonline(vso1_A)
  call jsonline(vso2_0)
  call jsonline(rso_0)
  call jsonline(rso_A)
  call jsonline(aso_0)
  read(2,'()')
  read(2,'()')
  call jsonline(wso1_0)
  call jsonline(wso2_0)
  read(2,'()')
  read(2,'()')
  call jsonline(x)
  if (k == 1) w1_0 = x
  call jsonline(x)
  if (k == 1) w1_A = x
  call jsonline(x)
  if (k == 2) w1_0 = x
  call jsonline(x)
  if (k == 2) w1_A = x
  call jsonline(w2_0)
  call jsonline(w2_A)
  read(2,'()')
  read(2,'()')
  call jsonline(d1_0)
  call jsonline(d1_asymm)
  call jsonline(d2_0)
  call jsonline(d2_A)
  call jsonline(d2_a2)
  call jsonline(d2_a3)
  call jsonline(d3_0)
  call jsonline(rd_0)
  call jsonline(rd_A)
  call jsonline(x)
  if (k == 1) ad_0 = x
  call jsonline(x)
  if (k == 1) ad_A = x
  call jsonline(x)
  if (k == 2) ad_0 = x
  call jsonline(x)
  if (k == 2) ad_A = x
  close(2)
  return
end subroutine kd03
! Copyright A.J. Koning 2021
subroutine jsonline(var)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read data from JSON line
!
! Author    : Arjan Koning
!
! 2024-07-18: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  implicit none
  character(len=132):: string        ! string
  integer           :: commaix       ! 
  integer           :: colix         ! 
  real              :: var           !
  read(2,'(a)') string  
  colix=index(string,':')
  commaix=index(string,',')
  if (commaix == 0) commaix=132
  read(string(colix+1:commaix-1),*) var
  return
end subroutine jsonline
! Copyright A.J. Koning 2021
