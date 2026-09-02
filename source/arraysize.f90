subroutine arraysize
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output actual memory use of selected large arrays
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2026-09-02: Runtime memory accounting with size() and storage_size()
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
  use iso_fortran_env, only: int64
!
! *** Declaration of local data
!
  implicit none
  integer        :: Nix
  integer        :: Zix
  integer(int64) :: localbytes
  integer(int64) :: localelems
  integer(int64) :: qrpa_e_elems
  integer(int64) :: qrpa_f_elems
  integer(int64) :: qrpa_fJP_elems
  integer(int64) :: qrpabytes
  integer(int64) :: totalbytes
!
! *************************** Array size *******************************
!
! The reported memory is the storage occupied by the array elements.
! It does not include executable code, local scalar variables, allocator
! overhead, runtime-library memory, or other process memory.
!
  totalbytes = 0_int64
  qrpa_e_elems = 0_int64
  qrpa_f_elems = 0_int64
  qrpa_fJP_elems = 0_int64
!
  write(*, '(/"  Memory use of selected large arrays"/)')
  open (unit=1, file='arrays.out', status='unknown')
  write(1, '(" Array",t27,"Elements",t44,"MiB")')
  write(1, '(" -------------------------------------------------------")')
!
! Static module arrays
!
  call report_array('gamexist',       size(gamexist,kind=int64),       storage_size(gamexist))
  call report_array('chanisoexist',   size(chanisoexist,kind=int64),   storage_size(chanisoexist))
  call report_array('bassign',        size(bassign,kind=int64),        storage_size(bassign))
  call report_array('ldtable',        size(ldtable,kind=int64),        storage_size(ldtable))
  call report_array('ldtableT',       size(ldtableT,kind=int64),       storage_size(ldtableT))
  call report_array('ldtableN',       size(ldtableN,kind=int64),       storage_size(ldtableN))
  call report_array('xspopph2',       size(xspopph2,kind=int64),       storage_size(xspopph2))
  call report_array('preeqpop',       size(preeqpop,kind=int64),       storage_size(preeqpop))
  call report_array('xspop',          size(xspop,kind=int64),          storage_size(xspop))
  call report_array('wemission2',     size(wemission2,kind=int64),     storage_size(wemission2))
  call report_array('xsrp',           size(xsrp,kind=int64),           storage_size(xsrp))
  call report_array('Tjlnex',         size(Tjlnex,kind=int64),         storage_size(Tjlnex))
  call report_array('fxsgamdischan',  size(fxsgamdischan,kind=int64),  storage_size(fxsgamdischan))
  call report_array('areaejlab',      size(areaejlab,kind=int64),      storage_size(areaejlab))
  call report_array('areareclab',     size(areareclab,kind=int64),     storage_size(areareclab))
  call report_array('xsastroex',      size(xsastroex,kind=int64),      storage_size(xsastroex))
  call report_array('sfactor',        size(sfactor,kind=int64),        storage_size(sfactor))
  call report_array('phdensjp',       size(phdensjp,kind=int64),       storage_size(phdensjp))
!
! Allocatable module arrays
!
  if (allocated(Erescue)) then
    call report_array('Erescue', size(Erescue,kind=int64), storage_size(Erescue))
  else
    call report_array('Erescue', 0_int64, storage_size(0.0_sgl))
  endif
!
  if (allocated(frescue)) then
    call report_array('frescue', size(frescue,kind=int64), storage_size(frescue))
  else
    call report_array('frescue', 0_int64, storage_size(0.0_sgl))
  endif
!
  if (allocated(phtable1)) then
    call report_array('phtable1', size(phtable1,kind=int64), storage_size(phtable1))
  else
    call report_array('phtable1', 0_int64, storage_size(0.0_sgl))
  endif
!
  if (allocated(phtable2)) then
    call report_array('phtable2', size(phtable2,kind=int64), storage_size(phtable2))
  else
    call report_array('phtable2', 0_int64, storage_size(0.0_sgl))
  endif
!
  if (allocated(ddxrec)) then
    call report_array('ddxrec', size(ddxrec,kind=int64), storage_size(ddxrec))
  else
    call report_array('ddxrec', 0_int64, storage_size(0.0_sgl))
  endif
!
! feedexcl is assumed to be allocatable after the memory-reduction change.
!
  if (allocated(feedexcl)) then
    call report_array('feedexcl', size(feedexcl,kind=int64), storage_size(feedexcl))
  else
    call report_array('feedexcl', 0_int64, storage_size(0.0_sgl))
  endif
!
! QRPA data are allocatable components of qrpa_type and must be counted
! nucleus by nucleus. The descriptor array qrpa itself is not included.
!
  do Nix = 0, numN
    do Zix = 0, numZ
      if (allocated(qrpa(Zix,Nix)%e)) &
        qrpa_e_elems = qrpa_e_elems + size(qrpa(Zix,Nix)%e,kind=int64)
      if (allocated(qrpa(Zix,Nix)%f)) &
        qrpa_f_elems = qrpa_f_elems + size(qrpa(Zix,Nix)%f,kind=int64)
      if (allocated(qrpa(Zix,Nix)%fJP)) &
        qrpa_fJP_elems = qrpa_fJP_elems + size(qrpa(Zix,Nix)%fJP,kind=int64)
    enddo
  enddo
!
  call report_array('qrpa%e',   qrpa_e_elems,   storage_size(0.0_sgl))
  call report_array('qrpa%f',   qrpa_f_elems,   storage_size(0.0_sgl))
  call report_array('qrpa%fJP', qrpa_fJP_elems, storage_size(0.0_sgl))
!
  qrpabytes = (qrpa_e_elems + qrpa_f_elems + qrpa_fJP_elems) * &
    int(storage_size(0.0_sgl),int64) / 8_int64
  write(1, '(" QRPA total",t27,i14,t44,f12.3)') &
    qrpa_e_elems + qrpa_f_elems + qrpa_fJP_elems, real(qrpabytes,dbl) / 1048576.d0
!
! Total of the module-array entries listed above.
!
  write(1, '(" -------------------------------------------------------")')
  write(1, '(" TOTAL listed module arrays",t44,f12.3," MiB")') real(totalbytes,dbl) / 1048576.d0
!
! Local automatic arrays cannot be queried with size() from this routine.
! Report their declared size separately; do not include them in totalbytes.
!
  write(1, '(/," Local automatic arrays (estimated from declarations)")')
  write(1, '(" Array",t27,"Elements",t44,"MiB")')
  write(1, '(" -------------------------------------------------------")')
!
! channels.f90:
! real(sgl) :: specexcl(0:numchantot,0:numpar,0:numex+1,0:numen)
!
  localelems = int(numchantot+1,int64) * int(numpar+1,int64) * &
    int(numex+2,int64) * int(numen+1,int64)
  localbytes = localelems * int(storage_size(0.0_sgl),int64) / 8_int64
  write(1, '(" specexcl (local)",t27,i14,t44,f12.3)') &
    localelems, real(localbytes,dbl) / 1048576.d0
!
  write(1, '()')
  close(1)
  call write_outfile('arrays.out',flagoutall)
  return
!
contains
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Report element count and memory of one module array and add it to the total.
!-----------------------------------------------------------------------------------------------------------------------------------
!
  subroutine report_array(name, nelem, nbits)
    implicit none
    character(len=*), intent(in) :: name
    integer(int64),   intent(in) :: nelem
    integer,          intent(in) :: nbits
    integer(int64)               :: nbytes
!
    nbytes = nelem * int(nbits,int64) / 8_int64
    totalbytes = totalbytes + nbytes
    write(1, '(" ",a24,t27,i14,t44,f12.3)') trim(name), nelem, real(nbytes,dbl) / 1048576.d0
    return
  end subroutine report_array
!
end subroutine arraysize
! Copyright A.J. Koning 2026
