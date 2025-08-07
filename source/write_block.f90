  subroutine write_header(indent,title,source,user,date,comment)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for title
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) title
  character(len=*) source
  character(len=*) user
  character(len=*) date
  character(len=*) comment
  integer          indent
  integer          id2
!
! ************* Write output block *************************************
!
  call write_char(indent,'header','')
  id2 = indent + 2
  call write_char(id2,'title',trim(title))
  call write_char(id2,'source',trim(source))
  if (len_trim(user) > 0) call write_char(id2,'user',trim(user))
  call write_char(id2,'date',trim(date))
  if (len_trim(comment) > 0) call write_char(id2,'format',trim(comment))
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_target(indent)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for target
! +---------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! ****************** Declarations **************************************
!
  implicit none
  integer          Zixx
  integer          Nixx
  integer          ilev
  integer          P
  integer          indent
  integer          id2
  real             J
  real             Ee
  real             tau0
!
! ************* Write output block *************************************
!
  id2 = indent + 2
  call write_char(indent,'target','')
  call write_nuc(id2,Ztarget,Atarget,targetnuclide)
  if (Liso > 0) then
    Zixx=parZ(k0)
    Nixx=parN(k0)
    ilev=Lisomer(Zixx, Nixx, Liso)
    Ee=edis(Zixx,Nixx,ilev)
    J=jdis(Zixx,Nixx,ilev)
    P=parlev(Zixx,Nixx,ilev)
    tau0=tau(Zixx,Nixx,ilev)
    call write_level(id2,Liso,ilev,Ee,J,P,tau0)
  endif
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_reaction(indent,reaction,Q,E,MF,MT)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for reaction
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) reaction
  integer          indent
  integer          id2
  integer          MF
  integer          MT
  real*8           Q
  real*8           E
!
! ************* Write output block *************************************
!
  id2 = indent + 2
  call write_char(indent,'reaction','')
  call write_char(id2,'type',trim(reaction))
  if (Q /= 0.D0 .or. E /= 0.D0) then
    call write_double(id2,'Q-value [MeV]',Q)
    call write_double(id2,'E-threshold [MeV]',E)
  endif
  if( MF > 0 .and. MT > 0) call write_endfMFMT(id2,MF,MT)
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_endfMFMT(indent,MF,MT)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for ENDF MF and MT numbers
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  integer          indent
  integer          MF
  integer          MT
!
! ************* Write output block *************************************
!
  call write_integer(indent,'ENDF_MF',MF)
  call write_integer(indent,'ENDF_MT',MT)
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_residual(indent,Z,A,nuc)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for residual nuclide
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) nuc
  integer          Z
  integer          A
  integer          indent
  integer          id2
!
! ************* Write output block *************************************
!
  call write_char(indent,'residual','')
  id2 = indent + 2
  call write_nuc(id2,Z,A,nuc)
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_nuc(indent,Z,A,nuc)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for nuclide
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) nuc
  integer          indent
  integer          Z
  integer          A
!
! ************* Write output block *************************************
!
  if (Z > 0) call write_integer(indent,'Z',Z)
  call write_integer(indent,'A',A)
  if (Z > 0) call write_char(indent,'nuclide',trim(nuc))
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_level(indent,I,L,E,J,P,tau)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for level
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  integer          indent
  integer          id2
  integer          I
  integer          L
  integer          P
  real             J
  real             E
  real             tau
!
! ************* Write output block *************************************
!
  call write_char(indent,'level','')
  id2 = indent + 2
  if ( L >= 0 ) then
    call write_integer(id2,'number',L)
    call write_real(id2,'energy [MeV]',E)
  endif
  if (J >= 0.) then
    call write_real(id2,'spin',J)
    call write_integer(id2,'parity',P)
  endif
  if ( I >= 0 ) call write_integer(id2,'isomer',I)
  if ( tau > 0.) call write_real(id2,'half-life [sec]',tau)
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_quantity(indent,quantity)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for quantity
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) quantity
  integer          indent
  integer          id2
!
! ************* Write output block *************************************
!
  id2 = indent + 2
  call write_char(indent,'quantity','')
  call write_char(id2,'type',trim(quantity))
  return
  end
!Copyright (C)  2025 A.J. Koning
  subroutine write_datablock(indent,Nc,Ne,col,un)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : July 25, 2025
! | Task  : Write output block for data block
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  integer            Ne
  integer            Nc
  integer            L
  integer            i
  integer            indent
  integer            id2
  integer            k
  integer            ibeg
  integer            width
  character(len=*)   col(Nc)
  character(len=*)   un(Nc)
  character(len=15)  word
  character(len=3000) obs
  character(len=3000) units
!
! ************* Write output block *************************************
!
  call write_char(indent,'datablock','')
  id2 = indent + 2
  call write_integer(id2,'columns',Nc)
  call write_integer(id2,'entries',Ne)
  width=15
  obs="##"
  k=2
  do i=1,Nc
    word=''
    L=len_trim(col(i))
    ibeg=max(width/2-L/2,0)
    word(ibeg+1:ibeg+L)=col(i)(1:L)
    obs=obs(1:k)//word
    k=k+width
  enddo
  write(1, '(a)') trim(obs)
  units="##"
  k=2
  do i=1,Nc
    word=''
    L=len_trim(un(i))
    ibeg=max(width/2-L/2,0)
    word(ibeg+1:ibeg+L)=un(i)(1:L)
    word(ibeg:ibeg)='['
    word(ibeg+L+1:ibeg+L+1)=']'
    units=units(1:k)//word
    k=k+width
  enddo
  write(1, '(a)') trim(units)
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_real(indent,key,x)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for real number
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=80) fmt
  integer           indent
  real              x
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",es13.6)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",es13.6)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),x
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_double(indent,key,x)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for double precision number
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=80) fmt
  integer           indent
  real*8            x
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",es13.6)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",es13.6)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),x
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_integer(indent,key,k)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for integer
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=80) fmt
  integer           indent
  integer           k
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",i0)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",i0)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),k
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_char(indent,key,word)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for character
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=*)  word
  character(len=80) fmt
  integer           indent
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",a)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",a)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),trim(word)
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_outfile(tfile,flagout)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : August 8, 2023
! | Task  : Write specific output file to main output file
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  logical            :: flagout
  logical            :: lexist
  character(len=*)   :: tfile
  character(len=800) :: substring
  character(len=800) :: string
  integer            :: istat
!
! ************* Write output file **************************************
!
  inquire (file = tfile, exist = lexist)
  if (lexist) then
    if (flagout) then
      write(*, '()') 
      open (unit = 1, file = tfile, status = 'old')
      do              
        read(1, '(a)', iostat = istat) string
        if (istat /= 0) exit
        if (index(string,'# header') > 0) then
          do
            read(1,'(a)') substring 
            if (substring(1:4) /= '#   ' .and. substring(1:8) /= '# target') then
              write(*, '(1x, a)') trim(substring)
              exit
            endif
          enddo
        else
          write(*, '(1x, a)') trim(string)
        endif
      enddo           
      close (unit = 1)  
    else            
      write(*,'("file: ",a)') trim(tfile) 
    endif           
  endif           
  return
  end
!Copyright (C)  2023 A.J. Koning
