program talys
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main program
!
! Author    : Arjan Koning, Stephane Hilaire and Stephane Goriely
!
! 2023-12-30: Original code
! 2024-10-25: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
!   |-------------------------------------------------------|
!   |                 TALYS-2.04                            |
!   |                 Arjan Koning                          |
!   |                 Stephane Hilaire                      |
!   |                 Stephane Goriely                      |
!   |                                                       |
!   | Email: A.Koning@iaea.org                              |
!   |-------------------------------------------------------|
!
! MIT License
!
! Copyright (c) 2024 Arjan Koning
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Variables for main input
!   flagnatural    ! flag for calculation of natural element
!
! ************* Input, initialization and reaction models **************
!
! machine      : subroutine for machine dependent statements
! constants    : subroutine for constants and initialization
! talysinput   : subroutine for user input and defaults
! talysinitial : subroutine for initialization of nuclear structure
! talysreaction: subroutine with reaction models
! natural      : subroutine for calculation of natural element
!
  call machine
  call constants
  call talysinput
  call talysinitial
  call talysreaction
  if (flagnatural) call natural
end program talys
! Copyright A.J. Koning 2024
