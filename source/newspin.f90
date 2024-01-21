function newspin(A,e,n,iglob,f)
  implicit none
  ! Function to calculate spin cut-off parameter used
  ! for residual nucleus formed after pre-equilibrium emission.
  ! Parameters :
  ! ------------
  ! E : incident energy in MeV
  ! A : target nucleon number.
  ! n : exciton number.
  ! iglob : use 1 fo global parameterization, 0 for local parameterization (one for each nucleus).
  ! f : parameter to give in case of local parameterization (not yet used).


  real*8 :: e
  integer :: a,n,iglob ! iglob=1 : global, iglob=0 : local
  real*8 :: newspin,fa,fb,fc,a13,f(3)

  a13=real(a,8)**(1.d0/3.d0)
  !write(899,*)"a13,e",a13,e
  if(iglob == 1) then
     fa= -0.680858d0  * a13 + 6.37643d0
     fb=  0.0818127d0 * a13 - 0.557189d0
     fc=  1.14442d0   * a13 - 7.1285d0
!    write(899,*)fa,fb,fc
  else
     fa=f(1)
     fb=f(2)
     fc=f(3)
  endif
  !write(899,*)fa*e**.25d0, fb*e**.75d0,fc + .5d0
  newspin = 0.7d0 * ( fa*e**.25d0 + fb*e**.75d0 + fc + .5d0) &
       * real(n,8)**.44d0

end function newspin

function spin_wigner(j,s)
  ! Funtion used to define spin distribution R(J) as a function of the spin cut-off.
  ! It is normalized such as sum_J R(J) = 1 (thus it is 2J+1 times the R(J) defined in preeqinit.f).
  ! Parameters:
  ! -----------
  ! J : spin
  ! s : spin cut-off
  implicit none
  real*8, intent(in) :: j,s
  real*8 :: spin_wigner
  spin_wigner=(2.d0*j+1.d0)**2.d0/(s**3.d0*sqrt(4.d0*atan(1.d0)))*&
       exp(-(j+0.5d0)**2.d0/(s**2.d0))
end function spin_wigner
