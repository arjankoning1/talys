      INTEGER*4 FUNCTION U_Valid(I_Z,I_A)
      IMPLICIT NONE
      INTEGER*4 I_Z
      INTEGER*4 I_A
      INTEGER*4  Ivalid
      Ivalid = 1
C     '   If I_A / I_Z < 210.E0/90.E0
      IF (  I_A / I_Z .LT. 172.E0 / 80.E0 .OR. I_A / I_Z .GT.
     *250.E0/90.E0  ) THEN
      Ivalid = 0
      End If
      IF (  I_Z .LT. 76 .OR. I_Z .GT. 120  ) THEN
      Ivalid = 0
      End If
      U_Valid = Ivalid
      END
C     '
C     '
      REAL*4 FUNCTION U_Delta_S0(I_Z,I_A)
      IMPLICIT NONE
      INTEGER*4 I_Z
      INTEGER*4 I_A
C     ' I_Z and I_A refer to the fissioning nucleus90 22
      REAL*4  Delta
      Delta = 0.3
      IF (  I_Z .EQ. 90 .AND. I_A .EQ. 228  )  Delta = 0.70
C     'N
      IF (  I_Z .EQ. 90 .AND. I_A .EQ. 230  )  Delta = 0.6
C     'N
      IF (  I_Z .EQ. 90 .AND. I_A .EQ. 233  )  Delta = 0.3
C     '
      IF (  I_Z .EQ. 91 .AND. I_A .EQ. 228  )  Delta = 0.65
C     '
      IF (  I_Z .EQ. 92 .AND. I_A .EQ. 233  )  Delta = 0.5
C     'N
      IF (  I_Z .EQ. 92 .AND. I_A .EQ. 234  )  Delta = 0.6
C     'N
      IF (  I_Z .EQ. 92 .AND. I_A .EQ. 235  )  Delta = 0.3
      IF (  I_Z .EQ. 92 .AND. I_A .EQ. 236  )  Delta = 0.3
C     'N
      IF (  I_Z .EQ. 92 .AND. I_A .EQ. 237  )  Delta = 0.3
      IF (  I_Z .EQ. 92 .AND. I_A .EQ. 238  )  Delta = 0.3
      IF (  I_Z .EQ. 92 .AND. I_A .EQ. 239  )  Delta = 0.1
C     '
      IF (  I_Z .EQ. 93 .AND. I_A .EQ. 238  )  Delta = -0.1
C     'N
C     '
      IF (  I_Z .EQ. 94 .AND. I_A .EQ. 240  )  Delta = -0.1
C     'N
      IF (  I_Z .EQ. 94 .AND. I_A .EQ. 241  )  Delta = -0.5
C     'N
      IF (  I_Z .EQ. 94 .AND. I_A .EQ. 242  )  Delta = -0.15
C     'N
      IF (  I_Z .EQ. 94 .AND. I_A .EQ. 243  )  Delta = -0.45
C     'N
      IF (  I_Z .EQ. 94  )  Delta = 0.25
C     '
      IF (  I_Z .EQ. 95 .AND. I_A .EQ. 242  )  Delta = -0.35
C     'N
C     '
      IF (  I_Z .EQ. 95 .AND. I_A .EQ. 243  )  Delta = -0.1
C     'N
C     '
      IF (  I_Z .EQ. 95 .AND. I_A .EQ. 244  )  Delta = -0.1
C     '
      IF (  I_Z .EQ. 96 .AND. I_A .EQ. 244  )  Delta = 0
C     'N
      IF (  I_Z .EQ. 96 .AND. I_A .EQ. 246  )  Delta = -0.2
C     'N
      U_Delta_S0 = Delta
      END
C     '
C     '
      REAL*4 FUNCTION Getyield(E_rel,E_ref,T_low,T_high)
      IMPLICIT NONE
      REAL*4 E_rel
      REAL*4 E_ref
      REAL*4 T_low
      REAL*4 T_high
C     /' Erel: Energy relative to the barrier '/
C     /' T_low: Effective temperature below barrier '/
C     /' T_high: Effective temperature above barrier '/
      REAL*4  Exp1
      REAL*4  Yield
CAK
      REAL*4 expo1,expo2
C     '
      Exp1 = E_rel/T_low - E_ref/0.4
C     ' energy far below barrier
C     ' Subtraction of E_ref/0.4 to prevent numerical problems.
CAK
       expo1=E_rel / T_high - E_ref/0.4
       expo2=-E_rel/ (T_high*T_low/(T_high-T_low) )
C     ' energy far below barrier
       If(Exp1.lt.-50..or.expo1.gt.80..or.expo2.gt.80.) Then
C     IF (  Exp1 .LT. -50  ) THEN
CAK end
      Yield = 0
      Else
CAK
      Yield = Exp(expo1) * 1.E0 /(1.E0 + exp(expo2) )
C     Yield = Exp(E_rel / T_high - E_ref/0.4) * 1.E0 / (1.E0 +
C    *exp(-E_rel/ (T_high*T_low/(T_high-T_low) ) ) )
CAK end
      End If
C     '   print  E_rel,T_high,E_ref,Yield
      Getyield = Max(Yield,0.0)
C     '
      END
C     '
C     '
      REAL*4 FUNCTION F1(Z_S_A)
      IMPLICIT NONE
      REAL*4 Z_S_A
C     /' Fit to the lower part of the data '/
      REAL*4  Result
      Result = exp(-9.05E0 + 4.58E0 * Log(Z_S_A/2.3E0))
      F1 = Result
      END
      REAL*4 FUNCTION F2(Z_S_A)
      IMPLICIT NONE
      REAL*4 Z_S_A
C     /' Fit to the upper part of the data '/
      REAL*4  Result
      Result = exp(12.08E0 - 3.27E0 * Log(Z_S_A/2.3E0))
      F2 = Result
      END
C     '
      REAL*4 FUNCTION Masscurv(Z,A,RL,kappa)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      REAL*4 RL
      REAL*4 kappa
C     /'  Fit to  Data of Fig. 7 of                                             '/
C     /'  "Shell effect in the symmetric-modal fission of pre-actinide nuclei"  '/
C     /'  S. I. Mulgin,K.-H. Schmidt,A. Grewe,S. V. Zhdanov                  '/
C     /'  Nucl. Phys. A 640 (1998) 375
C     /' (From fit of the width of the mass distributions.) '/                                         '/
      REAL*4  RI,Result1,Result2,Result
      REAL*4  Z_square_over_A
      REAL*4  ZsqrA
      REAL*4  c_rot
      DATA c_rot/600.0/
      REAL*4 F1
      REAL*4 F2
C     '
      Z_square_over_A = Z**2/A
      RI = (A - 2*Z)/A
      ZsqrA = Z_square_over_A * (1.E0 - kappa * RI**2) / (1.E0 - kappa
     ** ((226.E0 - 2.E0*91.E0)/226.E0)**2) + c_rot * RL**2 /
     *A**(7.0/3.0)
C     ' Hasse & Myers
C     '      + 0.0017 * RL^2
C     '
      Result1 = F1(ZsqrA)
      Result2 = F2(ZsqrA)
      Result = Min(Result1,Result2)
      Masscurv = Result
C     '
      END
C     '
      REAL*4 FUNCTION Masscurv1(Z,A,RL,kappa)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      REAL*4 RL
      REAL*4 kappa
C     /'  Fit to  Data of Fig. 7 of                                             '/
C     /'  "Shell effect in the symmetric-modal fission of pre-actinide nuclei"  '/
C     /'  S. I. Mulgin,K.-H. Schmidt,A. Grewe,S. V. Zhdanov                  '/
C     /'  Nucl. Phys. A 640 (1998) 375
C     /' (The left part assumed to be valid for the yields of the fission channels.) '/                                         '/
      REAL*4  RI,Result1,Result2,Result
C     '    Dim As Single A,A_central,Z
      REAL*4  Z_square_over_A
      REAL*4  ZsqrA
      REAL*4  c_rot
      DATA c_rot/600.0/
      REAL*4 F1
      REAL*4 F2
C     '
C     'A_central = -28.8156 + Z * 2.86587  ' Stability line for heavy nuclei
C     '
      Z_square_over_A = Z**2/A
      RI = (A - 2*Z)/A
      ZsqrA = Z_square_over_A * (1.E0 - kappa * RI**2) / (1.E0 - kappa
     ** ((226.E0 - 2.E0*91.E0)/226.E0)**2) + c_rot * RL**2 /
     *A**(7.0/3.0)
C     ' Hasse & Myers
C     '      + 0.0017 * RL^2
C     '
      IF (  ZsqrA .LT. 36.0  ) THEN
C     ' adjusted to Y(S2) in light nuclei (80<Z<92)
      ZsqrA = ZsqrA + 0.9 * (36.0 - ZsqrA)
      End If
C     '
      Result1 = F1(ZsqrA)
C     '  Result2 = F2(ZsqrA)
C     '  Result = Min(Result1,Result2)
      Masscurv1 = Result1
C     '
      END
C     '
C     '
      REAL*4 FUNCTION De_Saddle_Scission(Z_square_over_Athird,
     *ESHIFTSASCI)
      IMPLICIT NONE
      REAL*4 Z_square_over_Athird
      REAL*4 ESHIFTSASCI
C     /' Energy release between saddle and scission '/
C     /' M. Asghar,R. W. Hasse,J. Physique C 6 (1984) 455 '/
      REAL*4  Result
      Result = (31.E0 - 11.E0) / (1550.E0 - 1300.E0) *
     *(Z_square_over_Athird - 1300.E0 + ESHIFTSASCI) + 11.E0
C     ' This formula with ESHIFTSASCI = 0 is the parameterisation of the results
C     ' of Ashgar and Hasse,JPC 6 (1984) 455,see
C     ' F. Rejmund,A. V. Ignatyuk,A. R. Junghans,K.-H. Schmidt
C     ' Nucl. Phys. A 678 (2000) 215
      Result = max(Result,0.0)
      De_Saddle_Scission = Result
      END
C     '
C     '
      REAL*4 FUNCTION TEgidy(A,DU,Fred)
      IMPLICIT NONE
      REAL*4 A
      REAL*4 DU
      REAL*4 Fred
C     /' Temperature parameter of the constant-temperature formula for the
C     nuclear level density.
C     Input parameters: A = Mass number of nucleus
C     DU = Shell effect (corrected for pairing:P=0 for odd-A nuclei)
C     From "Correlations between the nuclear level density parameters"
C     Dorel Bucurescu,Till von Egidy
C     Phys. Rev. C 72 (2005) 067304    and
C     "Systematics of nuclear level density parameters"
C     Dorel Bucurescu,Till von Egidy
C     J. Phys. G: Nucl. Part. Phys. 31 (2005) S1675 and
C     "Systematics of nuclear level density parameters"
C     Till von Egidy,Dorel Bucurescu
C     Phys. Rev. C 72 (2005) 044311 '/
      REAL*4  Temp_smooth,Temp,T_Fac
C     ' Temp_smooth = 17.45E0 / (A^0.666667E0)
C     ' Temp = (17.45E0 - 0.51E0 * DU + 0.051 * DU^2) / (A^0.666667E0)
      Temp_smooth = 1.0 / (0.0570 * A**0.6666667)
      Temp = 1.0 / ( (0.0570 + 0.00193*DU) * A**0.6666667)
C     ' from  PRC 80 (2009) 054310
      T_Fac = Temp / Temp_smooth
      Temp = Temp * Fred
C     /' (For influence of deformation) '/
      TEgidy = Temp
      END
C     '
C     '
      REAL*4 FUNCTION TRusanov(E,A)
      IMPLICIT NONE
      REAL*4 E
      REAL*4 A
C     /' Fermi-gas level density,parameterisation of Rusanov et al. '/
      IF (  E >0  ) THEN
      TRusanov = SQRT(E / (0.094E0 * A) )
      Else
      TRusanov = 0.0
      End If
      END
C     '
      REAL*4 FUNCTION LyMass(Z,A,beta)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      REAL*4 beta
C     '
C     /' liquid-drop mass,Myers & Swiatecki,Lysekil,1967  '/
C     /' pure liquid drop,without pairing and shell effects '/
C     '
C     /' On input:    Z     nuclear charge of nucleus        '/
C     /'              N     number of neutrons in nucleus    '/
C     /'              beta  deformation of nucleus           '/
C     /' On output:   binding energy of nucleus              '/
C     '
      REAL*4  pi
      PARAMETER (pi=3.14159)
      REAL*4  N
      REAL*4  alpha
      REAL*4  XCOM,XVS,XE,EL
C     '
      N = A - Z
      alpha = SQRT(5.E0/(4.E0*pi)) * beta
      XCOM = 1.E0 - 1.7826E0 * ((A - 2.E0*Z)/A)**2
C     /' factor for asymmetry dependence of surface and volume term '/
      XVS = - XCOM * (15.4941E0*A                    -
     *17.9439E0*A**(2.E0/3.E0)*(1.E0+0.4E0*Alpha**2))
C     /' sum of volume and surface energy '/
      XE = Z**2 * (0.7053E0/A**(1.E0/3.E0)*(1.E0-0.2E0*Alpha**2)
     *           - 1.1529E0/A)
      EL = XVS + XE
C     /'   EL = EL + LyPair(Z,A); '/
      LyMass = EL
      END
C     '
C     '
      REAL*4 FUNCTION LyPair(Z,A)
      IMPLICIT NONE
      INTEGER*4 Z
      INTEGER*4 A
C     /' Calculates pairing energy '/
C     /' odd-odd nucleus:   Lypair = 0 '/
C     /' even-odd nucleus:  Lypair = -12/sqr(A) '/
C     /' even-even nucleus: Lypair = -2*12/sqr(A) '/
      REAL*4  E_PAIR
C     '
      E_PAIR = - 12.E0 / SQRT(REAL(A)) * ( MOD((Z+1) , 2) + MOD((A-Z+1)
     *, 2))
C     '
      Lypair = E_PAIR
      END
C     '
C     '
      REAL*4 FUNCTION TFPair(Z,A)
      IMPLICIT NONE
      INTEGER*4 Z
      INTEGER*4 A
C     /' Pairing energy from Thomas-Fermi model of Myers and Swiatecki '/
C     /' Shifted that TFPair is zero for odd-odd nuclei '/
      INTEGER*4  N
      REAL*4  E_Pair
      N = A - Z
      IF (   MOD(Z,2)  .EQ. 0 .AND.  MOD(N,2)  .EQ. 0  ) THEN
C     /' even-even '/
      E_Pair = - 4.8E0 / Z**0.333333E0 - 4.8E0 / N**0.333333E0 + 6.6E0
     */ A**0.666666E0
      END IF
      IF (   MOD(Z,2)  .EQ. 0 .AND.  MOD(N,2)  .EQ. 1  ) THEN
C     /' even Z,odd N '/
      E_Pair = - 4.8E0 / Z**0.333333E0 + 6.6E0 / A**0.666666E0
      END IF
      IF (   MOD(Z,2)  .EQ. 1 .AND.  MOD(N,2)  .EQ. 0  ) THEN
C     /' odd Z,even N '/
      E_Pair = - 4.8E0 / N**0.333333E0 + 6.6E0 / A**0.666666E0
      END IF
      IF (   MOD(Z,2)  .EQ. 1 .AND.  MOD(N,2)  .EQ. 1  ) THEN
C     /' odd N,odd N '/
      E_Pair = 0.0
      END IF
      TFPair = E_Pair
      END
C     '
C     '
      REAL*4 FUNCTION Pmass(Z,A,beta)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      REAL*4 beta
C     /' Liquid-drop model of Pearson,2001 '/
      REAL*4  N,EA,BE
      REAL*4  avol
      DATA avol/-15.65/
      REAL*4  asf
      DATA asf/17.63/
      REAL*4  r0
      DATA r0/1.233/
      REAL*4  asym
      DATA asym/27.72/
      REAL*4  ass
      DATA ass/-25.60/
      REAL*4  alpha
      REAL*4  pi
      PARAMETER (pi=3.14159)
C     '
      N = A - Z
      alpha = SQRT(5.E0/(4.E0*pi)) * beta
      EA = avol + asf * A**(-0.333333)*(1.E0+0.4E0*Alpha**2) + 0.6E0 *
     *1.44E0 * Z**2 / (A**1.333333 * r0 )*(1.E0-0.2E0*Alpha**2) + (asym
     *+ ass * A**(-0.333333)) * (N-Z)**2 / A**2
      BE = EA * A
      Pmass = BE
      END
C     '
C     '
      REAL*4 FUNCTION FEDEFOP(Z,A,beta)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      REAL*4 beta
C     /' According to liquid-drop model of Pearson 2001 '/
      REAL*4  asf
      DATA asf/17.63/
      REAL*4  r0
      DATA r0/1.233/
      REAL*4  N,Alpha
      REAL*4  pi
      PARAMETER (pi=3.14159)
C     '
      N = A - Z
      alpha = SQRT(5.E0/(4.E0*pi)) * beta
      FEDEFOP = asf * A**(0.666667)*(0.4E0*Alpha**2) - 0.6E0 * 1.44E0 *
     *Z**2 / (A**0.333333 * r0 )*(0.2E0*Alpha**2)
      END
C     '
C     '
      REAL*4 FUNCTION FEDEFOLys(Z,A,beta)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      REAL*4 beta
      REAL*4 LYMASS
      FEDEFOLys = Lymass(Z,A,beta) - Lymass(Z,A,0.0)
      END
C     '
C     '
      REAL*4 FUNCTION LDMass(Z,A,beta)
      IMPLICIT NONE
CAK
      include "gef.cmb"
      REAL*4 Z
      REAL*4 A
      REAL*4 beta
      REAL*4  N,BEtab
      REAL*4 LYMASS
      REAL*4 FEDEFOLYS
      REAL*4 BEldmTF
      REAL*4 BEexp
      N = A - Z
CAK
      BEtab = BEldm(NINT(N),NINT(Z)) + 2.0 * 12.0 / SQRT(REAL(A)) -
C     BEtab = BEldmTF(NINT(N),NINT(Z)) + 2.0 * 12.0 / SQRT(REAL(A)) -
CAKend
     *0.00001433*Z**2.39
C     ' The values in BEtab are the negative binding energies!
C     ' Pairing in Thomas Fermi masses is zero for Z,N even !
      IF (  BEtab .EQ. 0.0  ) THEN
      BEtab = Lymass(Z,A,0.0)
C     '         Print "Warning: Binding energy of Z=";Z;",A=";A;" not in mass table,";                         " replaced by LYMASS"
C     '         Print "I_Mode = ";I_Mode
      End If
      LDMASS = BEtab + FEDEFOLys(Z,A,beta)
      END
C     '
      REAL*4 FUNCTION AME2012(IZ,IA)
      IMPLICIT NONE
CAK
      include "gef.cmb"
      INTEGER*4 IZ
      INTEGER*4 IA
C     ' Masses from the 2003 mass evaluation,complemented by TF masses
C     ' and Lysekil masses.
      REAL*4  BEexpval
      REAL*4  Z,A,N
      INTEGER*4  INeu
      REAL*4 LYPAIR
      REAL*4 U_SHELL
      REAL*4 LDMASS
      REAL*4 BEexp
      INeu = IA - IZ
      A = REAL(IA)
      Z = REAL(IZ)
      N = A - Z
CAK
      BEexpval = beldm(INeu,IZ)
C     BEexpval = BEexp(INeu,IZ)
CAK end
      IF (  BEexpval .GT. -1.E10  ) THEN
      AME2012 = BEexpval
      Else
      AME2012 = Ldmass(Z,A,0.0) + U_SHELL(IZ,IA) + Lypair(IZ,IA)
      End If
      END
C     '
      REAL*4 FUNCTION U_SHELL(Z,A)
      IMPLICIT NONE
CAK
      include "gef.cmb"
      INTEGER*4 Z
      INTEGER*4 A
      INTEGER*4  N
      REAL*4  Res
      REAL*4 ShellMO
      N = A - Z
CAK
      Res = ushel(N,Z)
C     Res = ShellMO(N,Z)
CAK end
      IF (  Res .GT. 0.0  )  Res = 0.3 * Res
C     ' KHS (12. Feb. 2012)
C     '      ' The positive shell effects for deformed nuclei seem to be too positive
C     ' This gives too many high-energetic prompt neutrons.
      U_SHELL = Res
      END
C     '
      REAL*4 FUNCTION U_SHELL_exp(IZ,IA)
      IMPLICIT NONE
      INTEGER*4 IZ
      INTEGER*4 IA
      REAL*4  Res
      REAL*4  Z,A
      REAL*4 LDMASS
      REAL*4 LYPAIR
      REAL*4 AME2012
      Z = REAL(IZ)
      A = REAL(IA)
C     '   Res = 2.0 * ( AME2012(IZ,IA) - Lypair(IZ,IA) - LDMass(Z,A,0.0) )    '          - 0.25 * ( AME2012(IZ,IA-1) - Lypair(IZ,IA-1) - LDMass(Z,A-1.0,0.0) )    '          - 0.25 * ( AME2012(IZ,IA+1) - Lypair(IZ,IA+1) - LDMass(Z,A+1.0,0.0) )    '          - 0.25 * ( AME2012(IZ+1,IA+1) - Lypair(IZ+1,IA+1) - LDMass(Z+1.0,A+1.0,0.0) )    '          - 0.25 * ( AME2012(IZ-1,IA-1) - Lypair(IZ-1,IA-1) - LDMass(Z-1.0,A-1.0,0.0) )
      Res = 0.5 * ( AME2012(IZ,IA) - Lypair(IZ,IA) - LDMass(Z,A,0.0) )
     *+ 0.125 * ( AME2012(IZ,IA-1) - Lypair(IZ,IA-1) - LDMass(Z,A-1.0,
     *0.0) ) + 0.125 * ( AME2012(IZ,IA+1) - Lypair(IZ,IA+1) - LDMass(Z,
     *A+1.0,0.0) ) + 0.125 * ( AME2012(IZ+1,IA+1) - Lypair(IZ+1,IA+1) -
     *LDMass(Z+1.0,A+1.0,0.0) ) + 0.125 * ( AME2012(IZ-1,IA-1) -
     *Lypair(IZ-1,IA-1) - LDMass(Z-1.0,A-1.0,0.0) )
      U_SHELL_exp = Res
      END
C     '
      REAL*4 FUNCTION U_SHELL_EO_exp(IZ,IA)
      IMPLICIT NONE
      INTEGER*4 IZ
      INTEGER*4 IA
C     ' Returns experimental shell and even-odd staggering
      REAL*4  Res
      REAL*4  Z,A
      REAL*4 LDMASS
      REAL*4 LYPAIR
      REAL*4 AME2012
      Z = REAL(IZ)
      A = REAL(IA)
      Res = AME2012(IZ,IA) - LDMass(Z,A,0.0)
      U_SHELL_EO_exp = Res
      END
C     '
C     '
C     '
      REAL*4 FUNCTION U_MASS(Z,A)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
C     /' LD + congruence energy + shell (no pairing) '/
      REAL*4  BE
      REAL*4 U_SHELL
      REAL*4 LDMASS
      IF (  Z .LT. 0 .OR. A .LT. 0  ) THEN
C     '       Print "U_Mass: Z,A",Z,A
      End If
      BE = Ldmass(Z,A,0.0) + U_SHELL(NINT(Z),NINT(A))
C     '    BE = AME2012(Cint(Z),Cint(A)) - Lypair(Z,A)
C     '    BE = Lymass(Z,A,0.0) + U_Shell(CInt(Z),CInt(A))
C     '    BE = Lymass(Z,A,0.0)
      U_MASS = BE
      END
C     '
C     '
      REAL*4 FUNCTION ECOUL(Z1,A1,beta1,Z2,A2,beta2,d)
      IMPLICIT NONE
      REAL*4 Z1
      REAL*4 A1
      REAL*4 beta1
      REAL*4 Z2
      REAL*4 A2
      REAL*4 beta2
      REAL*4 d
C     '
C     /' Coulomb potential between two nuclei                    '/
C     /' surfaces are in a distance of d                         '/
C     /' in a tip to tip configuration                           '/
C     '
C     /' approximate formulation                                 '/
C     /' On input: Z1      nuclear charge of first nucleus       '/
C     /'           A1      mass number of irst nucleus   '/
C     /'           beta1   deformation of first nucleus          '/
C     /'           Z2      nuclear charge of second nucleus      '/
C     /'           A2      mass number of second nucleus  '/
C     /'           beta2   deformation of second nucleus         '/
C     /'           d       distance of surfaces of the nuclei    '/
C     '
      REAL*4  N1,N2,recoul
      REAL*4  dtot
      REAL*4  r0
      DATA r0/1.16/
C     '
      N1 = A1 - Z1
      N2 = A2 - Z2
      dtot = r0 *( (Z1+N1)**0.3333333E0 * (1.E0+0.6666667E0*beta1)
     *        + (Z2+N2)**0.3333333E0 * (1.E0+0.6666667E0*beta2) ) + d
      REcoul = Z1 * Z2 * 1.44E0 / dtot
C     '
      ECOUL = REcoul
      END
C     '
C     '
      REAL*4 FUNCTION beta_light(Z,betaL0,betaL1)
      IMPLICIT NONE
      INTEGER*4 Z
      REAL*4 betaL0
      REAL*4 betaL1
C     /' Deformation of light fission fragment for S1 and S2 '/
C     /' Systematic correlation Z vs. beta for deformed shells '/
C     /' Z of fission fragment '/
      REAL*4  beta
      beta = (Z - betaL0) * betaL1/20.E0
      beta_light = beta
      END
C     '
C     '
      REAL*4 FUNCTION beta_heavy(Z,betaH0,betaH1)
      IMPLICIT NONE
      INTEGER*4 Z
      REAL*4 betaH0
      REAL*4 betaH1
C     /' Deformation of heavy fission fragment for S2 '/
C     /' Systematic correlation Z vs. beta for deformed shells '/
C     /' Z of fission fragment '/
      REAL*4  beta
      beta = (Z - betaH0) * betaH1/20.E0
      beta_heavy = beta
      END
C     '
C     '
C     '
      REAL*4 FUNCTION Z_equi(ZCN,A1,A2,beta1,beta2,d,Imode,POLARadd,
     *POLARfac)
      IMPLICIT NONE
      INTEGER*4 ZCN
      INTEGER*4 A1
      INTEGER*4 A2
      REAL*4 beta1
      REAL*4 beta2
      REAL*4 d
      INTEGER*4 Imode
      REAL*4 POLARadd
      REAL*4 POLARfac
C     /' Determines the minimum potential of the scission-point configuration
C     represented by two deformed nuclei divided by a tip distance d.
C     A1,A2,beta1,beta2,d are fixed,Z1 is searched for and returned on output.  '/
C     '
C     /' ZCN: Z of fissioning nucleus '/
C     /' A1: A of first fission fragment '/
C     /' A2: A of second fission fragment '/
C     /' beta1: deformation of first fission fragment '/
C     /' beta2: deformation of second fission fragment '/
C     /' d: tip distance '/
C     '
      REAL*4  RZ_equi
      REAL*4  RA1,RA2,RZCN,RACN
      REAL*4  Z1UCD,Z2UCD
      REAL*4  re1,re2,re3,eps1,eps2,DZ_Pol
C     /' help variables '/
      REAL*4 ECOUL
      REAL*4 LYMASS
C     '
      RA1 = REAL(A1)
      RA2 = REAL(A2)
      RZCN = REAL(ZCN)
      RACN = RA1 + RA2
      Z1UCD = RA1 / (RA1 + RA2) * RZCN
      Z2UCD = RZCN - Z1UCD
      re1 = LyMass( Z1UCD-1.E0,RA1,beta1 ) + LyMass( Z2UCD+1.E0,RA2,
     *beta2 ) + ECoul( Z1UCD-1.E0,RA1,beta1,Z2UCD+1.E0,RA2,beta2,d )
      re2 = LyMass( Z1UCD,RA1,beta1) + LyMass( Z2UCD,RA2,beta2) +
     *ECoul( Z1UCD,RA1,beta1,Z2UCD,RA2,beta2,d )
      re3 = LyMass( Z1UCD+1.E0,RA1,beta1 ) + LyMass( Z2UCD-1.E0,RA2,
     *beta2 ) + ECoul( Z1UCD+1.E0,RA1,beta1,Z2UCD-1.E0,RA2,beta2,d )
      eps2 = ( re1 - 2.E0*re2 + re3 ) / 2.E0
      eps1 = ( re3 - re1 ) / 2.E0
      DZ_Pol = -eps1 / ( 2.E0 * eps2 )
C     '
      IF (  DZ_Pol .GT. 2 .OR. DZ_Pol .LT. -2  )  DZ_Pol = 0
C     '
      IF (  Imode .GT. 0  ) THEN
C     /' Purely empirical enhancement of charge polarization '/
      DZ_POL = DZ_POL * POLARfac + POLARadd
      End If
C     '
      RZ_equi = Z1UCD + DZ_POL
      Z_equi = RZ_equi
      END
C     '
C     '
      SUBROUTINE Beta_opt_light(A1,A2,Z1,Z2,d,beta2_imposed,beta1_opt)
      IMPLICIT NONE
      REAL*4 A1
      REAL*4 A2
      REAL*4 Z1
      REAL*4 Z2
      REAL*4 d
      REAL*4 beta2_imposed
      REAL*4 beta1_opt
C     /' Determines the optimum deformation of the light fragment when the deformation of the
C     heavy fragment is imposed. '/
C     '
      REAL*4  beta1,dbeta1,beta1_prev,beta1_next
      REAL*4  Uguess,Uplus,Uminus,Uprev,Unext
      INTEGER*4  I
      REAL*4 ECOUL
      REAL*4 LYMASS
C     '
C     /' List('Beta_opt_light called with ');
C     List(A1,A2,Z1,Z2,d,beta2_imposed,beta1_opt);
C     DCL Byes Bit(1) aligned;
C     Call GPYES('Continue',Byes); '/
      beta1 = 0.5
      dbeta1 = 0.01
      Uguess = LyMass(Z1,A1,beta1) + Lymass(Z2,A2,beta2_imposed) +
     *ECoul(Z1,A1,beta1,Z2,A2,beta2_imposed,d)
      Uplus = LyMass(Z1,A1,beta1 + dbeta1) + Lymass(Z2,A2,
     *beta2_imposed) + ECoul(Z1,A1,beta1 + dbeta1,Z2,A2,beta2_imposed,
     *d)
      Uminus = LyMass(Z1,A1,beta1 - dbeta1) + Lymass(Z2,A2,
     *beta2_imposed) + ECoul(Z1,A1,beta1 - dbeta1,Z2,A2,beta2_imposed,
     *d)
      IF (  Uplus .GT. Uguess .AND. Uminus .GT. Uguess  ) THEN
      beta1_opt = beta1
      Else
      IF (  Uplus .LT. Uguess  )  dbeta1 = 0.01
      IF (  Uminus .LT. Uguess  )  dbeta1 = -0.01
      Unext = Uguess
      beta1_next = beta1
      DO I = 1 , 10000
      beta1_prev = beta1_next
      Uprev = Unext
      beta1_next = beta1_prev + dbeta1
      Unext = LyMass(Z1,A1,beta1_next) + Lymass(Z2,A2,beta2_imposed) +
     *ECoul(Z1,A1,beta1_next,Z2,A2,beta2_imposed,d)
      IF (  Unext .GE. Uprev  )  Exit
      END DO
      beta1_opt = beta1_prev
      END IF
C     '
      END
C     '
C     '
      SUBROUTINE Beta_Equi(A1,A2,Z1,Z2,d,beta1prev,beta2prev,beta1opt,
     *beta2opt)
      IMPLICIT NONE
      REAL*4 A1
      REAL*4 A2
      REAL*4 Z1
      REAL*4 Z2
      REAL*4 d
      REAL*4 beta1prev
      REAL*4 beta2prev
      REAL*4 beta1opt
      REAL*4 beta2opt
C     /' Determines the minimum potential of the scission-point configuration
C     represented by two deformed nuclei,divided by a tip distance d.
C     A1,A2,Z1,Z2,d are fixed,beta1 and beta2 are searched for and returned on output '/
C     '
      INTEGER*4  B_analytical
      DATA B_analytical/0/
C     ' Switch to use the analytical approximation
C     ' that replaces the long numerical calculation.
      REAL*4  x,y,xcoul
      REAL*4  xcoul236U
      DATA xcoul236U/1369.64/
C     '
      REAL*4  beta1,beta2
C     '
C     '      Dim As Double U,Uprev,Ulast,Ubest,Uopt
      REAL*4  U,Uprev,Ulast,Ubest,Uopt
C     '
C     '      Dim As Double sbeta1,sbeta2
      REAL*4  sbeta1,sbeta2
C     '
      INTEGER*4  N,N1,N2,Nopt
C     '
C     '      Dim As Double eps = 5.E-4
      REAL*4  eps
      DATA eps/5.E-4/
C     '
      INTEGER*4  I
      REAL*4 LYMASS
      REAL*4 ECOUL
C     '
      IF (  B_analytical .EQ. 0  ) THEN
C     ' Numerical algorithm
C     '
      beta1 = beta1prev
      beta2 = beta2prev
CAK
      sbeta1=0.
      sbeta2=0.
      Nopt=0
CAK end
      Uprev = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECoul(Z1,A1,
     *beta1,Z2,A2,beta2,d)
      Uopt = Uprev
C     '
C     /' Test slope of variation of U '/
      beta1 = beta1prev + eps
      U = 1.E30
C     '
      beta2 = beta2prev
C     '     For beta2 = beta2prev to 0 Step -eps
      DO I = 1 , Int(beta2prev/eps)
      beta2 = beta2 - eps
      Ulast = U
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECoul(Z1,A1,beta1,
     *Z2,A2,beta2,d)
      IF (  U .GT. Ulast  ) THEN
      Exit
      Else
      Ubest = U
      END IF
      END DO
      IF (  Ubest .LT. Uopt  ) THEN
      Uopt = Ubest
      sbeta1 = eps
      sbeta2 = -eps
      END IF
C     '
      U = 1.E30
      beta2 = beta2prev
C     '   For beta2 = beta2prev To 1 Step eps
      DO I = 1 , Int((1 - beta2prev)/eps)
      beta2 = beta2 + eps
      Ulast = U
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECoul(Z1,A1,beta1,
     *Z2,A2,beta2,d)
      IF (  U .GT. Ulast  ) THEN
      Exit
      Else
      Ubest = U
      END IF
      END DO
      IF (  Ubest .LT. Uopt  ) THEN
      Uopt = Ubest
      sbeta1 = eps
      sbeta2 = eps
      End If
C     '
      beta1 = beta1prev - eps
      U = 1.E30
      beta2 = beta2prev
C     '   For beta2 = beta2prev To 0 Step -eps
      DO I = 1 , Int(beta2prev/eps)
      beta2 = beta2 - eps
      Ulast = U
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECoul(Z1,A1,beta1,
     *Z2,A2,beta2,d)
      IF (  U .GT. Ulast  ) THEN
      Exit
      Else
      Ubest = U
      End If
      END DO
      IF (  Ubest .LT. Uopt  ) THEN
      Uopt = Ubest
      sbeta1 = -eps
      sbeta2 = -eps
      END IF
C     '
      U = 1.E30
      beta2 = beta2prev
C     '   For beta2 = beta2prev To 1 Step eps
      DO I = 1 , Int((1-beta2prev)/eps)
      beta2 = beta2 + eps
      Ulast = U
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECoul(Z1,A1,beta1,
     *Z2,A2,beta2,d)
      IF (  U .GT. Ulast  ) THEN
      Exit
      Else
      Ubest = U
      END IF
      END DO
      IF (  Ubest .LT. Uopt  ) THEN
      Uopt = Ubest
      sbeta1 = -eps
      sbeta2 = eps
      END IF
C     '
C     '
      Ubest = Lymass(Z1,A1,beta1prev) + Lymass(Z2,A2,beta2prev) +
     *ECoul(Z1,A1,beta1prev,Z2,A2,beta2prev,d)
      U = Lymass(Z1,A1,beta1prev+REAL(sbeta1)) + Lymass(Z2,A2,
     *beta2prev+REAL(sbeta2)) + ECoul(Z1,A1,beta1prev+sbeta1,Z2,A2,
     *beta2prev+REAL(sbeta2),d)
C     '
C     '   L1:
      DO N = 1 , 1000
C     '
C     '   L2:
      DO N1 = 1 , N
      N2 = N-N1
      beta1 = beta1prev + sbeta1*N1
      beta2 = beta2prev + sbeta2*N2
      U = LyMass(Z1,A1,beta1) + LyMass(Z2,A2,beta2) + ECoul(Z1,A1,beta1,
     *Z2,A2,beta2,d)
      IF (  U .LT. Ubest  ) THEN
      Ubest = U
      beta1opt = beta1
      beta2opt = beta2
      Nopt = N
      END IF
      END DO
      IF (  N-Nopt .GT. 2  )  Exit
      END DO
C     '
C     '
      Else
C     ' Analytical approximation
C     ' Must be adapted if the relevant parameters of GEF are modified!
      xcoul = (Z1 + Z2)**2 / (A1 + A2)**(1.0/3.0)
      x = (Z1 / (Z1 + Z2))**(xcoul/xcoul236U)
      y = 1.2512E-4 + 0.00122851*x - 0.00267707*x**2 + 0.00372901*x**3
     *- 0.00219903*x**4
      beta1opt = y * xcoul
C     '
      x = (Z2 / (Z1 + Z2))**(xcoul/xcoul236U)
      y = 1.2512E-4 + 0.00122851*x - 0.00267707*x**2 + 0.00372901*x**3
     *- 0.00219903*x**4
      beta2opt = y * xcoul
C     '
      End If
C     '
      END
C     '
      REAL*4 FUNCTION U_Ired(Z,A)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
C     ' Effective moment of inertia by pairing with correction for excitation energy
      REAL*4  I_rigid_spher,IfragEff
C     '
      REAL*4 U_SHELL
C     '
      I_rigid_spher = 1.16E0**2 * A**1.6667E0 / 103.8415E0
C     '   IfragEff = I_rigid_spher + 0.003 * A^(4.0/3.0) * U_shell(Cint(Z),Cint(A))
C     '   IfragEff = I_rigid_spher + 0.005 * A^(4.0/3.0) * U_shell(Cint(Z),Cint(A))
C     ' reduction due to shell (Deleplanque et al. PRC 69 (2004) 044309)
      IfragEff = 0.45 * I_rigid_spher
C     ' Effect of superfluidity
C     '   IfragEff = 0.65 * IfragEff   ' Average effect of superfluidity and deformation
C     '
      U_Ired = IfragEff
      END
C     '
      REAL*4 FUNCTION U_IredFF(Z,A)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
C     ' Effective moment of inertia by pairing with correction for excitation energy
C     ' of final fission fragments
C     '
      REAL*4 U_Ired
      REAL*4 U_I_Shell
C     '
      U_IredFF = U_Ired(Z,A) * U_I_Shell(Z,A)
      END
C     '
      REAL*4 FUNCTION U_I_Shell(Z,A)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      INTEGER*4  N_shells(6)
C     ' Shell effect on the effective moment of inertia
      INTEGER*4  I
      REAL*4  dNmin,dZmin,dNsubmin
      REAL*4  Inv_add
      DATA Inv_add/0/
      REAL*4  I_inv_add_Z
      DATA I_inv_add_Z/0/
      REAL*4  I_inv_add_N
      DATA I_inv_add_N/0/
      REAL*4  I_inv_add_Nsub
      DATA I_inv_add_Nsub/0/
      N_shells(1) = 20
      N_shells(2) = 28
      N_shells(3) = 50
      N_shells(4) = 82
      N_shells(5) = 126
      N_shells(6) = 56
      dNmin = 100
      dZmin = 100
      dNsubmin = 100
      DO I = 1 , 5
      dZmin = Min(dZmin,Abs(N_shells(I) - Z))
      END DO
C     '
      DO I = 1 , 5
      dNmin = Min(dNmin,Abs(N_shells(I) - (A-Z)))
      END DO
C     '
      dNsubmin = Abs(N_shells(6) - (A-Z))
C     '
C     ' Effect of shells:
      IF (  dZmin .LT. 10.0  ) THEN
C     '        I_inv_add_Z = 0.33 * (6.0 * sqr(A/140.) - dZmin) * sqr(140./A)
      I_inv_add_Z = 0.33 * (6.0 * SQRT(A/140.0) - dZmin) *
     *(140.0/A)**1.5
C     ' A^(-1/3) dependence: "A simple phenomenology for 2gamma+ states",
C     ' N. V. Zamfir,D. Bucurescu,R. F. Casten,M. Ivascu,
C     ' Phys. Lett. B 241 (1990) 463
      I_inv_add_Z = Max(I_inv_add_Z,0.0)
      End If
      IF (  dNmin .LT. 10.0  ) THEN
C     '        I_inv_add_N = 0.42 * (8.0 * sqr(A/140.) - dNmin) * sqr(140./A)
      I_inv_add_N = 0.42 * (8.0 * SQRT(A/140.0) - dNmin) *
     *(140.0/A)**1.5
      I_inv_add_N = Max(I_inv_add_N,0.0)
      End If
      IF (  DNsubmin .LT. 6.0  ) THEN
C     '    I_inv_add_Nsub = 1.7 * (4.0 - dNsubmin) * (1.0 - 0.32 * Abs(40.0-Z))
      I_inv_add_Nsub = 1.7 * (4.0 - dNsubmin) * (1.0 - 0.18 *
     *Abs(40.0-Z))
C     ' N = 56 subshell only around Z = 40
      I_inv_add_Nsub = Max(I_inv_add_Nsub,0.0)
      End If
      U_I_shell = 1.0 / (1.0 + Max(I_inv_add_N,I_inv_add_Nsub) +
     *I_inv_add_Z)
C     'Print "*",I_inv_add_Z,I_inv_add_N,I_inv_add_Nsub,1.0 / (1.0 + Max(I_inv_add_N,I_inv_add_Nsub) + I_inv_add_Z)
      END
C     '
C     '
      REAL*4 FUNCTION U_alev_ld(Z,A)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
C     '  U_alev_ld = 0.073 * A + 0.095 * A^0.666667  'Ignatyuk (1970's)
      U_alev_ld = 0.078 * A + 0.115 * A**0.6666667
C     ' Ignatyuk (Bologna 2000)
C     '  U_alev_ld = 0.089 * A    ' only volume term
      END
C     '
      REAL*4 FUNCTION U_Temp(Z,A,E,Ishell,Ipair,Tscale,Econd)
      IMPLICIT NONE
      REAL*4 Z
      REAL*4 A
      REAL*4 E
      INTEGER*4 Ishell
      INTEGER*4 Ipair
      REAL*4 Tscale
      REAL*4 Econd
C     ' Temperature (modified Gilbert-Cameron composite level density)
C     ' KHS (10. 2. 2012)
      REAL*4  alev
      REAL*4  Eeff0,Eeff1,Rho0,Rho1,TCT,TFG
      REAL*4  fgamma
      DATA fgamma/0.055/
      REAL*4  RShell,RPair,Res
      REAL*4 U_ALEV_LD
      REAL*4 U_SHELL
      REAL*4 LYPAIR
      REAL*4 TEGIDY
C     ' Used global parameters: Tscale
C     '   alev = U_alev_ld(Z,A) * 1.1   ' Factor adjusted to high-energy prompt neutrons in U235(nth,f)
C     '  alev = U_alev_ld(Z,A) * 0.8  ' " with the correction for non-constant T (FG range)
      alev = U_alev_ld(Z,A)
C     '
      IF (  Ishell .EQ. 1  ) THEN
      RShell = U_Shell(NINT(Z),NINT(A))
      Else
      RShell = 0.0
      End If
      TCT = TEgidy(A,RShell,Tscale)
C     '
      IF (  Ipair .EQ. 1  ) THEN
      RPair = Lypair(NINT(Z),NINT(A))
      Else
      Rpair = 0.0
      End If
      Eeff0 = E - Econd + RPair + Rshell*(1.0 - exp(-fgamma * E))
C     '
      IF (  Eeff0 .GT. 0.5  ) THEN
      Eeff1 = Eeff0 + 0.1
      Rho0 = 1.E0/Eeff0**1.25 * exp(2.E0 * SQRT(alev * Eeff0))
      Rho1 = 1.E0/Eeff1**1.25 * exp(2.E0 * SQRT(alev * Eeff1))
C     '         Rho0 = 1.E0/Eeff0 * exp(2.E0 * sqr(alev * Eeff0))
C     '         Rho1 = 1.E0/Eeff1 * exp(2.E0 * sqr(alev * Eeff1))
      TFG = 0.1E0 / (log(Rho1) - log(Rho0))
      Else
      TFG = 0.0
      End If
      Res = TCT
      IF (  TFG .GT. Res  )  Res = TFG
C     '
C     ' If Res > 1.4 Then Res = 1.4
C     '
      U_Temp = Res
      END
C     '
      REAL*4 FUNCTION U_Even_Odd(I_Channel,PEO)
      IMPLICIT NONE
      INTEGER*4 I_Channel
      REAL*4 PEO
C     ' Creates even-odd fluctuations
      REAL*4  R
      IF (   MOD(I_Channel,2)  .EQ. 0  ) THEN
      R = 1.0 + PEO
      Else
      R = 1.0 - PEO
      End If
      U_Even_Odd = R
      END
C     '
C     '
      REAL*4 FUNCTION BFTF(RZ,RA,I_Switch)
      IMPLICIT NONE
      REAL*4 RZ
      REAL*4 RA
      INTEGER*4 I_Switch
C     /' Fission barriers from Myers and Swiatecki,Thomas-Fermi model '/
C     /'  I_Switch: 0: liquid-drop; 1: with shells and pairing,
C     2: averaged over pairing,3: with shell and pairing + pairing gap at barrier '/
C     ' 4: liquid-drop + g.s. shell,no Z correction
      REAL*4  RN,RI,Rkappa,RS,RF,RX
      REAL*4  RX0
      DATA RX0/48.5428/
      REAL*4  RX1
      DATA RX1/34.15/
      REAL*4  RB
      INTEGER*4  IZ,IA
      REAL*4 U_SHELL
      REAL*4 U_SHELL_EXP
      REAL*4 U_SHELL_EO_EXP
      REAL*4 LYPAIR
C     '
      IZ = NINT(RZ)
      IA = NINT(RA)
      RN = RA - RZ
      RI = (RN-RZ) / RA
      Rkappa = 1.9E0 + (RZ - 80.E0) / 75.E0
      RS = RA**0.666667E0 * (1.E0 - Rkappa * RI**2)
      RX = RZ**2 / (RA * (1.E0 - Rkappa * RI**2))
      IF (  RX .LT. 30  ) THEN
C     /' out of range '/
      RF = 1.E10
      End If
      IF (  RX .GT. RX0  ) THEN
C     /' out of range '/
      RF = 0.0
      End If
      IF (  RX .LT. RX1 .AND. RX .GT. 30  ) THEN
      RF = 0.595553E0 - 0.124136E0 * (RX - RX1)
      End If
      IF (  RX .GE. RX1 .AND. RX .LE. RX0  ) THEN
      RF = 0.000199749 * (RX0 - RX)**3
      End If
      RB = RF * RS
C     '
      Select CASE( I_Switch)
      CASE( 0)
      BFTF = RB
      CASE( 1)
C     ' including even-odd staggering due to increased pairing strength at barrier
C     ' Tentative modification from comparison with experimental fission barriers
C     ' (shell correction at the barrier?)
      IF (  RZ .GT. 86.5  )  RB = RB - 0.15 * (RZ - 86.5)
C     '    If RZ > 90 Then RB = RB + 0.3 * (RZ - 90.0)
C     '    If RZ > 98 Then RB = RB - 0.15 * (RZ - 98.0)
      IF (  RZ .GT. 90  )  RB = RB + 0.35 * (RZ - 90.0)
      IF (  RZ .GT. 93  )  RB = RB + 0.15 * (RZ - 93.0)
      IF (  RZ .GT. 95  )  RB = RB - 0.25 * (RZ - 95.0)
C     '    BFTF = RB - U_Shell(IZ,IA)
C     '    BFTF = RB - U_Shell_exp(IZ,IA)
      BFTF = RB - U_Shell_EO_exp(IZ,IA) + Lypair(IZ,IA) * 14.0/12.0
      CASE( 2)
C     ' averaged over even-odd staggering
      IF (  RZ .GT. 86.5  )  RB = RB - 0.15 * (RZ - 86.5)
      IF (  RZ .GT. 90  )  RB = RB + 0.35 * (RZ - 90.0)
      IF (  RZ .GT. 93  )  RB = RB + 0.15 * (RZ - 93.0)
      IF (  RZ .GT. 95  )  RB = RB - 0.25 * (RZ - 95.0)
      BFTF = RB - U_Shell_exp(IZ,IA)
      CASE( 3)
C     ' like Case 1 + pairing gap at barrier
      IF (  RZ .GT. 86.5  )  RB = RB - 0.15 * (RZ - 86.5)
      IF (  RZ .GT. 90  )  RB = RB + 0.35 * (RZ - 90.0)
      IF (  RZ .GT. 93  )  RB = RB + 0.15 * (RZ - 93.0)
      IF (  RZ .GT. 95  )  RB = RB - 0.25 * (RZ - 95.0)
      BFTF = RB - U_Shell_EO_exp(IZ,IA)
      CASE( 4)
C     ' like case 3 but without Z correction
C     ' This is the direct description from the topographic theorem.
      BFTF = RB - U_Shell_exp(IZ,IA)
      Case DEFAULT
C     '         Print "Undefined option in BFTF"
C     '         Sleep
      End Select
C     /'  If I_Switch = 0 Then
C     BFTF = RB
C     Else
C     ' Tentative modification from comparison with experimental fission barriers
C     ' (shell correction at the barrier?)
C     If RZ > 86.5 Then RB = RB - 0.15 * (RZ - 86.5)
C     '    If RZ > 90 Then RB = RB + 0.3 * (RZ - 90.0)
C     '    If RZ > 98 Then RB = RB - 0.15 * (RZ - 98.0)
C     If RZ > 90 Then RB = RB + 0.35 * (RZ - 90.0)
C     If RZ > 93 Then RB = RB + 0.15 * (RZ - 93.0)
C     If RZ > 95 Then RB = RB - 0.25 * (RZ - 95.0)
C     '
C     '    BFTF = RB - U_Shell(IZ,IA)
C     '    BFTF = RB - U_Shell_exp(IZ,IA)
C     BFTF = RB - U_Shell_EO_exp(IZ,IA) + Lypair(IZ,IA) * 14.0/12.0
C     End If '/
      END
C     '
      REAL*4 FUNCTION BFTFA(RZ,RA,I_Switch)
      IMPLICIT NONE
      REAL*4 RZ
      REAL*4 RA
      INTEGER*4 I_Switch
C     /' inner barrier height '/
      REAL*4  EA,BF0,Z4A,Z3A,DB
      REAL*4  coeff
      DATA coeff/0.5/
      REAL*4 BFTF
      BF0 = BFTF(RZ,RA,I_Switch)
C     ' Z4A = RZ^4 / RA
C     '  EB - EA from fit to Smirenkin barriers:
C     '  V. M. Kupriyanov,K. K. Istekov,B. I. Fursov,G. N. Smirenkin
C     '  Sov. J. Nucl. Phys. 32 (1980) 184
C     '  DB = -10.3517 + 1.6027E-5 * Z4A + 5.4945E-11 * Z4A^2  ' EA - EB
C     '
C     '  EB - EA from fit to data from Dahlinger et al. (KHS,21. Dec. 2012)
      Z3A = RZ**3 / RA
      DB = -(5.40101 - 0.00666175*Z3A + 1.52531E-6*Z3A**2)
      IF (  DB .GT. 0.0  ) THEN
      EA = BF0 - DB
      Else
      EA = BF0
      End If
      BFTFA = EA
      END
C     '
      REAL*4 FUNCTION BFTFB(RZ,RA,I_Switch)
      IMPLICIT NONE
      REAL*4 RZ
      REAL*4 RA
      INTEGER*4 I_Switch
C     /' outer barrier height '/
      REAL*4  EB,BF0,Z4A,Z3A,DB
      REAL*4  coeff
      DATA coeff/0.5/
      REAL*4 BFTF
      BF0 = BFTF(RZ,RA,I_Switch)
C     ' Z4A = RZ^4 / RA
C     '  EB - EA from fit to Smirenkin barriers:
C     '  V. M. Kupriyanov,K. K. Istekov,B. I. Fursov,G. N. Smirenkin
C     '  Sov. J. Nucl. Phys. 32 (1980) 184
C     '   DB = -10.3517 + 1.6027E-5 * Z4A + 5.4945E-11 * Z4A^2  ' EA - EB
C     '
C     '  EB - EA from fit to data from Dahlinger et al. (KHS,21. Dec. 2012)
      Z3A = RZ**3 / RA
      DB = -(5.40101 - 0.00666175*Z3A + 1.52531E-6*Z3A**2)
      IF (  DB .LT. 0.0  ) THEN
      EB = BF0 + DB
      Else
      EB = BF0
      End If
      BFTFB = EB
      END
C     '
C     '
C     '
C     /' Utility functions '/
C     '
C     '
      REAL*4 FUNCTION Gaussintegral(R_x,R_sigma)
      IMPLICIT NONE
      REAL*4 R_x
      REAL*4 R_sigma
C     /' Smoothed step function. Grows from 0 to 1 around R_x
C     with a Gauss-integral function with given sigma'/
      REAL*4  R_ret
C     ' Note: The variable R_sigma = standard deviation / sqr(2) !
      REAL*4 ERF
      R_ret = 0.5E0 + 0.5E0 * Erf(R_x / R_sigma)
      Gaussintegral = R_ret
      END
C     '
      REAL*4 FUNCTION U_Box(x,sigma,length)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 sigma
      REAL*4 length
      REAL*4  y
C     ' Note: The variable sigma = standard deviation / sqr(2) !
      REAL*4 GAUSSINTEGRAL
      y = Gaussintegral(x+0.5*length,sigma) -
     *Gaussintegral(x-0.5*length,sigma)
      U_Box = y/length
      END
C     '
      REAL*4 FUNCTION U_Box2(x,sigma1,sigma2,length)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 sigma1
      REAL*4 sigma2
      REAL*4 length
      REAL*4  y
C     ' Note: The variable sigma = standard deviation / sqr(2) !
      REAL*4 GAUSSINTEGRAL
      y = Gaussintegral(x+0.5*length,sigma2) -
     *Gaussintegral(x-0.5*length,sigma1)
      U_Box2 = y/length
      END
C     '
      REAL*4 FUNCTION U_Gauss(x,sigma)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 sigma
      REAL*4  y
      REAL*4  pi
      PARAMETER (pi=3.14159)
C     '
      y = 1.0 / (SQRT(2.0 * pi) * sigma) * exp(-x**2/ ( 2.0 * sigma**2
     *) )
      U_Gauss = y
      END
C     '
      REAL*4 FUNCTION U_Gauss_mod(x,sigma)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 sigma
C     ' Gaussian with Sheppard correction
      REAL*4  y
      REAL*4  sigma_mod
      REAL*4  pi
      PARAMETER (pi=3.14159)
      sigma_mod = SQRT(sigma**2 + 1./12.)
C     '
      y = 1.0 / (SQRT(2.0 * pi) * sigma_mod) * exp(-x**2/ ( 2.0 *
     *sigma_mod**2 ) )
      U_Gauss_mod = y
      END
C     '
      REAL*4 FUNCTION U_LinGauss(x,R_Sigma)
      IMPLICIT NONE
      REAL*4 x
      REAL*4 R_Sigma
C     /' Gaussian times a linear function '/
C     /' Not normalized! '/
      REAL*4  R_Res
      IF (  R_Sigma .GT. 0.0  ) THEN
      R_Res = x * exp(-x**2/(2.0 * R_Sigma**2))
      Else
      R_Res = 0.0
      End If
      U_LinGauss = R_Res
      END
