C     Output of FBtoFO from GEFSUB.BAS
CAK   PROGRAM MAIN
CAK   IMPLICIT NONE
C     '
C     '
C     '    Copyright 2009,2010,2011,2012,2013,2014,2015,2016:
C     '       Dr. Karl-Heinz Schmidt,Rheinstrasse 4,64390 Erzhausen,Germany
C     '       and
C     '       Dr. Beatriz Jurado,Centre d'Etudes Nucleaires de Bordeaux-Gradignan,
C     '       Chemin du Solarium,Le Haut Vigneau,BP 120,33175 Gradignan,Cedex,
C     '       France
C     '
C     '    This program is free software: you can redistribute it and/or modify
C     '    it under the terms of the GNU General Public License as published by
C     '    the Free Software Foundation,either version 3 of the License,or
C     '    (at your option) any later version.
C     '
C     '    This program is distributed in the hope that it will be useful,
C     '    but WITHOUT ANY WARRANTY; without even the implied warranty of
C     '    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     '    GNU General Public License for more details.
C     '
C     '    You should have received a copy of the GNU General Public License
C     '    along with this program.  If not,see <http://www.gnu.org/licenses/>.
C     '
C     '
C     /' Documentation: '/
C     /' (1) K.-H. Schmidt and B. Jurado,Contribution to '/
C     /'     ESNT Workshop "The scission process",Saclay (France),April 12-16,2010 '/
C     /' (2) B. Jurado and K.-H. Schmidt,Contribution to '/
C     /'     Seminar an fission,Gent (Belgium),May 17-20,2010 '/
C     /' (3) B. Jurado and K.-H. Schmidt,Contribution to '/
C     /'     EFNUDAT Workshop,Paris (France),May 25-27,2010 '/
C     /' (4) K.-H. Schmidt and B. Jurado,Contribution to '/
C     /'     EFNUDAT Workshop,Paris (France),May 25-27,2010 '/
C     /' (5) K.-H. Schmidt and B. Jurado,'/
C     /'     Final Report to EFNUDAT,October,2010 '/
C     /' (6) K.-H. Schmidt and B. Jurado,Phys. Rev. Lett. 104 (2010) 21250 '/
C     /' (7) K.-H. Schmidt and B. Jurado,Phys. Rev. C 82 (2011) 014607 '/
C     /' (8) K.-H. Schmidt and B. Jurado,Phys. Rev. C 83 (2011) 061601 '/
C     /' (9) K.-H. Schmidt and B. Jurado,arXiv:1007.0741v1[nucl-th] (2010) '/
C     /' (10) K.-H. Schmidt and B. Jurado,JEF/DOC 1423,NEA of OECD,2012 '/
C     /' (11) K.-H. Schmidt and B. Jurado,Phys. Rev. C 86 (2012) 044322 '/
C     /' (12) K.-H. Schmidt,B. Jurado,Ch. Amouroux,JEFF-Report 24,NEA of OECD,2014 '/
C     /' (13) B. Jurado,K.-H. Schmidt,J. Phys. G: Nucl. Part. Phys. 42 (2015) 055101 '/
C     /' (14) K.-H. Schmidt,B. Jurado,Eur. Phys. J. A 51 (2015) 176 '/
C     /' (15) K.-H. Schmidt,B. Jurado,C. Amouroux,C. Schmitt,Nucl. Data Sheets 131 (2016) 107 '/
C     '
C     '
C     /' Further documentation and the newest version of the GEF code are '/
C     /' available from                                                   '/
C     /' http://www.cenbg.in2p3.fr/GEF and http://www.khs-erzhausen.de/ . '/
C     '
C     '
C     '    The development of the GEF code has been supported by the European Union,
C     '    EURATOM 6 in the Framework Program "European Facilities for Nuclear Data
C     '    Measurements" (EFNUDAT),contract number FP6-036434,the Framework
C     '    Program "European Research Infrastructure for Nuclear Data Applications
C     '    (ERINDA),contract number FP7-269499,and by the OECD Nuclear Energy Agency.
C     '
C     '
C     ' Technical remark: The code contains commented sections with
C     ' produce a deterministic version of the GEF code as a subroutine from this
C     ' source with a dedicated pre-processor,named ExtractSub.bas.
C     ' are used in the deterministic GEF version.
C     ' ExtractSub.bas produces the deterministic code GEFSUB.bas in FreeBASIC.
C     ' The translator FBtoFO.bas can be used to translate this code into FORTRAN99.
C     ' In order to enable this translation,GEF.bas contains already some additional
C     ' that are only used by the translator to provide some specific
C     ' FORTRAN statements that cannot be produced automatically.
C     '
C     '
C     /' K.-H. Schmidt / B. Jurado,07/Feb./2009 '/
C     /' SEFI9 is taken as a basis and extended by new features in SEFI14 (May 2010),KHS '/
C     /' Several improvements (even-odd effect,charge polarization etc. (June 2010),KHS '/
C     /' SEFI15 converted from PL/I to FreeBASIC on 04/July/2010,KHS '/
C     /' Error in LyPair corrected (26/July/2010) KHS '/
C     /' Indices corrected in U_Shell inside Eva (1/Aug/2010) KHS '/
C     /' Major developments,sigma_E*(scission),sigma_Z(A) etc. (14/Aug/2010) KHS '/
C     /' Macroscopic masses from Thomas-Fermi model (Myers & Swiatecki) (17/Aug/2010) KHS '/
C     /' 3 reference options for energy input (4/Sept/2010) KHS '/
C     /' Graphic output of mass distribution added (if there are problems with the X11
C     installation on LINUX,the graphics may be suppressed by simply commenting
C     the line  -> #Include Once "DCLPlotting.bas" <- )  (5/Sept/2010) KHS '/
C     /' Comparison with ENDF compilation in graphic output (12/Sept/2010) KHS '/
C     /' Super-long fission channel included (14/Sept/2010) KHS '/
C     /' Overlap of S1 and S2 fission channels in both fragments included (18/Oct/2010) KHS '/
C     /' Output of neutron multiplicity distribution added (20/Oct/2010) KHS '/
C     /' Decreasing curvature of shells with increasing E* (24/Oct/2010) KHS '/
C     /' Angular momenta of fission fragments added (18/Dec/2010) KHS '/
C     /' CN angular momentum considered (14/Jan/2011) KHS '/
C     /' Numerical stability improved (28/Jan/2011) KHS '/
C     /' Output in ENDF format (optional) (4/Feb/2011) KHS '/
C     /' Input list from file supported (31/Jan/2011) KHS '/
C     /' Multiprocessing supported (5/Feb/2011) KHS '/
C     /' Polarization for symmetric fission channel improved (12/Feb/2011) KHS '/
C     /' GUI for input (23/Feb/2011) KHS '/
C     /' Calculation of fission-fragment angular momentum refined (5/May/2011) KHS '/
C     /' Neutron inverse cross section modified (5/August/2011) KHS '/
C     /' Even-odd staggering in neutron emission improved (5/August/2011) KHS '/
C     /' Slight modifications in angular-momentum distributions (19/October/2011) KHS '/
C     /' Gamma emission added (23/November/2011) KHS '/
C     /' TKE added (24/November/2011) KHS '/
C     /' Neutron spectrum added (25/November/2011) KHS '/
C     /' Neutron-gamma competition added (4/December/2011) KHS '/
C     /' Composite level density (Egidy + Ignatyuk) refined (31/January/2012) KHS '/
C     /' Treatment of GDR refined (14/February/2012) KHS '/
C     /' Deformation of S3 channel changed (27/February/2012) KHS '/
C     /' Z=44 deformed shell added (supports S1 around Pu) (27/February/2012) KHS '/
C     /' Uncertainties from perturbed fission yields (3/March/2012) KHS '/
C     /' Validity range extended to Z=120 (with a warning message) (8/March/2012) KHS '/
C     /' Neutron emission during fragment acceleration (18/March/2012) KHS '/
C     /' Several optimizations (15/April/2012) KHS '/
C     /' TF masses of Myers & Swiatecki corrected (pairing shift) (29/May/2012) KHS '/
C     /' Correction on intrinsic excitation energy (05/June/2012) KHS '/
C     /' Correction on gamma emission (25/September/2012) KHS '/
C     /' Free choice of listmode values (13/October/2012) KHS '/
C     /' Excitation-energy distribution from file (14/October/2012) KHS '/
C     /' Transfer of input from GUI corrected (02/November/2012) KHS '/
C     /' Parameters of perturbed calculations modified (02/November/2012) KHS '/
C     /' Model parameters je-adjusted (08/November/2012) KHS '/
C     /' Input options for isomeric target nuclei (09/November/2012) KHS '/
C     /' Random initialisation of the random generator (26/November/2012) KHS '/
C     /' Covariance matrix for Z,Apre,Apost,ZApre,ZApost (06/December/2012) KHS '/
C     /' Output file in XML format (06/December/2012) KHS '/
C     /' Multi-chance fission supported (13/December/2012) KHS '/
C     /' Pre-compound emission for (n,f) included (13/December/2012) KHS '/
C     /' Some technical corrections and modifications (17/December/2012) KHS '/
C     /' Transition from asymmetric to symmetric fission around Fm improved (19/December/2012) KHS '/
C     /' Influence of S2 channel on S1 channel in the other fragment included (20/December/2012) KHS'/
C     /' List-mode output of pre-fission neutron energies (21/December/2012) KHS '/
C     /' Parameterisation for EB-EA from fit to data in Dahlinger et al. (21/December/2012) KHS '/
C     /' Fission channel at Z=42 added (seen around Pu in light fragment and around Hg) (23/Dec./2012) KHS '/
C     /' Gamma-n / Gamma-f according to Moretto (IAEA Rochester) (23/December/2012) KHS '/
C     /' Pre-compound neutron energies modified (24/December/2012) KHS '/
C     /' Some technical revisions to avoid crashes in covariances (26/December/2012) KHS '/
C     /' Influence of shells on yrast line from Deleplanque et al. (26/December/2012) KHS '/
C     /' Fission threshold in multi-chance fission modified (30/December/2012) KHS '/
C     /' Output of energies at fission for multi-chance fission (30/December/2012) KHS '/
C     /' Several revisions (15/January/2013) KHS '/
C     /' New optical model fit (3/February/2013) KHS '/
C     /' Gamma-f / Gamma-n modified (3/February/2013) KHS '/
C     /' Handling for reading input from file corrected (5/February/2013) KHS '/
C     /' Data transfer from GUI corrected (6/February/2013) KHS '/
C     /' Input dialog re-organized (6/February/2013) KHS '/
C     /' Mass-dependent deformation and charge polarization revised (9/February/2013) KHS '/
C     /' Calculation of combined fission channels S12,S22 revised (10/February/2013) KHS '/
C     /' Extension of validity range to heavier nuclei (10/February/2013) KHS '/
C     /' Improved description of prompt-neutron spectra (21/February/2013) KHS '/
C     /' Several technical corrections and developments (April-May/2013) KHS '/
C     /' Pre-fission emission of protons considered (14/Mai,2013) KHS '/
C     /' Structure of energy list in input file modified (26/May/2013) KHS '/
C     /' Option "local fit" added (26/May/2013) KHS '/
C     /' Neutron evaporation corrected (more realistic even-odd effect in isotonic distr.) (21/June(2013) KHS '/
C     /' Even-odd effect in TKE added (23/June/2013) KHS '/
C     /' Calculation of Z-A-covariance matrix corrected (16/July/2013) KHS '/
C     /' Output of multi-variant distributions corrected for input from file (24/July/2013) KHS '/
C     /' New global fit (most important: Energy gain from saddle to scission reduced) (15/September/2013) KHS '/
C     /' Even-odd effect in neutron number of fragments modified (17/September/2013) KHS '/
C     /' Curvatures of fission valleys adjusted to experimental shells around 132Sn (18/September/2013) KHS '/
C     /' Width of S0 corrected: Fit of Rusanov (18/September/2013) KHS '/
C     /' Random generator Box with asymmetric diffuseness for S2 (18/September/2013) KHS '/
C     /' Gaussian random generator revised (20/September/2013) KHS '/
C     /' Mass shift of fission channels with E* modified (22/September/2013) KHS '/
C     /' Energy dependence of S1 position corrected,slightly modified parameters (25/September/2013) KHS '/
C     /' Initial angular momentum introduced as an input parameter (02/October/2013) KHS '/
C     /' Calculation of prompt-neutron emission improved,some model parameters modified (12/October/2013) KHS '/
C     /' Technical error,causing incomplete covariance matrices on output corrected  (17/October/2013) KHS '/
C     /' Washing of shell effects considerd in shape fluctuations (26/October/2013) KHS '/
C     /' Post-scission neutrons added to list-mode output (8/November/2013) KHS '/
C     /' Fission-gamma competition refined (10/November/2013) KHS '/
C     /' New global fit of model parameters (18/November/2013) KHS '/
C     /' Multi-chance fission modified (18/November/2013) KHS '/
C     /' A numerical instability removed (20/November/2013) KHS '/
C     /' Calculation of multi-chance fission: corrected and modified (28/November/2013) KHS '/
C     /' Spurious even-odd effect in fission probabilities (from Moeller's shells) removed (4/December/2013) KHS '/
C     /' More precise calculation of E* in n-induced fission (6/December/2013) KHS '/
C     /' Pre-fission gamma strength increased (adjusted to Pf(E*) of 238U (22/December/2013) KHS '/
C     /' Kinetic energy of prompt neutrons in cm system (fragment-neutron) (10/Janurary/2014) KHS '/
C     /' Binnig of prompt-neutron spectra corrected (shift by 50 keV) (10/January/2014) KHS '/
C     /' Description of fission barriers modified (15/January/2014) KHS '/
C     /' Description of energy-dependent fission probability modified (17/January/2014) KHS '/
C     /' E* at scission modified for fission below Bf (4/February/2014) KHS '/
C     /' Influence of Z=44 shell on deformation of light fragment (4/February/2014) KHS '/
C     /' Tunneling for S1 enhanced (16/February/2014) KHS '/
C     /' New global fit of fission channels (strength and position) (19/March/2014) KHS '/
C     /' Normalization of distributions for deterministic version corrected (23/May/2014) KHS '/
C     /' Calculation of "corrected sample variance" revised (05/September/2014) KHS '/
C     /' Uncertainties calculated with "corrected sample variance" (05/September/2014) KHS '/
C     /' Shape of "S22" fission channel modified (08/December/2014) KHS '/
C     /' Binsize for prompt neutrons set to 1 keV (09/December/2014) KHS '/
C     /' Binsize for prompt gammas set to 1 keV (09/December/2014) KHS '/
C     /' Fragment excitation energies in listmode (10/December/2014) KHS '/
C     /' Shell in symmetric fission channel included in GEFSUB (12/December/2014) KHS '/
C     /' Ground-state spin of fragments is taken into account (22/January/2015) KHS '/
C     /' Pre-fission gamma emission included (20/February/2015) KHS '/
C     /' VMI model for E2 gammas developed and implemented (15/March/2015) KHS '/
C     /' Binning of spectrum EN corrected (28/March/2015) KHS '/
C     /' Energetics of symmetric fission channel revised (5/April/2015) KHS '/
C     /' Output of isomeric yields extended (# of events added) (10/April/2015) KHS '/
C     /' Control-output for progress of calculation (17/April/2015) KHS '/
C     /' Angular-momentum dependence of pairing gap in eva (7/May/2015) KHS '/
C     /' Angular-momentum dependent Gf: influence of fissility and temperature (13/May/2015) KHS '/
C     /' Influence of angular momentum on mass width by increasing fissiliy (14/May/2015) KHS '/
C     /' Mass distribution truncated at zero and A_CN (14/May/2015) KHS '/
C     /' Yields of fission modes also for multi-chance fission (7/August/2015) KHS '/
C     /' Tentative description of structure in FF distributions around A_CN = 180-200 (8/August/2015) KHS '/
C     /' Improved description of asymmetric fission in light nuclei (A_CN < 220) (9/August/2015) KHS '/
C     /' Extension of the validity range to lighter nuclei (9/August/2015) KHS '/
C     /' Adaptions of ERF function to new FreeBASIC compiler (12/Sept(2015) KHS '/
C     /' Allowance for calculations with very large number of events (15/Sept(2015) KHS '/
C     /' Slight adjustement of parameters to "repair" deterioation for low E* (21/Sept/2015) KHS '/
C     /' Correction of some formatting issues (03/Oct/2015) KHS '/
C     /' Calculation of covariances for independent yields corrected (6/Nov/2015) KHS '/
C     /' Calculation of covariances between two systems corrected (6/Nov/2015) KHS '/
C     /' Output of correlations matrices added (6/Nov/2015) KHS '/
C     /' Espectrum.in can now also provide the CN spin (13/Nov/2015) KHS '/
C     /' Output for parallel computing better protected (1/March/2016) KHS '/
C     /' Output of random files provided (4/March/2016) KHS '/
C     /' Prompt-neutron spectrum with variable bin size added (6/April/2016) KHS '/
C     /' Overflow problem in tanh and coth corrected (11/April/2016) KHS '/
C     /' Experimental masses in evaporation routine (not yet finished) (23/April/2016) KHS '/
C     /' Subscriptrange problem in output of neutron energies corrected (29/April/2016) KHS '/
C     /' Gamma spectra with conditions on pre-neutron mass (29/Mai/2016) KHS '/
C     /' Relation between energies: experimental masses used (07/Sept/2016) KHS '/
C     /' Partial support of p-induced fission (up to Ep = 30 MeV)
C     (pre-equilibrium emission not yet implemented,preliminary l-distribution) (07/Sept/2016) KHS '/
C     /' Model parameters of fission channels optimized (20/Sept/2016) KHS '/
C     /' Neutrons emitted between outer saddle and scission added (21/Sept/2016) KHS '/
C     /' Output of delayed-neutron multiplicities added (25/Sept/2016) KHS '/
C     /' New adjustment of model parameters (13/Oct/2016) KHS '/
C     /' Output of delayed-neutron emitters added (13/Oct/2016) KHS '/
C     /' Uncertainty range of angular momentum added (24/Oct/2016) KHS '/
C     /' Uncertainties of prompt-gamma and prompt-neutron characteristics added (24/Oct/2016) KHS '/
C     /' Uncertainty of TKE (pre and post) added (06/Nov/2016) KHS '/
C     /' Normalization of Maxwell distribution corrected (07/Nov/2016) KHS '/
C     '
C     /' FreeBASIC is available from http://www.freebasic.net/ '/
C     /' FreeBASIC runs on Windows,Linux,and DOS. '/
C     /' FreeBASIC compiles a binary code that uses the C run-time library. '/
C     '
C     '
C     '  #Include "utilities.bi"
C     '
C     '
C     /' Functions and subroutines '/
C     '
CAK
CAK   CALL GEFSUB(92,236,6.,0.0)
C     '
C     /'
C     Dim As Single Zsum
C     '
C     'Print
C     'Print "Z,A,Yield"
C     For J = 10 To 140
C     Zsum = 0
C     For I = 10 To 190
C     Zsum = Zsum + NZpre(I,J)
C     If NZPRE(I,J) > 0.0 Then
C     '      Print J,I+J,NZPRE(I,J)*200
C     End If
C     Next
C     Next
C     '
C     '
C     'Print
C     'Print "Z yields"
C     For J = 10 To 140
C     Zsum = 0
C     For I = 10 To 190
C     Zsum = Zsum + NZpre(I,J)
C     Next
C     If Zsum > 0.0 Then
C     '    Print J,Zsum * 200
C     End If
C     Next '/
C     '
C     '
C     /'
C     Dim As Single Asum
C     'Print
C     'Print "N yields"
C     For I = 10 To 140
C     Asum = 0
C     For J = 10 To 190
C     Asum = Asum + NZpre(I,J)
C     '   If NZPRE(I,J) > 0.001 Then
C     '     Print J,I+J,NZPRE(I,J)*200
C     '   End If
C     Next
C     '  Print I,Asum * 200
C     Next   '/
C     '
CAK   End
C     '
C     '
      SUBROUTINE GEFSUB(P_Z_CN,P_A_CN,P_E_EXC,P_J_CN)
      IMPLICIT NONE
CAK
      include "gef.cmb"
      INTEGER*4 P_Z_CN
      INTEGER*4 P_A_CN
      REAL*4 P_E_EXC
      REAL*4 P_J_CN
C     /' Input parameters: '/
C     /' Atomic number,mass number,excitation energy/MeV,spin/h_bar of CN '/
C     /' Results are stored in external arrays. '/
C     '
CAK
CAK   INCLUDE "GEFSUBdcl2.FOR"
      include "gefsubdcl2.cmb"
CAK
      Beta = 0.
      EPART = 0.
CAK
      xDelta_S0 = U_Delta_S0(P_Z_CN,P_A_CN)
C     ' default values
C     '
C     ' Use nominal parameter values:
      P_DZ_Mean_S1 = xP_DZ_Mean_S1
      P_DZ_Mean_S2 = xP_DZ_Mean_S2
      P_DZ_Mean_S3 = xP_DZ_Mean_S3
      P_DZ_Mean_S4 = xP_DZ_Mean_S4
      P_Z_Curv_S1 = xP_Z_Curv_S1
      P_Z_Curv_S2 = xP_Z_Curv_S2
      P_A_Width_S2 = xP_A_Width_S2
      P_Z_Curv_S3 = xP_Z_Curv_S3
      P_Z_Curv_S4 = xP_Z_Curv_S4
      Delta_S0 = xDelta_S0
      P_Shell_S1 = xP_Shell_S1
      P_Shell_S2 = xP_Shell_S2
      P_Shell_S3 = xP_Shell_S3
      P_Shell_S4 = xP_Shell_S4
      T_low_S1 = xT_low_S1
      T_low_S2 = xT_low_S2
      T_low_S3 = xT_low_S3
      T_low_S4 = xT_low_S4
      T_low_SL = xT_low_SL
      P_att_pol = xP_att_pol
      HOMPOL = xHOMPOL
      POLARadd = xPOLARadd
      Jscaling = xJscaling
C     '
      R_E_exc_used = P_E_exc
      I_A_CN = P_A_CN
      I_Z_CN = P_Z_CN
C     '
C     /' Central Z values of fission modes '/
C     '
C     /' Fit to positions of fission channels (Boeckstiegel et al.,2008) '/
C     /' P_DZ_Mean_S1 and P_DZ_Mean_S2 allow for slight adjustments '/
C     '    Scope
      R_Z_mod = I_Z_CN
      ZC_Mode_0 = R_Z_mod * 0.5E0
C     /' Central Z value of SL mode '/
      ZC_Mode_1 = (53.0E0 - 51.5E0) / (1.56E0 - 1.50E0) *
     *(R_Z_mod**1.3E0 / I_A_CN - 1.50E0) + 51.5E0 + P_DZ_Mean_S1
      ZC_Mode_2 = (55.8E0 - 54.5E0) / (1.56E0 - 1.50E0) *
     *(R_Z_mod**1.3E0 / I_A_CN - 1.50E0) + 54.5E0 + P_DZ_Mean_S2
      ZC_Mode_3 = ZC_Mode_2 + 4.5E0 + P_DZ_Mean_S3
C     '  ZC_Mode_4 = 38.5 + P_DZ_Mean_S4  ' structure in nuclei with A around 190 for 201Tl
C     '  ZC_Mode_4 = 35.5 + P_DZ_Mean_S4  ' for 180Hg  ( 36.2 for 208Po )
C     '
C     ' Do not delete these lines (,because this is a very good fit!):
C     '    ZC_Mode_4 = 38.5 + (I_A_CN-I_Z_CN-110)*0.12 - (I_A_CN-I_Z_CN-110)^2 * 0.009   '                - (I_Z_CN-77)*0.34 + P_DZ_Mean_S4
C     '
      ZC_Mode_4 = 38.5 + (I_A_CN-I_Z_CN-110)*0.12 -
     *(I_A_CN-I_Z_CN-110)**2 * 0.009 - (I_Z_CN-77)*0.34 + P_DZ_Mean_S4
C     ' assumption: mode position moves with Z and A (adjusted to exp. data
C     ' of Itkis and Andreyev et al.
C     '
C     '    End Scope
C     '
C     '
C     '
      I_N_CN = I_A_CN - I_Z_CN
CAK
      if (I_N_CN.gt.203) then
        write(*,'(" TALYS-error: GEF Z=",i3," N=",i3," A=",i3,
     +    " out of bounds")') I_Z_CN, I_N_CN, I_A_CN
        stop
      endif
C     /' Mean deformation at scission as a function of mass '/
C     '
C     /' Mode 0: liquid drop and mode 4: Z = 38 '/
      beta1_prev = 0.3
      beta2_prev = 0.3
      beta1_opt = beta1_prev
      beta2_opt = beta2_prev
      DO I = 10 , I_Z_CN - 10
      IZ1 = I
      Z1 = REAL(IZ1)
      IZ2 = I_Z_CN - IZ1
      Z2 = REAL(IZ2)
      A1 = Z1 / REAL(I_Z_CN) * REAL(I_A_CN)
      A2 = I_A_CN - A1
C     '
      CALL Beta_Equi(A1,A2,Z1,Z2,dneck,beta1_prev,beta2_prev,beta1_opt,
     *beta2_opt)
C     '
C     'Print "Mode 0,Z1,Z2,beta1,beta2 ";Z1;" ";Z2;" ";beta1_opt,beta2_opt
C     'Print Z1;" ";Z2;" ";beta1_opt,beta2_opt
      Beta(0,1,IZ1) = beta1_opt
C     /' "light" fragment '/
      Beta(4,1,IZ1) = beta1_opt
      Beta(0,2,IZ2) = beta2_opt
C     /' "heavy" fragment '/
      Beta(4,2,IZ2) = beta2_opt
      beta1_prev = beta1_opt
      beta2_prev = beta2_opt
      E_defo = Lymass(Z1,A1,beta1_opt) - Lymass(Z1,A1,0.0)
      Edefo(0,1,IZ1) = E_defo
C     /' "light" fragment '/
      Edefo(4,1,IZ1) = E_defo
      E_defo = Lymass(Z2,A2,beta2_opt) - Lymass(Z2,A2,0.0)
      Edefo(0,2,IZ2) = E_defo
C     /' "heavy" fragment '/
      Edefo(4,2,IZ2) = E_defo
      END DO
C     '
C     /' Mode 1: deformed shells (light) and spherical (heavy) '/
      DO I = 10 , I_Z_CN - 10
      Z1 = I
      Z2 = I_Z_CN - Z1
      A1 = (Z1 - 0.5E0) / REAL(I_Z_CN) * REAL(I_A_CN)
C     /' polarization roughly considered '/
      A2 = I_A_CN - A1
      IF (  I_Z_CN * 0.5 .LT. ZC_Mode_1  ) THEN
C     ' Beta_opt_light(A1,A2,Z1,Z2,dneck,0,rbeta_ld)
C     /' nu_mean of Cf requires shells in the light fragment: '/
      rbeta = beta_light(I,betaL0,betaL1) - 0.1
C     ' smaller than general deformation of light fragment
C     '        (less neck influence due to spherical heavy fragment)
      IF (  rbeta .LT. 0  )  rbeta = 0
      Else
      rbeta = beta_heavy(I,betaH0,betaH1)
C     ' equal to S2 channel
      IF (  rbeta .LT. 0  )  rbeta = 0
      End If
      Beta(1,1,I) = rbeta
C     /' "light" fragment '/
      E_defo = Lymass(Z1,A1,rbeta) - Lymass(Z1,A1,0.0)
      Edefo(1,1,I) = E_defo
C     /' "light" fragment '/
      END DO
C     '
      DO I = 10 , I_Z_CN - 10
      rbeta = 0
      Beta(1,2,I) = rbeta
      Edefo(1,2,I) = 0
C     /' "heavy" fragment (at S1 shell) '/
      END DO
C     '
C     /' Mode 2: deformed shells (light and heavy) '/
      DO I = 10 , I_Z_CN - 10
      Z1 = I
      Z2 = I_Z_CN - Z1
      A1 = (Z1 - 0.5E0) / REAL(I_Z_CN) * REAL(I_A_CN)
C     /' polarization roughly considered '/
      A2 = I_A_CN - A1
      IF (  I_Z_CN * 0.5 .LT. ZC_Mode_2  ) THEN
C     ' Beta_opt_light(A1,A2,Z1,Z2,dneck,beta_heavy(Z2),rbeta_ld)
      rbeta = beta_light(I,betaL0,betaL1)
C     ' general deformation of light fragment
      IF (  rbeta .LT. 0  )  rbeta = 0
C     ' negative values replaced by 0
      Else
      rbeta = beta_heavy(I,betaH0,betaH1)
C     ' equal to S2 channel
      End If
      Beta(2,1,I) = rbeta
      E_defo = Lymass(Z1,A1,rbeta) - Lymass(Z1,A1,0.0)
      Edefo(2,1,I) = E_defo
      END DO
      DO I = 10 , I_Z_CN - 10
      rbeta = beta_heavy(I,betaH0,betaH1)
C     /' "heavy" fragment (at S2 shell)'/
      IF (  rbeta .LT. 0  )  rbeta = 0
C     ' negative values replaced by 0
      Beta(2,2,I) = rbeta
      Z1 = I
      A1 = (Z1 + 0.5E0) / I_Z_CN * I_A_CN
C     /' polarization roughly considered '/
      E_defo = Lymass(Z1,A1,rbeta) - Lymass(Z1,A1,0.0)
      Edefo(2,2,I) = E_defo
      END DO
C     '
C     /' Mode 3 '/
      DO I = 10 , I_Z_CN - 10
      Z1 = I
      Z2 = I_Z_CN - Z1
      A1 = (Z1 - 0.5E0) / REAL(I_Z_CN) * REAL(I_A_CN)
C     /' polarization roughly considered '/
      A2 = I_A_CN - A1
      rbeta = beta_light(I,betaL0,betaL1)
      rbeta = Max(rbeta-0.10,0.0)
C     /' for low nu-bar of lightest fragments '/
C     '  Beta_opt_light(A1,A2,Z1,Z2,dneck,beta_heavy(Z2,betaH0,betaH1),rbeta)
      Beta(3,1,I) = rbeta
      E_defo = Lymass(Z1,A1,rbeta) - Lymass(Z1,A1,0.0)
      Edefo(3,1,I) = E_defo
      END DO
      DO I = 10 , I_Z_CN - 10
      rbeta = beta_heavy(I,betaH0,betaH1) + 0.2
C     /' for high nu-bar of heaviest fragments '/
      IF (  rbeta .LT. 0  )  rbeta = 0
      Beta(3,2,I) = rbeta
      Z1 = I
      A1 = (Z1 + 0.5E0) / REAL(I_Z_CN) * REAL(I_A_CN)
C     /' polarization roughly considered '/
      E_defo = Lymass(Z1,A1,rbeta) - Lymass(Z1,A1,0.0)
      Edefo(3,2,I) = E_defo
      END DO
C     '
C     /' Mode 5: (Channel ST1 in both fragments) '/
      DO I = 10 , I_Z_CN - 10
      Z1 = I
      Z2 = I_Z_CN - Z1
      rbeta = Beta(1,2,I)
      IF (  rbeta .LT. 0  )  rbeta = 0
      Beta(5,1,Int(Z1)) = rbeta
      Beta(5,2,Int(Z1)) = rbeta
      END DO
C     '
C     /' Mode 6: (Channel ST2 in both fragments) '/
      DO I = 10 , I_Z_CN - 10
      Z1 = I
      Z2 = I_Z_CN - Z1
      rbeta = Beta(2,2,I)
      IF (  rbeta .LT. 0  )  rbeta = 0
      Beta(6,1,Int(Z1)) = rbeta
      Beta(6,2,Int(Z1)) = rbeta
      END DO
C     '
C     '
C     /' Mean Z as a function of mass '/
C     '
C     /' Mode 0 '/
      DO I = 10 , I_A_CN - 10
      ZUCD = REAL(I) / REAL(I_A_CN) * REAL(I_Z_CN)
      beta1 = Beta(0,1,Int(ZUCD + 0.5))
      beta2 = Beta(0,2,Int(I_Z_CN - ZUCD + 0.5))
      Z1 = Z_equi(I_Z_CN,I,I_A_CN - I,beta1,beta2,dneck,0,0.0,1.0)
      Zmean(0,1,I) = Z1
      Zshift(0,1,I) = Z1 - ZUCD
      Zmean(0,2,I_A_CN - I) = I_Z_CN - Z1
      Zshift(0,2,I_A_CN - I) = ZUCD - Z1
      END DO
C     '
C     /' Mode 1 '/
      DO I = 10 , I_A_CN - 10
      ZUCD = REAL(I) / REAL(I_A_CN) * REAL(I_Z_CN)
      Z = ZUCD + ZPOL1
C     /' Charge polarisation is considered in a crude way '/
      beta1 = Beta(1,1,NINT(Z))
C     /' "light" fragment '/
      Z = ZUCD - ZPOL1
      beta2 = Beta(1,2,NINT(I_Z_CN-Z))
C     /' "heavy" fragment  at S1 shell '/
      IF (  REAL(I_Z_CN) * 0.5 .LT. ZC_Mode_1  ) THEN
      Z1 = Z_equi(I_Z_CN,I,I_A_CN - I,beta1,beta2,dneck,1,POLARadd,
     *POLARfac)
      Else
      Z1 = Z_equi(I_Z_CN,I,I_A_CN - I,beta1,beta2,dneck,1,0.0,0.0)
      End If
      Z1 = Z1 + ZPOL1
C     /' Charge polarization by shell '/
C     '
      IF (  I_Z_CN - Z1 .LT. 50 .AND. (I_Z_CN - Z1) .GT. Z1  ) THEN
      Z1 = I_Z_CN - 50
C     /' Z of mean heavy fragment not below 50 '/
      END IF
C     '
      Zmean(1,1,I) = Z1
      Zshift(1,1,I) = Z1 - ZUCD
C     ' neutron-deficient
      Zmean(1,2,I_A_CN - I) = I_Z_CN - Z1
      Zshift(1,2,I_A_CN - I) = ZUCD - Z1
C     ' neutron rich at shell
      END DO
C     '
C     /' Mode 2 '/
      DO I = 10 , I_A_CN - 10
      ZUCD = REAL(I) / REAL(I_A_CN) * REAL(I_Z_CN)
      Z = ZUCD
C     /' Charge polarisation is here neglected '/
      beta1 = Beta(2,1,NINT(Z))
      beta2 = Beta(2,2,NINT(I_Z_CN-Z))
      IF (  REAL(I_Z_CN) * 0.5 .LT. ZC_Mode_2  ) THEN
      Z1 = Z_equi(I_Z_CN,I,I_A_CN-I,beta1,beta2,dneck,2,POLARadd,
     *POLARfac)
      Else
      Z1 = Z_equi(I_Z_CN,I,I_A_CN-I,beta1,beta2,dneck,2,0.0,0.0)
      End If
C     '
      Zmean(2,1,I) = Z1
      Zshift(2,1,I) = Z1 - ZUCD
C     ' neutron deficieint
      Zmean(2,2,I_A_CN - I) = I_Z_CN - Z1
      Zshift(2,2,I_A_CN - I) = ZUCD - Z1
C     ' neutron rich at shell
      END DO
C     '
C     /' Mode 3 '/
      DO I = 10 , I_A_CN - 10
      ZUCD = REAL(I) / REAL(I_A_CN) * REAL(I_Z_CN)
      Z = ZUCD
C     /' Charge polarisation is here neglected '/
      beta1 = Beta(3,1,NINT(Z))
      beta2 = Beta(3,2,NINT(I_Z_CN-Z))
      Z1 = Z_equi(I_Z_CN,I,I_A_CN - I,beta1,beta2,dneck,3,POLARadd,
     *POLARfac)
      Zmean(3,1,I) = Z1
      Zshift(3,1,I) = Z1 - ZUCD
      Zmean(3,2,I_A_CN - I) = I_Z_CN - Z1
      Zshift(3,2,I_A_CN - I) = ZUCD - Z1
      END DO
C     '
C     /' Mode 4 (assumed to be equal to mode 0) '/
      DO I = 10 , I_A_CN - 10
      Zmean(4,1,I) = Zmean(0,1,I)
      Zshift(4,1,I) = Zshift(0,1,I)
      Zmean(4,2,I_A_CN - I) = Zmean(0,2,I_A_CN - I)
      Zshift(4,2,I_A_CN - I) = Zshift(0,2,I_A_CN - I)
      END DO
C     '
C     '
C     /' General relations between Z and A of fission channels '/
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_0 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(0,2,NINT(RA))
      END DO
      AC_Mode_0 = (ZC_Mode_0 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
C     /' mean position in mass '/
      NC_Mode_0 = AC_Mode_0 - ZC_Mode_0
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_1 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(1,2,NINT(RA))
      END DO
      AC_Mode_1 = (ZC_Mode_1 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_1 = AC_Mode_1 - ZC_Mode_1
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_2 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(2,2,NINT(RA))
      END DO
      AC_Mode_2 = (ZC_Mode_2 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_2 = AC_Mode_2 - ZC_Mode_2
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_3 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(3,2,NINT(RA))
      END DO
      AC_Mode_3 = (ZC_Mode_3 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_3 = AC_Mode_3 - ZC_Mode_3
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_4 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(4,2,ABS(NINT(RA)))
      END DO
      AC_Mode_4 = (ZC_Mode_4 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_4 = AC_Mode_4 - ZC_Mode_4
C     '
C     '
C     /' Potential curvatures of fission modes '/
C     '
C     ' For the width of the mass distribution (potential between saddle and scission):
C     ' Print Spin_pre_fission,P_I_rms_CN
      R_Z_Curv_S0 = 8.E0 / REAL(I_Z_CN)**2 * Masscurv(REAL(I_Z_CN),
     *REAL(I_A_CN),Spin_pre_fission,kappa)
C     ' For the yields of the fission channels (potential near saddle):
      R_Z_Curv1_S0 = 8.E0 / REAL(I_Z_CN)**2 * Masscurv1(REAL(I_Z_CN),
     *REAL(I_A_CN),0.0,kappa)
      R_A_Curv1_S0 = 8.E0 / REAL(I_A_CN)**2 * Masscurv1(REAL(I_Z_CN),
     *REAL(I_A_CN),0.0,kappa)
C     '
C     '
C     /' Energy transformation '/
C     '
      Select CASE( Emode)
      CASE( 0)
C     ' Energy above outer barrier given
      R_E_exc_Eb = R_E_exc_used
      R_E_exc_GS = R_E_exc_used + BFTFB(REAL(I_Z_CN),REAL(I_A_CN),1)
      CASE( 1,3,-1)
C     ' Energy above ground state given
      R_E_exc_Eb = R_E_exc_used - BFTFB(REAL(I_Z_CN),REAL(I_A_CN),1)
      R_E_exc_GS = R_E_exc_used
      CASE( 2)
C     ' kinetic energy of neutron given (SN = neutron separation energy)
C     '    SN = (U_Mass(Csng(I_Z_CN),Csng(I_A_CN-1)) + Lypair(I_Z_CN,I_A_CN-1))     '       -(U_Mass(Csng(I_Z_CN),Csng(I_A_CN)) + Lypair(I_Z_CN,I_A_CN))
C     '    R_E_exc_GS = R_E_exc_used + SN
      SN = AME2012(I_Z_CN,I_A_CN-1) - AME2012(I_Z_CN,I_A_CN)
      R_E_exc_GS = R_E_exc_used * ((P_A_CN-1) / P_A_CN) + SN
C     '           target                      CN
      R_E_exc_Eb = R_E_exc_GS - BFTFB(REAL(I_Z_CN),REAL(I_A_CN),1)
      CASE( 12)
C     ' kinetic energy of proton given (Sprot = proton separation energy)
C     '    Sprot = (U_Mass(Csng(I_Z_CN-1),Csng(I_A_CN-1)) + Lypair(I_Z_CN-1,I_A_CN-1))     '       -(U_Mass(Csng(I_Z_CN),Csng(I_A_CN)) + Lypair(I_Z_CN,I_A_CN))
C     '    R_E_exc_GS = R_E_exc_used + Sprot
      Sprot = AME2012(I_Z_CN-1,I_A_CN-1) - AME2012(I_Z_CN,I_A_CN)
      R_E_exc_GS = R_E_exc_used * ((P_A_CN-1) / P_A_CN) + Sprot
      R_E_exc_Eb = R_E_exc_GS - BFTFB(REAL(I_Z_CN),REAL(I_A_CN),1)
      End Select
C     '
C     '
C     /' Fission barriers -> global parameters '/
C     '
      B_F = BFTF(REAL(I_Z_CN),REAL(I_A_CN),1)
      B_F_ld = BFTF(REAL(I_Z_CN),REAL(I_A_CN),0)
      E_B = BFTFB(REAL(I_Z_CN),REAL(I_A_CN),1)
      E_B_ld = BFTFB(REAL(I_Z_CN),REAL(I_A_CN),0)
C     '
C     '
C     /' Barriers and excitation energies of the fission modes '/
C     '
      E_exc_S0_prov = R_E_exc_Eb
C     '
C     '
C     /' Additional influence of N=82 assumed '/
      Delta_NZ_Pol = 82.E0/50.E0 - REAL(I_N_CN)/REAL(I_Z_CN)
      R_Shell_S1_eff = P_Shell_S1 * (1.E0 - P_Att_Pol *
     *Abs(Delta_NZ_Pol))
C     '
C     '
C     /' In Pu,the Z=50 shell meets Z=44 in the light fragment. '/
C     /' A deformed shell at Z=44 is assumed to explain the enhancement        of the S1 channel around Pu '/
C     /' This very same shell automatically produces the double-humped '/
C     /' mass distribution in 180Hg '/
      S1_enhance = P_Shell_SL4 + (REAL(I_Z_CN) - ZC_Mode_1 -
     *ZC_Mode_4L)**2 * P_Z_Curv_SL4
C     'Print "ZC_Mode_1,ZC_Mode_4",ZC_Mode_1,ZC_Mode_4
C     'Print "Delta-Z S1-S4,S1_enhance",I_Z_CN-ZC_Mode_1 - ZC_Mode_4L,S1_enhance
      IF (  S1_enhance .GT. 0  )  S1_enhance = 0
      R_Shell_S1_eff = R_Shell_S1_eff + S1_enhance
C     '
C     /' The high TKE of S1 in 242Pu(sf) (and neighbours) is obtained by assuming '/
C     /' that the Z=44 shell reduces the deformation of the light fragment. '/
      DO I = 10 , I_Z_CN - 10
      Z1 = I
      A1 = (Z1 - 0.5E0) / REAL(I_Z_CN) * REAL(I_A_CN)
C     /' polarization roughly considered '/
C     '      Beta(1,1,Z1) = Beta(1,1,Z1) + 0.15 * S1_enhance   /' "light" fragment '/
      Beta(1,1,I) = exp(S1_enhance) * Beta(1,1,I) +
     *(1.E0-exp(S1_enhance)) * (Beta(1,1,I)-0.25)
      Beta(1,1,I) = Max(Beta(1,1,I),0.0)
      E_defo = Lymass(Z1,A1,Beta(1,1,I)) - Lymass(Z1,A1,0.0)
      Edefo(1,1,I) = E_defo
C     /' "light" fragment '/
      END DO
C     '
C     ' Influence of S2 shell in complementary fragment
C     ' May be called "S12 fission channel"
      T_Asym_Mode_2 = 0.5
      SigZ_Mode_2 = SQRT(0.5E0 * T_Asym_Mode_2/(P_Z_Curv_S2))
      SigA_Mode_2 = SigZ_Mode_2 * REAL(I_A_CN) / REAL(I_Z_CN)
      S1_enhance = P_Shell_S2 * U_Box(REAL(P_A_CN) - AC_Mode_2 -
     *AC_Mode_1,SigA_Mode_2,P_A_Width_S2) *P_A_Width_S2
      IF (  S1_enhance .LT. 0.01  ) THEN
      R_Shell_S1_eff = R_Shell_S1_eff + S1_enhance
      End If
C     ' Modify deformation of complementary fragment in corresponding analyzer
C     '
C     ' Overlap of S2 and shell in light fragment
      R_Shell_S2_eff = P_Shell_S2
C     '   S2_enhance = P_Shell_S4 +  '             (Csng(I_Z_CN) - ZC_Mode_2 - ZC_Mode_4)^2 * P_Z_Curv_S4
C     '   If S2_enhance > 0 Then S2_enhance = 0
C     '   R_Shell_S2_eff = R_Shell_S2_eff + S2_enhance
C     '
C     ' Overlap of S3 and shell in light fragment
      R_Shell_S3_eff = P_Shell_S3 * (1.E0 - PZ_S3_olap_curv          *
     *(REAL(I_Z_CN) - ZC_Mode_3 - PZ_S3_olap_pos)**2)
C     '        * (Csng(I_Z_CN) - 60.5E0 - PZ_S3_olap_pos)^2)
      R_Shell_S3_eff = Min(R_Shell_S3_eff,0.0)
C     '
C     '   R_Shell_S4_eff = 2.0 * (P_Shell_S4 + P_Z_Curv_S4*(ZC_Mode_4 - ZC_Mode_0)^2)
      R_Shell_S4_eff = 2.0 * (P_Shell_S4 + P_Z_Curv_S4 * (ZC_Mode_4 -
     *ZC_Mode_0)**2)
C     ' overlap of S4 in both fragments
      IF (  R_Shell_S4_eff .GT. P_Shell_S4  )  R_Shell_S4_eff =
     *P_Shell_S4
C     ' no overlap at large distance
C     '
      E_ld_S1 = R_A_Curv1_S0 * (REAL(I_A_CN)/REAL(I_Z_CN)*(ZC_MODE_1 -
     *ZC_MODE_0) )**2
      B_S1 = E_ld_S1 + R_Shell_S1_eff
      E_exc_S1_prov = E_Exc_S0_prov - B_S1
C     '
      E_ld_S2 = R_A_Curv1_S0 * (REAL(I_A_CN)/REAL(I_Z_CN)*(ZC_MODE_2 -
     *ZC_MODE_0) )**2
      B_S2 = E_ld_S2 + R_Shell_S2_eff
      E_exc_S2_prov = E_Exc_S0_prov - B_S2
C     '
      E_ld_S3 = R_A_Curv1_S0 * (REAL(I_A_CN)/REAL(I_Z_CN)*(ZC_MODE_3 -
     *ZC_MODE_0) )**2
      B_S3 = E_ld_S3 + R_Shell_S3_eff
      E_exc_S3_prov = E_Exc_S0_prov - B_S3
C     '
      IF (  I_A_CN .LT. 220  ) THEN
C     ' Only here S4 is close enough to symmetry to have a chance
      E_ld_S4 = R_A_Curv1_S0 * (REAL(I_A_CN)/REAL(I_Z_CN)*(ZC_MODE_4 -
     *ZC_MODE_0) )**2
      B_S4 = E_ld_S4 + R_Shell_S4_eff
      E_exc_S4_prov = E_Exc_S0_prov - B_S4
      Else
      B_S4 = 9999
      E_exc_S4_prov = - 9999
      End If
C     '
C     /' Mode 11 (overlap of channel 1 in light and heavy fragment '/
C     /' Potential depth with respect to liquid-drop potential: B_S11 '/
      B_S11 = 2.E0 * (R_Shell_S1_eff + De_Defo_S1              +
     *P_Z_Curv_S1 * (ZC_Mode_1 - ZC_Mode_0)**2 ) - De_Defo_S1
C     '
C     '
C     /' Lowering of effective barrier by lower ZPM due to larger width in
C     partial overlap region (shells in light and heavy fragment) '/
      DES11ZPM = Level_S11 * Min(Abs(ZC_Mode_1 - ZC_Mode_0),
     *4.E0*P_Z_Curv_S1)
C     ' Print B_S11,DES11ZPM,ZC_Mode_1-ZC_Mode_0
C     '
      B_S11 = B_S11 + DES11ZPM
C     '
C     '  If B_S11 > R_Shell_S1_eff + 0.5E0 Then
C     '   If B_S11 > R_Shell_S1_eff + Level_S11 Then
C     '     B_S11 = 100   ' S1 and S11 are exclusive
C     '   Else
C     '     B_S11 = Min(B_S11,R_Shell_S1_eff)
C     '   End If
C     '
C     '
      E_exc_S11_prov = E_Exc_S0_prov - B_S11
C     '
C     /' Mode 22 (overlap of channel 2 in light and heavy fragment '/
C     /' Potential depth with respect to liquid-drop potential: B_S22 '/
C     '
C     '   B_S22 = 2.E0 * (E_ld_S2 + P_Shell_S2)  '       + 2.E0 * P_Z_Curv_S2 * (ZC_Mode_2 - ZC_Mode_0)^2   /' Parabola '/
C     'Print E_ld_S2,P_Shell_S2,P_Z_Curv_S2,ZC_Mode_2,ZC_Mode_0
      B_S22 = 2.E0 * R_Shell_S2_eff * U_Box(REAL(P_A_CN)/2.0 -
     *AC_Mode_2,SigA_Mode_2,P_A_Width_S2) * P_A_Width_S2
C     ' The integral of U_Box is normalized,not the height!
C     '    If Abs((P_A_CN/2.E0) - AC_Mode_2) > P_A_Width_S2 Then B_S22 = 9999
      IF (  P_A_CN .LT. 226  )  B_S22 = 9999
C     '
      E_exc_S22_prov = E_Exc_S0_prov - B_S22
C     '
C     '
      E_Min_Barr = Min(0.0,B_S1)
      E_Min_Barr = Min(E_Min_Barr,B_S2)
      E_Min_Barr = Min(E_Min_Barr,B_S3)
      E_Min_Barr = Min(E_Min_Barr,B_S4)
      E_Min_Barr = Min(E_Min_Barr,B_S11)
      E_Min_Barr = Min(E_Min_Barr,B_S22)
C     '
C     /' Energy minus the height of the respective fission saddle '/
      E_exc_S0 = E_exc_S0_prov + E_Min_Barr - Delta_S0
      E_exc_S1 = E_exc_S1_prov + E_Min_Barr
      E_exc_S2 = E_exc_S2_prov + E_Min_Barr
      E_exc_S3 = E_exc_S3_prov + E_Min_Barr
      E_exc_S4 = E_exc_S4_prov + E_Min_Barr
      E_exc_S11 = E_exc_S11_prov + E_Min_Barr
      E_exc_S22 = E_exc_S22_prov + E_Min_Barr
C     '
C     /' Energy above the lowest fission saddle '/
      E_exc_Barr = Max(E_Exc_S0,E_Exc_S1)
      E_exc_Barr = Max(E_exc_Barr,E_Exc_S2)
      E_exc_Barr = Max(E_exc_Barr,E_Exc_S3)
      E_exc_Barr = Max(E_exc_Barr,E_Exc_S4)
      E_exc_Barr = Max(E_exc_Barr,E_exc_S11)
      E_exc_Barr = Max(E_exc_Barr,E_exc_S22)
C     '
C     '
C     /' Collective temperature used for calculating the widths
C     in mass asymmetry and charge polarization '/
C     '
      IF (  E_Exc_S0 .LT. 0  ) THEN
      E_tunn = -E_Exc_S0
      Else
      E_tunn = 0
      END IF
      R_E_exc_eff = Max(0.1,E_Exc_S0)
C     '  T_Coll_Mode_0 = TFCOLL * R_E_exc_eff + _  /' empirical,replaced by TRusanov '/
      T_Coll_Mode_0 = TCOLLFRAC * (De_Saddle_Scission(REAL(I_Z_CN)**2 /
     *           REAL(I_A_CN)**0.33333E0,ESHIFTSASCI_coll) - E_tunn)
      T_Coll_Mode_0 = Max(T_Coll_Mode_0,0.0)
C     '
C     ' Print "De_SS,E_tunn,T_Coll ";De_Saddle_Scission(I_Z_CN^2/I_A_CN^0.3333,ESHIFTSASCI_coll),E_tunn,T_Coll_Mode_0
C     '
C     /' Temperature description fitting to the empirical systematics of Rusanov et al. '/
C     /' Here from Ye. N. Gruzintsev et al.,Z. Phys. A 323 (1986) 307 '/
C     /' Empirical description of the nuclear temperature according to the '/
C     /' Fermi-gas description. Should be valid at higher excitation energies '/
      T_Rusanov = TRusanov(R_E_exc_eff,REAL(I_A_CN))
C     '  Print "Temperatures,(GEF,Total,Rusanov): ",T_Coll_Mode_0,TFCOLL * R_E_exc_eff,T_Rusanov
      T_Coll_Mode_0 = Max(T_Coll_Mode_0,T_Rusanov)
C     /' Transition vom const. temp. to Fermi gas occurs around 20 MeV by MAX function '/
C     '    T_Pol_Mode_0 = T_Pol_Red * T_Coll_Mode_0
C     '
C     ' Application of the statistical model,intrinsic temperature at saddle
      T_Pol_Mode_0 = U_Temp(0.5 * REAL(I_Z_CN),0.5 *REAL(I_A_CN),
     *R_E_exc_eff,0,0,Tscale,Econd)
      T_Asym_Mode_0 = SQRT(T_Coll_Mode_0**2 + (6E0*TCOLLMIN)**2)
C     '
      E_pot_scission = (De_Saddle_Scission(REAL(I_Z_CN)**2 /
     *    REAL(I_A_CN)**0.33333E0,ESHIFTSASCI_intr) - E_tunn)
C     '
C     /' Suppression of S1 fission channel due to reduced pairing in 132Sn '/
C     /' At very low excitation energy on the fission path,the binding energy at the
C     S1 fission channel does not profit as much from pairing as SL and S2,
C     because pairing is reduced in magic nuclei. This leads to a reduction of
C     the yield in S1 in the case that the fully paired ground-state configuration
C     is populated on the fission path with a considerable probability. '/
C     '   EeffS2 = Max(E_exc_S2,0.0) + EDISSFRAC * E_pot_scission - 2.3E0
C     '   EeffS2 = Max(0.0,EeffS2)
C     /' -2.3 MeV,because fission channels are assumed to be chosen before scission '/
C     '
C     '   If EeffS2 < ETHRESHSUPPS1 + 2.E0 * ESIGSUPPS1 Then
C     '     E_exc_S1 = E_exc_S1 -  '        0.5E0 * 4.E0 * 12.E0 / Sqr(132.E0) * Gaussintegral(ETHRESHSUPPS1 - EeffS2,ESIGSUPPS1)
C     '   EndIf
C     '
      T_low_S1_used = T_low_S1
C     '
      T_Coll_Mode_1 = TFCOLL * Max(E_exc_S1,0.E0) + TCOLLFRAC *
     *(De_Saddle_Scission(I_Z_CN**2 / I_A_CN**0.33333E0,
     *ESHIFTSASCI_coll) - E_tunn)
      T_Coll_Mode_1 = Max(T_Coll_mode_1,0.0)
C     '    T_Pol_Mode_1 = T_Pol_Red * T_Coll_Mode_1
      T_Pol_Mode_1 = T_Pol_Mode_0
      T_Asym_Mode_1 = SQRT(T_Coll_Mode_1**2 + (4.0*TCOLLMIN)**2)
C     ' TCOLLMIN for ZPM
C     '
      T_Coll_Mode_2 = TFCOLL * Max(E_exc_S2,0.E0) + TCOLLFRAC *
     *(De_Saddle_Scission(REAL(I_Z_CN)**2 /
     *REAL(I_A_CN)**0.33333E0,ESHIFTSASCI_coll) - E_tunn)
      T_Coll_Mode_2 = Max(T_Coll_mode_2,0.0)
C     '    T_Pol_Mode_2 = T_Pol_Red * T_Coll_Mode_2
      T_Pol_Mode_2 = T_Pol_Mode_0
      T_Asym_Mode_2 = SQRT(T_Coll_Mode_2**2 + TCOLLMIN**2)
C     '
      T_Coll_Mode_3 = TFCOLL * Max(E_exc_S3,0.E0) + TCOLLFRAC *
     *(De_Saddle_Scission(REAL(I_Z_CN)**2 /
     *REAL(I_A_CN)**0.33333E0,ESHIFTSASCI_coll) - E_tunn)
      T_Coll_Mode_3 = Max(T_Coll_mode_3,0.0)
C     '    T_Pol_Mode_3 = T_Pol_Red * T_Coll_Mode_3
      T_Pol_Mode_3 = T_Pol_Mode_0
      T_Asym_Mode_3 = SQRT(T_Coll_Mode_3**2 + TCOLLMIN**2)
C     '
      T_Coll_Mode_4 = TFCOLL * Max(E_exc_S4,0.E0) + TCOLLFRAC *
     *(De_Saddle_Scission(REAL(I_Z_CN)**2 /
     *REAL(I_A_CN)**0.33333E0,ESHIFTSASCI_coll) - E_tunn)
      T_Coll_Mode_4 = Max(T_Coll_mode_4,0.0)
C     '    T_Pol_Mode_4 = T_Pol_Red * T_Coll_Mode_4
      T_Pol_Mode_4 = T_Pol_Mode_0
      T_Asym_Mode_4 = SQRT(T_Coll_Mode_4**2 + 4.0*TCOLLMIN**2)
C     ' ZPM like S1
C     '
C     /' Stiffness in polarization '/
C     '
      RZ = REAL(I_Z_CN) * 0.5E0
      RA = REAL(I_A_CN) * 0.5E0
      beta1 = Beta(0,1,NINT(RZ))
      beta2 = Beta(0,2,NINT(RZ))
      R_Pol_Curv_S0 = ( LyMass( RZ - 1.E0,RA,beta1 ) +
     *LyMass( RZ + 1.0E0,RA,beta2 ) +              LyMass( RZ + 1.0E0,
     *RA,beta1 ) +              LyMass( RZ - 1.0E0,RA,beta2 ) +
     *     ecoul( RZ - 1.0E0,RA,beta1,RZ + 1.0E0,RA,beta2,dneck) +
     *        ecoul( RZ + 1.0E0,RA,beta1,RZ - 1.0E0,RA,beta2,dneck) -
     *       2.0E0*ecoul( RZ,RA,beta1,RZ,RA,beta2,dneck) -
     *2.0E0*LyMass( RZ,RA,beta1 ) -          2.0E0*LyMass( RZ,RA,beta2)
     *) * 0.5E0
C     '
      P_Pol_Curv_S0 = R_Pol_Curv_S0
C     '
      R_Pol_Curv_S1 = R_Pol_Curv_S0
      R_Pol_Curv_S2 = R_Pol_Curv_S0
      R_Pol_Curv_S3 = R_Pol_Curv_S0
      R_Pol_Curv_S4 = R_Pol_Curv_S0
C     '
C     '
C     '
C     /' Mean values and standard deviations of fission modes '/
C     '
      SIGZ_Mode_0 = SQRT(0.5E0 * T_Asym_Mode_0/R_Z_Curv_S0)
      IF (  T_Pol_Mode_0 .GT. 1.E-2  ) THEN
      SigPol_Mode_0 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S0 /
     *          Tanh(HOMPOL/(2.E0 * T_Pol_Mode_0)))
      Else
      SigPol_Mode_0 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S0)
C     /' including influence of zero-point motion '/
      END IF
C     '
C     '
      R_E_intr_S1 = Max(E_Exc_S1+Lypair(I_Z_CN,I_A_CN),0.0)
      R_Att(1) = exp(-R_E_intr_S1/Shell_fading)
      R_Att(5) = R_Att(1)
      R_Att_Sad(1) = exp(-R_E_intr_S1/Shell_fading)
      R_Att_Sad(5) = R_Att_Sad(1)
      SIGZ_Mode_1 = SQRT(0.5E0 *
     *T_Asym_Mode_1/(P_Z_Curv_S1*SQRT(R_Att(1))))
      IF (  T_Pol_Mode_1 .GT. 1.E-2  ) THEN
      SigPol_Mode_1 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S1 /
     *          Tanh(HOMPOL/(2.E0 * T_Pol_Mode_1)))
      Else
      SigPol_Mode_1 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S1)
      END IF
C     '
      R_E_intr_S2 = Max(E_Exc_S2+Lypair(I_Z_CN,I_A_CN),0.0)
      R_Att(2) = exp(-R_E_intr_S2/Shell_fading)
      R_Att(6) = R_Att(2)
      R_Att_Sad(2) = exp(-R_E_intr_S2/Shell_fading)
      R_Att_Sad(6) = R_Att_Sad(2)
      SIGZ_Mode_2 = SQRT(0.5E0 *
     *T_Asym_Mode_2/(P_Z_Curv_S2*SQRT(R_Att(2))))
      IF (  T_Pol_Mode_2 .GT. 1.E-2  ) THEN
      SigPol_Mode_2 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S2 /
     *          Tanh(HOMPOL/(2.E0 * T_Pol_Mode_2)))
      Else
      SigPol_Mode_2 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S2)
      End If
C     '
      R_E_intr_S3 = Max(E_exc_S3+Lypair(I_Z_CN,I_A_CN),0.0)
      R_Att(3) = exp(-R_E_intr_S3/Shell_fading)
      R_Att_Sad(3) = exp(-R_E_intr_S3/Shell_fading)
      SIGZ_Mode_3 = SQRT(0.5E0 *
     *T_Asym_Mode_3/(P_Z_Curv_S3*SQRT(R_Att(3))))
      IF (  T_Pol_Mode_3 .GT. 1.E-2  ) THEN
      SigPol_Mode_3 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S3 /
     *          Tanh(HOMPOL/(2.E0 * T_Pol_Mode_3)))
      Else
      SigPol_Mode_3 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S3)
      End if
C     '
      R_E_intr_S4 = Max(E_exc_S4+Lypair(I_Z_CN,I_A_CN),0.0)
      R_Att(4) = exp(-R_E_intr_S4/Shell_fading)
      R_Att_Sad(4) = exp(-R_E_intr_S4/Shell_fading)
      SIGZ_Mode_4 = SQRT(0.5E0 *
     *T_Asym_Mode_4/(P_Z_Curv_S4*SQRT(R_Att(4))))
      IF (  T_Pol_Mode_4 .GT. 1.E-2  ) THEN
      SigPol_Mode_4 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S4 /
     *          Tanh(HOMPOL/(2.E0 * T_Pol_Mode_4)))
      Else
      SigPol_Mode_4 = SQRT(0.25E0 * HOMPOL / R_Pol_Curv_S4)
      End if
C     '
C     '
C     '
C     /' Energy-dependent shift of fission channels '/
C     '    Scope
      P_Z_Curv_S1_eff = P_Z_Curv_S1 * P_Z_Curvmod_S1
      P_Z_Curv_S2_eff = P_Z_Curv_S2 * P_Z_Curvmod_S2
      P_Z_Curv_S3_eff = P_Z_Curv_S3 * P_Z_Curvmod_S3
      P_Z_Curv_S4_eff = P_Z_Curv_S4 * P_Z_Curvmod_S4
C     '
      DZ_S1 = ZC_Mode_1 * (P_Z_Curv_S1_eff*R_Att(1) / (R_Z_Curv_S0 +
     *P_Z_Curv_S1_eff*R_Att(1))             - (P_Z_Curv_S1_eff /
     *(R_Z_Curv_S0 + P_Z_Curv_S1_eff) ) )
      DZ_S2 = ZC_Mode_2 * (P_Z_Curv_S2_eff*R_Att(2) / (R_Z_Curv_S0 +
     *P_Z_Curv_S2_eff*R_Att(2))              - (P_Z_Curv_S2_eff /
     *(R_Z_Curv_S0 + P_Z_Curv_S2_eff) ) )
      DZ_S3 = ZC_Mode_3 * (P_Z_Curv_S3_eff*R_Att(3) / (R_Z_Curv_S0 +
     *P_Z_Curv_S3_eff*R_Att(3))              - (P_Z_Curv_S3_eff /
     *(R_Z_Curv_S0 + P_Z_Curv_S3_eff) ) )
      DZ_S4 = SIGN(1.0,ZC_Mode_4 - ZC_Mode_0) * ZC_Mode_4 *
     *(P_Z_Curv_S4_eff*R_Att(4) / (R_Z_Curv_S0 +
     *P_Z_Curv_S4_eff*R_Att(4))              - (P_Z_Curv_S4_eff /
     *(R_Z_Curv_S0 + P_Z_Curv_S4_eff) ) )
C     '
C     ' Empirical shift of S2 channel at low excitation energy at scission
C     ' for better reproduction of 238U(s,f) and some data for Th isotopes.
C     ' Does not solve the problem of 229Th(nth,f).
      EtotS2 = Max(E_Exc_S2 + EDISSFRAC * E_pot_scission,0.0)
      IF (  EtotS2 .LT. 5.E0  ) THEN
      DZ_S2 = DZ_S2 + (5.E0 - EtotS2) * 0.1
      End If
C     '
C     '   DZ_S1 = 0
C     '   DZ_S2 = 0
C     '   DZ_S3 = 0
C     '   DZ_S4 = 0
C     '
C     '
      P_Z_Mean_S0 = ZC_Mode_0
      ZC_Mode_1 = ZC_Mode_1 + DZ_S1
      P_Z_Mean_S1 = ZC_Mode_1
C     /' Copy to global parameter '/
      ZC_Mode_2 = ZC_Mode_2 + DZ_S2
      P_Z_Mean_S2 = ZC_Mode_2
C     /'             "            '/
      ZC_Mode_3 = ZC_Mode_3 + DZ_S3
      P_Z_Mean_S3 = ZC_Mode_3
C     '   ZC_Mode_4 = ZC_Mode_4 + DZ_S4
C     ' shift is very small,because S4 exists only close to symmetry
      P_Z_Mean_S4 = ZC_Mode_4
C     '    End Scope
C     '
C     /' Energy dependence of charge polarization '/
C     /' Due to washing out of shells '/
C     '
      DO I = 10 , I_A_CN - 10
C     ' mass number
      DO J = 1 , 4
C     ' fission channel
      DO K = 1 , 2
C     ' light - heavy group
      Zshift(J,K,I) = Zshift(0,K,I) + (Zshift(J,K,I) - Zshift(0,K,
     *I))*R_Att(J)
      END DO
      END DO
      END DO
C     '
C     '
C     /' Energy dependence of shell-induced deformation '/
C     /' Due to washing out of shells '/
C     /' (Under development) '/
C     /'For I = 10 To I_Z_CN - 10  ' mass number
C     For J = 1 To 4           ' fission channel
C     For K = 1 To 2         ' light - heavy group
C     beta(J,K,I) = beta(0,K,I) + (beta(J,K,I) - beta(0,K,I))*R_Att_Sad(J)
C     if beta(J,K,I) < 0 Then
C     beta(J,K,I) = 0
C     End If
C     Z1 = I
C     Z2 = I_Z_CN - Z1
C     A1 = Z1 / Csng(I_Z_CN) * Csng(I_A_CN)
C     A2 = I_A_CN - A1
C     E_defo = Lymass(Z1,A1,beta(J,K,I)) - Lymass(Z1,A1,0.0)
C     Edefo(J,K,I) = E_defo
C     Next
C     Next
C     Next  '/
C     '
C     '
C     '
C     '
C     /' General relations between Z and A of fission channels '/
C     /' 2nd iteration '/
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_0 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(0,2,NINT(RA))
      END DO
      AC_Mode_0 = (ZC_Mode_0 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
C     /' mean position in mass '/
      NC_Mode_0 = AC_Mode_0 - ZC_Mode_0
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_1 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(1,2,NINT(RA))
      END DO
      AC_Mode_1 = (ZC_Mode_1 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_1 = AC_Mode_1 - ZC_Mode_1
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_2 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(2,2,NINT(RA))
      END DO
      AC_Mode_2 = (ZC_Mode_2 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_2 = AC_Mode_2 - ZC_Mode_2
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_3 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(3,2,NINT(RA))
      END DO
      AC_Mode_3 = (ZC_Mode_3 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_3 = AC_Mode_3 - ZC_Mode_3
C     '
      RZpol = 0
      DO I = 1 , 3
      RA = (ZC_Mode_4 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      RZpol = Zshift(4,2,ABS(NINT(RA)))
      END DO
      AC_Mode_4 = (ZC_Mode_4 - RZPol) * REAL(I_A_CN) / REAL(I_Z_CN)
      NC_Mode_4 = AC_Mode_4 - ZC_Mode_4
C     '
C     '
C     '
C     /' Yields of the fission modes '/
C     '
      Yield_Mode_0 = Getyield(E_exc_S0,E_exc_S0,T_low_SL,
     *TEgidy(REAL(I_A_CN),0.E0,Tscale))
C     '
      Yield_Mode_1 = Getyield(E_exc_S1,E_exc_S0,T_low_S1_used,
     *TEgidy(REAL(I_A_CN),R_Shell_S1_eff + dE_Defo_S1,Tscale))
C     /'  - Getyield(E_exc_S0 - E_ld_S1,T_low,T_high) '/
C     '
      Yield_Mode_2 = Getyield(E_exc_S2,E_exc_S0,T_low_S2,
     *TEgidy(REAL(I_A_CN),R_Shell_S2_eff + dE_Defo_S2,Tscale))
C     /'  - Getyield(E_exc_S0 - E_ld_S2,T_low,T_high) '/
C     '
      Yield_Mode_3 = Getyield(E_exc_S3,E_exc_S0,T_low_S3,
     *TEgidy(REAL(I_A_CN),R_Shell_S3_eff + dE_Defo_S3,Tscale))
C     /'  - Getyield(E_exc_S0 - E_ld_S3,T_low,T_high) '/
C     '
      Yield_Mode_4 = Getyield(E_exc_S4,E_exc_S0,T_low_S4,
     *TEgidy(REAL(I_A_CN),R_Shell_S4_eff + dE_Defo_S4,Tscale))
C     /'   - Getyield(E_exc_S0 - E_ld_S4,T_low,T_high) '/
C     '
C     'Print TEgidy(Csng(I_A_CN),0.E0,Tscale),TEgidy(Csng(I_A_CN),R_Shell_S2_eff + dE_Defo_S2,Tscale),de_Defo_S2
C     'sleep
C     '
      IF (  B_S11 .GT. B_S1  ) THEN
      Yield_Mode_11 = 0.0
      Else
      Yield_Mode_11 = Getyield(E_exc_S11,E_exc_S0,T_low_S11,
     *TEgidy(REAL(I_A_CN),R_Shell_S1_eff + 2.E0 * dE_Defo_S1,Tscale))
      End If
C     '
      IF (  B_S22 .GT. B_S2  ) THEN
      Yield_Mode_22 = 0.0
      Else
      Yield_Mode_22 = Getyield(E_exc_S22,E_exc_S0,T_low_S2,
     *TEgidy(REAL(I_A_CN),R_Shell_S2_eff,Tscale))
      End If
C     '
C     '
      Yield_Norm = Yield_Mode_0 + Yield_Mode_1 + Yield_Mode_2 +
     *Yield_Mode_3 + Yield_Mode_4 + Yield_Mode_11 + Yield_Mode_22
      Yield_Mode_0 = Yield_Mode_0 / Yield_Norm
      Yield_Mode_1 = Yield_Mode_1 / Yield_Norm
      Yield_Mode_2 = Yield_Mode_2 / Yield_Norm
      Yield_Mode_3 = Yield_Mode_3 / Yield_Norm
      Yield_Mode_4 = Yield_Mode_4 / Yield_Norm
      Yield_Mode_11 = Yield_Mode_11 / Yield_Norm
      Yield_Mode_22 = Yield_Mode_22 / Yield_Norm
C     '
C     '
C     /' Mass widhts of the fission channels '/
C     '
      SigA_Mode_0 = SigZ_Mode_0 * REAL(I_A_CN) / REAL(I_Z_CN)
C     /' width in mass '/
      SigA_Mode_1 = SigZ_Mode_1 * REAL(I_A_CN) / REAL(I_Z_CN)
      SigA_Mode_1 = Min(SigA_Mode_1,SigA_Mode_0)
C     ' not broader than liquid-drop
      SigA_Mode_2 = SigZ_Mode_2 * REAL(I_A_CN) / REAL(I_Z_CN)
      SigA_Mode_2 = Min(SigA_Mode_2,SigA_Mode_0)
C     ' not broader than liquid-drop
      SigA_Mode_3 = SigZ_Mode_3 * REAL(I_A_CN) / REAL(I_Z_CN)
      SigA_Mode_3 = Min(SigA_Mode_3,SigA_Mode_0)
      SigA_Mode_4 = SigZ_mode_4 * REAL(I_A_CN) / REAL(I_Z_CN)
      SigA_Mode_4 = Min(SigA_Mode_4,SigA_Mode_0)
      SigA_Mode_11 = SigZ_Mode_1 * SQRT(2.E0) * REAL(I_A_CN) /
     *REAL(I_Z_CN)
      SigA_Mode_11 = Min(SigA_Mode_11,SigA_Mode_0)
      SigA_Mode_22 = SigZ_Mode_2 * SQRT(2.E0) * REAL(I_A_CN) /
     *REAL(I_Z_CN)
      SigA_Mode_22 = Min(SigA_Mode_22,SigA_Mode_0)
C     '
C     '
C     '
C     /' Shell effects of different fission channels '/
C     /' This is the "real" microscopic shell effect,not the effective shell-correction energy '/
C     /' EShell acts on the level density and determines the T parameter '/
C     '
      DO I = 1 , I_A_CN - 1
      DO J = 0 , 4
      EShell(J,1,I) = 0
C     /' Shells in "light" fragment assumed to be zero '/
      END DO
      DU0 = 0
      EShell(0,2,I) = 0
C     /' Shell = 0 in symmetric mode '/
      DU1 = R_Shell_S1_eff + dE_Defo_S1
C     /' + R_A_Curv1_S1 * (AC_Mode_1 - Float(I,6))**2; '/
      DU1 = MIN(DU1,0.E0)
C     /' Technical limit '/
      EShell(1,2,I) = DU1
C     '
      DU2 = R_Shell_S2_eff + dE_Defo_S2
C     /' + R_A_Curv1_S2 * (AC_Mode_2 - Float(I,6))**2; '/
      DU2 = Min(DU2,0.E0)
C     /' Technical limit '/
      EShell(2,2,I) = DU2
C     '
      DU3 = R_Shell_S3_eff + dE_Defo_S3
C     /' + R_A_Curv1_S3 * (AC_Mode_3 - Float(I,6))**2; '/
      DU3 = Min(DU3,0.E0)
C     /' Technical limit '/
      EShell(3,2,I) = DU3
C     '
      DU4 = R_Shell_S4_eff + dE_Defo_S4
C     /' + R_A_Curv1_S4 * (AC_Mode_4 - Float(I,6))**2; '/
      DU4 = Min(DU4,0.E0)
C     /' Technical limit '/
      EShell(4,2,I) = DU4
C     '
      END DO
C     '
C     '
C     /' Intrinsic temperatures of fragments at scission '/
C     '
C     /' Mean values '/
      T_intr_Mode_0 = TEgidy(AC_Mode_0,0.0,0.8)
      T_intr_Mode_1_heavy = TEgidy(AC_Mode_1,R_Shell_S1_eff +
     *dE_Defo_S1,Tscale)
      T_intr_Mode_1_light = TEgidy(REAL(I_A_CN) - AC_Mode_1,0.0,Tscale)
      T_intr_Mode_2_heavy = TEgidy(AC_Mode_2,R_Shell_S2_eff +
     *dE_Defo_S2,Tscale)
      T_intr_Mode_2_light = TEgidy(REAL(I_A_CN) - AC_Mode_2,0.0,Tscale)
      T_intr_Mode_3_heavy = TEgidy(AC_Mode_3,R_Shell_S3_eff +
     *dE_Defo_S3,Tscale)
      T_intr_Mode_3_light = TEgidy(REAL(I_A_CN) - AC_Mode_3,0.0,Tscale)
      T_intr_Mode_4_heavy = TEgidy(AC_Mode_4,R_Shell_S4_eff +
     *dE_Defo_S4,Tscale)
      T_intr_Mode_4_light = TEgidy(REAL(I_A_CN) - AC_Mode_4,0.0,Tscale)
C     '
C     '
C     /' Mass-dependent values of individual fragments '/
C     /' Mode 0 '/
      DO I = 1 , I_A_CN - 1
      T = TEgidy(REAL(I),EShell(0,1,I),Tscale)
      Temp(0,1,I) = T
C     /' "light" fragment at freeze-out (somewhere before scission) '/
      T = TEgidy(REAL(I),EShell(0,2,I),Tscale)
      Temp(0,2,I) = T
C     /' "heavy" fragment at freeze-out (somewhere before scission) '/
C     '
      T = TEgidy(REAL(I),0.0,1.0)
      TempFF(0,1,I) = T
C     ' FF in their ground state
      TempFF(0,2,I) = T
C     ' FF in their ground state
      END DO
C     '
C     /' Mode 1 '/
      DO I = 1 , I_A_CN - 1
      T = TEgidy(REAL(I),EShell(1,1,I),Tscale)
      Temp(1,1,I) = T
C     /' "light" fragment '/
      T = TEgidy(REAL(I),EShell(1,2,I),Tscale)
      Temp(1,2,I) = T
C     /' "heavy" fragment '/
C     '
      T = TEgidy(REAL(I),0.0,1.0)
      TempFF(1,1,I) = T
C     ' FF in their ground state
      TempFF(1,2,I) = T
C     ' FF in their ground state
      END DO
C     '
C     /' Mode 2 '/
      DO I = 1 , I_A_CN - 1
      T = TEgidy(REAL(I),EShell(2,1,I),Tscale)
      Temp(2,1,I) = T
C     /' "light" fragment '/
      T = TEgidy(REAL(I),EShell(2,2,I),Tscale)
      Temp(2,2,I) = T
C     /' "heavy" fragment '/
C     '
C     /' The next section is introduced,because energy sorting is not strong enough,
C     when shells are only introduced in the heavy fragment.
C     Ad hoc assumption: For Mode 2 there are shells in both fragments of about
C     equal size. Technically,we neglect the shells in both fragments.
C     This has about the same effect for the energy sorting. '/
      T = TEgidy(REAL(I),0.0,Tscale)
C     ' FF at scssion
      Temp(2,1,I) = T
C     /' "light" fragment '/
      T = TEgidy(REAL(I),0.0,Tscale)
C     ' FF at scission
      Temp(2,2,I) = T
C     /' "heavy" fragment '/
C     '
      T = TEgidy(REAL(I),0.0,1.0)
C     ' shell effect neglected
      TempFF(2,1,I) = T
C     ' FFs in their ground state
      TempFF(2,2,I) = T
C     ' FFs in their ground state
      END DO
C     '
C     /' Mode 3 '/
      DO I = 1 , I_A_CN -1
      T = TEgidy(REAL(I),0.0,Tscale)
      Temp(3,1,I) = T
      T = TEgidy(REAL(I),0.0,Tscale)
      Temp(3,2,I) = T
C     '
      T = TEgidy(REAL(I),0.0,1.0)
      TempFF(3,1,I) = T
C     ' FF in their ground state
      TempFF(3,2,I) = T
C     ' FF in their ground state
      END DO
C     '
C     /' Mode 4 '/
      DO I = 1 , I_A_CN -1
      T = TEgidy(REAL(I),0.0,Tscale)
      Temp(4,1,I) = T
      T = TEgidy(REAL(I),0.0,Tscale)
      Temp(4,2,I) = T
C     '
      T = TEgidy(REAL(I),0.0,1.0)
      TempFF(4,1,I) = T
C     ' FF in their ground state
      TempFF(4,2,I) = T
C     ' FF in their ground state
      END DO
C     '
C     '
C     /'** Intrinsic excitation energy at saddle and at scission as well as   **'/
C     /'** Even-odd effect in proton and neutron number for each fission mode **'/
      DO I_Mode = 0 , 6
      E_coll_saddle(I_Mode) = 0
      IF (  I_Mode .EQ. 0  )  Etot = E_exc_S0
      IF (  I_Mode .EQ. 1  )  Etot = E_exc_S1
      IF (  I_Mode .EQ. 2  )  Etot = E_exc_S2
      IF (  I_Mode .EQ. 3  )  Etot = E_exc_S3
      IF (  I_Mode .EQ. 4  )  Etot = E_exc_S4
      IF (  I_Mode .EQ. 5  )  Etot = E_exc_S11
      IF (  I_Mode .EQ. 6  )  Etot = E_exc_S22
C     '
      IF (   MOD(I_Z_CN,2)  +  MOD(I_N_CN,2)  .EQ. 0  ) THEN
C     /' Even-even CN '/
      IF (  Etot .GT. 0 .AND. Etot .LT. 2.E0 * 14.E0/SQRT(REAL(I_A_CN))
     * ) THEN
      E_coll_saddle(I_Mode) = Etot
      Etot = 0
C     /' Excitation below the pairing gap in even-even CN goes into collective excitations '/
      End If
      End If
C     '
C     '    If I_Z_CN Mod 2 + I_N_CN Mod 2 = 0 Then    ' even-even
C     '      Ediff = Min(Etot,14.0/sqr(Csng(I_A_CN)))
C     '    End If
C     '    If I_Z_CN Mod 2 + I_N_CN Mod 2 = 1 Then    ' even-odd or odd-even
C     '       Ediff = Min(Etot,2.0 * 14.0/sqr(Csng(I_A_CN)))
C     '    End If
C     '    Ediff = Max(Ediff,0.0)
C     '    Etot = Etot - Ediff
C     '
C     '
      IF (  Etot .LT. 0  ) THEN
      E_tunn = -Etot
      Else
      E_tunn = 0
      END IF
      Etot = Max(Etot,0.0)
C     '
      E_pot_scission = (De_Saddle_Scission(REAL(I_Z_CN)**2 /
     *    REAL(I_A_CN)**0.33333E0,ESHIFTSASCI_intr) )
      Etot = Etot + EDISSFRAC * (E_pot_scission - E_tunn)
C     /' All excitation energy at saddle and part of the potential-energy gain to scission
C     go into intrinsic excitation energy at scission '/
C     '
C     '
C     '
C     '
      IF (  I_Mode .EQ. 2  ) THEN
      EINTR_SCISSION = Etot
C     /' (For Mode 2) Global parameter '/
      End If
C     '
      DO IA1 = 40 , I_A_CN - 40
C     '
      IA2 = I_A_CN - IA1
      IF (  I_Mode .LE. 4  ) THEN
      T1 = Temp(I_Mode,1,IA1)
      T2 = Temp(I_Mode,2,IA2)
      End If
      IF (  I_Mode .EQ. 5  ) THEN
      T1 = Temp(1,2,IA1)
      T2 = Temp(1,2,IA2)
      End If
      IF (  I_Mode .EQ. 6  ) THEN
      T1 = Temp(2,2,IA1)
      T2 = Temp(2,2,IA2)
      End If
      DT = ABS(T2 - T1)
C     '
C     /' Even-odd effect '/
      IF (   MOD(I_Z_CN,2)  .EQ. 0  ) THEN
CAK   Rincr1P = Exp(-Etot/PZ_EO_symm)
      Rincr1P = Exp(min(-Etot/PZ_EO_symm,80.))
      Else
      Rincr1P = 0
      End If
      IF (   MOD(I_N_CN,2)  .EQ. 0  ) THEN
CAK   Rincr1N = Exp(-Etot/PN_EO_symm)
      Rincr1N = Exp(min(-Etot/PN_EO_symm,80.))
      Else
      Rincr1N = 0
      End If
      PEOZ(I_Mode,1,IA1) = Rincr1P
      PEOZ(I_Mode,2,IA2) = Rincr1P
      PEON(I_Mode,1,IA1) = Rincr1N
      PEON(I_Mode,2,IA2) = Rincr1N
C     '
      Rincr2 = Gaussintegral(DT/Etot-R_EO_Thresh,
     *R_EO_Sigma*(DT+0.0001))
C     /' even-odd effect due to asymmetry '/
      Rincr2P = (R_EO_MAX - Rincr1P) * Rincr2
      Rincr2N = (R_EO_MAX - Rincr1N) * Rincr2
C     '
      IF (  IA1 .LT. IA2  ) THEN
C     ' A1 is lighter
      PEOZ(I_Mode,1,IA1) = PEOZ(I_Mode,1,IA1) + Rincr2P
      IF (   MOD(I_Z_CN,2)  .EQ. 0  ) THEN
      PEOZ(I_Mode,2,IA2) = PEOZ(I_Mode,2,IA2) + Rincr2P
      Else
      PEOZ(I_Mode,2,IA2) = PEOZ(I_Mode,2,IA2) - Rincr2P
      End if
      PEON(I_Mode,1,IA1) = PEON(I_Mode,1,IA1) + Rincr2N
      IF (   MOD(I_N_CN,2)  .EQ. 0  ) THEN
      PEON(I_Mode,2,IA2) = PEON(I_Mode,2,IA2) + Rincr2N
      Else
      PEON(I_Mode,2,IA2) = PEON(I_Mode,2,IA2) - Rincr2N
      End if
      Else
      PEOZ(I_Mode,1,IA1) = PEOZ(I_Mode,2,IA1)
      PEON(I_Mode,1,IA1) = PEON(I_Mode,2,IA1)
      PEOZ(I_Mode,2,IA2) = PEOZ(I_Mode,1,IA2)
      PEON(I_Mode,2,IA2) = PEON(I_Mode,1,IA2)
      End If
C     '
C     '
C     /'  Else
C     PEOZ(I_Mode,2,IA2) =                PEOZ(I_Mode,1,IA2) + Rincr2P
C     IF I_Z_CN Mod 2 = 0 Then
C     PEOZ(I_Mode,1,IA1) =                 PEOZ(I_Mode,1,IA1) + Rincr2P
C     Else
C     PEOZ(I_Mode,1,IA1) =                 PEOZ(I_Mode,1,IA1) - Rincr2P
C     End if
C     PEON(I_Mode,2,IA2) =              PEON(I_Mode,2,IA2) + Rincr2N
C     IF I_N_CN Mod 1 = 0 Then
C     PEON(I_Mode,1,IA1) =                 PEON(I_Mode,1,IA1) + Rincr2N
C     Else
C     PEON(I_Mode,1,IA1) =                 PEON(I_Mode,1,IA1) - Rincr2N
C     End if
C     End If  '/
C     '
      PEOZ(I_Mode,1,IA1) = PEOZ(I_Mode,1,IA1) * EOscale
      PEOZ(I_Mode,2,IA2) = PEOZ(I_Mode,2,IA2) * EOscale
      PEON(I_Mode,1,IA1) = PEON(I_Mode,1,IA1) * EOscale
      PEON(I_Mode,2,IA2) = PEON(I_Mode,2,IA2) * EOscale
C     '
C     /' Energy sorting '/
C     /' E1 = Etot * Gaussintegral(T2-T1,0.03); '/
      IF (  Abs(T1-T2) .LT. 1.E-6  ) THEN
      E1 = 0.5E0 * Etot
      Else
      E1ES = Csort * T1 * T2 / ( Abs(T1 - T2) )
      E1ES = Min(E1ES,0.5E0*Etot)
C     /' Asymptotic value after "complete" energy sorting '/
      E1FG = Etot * IA1 / I_A_CN
C     /' in Fermi-gas regime '/
      IF (  Etot .LT. 13  )  E1 = E1ES
C     ' complete energy sorting
      IF (  Etot .GE. 13 .AND. Etot .LE. 20  ) THEN
C     ' transition region
      E1 = E1ES + (Etot-13)/7*(E1FG-E1ES)
      End If
      IF (  Etot .GT. 20  )  E1 = E1FG
C     ' Fermi-gas regime
      End If
      E2 = Etot - E1
      EPART(I_Mode,1,IA1) = Max(E1,0.0)
C     /' Mean E* in light fragment '/
      EPART(I_Mode,2,IA2) = Max(E2,0.0)
C     /' Mean E* in heavy fragment '/
      END DO
      END DO
C     '
C     '
C     /'** RMS angular momentum of fission fragments **'/
C     /' Following Naik et al.,EPJ A 31 (2007) 195 and  '/
C     /' S. G. Kadmensky,Phys. At. Nucl. 71 (2008) 1193 '/
C     '
C     '   Scope
      Spin_CN = P_J_CN
      P_I_rms_CN = P_J_CN
      Spin_pre_fission = SPIN_CN
C     ' target or CN ground-state spin
C     '
      DO IZ1 = 10 , I_Z_CN - 10
      AUCD = Int(REAL(IZ1) * REAL(I_A_CN) / REAL(I_Z_CN))
      DO IA1 = Int(AUCD - 15) , Int(AUCD + 15)
      IN1 = IA1 - IZ1
      IF (  IA1 - IZ1 .GE. 10  ) THEN
C     /' Rigid momentum of inertia for spherical nucleus '/
      I_rigid_spher = 1.16E0**2 * REAL(IA1)**1.6667E0 / 103.8415E0
C     /' unit: hbar^2/MeV '/
      DO I_Mode = 0 , 6
C     '
C     /' First (normally light) fission fragment: '/
C     '
      beta1 = Beta(I_Mode,1,IZ1)
      alph = beta1 / SQRT(4.E0 * pi / 5.E0)
      I_rigid = I_rigid_spher * (1.E0 + 0.5E0*alph + 9.E0/7.E0*alph**2)
C     /' From Hasse & Myers,Geometrical Relationships ... '/
      E_exc = EPART(I_Mode,1,IA1)
      IF (  E_exc .LT. 0  )  E_exc = 0
      T = U_Temp(REAL(IZ1),REAL(IA1),E_exc,1,1,Tscale,Econd)
C     '   T = sqr(T^2 + 0.8^2)       ' For ZPM
C     '   T = T_orbital
C     '   T =  sqr(T^2 + T_orbital^2)
      IF (  T_orbital .GT. 0.1  ) THEN
      T = T_orbital / tanh(T_orbital/T)
C     ' T_orbital represents the ZPM
      End If
      I_eff = I_rigid * (1.E0 - 0.8E0 * exp(-0.693E0 * E_exc / 5.E0))
      J_rms = SQRT(2.E0 * I_eff * T)
C     '
      J_rms = J_rms * Jscaling
C     '
      IF (   MOD(IZ1,2)  .EQ. 1 .OR.  MOD(IN1,2)  .EQ. 1  )  J_rms =
     *J_rms + Spin_odd * (REAL(IA1)/140.0)**0.66667
C     '                * Max(0,1 - (E_exc-1)/9) /' empirical '/
C     /' Additional angular momentum of unpaired proton. '/
C     /' See also Tomar et al.,Pramana 68 (2007) 111 '/
C     '
C     ' Print Z1,I_Mode,beta1,T,E_exc,Spin_CN
C     ' Print " ",I_rigid_spher,I_rigid,I_eff,J_rms
C     '
      J_rms = SQRT(J_rms**2 + (IA1/I_A_CN * Spin_pre_fission)**2)
C     '
      SpinRMSNZ(I_Mode,1,IA1-IZ1,IZ1) = J_rms
C     '
C     '     Print
C     '     Print IA1,T,E_exc,I_rigid_spher,I_rigid,I_eff,J_rms
C     '
C     /' Second (normally heavy) fission fragment: '/
C     '
      beta2 = Beta(I_Mode,2,IZ1)
      alph = beta2 / SQRT(4.E0 * pi / 5.E0)
      I_rigid = I_rigid_spher * (1.E0 + 0.5E0*alph + 9.E0/7.E0*alph**2)
C     /' From Hasse & Myers,Geometrical Relationships ... '/
      E_exc = EPART(I_Mode,2,IA1)
      IF (  E_exc .LT. 0  )  E_exc = 0
      T = U_Temp(REAL(IZ1),REAL(IA1),E_exc,1,1,Tscale,Econd)
C     '    T = sqr(T^2 + 0.8^2)       ' For ZPM
C     '    T = T_orbital
C     '    T =  sqr(T^2 + T_orbital^2)
      IF (  T_orbital .GT. 0.1  ) THEN
      T = T_orbital / tanh(T_orbital/T)
C     ' T_orbital represents the ZPM
      End If
      I_eff = I_rigid * (1.E0 - 0.8E0 * exp(-0.693E0 * E_exc / 5.E0))
      J_rms = SQRT(2.E0 * I_eff * T)
C     '
      J_rms = J_rms * Jscaling
C     '
      IF (   MOD(IZ1,2)  .EQ. 1 .OR.  MOD(IN1,2)  .EQ. 1  )  J_rms =
     *J_rms + Spin_odd * (REAL(IA1)/140.0)**0.66667
C     '                 * Max(0,1 - (E_exc-1)/9) /' empirical '/
C     /' Additional angular momentum of unpaired proton. '/
C     /' See also Tomar et al.,Pramana 68 (2007) 111 '/
C     '
      J_rms = SQRT(J_rms**2 + (IA1/I_A_CN * Spin_pre_fission)**2)
C     '
      SpinRMSNZ(I_Mode,2,IA1-IZ1,IZ1) = J_rms
C     '
C     '      Print IA1,T,E_exc,I_rigid_spher,I_rigid,I_eff,J_rms
C     '
      END DO
      ENd If
      END DO
      END DO
C     '   End Scope
C     '
C     ' ****************************************************************
C     ' *** Filling arrays with results in the folding mode (GEFSUB) ***
C     ' ****************************************************************
C     '
      DO I = 10 , I_A_CN - P_Z_CN - 10
      DO J = 10 , P_Z_CN - 10
      DO K = 0 , 6
      NZMPRE(K,I,J) = 0.0
      END DO
      END DO
      END DO
C     '
C     ' Mode 0
      DO I = 20 , I_A_CN - 20
      Ic = I_A_CN - I
      R_Help = Yield_Mode_0 * (U_Gauss_mod(AC_Mode_0 - REAL(I),
     *SigA_Mode_0)                  + U_Gauss_mod(AC_Mode_0 - REAL(Ic),
     *SigA_Mode_0))
C     ' Mass yield
      IF (  I .LT. Ic  ) THEN
      Zs = ZShift(0,1,I)
      Else
      Zs = -ZShift(0,1,Ic)
      End If
      DO J = 10 , P_Z_CN - 10
      Jc = P_Z_CN - J
      IF (  I-J .GE. 0 .AND. Ic-Jc .GE. 0 .AND. I-J .LE. 200 .AND.
     *Ic-Jc .LE. 200  ) THEN
      NZMPRE(0,I-J,J) = R_Help *
     *U_Gauss_mod(REAL(P_Z_CN)/REAL(I_A_CN)*REAL(I) + Zs - REAL(J),
     *SigPol_Mode_0) * U_Even_Odd(J,PEOZ(0,1,I)) * U_Even_Odd(I-J,
     *PEON(0,1,I))
      End If
      END DO
      END DO
C     '
C     ' Mode 1
      DO I = 20 , I_A_CN - 20
      Ic = I_A_CN - I
      R_Help = Yield_Mode_1 * (U_Gauss_mod(AC_Mode_1 - REAL(I),
     *SigA_Mode_1)                + U_Gauss_mod(AC_Mode_1 - REAL(Ic),
     *SigA_Mode_1))
C     ' Mass yield
      IF (  I .LT. Ic  ) THEN
      Zs = ZShift(1,1,I)
      Else
      Zs = -ZShift(1,1,Ic)
      End If
      DO J = 10 , P_Z_CN - 10
      Jc = P_Z_CN - J
      IF (  I-J .GE. 0 .AND. Ic-Jc .GE. 0 .AND. I-J .LE. 200 .AND.
     *Ic-Jc .LE. 200  ) THEN
      NZMPRE(1,I-J,J) = R_Help *
     *U_Gauss_mod(REAL(P_Z_CN)/REAL(I_A_CN)*REAL(I) + Zs - REAL(J),
     *SigPol_Mode_1)* U_Even_Odd(J,PEOZ(1,1,I)) * U_Even_Odd(I-J,PEON(1,
     *1,I))
      End If
      END DO
      END DO
C     '
C     ' Mode 2
      DO I = 20 , I_A_CN - 20
      Ic = I_A_CN - I
      R_Help = Yield_Mode_2 * (U_Box2(AC_Mode_2 - REAL(I),
     *SQRT(2.0)*S2leftmod*SigA_Mode_2,SQRT(2.0)*SigA_Mode_2,
     *P_A_Width_S2) +             U_Box2(AC_Mode_2 - REAL(Ic),
     *SQRT(2.0)*S2leftmod*SigA_Mode_2,SQRT(2.0)*SigA_Mode_2,
     *P_A_Width_S2))
      IF (  I .LT. Ic  ) THEN
      Zs = ZShift(2,1,I)
      Else
      Zs = -ZShift(2,1,Ic)
      End If
      DO J = 10 , P_Z_CN - 10
      Jc = P_Z_CN - J
      IF (  I-J .GE. 0 .AND. Ic-Jc .GE. 0 .AND. I-J .LE. 200 .AND.
     *Ic-Jc .LE. 200  ) THEN
      R_Cut1 = R_Help
      R_Cut2 = R_Help
      IF (  J .GT. Jc  ) THEN
      R_Cut1 = R_Help * Gaussintegral(REAL(J)-ZTRUNC50,
     *FTRUNC50*SigZ_Mode_2)
      Else
      R_Cut2 = R_Help * Gaussintegral(REAL(J)-ZTRUNC50,
     *FTRUNC50*SigZ_Mode_2)
      End If
      NZMPRE(2,I-J,J) = R_Help *
     *U_Gauss_mod(REAL(P_Z_CN)/REAL(I_A_CN)*REAL(I) + Zs - REAL(J),
     *SigPol_Mode_2) * U_Even_Odd(J,PEOZ(2,1,I)) * U_Even_Odd(I-J,
     *PEON(2,1,I))
      End If
      END DO
      END DO
C     '
C     ' Mode 3
      DO I = 20 , I_A_CN - 20
      Ic = I_A_CN - I
      R_Help = Yield_Mode_3 * (U_Gauss_mod(AC_Mode_3 - REAL(I),
     *SigA_Mode_3) +                     U_Gauss_mod(AC_Mode_3 -
     *REAL(Ic),SigA_Mode_3))
C     ' Mass yield
      IF (  I .LT. Ic  ) THEN
      Zs = ZShift(3,1,I)
      Else
      Zs = -ZShift(3,1,Ic)
      End If
      DO J = 10 , P_Z_CN - 10
      Jc = P_Z_CN - J
      IF (  I-J .GE. 0 .AND. Ic-Jc .GE. 0 .AND. I-J .LE. 200 .AND.
     *Ic-Jc .LE. 200  ) THEN
      NZMPRE(3,I-J,J) = R_Help *
     *U_Gauss_mod(REAL(P_Z_CN)/REAL(I_A_CN)*REAL(I) + Zs - REAL(J),
     *SigPol_Mode_3) * U_Even_Odd(J,PEOZ(3,1,I)) * U_Even_Odd(I-J,
     *PEON(3,1,I))
      End If
      END DO
      END DO
C     '
C     ' Mode 4
      DO I = 20 , I_A_CN - 20
      Ic = I_A_CN - I
      R_Help = Yield_Mode_4 * (U_Gauss_mod(AC_Mode_4 - REAL(I),
     *SigA_Mode_4) +                     U_Gauss_mod(AC_Mode_4 -
     *REAL(Ic),SigA_Mode_4))
C     ' Mass yield
      IF (  I .LT. Ic  ) THEN
      Zs = ZShift(3,1,I)
      Else
      Zs = -ZShift(3,1,Ic)
      End If
      DO J = 10 , P_Z_CN - 10
      Jc = P_Z_CN - J
      IF (  I-J .GE. 0 .AND. Ic-Jc .GE. 0 .AND. I-J .LE. 200 .AND.
     *Ic-Jc .LE. 200  ) THEN
      NZMPRE(4,I-J,J) = R_Help *
     *U_Gauss_mod(REAL(P_Z_CN)/REAL(I_A_CN)*REAL(I) + Zs - REAL(J),
     *SigPol_Mode_4) * U_Even_Odd(J,PEOZ(4,1,I)) * U_Even_Odd(I-J,
     *PEON(4,1,I))
      End If
      END DO
      END DO
C     '
C     ' Mode 11
      DO I = 20 , I_A_CN - 20
      Ic = I_A_CN - I
      R_Help = Yield_Mode_11 * (U_Gauss_mod(AC_Mode_0 - REAL(I),
     *SigA_Mode_11) +                     U_Gauss_mod(AC_Mode_0 -
     *REAL(Ic),SigA_Mode_11))
C     ' Mass yield
      DO J = 10 , P_Z_CN - 10
      Jc = P_Z_CN - J
      IF (  I-J .GE. 0 .AND. Ic-Jc .GE. 0 .AND. I-J .LE. 200 .AND.
     *Ic-Jc .LE. 200  ) THEN
      NZMPRE(5,I-J,J) = R_Help *
     *U_Gauss_mod(REAL(P_Z_CN)/REAL(I_A_CN)*REAL(I) - REAL(J),
     *SigPol_Mode_0) * U_Even_Odd(J,PEOZ(5,1,I)) * U_Even_Odd(I-J,
     *PEON(5,1,I))
      End If
      END DO
      END DO
C     '
C     ' Mode 22
      DO I = 20 , I_A_CN - 20
      Ic = I_A_CN - I
      R_Help = Yield_Mode_22 * (U_Gauss_mod(AC_Mode_0 - REAL(I),
     *SigA_Mode_22) +                     U_Gauss_mod(AC_Mode_0 -
     *REAL(Ic),SigA_Mode_22))
C     ' Mass yield
      DO J = 10 , P_Z_CN - 10
      Jc = P_Z_CN - J
      IF (  I-J .GE. 0 .AND. Ic-Jc .GE. 0 .AND. I-J .LE. 200 .AND.
     *Ic-Jc .LE. 200  ) THEN
      NZMPRE(6,I-J,J) = R_Help *
     *U_Gauss_mod(REAL(P_Z_CN)/REAL(I_A_CN)*REAL(I) - REAL(J),
     *SigPol_Mode_0) * U_Even_Odd(J,PEOZ(6,1,I)) * U_Even_Odd(I-J,
     *PEON(6,1,I))
      End If
      END DO
      END DO
C     '
C     '
C     ' Normalization
      R_Sum = 0
      DO I = 10 , (I_A_CN - P_Z_CN) - 10
      DO J = 10 , P_Z_CN - 10
      NZPRE(I,J) = 0
      DO K = 0 , 6
      IF (  NZMPRE(K,I,J) .GT. 0  ) THEN
      R_Sum = R_Sum + NZMPRE(K,I,J)
      NZPRE(I,J) = NZPRE(I,J) + NZMPRE(K,I,J)
C     ' sum of all modes
      End If
      END DO
      END DO
      END DO
C     ' Print R_Sum
      DO I = 10 , (I_A_CN - P_Z_CN) - 10
      DO J = 10 , P_Z_CN - 10
      NZPRE(I,J) = NZPRE(I,J) / R_Sum
      DO K = 0 , 6
      NZMPRE(K,I,J) = NZMPRE(K,I,J) / R_Sum
      END DO
      END DO
      END DO
C     '
C     ' Calculate and store distributions of fragment excitation energy and spin
C     '
      N_cases = 0
      DO N_index = 10 , (I_A_CN - P_Z_CN) - 10
C     ' Neutron number
      DO Z_index = 10 , P_Z_CN - 10
C     ' Atomic number
      DO M_index = 0 , 6
C     ' Fission channel
      IF (  NZMPRE(M_index,N_index,Z_index) .GT. Ymin  ) THEN
      N_cases = N_cases + 1
      IF (  N_cases .EQ. Ubound(NZMkey,1)  ) THEN
C     '           Print "Upper bound of NZkey reached"
C     '           Print "Result will be incomplete"
      End If
      NZMkey(N_cases,1) = M_index
C     ' Fission mode
      NZMkey(N_cases,2) = N_index
C     ' Neutron number of fragment
      NZMkey(N_cases,3) = Z_index
C     ' Atomic number of fragment
      End If
      END DO
      END DO
      END DO
C     ' Print "N_cases  ",N_cases
CAK   WRITE (*,*) "N_cases ",N_cases
C     '
      DO K = 1 , N_cases
      M_index = NZMkey(K,1)
C     ' fission mode
      N_index = NZMkey(K,2)
C     ' neutron number
      Z_index = NZMkey(K,3)
C     ' atomic number
      A_index = N_index + Z_index
C     '
C     ' Yield
      Ytab(K) = NZMpre(M_index,N_index,Z_index)
C     '
C     ' Angular momentum:
      DO I = 1 , 100
      IF (  M_index .LE. 4  ) THEN
      IF (  Z_index .LT. 0.5 * P_Z_CN  ) THEN
      Jtab(K,I) = U_LinGauss(REAL(I),SpinRMSNZ(M_index,1,N_index,
     *Z_index)/SQRT(2.0))
      Else
      Jtab(K,I) = U_LinGauss(REAL(I),SpinRMSNZ(M_index,2,N_index,
     *Z_index)/SQRT(2.0))
      End If
      End If
      IF (  M_index .EQ. 5  ) THEN
      Jtab(K,I) = U_LinGauss(REAL(I),SpinRMSNZ(1,2,N_index,
     *Z_index)/SQRT(2.0))
      End If
      IF (  M_index .EQ. 6  ) THEN
      Jtab(K,I) = U_LinGauss(REAL(I),SpinRMSNZ(2,2,N_index,
     *Z_index)/SQRT(2.0))
      End If
      END DO
C     '
C     ' Normalize numerically (due to non-continuous values)
C     '   Scope
      Rint = 0
      DO I = 1 , 100
      Rint = Rint + Jtab(K,I)
      END DO
      IF (  Rint .GT. 0  ) THEN
      DO I = 1 , 100
      Jtab(K,I) = Jtab(K,I) / Rint
      END DO
      End If
C     '   End Scope
C     '
C     '
C     ' Excitation energy:
C     ' 1. Deformation energy at scission
      IF (  M_index .EQ. 0  ) THEN
      IF (  Z_index .LT. 0.5 * P_Z_CN  ) THEN
      Eexc_mean = Edefo(M_index,1,Z_index)
      Eexc_sigma = ( Lymass(REAL(Z_index),REAL(A_index),beta(M_index,1,
     *Z_index) + SIGDEFO_0) -             Lymass(REAL(Z_index),
     *REAL(A_index),beta(M_index,1,Z_index) ))
      Else
      Eexc_mean = Edefo(M_index,2,Z_index)
      Eexc_sigma = ( Lymass(REAL(Z_index),REAL(A_index),beta(M_index,2,
     *Z_index) + SIGDEFO_0) -             Lymass(REAL(Z_index),
     *REAL(A_index),beta(M_index,2,Z_index) ))
      End If
      End If
      IF (  M_index .GT. 0 .AND. M_index .LE. 4  ) THEN
      IF (  Z_index .LT. 0.5 * P_Z_CN  ) THEN
      Eexc_mean = Edefo(M_index,1,Z_index)
      RS = SIGDEFO/SQRT(R_Att_Sad(M_index))
      Eexc_sigma = ( Lymass(REAL(Z_index),REAL(A_index),beta(M_index,1,
     *Z_index) + RS) -             Lymass(REAL(Z_index),REAL(A_index),
     *beta(M_index,1,Z_index) ))
      Else
      Eexc_mean = Edefo(M_index,2,Z_index)
      RS = SIGDEFO/SQRT(R_Att_Sad(M_index))
      Eexc_sigma = ( Lymass(REAL(Z_index),REAL(A_index),beta(M_index,2,
     *Z_index) + RS) -             Lymass(REAL(Z_index),REAL(A_index),
     *beta(M_index,2,Z_index) ))
      End If
      End If
      IF (  M_index .EQ. 5  ) THEN
      Eexc_mean = Edefo(1,2,Z_index)
      RS = SIGDEFO/SQRT(R_Att_Sad(M_index))
      Eexc_sigma = ( Lymass(REAL(Z_index),REAL(A_index),beta(1,2,
     *Z_index) + RS) -             Lymass(REAL(Z_index),REAL(A_index),
     *beta(1,2,Z_index) ))
      End If
      IF (  M_index .EQ. 6  ) THEN
      Eexc_mean = Edefo(2,2,Z_index)
      RS = SIGDEFO/SQRT(R_Att_Sad(M_index))
      Eexc_sigma = ( Lymass(REAL(Z_index),REAL(A_index),beta(2,2,
     *Z_index) + RS) -             Lymass(REAL(Z_index),REAL(A_index),
     *beta(2,2,Z_index) ))
      End If
      Eexc_mean = Max(Eexc_mean,0.0)
C     '
C     ' 2. Intrinsic excitation energy at scission
      IF (  Z_index .LT. 0.5 * REAL(P_Z_CN)  ) THEN
      Eexc_intr = EPART(M_index,1,A_index)
      Else
      Eexc_intr = EPART(M_index,2,A_index)
      End If
      IF (  M_index .EQ. 0  ) THEN
C     ' add shell and pairing of final fragment
      Eexc_intr = Eexc_intr - AME2012(Z_index,A_index) +
     *LDMass(REAL(Z_index),REAL(A_index),0.) - 2.0 * 12.0 /
     *SQRT(REAL(A_index))
C     ' general shift
      End If
      Eexc_intr = Max(Eexc_intr,0.0)
      Eexc_mean = Eexc_mean + Eexc_intr
      Eexc_sigma = SQRT(Eexc_sigma**2 + (EexcSIGrel * Eexc_intr)**2)
C     '
C     ' 3. Pairing staggering
      Eexc_mean = Eexc_mean - Lypair(Z_index,A_index)
C     '
C     ' 4. Collective energy
      Eexc_coll = 0.5 * ECOLLFRAC * (De_Saddle_Scission(REAL(P_Z_CN)**2
     */      REAL(I_A_CN)**0.33333E0,ESHIFTSASCI_coll) - E_tunn)
      Eexc_coll = Max(Eexc_coll,0.0)
      Eexc_sigma = SQRT(Eexc_sigma**2 + 0.5*(EexcSIGrel*Eexc_coll)**2)
      Eexc_mean = Eexc_mean + Eexc_coll + 0.5 * E_coll_saddle(M_index)
C     '
C     ' 5. Total excitation energy distribution of fragments (all contributions summed up)
C     ' This is the value above the yrast line. Erot must be added!
      DO I = 0 , 1000
C     ' 100 keV bins up to 100 MeV
      Etab(K,I) = exp(-(0.1*REAL(I)-Eexc_mean)**2/(2.0 * Eexc_sigma))
      END DO
CAK
      Emean(K)=Eexc_mean
      dEmean(K)=Eexc_sigma
CAK   print*,Z_index,A_index,Eexc_mean,Eexc_sigma,
CAK  +  Eexc_intr,Lypair(Z_index,A_index),Eexc_coll,
CAK  +  0.5 * E_coll_saddle(M_index),Ytab(K),M_index
CAK end
C     '
C     ' Normalize excitation-energy distribution
C     '   Scope
      RintE = 0
      DO I = 0 , 1000
      RintE = RintE + Etab(K,I)
      END DO
      IF (  RintE .GT. 0  ) THEN
      DO I = 0 , 1000
      Etab(K,I) = Etab(K,I) / RintE
      END DO
      End If
C     '   End Scope
C     '
      END DO
C     '
C     '
C     '
CAK
      return
      End
