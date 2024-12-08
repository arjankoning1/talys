subroutine checkkeyword
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in keywords
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2023-12-08: Latest version
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! Variables for reading input lines
!   inline           ! input line
!   nlines           ! number of input lines
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numkey=409        ! number of keywords
  integer            :: i                 ! counter
  integer            :: j                 ! counter
  character(len=132) :: key               ! keyword
  character(len=132) :: keyword(numkey)   ! keyword
  character(len=132) :: word(40)          ! words on input line
!
! Although it is difficult to prevent the user from all possible input errors, we can check for the use of wrong keywords
! and for unphysical values for most of the input variables.
!
! *********************** Check for wrong keywords *********************
!
! TALYS will stop if a keyword is incorrect
!
  data (keyword(i), i = 1, numkey) / ' ', 'a', 'aadjust', 'abundance', 'adddiscrete', 'addelastic', &
 &  'adepthcor', 'alimit', 'alphald', 'alphaomp', 'anfit', 'angles', 'anglescont', 'anglesrec', 'aradialcor', 'area', 'astro', &
 &  'astroe', 'astroex', 'astrogs', 'astrot', 'asys', 'autorot', 'avadjust', 'avadjustf', 'avdadjust', 'avdadjustf', &
 &  'avsoadjust', 'avsoadjustf', 'awadjust', 'awadjustf', 'awdadjust', 'awdadjustf', 'awsoadjust', 'awsoadjustf', 'axtype', &
 &  'bdamp', 'bdampadjust', 'best', 'bestbranch', 'bestend', 'bestpath', 'beta2', 'betafiscor', 'betafiscoradjust', 'betald', &
 &  'bins', 'block', &
 &  'branch', 'breakupmodel', 'cbarrier', 'cbreak', 'cfermi', 'cglobal', 'channelenergy', 'channels', 'cknock', 'class2', &
 &  'class2file', 'class2width', 'cnubar1', 'cnubar2', 'colenhance', 'colldamp', 'components', 'compound', 'core', &
 &  'coulomb', 'cpang', 'cstrip', 'ctable', 'ctableadjust', 'ctmglobal', 'd0', 'd1adjust', 'd2adjust', 'd3adjust', 'ddxmode', &
 &  'deformfile', 'deltaw', 'densfile',  &
 &  'deuteronomp', 'disctable', 'dispersion', 'dnfit', 'e0', 'e0adjust', 'e1file', 'eback', 'ebeam', 'eciscalc', 'eciscompound', &
 &  'ecisdwba', 'ecissave', 'egr', 'egradjust', 'ejectiles', 'ejoin', 'electronconv', 'element', 'elow', 'elwidth', &
 &  'emsdmin', 'endf', 'endfdetail', 'endfecis', 'energy', 'epr', 'epradjust', 'equidistant', 'equispec', &
 &  'estop', 'esurf', 'etable', &
 &  'etableadjust', 'exmatch', 'exmatchadjust', 'expmass', 'ffevaporation', 'ffmodel', 'ffspin', 'fileangle', 'filechannels', &
 &  'fileddxa', 'fileddxe', 'filedensity', 'filediscrete', 'fileelastic', 'filefission', 'filegamdis', 'filepsf', 'filerecoil', &
 &  'fileresidual', 'filespectrum', 'filetotal', 'fisbar', 'fisbaradjust', 'fisfeed', 'fishw', 'fishwadjust', 'fismodel', &
 &  'fismodelalt', 'fiso', 'fisom', 'fispartdamp', 'fission', 'fit', 'format', 'fsadjust', 'ftable', 'ftableadjust', 'fullhf', & 
 &  'fymodel', 'g', 'gadjust', 'gamgam', 'gamgamadjust', 'gamgamfit', 'gammald', 'gammashell1', &
 &  'gammashell2', 'gammax', 'gefran', 'ggr', 'ggradjust', 'giantresonance', 'gn', 'gnadjust', 'gnfit', 'gnorm', 'gp', 'gpadjust', &
 &  'gpr', 'gpradjust', 'group', 'gshell', 'hbstate', 'hbtransfile', 'ibeam', 'incadjust', &
 &  'inccalc', 'integral', 'isomer', 'jlmmode', 'jlmomp', 'kph', 'krotconstant', 'kvibmodel', 'labddx', &
 &  'ldmodel', 'ldmodelcn', 'ldmodelracap', 'levelfile', 'liso', 'localomp', 'ltarget', 'lurr', 'lv1adjust', 'lvadjust', &
 &  'lvsoadjust', 'lw1adjust', 'lwadjust', 'lwsoadjust', 'm1file', 'm2constant', 'm2limit', 'm2shift', &
 &  'macsfit', 'mass', 'massdir', &
 &  'massdis', 'massexcess', 'massmodel', 'massnucleus', 'maxband', 'maxchannel', 'maxenrec', 'maxlevelsbin', 'maxlevelsres', &
 &  'maxlevelstar', 'maxn', 'maxnrp', 'maxrot', 'maxz', 'maxzrp', 'micro', 'mpreeqmode', 'msdbins', 'multipreeq', 'nafit', &
 &  'ndfit', 'nffit', 'ngfit', 'nlevels', 'nlow', 'nnfit', 'nonthermlev', 'ntop', 'nulldev', 'ompenergyfile', 'omponly', &
 &  'onestep', 'optmod', 'optmodall', 'optmodfilen', 'optmodfilep', 'outall', 'outangle', 'outbasic', 'outbinspectra', 'outcheck', &
 &  'outdecay', 'outdensity', 'outdirect', 'outdiscrete', 'outdwba', 'outecis', 'outexcitation', &
 &  'outfission', 'outfy', 'outgamdis', 'outgamma', 'outinverse', 'outkd', 'outlegendre', 'outlevels', 'outmain', 'outomp', &
 &  'outpopulation', 'outpreequilibrium', 'outspectra', 'outtransenergy', 'pair', &
 &  'pairconstant', 'pairmodel', 'parity', 'partable', 'pfnsmodel', 'pglobal', 'phmodel', 'pnfit', 'popeps', 'popmev', &
 &  'preeqcomplex', 'preeqmode', 'preeqspin', 'preeqsurface', 'preequilibrium', 'production', 'projectile', 'pruitt', 'pruittset', &
 &  'pseudoresonances', 'psfglobal', 'pshift', 'pshiftadjust', 'pshiftconstant', 'ptable', 'ptableadjust', 'racap', & 
 &  'radialfile', 'radialmodel', 'radiounit', 'rcadjust', 'rclass2mom', 'reaction', 'recoil', 'recoilaverage', 'relativistic', &
 &  'rescuefile', 'reslib', 'resonance', 'rfiseps', 'rgamma', 'rho', 'riplomp', 'riplrisk', 'risomer', 'rnunu', 'rnupi', &
 &  'rotational', 'rpevap', 'rpinu', 'rpipi', 'rprime', 'rspincut', 'rspincutff', 'rspincutpreeq', 'rtransmom', 'rvadjust', &
 &  'rvadjustf', 'rvdadjust', &
 &  'rvdadjustf', 'rvsoadjust', 'rvsoadjustf', 'rwadjust', 'rwadjustf', 'rwdadjust', 'rwdadjustf', 'rwsoadjust', 'rwsoadjustf', &
 &  's2adjust', 'sacs', 'segment', 'sfexp', 'sfth', 'sgr', 'sgradjust', 'shellmodel', 'skipcn', 'soswitch', 'soukho', 'source', &
 &  'spherical', 'spincutmodel', 'spr', 'spradjust', 'statepot', 'strength', 'strengthm1', 'strucpath', 'sysreaction', &
 &  't', 'tadjust', 'tcool', 'tirrad', 'tjadjust', 'tmadjust', 'transeps', 'transpower', 'tres', 'twocomponent', 'ufermi', &
 &  'upbend', 'upbendc', 'upbendcadjust', 'upbende','upbendeadjust', 'upbendf', 'upbendfadjust', 'urr', 'urrnjoy', 'user', & 
 &  'v1adjust', 'v2adjust', 'v3adjust', 'v4adjust', 'vfiscor', 'vfiscoradjust', 'vinfadjust', 'vso1adjust', 'vso2adjust', &
 &  'w1adjust', 'w2adjust', 'w3adjust', 'w4adjust', 'wfcfactor', 'wso1adjust', 'wso2adjust', 'widthfluc', 'widthmode', &
 &  'wtable', 'wtableadjust', 'xsalphatherm', 'xscaptherm', 'xseps', 'xsptherm', 'yieldfile', 'yieldunit'/
!
! A keyword can be de-activated by putting a # in front of it.
! All first words of the input lines are checked against the list of keywords.
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified.
!
Loop1:  do i = 1, nlines
    call getkeywords(inline(i), word)
    key = word(1)
    if (key(1:1) == '#') cycle
    do j = 1, numkey
      if (keyword(j) == key) cycle Loop1
    enddo
    write(*, '(/" TALYS-error: Wrong keyword: ", a20)') key
    stop
  enddo Loop1
  return
end subroutine checkkeyword
! Copyright A.J. Koning 2021
