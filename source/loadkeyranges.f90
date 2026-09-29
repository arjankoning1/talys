subroutine loadkeyranges
  ! 
  ! Assign key word limits to limit arrays
  ! 
  use A0_talys_mod
  use A1_error_handling_mod

  ! Global Data
  ! 
  ! path        ! Full path to TALYS structure directory
  ! xyzLIM      ! Array to hold limits for key xyz
  ! 
  ! Local Data
  !   
  logical            :: lexist       ! Boolean for existence of file
  integer            :: istat        ! File IO state
  character(len=132) :: word(40)     ! To hold key-value pairs per line in range file
  character(len=132) :: line         ! To hold line from file temporarily
  character(len=132) :: key          ! To hold key per line
  character(len=132) :: lowerLim     ! To hold lower limit value (as string) per line
  character(len=132) :: upperLim     ! To hold upper limit value (as string) per line
  character(len=132) :: rangeFile    ! Path to keyword range file
  integer            :: ix           ! String index var
  
  ix=index(path, '/structure/')
  rangeFile=trim(path(1:ix))//'misc/key.ranges'

  ! Verify key range file can be accessed
  inquire(file=rangeFile, exist=lexist)
  if (lexist) then
    open(3, file=rangeFile, status='old', iostat=istat)
    if (istat/=0) call read_error(rangeFile, istat)
    do
      read(3, '(a132)', iostat=istat) line
      if (istat==-1) exit
      if (istat/=0) call read_error(rangeFile, istat)
      if (line(1:1)=='#') cycle
      call getkeywords(line, word)
      key=trim(word(1))
      lowerLim=trim(word(2))
      upperLim=trim(word(3))

      ! Assign numerical limits to associated LIM vars in A0_talys_mod
      if (key == "maxlevelstar") then
        read(lowerLim, *) nlevmaxLIM(1)
        nlevmaxLIM(2) = numlev
        cycle
      endif
      if (key == "maxlevelsres") then
        read(lowerLim, *) nlevmaxresLIM(1)
        nlevmaxresLIM(2) = numlev
        cycle
      endif
      if (key == "maxlevelsbin") then
        read(lowerLim, *) nlevbinLIM(1)
        nlevbinLIM(2) = numlev
        cycle
      endif
      if (key == "nlevels") then
        read(lowerLim, *) nlevLIM(1)
        nlevLIM(2) = numlev
        cycle
      endif
      if (key == "massexcess") then
        read(lowerLim, *) massexcessLIM(1)
        read(upperLim, *) massexcessLIM(2)
        cycle
      endif
      if (key == "Ltarget") then
        read(lowerLim, *) LtargetLIM(1)
        LtargetLIM(2) = numlev
        cycle
      endif
      if (key == "Liso") then
        read(lowerLim, *) LisoinpLIM(1)
        read(upperLim, *) LisoinpLIM(2)
        cycle
      endif
      if (key == "isomer") then
        read(lowerLim, *) isomerLIM(1)
        read(upperLim, *) isomerLIM(2)
        cycle
      endif
      if (key == "core") then
        read(lowerLim, *) coreLIM(1)
        read(upperLim, *) coreLIM(2)
        cycle
      endif
      if (key == "transpower") then
        read(lowerLim, *) transpowerLIM(1)
        read(upperLim, *) transpowerLIM(2)
        cycle
      endif
      if (key == "transeps") then
        read(lowerLim, *) transepsLIM(1)
        read(upperLim, *) transepsLIM(2)
        cycle
      endif
      if (key == "xseps") then
        read(lowerLim, *) xsepsLIM(1)
        read(upperLim, *) xsepsLIM(2)
        cycle
      endif
      if (key == "popeps") then
        read(lowerLim, *) popepsLIM(1)
        read(upperLim, *) popepsLIM(2)
        cycle
      endif
      if (key == "Rfiseps") then
        read(lowerLim, *) RfisepsLIM(1)
        read(upperLim, *) RfisepsLIM(2)
        cycle
      endif
      if (key == "Elow") then
        read(lowerLim, *) eninclowLIM(1)
        read(upperLim, *) eninclowLIM(2)
        cycle
      endif
      if (key == "angles") then
        read(lowerLim, *) nangleLIM(1)
        nangleLIM(2) = numang
        cycle
      endif
      if (key == "anglescont") then
        read(lowerLim, *) nanglecontLIM(1)
        nanglecontLIM(2) = numangcont
        cycle
      endif
      if (key == "anglesrec") then
        read(lowerLim, *) nanglerecLIM(1)
        nanglerecLIM(2) = numangrec
        cycle
      endif
      if (key == "maxenrec") then
        read(lowerLim, *) maxenrecLIM(1)
        maxenrecLIM(2) = numenrec
        cycle
      endif
      if (key == "maxchannel") then
        read(lowerLim, *) maxchannelLIM(1)
        read(upperLim, *) maxchannelLIM(2)
        cycle
      endif
      if (key == "massmodel") then
        read(lowerLim, *) massmodelLIM(1)
        read(upperLim, *) massmodelLIM(2)
        cycle
      endif
      if (key == "disctable") then
        read(lowerLim, *) disctableLIM(1)
        read(upperLim, *) disctableLIM(2)
        cycle
      endif
      if (key == "astroT") then
        read(lowerLim, *) astroT9LIM(1)
        read(upperLim, *) astroT9LIM(2)
        cycle
      endif
      if (key == "astroE") then
        read(lowerLim, *) astroELIM(1)
        read(upperLim, *) astroELIM(2)
        cycle
      endif
      if (key == "nonthermlev") then
        read(lowerLim, *) nonthermlevLIM(1)
        nonthermlevLIM(2) = numlev
        cycle
      endif
      if (key == "Ebeam") then
        read(lowerLim, *) EbeamLIM(1)
        EbeamLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "Eback") then
        read(lowerLim, *) EbackLIM(1)
        EbackLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "Ebeam") then
        EbeamLIM(1) = Eback
        EbeamLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "Ibeam") then
        read(lowerLim, *) IbeamLIM(1)
        read(upperLim, *) IbeamLIM(2)
        cycle
      endif
      if (key == "Area") then
        read(lowerLim, *) AreaLIM(1)
        read(upperLim, *) AreaLIM(2)
        cycle
      endif
      if (key == "Tirrad") then
        read(lowerLim, *) TirradLIM(1)
        read(upperLim, *) TirradLIM(2)
        cycle
      endif
      if (key == "Tcool") then
        read(lowerLim, *) TcoolLIM(1)
        read(upperLim, *) TcoolLIM(2)
        cycle
      endif
      if (key == "rhotarget") then
        read(lowerLim, *) rhotargetLIM(1)
        read(upperLim, *) rhotargetLIM(2)
        cycle
      endif
      if (key == "Tres") then
        read(lowerLim, *) TresLIM(1)
        read(upperLim, *) TresLIM(2)
        cycle
      endif
      if (key == "grescue") then
        read(lowerLim, *) grescueLIM(1)
        read(upperLim, *) grescueLIM(2)
        cycle
      endif
      if (key == "alphaomp") then
        read(lowerLim, *) alphaompLIM(1)
        read(upperLim, *) alphaompLIM(2)
        cycle
      endif
      if (key == "deuteronomp") then
        read(lowerLim, *) deuteronompLIM(1)
        read(upperLim, *) deuteronompLIM(2)
        cycle
      endif
      if (key == "radialmodel") then
        read(lowerLim, *) radialmodelLIM(1)
        read(upperLim, *) radialmodelLIM(2)
        cycle
      endif
      if (key == "rvadjust") then
        read(lowerLim, *) rvadjustLIM(1)
        read(upperLim, *) rvadjustLIM(2)
        cycle
      endif
      if (key == "avadjust") then
        read(lowerLim, *) avadjustLIM(1)
        read(upperLim, *) avadjustLIM(2)
        cycle
      endif
      if (key == "v1adjust") then
        read(lowerLim, *) v1adjustLIM(1)
        read(upperLim, *) v1adjustLIM(2)
        cycle
      endif
      if (key == "v2adjust") then
        read(lowerLim, *) v2adjustLIM(1)
        read(upperLim, *) v2adjustLIM(2)
        cycle
      endif
      if (key == "v3adjust") then
        read(lowerLim, *) v3adjustLIM(1)
        read(upperLim, *) v3adjustLIM(2)
        cycle
      endif
      if (key == "v4adjust") then
        read(lowerLim, *) v4adjustLIM(1)
        read(upperLim, *) v4adjustLIM(2)
        cycle
      endif
      if (key == "rwadjust") then
        read(lowerLim, *) rwadjustLIM(1)
        read(upperLim, *) rwadjustLIM(2)
        cycle
      endif
      if (key == "awadjust") then
        read(lowerLim, *) awadjustLIM(1)
        read(upperLim, *) awadjustLIM(2)
        cycle
      endif
      if (key == "w1adjust") then
        read(lowerLim, *) w1adjustLIM(1)
        read(upperLim, *) w1adjustLIM(2)
        cycle
      endif
      if (key == "w2adjust") then
        read(lowerLim, *) w2adjustLIM(1)
        read(upperLim, *) w2adjustLIM(2)
        cycle
      endif
      if (key == "w3adjust") then
        read(lowerLim, *) w3adjustLIM(1)
        read(upperLim, *) w3adjustLIM(2)
        cycle
      endif
      if (key == "w4adjust") then
        read(lowerLim, *) w4adjustLIM(1)
        read(upperLim, *) w4adjustLIM(2)
        cycle
      endif
      if (key == "rvdadjust") then
        read(lowerLim, *) rvdadjustLIM(1)
        read(upperLim, *) rvdadjustLIM(2)
        cycle
      endif
      if (key == "avdadjust") then
        read(lowerLim, *) avdadjustLIM(1)
        read(upperLim, *) avdadjustLIM(2)
        cycle
      endif
      if (key == "d1adjust") then
        read(lowerLim, *) d1adjustLIM(1)
        read(upperLim, *) d1adjustLIM(2)
        cycle
      endif
      if (key == "d2adjust") then
        read(lowerLim, *) d2adjustLIM(1)
        read(upperLim, *) d2adjustLIM(2)
        cycle
      endif
      if (key == "d3adjust") then
        read(lowerLim, *) d3adjustLIM(1)
        read(upperLim, *) d3adjustLIM(2)
        cycle
      endif
      if (key == "rwdadjust") then
        read(lowerLim, *) rwdadjustLIM(1)
        read(upperLim, *) rwdadjustLIM(2)
        cycle
      endif
      if (key == "awdadjust") then
        read(lowerLim, *) awdadjustLIM(1)
        read(upperLim, *) awdadjustLIM(2)
        cycle
      endif
      if (key == "rvsoadjust") then
        read(lowerLim, *) rvsoadjustLIM(1)
        read(upperLim, *) rvsoadjustLIM(2)
        cycle
      endif
      if (key == "avsoadjust") then
        read(lowerLim, *) avsoadjustLIM(1)
        read(upperLim, *) avsoadjustLIM(2)
        cycle
      endif
      if (key == "vso1adjust") then
        read(lowerLim, *) vso1adjustLIM(1)
        read(upperLim, *) vso1adjustLIM(2)
        cycle
      endif
      if (key == "vso2adjust") then
        read(lowerLim, *) vso2adjustLIM(1)
        read(upperLim, *) vso2adjustLIM(2)
        cycle
      endif
      if (key == "rwsoadjust") then
        read(lowerLim, *) rwsoadjustLIM(1)
        read(upperLim, *) rwsoadjustLIM(2)
        cycle
      endif
      if (key == "awsoadjust") then
        read(lowerLim, *) awsoadjustLIM(1)
        read(upperLim, *) awsoadjustLIM(2)
        cycle
      endif
      if (key == "wso1adjust") then
        read(lowerLim, *) wso1adjustLIM(1)
        read(upperLim, *) wso1adjustLIM(2)
        cycle
      endif
      if (key == "wso2adjust") then
        read(lowerLim, *) wso2adjustLIM(1)
        read(upperLim, *) wso2adjustLIM(2)
        cycle
      endif
      if (key == "rcadjust") then
        read(lowerLim, *) rcadjustLIM(1)
        read(upperLim, *) rcadjustLIM(2)
        cycle
      endif
      if (key == "ecisstep") then
        read(lowerLim, *) ecisstepLIM(1)
        read(upperLim, *) ecisstepLIM(2)
        cycle
      endif
      if (key == "ompadjustE1") then
        read(lowerLim, *) ompadjustE1LIM(1)
        ompadjustE1LIM(2) = Emaxtalys
        cycle
      endif
      if (key == "ompadjustE2") then
        read(lowerLim, *) ompadjustE2LIM(1)
        ompadjustE2LIM(2) = Emaxtalys
        cycle
      endif
      if (key == "ompadjustD") then
        read(lowerLim, *) ompadjustDLIM(1)
        read(upperLim, *) ompadjustDLIM(2)
        cycle
      endif
      if (key == "ompadjusts") then
        read(lowerLim, *) ompadjustsLIM(1)
        read(upperLim, *) ompadjustsLIM(2)
        cycle
      endif
      if (key == "jlmmode") then
        read(lowerLim, *) jlmmodeLIM(1)
        read(upperLim, *) jlmmodeLIM(2)
        cycle
      endif
      if (key == "lvadjust") then
        read(lowerLim, *) lvadjustLIM(1)
        read(upperLim, *) lvadjustLIM(2)
        cycle
      endif
      if (key == "lwadjust") then
        read(lowerLim, *) lwadjustLIM(1)
        read(upperLim, *) lwadjustLIM(2)
        cycle
      endif
      if (key == "lv1adjust") then
        read(lowerLim, *) lv1adjustLIM(1)
        read(upperLim, *) lv1adjustLIM(2)
        cycle
      endif
      if (key == "lw1adjust") then
        read(lowerLim, *) lw1adjustLIM(1)
        read(upperLim, *) lw1adjustLIM(2)
        cycle
      endif
      if (key == "lvsoadjust") then
        read(lowerLim, *) lvsoadjustLIM(1)
        read(upperLim, *) lvsoadjustLIM(2)
        cycle
      endif
      if (key == "lwsoadjust") then
        read(lowerLim, *) lwsoadjustLIM(1)
        read(upperLim, *) lwsoadjustLIM(2)
        cycle
      endif
      if (key == "aradialcor") then
        read(lowerLim, *) aradialcorLIM(1)
        read(upperLim, *) aradialcorLIM(2)
        cycle
      endif
      if (key == "adepthcor") then
        read(lowerLim, *) adepthcorLIM(1)
        read(upperLim, *) adepthcorLIM(2)
        cycle
      endif
      if (key == "soswitch") then
        read(lowerLim, *) soswitchLIM(1)
        read(upperLim, *) soswitchLIM(2)
        cycle
      endif
      if (key == "Ejoin") then
        read(lowerLim, *) EjoinLIM(1)
        EjoinLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "Vinfadjust") then
        read(lowerLim, *) VinfadjustLIM(1)
        read(upperLim, *) VinfadjustLIM(2)
        cycle
      endif
      if (key == "pruittset") then
        read(lowerLim, *) pruittsetLIM(1)
        read(upperLim, *) pruittsetLIM(2)
        cycle
      endif
      if (key == "maxband") then
        read(lowerLim, *) maxbandLIM(1)
        read(upperLim, *) maxbandLIM(2)
        cycle
      endif
      if (key == "maxrot") then
        read(lowerLim, *) maxrotLIM(1)
        read(upperLim, *) maxrotLIM(2)
        cycle
      endif
      if (key == "ewfc") then
        read(lowerLim, *) ewfcLIM(1)
        read(upperLim, *) ewfcLIM(2)
        cycle
      endif
      if (key == "eurr") then
        read(lowerLim, *) eurrLIM(1)
        read(upperLim, *) eurrLIM(2)
        cycle
      endif
      if (key == "wmode") then
        read(lowerLim, *) wmodeLIM(1)
        read(upperLim, *) wmodeLIM(2)
        cycle
      endif
      if (key == "wfcfactor") then
        read(lowerLim, *) wfcfactorLIM(1)
        read(upperLim, *) wfcfactorLIM(2)
        cycle
      endif
      if (key == "lurr") then
        read(lowerLim, *) lurrLIM(1)
        lurrLIM(2) = numl
        cycle
      endif
      if (key == "gammax") then
        read(lowerLim, *) gammaxLIM(1)
        read(upperLim, *) gammaxLIM(2)
        cycle
      endif
      if (key == "strength") then
        read(lowerLim, *) strengthLIM(1)
        read(upperLim, *) strengthLIM(2)
        cycle
      endif
      if (key == "etable") then
        read(lowerLim, *) etableLIM(1)
        read(upperLim, *) etableLIM(2)
        cycle
      endif
      if (key == "ftable") then
        read(lowerLim, *) ftableLIM(1)
        read(upperLim, *) ftableLIM(2)
        cycle
      endif
      if (key == "wtable") then
        read(lowerLim, *) wtableLIM(1)
        read(upperLim, *) wtableLIM(2)
        cycle
      endif
      if (key == "etableadjust") then
        read(lowerLim, *) etableadjustLIM(1)
        read(upperLim, *) etableadjustLIM(2)
        cycle
      endif
      if (key == "ftableadjust") then
        read(lowerLim, *) ftableadjustLIM(1)
        read(upperLim, *) ftableadjustLIM(2)
        cycle
      endif
      if (key == "wtableadjust") then
        read(lowerLim, *) wtableadjustLIM(1)
        read(upperLim, *) wtableadjustLIM(2)
        cycle
      endif
      if (key == "egr") then
        read(lowerLim, *) egrLIM(1)
        read(upperLim, *) egrLIM(2)
        cycle
      endif
      if (key == "ggr") then
        read(lowerLim, *) ggrLIM(1)
        read(upperLim, *) ggrLIM(2)
        cycle
      endif
      if (key == "sgr") then
        read(lowerLim, *) sgrLIM(1)
        read(upperLim, *) sgrLIM(2)
        cycle
      endif
      if (key == "epr") then
        read(lowerLim, *) eprLIM(1)
        read(upperLim, *) eprLIM(2)
        cycle
      endif
      if (key == "gpr") then
        read(lowerLim, *) gprLIM(1)
        read(upperLim, *) gprLIM(2)
        cycle
      endif
      if (key == "spr") then
        read(lowerLim, *) tprLIM(1)
        read(upperLim, *) tprLIM(2)
        cycle
      endif
      if (key == "egradjust") then
        read(lowerLim, *) egradjustLIM(1)
        read(upperLim, *) egradjustLIM(2)
        cycle
      endif
      if (key == "ggradjust") then
        read(lowerLim, *) ggradjustLIM(1)
        read(upperLim, *) ggradjustLIM(2)
        cycle
      endif
      if (key == "sgradjust") then
        read(lowerLim, *) sgradjustLIM(1)
        read(upperLim, *) sgradjustLIM(2)
        cycle
      endif
      if (key == "epradjust") then
        read(lowerLim, *) epradjustLIM(1)
        read(upperLim, *) epradjustLIM(2)
        cycle
      endif
      if (key == "gpradjust") then
        read(lowerLim, *) gpradjustLIM(1)
        read(upperLim, *) gpradjustLIM(2)
        cycle
      endif
      if (key == "spradjust") then
        read(lowerLim, *) tpradjustLIM(1)
        read(upperLim, *) tpradjustLIM(2)
        cycle
      endif
      if (key == "upbendc") then
        read(lowerLim, *) upbendcLIM(1)
        read(upperLim, *) upbendcLIM(2)
        cycle
      endif
      if (key == "upbende") then
        read(lowerLim, *) upbendeLIM(1)
        read(upperLim, *) upbendeLIM(2)
        cycle
      endif
      if (key == "upbendf") then
        read(lowerLim, *) upbendfLIM(1)
        read(upperLim, *) upbendfLIM(2)
        cycle
      endif
      if (key == "upbendcadjust") then
        read(lowerLim, *) upbendcadjustLIM(1)
        read(upperLim, *) upbendcadjustLIM(2)
        cycle
      endif
      if (key == "upbendeadjust") then
        read(lowerLim, *) upbendeadjustLIM(1)
        read(upperLim, *) upbendeadjustLIM(2)
        cycle
      endif
      if (key == "upbendfadjust") then
        read(lowerLim, *) upbendfadjustLIM(1)
        read(upperLim, *) upbendfadjustLIM(2)
        cycle
      endif
      if (key == "gamgam") then
        read(lowerLim, *) gamgamLIM(1)
        read(upperLim, *) gamgamLIM(2)
        cycle
      endif
      if (key == "D0") then
        read(lowerLim, *) D0LIM(1)
        read(upperLim, *) D0LIM(2)
        cycle
      endif
      if (key == "fiso") then
        read(lowerLim, *) fisoLIM(1)
        read(upperLim, *) fisoLIM(2)
        cycle
      endif
      if (key == "fisom") then
        read(lowerLim, *) fisomLIM(1)
        read(upperLim, *) fisomLIM(2)
        cycle
      endif
      if (key == "gamgamadjust") then
        read(lowerLim, *) gamgamadjustLIM(1)
        read(upperLim, *) gamgamadjustLIM(2)
        cycle
      endif
      if (key == "RprimeU") then
        read(lowerLim, *) RprimeULIM(1)
        read(upperLim, *) RprimeULIM(2)
        cycle
      endif
      if (key == "ldmodelracap") then
        read(lowerLim, *) ldmodelracapLIM(1)
        read(upperLim, *) ldmodelracapLIM(2)
        cycle
      endif
      if (key == "levinger") then
        read(lowerLim, *) levingerLIM(1)
        read(upperLim, *) levingerLIM(2)
        cycle
      endif
      if (key == "sfth") then
        read(lowerLim, *) spectfacthLIM(1)
        read(upperLim, *) spectfacthLIM(2)
        cycle
      endif
      if (key == "sfexp") then
        read(lowerLim, *) spectfacexpLIM(1)
        read(upperLim, *) spectfacexpLIM(2)
        cycle
      endif
      if (key == "epreeq") then
        read(lowerLim, *) epreeqLIM(1)
        epreeqLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "preeqmode") then
        read(lowerLim, *) preeqmodeLIM(1)
        read(upperLim, *) preeqmodeLIM(2)
        cycle
      endif
      if (key == "mpreeqmode") then
        read(lowerLim, *) mpreeqmodeLIM(1)
        read(upperLim, *) mpreeqmodeLIM(2)
        cycle
      endif
      if (key == "breakupmodel") then
        read(lowerLim, *) breakupmodelLIM(1)
        read(upperLim, *) breakupmodelLIM(2)
        cycle
      endif
      if (key == "phmodel") then
        read(lowerLim, *) phmodelLIM(1)
        read(upperLim, *) phmodelLIM(2)
        cycle
      endif
      if (key == "pairmodel") then
        read(lowerLim, *) pairmodelLIM(1)
        read(upperLim, *) pairmodelLIM(2)
        cycle
      endif
      if (key == "pespinmodel") then
        read(lowerLim, *) pespinmodelLIM(1)
        read(upperLim, *) pespinmodelLIM(2)
        cycle
      endif
      if (key == "emulpre") then
        read(lowerLim, *) emulpreLIM(1)
        emulpreLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "M2constant") then
        read(lowerLim, *) M2constantLIM(1)
        read(upperLim, *) M2constantLIM(2)
        cycle
      endif
      if (key == "M2limit") then
        read(lowerLim, *) M2limitLIM(1)
        read(upperLim, *) M2limitLIM(2)
        cycle
      endif
      if (key == "M2shift") then
        read(lowerLim, *) M2shiftLIM(1)
        read(upperLim, *) M2shiftLIM(2)
        cycle
      endif
      if (key == "Rpipi") then
        read(lowerLim, *) RpipiLIM(1)
        read(upperLim, *) RpipiLIM(2)
        cycle
      endif
      if (key == "Rnunu") then
        read(lowerLim, *) RnunuLIM(1)
        read(upperLim, *) RnunuLIM(2)
        cycle
      endif
      if (key == "Rpinu") then
        read(lowerLim, *) RpinuLIM(1)
        read(upperLim, *) RpinuLIM(2)
        cycle
      endif
      if (key == "Rnupi") then
        read(lowerLim, *) RnupiLIM(1)
        read(upperLim, *) RnupiLIM(2)
        cycle
      endif
      if (key == "Rgamma") then
        read(lowerLim, *) RgammaLIM(1)
        read(upperLim, *) RgammaLIM(2)
        cycle
      endif
      if (key == "Esurf") then
        read(lowerLim, *) Esurf0LIM(1)
        read(upperLim, *) Esurf0LIM(2)
        cycle
      endif
      if (key == "msdbins") then
        read(lowerLim, *) msdbinsLIM(1)
        msdbinsLIM(2) = numenmsd/2 - 1
        cycle
      endif
      if (key == "Emsdmin") then
        read(lowerLim, *) EmsdminLIM(1)
        EmsdminLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "elwidth") then
        read(lowerLim, *) elwidthLIM(1)
        read(upperLim, *) elwidthLIM(2)
        cycle
      endif
      if (key == "xscaptherm") then
        read(lowerLim, *) xscapthermLIM(1)
        read(upperLim, *) xscapthermLIM(2)
        cycle
      endif
      if (key == "xsptherm") then
        read(lowerLim, *) xspthermLIM(1)
        read(upperLim, *) xspthermLIM(2)
        cycle
      endif
      if (key == "xsalphatherm") then
        read(lowerLim, *) xsalphathermLIM(1)
        read(upperLim, *) xsalphathermLIM(2)
        cycle
      endif
      if (key == "Cstrip") then
        read(lowerLim, *) CstripLIM(1)
        read(upperLim, *) CstripLIM(2)
        cycle
      endif
      if (key == "Cknock") then
        read(lowerLim, *) CknockLIM(1)
        read(upperLim, *) CknockLIM(2)
        cycle
      endif
      if (key == "Cbreak") then
        read(lowerLim, *) CbreakLIM(1)
        read(upperLim, *) CbreakLIM(2)
        cycle
      endif
      if (key == "GMRadjustE") then
        read(lowerLim, *) GMRadjustELIM(1)
        read(upperLim, *) GMRadjustELIM(2)
        cycle
      endif
      if (key == "GQRadjustE") then
        read(lowerLim, *) GQRadjustELIM(1)
        read(upperLim, *) GQRadjustELIM(2)
        cycle
      endif
      if (key == "LEORadjustE") then
        read(lowerLim, *) LEORadjustELIM(1)
        read(upperLim, *) LEORadjustELIM(2)
        cycle
      endif
      if (key == "HEORadjustE") then
        read(lowerLim, *) HEORadjustELIM(1)
        read(upperLim, *) HEORadjustELIM(2)
        cycle
      endif
      if (key == "GMRadjustG") then
        read(lowerLim, *) GMRadjustGLIM(1)
        read(upperLim, *) GMRadjustGLIM(2)
        cycle
      endif
      if (key == "GQRadjustG") then
        read(lowerLim, *) GQRadjustGLIM(1)
        read(upperLim, *) GQRadjustGLIM(2)
        cycle
      endif
      if (key == "LEORadjustG") then
        read(lowerLim, *) LEORadjustGLIM(1)
        read(upperLim, *) LEORadjustGLIM(2)
        cycle
      endif
      if (key == "HEORadjustG") then
        read(lowerLim, *) HEORadjustGLIM(1)
        read(upperLim, *) HEORadjustGLIM(2)
        cycle
      endif
      if (key == "GMRadjustD") then
        read(lowerLim, *) GMRadjustDLIM(1)
        read(upperLim, *) GMRadjustDLIM(2)
        cycle
      endif
      if (key == "GQRadjustD") then
        read(lowerLim, *) GQRadjustDLIM(1)
        read(upperLim, *) GQRadjustDLIM(2)
        cycle
      endif
      if (key == "LEORadjustD") then
        read(lowerLim, *) LEORadjustDLIM(1)
        read(upperLim, *) LEORadjustDLIM(2)
        cycle
      endif
      if (key == "HEORadjustD") then
        read(lowerLim, *) HEORadjustDLIM(1)
        read(upperLim, *) HEORadjustDLIM(2)
        cycle
      endif
      if (key == "spincutmodel") then
        read(lowerLim, *) spincutmodelLIM(1)
        read(upperLim, *) spincutmodelLIM(2)
        cycle
      endif
      if (key == "shellmodel") then
        read(lowerLim, *) shellmodelLIM(1)
        read(upperLim, *) shellmodelLIM(2)
        cycle
      endif
      if (key == "kvibmodel") then
        read(lowerLim, *) kvibmodelLIM(1)
        read(upperLim, *) kvibmodelLIM(2)
        cycle
      endif
      if (key == "ldmodelcn") then
        read(lowerLim, *) ldmodelCNLIM(1)
        read(upperLim, *) ldmodelCNLIM(2)
        cycle
      endif
      if (key == "ldmodel") then
        read(lowerLim, *) ldmodelLIM(1)
        read(upperLim, *) ldmodelLIM(2)
        cycle
      endif
      if (key == "a") then
        read(lowerLim, *) alevLIM(1)
        read(upperLim, *) alevLIM(2)
        cycle
      endif
      if (key == "alimit") then
        read(lowerLim, *) alimitLIM(1)
        read(upperLim, *) alimitLIM(2)
        cycle
      endif
      if (key == "gammald") then
        read(lowerLim, *) gammaldLIM(1)
        read(upperLim, *) gammaldLIM(2)
        cycle
      endif
      if (key == "risomer") then
        read(lowerLim, *) RisomerLIM(1)
        read(upperLim, *) RisomerLIM(2)
        cycle
      endif
      if (key == "deltaW") then
        read(lowerLim, *) deltaWLIM(1)
        read(upperLim, *) deltaWLIM(2)
        cycle
      endif
      if (key == "Nlow") then
        read(lowerLim, *) NlowLIM(1)
        read(upperLim, *) NlowLIM(2)
        cycle
      endif
      if (key == "Ntop") then
        read(lowerLim, *) NtopLIM(1)
        read(upperLim, *) NtopLIM(2)
        cycle
      endif
      if (key == "E0") then
        read(lowerLim, *) E0LIM(1)
        read(upperLim, *) E0LIM(2)
        cycle
      endif
      if (key == "beta2") then
        read(lowerLim, *) beta2LIM(1)
        read(upperLim, *) beta2LIM(2)
        cycle
      endif
      if (key == "s2adjust") then
        read(lowerLim, *) s2adjustLIM(1)
        read(upperLim, *) s2adjustLIM(2)
        cycle
      endif
      if (key == "Krotconstant") then
        read(lowerLim, *) KrotconstantLIM(1)
        read(upperLim, *) KrotconstantLIM(2)
        cycle
      endif
      if (key == "Ufermi") then
        read(lowerLim, *) UfermiLIM(1)
        read(upperLim, *) UfermiLIM(2)
        cycle
      endif
      if (key == "cfermi") then
        read(lowerLim, *) cfermiLIM(1)
        read(upperLim, *) cfermiLIM(2)
        cycle
      endif
      if (key == "T") then
        read(lowerLim, *) TLIM(1)
        read(upperLim, *) TLIM(2)
        cycle
      endif
      if (key == "Exmatch") then
        read(lowerLim, *) ExmatchLIM(1)
        read(upperLim, *) ExmatchLIM(2)
        cycle
      endif
      if (key == "Tadjust") then
        read(lowerLim, *) TadjustLIM(1)
        read(upperLim, *) TadjustLIM(2)
        cycle
      endif
      if (key == "E0adjust") then
        read(lowerLim, *) E0adjustLIM(1)
        read(upperLim, *) E0adjustLIM(2)
        cycle
      endif
      if (key == "Exmatchadjust") then
        read(lowerLim, *) ExmatchadjustLIM(1)
        read(upperLim, *) ExmatchadjustLIM(2)
        cycle
      endif
      if (key == "Pshift") then
        read(lowerLim, *) PshiftLIM(1)
        read(upperLim, *) PshiftLIM(2)
        cycle
      endif
      if (key == "Pshiftadjust") then
        read(lowerLim, *) PshiftadjustLIM(1)
        read(upperLim, *) PshiftadjustLIM(2)
        cycle
      endif
      if (key == "ctable") then
        read(lowerLim, *) ctableLIM(1)
        read(upperLim, *) ctableLIM(2)
        cycle
      endif
      if (key == "ptable") then
        read(lowerLim, *) ptableLIM(1)
        read(upperLim, *) ptableLIM(2)
        cycle
      endif
      if (key == "ctableadjust") then
        read(lowerLim, *) ctableadjustLIM(1)
        read(upperLim, *) ctableadjustLIM(2)
        cycle
      endif
      if (key == "ptableadjust") then
        read(lowerLim, *) ptableadjustLIM(1)
        read(upperLim, *) ptableadjustLIM(2)
        cycle
      endif
      if (key == "aadjust") then
        read(lowerLim, *) aadjustLIM(1)
        read(upperLim, *) aadjustLIM(2)
        cycle
      endif
      if (key == "gadjust") then
        read(lowerLim, *) gadjustLIM(1)
        read(upperLim, *) gadjustLIM(2)
        cycle
      endif
      if (key == "gnadjust") then
        read(lowerLim, *) gnadjustLIM(1)
        read(upperLim, *) gnadjustLIM(2)
        cycle
      endif
      if (key == "gpadjust") then
        read(lowerLim, *) gpadjustLIM(1)
        read(upperLim, *) gpadjustLIM(2)
        cycle
      endif
      if (key == "pair") then
        read(lowerLim, *) pairLIM(1)
        read(upperLim, *) pairLIM(2)
        cycle
      endif
      if (key == "g") then
        read(lowerLim, *) gLIM(1)
        read(upperLim, *) gLIM(2)
        cycle
      endif
      if (key == "gn") then
        read(lowerLim, *) gnLIM(1)
        read(upperLim, *) gnLIM(2)
        cycle
      endif
      if (key == "gp") then
        read(lowerLim, *) gpLIM(1)
        read(upperLim, *) gpLIM(2)
        cycle
      endif
      if (key == "alphald") then
        read(lowerLim, *) alphaldLIM(1)
        read(upperLim, *) alphaldLIM(2)
        cycle
      endif
      if (key == "betald") then
        read(lowerLim, *) betaldLIM(1)
        read(upperLim, *) betaldLIM(2)
        cycle
      endif
      if (key == "gammashell1") then
        read(lowerLim, *) gammashell1LIM(1)
        read(upperLim, *) gammashell1LIM(2)
        cycle
      endif
      if (key == "Pshiftconstant") then
        read(lowerLim, *) PshiftconstantLIM(1)
        read(upperLim, *) PshiftconstantLIM(2)
        cycle
      endif
      if (key == "cglobal") then
        read(lowerLim, *) cglobalLIM(1)
        read(upperLim, *) cglobalLIM(2)
        cycle
      endif
      if (key == "pglobal") then
        read(lowerLim, *) pglobalLIM(1)
        read(upperLim, *) pglobalLIM(2)
        cycle
      endif
      if (key == "gammashell2") then
        read(lowerLim, *) gammashell2LIM(1)
        read(upperLim, *) gammashell2LIM(2)
        cycle
      endif
      if (key == "pairconstant") then
        read(lowerLim, *) pairconstantLIM(1)
        read(upperLim, *) pairconstantLIM(2)
        cycle
      endif
      if (key == "Kph") then
        read(lowerLim, *) KphLIM(1)
        read(upperLim, *) KphLIM(2)
        cycle
      endif
      if (key == "Rspincut") then
        read(lowerLim, *) RspincutLIM(1)
        read(upperLim, *) RspincutLIM(2)
        cycle
      endif
      if (key == "Rspincutpreeq") then
        read(lowerLim, *) RspincutpreeqLIM(1)
        read(upperLim, *) RspincutpreeqLIM(2)
        cycle
      endif
      if (key == "Rspincutff") then
        read(lowerLim, *) RspincutffLIM(1)
        read(upperLim, *) RspincutffLIM(2)
        cycle
      endif
      if (key == "fismodel") then
        read(lowerLim, *) fismodelLIM(1)
        read(upperLim, *) fismodelLIM(2)
        cycle
      endif
      if (key == "fismodelalt") then
        read(lowerLim, *) fismodelaltLIM(1)
        read(upperLim, *) fismodelaltLIM(2)
        cycle
      endif
      if (key == "fymodel") then
        read(lowerLim, *) fymodelLIM(1)
        read(upperLim, *) fymodelLIM(2)
        cycle
      endif
      if (key == "ffmodel") then
        read(lowerLim, *) ffmodelLIM(1)
        read(upperLim, *) ffmodelLIM(2)
        cycle
      endif
      if (key == "pfnsmodel") then
        read(lowerLim, *) pfnsmodelLIM(1)
        read(upperLim, *) pfnsmodelLIM(2)
        cycle
      endif
      if (key == "gefran") then
        read(lowerLim, *) gefranLIM(1)
        read(upperLim, *) gefranLIM(2)
        cycle
      endif
      if (key == "Cnubar1") then
        read(lowerLim, *) Cnubar1LIM(1)
        read(upperLim, *) Cnubar1LIM(2)
        cycle
      endif
      if (key == "Cnubar2") then
        read(lowerLim, *) Cnubar2LIM(1)
        read(upperLim, *) Cnubar2LIM(2)
        cycle
      endif
      if (key == "Tmadjust") then
        read(lowerLim, *) TmadjustLIM(1)
        read(upperLim, *) TmadjustLIM(2)
        cycle
      endif
      if (key == "Fsadjust") then
        read(lowerLim, *) FsadjustLIM(1)
        read(upperLim, *) FsadjustLIM(2)
        cycle
      endif
      if (key == "Cbarrier") then
        read(lowerLim, *) CbarrierLIM(1)
        read(upperLim, *) CbarrierLIM(2)
        cycle
      endif
      if (key == "axtype") then
        read(lowerLim, *) axtypeLIM(1)
        read(upperLim, *) axtypeLIM(2)
        cycle
      endif
      if (key == "fisbar") then
        read(lowerLim, *) fbarrierLIM(1)
        read(upperLim, *) fbarrierLIM(2)
        cycle
      endif
      if (key == "fisbaradjust") then
        read(lowerLim, *) fbaradjustLIM(1)
        read(upperLim, *) fbaradjustLIM(2)
        cycle
      endif
      if (key == "fishw") then
        read(lowerLim, *) fwidthLIM(1)
        read(upperLim, *) fwidthLIM(2)
        cycle
      endif
      if (key == "fishwadjust") then
        read(lowerLim, *) fwidthadjustLIM(1)
        read(upperLim, *) fwidthadjustLIM(2)
        cycle
      endif
      if (key == "bdamp") then
        read(lowerLim, *) bdampLIM(1)
        read(upperLim, *) bdampLIM(2)
        cycle
      endif
      if (key == "bdampadjust") then
        read(lowerLim, *) bdampadjustLIM(1)
        read(upperLim, *) bdampadjustLIM(2)
        cycle
      endif
      if (key == "Rtransmom") then
        read(lowerLim, *) RtransmomLIM(1)
        read(upperLim, *) RtransmomLIM(2)
        cycle
      endif
      if (key == "Rclass2mom") then
        read(lowerLim, *) Rclass2momLIM(1)
        read(upperLim, *) Rclass2momLIM(2)
        cycle
      endif
      if (key == "betafiscor") then
        read(lowerLim, *) betafiscorLIM(1)
        read(upperLim, *) betafiscorLIM(2)
        cycle
      endif
      if (key == "rmiufiscor") then
        read(lowerLim, *) rmiufiscorLIM(1)
        read(upperLim, *) rmiufiscorLIM(2)
        cycle
      endif
      if (key == "vfiscor") then
        read(lowerLim, *) vfiscorLIM(1)
        read(upperLim, *) vfiscorLIM(2)
        cycle
      endif
      if (key == "betafiscoradjust") then
        read(lowerLim, *) betafiscoradjustLIM(1)
        read(upperLim, *) betafiscoradjustLIM(2)
        cycle
      endif
      if (key == "vfiscoradjust") then
        read(lowerLim, *) vfiscoradjustLIM(1)
        read(upperLim, *) vfiscoradjustLIM(2)
        cycle
      endif
      if (key == "rmiufiscoradjust") then
        read(lowerLim, *) rmiufiscoradjustLIM(1)
        read(upperLim, *) rmiufiscoradjustLIM(2)
        cycle
      endif
      if (key == "eadd") then
        read(lowerLim, *) eaddLIM(1)
        eaddLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "eaddel") then
        read(lowerLim, *) eaddelLIM(1)
        eaddelLIM(2) = Emaxtalys
        cycle
      endif
      if (key == "ddxmode") then
        read(lowerLim, *) ddxmodeLIM(1)
        read(upperLim, *) ddxmodeLIM(2)
        cycle
      endif
      if (key == "fileddxe") then
        read(lowerLim, *) fileddxeLIM(1)
        fileddxeLIM(2) = enincmax
        cycle
      endif
      if (key == "fileddxa") then
        read(lowerLim, *) fileddxaLIM(1)
        read(upperLim, *) fileddxaLIM(2)
        cycle
      endif
      if (key == "D") then
        read(lowerLim, *) DLIM(1)
        read(upperLim, *) DLIM(2)
        cycle
      endif
      if (key == "Emaxpseudores") then
        read(lowerLim, *) EmaxpseudoresLIM(1)
        read(upperLim, *) EmaxpseudoresLIM(2)
        cycle
      end if
      if (key == "pseudoreswidth") then
        read(lowerLim, *) pseudoreswidthLIM(1)
        read(upperLim, *) pseudoreswidthLIM(2)
        cycle
      end if
      if (key == "pseudoresfade") then
        read(lowerLim, *) pseudoresfadeLIM(1)
        read(upperLim, *) pseudoresfadeLIM(2)
        cycle
      end if
      if (key == "class2width") then
        read(lowerLim, *) widthc2LIM(1)
        read(upperLim, *) widthc2LIM(2)
        cycle
      end if
      ! Terminate if there is no match
      write(*, '( " TALYS-error: The following key in key.ranges can", &
      & "not be found in loadkeyranges.f90: ", a30)') key
      stop
    end do
    ! For keys whose limits are derived from others, assign here
    massnucleusLIM(1) = real(A) - 1.
    massnucleusLIM(2) = real(A) + 1.
    close(unit = 3)
  end if
end subroutine loadkeyranges