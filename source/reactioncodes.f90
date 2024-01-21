subroutine reactioncodes
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Relate reaction strings, MT numbers and TALYS output files
!
! Author    : Arjan Koning
!
! 2023-08-08: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_talys_mod
!
! *** Declaration of local data
!
  integer :: i
!
! ************************** Initialization ****************************
!
! Relate MT numbers to ejectiles, also used for TALYS exclusive output files
!
  MTchannel = -1
  MTchannel(4) = 100000
  MTchannel(11) = 201000
  MTchannel(16) = 200000
  MTchannel(17) = 300000
  MTchannel(18) = -99
  MTchannel(22) = 100001
  MTchannel(23) = 100003
  MTchannel(24) = 200001
  MTchannel(25) = 300001
  MTchannel(28) = 110000
  MTchannel(29) = 100002
  MTchannel(30) = 200002
  MTchannel(32) = 101000
  MTchannel(33) = 100100
  MTchannel(34) = 100010
  MTchannel(35) = 101002
  MTchannel(36) = 100102
  MTchannel(37) = 400000
  MTchannel(41) = 210000
  MTchannel(42) = 310000
  MTchannel(44) = 120000
  MTchannel(45) = 110001
  MTchannel(91) = 100000
  MTchannel(102) = 000000
  MTchannel(103) = 010000
  MTchannel(104) = 001000
  MTchannel(105) = 000100
  MTchannel(106) = 000010
  MTchannel(107) = 000001
  MTchannel(108) = 000002
  MTchannel(109) = 000003
  MTchannel(111) = 020000
  MTchannel(112) = 010001
  MTchannel(113) = 000102
  MTchannel(114) = 001002
  MTchannel(115) = 011000
  MTchannel(116) = 010100
  MTchannel(117) = 001001
  MTchannel(152) = 500000
  MTchannel(153) = 600000
  MTchannel(154) = 200100
  MTchannel(155) = 000101
  MTchannel(156) = 410000
  MTchannel(157) = 301000
  MTchannel(158) = 101001
  MTchannel(159) = 210001
  MTchannel(160) = 700000
  MTchannel(161) = 800000
  MTchannel(162) = 510000
  MTchannel(163) = 610000
  MTchannel(164) = 710000
  MTchannel(165) = 400001
  MTchannel(166) = 500001
  MTchannel(167) = 600001
  MTchannel(168) = 700001
  MTchannel(169) = 401000
  MTchannel(170) = 501000
  MTchannel(171) = 601000
  MTchannel(172) = 300100
  MTchannel(173) = 400100
  MTchannel(174) = 500100
  MTchannel(175) = 600100
  MTchannel(176) = 200010
  MTchannel(177) = 300010
  MTchannel(178) = 400010
  MTchannel(179) = 320000
  MTchannel(180) = 300002
  MTchannel(181) = 310001
  MTchannel(182) = 001100
  MTchannel(183) = 111000
  MTchannel(184) = 110100
  MTchannel(185) = 101100
  MTchannel(186) = 110010
  MTchannel(187) = 101010
  MTchannel(188) = 100110
  MTchannel(189) = 100101
  MTchannel(190) = 220000
  MTchannel(191) = 010010
  MTchannel(192) = 001010
  MTchannel(193) = 000011
  MTchannel(194) = 420000
  MTchannel(195) = 400002
  MTchannel(196) = 410001
  MTchannel(197) = 030000
  MTchannel(198) = 130000
  MTchannel(199) = 320001
  MTchannel(200) = 520000
  MTchannel(649) = 010000
  MTchannel(699) = 001000
  MTchannel(749) = 000100
  MTchannel(799) = 000010
  MTchannel(849) = 000001
!
! Reaction strings
!
  reactionstring=''
  reactionstring(1)='(n,tot)'
  reactionstring(2)='(n,el)'
  reactionstring(3)='(n,non)'
  reactionstring(4)="(n,n')"
  reactionstring(5)='(n,x)'
  reactionstring(11)='(n,2nd)'
  reactionstring(16)='(n,2n)'
  reactionstring(17)='(n,3n)'
  reactionstring(18)='(n,f)'
  reactionstring(22)='(n,na)'
  reactionstring(23)='(n,n3a)'
  reactionstring(24)='(n,2na)'
  reactionstring(25)='(n,3na)'
  reactionstring(28)='(n,np)'
  reactionstring(29)='(n,n2a)'
  reactionstring(30)='(n,2n2a)'
  reactionstring(32)='(n,nd)'
  reactionstring(33)='(n,nt)'
  reactionstring(34)='(n,nh)'
  reactionstring(35)='(n,nd2a)'
  reactionstring(36)='(n,nt2a)'
  reactionstring(37)='(n,4n)'
  reactionstring(41)='(n,2np)'
  reactionstring(42)='(n,3np)'
  reactionstring(44)='(n,n2p)'
  reactionstring(45)='(n,npa)'
  do i=0,40
    reactionstring(50+i)="(n,n'_01) "
    write(reactionstring(50+i)(7:8),'(i2.2)') i
  enddo
  reactionstring(91)="(n,n'_con)"
  reactionstring(101)='(n,abs)'
  reactionstring(102)='(n,g)'
  reactionstring(103)='(n,p)'
  reactionstring(104)='(n,d)'
  reactionstring(105)='(n,t)'
  reactionstring(106)='(n,h)'
  reactionstring(107)='(n,a)'
  reactionstring(108)='(n,2a)'
  reactionstring(109)='(n,3a)'
  reactionstring(111)='(n,2p)'
  reactionstring(112)='(n,pa)'
  reactionstring(113)='(n,t2a)'
  reactionstring(114)='(n,d2a)'
  reactionstring(115)='(n,pd)'
  reactionstring(116)='(n,pt)'
  reactionstring(117)='(n,da)'
  reactionstring(152)='(n,5n)'
  reactionstring(153)='(n,6n)'
  reactionstring(154)='(n,2nt)'
  reactionstring(155)='(n,ta)'
  reactionstring(156)='(n,4np)'
  reactionstring(157)='(n,3nd)'
  reactionstring(158)='(n,nda)'
  reactionstring(159)='(n,2npa)'
  reactionstring(160)='(n,7n)'
  reactionstring(161)='(n,8n)'
  reactionstring(162)='(n,5np)'
  reactionstring(163)='(n,6np)'
  reactionstring(164)='(n,7np)'
  reactionstring(165)='(n,4na)'
  reactionstring(166)='(n,5na)'
  reactionstring(167)='(n,6na)'
  reactionstring(168)='(n,7na)'
  reactionstring(169)='(n,4nd)'
  reactionstring(170)='(n,5nd)'
  reactionstring(171)='(n,6nd)'
  reactionstring(172)='(n,3nt)'
  reactionstring(173)='(n,4nt)'
  reactionstring(174)='(n,5nt)'
  reactionstring(175)='(n,6nt)'
  reactionstring(176)='(n,2nh)'
  reactionstring(177)='(n,3nh)'
  reactionstring(178)='(n,4nh)'
  reactionstring(179)='(n,3n2p)'
  reactionstring(180)='(n,3n2a)'
  reactionstring(181)='(n,3npa)'
  reactionstring(182)='(n,dt)'
  reactionstring(183)='(n,npd)'
  reactionstring(184)='(n,npt)'
  reactionstring(185)='(n,ndt)'
  reactionstring(186)='(n,nph)'
  reactionstring(187)='(n,ndh)'
  reactionstring(188)='(n,nth)'
  reactionstring(189)='(n,nta)'
  reactionstring(190)='(n,2n2p)'
  reactionstring(191)='(n,ph)'
  reactionstring(192)='(n,dh)'
  reactionstring(193)='(n,ha)'
  reactionstring(194)='(n,4n2p)'
  reactionstring(195)='(n,4n2a)'
  reactionstring(196)='(n,4npa)'
  reactionstring(197)='(n,3p)'
  reactionstring(198)='(n,n3p)'
  reactionstring(199)='(n,3n2pa)'
  reactionstring(200)='(n,5n2p)'
  reactionstring(201)='(n,xn)'
  reactionstring(202)='(n,xg)'
  reactionstring(203)='(n,xp)'
  reactionstring(204)='(n,xd)'
  reactionstring(205)='(n,xt)'
  reactionstring(206)='(n,xh)'
  reactionstring(207)='(n,xa)'
  do i=0,48
    reactionstring(600+i)="(n,p_00)"
    write(reactionstring(600+i)(6:7),'(i2.2)') i
    reactionstring(650+i)="(n,d_00)"
    write(reactionstring(650+i)(6:7),'(i2.2)') i
    reactionstring(700+i)="(n,t_00)"
    write(reactionstring(700+i)(6:7),'(i2.2)') i
    reactionstring(750+i)="(n,h_00)"
    write(reactionstring(750+i)(6:7),'(i2.2)') i
    reactionstring(800+i)="(n,a_00)"
    write(reactionstring(800+i)(6:7),'(i2.2)') i
  enddo
  reactionstring(649)='(n,p_con)'
  reactionstring(699)='(n,d_con)'
  reactionstring(749)='(n,t_con)'
  reactionstring(799)='(n,h_con)'
  reactionstring(849)='(n,a_con)'
  do i = 1, 849
    if (reactionstring(i) /= '') reactionstring(i)(2:2) = parsym(k0)
  enddo
  return
end subroutine reactioncodes
! Copyright A.J. Koning 2023
