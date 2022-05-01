*DECK CUNHJ
      SUBROUTINE CUNHJ (Z, FNU, IPMTR, TOL, PHI, ARG, ZETA1, ZETA2,
     +   ASUM, BSUM)
C***BEGIN PROLOGUE  CUNHJ
C***SUBSIDIARY
C***PURPOSE  Subsidiary to CBESI and CBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNHJ-A, ZUNHJ-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         CUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***SEE ALSO  CBESI, CBESK
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CUNHJ
      COMPLEX ARG, ASUM, BSUM, CFNU, CONE, CR, CZERO, DR, P, PHI,
     * PRZTH, PTFN, RFN13, RTZTA, RZTH, SUMA, SUMB, TFN, T2, UP, W, W2,
     * Z, ZA, ZB, ZC, ZETA, ZETA1, ZETA2, ZTH
      REAL ALFA, ANG, AP, AR, ATOL, AW2, AZTH, BETA, BR, BTOL, C, EX1,
     * EX2, FNU, FN13, FN23, GAMA, HPI, PI, PP, RFNU, RFNU2, THPI, TOL,
     * WI, WR, ZCI, ZCR, ZETAI, ZETAR, ZTHI, ZTHR, ASUMR, ASUMI, BSUMR,
     * BSUMI, TEST, TSTR, TSTI, AC, R1MACH
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     * LRP1, L1, L2, M
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     * AP(30), P(30), UP(14), CR(14), DR(14)
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2     1.00000000000000000E+00,     1.04166666666666667E-01,
     3     8.35503472222222222E-02,     1.28226574556327160E-01,
     4     2.91849026464140464E-01,     8.81627267443757652E-01,
     5     3.32140828186276754E+00,     1.49957629868625547E+01,
     6     7.89230130115865181E+01,     4.74451538868264323E+02,
     7     3.20749009089066193E+03,     2.40865496408740049E+04,
     8     1.98923119169509794E+05,     1.79190200777534383E+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2     1.00000000000000000E+00,    -1.45833333333333333E-01,
     3    -9.87413194444444444E-02,    -1.43312053915895062E-01,
     4    -3.17227202678413548E-01,    -9.42429147957120249E-01,
     5    -3.51120304082635426E+00,    -1.57272636203680451E+01,
     6    -8.22814390971859444E+01,    -4.92355370523670524E+02,
     7    -3.31621856854797251E+03,    -2.48276742452085896E+04,
     8    -2.04526587315129788E+05,    -1.83844491706820990E+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000E+00,    -2.08333333333333333E-01,
     4     1.25000000000000000E-01,     3.34201388888888889E-01,
     5    -4.01041666666666667E-01,     7.03125000000000000E-02,
     6    -1.02581259645061728E+00,     1.84646267361111111E+00,
     7    -8.91210937500000000E-01,     7.32421875000000000E-02,
     8     4.66958442342624743E+00,    -1.12070026162229938E+01,
     9     8.78912353515625000E+00,    -2.36408691406250000E+00,
     A     1.12152099609375000E-01,    -2.82120725582002449E+01,
     B     8.46362176746007346E+01,    -9.18182415432400174E+01,
     C     4.25349987453884549E+01,    -7.36879435947963170E+00,
     D     2.27108001708984375E-01,     2.12570130039217123E+02,
     E    -7.65252468141181642E+02,     1.05999045252799988E+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541E+02,     2.18190511744211590E+02,
     4    -2.64914304869515555E+01,     5.72501420974731445E-01,
     5    -1.91945766231840700E+03,     8.06172218173730938E+03,
     6    -1.35865500064341374E+04,     1.16553933368645332E+04,
     7    -5.30564697861340311E+03,     1.20090291321635246E+03,
     8    -1.08090919788394656E+02,     1.72772750258445740E+00,
     9     2.02042913309661486E+04,    -9.69805983886375135E+04,
     A     1.92547001232531532E+05,    -2.03400177280415534E+05,
     B     1.22200464983017460E+05,    -4.11926549688975513E+04,
     C     7.10951430248936372E+03,    -4.93915304773088012E+02,
     D     6.07404200127348304E+00,    -2.42919187900551333E+05,
     E     1.31176361466297720E+06,    -2.99801591853810675E+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400E+06,    -2.81356322658653411E+06,
     4     1.26836527332162478E+06,    -3.31645172484563578E+05,
     5     4.52187689813627263E+04,    -2.49983048181120962E+03,
     6     2.43805296995560639E+01,     3.28446985307203782E+06,
     7    -1.97068191184322269E+07,     5.09526024926646422E+07,
     8    -7.41051482115326577E+07,     6.63445122747290267E+07,
     9    -3.75671766607633513E+07,     1.32887671664218183E+07,
     A    -2.78561812808645469E+06,     3.08186404612662398E+05,
     B    -1.38860897537170405E+04,     1.10017140269246738E+02,
     C    -4.93292536645099620E+07,     3.25573074185765749E+08,
     D    -9.39462359681578403E+08,     1.55359689957058006E+09,
     E    -1.62108055210833708E+09,     1.10684281682301447E+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309E+08,     1.42062907797533095E+08,
     4    -2.44740627257387285E+07,     2.24376817792244943E+06,
     5    -8.40054336030240853E+04,     5.51335896122020586E+02,
     6     8.14789096118312115E+08,    -5.86648149205184723E+09,
     7     1.86882075092958249E+10,    -3.46320433881587779E+10,
     8     4.12801855797539740E+10,    -3.30265997498007231E+10,
     9     1.79542137311556001E+10,    -6.56329379261928433E+09,
     A     1.55927986487925751E+09,    -2.25105661889415278E+08,
     B     1.73951075539781645E+07,    -5.49842327572288687E+05,
     C     3.03809051092238427E+03,    -1.46792612476956167E+10,
     D     1.14498237732025810E+11,    -3.99096175224466498E+11,
     E     8.19218669548577329E+11,    -1.09837515608122331E+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105)/
     2     1.00815810686538209E+12,    -6.45364869245376503E+11,
     3     2.87900649906150589E+11,    -8.78670721780232657E+10,
     4     1.76347306068349694E+10,    -2.16716498322379509E+09,
     5     1.43157876718888981E+08,    -3.87183344257261262E+06,
     6     1.82577554742931747E+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3     ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4    -4.44444444444444444E-03,    -9.22077922077922078E-04,
     5    -8.84892884892884893E-05,     1.65927687832449737E-04,
     6     2.46691372741792910E-04,     2.65995589346254780E-04,
     7     2.61824297061500945E-04,     2.48730437344655609E-04,
     8     2.32721040083232098E-04,     2.16362485712365082E-04,
     9     2.00738858762752355E-04,     1.86267636637545172E-04,
     A     1.73060775917876493E-04,     1.61091705929015752E-04,
     B     1.50274774160908134E-04,     1.40503497391269794E-04,
     C     1.31668816545922806E-04,     1.23667445598253261E-04,
     D     1.16405271474737902E-04,     1.09798298372713369E-04,
     E     1.03772410422992823E-04,     9.82626078369363448E-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3     ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4     9.32120517249503256E-05,     8.85710852478711718E-05,
     5     8.42963105715700223E-05,     8.03497548407791151E-05,
     6     7.66981345359207388E-05,     7.33122157481777809E-05,
     7     7.01662625163141333E-05,     6.72375633790160292E-05,
     8     6.93735541354588974E-04,     2.32241745182921654E-04,
     9    -1.41986273556691197E-05,    -1.16444931672048640E-04,
     A    -1.50803558053048762E-04,    -1.55121924918096223E-04,
     B    -1.46809756646465549E-04,    -1.33815503867491367E-04,
     C    -1.19744975684254051E-04,    -1.06184319207974020E-04,
     D    -9.37699549891194492E-05,    -8.26923045588193274E-05,
     E    -7.29374348155221211E-05,    -6.44042357721016283E-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3     ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4    -5.69611566009369048E-05,    -5.04731044303561628E-05,
     5    -4.48134868008882786E-05,    -3.98688727717598864E-05,
     6    -3.55400532972042498E-05,    -3.17414256609022480E-05,
     7    -2.83996793904174811E-05,    -2.54522720634870566E-05,
     8    -2.28459297164724555E-05,    -2.05352753106480604E-05,
     9    -1.84816217627666085E-05,    -1.66519330021393806E-05,
     A    -1.50179412980119482E-05,    -1.35554031379040526E-05,
     B    -1.22434746473858131E-05,    -1.10641884811308169E-05,
     C    -3.54211971457743841E-04,    -1.56161263945159416E-04,
     D     3.04465503594936410E-05,     1.30198655773242693E-04,
     E     1.67471106699712269E-04,     1.70222587683592569E-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3     ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4     1.56501427608594704E-04,     1.36339170977445120E-04,
     5     1.14886692029825128E-04,     9.45869093034688111E-05,
     6     7.64498419250898258E-05,     6.07570334965197354E-05,
     7     4.74394299290508799E-05,     3.62757512005344297E-05,
     8     2.69939714979224901E-05,     1.93210938247939253E-05,
     9     1.30056674793963203E-05,     7.82620866744496661E-06,
     A     3.59257485819351583E-06,     1.44040049814251817E-07,
     B    -2.65396769697939116E-06,    -4.91346867098485910E-06,
     C    -6.72739296091248287E-06,    -8.17269379678657923E-06,
     D    -9.31304715093561232E-06,    -1.02011418798016441E-05,
     E    -1.08805962510592880E-05,    -1.13875481509603555E-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4    -1.17519675674556414E-05,    -1.19987364870944141E-05,
     5     3.78194199201772914E-04,     2.02471952761816167E-04,
     6    -6.37938506318862408E-05,    -2.38598230603005903E-04,
     7    -3.10916256027361568E-04,    -3.13680115247576316E-04,
     8    -2.78950273791323387E-04,    -2.28564082619141374E-04,
     9    -1.75245280340846749E-04,    -1.25544063060690348E-04,
     A    -8.22982872820208365E-05,    -4.62860730588116458E-05,
     B    -1.72334302366962267E-05,     5.60690482304602267E-06,
     C     2.31395443148286800E-05,     3.62642745856793957E-05,
     D     4.58006124490188752E-05,     5.24595294959114050E-05,
     E     5.68396208545815266E-05,     5.94349820393104052E-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4     6.06478527578421742E-05,     6.08023907788436497E-05,
     5     6.01577894539460388E-05,     5.89199657344698500E-05,
     6     5.72515823777593053E-05,     5.52804375585852577E-05,
     7     5.31063773802880170E-05,     5.08069302012325706E-05,
     8     4.84418647620094842E-05,     4.60568581607475370E-05,
     9    -6.91141397288294174E-04,    -4.29976633058871912E-04,
     A     1.83067735980039018E-04,     6.60088147542014144E-04,
     B     8.75964969951185931E-04,     8.77335235958235514E-04,
     C     7.49369585378990637E-04,     5.63832329756980918E-04,
     D     3.68059319971443156E-04,     1.88464535514455599E-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4     3.70663057664904149E-05,    -8.28520220232137023E-05,
     5    -1.72751952869172998E-04,    -2.36314873605872983E-04,
     6    -2.77966150694906658E-04,    -3.02079514155456919E-04,
     7    -3.12594712643820127E-04,    -3.12872558758067163E-04,
     8    -3.05678038466324377E-04,    -2.93226470614557331E-04,
     9    -2.77255655582934777E-04,    -2.59103928467031709E-04,
     A    -2.39784014396480342E-04,    -2.20048260045422848E-04,
     B    -2.00443911094971498E-04,    -1.81358692210970687E-04,
     C    -1.63057674478657464E-04,    -1.45712672175205844E-04,
     D    -1.29425421983924587E-04,    -1.14245691942445952E-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4     1.92821964248775885E-03,     1.35592576302022234E-03,
     5    -7.17858090421302995E-04,    -2.58084802575270346E-03,
     6    -3.49271130826168475E-03,    -3.46986299340960628E-03,
     7    -2.82285233351310182E-03,    -1.88103076404891354E-03,
     8    -8.89531718383947600E-04,     3.87912102631035228E-06,
     9     7.28688540119691412E-04,     1.26566373053457758E-03,
     A     1.62518158372674427E-03,     1.83203153216373172E-03,
     B     1.91588388990527909E-03,     1.90588846755546138E-03,
     C     1.82798982421825727E-03,     1.70389506421121530E-03,
     D     1.55097127171097686E-03,     1.38261421852276159E-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2     1.20881424230064774E-03,     1.03676532638344962E-03,
     3     8.71437918068619115E-04,     7.16080155297701002E-04,
     4     5.72637002558129372E-04,     4.42089819465802277E-04,
     5     3.24724948503090564E-04,     2.20342042730246599E-04,
     6     1.28412898401353882E-04,     4.82005924552095464E-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3     BETA(19), BETA(20), BETA(21), BETA(22)/
     4     1.79988721413553309E-02,     5.59964911064388073E-03,
     5     2.88501402231132779E-03,     1.80096606761053941E-03,
     6     1.24753110589199202E-03,     9.22878876572938311E-04,
     7     7.14430421727287357E-04,     5.71787281789704872E-04,
     8     4.69431007606481533E-04,     3.93232835462916638E-04,
     9     3.34818889318297664E-04,     2.88952148495751517E-04,
     A     2.52211615549573284E-04,     2.22280580798883327E-04,
     B     1.97541838033062524E-04,     1.76836855019718004E-04,
     C     1.59316899661821081E-04,     1.44347930197333986E-04,
     D     1.31448068119965379E-04,     1.20245444949302884E-04,
     E     1.10449144504599392E-04,     1.01828770740567258E-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3     BETA(41), BETA(42), BETA(43), BETA(44)/
     4     9.41998224204237509E-05,     8.74130545753834437E-05,
     5     8.13466262162801467E-05,     7.59002269646219339E-05,
     6     7.09906300634153481E-05,     6.65482874842468183E-05,
     7     6.25146958969275078E-05,     5.88403394426251749E-05,
     8    -1.49282953213429172E-03,    -8.78204709546389328E-04,
     9    -5.02916549572034614E-04,    -2.94822138512746025E-04,
     A    -1.75463996970782828E-04,    -1.04008550460816434E-04,
     B    -5.96141953046457895E-05,    -3.12038929076098340E-05,
     C    -1.26089735980230047E-05,    -2.42892608575730389E-07,
     D     8.05996165414273571E-06,     1.36507009262147391E-05,
     E     1.73964125472926261E-05,     1.98672978842133780E-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3     BETA(63), BETA(64), BETA(65), BETA(66)/
     4     2.14463263790822639E-05,     2.23954659232456514E-05,
     5     2.28967783814712629E-05,     2.30785389811177817E-05,
     6     2.30321976080909144E-05,     2.28236073720348722E-05,
     7     2.25005881105292418E-05,     2.20981015361991429E-05,
     8     2.16418427448103905E-05,     2.11507649256220843E-05,
     9     2.06388749782170737E-05,     2.01165241997081666E-05,
     A     1.95913450141179244E-05,     1.90689367910436740E-05,
     B     1.85533719641636667E-05,     1.80475722259674218E-05,
     C     5.52213076721292790E-04,     4.47932581552384646E-04,
     D     2.79520653992020589E-04,     1.52468156198446602E-04,
     E     6.93271105657043598E-05,     1.76258683069991397E-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3     BETA(85), BETA(86), BETA(87), BETA(88)/
     4    -1.35744996343269136E-05,    -3.17972413350427135E-05,
     5    -4.18861861696693365E-05,    -4.69004889379141029E-05,
     6    -4.87665447413787352E-05,    -4.87010031186735069E-05,
     7    -4.74755620890086638E-05,    -4.55813058138628452E-05,
     8    -4.33309644511266036E-05,    -4.09230193157750364E-05,
     9    -3.84822638603221274E-05,    -3.60857167535410501E-05,
     A    -3.37793306123367417E-05,    -3.15888560772109621E-05,
     B    -2.95269561750807315E-05,    -2.75978914828335759E-05,
     C    -2.58006174666883713E-05,    -2.41308356761280200E-05,
     D    -2.25823509518346033E-05,    -2.11479656768912971E-05,
     E    -1.98200638885294927E-05,    -1.85909870801065077E-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4    -1.74532699844210224E-05,    -1.63997823854497997E-05,
     5    -4.74617796559959808E-04,    -4.77864567147321487E-04,
     6    -3.20390228067037603E-04,    -1.61105016119962282E-04,
     7    -4.25778101285435204E-05,     3.44571294294967503E-05,
     8     7.97092684075674924E-05,     1.03138236708272200E-04,
     9     1.12466775262204158E-04,     1.13103642108481389E-04,
     A     1.08651634848774268E-04,     1.01437951597661973E-04,
     B     9.29298396593363896E-05,     8.40293133016089978E-05,
     C     7.52727991349134062E-05,     6.69632521975730872E-05,
     D     5.92564547323194704E-05,     5.22169308826975567E-05,
     E     4.58539485165360646E-05,     4.01445513891486808E-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4     3.50481730031328081E-05,     3.05157995034346659E-05,
     5     2.64956119950516039E-05,     2.29363633690998152E-05,
     6     1.97893056664021636E-05,     1.70091984636412623E-05,
     7     1.45547428261524004E-05,     1.23886640995878413E-05,
     8     1.04775876076583236E-05,     8.79179954978479373E-06,
     9     7.36465810572578444E-04,     8.72790805146193976E-04,
     A     6.22614862573135066E-04,     2.85998154194304147E-04,
     B     3.84737672879366102E-06,    -1.87906003636971558E-04,
     C    -2.97603646594554535E-04,    -3.45998126832656348E-04,
     D    -3.53382470916037712E-04,    -3.35715635775048757E-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4    -3.04321124789039809E-04,    -2.66722723047612821E-04,
     5    -2.27654214122819527E-04,    -1.89922611854562356E-04,
     6    -1.55058918599093870E-04,    -1.23778240761873630E-04,
     7    -9.62926147717644187E-05,    -7.25178327714425337E-05,
     8    -5.22070028895633801E-05,    -3.50347750511900522E-05,
     9    -2.06489761035551757E-05,    -8.70106096849767054E-06,
     A     1.13698686675100290E-06,     9.16426474122778849E-06,
     B     1.56477785428872620E-05,     2.08223629482466847E-05,
     C     2.48923381004595156E-05,     2.80340509574146325E-05,
     D     3.03987774629861915E-05,     3.21156731406700616E-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4    -1.80182191963885708E-03,    -2.43402962938042533E-03,
     5    -1.83422663549856802E-03,    -7.62204596354009765E-04,
     6     2.39079475256927218E-04,     9.49266117176881141E-04,
     7     1.34467449701540359E-03,     1.48457495259449178E-03,
     8     1.44732339830617591E-03,     1.30268261285657186E-03,
     9     1.10351597375642682E-03,     8.86047440419791759E-04,
     A     6.73073208165665473E-04,     4.77603872856582378E-04,
     B     3.05991926358789362E-04,     1.60315694594721630E-04,
     C     4.00749555270613286E-05,    -5.66607461635251611E-05,
     D    -1.32506186772982638E-04,    -1.90296187989614057E-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4    -2.32811450376937408E-04,    -2.62628811464668841E-04,
     5    -2.82050469867598672E-04,    -2.93081563192861167E-04,
     6    -2.97435962176316616E-04,    -2.96557334239348078E-04,
     7    -2.91647363312090861E-04,    -2.83696203837734166E-04,
     8    -2.73512317095673346E-04,    -2.61750155806768580E-04,
     9     6.38585891212050914E-03,     9.62374215806377941E-03,
     A     7.61878061207001043E-03,     2.83219055545628054E-03,
     B    -2.09841352012720090E-03,    -5.73826764216626498E-03,
     C    -7.70804244495414620E-03,    -8.21011692264844401E-03,
     D    -7.65824520346905413E-03,    -6.47209729391045177E-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4    -4.99132412004966473E-03,    -3.45612289713133280E-03,
     5    -2.01785580014170775E-03,    -7.59430686781961401E-04,
     6     2.84173631523859138E-04,     1.10891667586337403E-03,
     7     1.72901493872728771E-03,     2.16812590802684701E-03,
     8     2.45357710494539735E-03,     2.61281821058334862E-03,
     9     2.67141039656276912E-03,     2.65203073395980430E-03,
     A     2.57411652877287315E-03,     2.45389126236094427E-03,
     B     2.30460058071795494E-03,     2.13684837686712662E-03,
     C     1.95896528478870911E-03,     1.77737008679454412E-03,
     D     1.59690280765839059E-03,     1.42111975664438546E-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3     GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4     6.29960524947436582E-01,     2.51984209978974633E-01,
     5     1.54790300415655846E-01,     1.10713062416159013E-01,
     6     8.57309395527394825E-02,     6.97161316958684292E-02,
     7     5.86085671893713576E-02,     5.04698873536310685E-02,
     8     4.42600580689154809E-02,     3.93720661543509966E-02,
     9     3.54283195924455368E-02,     3.21818857502098231E-02,
     A     2.94646240791157679E-02,     2.71581677112934479E-02,
     B     2.51768272973861779E-02,     2.34570755306078891E-02,
     C     2.19508390134907203E-02,     2.06210828235646240E-02,
     D     1.94388240897880846E-02,     1.83810633800683158E-02,
     E     1.74293213231963172E-02,     1.65685837786612353E-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1     GAMA(29), GAMA(30)/
     2     1.57865285987918445E-02,     1.50729501494095594E-02,
     3     1.44193250839954639E-02,     1.38184805735341786E-02,
     4     1.32643378994276568E-02,     1.27517121970498651E-02,
     5     1.22761545318762767E-02,     1.18338262398482403E-02/
      DATA EX1, EX2, HPI, PI, THPI /
     1     3.33333333333333333E-01,     6.66666666666666667E-01,
     2     1.57079632679489662E+00,     3.14159265358979324E+00,
     3     4.71238898038468986E+00/
      DATA CZERO, CONE / (0.0E0,0.0E0), (1.0E0,0.0E0) /
C***FIRST EXECUTABLE STATEMENT  CUNHJ
      RFNU = 1.0E0/FNU
C     ZB = Z*CMPLX(RFNU,0.0E0)
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TSTR = REAL(Z)
      TSTI = AIMAG(Z)
      TEST = R1MACH(1)*1.0E+3
      AC = FNU*TEST
      IF (ABS(TSTR).GT.AC .OR. ABS(TSTI).GT.AC) GO TO 15
      AC = 2.0E0*ABS(ALOG(TEST))+FNU
      ZETA1 = CMPLX(AC,0.0E0)
      ZETA2 = CMPLX(FNU,0.0E0)
      PHI=CONE
      ARG=CONE
      RETURN
   15 CONTINUE
      ZB = Z*CMPLX(RFNU,0.0E0)
      RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = CMPLX(1.0E0/FN13,0.0E0)
      W2 = CONE - ZB*ZB
      AW2 = ABS(W2)
      IF (AW2.GT.0.25E0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR ABS(W2).LE.0.25E0
C-----------------------------------------------------------------------
      K = 1
      P(1) = CONE
      SUMA = CMPLX(GAMA(1),0.0E0)
      AP(1) = 1.0E0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        P(K) = P(K-1)*W2
        SUMA = SUMA + P(K)*CMPLX(GAMA(K),0.0E0)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETA = W2*SUMA
      ARG = ZETA*CMPLX(FN23,0.0E0)
      ZA = CSQRT(SUMA)
      ZETA2 = CSQRT(W2)*CMPLX(FNU,0.0E0)
      ZETA1 = ZETA2*(CONE+ZETA*ZA*CMPLX(EX2,0.0E0))
      ZA = ZA + ZA
      PHI = CSQRT(ZA)*RFN13
      IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
      SUMB = CZERO
      DO 30 K=1,KMAX
        SUMB = SUMB + P(K)*CMPLX(BETA(K),0.0E0)
   30 CONTINUE
      ASUM = CZERO
      BSUM = SUMB
      L1 = 0
      L2 = 30
      BTOL = TOL*ABS(BSUM)
      ATOL = TOL
      PP = 1.0E0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMA = CZERO
        DO 40 K=1,KMAX
          M = L1 + K
          SUMA = SUMA + P(K)*CMPLX(ALFA(M),0.0E0)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUM = ASUM + SUMA*CMPLX(PP,0.0E0)
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMB = CZERO
        DO 70 K=1,KMAX
          M = L2 + K
          SUMB = SUMB + P(K)*CMPLX(BETA(M),0.0E0)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUM = BSUM + SUMB*CMPLX(PP,0.0E0)
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUM = ASUM + CONE
      PP = RFNU*REAL(RFN13)
      BSUM = BSUM*CMPLX(PP,0.0E0)
  120 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     ABS(W2).GT.0.25E0
C-----------------------------------------------------------------------
  130 CONTINUE
      W = CSQRT(W2)
      WR = REAL(W)
      WI = AIMAG(W)
      IF (WR.LT.0.0E0) WR = 0.0E0
      IF (WI.LT.0.0E0) WI = 0.0E0
      W = CMPLX(WR,WI)
      ZA = (CONE+W)/ZB
      ZC = CLOG(ZA)
      ZCR = REAL(ZC)
      ZCI = AIMAG(ZC)
      IF (ZCI.LT.0.0E0) ZCI = 0.0E0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0E0) ZCR = 0.0E0
      ZC = CMPLX(ZCR,ZCI)
      ZTH = (ZC-W)*CMPLX(1.5E0,0.0E0)
      CFNU = CMPLX(FNU,0.0E0)
      ZETA1 = ZC*CFNU
      ZETA2 = W*CFNU
      AZTH = ABS(ZTH)
      ZTHR = REAL(ZTH)
      ZTHI = AIMAG(ZTH)
      ANG = THPI
      IF (ZTHR.GE.0.0E0 .AND. ZTHI.LT.0.0E0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0E0) GO TO 140
      ANG = ATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0E0) ANG = ANG + PI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*COS(ANG)
      ZETAI = PP*SIN(ANG)
      IF (ZETAI.LT.0.0E0) ZETAI = 0.0E0
      ZETA = CMPLX(ZETAR,ZETAI)
      ARG = ZETA*CMPLX(FN23,0.0E0)
      RTZTA = ZTH/ZETA
      ZA = RTZTA/W
      PHI = CSQRT(ZA+ZA)*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      TFN = CMPLX(RFNU,0.0E0)/W
      RZTH = CMPLX(RFNU,0.0E0)/ZTH
      ZC = RZTH*CMPLX(AR(2),0.0E0)
      T2 = CONE/W2
      UP(2) = (T2*CMPLX(C(2),0.0E0)+CMPLX(C(3),0.0E0))*TFN
      BSUM = UP(2) + ZC
      ASUM = CZERO
      IF (RFNU.LT.TOL) GO TO 220
      PRZTH = RZTH
      PTFN = TFN
      UP(1) = CONE
      PP = 1.0E0
      BSUMR = REAL(BSUM)
      BSUMI = AIMAG(BSUM)
      BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZA = CMPLX(C(L),0.0E0)
          DO 150 J=2,KP1
            L = L + 1
            ZA = ZA*T2 + CMPLX(C(L),0.0E0)
  150     CONTINUE
          PTFN = PTFN*TFN
          UP(KP1) = PTFN*ZA
          CR(KS) = PRZTH*CMPLX(BR(KS+1),0.0E0)
          PRZTH = PRZTH*RZTH
          DR(KS) = PRZTH*CMPLX(AR(KS+2),0.0E0)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMA = UP(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMA = SUMA + CR(JR)*UP(JU)
  170   CONTINUE
        ASUM = ASUM + SUMA
        ASUMR = REAL(ASUM)
        ASUMI = AIMAG(ASUM)
        TEST = ABS(ASUMR) + ABS(ASUMI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMB = UP(LR+2) + UP(LRP1)*ZC
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMB = SUMB + DR(JR)*UP(JU)
  190   CONTINUE
        BSUM = BSUM + SUMB
        BSUMR = REAL(BSUM)
        BSUMI = AIMAG(BSUM)
        TEST = ABS(BSUMR) + ABS(BSUMI)
        IF (PP.LT.BTOL .AND. TEST.LT.TOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUM = ASUM + CONE
      BSUM = -BSUM*RFN13/RTZTA
      GO TO 120
      END