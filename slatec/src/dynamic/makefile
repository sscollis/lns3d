FC = gfortran
FFLAGS= -O2 -fPIC -dynamiclib -std=legacy -Wno-array-bounds
OBJ = \
aaaaaa.o  cnbdi.o   dcv.o     dqc25f.o  hpperm.o  qc25s.o       spoco.o \
acosh.o   cnbfa.o   ddaini.o  dqc25s.o  hpsort.o  qcheb.o       spodi.o \
ai.o      cnbfs.o   ddajac.o  dqcheb.o  hqr.o     qelg.o        spofa.o \
aie.o     cnbir.o   ddanrm.o  dqdota.o  hqr2.o    qform.o       spofs.o \
albeta.o  cnbsl.o   ddaslv.o  dqdoti.o  hstart.o  qk15.o        spoir.o \
algams.o  combak.o  ddassl.o  dqelg.o   hstcrt.o  qk15i.o       spopt.o \
ali.o     comhes.o  ddastp.o  dqform.o  hstcs1.o  qk15w.o       sposl.o \
alngam.o  comlr.o   ddatrp.o  dqk15.o   hstcsp.o  qk21.o        sppco.o \
alnrel.o  comlr2.o  ddaws.o   dqk15i.o  hstcyl.o  qk31.o        sppdi.o \
asinh.o   compb.o   ddawts.o  dqk15w.o  hstplr.o  qk41.o        spperm.o \
asyik.o   comqr.o   ddcor.o   dqk21.o   hstssp.o  qk51.o        sppfa.o \
asyjy.o   comqr2.o  ddcst.o   dqk31.o   htrib3.o  qk61.o        sppsl.o \
atanh.o   cortb.o   ddeabm.o  dqk41.o   htribk.o  qmomo.o       spsort.o \
avint.o   corth.o   ddebdf.o  dqk51.o   htrid3.o  qnc79.o       sptsl.o \
bakvec.o  cosdg.o   dderkf.o  dqk61.o   htridi.o  qng.o         sqrdc.o \
balanc.o  cosgen.o  ddes.o    dqmomo.o  hvnrm.o   qpdoc.o       sqrsl.o \
balbak.o  cosqb.o   ddntl.o   dqnc79.o  hw3crt.o  qpsrt.o       sreadp.o \
bandr.o   cosqb1.o  ddntp.o   dqng.o    hwscrt.o  qrfac.o       srlcal.o \
bandv.o   cosqf.o   ddoglg.o  dqpsrt.o  hwscs1.o  qrsolv.o      srot.o \
bcrh.o    cosqf1.o  ddot.o    dqrdc.o   hwscsp.o  qs2i1d.o      srotg.o \
bdiff.o   cosqi.o   ddpsc.o   dqrfac.o  hwscyl.o  qs2i1r.o      srotm.o \
besi.o    cost.o    ddpst.o   dqrsl.o   hwsplr.o  qwgtc.o       srotmg.o \
besi0.o   costi.o   ddriv1.o  dqrslv.o  hwsss1.o  qwgtf.o       ss2lt.o \
besi0e.o  cot.o     ddriv2.o  dqwgtc.o  hwsssp.o  qwgts.o       ss2y.o \
besi1.o   cpadd.o   ddriv3.o  dqwgtf.o  i1mach.o  qzhes.o       ssbmv.o \
besi1e.o  cpbco.o   ddscl.o   dqwgts.o  i1merg.o  qzit.o        sscal.o \
besj.o    cpbdi.o   ddstp.o   drc.o     icamax.o  qzval.o       ssd2s.o \
besj0.o   cpbfa.o   ddzro.o   drc3jj.o  icopy.o   qzvec.o       ssdbcg.o \
besj1.o   cpbsl.o   de1.o     drc3jm.o  idamax.o  r1mach.o      ssdcg.o \
besk.o    cpevl.o   deabm.o   drc6j.o   idloc.o   r1mpyq.o      ssdcgn.o \
besk0.o   cpevlr.o  debdf.o   drd.o     imtql1.o  r1updt.o      ssdcgs.o \
besk0e.o  cpoco.o   defc.o    dreadp.o  imtql2.o  r9aimp.o      ssdgmr.o \
besk1.o   cpodi.o   defcmn.o  dreort.o  imtqlv.o  r9atn1.o      ssdi.o \
besk1e.o  cpofa.o   defe4.o   drf.o     indxa.o   r9chu.o       ssdomn.o \
beskes.o  cpofs.o   defehl.o  drj.o     indxb.o   r9gmic.o      ssds.o \
besknu.o  cpoir.o   defer.o   indxc.o   r9gmit.o      ssdscl.o \
besks.o   cposl.o   dei.o     drkfs.o   initds.o  r9knus.o      ssgs.o \
besy.o    cppco.o   denorm.o  drlcal.o  inits.o   r9lgic.o      ssiccg.o \
besy0.o   cppdi.o   derf.o    drot.o    intrv.o   r9lgit.o      ssico.o \
besy1.o   cppfa.o   derfc.o   drotg.o   intyd.o   r9lgmc.o      ssics.o \
besynu.o  cppsl.o   derkf.o   drotm.o   invit.o   r9ln2r.o      ssidi.o \
beta.o    cpqr79.o  derkfs.o  drotmg.o  inxca.o   r9pak.o       ssiev.o \
betai.o   cproc.o   des.o     drsco.o   inxcb.o   r9upak.o      ssifa.o \
bfqad.o   cprocp.o  ds2lt.o   inxcc.o   radb2.o       ssilur.o \
bi.o      cprod.o   dexint.o  ds2y.o    iploc.o   radb3.o       ssilus.o \
bie.o     cprodp.o  dexprl.o  dsbmv.o   ipperm.o  radb4.o       ssisl.o \
binom.o   cpsi.o    dfac.o    dscal.o   ipsort.o  radb5.o       ssjac.o \
bint4.o   cptsl.o   dfc.o     dsd2s.o   isamax.o  radbg.o       ssli.o \
bintk.o   cpzero.o  dfcmn.o   dsdbcg.o  isdbcg.o  radf2.o       ssli2.o \
bisect.o  cqrdc.o   dfdjc1.o  dsdcg.o   isdcg.o   radf3.o       ssllti.o \
bkias.o   cqrsl.o   dfdjc3.o  dsdcgn.o  isdcgn.o  radf4.o       sslubc.o \
bkisr.o   crati.o   dfehl.o   dsdcgs.o  isdcgs.o  radf5.o       sslucn.o \
bksol.o   crotg.o   dfspvd.o  dsdgmr.o  isdgmr.o  radfg.o       sslucs.o \
blktr1.o  cs1s2.o   dfspvn.o  dsdi.o    isdir.o   rand.o        sslugm.o \
blktri.o  cscal.o   dfulmt.o  dsdomn.o  isdomn.o  ratqr.o       sslui.o \
bndacc.o  cscale.o  dfzero.o  dsdot.o   isort.o   rc.o          sslui2.o \
bndsol.o  cseri.o   dgami.o   dsds.o    issbcg.o  rc3jj.o       sslui4.o \
bnfac.o   csevl.o   dgamic.o  dsdscl.o  isscg.o   rc3jm.o       ssluom.o \
bnslv.o   cshch.o   dgamit.o  dsgs.o    isscgn.o  rc6j.o        ssluti.o \
bqr.o     csico.o   dgamlm.o  dsiccg.o  isscgs.o  rd.o          ssmmi2.o \
bsgq8.o   csidi.o   dgamln.o  dsico.o   issgmr.o  rebak.o       ssmmti.o \
bskin.o   csifa.o   dgamma.o  dsics.o   issir.o   rebakb.o      ssmtv.o \
bspdoc.o  csinh.o   dgamr.o   dsidi.o   issomn.o  reduc.o       ssmv.o \
bspdr.o   csisl.o   dgamrn.o  dsifa.o   iswap.o   reduc2.o      ssort.o \
bspev.o   cspco.o   dgaus8.o  dsilur.o  ivout.o   reort.o       sspco.o \
bsplvd.o  cspdi.o   dgbco.o   dsilus.o  j4save.o  rf.o          sspdi.o \
bsplvn.o  cspfa.o   dgbdi.o   dsindg.o  jairy.o   rfftb.o       sspev.o \
bsppp.o   cspsl.o   dgbfa.o   dsisl.o   la05ad.o  rfftb1.o      sspfa.o \
bspvd.o   csroot.o  dgbmv.o   dsjac.o   la05as.o  rfftf.o       sspmv.o \
bspvn.o   csrot.o   dgbsl.o   dsli.o    la05bd.o  rfftf1.o      sspr.o \
bsqad.o   csscal.o  dgeco.o   dsli2.o   la05bs.o  rffti.o       sspr2.o \
bsrh.o    csvdc.o   dgedi.o   dsllti.o  la05cd.o  rffti1.o      sspsl.o \
bvalu.o   cswap.o   dgefa.o   dslubc.o  la05cs.o  rg.o          ssvdc.o \
csymm.o   dgefs.o   dslucn.o  la05ed.o  rgauss.o      sswap.o \
csyr2k.o  dgemm.o   dslucs.o  la05es.o  rgg.o         ssymm.o \
csyrk.o   dgemv.o   dslugm.o  llsia.o   rj.o          ssymv.o \
c0lgmc.o  ctan.o    dger.o    dslui.o   lmpar.o   ssyr.o \
c1merg.o  ctanh.o   dgesl.o   dslui2.o  lpdp.o    rpqr79.o      ssyr2.o \
c9lgmc.o  ctbmv.o   dglss.o   dslui4.o  lsame.o   rpzero.o      ssyr2k.o \
c9ln2r.o  ctbsv.o   dgmres.o  dsluom.o  lsei.o    rs.o          ssyrk.o \
cacai.o   ctpmv.o   dgtsl.o   dsluti.o  lsi.o     rsb.o         stbmv.o \
cacon.o   ctpsv.o   dh12.o    dslvs.o   lsod.o    rsco.o        stbsv.o \
cacos.o   ctrco.o   dhels.o   dsmmi2.o  lssods.o  rsg.o         steps.o \
cacosh.o  ctrdi.o   dheqr.o   dsmmti.o  lssuds.o  rsgab.o       stin.o \
cairy.o   ctrmm.o   dhfti.o   dsmtv.o   macon.o   rsgba.o       stod.o \
carg.o    ctrmv.o   dhkseq.o  dsmv.o    rsp.o         stor1.o \
casin.o   ctrsl.o   dhstrt.o  dsort.o   mc20ad.o  rst.o         stout.o \
casinh.o  ctrsm.o   dhvnrm.o  dsos.o    mc20as.o  rt.o          stpmv.o \
casyi.o   ctrsv.o   dintp.o   dsoseq.o  mgsbv.o   runif.o       stpsv.o \
catan.o   cuchk.o   dintrv.o  dsossl.o  minfit.o  rwupdt.o      strco.o \
catan2.o  cunhj.o   dintyd.o  dspco.o   minso4.o  s1merg.o      strdi.o \
catanh.o  cuni1.o   dir.o     dspdi.o   minsol.o  sasum.o       strmm.o \
caxpy.o   cuni2.o   djairy.o  dspenc.o  mpadd.o   saxpy.o       strmv.o \
cbabk2.o  cunik.o   dlbeta.o  dspfa.o   mpadd2.o  sbcg.o        strsl.o \
cbal.o    cunk1.o   dlgams.o  dsplp.o   mpadd3.o  sbhin.o       strsm.o \
cbesh.o   cunk2.o   dli.o     dspmv.o   mpblas.o  sbocls.o      strsv.o \
cbesi.o   cuoik.o   dllsia.o  dspr.o    mpcdm.o   sbols.o       stway.o \
cbesj.o   cv.o      dllti2.o  dspr2.o   mpchk.o   sbolsm.o      suds.o \
cbesk.o   cwrsk.o   dlngam.o  dspsl.o   mpcmd.o   scasum.o      svco.o \
cbesy.o   d1mach.o  dlnrel.o  dsteps.o  mpdivi.o  scg.o         svd.o \
cbeta.o   d1merg.o  dlpdoc.o  dstod.o   mperr.o   scgn.o        svecs.o \
cbinu.o   d1mpyq.o  dlpdp.o   dstor1.o  mpmaxr.o  scgs.o        svout.o \
cbiry.o   d1updt.o  dlsei.o   dstway.o  mpmlp.o   schdc.o       swritp.o \
cbknu.o   d9aimp.o  dlsi.o    dsuds.o   mpmul.o   schdd.o       sxlcal.o \
cblkt1.o  d9atn1.o  dlsod.o   dsvco.o   mpmul2.o  schex.o       tevlc.o \
cblktr.o  d9b0mp.o  dlssud.o  dsvdc.o   mpmuli.o  schkw.o       tevls.o \
cbrt.o    d9b1mp.o  dmacon.o  dswap.o   mpnzr.o   schud.o       tinvit.o \
cbuni.o   d9chu.o   dmgsbv.o  dsymm.o   mpovfl.o  sclosm.o      tql1.o \
cbunk.o   d9gmic.o  dmout.o   dsymv.o   mpstr.o   scnrm2.o      tql2.o \
ccbrt.o   d9gmit.o  dmpar.o   dsyr.o    mpunfl.o  scoef.o       tqlrat.o \
cchdc.o   d9knus.o  dnbco.o   dsyr2.o   numxer.o  scopy.o       trbak1.o \
cchdd.o   d9lgic.o  dnbdi.o   dsyr2k.o  ohtrol.o  scopym.o      trbak3.o \
cchex.o   d9lgit.o  dnbfa.o   dsyrk.o   ohtror.o  scov.o        tred1.o \
cchud.o   d9lgmc.o  dnbfs.o   dtbmv.o   ortbak.o  scpplt.o      tred2.o \
ccmpb.o   d9ln2r.o  dnbsl.o   dtbsv.o   orthes.o  sdaini.o      tred3.o \
ccopy.o   d9pak.o   dnls1.o   dtin.o    ortho4.o  sdajac.o      tri3.o \
ccosh.o   d9upak.o  dnls1e.o  dtout.o   orthog.o  sdanrm.o      tridib.o \
ccot.o    dacosh.o  dnrm2.o   dtpmv.o   orthol.o  sdaslv.o      tridq.o \
cdcdot.o  dai.o     dnsq.o    dtpsv.o   orthor.o  sdassl.o      tris4.o \
cdcor.o   daie.o    dnsqe.o   dtrco.o   ortran.o  sdastp.o      trisp.o \
cdcst.o   dasinh.o  dogleg.o  dtrdi.o   passb.o   sdatrp.o      trix.o \
cdiv.o    dasum.o   dohtrl.o  dtrmm.o   passb2.o  sdawts.o      tsturm.o \
cdntl.o   dasyik.o  domn.o    dtrmv.o   passb3.o  sdcor.o       u11ls.o \
cdntp.o   dasyjy.o  dorth.o   dtrsl.o   passb4.o  sdcst.o       u11us.o \
cdotc.o   datanh.o  dorthr.o  dtrsm.o   passb5.o  sdntl.o       u12ls.o \
cdotu.o   davint.o  dp1vlu.o  dtrsv.o   passf.o   sdntp.o       u12us.o \
cdpsc.o   daws.o    dpbco.o   du11ls.o  passf2.o  sdot.o        ulsia.o \
cdpst.o   daxpy.o   dpbdi.o   du11us.o  passf3.o  sdpsc.o       usrmat.o \
cdriv1.o  dbcg.o    dpbfa.o   du12ls.o  passf4.o  sdpst.o       vnwrms.o \
cdriv2.o  dbdiff.o  dpbsl.o   du12us.o  passf5.o  sdriv1.o      wnlit.o \
cdriv3.o  dbesi.o   dpchbs.o  dulsia.o  pchbs.o   sdriv2.o      wnlsm.o \
cdscl.o   dbesi0.o  dpchce.o  dusrmt.o  pchce.o   sdriv3.o      wnlt1.o \
cdstp.o   dbesi1.o  dpchci.o  dvecs.o   pchci.o   sdscl.o       wnlt2.o \
cdzro.o   dbesj.o   dpchcm.o  dvnrms.o  pchcm.o   sdsdot.o      wnlt3.o \
cexprl.o  dbesj0.o  dpchcs.o  dvout.o   pchcs.o   sdstp.o       wnnls.o \
cfftb.o   dbesj1.o  dpchdf.o  dwnlit.o  pchdf.o   sdzro.o       xadd.o \
cfftb1.o  dbesk.o   dpchfd.o  dwnlsm.o  pchdoc.o  sepeli.o      xadj.o \
cfftf.o   dbesk0.o  dpchfe.o  dwnlt1.o  pchfd.o   sepx4.o       xc210.o \
cfftf1.o  dbesk1.o  dpchia.o  dwnlt2.o  pchfe.o   sgbco.o       xcon.o \
cffti.o   dbesks.o  dpchic.o  dwnlt3.o  pchia.o   sgbdi.o       xerbla.o \
cffti1.o  dbesy.o   dpchid.o  dwnnls.o  pchic.o   sgbfa.o       xerclr.o \
cfod.o    dbesy0.o  dpchim.o  dwritp.o  pchid.o   sgbmv.o       xercnt.o \
cg.o      dbesy1.o  dpchkt.o  dwupdt.o  pchim.o   sgbsl.o       xerdmp.o \
cgamma.o  dbeta.o   dpchng.o  dx.o      pchkt.o   sgeco.o       xerhlt.o \
cgamr.o   dbetai.o  dpchsp.o  dx4.o     pchngs.o  sgedi.o       xermax.o \
cgbco.o   dbfqad.o  dpchst.o  dxadd.o   pchsp.o   sgeev.o       xermsg.o \
cgbdi.o   dbhin.o   dpchsw.o  dxadj.o   pchst.o   sgefa.o       xerprn.o \
cgbfa.o   dbi.o     dpcoef.o  dxc210.o  pchsw.o   sgefs.o       xersve.o \
cgbmv.o   dbie.o    dpfqad.o  dxcon.o   pcoef.o   sgeir.o       xgetf.o \
cgbsl.o   dbinom.o  dpigmr.o  dxlcal.o  pfqad.o   sgemm.o       xgetua.o \
cgeco.o   dbint4.o  dpincw.o  dxlegf.o  pgsf.o    sgemv.o       xgetun.o \
cgedi.o   dbintk.o  dpinit.o  dxnrmp.o  pimach.o  sger.o        xlegf.o \
cgeev.o   dbkias.o  dpintm.o  dxpmu.o   pinitm.o  sgesl.o       xnrmp.o \
cgefa.o   dbkisr.o  dpjac.o   dxpmup.o  pjac.o    sglss.o       xpmu.o \
cgefs.o   dbksol.o  dplint.o  dxpnrm.o  pnnzrs.o  sgmres.o      xpmup.o \
cgeir.o   dbndac.o  dplpce.o  dxpqnu.o  poch.o    sgtsl.o       xpnrm.o \
cgemm.o   dbndsl.o  dplpdm.o  dxpsi.o   poch1.o   shels.o       xpqnu.o \
cgemv.o   dbnfac.o  dplpfe.o  dxqmu.o   pois3d.o  sheqr.o       xpsi.o \
cgerc.o   dbnslv.o  dplpfl.o  dxqnu.o   poisd2.o  sindg.o       xqmu.o \
cgeru.o   dbocls.o  dplpmn.o  dxred.o   poisn2.o  sinqb.o       xqnu.o \
cgesl.o   dbols.o   dplpmu.o  dxset.o   poisp2.o  sinqf.o       xred.o \
cgtsl.o   dbolsm.o  dplpup.o  dy.o      poistg.o  sinqi.o       xset.o \
ch.o      dbsgq8.o  dpnnzr.o  dy4.o     polcof.o  sint.o        xsetf.o \
chbmv.o   dbsi0e.o  dpoch.o   dyairy.o  polfit.o  sinti.o       xsetua.o \
chemm.o   dbsi1e.o  dpoch1.o  e1.o      polint.o  sintrp.o      xsetun.o \
chemv.o   dbsk0e.o  dpoco.o   efc.o     polyvl.o  sir.o         yairy.o \
cher.o    dbsk1e.o  dpodi.o   efcmn.o   pos3d1.o  zabs.o \
cher2.o   dbskes.o  dpofa.o   ei.o      postg2.o  sllti2.o      zacai.o \
cher2k.o  dbskin.o  dpofs.o   eisdoc.o  ppadd.o   slpdoc.o      zacon.o \
cherk.o   dbsknu.o  dpolcf.o  elmbak.o  ppgq8.o   slvs.o        zairy.o \
chfcm.o   dbspdr.o  dpolft.o  elmhes.o  ppgsf.o   smout.o       zasyi.o \
chfdv.o   dbspev.o  dpolvl.o  eltran.o  pppsf.o   snbco.o       zbesh.o \
chfev.o   dbsppp.o  dpopt.o   enorm.o   ppqad.o   snbdi.o       zbesi.o \
chfie.o   dbspvd.o  dposl.o   erf.o     ppsgf.o   snbfa.o       zbesj.o \
chico.o   dbspvn.o  dppco.o   erfc.o    ppspf.o   snbfs.o       zbesk.o \
chidi.o   dbsqad.o  dppdi.o   ppval.o   snbir.o       zbesy.o \
chiev.o   dbsynu.o  dpperm.o  exint.o   proc.o    snbsl.o       zbinu.o \
chifa.o   dbvalu.o  dppfa.o   exprel.o  procp.o   snls1.o       zbiry.o \
chisl.o   dppgq8.o  ezfft1.o  prod.o    snls1e.o      zbknu.o \
chkder.o  dppqad.o  ezfftb.o  prodp.o   snrm2.o       zbuni.o \
chkpr4.o  dppsl.o   ezfftf.o  prvec.o   snsq.o        zbunk.o \
chkprm.o  dcbrt.o   dppval.o  ezffti.o  prwpge.o  snsqe.o       zdiv.o \
chksn4.o  dcdot.o   dprvec.o  fac.o     prwvir.o  sods.o        zexp.o \
chksng.o  dcfod.o   dprwpg.o  fc.o      psgf.o    somn.o        zkscl.o \
chpco.o   dcg.o     dprwvr.o  fcmn.o    psi.o     sopenm.o      zlog.o \
chpdi.o   dcgn.o    dpsi.o    fdjac1.o  psifn.o   sorth.o       zmlri.o \
chpfa.o   dcgs.o    dpsifn.o  fdjac3.o  psixn.o   sos.o         zmlt.o \
chpmv.o   dchdc.o   dpsixn.o  fdump.o   pvalue.o  soseqs.o      zrati.o \
chpr.o    dchdd.o   dpsort.o  fftdoc.o  pythag.o  sossol.o      zs1s2.o \
chpr2.o   dchex.o   dptsl.o   figi.o    qag.o     spbco.o       zseri.o \
chpsl.o   dchfcm.o  dqag.o    figi2.o   qage.o    spbdi.o       zshch.o \
chu.o     dchfdv.o  dqage.o   fulmat.o  qagi.o    spbfa.o       zsqrt.o \
cinvit.o  dchfev.o  dqagi.o   fundoc.o  qagie.o   spbsl.o       zuchk.o \
ckscl.o   dchfie.o  dqagie.o  fzero.o   qagp.o    speli4.o      zunhj.o \
clbeta.o  dchkw.o   dqagp.o   gami.o    qagpe.o   spelip.o      zuni1.o \
clngam.o  dchu.o    dqagpe.o  gamic.o   qags.o    spenc.o       zuni2.o \
clnrel.o  dchud.o   dqags.o   gamit.o   qagse.o   spigmr.o      zunik.o \
clog10.o  dckder.o  dqagse.o  gamlim.o  qawc.o    spincw.o      zunk1.o \
cmgnbn.o  dcoef.o   dqawc.o   gamln.o   qawce.o   spinit.o      zunk2.o \
cmlri.o   dcopy.o   dqawce.o  gamma.o   qawf.o    splp.o        zuoik.o \
cmpcsg.o  dcopym.o  dqawf.o   gamr.o    qawfe.o   splpce.o      zwrsk.o \
cmposd.o  dcosdg.o  dqawfe.o  gamrn.o   qawo.o    splpdm.o \
cmposn.o  dcot.o    dqawo.o   gaus8.o   qawoe.o   splpfe.o \
cmposp.o  dcov.o    dqawoe.o  genbun.o  qaws.o    splpfl.o \
cmptr3.o  dcpplt.o  dqaws.o   h12.o     qawse.o   splpmn.o \
cmptrx.o  dcscal.o  dqawse.o  hfti.o    qc25c.o   splpmu.o \
cnbco.o   dcsevl.o  dqc25c.o  hkseq.o   qc25f.o   splpup.o

#
#  Dynamic library target
#

libslatec.so: $(OBJ)
	@echo $(FC)
	$(FC) -dynamiclib -Xlinker -rpath -Xlinker ../../lib -o $@ $(OBJ)

${OBJ}: %.o: ../%.f
	${FC} ${FFLAGS} -c -o $@ $<

clean:
	rm -f *.o libslatec.so
