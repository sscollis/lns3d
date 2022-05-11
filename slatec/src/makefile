#
# Set the installation destination for SLATEC to be up one level
#
FC = gfortan

ifndef DEST
  DEST = ..
endif

MAN = \
aaaaaa.1 cnbdi.1  dcv.1    dqc25f.1 hpperm.1 qc25s.1      spoco.1 \
acosh.1  cnbfa.1  ddaini.1 dqc25s.1 hpsort.1 qcheb.1      spodi.1 \
ai.1     cnbfs.1  ddajac.1 dqcheb.1 hqr.1    qelg.1       spofa.1 \
aie.1    cnbir.1  ddanrm.1 dqdota.1 hqr2.1   qform.1      spofs.1 \
albeta.1 cnbsl.1  ddaslv.1 dqdoti.1 hstart.1 qk15.1       spoir.1 \
algams.1 combak.1 ddassl.1 dqelg.1  hstcrt.1 qk15i.1      spopt.1 \
ali.1    comhes.1 ddastp.1 dqform.1 hstcs1.1 qk15w.1      sposl.1 \
alngam.1 comlr.1  ddatrp.1 dqk15.1  hstcsp.1 qk21.1       sppco.1 \
alnrel.1 comlr2.1 ddaws.1  dqk15i.1 hstcyl.1 qk31.1       sppdi.1 \
asinh.1  compb.1  ddawts.1 dqk15w.1 hstplr.1 qk41.1       spperm.1 \
asyik.1  comqr.1  ddcor.1  dqk21.1  hstssp.1 qk51.1       sppfa.1 \
asyjy.1  comqr2.1 ddcst.1  dqk31.1  htrib3.1 qk61.1       sppsl.1 \
atanh.1  cortb.1  ddeabm.1 dqk41.1  htribk.1 qmomo.1      spsort.1 \
avint.1  corth.1  ddebdf.1 dqk51.1  htrid3.1 qnc79.1      sptsl.1 \
bakvec.1 cosdg.1  dderkf.1 dqk61.1  htridi.1 qng.1        sqrdc.1 \
balanc.1 cosgen.1 ddes.1   dqmomo.1 hvnrm.1  qpdoc.1      sqrsl.1 \
balbak.1 cosqb.1  ddntl.1  dqnc79.1 hw3crt.1 qpsrt.1      sreadp.1 \
bandr.1  cosqb1.1 ddntp.1  dqng.1   hwscrt.1 qrfac.1      srlcal.1 \
bandv.1  cosqf.1  ddoglg.1 dqpsrt.1 hwscs1.1 qrsolv.1     srot.1 \
bcrh.1   cosqf1.1 ddot.1   dqrdc.1  hwscsp.1 qs2i1d.1     srotg.1 \
bdiff.1  cosqi.1  ddpsc.1  dqrfac.1 hwscyl.1 qs2i1r.1     srotm.1 \
besi.1   cost.1   ddpst.1  dqrsl.1  hwsplr.1 qwgtc.1      srotmg.1 \
besi0.1  costi.1  ddriv1.1 dqrslv.1 hwsss1.1 qwgtf.1      ss2lt.1 \
besi0e.1 cot.1    ddriv2.1 dqwgtc.1 hwsssp.1 qwgts.1      ss2y.1 \
besi1.1  cpadd.1  ddriv3.1 dqwgtf.1 i1mach.1 qzhes.1      ssbmv.1 \
besi1e.1 cpbco.1  ddscl.1  dqwgts.1 i1merg.1 qzit.1       sscal.1 \
besj.1   cpbdi.1  ddstp.1  drc.1    icamax.1 qzval.1      ssd2s.1 \
besj0.1  cpbfa.1  ddzro.1  drc3jj.1 icopy.1  qzvec.1      ssdbcg.1 \
besj1.1  cpbsl.1  de1.1    drc3jm.1 idamax.1 r1mach.1     ssdcg.1 \
besk.1   cpevl.1  deabm.1  drc6j.1  idloc.1  r1mpyq.1     ssdcgn.1 \
besk0.1  cpevlr.1 debdf.1  drd.1    imtql1.1 r1updt.1     ssdcgs.1 \
besk0e.1 cpoco.1  defc.1   dreadp.1 imtql2.1 r9aimp.1     ssdgmr.1 \
besk1.1  cpodi.1  defcmn.1 dreort.1 imtqlv.1 r9atn1.1     ssdi.1 \
besk1e.1 cpofa.1  defe4.1  drf.1    indxa.1  r9chu.1      ssdomn.1 \
beskes.1 cpofs.1  defehl.1 drj.1    indxb.1  r9gmic.1     ssds.1 \
besknu.1 cpoir.1  defer.1  indxc.1  r9gmit.1     ssdscl.1 \
besks.1  cposl.1  dei.1    drkfs.1  initds.1 r9knus.1     ssgs.1 \
besy.1   cppco.1  denorm.1 drlcal.1 inits.1  r9lgic.1     ssiccg.1 \
besy0.1  cppdi.1  derf.1   drot.1   intrv.1  r9lgit.1     ssico.1 \
besy1.1  cppfa.1  derfc.1  drotg.1  intyd.1  r9lgmc.1     ssics.1 \
besynu.1 cppsl.1  derkf.1  drotm.1  invit.1  r9ln2r.1     ssidi.1 \
beta.1   cpqr79.1 derkfs.1 drotmg.1 inxca.1  r9pak.1      ssiev.1 \
betai.1  cproc.1  des.1    drsco.1  inxcb.1  r9upak.1     ssifa.1 \
bfqad.1  cprocp.1 ds2lt.1  inxcc.1  radb2.1      ssilur.1 \
bi.1     cprod.1  dexint.1 ds2y.1   iploc.1  radb3.1      ssilus.1 \
bie.1    cprodp.1 dexprl.1 dsbmv.1  ipperm.1 radb4.1      ssisl.1 \
binom.1  cpsi.1   dfac.1   dscal.1  ipsort.1 radb5.1      ssjac.1 \
bint4.1  cptsl.1  dfc.1    dsd2s.1  isamax.1 radbg.1      ssli.1 \
bintk.1  cpzero.1 dfcmn.1  dsdbcg.1 isdbcg.1 radf2.1      ssli2.1 \
bisect.1 cqrdc.1  dfdjc1.1 dsdcg.1  isdcg.1  radf3.1      ssllti.1 \
bkias.1  cqrsl.1  dfdjc3.1 dsdcgn.1 isdcgn.1 radf4.1      sslubc.1 \
bkisr.1  crati.1  dfehl.1  dsdcgs.1 isdcgs.1 radf5.1      sslucn.1 \
bksol.1  crotg.1  dfspvd.1 dsdgmr.1 isdgmr.1 radfg.1      sslucs.1 \
blktr1.1 cs1s2.1  dfspvn.1 dsdi.1   isdir.1  rand.1       sslugm.1 \
blktri.1 cscal.1  dfulmt.1 dsdomn.1 isdomn.1 ratqr.1      sslui.1 \
bndacc.1 cscale.1 dfzero.1 dsdot.1  isort.1  rc.1         sslui2.1 \
bndsol.1 cseri.1  dgami.1  dsds.1   issbcg.1 rc3jj.1      sslui4.1 \
bnfac.1  csevl.1  dgamic.1 dsdscl.1 isscg.1  rc3jm.1      ssluom.1 \
bnslv.1  cshch.1  dgamit.1 dsgs.1   isscgn.1 rc6j.1       ssluti.1 \
bqr.1    csico.1  dgamlm.1 dsiccg.1 isscgs.1 rd.1         ssmmi2.1 \
bsgq8.1  csidi.1  dgamln.1 dsico.1  issgmr.1 rebak.1      ssmmti.1 \
bskin.1  csifa.1  dgamma.1 dsics.1  issir.1  rebakb.1     ssmtv.1 \
bspdoc.1 csinh.1  dgamr.1  dsidi.1  issomn.1 reduc.1      ssmv.1 \
bspdr.1  csisl.1  dgamrn.1 dsifa.1  iswap.1  reduc2.1     ssort.1 \
bspev.1  cspco.1  dgaus8.1 dsilur.1 ivout.1  reort.1      sspco.1 \
bsplvd.1 cspdi.1  dgbco.1  dsilus.1 j4save.1 rf.1         sspdi.1 \
bsplvn.1 cspfa.1  dgbdi.1  dsindg.1 jairy.1  rfftb.1      sspev.1 \
bsppp.1  cspsl.1  dgbfa.1  dsisl.1  la05ad.1 rfftb1.1     sspfa.1 \
bspvd.1  csroot.1 dgbmv.1  dsjac.1  la05as.1 rfftf.1      sspmv.1 \
bspvn.1  csrot.1  dgbsl.1  dsli.1   la05bd.1 rfftf1.1     sspr.1 \
bsqad.1  csscal.1 dgeco.1  dsli2.1  la05bs.1 rffti.1      sspr2.1 \
bsrh.1   csvdc.1  dgedi.1  dsllti.1 la05cd.1 rffti1.1     sspsl.1 \
bvalu.1  cswap.1  dgefa.1  dslubc.1 la05cs.1 rg.1         ssvdc.1 \
csymm.1  dgefs.1  dslucn.1 la05ed.1 rgauss.1     sswap.1 \
csyr2k.1 dgemm.1  dslucs.1 la05es.1 rgg.1        ssymm.1 \
csyrk.1  dgemv.1  dslugm.1 llsia.1  rj.1         ssymv.1 \
c0lgmc.1 ctan.1   dger.1   dslui.1  lmpar.1  ssyr.1 \
c1merg.1 ctanh.1  dgesl.1  dslui2.1 lpdp.1   rpqr79.1     ssyr2.1 \
c9lgmc.1 ctbmv.1  dglss.1  dslui4.1 lsame.1  rpzero.1     ssyr2k.1 \
c9ln2r.1 ctbsv.1  dgmres.1 dsluom.1 lsei.1   rs.1         ssyrk.1 \
cacai.1  ctpmv.1  dgtsl.1  dsluti.1 lsi.1    rsb.1        stbmv.1 \
cacon.1  ctpsv.1  dh12.1   dslvs.1  lsod.1   rsco.1       stbsv.1 \
cacos.1  ctrco.1  dhels.1  dsmmi2.1 lssods.1 rsg.1        steps.1 \
cacosh.1 ctrdi.1  dheqr.1  dsmmti.1 lssuds.1 rsgab.1      stin.1 \
cairy.1  ctrmm.1  dhfti.1  dsmtv.1  macon.1  rsgba.1      stod.1 \
carg.1   ctrmv.1  dhkseq.1 dsmv.1   rsp.1        stor1.1 \
casin.1  ctrsl.1  dhstrt.1 dsort.1  mc20ad.1 rst.1        stout.1 \
casinh.1 ctrsm.1  dhvnrm.1 dsos.1   mc20as.1 rt.1         stpmv.1 \
casyi.1  ctrsv.1  dintp.1  dsoseq.1 mgsbv.1  runif.1      stpsv.1 \
catan.1  cuchk.1  dintrv.1 dsossl.1 minfit.1 rwupdt.1     strco.1 \
catan2.1 cunhj.1  dintyd.1 dspco.1  minso4.1 s1merg.1     strdi.1 \
catanh.1 cuni1.1  dir.1    dspdi.1  minsol.1 sasum.1      strmm.1 \
caxpy.1  cuni2.1  djairy.1 dspenc.1 mpadd.1  saxpy.1      strmv.1 \
cbabk2.1 cunik.1  dlbeta.1 dspfa.1  mpadd2.1 sbcg.1       strsl.1 \
cbal.1   cunk1.1  dlgams.1 dsplp.1  mpadd3.1 sbhin.1      strsm.1 \
cbesh.1  cunk2.1  dli.1    dspmv.1  mpblas.1 sbocls.1     strsv.1 \
cbesi.1  cuoik.1  dllsia.1 dspr.1   mpcdm.1  sbols.1      stway.1 \
cbesj.1  cv.1     dllti2.1 dspr2.1  mpchk.1  sbolsm.1     suds.1 \
cbesk.1  cwrsk.1  dlngam.1 dspsl.1  mpcmd.1  scasum.1     svco.1 \
cbesy.1  d1mach.1 dlnrel.1 dsteps.1 mpdivi.1 scg.1        svd.1 \
cbeta.1  d1merg.1 dlpdoc.1 dstod.1  mperr.1  scgn.1       svecs.1 \
cbinu.1  d1mpyq.1 dlpdp.1  dstor1.1 mpmaxr.1 scgs.1       svout.1 \
cbiry.1  d1updt.1 dlsei.1  dstway.1 mpmlp.1  schdc.1      swritp.1 \
cbknu.1  d9aimp.1 dlsi.1   dsuds.1  mpmul.1  schdd.1      sxlcal.1 \
cblkt1.1 d9atn1.1 dlsod.1  dsvco.1  mpmul2.1 schex.1      tevlc.1 \
cblktr.1 d9b0mp.1 dlssud.1 dsvdc.1  mpmuli.1 schkw.1      tevls.1 \
cbrt.1   d9b1mp.1 dmacon.1 dswap.1  mpnzr.1  schud.1      tinvit.1 \
cbuni.1  d9chu.1  dmgsbv.1 dsymm.1  mpovfl.1 sclosm.1     tql1.1 \
cbunk.1  d9gmic.1 dmout.1  dsymv.1  mpstr.1  scnrm2.1     tql2.1 \
ccbrt.1  d9gmit.1 dmpar.1  dsyr.1   mpunfl.1 scoef.1      tqlrat.1 \
cchdc.1  d9knus.1 dnbco.1  dsyr2.1  numxer.1 scopy.1      trbak1.1 \
cchdd.1  d9lgic.1 dnbdi.1  dsyr2k.1 ohtrol.1 scopym.1     trbak3.1 \
cchex.1  d9lgit.1 dnbfa.1  dsyrk.1  ohtror.1 scov.1       tred1.1 \
cchud.1  d9lgmc.1 dnbfs.1  dtbmv.1  ortbak.1 scpplt.1     tred2.1 \
ccmpb.1  d9ln2r.1 dnbsl.1  dtbsv.1  orthes.1 sdaini.1     tred3.1 \
ccopy.1  d9pak.1  dnls1.1  dtin.1   ortho4.1 sdajac.1     tri3.1 \
ccosh.1  d9upak.1 dnls1e.1 dtout.1  orthog.1 sdanrm.1     tridib.1 \
ccot.1   dacosh.1 dnrm2.1  dtpmv.1  orthol.1 sdaslv.1     tridq.1 \
cdcdot.1 dai.1    dnsq.1   dtpsv.1  orthor.1 sdassl.1     tris4.1 \
cdcor.1  daie.1   dnsqe.1  dtrco.1  ortran.1 sdastp.1     trisp.1 \
cdcst.1  dasinh.1 dogleg.1 dtrdi.1  passb.1  sdatrp.1     trix.1 \
cdiv.1   dasum.1  dohtrl.1 dtrmm.1  passb2.1 sdawts.1     tsturm.1 \
cdntl.1  dasyik.1 domn.1   dtrmv.1  passb3.1 sdcor.1      u11ls.1 \
cdntp.1  dasyjy.1 dorth.1  dtrsl.1  passb4.1 sdcst.1      u11us.1 \
cdotc.1  datanh.1 dorthr.1 dtrsm.1  passb5.1 sdntl.1      u12ls.1 \
cdotu.1  davint.1 dp1vlu.1 dtrsv.1  passf.1  sdntp.1      u12us.1 \
cdpsc.1  daws.1   dpbco.1  du11ls.1 passf2.1 sdot.1       ulsia.1 \
cdpst.1  daxpy.1  dpbdi.1  du11us.1 passf3.1 sdpsc.1      usrmat.1 \
cdriv1.1 dbcg.1   dpbfa.1  du12ls.1 passf4.1 sdpst.1      vnwrms.1 \
cdriv2.1 dbdiff.1 dpbsl.1  du12us.1 passf5.1 sdriv1.1     wnlit.1 \
cdriv3.1 dbesi.1  dpchbs.1 dulsia.1 pchbs.1  sdriv2.1     wnlsm.1 \
cdscl.1  dbesi0.1 dpchce.1 dusrmt.1 pchce.1  sdriv3.1     wnlt1.1 \
cdstp.1  dbesi1.1 dpchci.1 dvecs.1  pchci.1  sdscl.1      wnlt2.1 \
cdzro.1  dbesj.1  dpchcm.1 dvnrms.1 pchcm.1  sdsdot.1     wnlt3.1 \
cexprl.1 dbesj0.1 dpchcs.1 dvout.1  pchcs.1  sdstp.1      wnnls.1 \
cfftb.1  dbesj1.1 dpchdf.1 dwnlit.1 pchdf.1  sdzro.1      xadd.1 \
cfftb1.1 dbesk.1  dpchfd.1 dwnlsm.1 pchdoc.1 sepeli.1     xadj.1 \
cfftf.1  dbesk0.1 dpchfe.1 dwnlt1.1 pchfd.1  sepx4.1      xc210.1 \
cfftf1.1 dbesk1.1 dpchia.1 dwnlt2.1 pchfe.1  sgbco.1      xcon.1 \
cffti.1  dbesks.1 dpchic.1 dwnlt3.1 pchia.1  sgbdi.1      xerbla.1 \
cffti1.1 dbesy.1  dpchid.1 dwnnls.1 pchic.1  sgbfa.1      xerclr.1 \
cfod.1   dbesy0.1 dpchim.1 dwritp.1 pchid.1  sgbmv.1      xercnt.1 \
cg.1     dbesy1.1 dpchkt.1 dwupdt.1 pchim.1  sgbsl.1      xerdmp.1 \
cgamma.1 dbeta.1  dpchng.1 dx.1     pchkt.1  sgeco.1      xerhlt.1 \
cgamr.1  dbetai.1 dpchsp.1 dx4.1    pchngs.1 sgedi.1      xermax.1 \
cgbco.1  dbfqad.1 dpchst.1 dxadd.1  pchsp.1  sgeev.1      xermsg.1 \
cgbdi.1  dbhin.1  dpchsw.1 dxadj.1  pchst.1  sgefa.1      xerprn.1 \
cgbfa.1  dbi.1    dpcoef.1 dxc210.1 pchsw.1  sgefs.1      xersve.1 \
cgbmv.1  dbie.1   dpfqad.1 dxcon.1  pcoef.1  sgeir.1      xgetf.1 \
cgbsl.1  dbinom.1 dpigmr.1 dxlcal.1 pfqad.1  sgemm.1      xgetua.1 \
cgeco.1  dbint4.1 dpincw.1 dxlegf.1 pgsf.1   sgemv.1      xgetun.1 \
cgedi.1  dbintk.1 dpinit.1 dxnrmp.1 pimach.1 sger.1       xlegf.1 \
cgeev.1  dbkias.1 dpintm.1 dxpmu.1  pinitm.1 sgesl.1      xnrmp.1 \
cgefa.1  dbkisr.1 dpjac.1  dxpmup.1 pjac.1   sglss.1      xpmu.1 \
cgefs.1  dbksol.1 dplint.1 dxpnrm.1 pnnzrs.1 sgmres.1     xpmup.1 \
cgeir.1  dbndac.1 dplpce.1 dxpqnu.1 poch.1   sgtsl.1      xpnrm.1 \
cgemm.1  dbndsl.1 dplpdm.1 dxpsi.1  poch1.1  shels.1      xpqnu.1 \
cgemv.1  dbnfac.1 dplpfe.1 dxqmu.1  pois3d.1 sheqr.1      xpsi.1 \
cgerc.1  dbnslv.1 dplpfl.1 dxqnu.1  poisd2.1 sindg.1      xqmu.1 \
cgeru.1  dbocls.1 dplpmn.1 dxred.1  poisn2.1 sinqb.1      xqnu.1 \
cgesl.1  dbols.1  dplpmu.1 dxset.1  poisp2.1 sinqf.1      xred.1 \
cgtsl.1  dbolsm.1 dplpup.1 dy.1     poistg.1 sinqi.1      xset.1 \
ch.1     dbsgq8.1 dpnnzr.1 dy4.1    polcof.1 sint.1       xsetf.1 \
chbmv.1  dbsi0e.1 dpoch.1  dyairy.1 polfit.1 sinti.1      xsetua.1 \
chemm.1  dbsi1e.1 dpoch1.1 e1.1     polint.1 sintrp.1     xsetun.1 \
chemv.1  dbsk0e.1 dpoco.1  efc.1    polyvl.1 sir.1        yairy.1 \
cher.1   dbsk1e.1 dpodi.1  efcmn.1  pos3d1.1 zabs.1 \
cher2.1  dbskes.1 dpofa.1  ei.1     postg2.1 sllti2.1     zacai.1 \
cher2k.1 dbskin.1 dpofs.1  eisdoc.1 ppadd.1  slpdoc.1     zacon.1 \
cherk.1  dbsknu.1 dpolcf.1 elmbak.1 ppgq8.1  slvs.1       zairy.1 \
chfcm.1  dbspdr.1 dpolft.1 elmhes.1 ppgsf.1  smout.1      zasyi.1 \
chfdv.1  dbspev.1 dpolvl.1 eltran.1 pppsf.1  snbco.1      zbesh.1 \
chfev.1  dbsppp.1 dpopt.1  enorm.1  ppqad.1  snbdi.1      zbesi.1 \
chfie.1  dbspvd.1 dposl.1  erf.1    ppsgf.1  snbfa.1      zbesj.1 \
chico.1  dbspvn.1 dppco.1  erfc.1   ppspf.1  snbfs.1      zbesk.1 \
chidi.1  dbsqad.1 dppdi.1  ppval.1  snbir.1      zbesy.1 \
chiev.1  dbsynu.1 dpperm.1 exint.1  proc.1   snbsl.1      zbinu.1 \
chifa.1  dbvalu.1 dppfa.1  exprel.1 procp.1  snls1.1      zbiry.1 \
chisl.1  dppgq8.1 ezfft1.1 prod.1   snls1e.1     zbknu.1 \
chkder.1 dppqad.1 ezfftb.1 prodp.1  snrm2.1      zbuni.1 \
chkpr4.1 dppsl.1  ezfftf.1 prvec.1  snsq.1       zbunk.1 \
chkprm.1 dcbrt.1  dppval.1 ezffti.1 prwpge.1 snsqe.1      zdiv.1 \
chksn4.1 dcdot.1  dprvec.1 fac.1    prwvir.1 sods.1       zexp.1 \
chksng.1 dcfod.1  dprwpg.1 fc.1     psgf.1   somn.1       zkscl.1 \
chpco.1  dcg.1    dprwvr.1 fcmn.1   psi.1    sopenm.1     zlog.1 \
chpdi.1  dcgn.1   dpsi.1   fdjac1.1 psifn.1  sorth.1      zmlri.1 \
chpfa.1  dcgs.1   dpsifn.1 fdjac3.1 psixn.1  sos.1        zmlt.1 \
chpmv.1  dchdc.1  dpsixn.1 fdump.1  pvalue.1 soseqs.1     zrati.1 \
chpr.1   dchdd.1  dpsort.1 fftdoc.1 pythag.1 sossol.1     zs1s2.1 \
chpr2.1  dchex.1  dptsl.1  figi.1   qag.1    spbco.1      zseri.1 \
chpsl.1  dchfcm.1 dqag.1   figi2.1  qage.1   spbdi.1      zshch.1 \
chu.1    dchfdv.1 dqage.1  fulmat.1 qagi.1   spbfa.1      zsqrt.1 \
cinvit.1 dchfev.1 dqagi.1  fundoc.1 qagie.1  spbsl.1      zuchk.1 \
ckscl.1  dchfie.1 dqagie.1 fzero.1  qagp.1   speli4.1     zunhj.1 \
clbeta.1 dchkw.1  dqagp.1  gami.1   qagpe.1  spelip.1     zuni1.1 \
clngam.1 dchu.1   dqagpe.1 gamic.1  qags.1   spenc.1      zuni2.1 \
clnrel.1 dchud.1  dqags.1  gamit.1  qagse.1  spigmr.1     zunik.1 \
clog10.1 dckder.1 dqagse.1 gamlim.1 qawc.1   spincw.1     zunk1.1 \
cmgnbn.1 dcoef.1  dqawc.1  gamln.1  qawce.1  spinit.1     zunk2.1 \
cmlri.1  dcopy.1  dqawce.1 gamma.1  qawf.1   splp.1       zuoik.1 \
cmpcsg.1 dcopym.1 dqawf.1  gamr.1   qawfe.1  splpce.1     zwrsk.1 \
cmposd.1 dcosdg.1 dqawfe.1 gamrn.1  qawo.1   splpdm.1 \
cmposn.1 dcot.1   dqawo.1  gaus8.1  qawoe.1  splpfe.1 \
cmposp.1 dcov.1   dqawoe.1 genbun.1 qaws.1   splpfl.1 \
cmptr3.1 dcpplt.1 dqaws.1  h12.1    qawse.1  splpmn.1 \
cmptrx.1 dcscal.1 dqawse.1 hfti.1   qc25c.1  splpmu.1 \
cnbco.1  dcsevl.1 dqc25c.1 hkseq.1  qc25f.1  splpup.1

#
#  Library targets and man pages.
#

#all: libslatec.a libslatec.so doc
all: libslatec.a doc

libslatec.a:
	cd static; make

libslatec.so:
	cd dynamic; make

doc: $(MAN)

${MAN}: %.1: %.f
	sed -f makedoc.sed $< > $@

install:
	@echo Installing SLATEC to $(DEST)
	mkdir -p $(DEST)/lib
	mkdir -p $(DEST)/man/man1
	cp static/libslatec.a $(DEST)/lib
	#cp dynamic/libslatec.so $(DEST)/lib
	cp *.1 $(DEST)/man/man1

clean:
	rm -f *.1; cd static; make clean; cd ../dynamic; make clean
