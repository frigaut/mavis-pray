func script1(void)
{
  require,"mavis_pray.i";

  // Study perf vs # mode
  nsamp = 10;
  nbiter = 100;
  nmov = [6,10,25,50,100,200];
  rseedv = random(nsamp);

  // nmov = [6,25,100];
  sst = est = nmov*0.;
  for (nmo=1;nmo<=nof(nmov);nmo++) {
    write,format="\n\nSimulation with %d modes\n\n",nmov(nmo);
    alls = array(0.,[3,2,2,nsamp]);
    for (n=1;n<=nsamp;n++) {
      include,"mavis_pray_conf.i",1;
      w = where(fit==1);
      nzer(w) = nmov(nmo);
      res=mavis_pray(,8,[0,-1.5,1.5],1000,0.,,disp=1,maxiter=nbiter, \
        rseed=rseedv(n),noinc=1);
      alls(,,n) = res;
    }
    alls;
    alls = median(alls,3);
    alls;
    sst(nmo) = alls(1,1); est(nmo) = alls(1,2);
    sst; est;
  }
  window,1;
  fma; limits,square=0; limits;
  plg,est,nmov,color="red"; plp,est,nmov,symbol="o",size=0.5,color="red";
  plg,sst,nmov; plp,sst,nmov,symbol="o",size=0.5;
  plmargin; range,0.,1.;
  pltitle,"End Strehl (red) vs # of fitted modes";
  xytitles,"Number of fitted modes","Strehl",[-0.015,0.];
}

func script2(void)
{
  require,"mavis_pray.i";

  // Study perf vs # mode
  nsamp = 25;
  nitv = [10,20,30,50,100,150,250,500];
  nmov = 100;
  rseedv = random(nsamp);

  // nmov = [6,25,100];
  sst = est = nitv*0.;
  for (ni=1;ni<=nof(nitv);ni++) {
    write,format="\n\nSimulation with %d iterations\n\n",nitv(ni);
    alls = array(0.,[3,2,2,nsamp]);
    for (n=1;n<=nsamp;n++) {
      include,"mavis_pray_conf.i",1;
      nit = nitv(ni);
      res=mavis_pray(,8,[0,-1.5,1.5],1000,0.,,disp=1,maxiter=nit, \
        rseed=rseedv(n),noinc=1);
      alls(,,n) = res;
    }
    alls;
    alls = median(alls,3);
    alls;
    sst(ni) = alls(1,1); est(ni) = alls(1,2);
    sst; est;
  }
  window,1;
  fma; limits,square=0; limits;
  plg,est,nitv,color="red"; plp,est,nitv,symbol="o",size=0.5,color="red";
  plg,sst,nitv; plp,sst,nitv,symbol="o",size=0.5;
  plmargin; range,0.,1.;
  pltitle,"End Strehl (red) vs # of iterations";
  xytitles,"Number of iterations","Strehl",[-0.015,0.];
}

func script3(void)
{
  require,"mavis_pray.i";

  // Study perf vs # mode
  nsamp = 25;
  nit = 50;
  nmo = 100;
  slov = [-1,-2,-2.5,-2.66,-3];
  rseedv = random(nsamp);

  sst = est = slov*0.;
  for (ns=1;ns<=nof(slov);ns++) {
    write,format="\n\nSimulation with slope = %.2f\n\n",slov(ns);
    alls = array(0.,[3,2,2,nsamp]);
    for (n=1;n<=nsamp;n++) {
      include,"mavis_pray_conf.i",1;
      ps_slope = slov(ns);
      res=mavis_pray(,8,[0,-1.5,1.5],1000,0.,,disp=1,maxiter=nit, \
        rseed=rseedv(n),noinc=1);
      alls(,,n) = res;
    }
    alls;
    alls = median(alls,3);
    alls;
    sst(ns) = alls(1,1); est(ns) = alls(1,2);
    sst; est;
  }
  window,1;
  fma; limits,square=0; limits;
  plg,est(::-1),slov(::-1),color="red"; plp,est(::-1),slov(::-1),symbol="o",size=0.5,color="red";
  plg,sst(::-1),slov(::-1); plp,sst(::-1),slov(::-1),symbol="o",size=0.5;
  plmargin; range,0.,1.;
  pltitle,"End Strehl (red) vs power spectrum slope";
  xytitles,"Power Spectrum slope","Strehl",[-0.015,0.];
}

status = script3();
