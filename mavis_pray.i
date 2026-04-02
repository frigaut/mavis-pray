/* mavis_pray.i
example of call:
Overall wrappers:
res=perfvsnit([100],4,[-1,0,1]*1.5,10000,0.1,rseed=random(),disp=0);
res=doitall([20,30,50,100],10,4,[-1,0,1]*1.8,1,0.01)

Configuration is defined in mavis_pray_config.i. Edit it and run
perfvsnit() as above. The DM phase window may need to be killed if
the number of DM is changed.


Example of a run:
> res=perfvsnit([50],8,[0,-1.5,1.5],100,0.001,rseed=random(),disp=1)
--------------------------------------------------------------------
Pray for MAVIS (mavis_pray and pray) - System Configuration
Number of optics: 3, Number of extra-focal pos.: 6, number of rotation: 2
Extra focal distances: [0,-1.5,1.5]
Optics conjug. altitude: [0,6000,13500]
Number of modes (DH) per optics: [50,50,50]
Fit optics?: [1,1,1,1,1,1,1,1,1,1,1,1,1]
Optics rotation, config 1: [0,0,0]
Optics rotation, config 2: [180,0,0]
Number of sources = 8 across 30" FoV, 59 total, hexagonal geometry
source flux = 100.0, RON= 0.001
--------------------------------------------------------------------
Strehl over FoV (rot=[0,0,0]): avg=0.342289 rms=0.034295
Strehl over FoV (rot=[180,0,0]): avg=0.375084 rms=0.038969
Using coeff = 0 as first guess
# Iter.   Time (ms)    Eval. Reject.       Obj. Func.           Grad.       Step
# ---------------------------------------------------------------------------------
      0     379.680       1       0   3.264587827870292e+10   5.958e+09   0.000e+00
     10    5817.125      14       1   1.053084155868230e+09   1.069e+09   1.000e+00
     20    9707.075      24       1   1.057027390264975e+08   2.385e+08   1.000e+00
     30   13553.921      34       1   4.101288959418561e+07   8.583e+07   1.000e+00
     40   17399.916      44       1   1.767627173208634e+07   5.198e+07   1.000e+00
     50   21251.018      54       1   9.551682008189257e+06   9.045e+07   1.000e+00
# Termination: too many iterations
Strehl over FoV (rot=[0,0,0]): avg=0.998477 rms=0.000507
Strehl over FoV (rot=[180,0,0]): avg=0.998474 rms=0.000481

Additional upgrades and questions:
- do I need to scale the rms according to the altitude (the diameter of the beam
on the particular optics). Probably. TODO.

*/

write,format="%s\n%s\n","Try running","res=mavis_pray(,8,[0,-1.5,1.5],1000,0.,,disp=1,maxiter=50,rseed=random())";
write,format="%s\n%s\n","Or","perfvsnit([50],8,[0,-1.5,1.5],100,0.001,rseed=random(),disp=1)";
write,format="%s\n%s\n","Or","res=doitall([10,20,50],10,8,[0.,-1.5,1.5],10000,0.,disp=1)";

system,"echo 'res=mavis_pray(,8,[0,-1.5,1.5],1000,0.,,disp=1,maxiter=50,rseed=random())' | wl-copy";

require,"pray.i";
require,"mavis_pray_lib.i";

func mavis_pray(coeff_offsets,ngrid,deltafoc,flux,ron,&strehlv,disp=,maxiter=,\
	rseed=,verbose=,noinc=)
{
  extern pray_data;

  if (!noinc) include,"mavis_pray_conf.i",1;

  // parameters (dynamic)
  if (maxiter==[])  maxiter = 100;
  if (verbose==[])  verbose = 0;
  if (deltafoc==[]) deltafoc = [-1.,0.,1.];
  if (!flux)        flux = 1000.;
  if (!ngrid)       ngrid = 4;
  if (usemodes==[]) error,"usemodes undefined (see config file mavis_pray_conf.i)";
  if (geometry==[]) error,"geometry undefined (see config file mavis_pray_conf.i)";
  if (fovshape==[]) error,"fovshape undefined (see config file mavis_pray_conf.i)";

  // Types check
  deltafoc = float(deltafoc);
  ngrid    = long(ngrid);
  flux     = float(flux);
  ron      = float(ron);

  nopt = numberof(alt); // number of optics in train
  nfoc = numberof(deltafoc); // number of extra focal distance in config
  nrot = dimsof(rotv)(0); // number of rotation sets in config

  if (fit==[]) fit = array(1,nopt);
  if (numberof(fit)!=nopt) error,"numberof(fit) != nopt";

  // build config extra-focal + rot from the entries above:
  config = array(config_struct,nfoc*nrot);
  config.foc = array(deltafoc,nrot)(*);
  tmp = [];
  for (i=1;i<=nrot;i++) grow,tmp,array(i,nfoc);
  config.roti = tmp;
  deltafoc_orig = deltafoc;
  deltafoc = config.foc; // we replace for simplicity in the rest of the code.
  nfoc = numberof(deltafoc);

  pupd = long(size/2/osampl);

  // window init
  win = [1,3,7,8,15];
  for (i=1;i<=numberof(win);i++) {
    if (window_exists(win(i))) {
      window,win(i);
      if (win(i)!=15) { limits; fma; }
    }
  }

  nw = clip(nopt,4,8); nw = 8;
  if (!window_exists(3)) plsplit,2,nw,win=3,style="nobox.gs",square=1,dpi=long(dpi_target*1.5),margin=-0.02,vp=[0.206+0.0115*(nw-3),0.656+0.0115*(nw-3),0.44,0.85];

  // define and fill pray data structure
  pray_data         = pray_struct();
  pray_data.teldiam = teldiam;
  pray_data.cobs    = cobs;
  pray_data.pupd    = pupd;
  pray_data.size    = size;
  pray_data.nzer    = &nzer;
  pray_data.alt     = &alt;

  // compute values of variables used later
  // noise variance
  variance = (ron?ron^2:0.0001);

  // target positions
  if (strpart(geometry,1:3)=="squ") {
    tmp = (indices(ngrid)-ngrid/2-1);
    tmp += (odd(ngrid)?0:0.5);
    tmp *= fullfield/ngrid;
    xpos = tmp(,,1)(*); ypos = tmp(,,2)(*);
  } else if (strpart(geometry,1:3)=="hex") {
    tmp = (indices(ngrid+4)-(ngrid+4)/2-1.);
    tmp(,,1) += 0.5*(indgen(ngrid+4)%2)(-,);
    tmp(,,1) += (odd(ngrid)?0:0.5);
    tmp(,,2) *= sin(60*pi/180.);
    tmp *= fullfield/(ngrid-0.8);
    w = where(abs(tmp(,,1))<(fullfield/2.)&abs(tmp(,,2))<(fullfield/2.));
    xpos = tmp(,,1)(w); ypos = tmp(,,2)(w);
  } else error,"geometry undefined";

  if (fovshape=="round") {
    if (gridpad==[]) gridpad = 0.;
    w = where(abs(xpos,ypos)<=(fullfield/2.+gridpad));
    xpos = xpos(w); ypos = ypos(w);
  }

  pray_data.xpos = &xpos;
  pray_data.ypos = &ypos;

  // configuration printout
  write,format="%s\n","--------------------------------------------------------------------";
  write,format="%s\n","Pray for MAVIS (mavis_pray and pray) - System Configuration";
  write,format="Number of optics: %d, Number of extra-focal pos.: %d, number of rotation: %d\n",nopt,numberof(deltafoc),dimsof(rotv)(0);
  write,format="%s","Extra focal distances: "; deltafoc_orig;
  write,format="%s","Optics conjug. altitude: "; alt;
  write,format="Number of modes (%s) per optics: ",strcase(1,usemodes); nzer;
  write,format="%s","Fit optics?: "; fit;
  perfrom = (initphase=="screens"?"Power spectrum":"Mode coefficients");
  if (initphase=="screens") perfrom += swrite(format=" (slope=%.3f)",ps_slope);
  write,format="Init phase perturbation from %s\n",perfrom;
  for (i=1;i<=nrot;i++) { write,format="Optics rotation, config %d: ",i; rotv(,i); }
  write,format="Number of sources = %d across %.0f\" FoV, %d total, %s geometry\n",ngrid,fullfield,numberof(xpos),geometry;
  write,format="source flux = %.1f, RON= %g\n",flux,ron;
  write,format="%s\n","Current random seed stored in extern last_random_seed, use as rseed to repeat current";
  write,format="%s\n","--------------------------------------------------------------------";

  // window,8; fma; plp,ypos,xpos; limits; limits,square=1; hitReturn;

  // pray init
  status = pray_init(pray_data);

  // init for Strehl estmation.
  airy = roll(abs(fft(*pray_data.ipupil,1))^2);
  airy = airy/sum(airy);
  peak_airy = max(airy);

  // create psfs
  coeff = cmax = cmin = [];
  nz12 = nzer(cum);
  extern last_random_seed;
  if (rseed!=[]) last_random_seed = rseed;
  random_seed,(rseed?rseed:0.3);
  // var = 0.;
  for (i=1;i<=nopt;i++) {
    c = weight(i)*(random(nzer(i))-0.5)/sqrt(indgen(nzer(i)));
    c(1) = 0.;
    cmx = 10*weight(i)/sqrt(indgen(nzer(i)))*fit(i);
    cmx(1) = 0.; // FIXME
    // if (i>=2) cmx *= 0;
    cmn = -cmx;
    grow,coeff,c;
    grow,cmax,cmx;
    grow,cmin,cmn;
  }
  // write,var,exp(-var/2.);
  coeff = coeff(*); cmax = cmax(*); cmin = cmin(*);
  extern truecoeff; truecoeff = coeff;
  if (coeff_offsets!=[]) coeff -= coeff_offsets;

  // fill mircube with phases
  if (initphase=="screens") {
	  for (i=1;i<=nopt;i++) {
	  	(*pray_data.mircube)(,,i) = make_phase_screens(*pray_data.ipupil,lambda, \
					nm_rmsv(i),ps_slope,rseed=rseed+i*0.01,remove_tt=1,remove_foc=1);
	  }
		// zero coefs:
		coeff *=0;
		truecoeff *= 0;
		cmin = cmin * 1e3;
		cmax = cmax * 1e3;
  }


  // compute true phase from truecoeffs for display comparison
  extern truecube;
  truecube = array(0.,[3,size,size,nopt]);
  def = *pray_data.def;
  cpt = 0;
  for (i=1;i<=nopt;i++) {
    zv = cpt+indgen(nzer(i));
    lpup = def(,,zv(2)) != 0;
    if (initphase=="screens") truecube(,,i) = (*pray_data.mircube)(,,i)*lpup; \
    else truecube(,,i) = def(,,zv)(,,+)*truecoeff(zv)(+);
    cpt += nzer(i);
  }
  def = [];

  res = [];
  for (i=1;i<=nfoc;i++) {
    // coeff(1) = deltafoc(i);
    // write,"Computing initial images";
    // write,config.roti(i),rotv(,config.roti(i));
    grow,res,&psfs_from_coeffs(pray_data,deltafoc(i),coeff,amp1,amp2, \
    rotv=rotv(,config.roti(i)),nodisp=1,fromscreens=(initphase=="screens"));
  }
  // OK so it's the implementation that is wrong.
  // probably we need to save the truecube and compare to it later.

  // simple square object ...
  object=array(float,[2,size,size]);
  object(size/2+1,size/2+1) = flux;
  // creating images
  ntarget = numberof(xpos);
  images = array(float,[4,size,size,ntarget,nfoc]);
  strehlv = array(0.,ntarget);
  for (n=1;n<=nfoc;n++) {
    for (i=1;i<=ntarget;i++) {
      // images are centred, object is centred res is rol (so eclat(res) is centred)
      // images(,,i,n) = fft_convolve(object,eclat((*res(n))(,,i)));
      images(,,i,n) = roll((*res(n))(,,i))*flux;
      if (deltafoc(n)==0) \
        strehlv(i) = max(images(,,i,n)/sum(images(,,i,n)))/peak_airy;
    }
    if (deltafoc(n)==0) {
      rotvstr = strjoin(swrite(format="%.0f",rotv(,config(n).roti)),",");
      write,format="\033[31mStrehl over FoV (rot=[%s]): avg=%.1f%%\033[0m rms=%.1f%%\n", \
        rotvstr,100*avg(strehlv),100*strehlv(rms);
      if (n==1) start_strehl = [avg(strehlv),strehlv(rms)];
      // FIXME, right now the start_strehl is for 180 while end_strehl is for 0
    }

    if (disp&&(n==1)) {
      // display for first defoc plane
      disp_im = build_bigim(roll(images(,,,n)),xpos,ypos,variance);
      extern bigim2; bigim2 = disp_im;
      if (window_exists(n)) window,n;
      else window,n,wait=1,dpi=dpi_target;
      fma; pli,disp_im; limits,square=1;
      pltitle,swrite(format="%.2f-focus images - data",deltafoc(n));
      pause,50;
    }
  }
  images = images + random_normal(dimsof(images))*sqrt(variance);

  pray_data.images = &images;

  res = pray(images,pray_data,deltafoc,variance,object,disp=disp,verbose=verbose,\
    threshold=threshold,nbiter=maxiter);

  // check estimation results:
  origcube = truecube;
  // first compute mircube for 0 focus:
  psfs_from_coeffs,pray_data,0,res,amp1,amp2,nodisp=1,fromscreens=0;
  // now pray_data.mircube should contains the estimated phase at foc=0.
  // subtract estimate from original phases:
  truecube = origcube - *pray_data.mircube;
  psfs = psfs_from_coeffs(pray_data,0,res*0,amp1,amp2,nodisp=1,fromscreens=1);
  disp_im = build_bigim(psfs,xpos,ypos,variance*0);
  window,1;
  fma; pli,disp_im; limits,square=1;
  pltitle,swrite(format="%s","Phase corrected image");

  strehlv = array(0.,ntarget);
  for (i=1;i<=ntarget;i++) {
    strehlv(i) = max(psfs(,,i)/sum(psfs(,,i)))/peak_airy;
  }
  write,format="\033[31mStrehl over FoV after compensation: avg=%.1f%%\033[0m rms=%.1f%%\n", \
    100*avg(strehlv),100*strehlv(rms);
  end_strehl = [avg(strehlv),strehlv(rms)];

  // compute strehl if all the fit=1 optics were perfectly corrected and the
  // fit=0 not corrected at all.
  tmp = nm_rmsv*(1-fit);
  st = exp(-(2*pi/lambda*sqrt(sum(tmp^2.))));
  write,format="Strehl due to fit=0 optics only (full original rms): %.1f%%\n",st*100;

  return [start_strehl,end_strehl]; // with rms, e.g. [[0.4492,0.0343104],[0.968047,0.0202398]]
}
// END OF MAVIS_PRAY

func perfvsnit(nitv,ngrid,deltafoc,flux,ron,rseed=,disp=)
{
  nitv = _(0,nitv);
  stvsnit = nitv*0.+1;
  st_start_vsnit = st_start_spatial_rms_vsnit = nitv*0.+1;
  st_end_vsnit = st_end_spatial_rms_vsnit = nitv*0.+1;
  for (nn=2;nn<=numberof(nitv);nn++) {
    strehls = mavis_pray(,ngrid,deltafoc,flux,ron,init_strehlv,disp=disp,maxiter=nitv(nn),rseed=rseed);
    st_start_vsnit(nn) = strehls(1,1);
    st_start_spatial_rms_vsnit(nn) = strehls(2,1);
    st_end_vsnit(nn) = strehls(1,2);
    st_end_spatial_rms_vsnit(nn) = strehls(2,2);
  }
  // fill in the uncorrected case
  st_start_vsnit(1) = strehls(1,1);
  st_start_spatial_rms_vsnit(1) = strehls(2,1);
  st_end_vsnit(1) = strehls(1,1);
  st_end_spatial_rms_vsnit(1) = strehls(2,1);

  cw = current_window();
  if (window_exists(15)) window,15;
  else window,15,wait=1,dpi=long(dpi_target_small);
  errnm = lambda/2/pi*sqrt(-log(st_end_vsnit));
  fma;
  plg,errnm,nitv,width=3;
  plp,errnm,nitv,symbol="o",size=0.5,width=3;
  plt,swrite(format="%.1f",errnm(0)),nitv(0)+2,errnm(0)+2,justify="LB",tosys=1,height=10;
  plmargin; range,0
  xytitles,"Number of iterations","Phase error [nm]",[-0.015,0.];
  pltitle,swrite(format="Pray: %dx%d grid, efd=%s, flux=%.0f, RON=%g",\
      ngrid,ngrid,print(deltafoc)(1),flux*1.,ron*1.);
  window,cw,wait=1;
  pause,50;
  return [nitv,errnm,st_start_vsnit,st_start_spatial_rms_vsnit,st_end_vsnit,st_end_spatial_rms_vsnit];
}

func doitall(nitv,nsamp,ngrid,deltafoc,flux,ron,disp=)
{
  allststart = allstend = array(_(0.,nitv*0.),nsamp);
  for (k=1;k<=nsamp;k++) {
    write,format="\ndoitall() iteration: %d/%d\n",k,nsamp;
    res = perfvsnit(nitv,ngrid,deltafoc,flux,ron,rseed=random(),disp=disp);
    allststart(,k) = res(,3); // start strehls vs nit
    allstend(,k) = res(,5); // end strehls vs nit
    nitv2 = res(,1); // nit vector, including 0
  }
  if (window_exists(1)) window,1;
  else window,1,wait=1,dpi=long(dpi_target_small);
  fma;
  w = where((abs(allstend(0,)-median(allstend(0,)))<3*allstend(0,rms)));
  write,format="Kept %d out of %d samples\n",numberof(w),nsamp;
  ststartavg = allststart(,w)(,avg); stendavg = allstend(,w)(,avg);
  ststartrms = allststart(,w)(,rms); stendrms = allstend(,w)(,rms);
  // write,format="Final Strehl [nm] = %.2fnm +/- %.2f\n",erravg(0),errrms(0);
  fma;
  plg,ststartavg,nitv2;
  plp,ststartavg,nitv2,symbol="o",size=0.5;
  pleb,ststartavg,nitv2,dy=ststartrms;
  plg,stendavg,nitv2,color="red";
  plp,stendavg,nitv2,symbol="o",size=0.5,color="red";
  pleb,stendavg,nitv2,dy=stendrms,color="red";
  plmargin; range,0,1;
  xytitles,"Number of iterations",swrite(format="Strehl @ %.0fnm",float(lambda)),[-0.015,0.];
  pltitle,"End Strehl (red) vs # of iterations";
  // pltitle,swrite(format="Pray: %dx%d grid, efd=%s, flux=%.0f, RON=%g",\
  //   ngrid,ngrid,print(deltafoc)(1),flux*1.,ron*1.);
}
