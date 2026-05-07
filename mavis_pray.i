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

require,"mavis_pray.h";
require,"mavis_pray_init.i";
require,"mavis_pray_lib.i";
require,"pray.i";

write,format="\n%s\n%s\n","Try running","random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=50,modes=\"dh\")";
system,"echo 'random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=50,modes=\"dh\")' | wl-copy";


func mavis_pray(coeff_offsets,ngrid,deltafoc,flux,ron,&strehlv,disp=,maxiter=,\
rseed=,verbose=,noinc=,modes=)
{
  extern pray_data;
  extern last_random_seed; // in case, to be able to repeat this random realisation

  if (!noinc) include,"mavis_pray_conf.i",1;

  //****************************************
  // parameters (dynamic) and default values
  //****************************************
  if (usemodes==[]) error,"usemodes undefined (see config file mavis_pray_conf.i)";
  if (geometry==[]) error,"geometry undefined (see config file mavis_pray_conf.i)";
  if (fovshape==[]) error,"fovshape undefined (see config file mavis_pray_conf.i)";
  if (maxiter==[])  maxiter = 100;
  if (verbose==[])  verbose = 0;
  if (deltafoc==[]) deltafoc = [-1.,0.,1.];
  if (!flux)        flux = 1000.;
  if (!ngrid)       ngrid = 4;
  if (!osampl)      osampl = 1;
  if (modes!=[])    usemodes = modes;
  // if (rseed==[])    rseed = 0.3; else
  last_random_seed = rseed;
  if (fit==[])      fit = array(1,nopt);
  if (nof(fit)!=nof(alt)) error,"sizeof(fit) != sizeof(alt)";

  if (rseed!=[]) random_seed,rseed;
  if (debug) tic;

  nopt = nof(alt);      // number of optics in train
  nfoc = nof(deltafoc); // number of extra focal distance in config
  nrot = dimsof(rotv)(0);    // number of rotation sets in config

  // Types check
  deltafoc = float(deltafoc);
  ngrid    = long(ngrid);
  flux     = float(flux);
  ron      = float(ron);

  //***********************
  // windows initialisation
  //***********************
  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"windows";
  status = init_windows();

  //**********************
  // other initialisations
  //**********************

  // build config extra-focal + rot from the entries above:
  // the idea is that we do not modify the other functions going through all
  // extra-focal distances, we just will apply for each EFD the rotation of the
  // system that applies (package initially developed w/o rotation feature)
  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"config";
  status = init_config(config);

  deltafoc_orig = deltafoc; // save original
  deltafoc = config.foc; // we replace for simplicity in the rest of the code.
  nfoc = nof(deltafoc); // update

  pupd = long(size/2/osampl);
  pupd = pupd/2*2; // we want it even
  centre = size/2+0.5;
  variance = (ron?ron^2:0.0001);   // noise variance (must be strickly positive in pray)

  // configuration printout
  if (debug) write,format="T=%.3fs -> configuration printout\n",tac();
  status = configuration_printout();

  // define and fill pray data structure
  pray_data         = pray_struct();
  pray_data.teldiam = teldiam;
  pray_data.cobs    = cobs;
  pray_data.pupd    = pupd;
  pray_data.size    = size;
  pray_data.centre  = centre;
  pray_data.nmod    = &nmod;
  pray_data.alt     = &alt;
  pray_data.active  = &active;
  pray_data.patch_diam = &(alt*0.);

  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"image centring";
  status = init_image_centring(pray_data);

  // target positions
  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"target positions";
  status = init_target_positions(geometry,fullfield,ngrid,gridpad,xpos,ypos);
  pray_data.xpos = &xpos; pray_data.ypos = &ypos;
  ntarget = nof(xpos);

  // Init ipupil
  pray_data.ipupil = &float(make_pupil(size,pupd,xc=centre,yc=centre,cobs=cobs));

  // init defs, etc used by pray, fills pray_data
  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"defs";
  status = init_defs(pray_data);

  // init masks (valid phase points at each optics)
  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"masks";
  status = init_masks(pray_data);

  // Airy pattern for Strehl estimation.
  airy = roll(abs(fft(*pray_data.ipupil,1))^2);
  peak_airy = max(airy/sum(airy));

  // create initial perturbation (coef/screens, PSFs)
  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"perturbations";
  status = init_perturbation(pray_data,coeff,cmin,cmax);

  // ok first let's do an estimate of Strehl for the deltafoc=0 and scale the
  // phase screens or coefficients from there:
  // (<- FIXME, doesn't work for coeff)
  if (strehl_normalise) {
    if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"strehl normalisation";
    write,format="Normalising Original Strehl to %.1f%% (%d iterations)\n",100*strehl_target,strehl_normalise;
    for (i=1;i<=strehl_normalise;i++) { // iterative, this works.
      status = strehl_normalisation(pray_data,coeff,config,rotv,peak_airy);
    }
  }

  if (debug) write,format="T=%.3fs -> Getting high order residuals\n",tac();
  status = get_high_order_residuals(pray_data,config);

  if (debug) write,format="T=%.3fs -> initialising %s\n",tac(),"images";
  strehlv = init_images(pray_data,config,object,start_strehl);

  if (debug&&(imask_radius_scaling!=[])) write,format="T=%.3fs -> %s\n",tac(),"Masking images to cut high frequencies";
  status = mask_images(pray_data,imask_radius_scaling);

  if (stop_after_init_images) return strehlv;

  //**************************************************************************
  // Call pray, which does the minimisation (with calls to vmlmb + pray_error)
  //**************************************************************************

  if (debug) write,format="T=%.3fs -> calling %s\n",tac(),"pray";
  res = pray(*pray_data.images,pray_data,deltafoc,variance,object,disp=disp,verbose=verbose,\
    threshold=threshold,nbiter=maxiter);

  //***************************************
  // Analyse and plot results
  //***************************************

  // check estimation results:
  origcube = *pray_data.truecube;
  // first compute mircube for 0 focus:
  compute_psfs,pray_data,0,res,amp1,amp2,nodisp=1,fromscreens=0;
  // now pray_data.mircube should contains the estimated phase at foc=0.
  // subtract estimate from original phases.
  // The PSFs are computed from "truecube"
  *pray_data.truecube = origcube - *pray_data.mircube;
  psfs = compute_psfs(pray_data,0,res*0,amp1,amp2,nodisp=1,fromscreens=1);
  disp_im = build_bigim(psfs,xpos,ypos,variance*0);
  // window,1;
  // fma; pli,disp_im; limits,square=1;
  // pltitle,swrite(format="%s","Phase corrected image");

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

  if (allof(active)) {
    write,format="%s\n","All active are set to one, no fitting to do.";
    return [start_strehl,end_strehl]; // with rms, e.g. [[0.4492,0.0343104],[0.968047,0.0202398]]
  }

  // now fit aberrations to actual DMs and compensate
  // strategy:
  // fitcube = mircube
  // project_to_dms(): we have mircube, for each non "active" optics OPT:
  // project_to_dms(): compute projection matrices
  // project_to_dms(): fit optics with modes -> coefs
  // project_to_dms(): projects coefs onto DMs, add phase to mircube DM plans
  // project_to_dms(): zero mircube(,,OPT)
  // *pray_data.truecube = origcube - *pray_data.mircube;
  // compute_psf

  return [start_strehl,end_strehl]; // with rms, e.g. [[0.4492,0.0343104],[0.968047,0.0202398]]
}
// END OF MAVIS_PRAY




// The following functions are mavis_pray wrappers to run it repeatedly
// or versus some changing parameters (number of modes etc)
func perfvsnit(nitv,ngrid,deltafoc,flux,ron,rseed=,disp=)
{
  nitv = _(0,nitv);
  stvsnit = nitv*0.+1;
  st_start_vsnit = st_start_spatial_rms_vsnit = nitv*0.+1;
  st_end_vsnit = st_end_spatial_rms_vsnit = nitv*0.+1;
  for (nn=2;nn<=nof(nitv);nn++) {
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
  if (window_exists(5)) window,5;
  else window,5,wait=1,dpi=long(dpi_target_small);
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
  if (window_exists(3)) window,3;
  else window,3,wait=1,dpi=long(dpi_target_small);
  fma;
  w = where((abs(allstend(0,)-median(allstend(0,)))<3*allstend(0,rms)));
  write,format="Kept %d out of %d samples\n",nof(w),nsamp;
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
