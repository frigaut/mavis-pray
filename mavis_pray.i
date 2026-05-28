/* mavis_pray.i
example of call:
Overall wrappers:
res=perfvsnit([100],4,[-1,0,1]*1.5,10000,0.1,rseed=random(),disp=0);
res=do_stats_vs_nit([20,30,50,100],10,4,[-1,0,1]*1.8,1,0.01)

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

write,format="\n%s\n%s\n","Try running","random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=50,modes=\"zer\")";
write,format="or \"%s\"\n","case=10; random_seed,0.75; do_stats,[50],20,8,[0.,-1.5,1.5],1e5,0,disp=1";
system,"echo 'random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=50,modes=\"dh\")' | wl-copy";

func mavis_pray(coeff_offsets,ngrid,deltafoc,flux,ron,&strehlv,disp=,maxiter=,\
rseed=,verbose=,noinc=,modes=,skip_proj=)
/* DOCUMENT func mavis_pray(coeff_offsets,ngrid,deltafoc,flux,ron,&strehlv,disp=,maxiter=,\
rseed=,verbose=,noinc=,modes=,skip_proj=)
 * mavis_pray simulates the MAVIS NCPA process. It does:
 * initliasation of the optics surface for all optics
 * computes the images
 * inits a number of system parameters
 * call pray() which does the minimisation
 * project on the active optics
 * output quality metrics, strehl avg , strehl rms
 * Configuration is read from executing mavis_pray_conf.i
 * PARAMETERS:
 * coeff_offsets: possible offsets for coefficients (rarely used, historical)
 * ngrid: side number of images in the output focal plane (e.g. 4-8)
 * deltafoc: extra focal distances, e.g. [0.,-1.5,1.5]. They don't have to be
 *   in increasing order. First one is the one for which images are displayed,
 *   so usually 0 to get in focus images.
 * flux: flua per individual image, in ADU (e.g. 100000)
 * ron: read out noise, in ADU
 * strehlv: output strehl (unused)
 * disp= set to display stuff
 * maxiter= maximum number of iterations in the minimisatioin process
 * rseed= random seed. can be useful to repeat similar runs. If not specified,
 *   will take a random seed
 * verbose= write out more stuff
 * noinc= if set, mavis_pray_conf.i will not be included. Can be useful to
 *   call again mavis_pray() after an initial call to change some parameters
 * modes= "zer" of "dh". Note that KL is broken at this time
 * skip_proj= skip the projection step, can be executed later using
 *   simple_projection_only(pray_data)
*/
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
  if (verbose==[])  verbose = 1;
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
  tic;

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
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"windows";
  status = init_windows();

  //**********************
  // other initialisations
  //**********************

  // build config extra-focal + rot from the entries above:
  // the idea is that we do not modify the other functions going through all
  // extra-focal distances, we just will apply for each EFD the rotation of the
  // system that applies (package initially developed w/o rotation feature)
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"config";
  status = init_config(config);

  deltafoc_orig = deltafoc; // save original
  deltafoc = config.foc; // we replace for simplicity in the rest of the code.
  nfoc = nof(deltafoc); // update

  pupd = long(size/2/osampl);
  pupd = pupd/2*2; // we want it even
  centre = size/2+0.5;
  variance = (ron?ron^2:0.0001);   // noise variance (must be strickly positive in pray)

  // configuration printout
  if (verbose) write,format="T=%.3fs -> \033[32mconfiguration printout\033[0m\n",tac();
  status = configuration_printout();

  // define and fill pray data structure
  pray_data         = pray_struct();
  pray_data.teldiam = teldiam;
  pray_data.cobs    = cobs;
  pray_data.pupd    = pupd;
  pray_data.size    = size;
  pray_data.ngrid    = ngrid;
  pray_data.centre  = centre;
  pray_data.nmod    = &nmod;
  pray_data.alt     = &alt;
  pray_data.active  = &active;
  pray_data.patch_diam = &(alt*0.);

  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"image centring";
  status = init_image_centring(pray_data);

  // target positions
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"target positions";
  status = init_target_positions(geometry,fullfield,ngrid,gridpad,xpos,ypos);
  pray_data.xpos = &xpos; pray_data.ypos = &ypos;
  ntarget = nof(xpos);

  // Init ipupil
  pray_data.ipupil = &float(make_pupil(size,pupd,xc=centre,yc=centre,cobs=cobs));

  // init defs, etc used by pray, fills pray_data
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"defs";
  status = init_defs(pray_data);

  // init masks (valid phase points at each optics)
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"masks";
  status = init_masks(pray_data);

  // Airy pattern for Strehl estimation.
  airy = roll(abs(fft(*pray_data.ipupil,1))^2);
  pray_data.peak_airy = max(airy/sum(airy));
  airy = [];

  // create initial perturbation (coef/screens, PSFs)
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"perturbations";
  status = init_perturbation(pray_data,coeff,cmin,cmax);

  // debug adding focus to check removal works
  // *pray_data.truecube = *pray_data.truecube+(*pray_data.focus)(,,-)*10;
  // remove average focus so that extra focal =0 images are in focus in average:
  if (verbose) write,format="T=%.3fs -> \033[32mRemoving average focus\033[0m\n",tac();
  foc = *pray_data.focus;
  for (i=1;i<=nopt;i++) {
    phase = (*pray_data.truecube)(,,i);
    wpup = where((*pray_data.pupil)(,,i));
    avg_focus = sum(phase(wpup)*foc(wpup))/sum(foc(wpup)*foc(wpup));
    (*pray_data.truecube)(,,i) -= avg_focus*foc;
  }

  // ok first let's do an estimate of Strehl for the deltafoc=0 and scale the
  // phase screens or coefficients from there:
  // (<- FIXME, doesn't work for coeff)
  if (strehl_normalise) {
    if (verbose) write,format="T=%.3fs -> \033[32mStrehl normalisation to %.1f%% (%d iterations)\033[0m\n",tac(),100*strehl_target,strehl_normalise;
    for (i=1;i<=strehl_normalise;i++) { // iterative, this works.
      status = strehl_normalisation(pray_data,coeff,config,rotv);
    }
  }

  if (skip_high_order!=1) {
    if (verbose) write,format="T=%.3fs -> \033[32mGetting high order residuals\033[0m\n",tac();
    status = get_high_order_residuals(pray_data,config);
  }

  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"images";
  strehlv_init = init_images(pray_data,config,object,start_strehl);

  if (debug&&(imask_radius_scaling!=[])) write,format="T=%.3fs -> %s\n",tac(),"Masking images to cut high frequencies";
  status = mask_images(pray_data,imask_radius_scaling);

  if (stop_after_init_images) return strehlv;

  //**************************************************************************
  // Call pray, which does the minimisation (with calls to vmlmb + pray_error)
  //**************************************************************************

  if (verbose) write,format="T=%.3fs -> \033[32mcalling %s\033[0m\n",tac(),"pray";
  res = pray(*pray_data.images,pray_data,deltafoc,variance,object,disp=disp,verbose=verbose,\
    threshold=threshold,nbiter=maxiter);
  pray_data.coeffs = &res;

  //***************************************
  // Analyse and plot results
  //***************************************

  // check estimation results:
  pray_data.origcube = &(*pray_data.truecube*1);
  // first compute mircube for 0 focus:
  compute_psfs,pray_data,0,res,amp1,amp2,nodisp=1,fromscreens=0;
  // now pray_data.mircube should contains the estimated phase at foc=0.
  // subtract estimate from original phases.
  // The PSFs are computed from "truecube" (when fromscreens==1)
  *pray_data.truecube = *pray_data.origcube - *pray_data.mircube;
  psfs = compute_psfs(pray_data,0,res*0,amp1,amp2,nodisp=1,fromscreens=1);
  disp_im = build_bigim(psfs,xpos,ypos,variance*0);

  strehlv = array(0.,ntarget);
  for (i=1;i<=ntarget;i++) { //FIXME this is only for rot = 0
    strehlv(i) = max(psfs(,,i)/sum(psfs(,,i)))/pray_data.peak_airy;
  }
  strehlv_corr = strehlv;
  write,format="\033[31mStrehl over FoV after compensation: avg=%.1f%%\033[0m rms=%.1f%%\n", \
    100*avg(strehlv),100*strehlv(rms);
  end_strehl = [avg(strehlv),strehlv(rms)];

  // compute strehl if all the fit=1 optics were perfectly corrected and the
  // fit=0 not corrected at all.
  tmp = nm_rmsv*(1-fit);
  st = exp(-(2*pi/lambda*sqrt(sum(tmp^2.))));
  write,format="Strehl due to fit=0 optics only (full original rms): %.1f%%\n",st*100;

  if (allof(active)) {
    write,format="%s\n","All optic active flags are set to one, no fitting to do.";
    strehlv_end = strehlv_corr;
  } else {
    strehlv_end = simple_projection_only(pray_data);
  }

  return [strehlv_init,strehlv_corr,strehlv_end]; // with rms, e.g. [[0.4492,0.0343104],[0.968047,0.0202398]]
}
// END OF MAVIS_PRAY


func simple_projection_only(pd)
{
  require,"projection.i";
  write,format="%s\n","\033[32mProjecting passive optics shape onto DMs\033[0m";
  compute_psfs,pd,0,*pd.coeffs,amp1,amp2,nodisp=1,fromscreens=0;
  // now pd.mircube should contains the estimated phase at foc=0.
  // subtract estimate from original phases.
  // The PSFs are computed from "truecube" (when fromscreens==1)
  pd = simple_project(pd);
  *pd.truecube = *pd.origcube - *pd.mircube;
  psfs = compute_psfs(pd,0,,amp1,amp2,nodisp=1,fromscreens=1);
  disp_im = build_bigim(psfs,*pd.xpos,*pd.ypos,0);
  window,1; plsys,1; pli,disp_im; pltitle_vp,"Final"; redraw;
  ntarget = nof(*pd.xpos);
  strehlv = array(0.,ntarget);
  airy = roll(abs(fft(*pd.ipupil,1))^2);
  for (i=1;i<=ntarget;i++) {
    strehlv(i) = max(psfs(,,i)/sum(psfs(,,i)))/pd.peak_airy;
  }
  end_strehlv = strehlv;
  write,format="\033[31mStrehl over FoV after DM projection: avg=%.1f%%\033[0m rms=%.1f%%\n", \
  100*avg(strehlv),100*strehlv(rms);
  end_strehl = [avg(strehlv),strehlv(rms)];
  return end_strehlv;
}


struct allstrehl_st {
  long nit;
  long nsamp;
  pointer strehls;
}

func do_stats(nitv,nsamp,ngrid,deltafoc,flux,ron,rseed=,disp=)
{
  // name = "do_stats_"+totxt(long(random()*100000));
  d=timestamp(); name = "do_stats_"+streplace(d,strfind(" ",d,n=10),"-");
  system,"mkdir "+name;
  system,"cp -p mavis_pray_conf.i "+name+"/";

  strehl_start = array(allstrehl_st(),nof(nitv)*nsamp);
  strehl_corr  = array(allstrehl_st(),nof(nitv)*nsamp);
  strehl_end   = array(allstrehl_st(),nof(nitv)*nsamp);
  info,strehl_start;
  ind = 1;
  rejected = 0;
  for (ns=1;ns<=nsamp;ns++) {
    write,format="\n\033[32mSample %d/%d\033[0m\n\n",ns,nsamp;
    for (nn=1;nn<=nof(nitv);nn++) {
      strehls = mavis_pray(,ngrid,deltafoc,flux,ron,init_strehlv,disp=disp,maxiter=nitv(nn),rseed=rseed);
      if (strehls(avg,2)<0.3) {
        write,format="%s\n","->>> Rejected run!";
        rejected = rejected+1;
        continue;
      }
      strehl_start(ind).nit = nitv(nn);
      strehl_start(ind).nsamp = ns;
      strehl_start(ind).strehls = &strehls(,1);
      strehl_corr(ind).nit = nitv(nn);
      strehl_corr(ind).nsamp = ns;
      strehl_corr(ind).strehls = &strehls(,2);
      strehl_end(ind).nit = nitv(nn);
      strehl_end(ind).nsamp = ns;
      strehl_end(ind).strehls = &strehls(,3);
      ind = ind+1;
    }
    plot_do_stats,strehl_start,strehl_corr,strehl_end,rejected;
  }
  strehl_start = strehl_start(1:ind-1);
  strehl_corr  = strehl_corr(1:ind-1);
  strehl_end   = strehl_end(1:ind-1);
  write,format="Rejected runs: %d\n",rejected;
  write,format="Saving data in folder %s/\n",name;
  f = createb(name+"/do_stats.dat");
  save,f,strehl_start,strehl_corr,strehl_end,lambda,case,ngrid,deltafoc,flux,ron,rejected;
  close,f;
}

func plot_do_stats(strehl_start,strehl_corr,strehl_end,rejected,binsize=,name=)
{
  if (binsize==[]) binsize=2.5;
  if (name) {
    f = openb(name+"/do_stats.dat");
    restore,f,strehl_start,strehl_corr,strehl_end,lambda,case,ngrid,deltafoc,flux,ron,rejected;
    close,f;
  }

  // Strehl histograms only for nit = max(nitv)
  if (max(strehl_start.nit)==0) return;
  w = wheremax(strehl_start.nit);
  nvalid = sum(strehl_start.nit!=0);
  nit = max(strehl_start.nit);
  ss = sc = se = [];
  cw = current_window();
  if (cw!=10) window,10,wait=1,style="clean.gs";
  fma; limits,square=0; limits;
  for (i=1;i<=nof(w);i++) {
    grow,ss,*(strehl_start(w(i)).strehls)*100;
    grow,sc,*(strehl_corr(w(i)).strehls)*100;
    grow,se,*(strehl_end(w(i)).strehls)*100;
  }
  hy = histo2(ss,hx,binsize=binsize);
  plh,hy,hx,color=torgb(tokyonight(3)),width=3;
  hy = histo2(sc,hx,binsize=binsize);
  plh,hy,hx,color=torgb(tokyonight(2)),width=3;
  hy = histo2(se,hx,binsize=binsize);
  plh,hy,hx,color=torgb(tokyonight(1)),width=3;
  pltitle,swrite(format="NCPA performance for %d iterations",nit);
  xytitles,swrite(format="Strehl@%dnm",long(lambda)),"Number in bin",[-0.012,0.008];
  y0 = y1 = limits()(4)*0.97; dy=limits()(4)*0.06; y1+=limits()(4)*0.012;
  plg,[y1,y1],[0.,2],width=20,color=torgb(tokyonight(3)); y1-=dy;
  plg,[y1,y1],[0.,2],width=20,color=torgb(tokyonight(2)); y1-=dy;
  plg,[y1,y1],[0.,2],width=20,color=torgb(tokyonight(1)); y1-=dy;
  plt,swrite(format="Strehl start median = %.1f%%",median(ss)),6,y0,tosys=1,height=12,color=torgb(tokyonight(3)),justify="LA"; y0-=dy;
  plt,swrite(format="Strehl fitted median = %.1f%%",median(sc)),6,y0,tosys=1,height=12,color=torgb(tokyonight(2)),justify="LA"; y0-=dy;
  plt,swrite(format="Strehl projected median = %.1f%%",median(se)),6,y0,tosys=1,height=12,color=torgb(tokyonight(1)),justify="LA"; y0-=dy;
  plmargin; range,0;

  write,format="%d rejected runs out of %d for this batch\n",rejected,nvalid;

  if (name) {
    pngcrop,name+"/do_stats";
    write,format="Plot saved in %s/do_stats.png\n",name;
  }
  hitReturn;

  /*
  // Strehl and strehl rms vs nit
  // find all unique nit:
  nitv = strehl_start.nit;
  nitv = nitv(sort(nitv));
  anit = [];
  grow,anit,nitv(1);
  w = where(nitv!=anit(0));
  while (nof(w)) {
    nitv = nitv(w);
    grow,anit,nitv(1);
    w = where(nitv!=anit(0));
  }
  nitv = anit;
  nnitv = nof(nitv);
  ss_avg = ss_rms = array(0.,nnitv);
  sc_avg = sc_rms = array(0.,nnitv);
  se_avg = se_rms = array(0.,nnitv);
  fma;
  for (ni=1;ni<=nnitv;ni++) {
    w = where(strehl_start.nit==nitv(ni));
    ss = sc = se = [];
    for (i=1;i<=nof(w);i++) {
      grow,ss,*(strehl_start(w(i)).strehls)*100;
      grow,sc,*(strehl_corr(w(i)).strehls)*100;
      grow,se,*(strehl_end(w(i)).strehls)*100;
    }
    // grow,ss,*(strehl_start(w(i)).strehls)*100;
    ss_avg(ni) = ss(avg); ss_rms(ni) = ss(rms);
    sc_avg(ni) = sc(avg); sc_rms(ni) = sc(rms);
    se_avg(ni) = se(avg); se_rms(ni) = se(rms);
  }
  ss_avg = _(ss_avg(1),ss_avg); ss_rms = _(ss_rms(1),ss_rms);
  sc_avg = _(ss_avg(1),sc_avg); sc_rms = _(ss_rms(1),sc_rms);
  se_avg = _(ss_avg(1),se_avg); se_rms = _(ss_rms(1),se_rms);
  nitv = _(0,nitv);
  plg,sc_avg,nitv,width=3,color=torgb(tokyonight(2));
  plp,sc_avg,nitv,symbol="o",size=0.7,width=3,color=torgb(tokyonight(2));
  pleb,sc_avg,nitv,dy=sc_rms,color=torgb(tokyonight(2));
  plg,se_avg,nitv,width=3,color=torgb(tokyonight(1));
  plp,se_avg,nitv,symbol="o",size=0.7,width=3,color=torgb(tokyonight(1));
  pleb,se_avg,nitv,dy=se_rms,color=torgb(tokyonight(1));
  plg,ss_avg,nitv,width=3,color=torgb(tokyonight(3));
  plp,ss_avg,nitv,symbol="o",size=0.7,width=3,color=torgb(tokyonight(3));
  pleb,ss_avg,nitv,dy=ss_rms,color=torgb(tokyonight(3));
  // pltitle,"End Strehl (all optics: red, DMs: green) vs # of iterations";
  pltitle,"Performance vs number of iterations";
  xytitles,"Number of iteration",swrite(format="Strehl@%dnm [%%]",long(lambda)),[-0.012,0.008];

  y0 = y1 = 17.; dy=6;
  plg,[y1,y1],[0.,2],width=20,color=torgb(tokyonight(3)); y1-=dy;
  plg,[y1,y1],[0.,2],width=20,color=torgb(tokyonight(2)); y1-=dy;
  plg,[y1,y1],[0.,2],width=20,color=torgb(tokyonight(1)); y1-=dy;
  plt,"Initial NCPA",5.5,y0-1,tosys=1,height=12,color=torgb(tokyonight(3)),justify="LA"; y0-=dy;
  plt,"Fitted w/ all optics",5.5,y0-1,tosys=1,height=12,color=torgb(tokyonight(2)),justify="LA"; y0-=dy;
  plt,"Projected on DMs",5.5,y0-1,tosys=1,height=12,color=torgb(tokyonight(1)),justify="LA"; y0-=dy;
  plmargin; range,0;

  write,format="%d rejected runs out of %d for this batch\n",rejected,nvalid;

  if (name) {
    pngcrop,name+"/s_vs_nit";
    write,format="Plot saved in %s/s_vs_nit.png\n",name;
  }

*/


  if (cw!=-1) window,cw;
}


// The following functions are mavis_pray wrappers to run it repeatedly
// or versus some changing parameters (number of modes etc)
func perfvsnit(nitv,ngrid,deltafoc,flux,ron,rseed=,disp=)
{
  error,"needs fixing";
  nitv = _(0,nitv);
  stvsnit = nitv*0.+1;
  st_start_vsnit = st_start_spatial_rms_vsnit = nitv*0.+1;
  st_end_vsnit = st_end_spatial_rms_vsnit = nitv*0.+1;
  st_proj_vsnit = st_proj_spatial_rms_vsnit = nitv*0.+1;
  allstrehls = [];
  for (nn=2;nn<=nof(nitv);nn++) {
    strehls = mavis_pray(,ngrid,deltafoc,flux,ron,init_strehlv,disp=disp,maxiter=nitv(nn),rseed=rseed);
    grow,allstrehls,strehls(,,-);
    st_start_vsnit(nn) = strehls(avg,1);
    st_start_spatial_rms_vsnit(nn) = strehls(rms,1);
    st_corr_vsnit(nn) = strehls(avg,2);
    st_corr_spatial_rms_vsnit(nn) = strehls(rms,2);
    st_end_vsnit(nn) = strehls(avg,3);
    st_end_spatial_rms_vsnit(nn) = strehls(rms,3);
  }
  // fill in the uncorrected case
  st_start_vsnit(1) = st_start_vsnit(2);
  st_start_spatial_rms_vsnit(1) = st_start_spatial_rms_vsnit(2);
  st_corr_vsnit(1) = st_start_vsnit(2);
  st_corr_spatial_rms_vsnit(1) = st_start_spatial_rms_vsnit(2);
  st_end_vsnit(1) = st_start_vsnit(2);
  st_end_spatial_rms_vsnit(1) = st_start_spatial_rms_vsnit(2);

  // cw = current_window();
  // if (window_exists(5)) window,5;
  // else window,5,wait=1,dpi=long(dpi_target_small);
  // errnm = lambda/2/pi*sqrt(-log(st_end_vsnit));
  // fma;
  // plg,errnm,nitv,width=3;
  // plp,errnm,nitv,symbol="o",size=0.5,width=3;
  // plt,swrite(format="%.1f",errnm(0)),nitv(0)+2,errnm(0)+2,justify="LB",tosys=1,height=10;
  // plmargin; range,0
  // xytitles,"Number of iterations","Phase error [nm]",[-0.015,0.];
  // pltitle,swrite(format="Pray: %dx%d grid, efd=%s, flux=%.0f, RON=%g",\
  //     ngrid,ngrid,print(deltafoc)(1),flux*1.,ron*1.);
  // window,cw,wait=1;
  // pause,50;
  return [nitv,st_start_vsnit,st_start_spatial_rms_vsnit,st_end_vsnit,\
    st_end_spatial_rms_vsnit,st_proj_vsnit,st_proj_spatial_rms_vsnit];
}

func do_stats_vs_nit(nitv,nsamp,ngrid,deltafoc,flux,ron,rseed,disp=)
{
  allststart = allstend = allstproj = array(_(0.,nitv*0.),nsamp);
  random_seed,rseed;
  for (k=1;k<=nsamp;k++) {
    write,format="\n\033[32mdo_stats_vs_nit() iteration: %d/%d\033[0m\n",k,nsamp;
    res = perfvsnit(nitv,ngrid,deltafoc,flux,ron,rseed=random(),disp=disp);
    allststart(,k) = res(,2); // start strehls vs nit
    allstend(,k) = res(,4); // end strehls vs nit
    allstproj(,k) = res(,6); // end strehls vs nit
    nitv2 = res(,1); // nit vector, including 0
  }
  if (window_exists(3)) window,3;
  else window,3,wait=1,dpi=long(dpi_target_small);
  fma;
  w = where((abs(allstend(0,)-median(allstend(0,)))<3*allstend(0,rms)));
  write,format="Kept %d out of %d samples\n",nof(w),nsamp;
  ststartavg = allststart(,w)(,avg); stendavg = allstend(,w)(,avg); stprojavg = allstproj(,w)(,avg);
  ststartrms = allststart(,w)(,rms); stendrms = allstend(,w)(,rms); stprojrms = allstproj(,w)(,rms);
  // write,format="Final Strehl [nm] = %.2fnm +/- %.2f\n",erravg(0),errrms(0);
  fma;

  plg,ststartavg,nitv2;
  plp,ststartavg,nitv2,symbol="o",size=0.5;
  pleb,ststartavg,nitv2,dy=ststartrms;

  plg,stprojavg,nitv2,color=torgb(tokyonight(2));
  plp,stprojavg,nitv2,symbol="o",size=0.5,color=torgb(tokyonight(2));
  pleb,stprojavg,nitv2,dy=stprojrms,color=torgb(tokyonight(2));

  plg,stendavg,nitv2,color=torgb(tokyonight(1));
  plp,stendavg,nitv2,symbol="o",size=0.5,color=torgb(tokyonight(1));
  pleb,stendavg,nitv2,dy=stendrms,color=torgb(tokyonight(1));

  plmargin; range,0,1;
  xytitles,"Number of iterations",swrite(format="Strehl @ %.0fnm",float(lambda)),[-0.015,0.];
  pltitle,"End Strehl (all optics: red, DMs: green) vs # of iterations";
  // pltitle,swrite(format="Pray: %dx%d grid, efd=%s, flux=%.0f, RON=%g",\
  //   ngrid,ngrid,print(deltafoc)(1),flux*1.,ron*1.);
}

/*
func projection_only(pd,cond,silent=)
{
  require,"projection.i";
  if (!silent) write,format="\033[32mProjecting passive optics shape onto DMs\033[0m (cond=%.3f)\n",cond;
  // this projects and update pd with new coefficients
  csaved = *pd.coeffs;
  pd = project_to_dms(pd,cond=cond);
  // same as above, this fills pd.mircube with phases at foc=0 and rot=0???
  compute_psfs,pd,0,*pd.coeffs,amp1,amp2,nodisp=1,fromscreens=0;
  *pd.truecube = *pd.origcube - *pd.mircube;
  psfs = compute_psfs(pd,0,,amp1,amp2,nodisp=1,fromscreens=1);
  disp_im = build_bigim(psfs,*pd.xpos,*pd.ypos,0);
  window,1; plsys,1; pli,disp_im; redraw;
  *pd.coeffs = csaved;
  // disp_im = build_bigim(psfs,xpos,ypos,variance*0);
  ntarget = nof(*pd.xpos);
  strehlv = array(0.,ntarget);
  airy = roll(abs(fft(*pd.ipupil,1))^2);
  for (i=1;i<=ntarget;i++) {
    strehlv(i) = max(psfs(,,i)/sum(psfs(,,i)))/pd.peak_airy;
  }
  write,format="\033[31mStrehl over FoV after DM projection: avg=%.1f%%\033[0m rms=%.1f%%\n", \
    100*avg(strehlv),100*strehlv(rms);
  end_strehl = [avg(strehlv),strehlv(rms)];
  return end_strehl;
}

func scan_cond(pd)
{
  ac=spanl(0.3,6,10);
  res=ac*0;
  for (i=1;i<=nof(ac);i++) {
    write,format="conditionning number=%.3f -> ",ac(i);
    res(i)=projection_only(pd,ac(i),silent=(i!=1))(1);
  }
  window,4;
  plot,res,ac;
  logxy,1,0;
  wbest = wheremax(res)
  write,format="Best performance = %.1f%% at conditioning number = %.3f\n",res(wbest)*100,ac(wbest);
  return ac(wbest);
}
*/
