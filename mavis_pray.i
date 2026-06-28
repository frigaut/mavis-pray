/* mavis_pray.i
example of call:
Overall wrappers:
case=10; random_seed,0.45; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=100,modes="zer")
case=10; random_seed,0.75; do_stats,[50],20,8,[0.,-1.5,1.5],1e5,0,disp=1;

Configuration is defined in mavis_pray_config.i. Edit it and run
functions as above (e.g. mavis_pray()).


Example of a run (also reference case):
> case=10; random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=100,modes="zer")
1.25
T=0.000s -> initialising windows
T=0.854s -> initialising config
T=0.854s -> initialising image centring
T=0.854s -> initialising target positions
T=0.854s -> configuration printout
--------------------------------------------------------------------------
Pray for MAVIS (mavis_pray and pray) - System Configuration
Number of optics: 10, Number of extra-focal pos.: 3, number of rotation: 2
Extra focal distances: [0,-1.5,1.5]
Optics conjug. altitude: [45500,13600,6000,1200,0,-1900,-4000,-12400,-23900,-29900]
Optics WFE [nm]: [10,25,25,10,25,47,9,11,6.9,48]
Number of modes (ZER) per optics: [100,100,100,100,100,100,100,100,100,100]
Fit optics?  : [1,1,1,1,1,1,1,1,1,1]
Active optics: [0,1,1,0,1,0,0,0,0,0]
Init phase perturbation from Power spectrum (slope=-2.500)
Optics rotation, config 1: [0,0,0,0,0,0,0,0,0,0]
Optics rotation, config 2: [180,180,180,180,180,180,180,90,0,0]
Number of sources = 8 across 30" FoV, 67 total, hexagonal geometry
source flux = 100000.0, RON= 1
Image centring: init:OFF  pray:OFF
Projection to DM(s) conditionning number: 1.2
--------------------------------------------------------------------------
T=0.854s -> initialising defs
T=13.276s -> initialising masks
T=13.296s -> initialising perturbations
T=13.319s -> Removing average focus
T=13.319s -> Strehl normalisation to 41.0% (10 iterations)
T=13.512s -> Getting high order residuals
Strehl over FoV (rot=[0,0,0,0,0,0,0,0,0,0]): avg=99.3% rms=0.1%
Strehl over FoV (rot=[180,180,180,180,180,180,180,90,0,0]): avg=99.3% rms=0.1%
Strehl over FoV (all rotations): avg=99.3% rms=0.1%
Strehl over FoV from high orders (fitting): avg=99.3% rms=0.1%
T=13.842s -> initialising images
Strehl over FoV (rot=[0,0,0,0,0,0,0,0,0,0]): avg=44.4% rms=15.6%
Strehl over FoV (rot=[180,180,180,180,180,180,180,90,0,0]): avg=37.5% rms=14.5%
Strehl over FoV (all rotations): avg=41.0% rms=15.4%
T=13.999s -> calling pray
Using coeff = 0 as first guess
Using the VMLMB method for minimization.
 ITER    EVAL     CPU [s]            FUNC             max(|G|)   STEPLEN
------  ------  ----------  -----------------------  ---------  ---------
# Iter.   Time (ms)    Eval. Reject.       Obj. Func.           Grad.       Step
# ---------------------------------------------------------------------------------
      0    1519.584       1       0   4.280392221656976e+10   1.461e+10   0.000e+00
     10   18743.320      12       0   2.692591830813855e+09   7.446e+09   1.000e+00
     20   35639.399      23       0   1.254362812538604e+09   3.295e+09   1.000e+00
     30   54099.775      35       0   6.395596149373661e+08   1.455e+09   1.000e+00
     40   72538.582      47       0   4.188744127318403e+08   1.124e+09   1.000e+00
     50   89661.632      58       0   3.031422851147903e+08   9.245e+08   1.000e+00
     60  108047.025      70       0   2.392147550769584e+08   4.854e+08   1.000e+00
     70  124838.690      81       0   1.854182988855661e+08   5.931e+08   1.000e+00
     80  144530.385      94       0   1.575390743933755e+08   4.256e+08   1.000e+00
     90  164363.523     107       0   1.246432861226766e+08   5.142e+08   5.541e-02
    100  181171.796     118       0   1.078448037931315e+08   2.821e+08   1.000e+00
# Termination: too many iterations
Strehl over FoV after compensation: avg=98.6% rms=0.3%
Strehl due to fit=0 optics only (full original rms): 100.0%
Projecting passive optics shape onto DMs
Strehl over FoV after DM projection: avg=76.5% rms=12.1%
> 
*/

/* Project solution onto a different star asterism:
case=10; random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=120,modes="zer"); 
coeff_offsets=*pray_data.coeffs; random_seed,0.75; fovshape="round"; geometry="square"; 
res=mavis_pray(coeff_offsets,16,[0.,-1.5,1.5],100000,1,,disp=1,noinc=1,maxiter=50,modes="zer"); 
This runs mavis_pray on the mavis geometry (or whatever is define din mavis_pray_conf.i),
then re-run by modifying the asterism (e.g. square with 16x16 here) and the fovshape.
Net result shows that perf doesn't have holes between the mavis stars:
Mavis asterism:
Strehl over FoV after compensation (all rotations): avg=98.1% rms=0.5%
Strehl over FoV after DM projection (all rotations): avg=70.9% rms=9.8%
16x16 square, round fov with same 3DM NCPA offsets:
Strehl over FoV after compensation (all rotations): avg=96.0% rms=0.5%
Strehl over FoV after DM projection (all rotations): avg=70.7% rms=8.0%
 */

require,"mavis_pray.h";
require,"mavis_pray_init.i";
require,"mavis_pray_lib.i";
require,"pray.i";

write,format="\n%s\n%s\n","Try running","case=10; random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=50,modes=\"zer\")";
write,format="(copied to cliboard) or\n%s\n","case=10; random_seed,0.75; do_stats,[50],20,8,[0.,-1.5,1.5],1e5,0,disp=1";
system,"echo 'case=10; random_seed,0.75; res=mavis_pray(,8,[0.,-1.5,1.5],100000,1,,disp=1,maxiter=50,modes=\"zer\")' | wl-copy";

func mavis_pray(coeff_offsets,ngrid,deltafoc,flux,ron,&strehlv,disp=,maxiter=,\
  rseed=,verbose=,noinc=,modes=,skip_proj=,proj_method=,config=)
/* DOCUMENT func mavis_pray(coeff_offsets,ngrid,deltafoc,flux,ron,&strehlv,disp=,maxiter=,\
rseed=,verbose=,noinc=,modes=,skip_proj=,proj_method=,config=)
 * mavis_pray simulates the MAVIS NCPA process. It does:
 * initialisation of the optics surface for all optics
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
 *   project_to_dms(pray_data)
 * proj_method= "simple" or "optimal", passed to project_to_dms() as
 *   method=. Defaults to projection_method (set in mavis_pray_conf.i,
 *   currently "simple") when not given here; see project_to_dms() and
 *   optimal_project_modal() for what each method does and for the
 *   dm_proj_cond/dm_proj_tikhonov regularization parameters.
 * config= config file. Default to mavis_pray_conf.i if not given
*/
{
  extern pray_data;
  extern last_random_seed; // in case, to be able to repeat this random realisation
  extern initialisation_timing, minimisation_timing, projection_timing;

  config_file = (config?config:"mavis_pray_conf.i");
  if (!noinc) include,config_file,1;
  tic,2;
  
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
  if (proj_method==[]) proj_method = projection_method;
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
  if (verbose&&disp) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"windows";
  if (disp) status = init_windows();

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

  // define and fill pray data structure
  pray_data         = pray_struct();
  pray_data.teldiam = teldiam;
  pray_data.cobs    = cobs;
  pray_data.pupd    = pupd;
  pray_data.size    = size;
  pray_data.ngrid   = ngrid;
  pray_data.centre  = centre;
  pray_data.nmod    = &nmod;
  pray_data.alt     = &alt;
  pray_data.active  = &active;
  pray_data.patch_diam = &(alt*0.);
  pray_data.config    = &config;
  pray_data.rotv      = &rotv;
  pray_data.disp      = (disp?1:0);
  pray_data.fovshape  = fovshape;
  pray_data.fullfield = fullfield;

  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"image centring";
  status = init_image_centring(pray_data);

  // target positions
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"target positions";
  status = init_target_positions(geometry,fullfield,ngrid,gridpad,fovshape,xpos,ypos);
  pray_data.xpos = &xpos; pray_data.ypos = &ypos;
  ntarget = nof(xpos);

  // configuration printout
  if (verbose) write,format="T=%.3fs -> \033[32mconfiguration printout\033[0m\n",tac();
  status = configuration_printout();


  // Init ipupil
  pray_data.ipupil = &float(make_pupil(size,pupd,xc=centre,yc=centre,cobs=cobs));

  // init defs, etc used by pray, fills pray_data
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"defs";
  status = init_defs(pray_data);

  // init masks (valid phase points at each optics)
  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"masks";
  status = init_masks(pray_data);

  // Airy pattern for Strehl estimation.
  extern workspace;
  workspace = fft_setup(dimsof(*pray_data.ipupil));
  airy = roll(abs(fft(*pray_data.ipupil,1,setup=workspace))^2);
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
    strehlv_ho = get_high_order_residuals(pray_data,config);
  }

  if (verbose) write,format="T=%.3fs -> \033[32minitialising %s\033[0m\n",tac(),"images";
  strehlv_init = init_images(pray_data,config,object,start_strehl,label="Init images: ");

  if (debug&&(imask_radius_scaling!=[])) write,format="T=%.3fs -> %s\n",tac(),"Masking images to cut high frequencies";
  status = mask_images(pray_data,imask_radius_scaling);

  if (stop_after_init_images) return strehlv;

  //**************************************************************************
  // Call pray, which does the minimisation (with calls to vmlmb + pray_error)
  //**************************************************************************

  initialisation_timing = tac(2); tic,2;
  if (coeff_offsets!=[]) {
    if (verbose) write,format="%s\n","\033[32mcoeff_offsets are set, skipping call to pray()\033[0m";
    res = coeff_offsets;
  } else {
    if (verbose) write,format="T=%.3fs -> \033[32mcalling %s\033[0m\n",tac(),"pray";
    res = pray(*pray_data.images,pray_data,deltafoc,variance,object,disp=disp,verbose=verbose,\
      threshold=threshold,nbiter=maxiter);
  }
  minimisation_timing = tac(2); tic,2;
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

  // we compute Strehls for *all* rotation configurations:
  psfs = []; bigim_done = 0;
  for (i=1;i<=nfoc;i++) {
    if (deltafoc(i)!=0) continue; // Strehl has meaning only for in-focus images
    grow,psfs,compute_psfs(pray_data,deltafoc(i),res*0,amp1,amp2, \
      rotv=rotv(,config.roti(i)),nodisp=1,fromscreens=1);
    // big im display only for first rot config:
    if (bigim_done==0) disp_im = build_bigim(psfs,xpos,ypos,variance*0);
    bigim_done = 1;
  }

  strehlv_corr = array(0.,dimsof(psfs)(0));
  for (i=1;i<=nof(strehlv_corr);i++) {
    strehlv_corr(i) = max(psfs(,,i)/sum(psfs(,,i)))/pray_data.peak_airy;
  }
  write,format="\033[31mStrehl over FoV after compensation (all rotations): avg=%.1f%%\033[0m rms=%.1f%%\n", \
    100*avg(strehlv_corr),100*strehlv_corr(rms);
  end_strehl = [avg(strehlv_corr),strehlv_corr(rms)];

  if (allof(active)) {
    write,format="%s\n","All optic active flags are set to one, no fitting to do.";
    strehlv_end = strehlv_corr;
  } else if (skip_proj) {
    write,format="%s\n","skip_proj is set: run project_to_dms(pray_data) later to project.";
    strehlv_end = strehlv_corr*0;
  } else {
    strehlv_end = project_to_dms(pray_data,method=proj_method);
  }
  
  projection_timing = tac(2);
  write,format="Initialisation took %5.1fs\n",initialisation_timing;
  if (!currentiter) error,"currentiter not defined or nill";
  write,format="Minimisation   took %5.1fs (%.0fms/it)\n",minimisation_timing,\
    minimisation_timing*1000/currentiter;
  write,format="Projection     took %5.1fs\n",projection_timing;


  return [strehlv_init,strehlv_corr,strehlv_end,strehlv_ho];
}
// END OF MAVIS_PRAY


func project_to_dms(pd,&phase_rms_max,method=,cond=,tikhonov=,report=,reset=)
/* DOCUMENT project_to_dms(pd,method=,cond=,tikhonov=,report=,reset=)
 * Projects the passive (non-fitted) optics shape onto the active DMs and
 * reports the resulting Strehl. pd must be a pray_struct as produced and
 * filled by mavis_pray(), e.g. global pray_data. Can be called either from
 * within mavis_pray() (when skip_proj is not set) or standalone afterwards
 * (e.g. when mavis_pray() was called with skip_proj=1).
 *
 * method= "optimal": joint multi-DM least-squares modal projection across
 *   all active optics simultaneously (optimal_project_modal()). Generally
 *   better FoV-averaged Strehl than the alternative below, at the cost of
 *   higher stroke demand on individual DMs -- see cond=/tikhonov=.
 *   "simple" (default here): each passive optic projects onto its single
 *   nearest active optic (simple_project(), real-space). Lower stroke
 *   demand, lower Strehl; kept available for comparison. Note that
 *   mavis_pray() itself defaults to projection_method from
 *   mavis_pray_conf.i (currently "simple") rather than this function's own
 *   default, via its proj_method= keyword.
 * cond=, tikhonov= passed through to optimal_project_modal() when
 *   method="optimal" (ignored otherwise); see optimal_project_modal() for
 *   their defaults (dm_proj_cond/dm_proj_tikhonov, or 15/1).
 * report= if set, print per-active-optic phase RMS/peak on mircube after
 *   projection (the physical DM stroke demand) -- works for both methods,
 *   so directly comparable. method="optimal" additionally prints its own
 *   per-DM modal coefficient RMS/peak from optimal_project_modal().
 * reset= projection mutates pd.coeffs and/or pd.mircube, so repeated calls
 *   on the same pd compound rather than compare cleanly. A pristine
 *   pre-projection snapshot of coeffs/mircube is saved automatically (on
 *   the global pray_data) the first time this is called for a given run;
 *   pass reset=1 to restore that snapshot before projecting, so you can
 *   A/B different method=/cond=/tikhonov= from the same starting point.
 *   A fresh mavis_pray() call (which rebuilds pray_data from scratch)
 *   clears the snapshot.
 */
{
  extern pray_data;
  require,"projection.i";

  if (method==[]) method = "optimal";

  if (*pray_data._coeffs0==[]) {
    tmp = (*pd.coeffs)*1;  pray_data._coeffs0  = &tmp;
    tmp = (*pd.mircube)*1; pray_data._mircube0 = &tmp;
  }
  if (reset) {
    *pd.coeffs  = *pray_data._coeffs0;
    *pd.mircube = *pray_data._mircube0;
  }

  deltafoc = *pd.deltafoc;
  nfoc     = nof(deltafoc);
  config   = *pd.config;
  rotv     = *pd.rotv;
  ngrid    = pd.ngrid;
  ntarget  = nof(*pd.xpos);
  disp     = pd.disp;

  if (method=="optimal") {
    write,format="%s\n","\033[32mProjecting passive optics shape onto DMs (optimal multi-DM modal fit)\033[0m";
    pd = optimal_project_modal(pd,cond=cond,tikhonov=tikhonov,report=report);
    compute_psfs,pd,0,*pd.coeffs,amp1,amp2,nodisp=1,fromscreens=0;
  } else {
    write,format="%s\n","\033[32mProjecting passive optics shape onto DMs\033[0m";
    compute_psfs,pd,0,*pd.coeffs,amp1,amp2,nodisp=1,fromscreens=0;
    // now pd.mircube should contains the estimated phase at foc=0.
    // subtract estimate from original phases.
    // The PSFs are computed from "truecube" (when fromscreens==1)
    pd = simple_project(pd);
  }
  *pd.truecube = *pd.origcube - *pd.mircube;

  if (report) {
    // phase (not coefficient) RMS/peak: comparable across both methods,
    // and the physically meaningful "DM stroke" quantity in either case.
    write,format="%s\n","Per-active-optic mircube (phase) RMS/peak after projection:";
    wact = where(*pd.active);
    phase_rms_max = array(0.,[2,nof(wact),2]);
    for (i=1;i<=nof(wact);i++) {
      no = wact(i);
      wm = where((*pd.pupil)(,,no));
      ph = ((*pd.mircube)(,,no))(*)(wm);
      phase_rms_max(i,) = [ph(rms),max(abs(ph))];
      write,format="  optic %d: rms=%.4g peak=%.4g\n",no,ph(rms),max(abs(ph));
    }
  }

  // we compute Strehls for *all* rotation configurations:
  zerocoeffs = (*pd.coeffs)*0;
  psfs = []; bigim_done = 0;
  for (i=1;i<=nfoc;i++) {
    if (deltafoc(i)!=0) continue; // Strehl has meaning only for in-focus images
    grow,psfs,compute_psfs(pd,deltafoc(i),zerocoeffs,amp1,amp2, \
      rotv=rotv(,config.roti(i)),nodisp=1,fromscreens=1);
    // big im display only for first rot config:
    if (bigim_done==0) disp_im = build_bigim(psfs,*pd.xpos,*pd.ypos,0);
    bigim_done = 1;
  }

  if (disp) {
    window, 1;
    plsys, 1;
    pli,disp_im;
    pltitle_vp, "Final";
    redraw;
  }
  strehlv = array(0.,dimsof(psfs)(0));
  for (i=1;i<=nof(strehlv);i++) {
    strehlv(i) = max(psfs(,,i)/sum(psfs(,,i)))/pd.peak_airy;
  }
  end_strehlv = strehlv;
  write,format="\033[31mStrehl over FoV after DM projection (all rotations): avg=%.1f%%\033[0m rms=%.1f%%\n", \
  100*avg(strehlv),100*strehlv(rms);
  end_strehl = [avg(strehlv),strehlv(rms)];

  window,3; fma;
  plot_strehl_contours,strehlv(1:ntarget),*pd.xpos,*pd.ypos,ngrid,\
    fovshape=pd.fovshape,fullfield=pd.fullfield,label="Final NCPA-compensated Strehl across FoV";

  return end_strehlv;
}


func do_stats(nitv,nsamp,ngrid,deltafoc,flux,ron,rseed=,disp=,batchname=,config=)
{
  if (batchname==[]) batchname=""; else batchname=batchname+"_";
  d=timestamp(); name = "do_stats_"+batchname+streplace(d,strfind(" ",d,n=10),"-");
  system,"mkdir "+name;
  config_file = (config?config:"mavis_pray_conf.i");
  system,"cp -p "+config_file+" "+name+"/";

  strehl_start = array(allstrehl_st(),nof(nitv)*nsamp);
  strehl_corr  = array(allstrehl_st(),nof(nitv)*nsamp);
  strehl_end   = array(allstrehl_st(),nof(nitv)*nsamp);
  strehl_ho    = array(allstrehl_st(),nof(nitv)*nsamp);
  ind = 1;
  rejected = 0;
  for (ns=1;ns<=nsamp;ns++) {
    write,format="\n\033[32mSample %d/%d\033[0m\n\n",ns,nsamp;
    for (nn=1;nn<=nof(nitv);nn++) {
      strehls = mavis_pray(,ngrid,deltafoc,flux,ron,init_strehlv,disp=disp,maxiter=nitv(nn),rseed=rseed,config=config);
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
      strehl_ho(ind).nit = nitv(nn);
      strehl_ho(ind).nsamp = ns;
      strehl_ho(ind).strehls = &strehls(,4);
      ind = ind+1;
    }
    plot_do_stats,strehl_start,strehl_corr,strehl_end,strehl_ho,rejected,disp=disp;
  }
  strehl_start = strehl_start(1:ind-1);
  strehl_corr  = strehl_corr(1:ind-1);
  strehl_end   = strehl_end(1:ind-1);
  write,format="Rejected runs: %d\n",rejected;
  write,format="Saving data in folder %s/\n",name;
  if (case==[]) case=-1;
  f = createb(name+"/do_stats.dat");
  save,f,strehl_start,strehl_corr,strehl_end,strehl_ho,lambda,case,xpos,ypos,ngrid,deltafoc,flux,ron,rejected;
  close,f;
}

func plot_do_stats(strehl_start,strehl_corr,strehl_end,strehl_ho,rejected,binsize=,name=,disp=)
{
  if (binsize==[]) binsize=1.0;
  if (name) {
    f = openb(name+"/do_stats.dat");
    restore,f,strehl_start,strehl_corr,strehl_end,strehl_ho,lambda,case,xpos,ypos,ngrid,deltafoc,flux,ron,rejected;
    close,f;
  }

  conf_file = findfiles(name+"/*conf*.i");
  include,conf_file(1),1; // for fovshape and fullfield.

  // Strehl histograms only for nit = max(nitv)
  if (max(strehl_start.nit)==0) return;
  w = wheremax(strehl_start.nit);
  nvalid = sum(strehl_start.nit!=0);
  nit = max(strehl_start.nit);
  ss = sc = se = sho = se4plc = [];
  ssin = scin = sein = shoin = []; // strehl inner FoV (30" diameter)

  // find points that are inside 30" diameter FoV:
  tmp = abs(xpos,ypos)<=1.03*max(abs(xpos));
  // above: 1.03 arbitrary, but otherwise for even ngrid, no edge point is included at all.
  // and replicate as we have accumulated strehl vector over different rotations:
  nrot = nof(*(strehl_start(1).strehls))/nof(tmp);
  in = [];
  for (i=1;i<=nrot;i++) grow,in,tmp;
  win = where(in);

  // compute overall metrics ad fill arrays for plots
  ssavg = ssrms = scavg = scrms = seavg = serms = shoavg = shorms = array(0.,nof(w));
  ssinavg = ssinrms = scinavg = scinrms = seinavg = seinrms = shoinavg = shoinrms = array(0.,nof(w));
  for (i=1;i<=nof(w);i++) {
    grow,ss,*(strehl_start(w(i)).strehls)*100;
    grow,sho,*(strehl_ho(w(i)).strehls)*100;
    grow,sc,*(strehl_corr(w(i)).strehls)*100;
    grow,se,*(strehl_end(w(i)).strehls)*100;
    grow,se4plc,(*(strehl_end(w(i)).strehls))(1:nof(xpos))(,-);
    grow,ssin,(*(strehl_start(w(i)).strehls))(win)*100;
    grow,shoin,(*(strehl_ho(w(i)).strehls))(win)*100;
    grow,scin,(*(strehl_corr(w(i)).strehls))(win)*100;
    grow,sein,(*(strehl_end(w(i)).strehls))(win)*100;

    ssavg(i)    = avg(*(strehl_start(w(i)).strehls)*100);
    shoavg(i)   = avg(*(strehl_ho(w(i)).strehls)*100);
    scavg(i)    = avg(*(strehl_corr(w(i)).strehls)*100);
    seavg(i)    = avg(*(strehl_end(w(i)).strehls)*100);
    ssinavg(i)  = avg((*(strehl_start(w(i)).strehls))(win)*100);
    shoinavg(i) = avg((*(strehl_ho(w(i)).strehls))(win)*100);
    scinavg(i)  = avg((*(strehl_corr(w(i)).strehls))(win)*100);
    seinavg(i)  = avg((*(strehl_end(w(i)).strehls))(win)*100);

    ssrms(i)    = (*(strehl_start(w(i)).strehls)*100)(rms);
    shorms(i)   = (*(strehl_ho(w(i)).strehls)*100)(rms);
    scrms(i)    = (*(strehl_corr(w(i)).strehls)*100)(rms);
    serms(i)    = (*(strehl_end(w(i)).strehls)*100)(rms);
    ssinrms(i)  = ((*(strehl_start(w(i)).strehls))(win)*100)(rms);
    shoinrms(i) = ((*(strehl_ho(w(i)).strehls))(win)*100)(rms);
    scinrms(i)  = ((*(strehl_corr(w(i)).strehls))(win)*100)(rms);
    seinrms(i)  = ((*(strehl_end(w(i)).strehls))(win)*100)(rms);
  }

  ssavg    = ssavg(avg);    ssrms    = ssrms(avg);
  shoavg   = shoavg(avg);   shorms   = shorms(avg);
  scavg    = scavg(avg);    scrms    = scrms(avg);
  seavg    = seavg(avg);    serms    = serms(avg);
  ssinavg  = ssinavg(avg);  ssinrms  = ssinrms(avg);
  shoinavg = shoinavg(avg); shoinrms = shoinrms(avg);
  scinavg  = scinavg(avg);  scinrms  = scinrms(avg);
  seinavg  = seinavg(avg);  seinrms  = seinrms(avg);
  
  write,format="%s\n","Average spatial Strehl, and average of the spatial RMS:";
  write,format="Initial    AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",ssavg,ssrms;
  write,format="High Order AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",shoavg,shorms;
  write,format="Fitted     AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",scavg,scrms;
  write,format="Projected  AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",seavg,serms;
  write,format="%s\n","Same within the circular Fov (MAVIS specs):";
  write,format="Initial    AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",ssinavg,ssinrms;
  write,format="High Order AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",shoinavg,shoinrms;
  write,format="Fitted     AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",scinavg,scinrms;
  write,format="Projected  AVG = %.1f%%, AVG spatial RMS = %.1f%%\n",seinavg,seinrms;

  if (disp) {
    cw = current_window();
    if (window_exists(10)==0) window,10,wait=1,style="clean.gs";
    else window,10;
    fma; limits,square=0; limits;
    if (colors==[]) colors = tokyonight;
    // start strehl (whole fov)
    hy = myhisto2(ss,hx,binsize=binsize);
    plh,hy,hx,color=torgb(colors(3)),width=2,type=3;
    // corrected strehl (whole fov)
    hy = myhisto2(sc,hx,binsize=binsize);
    plh,hy,hx,color=torgb(colors(2)),width=2,type=3;
    // end/fitted strehl (whole fov)
    hy = myhisto2(se,hx,binsize=binsize);
    plh,hy,hx,color=torgb(colors(1)),width=2,type=3;
    // high order (disk)
    hy = myhisto2(shoin,hx,binsize=binsize);
    plh,hy,hx,color=torgb(colors(5)),width=1;
    // start strehl (disk)
    hy = myhisto2(ssin,hx,binsize=binsize);
    plh,hy,hx,color=torgb(colors(3)),width=3;
    // corrected strehl (disk)
    hy = myhisto2(scin,hx,binsize=binsize);
    plh,hy,hx,color=torgb(colors(2)),width=3;
    // end/fitted strehl (disk)
    hy = myhisto2(sein,hx,binsize=binsize);
    plh,hy,hx,color=torgb(colors(1)),width=3;
    pltitle,swrite(format="NCPA performance for %d iterations",nit);
    xytitles,swrite(format="Strehl@%dnm",long(lambda)),"Number in bin",[-0.012,0.008];
    y0 = y1 = limits()(4)*0.97; dy=limits()(4)*0.06; y1+=limits()(4)*0.012;
    plg,[y1,y1],[0.,2],width=20,color=torgb(colors(3)); y1-=dy;
    plg,[y1,y1],[0.,2],width=20,color=torgb(colors(5)); y1-=dy;
    plg,[y1,y1],[0.,2],width=20,color=torgb(colors(2)); y1-=dy;
    plg,[y1,y1],[0.,2],width=20,color=torgb(colors(1)); y1-=dy;
    plt,swrite(format="Strehl start median = %.1f%% (inner: %.1f%%)",median(ss),median(ssin)),6,y0,tosys=1,height=10,color=torgb(colors(3)),justify="LA"; y0-=dy;
    plt,swrite(format="Strehl high order median = %.1f%% (inner: %.1f%%)",median(sho),median(shoin)),6,y0,tosys=1,height=10,color=torgb(colors(5)),justify="LA"; y0-=dy;
    plt,swrite(format="Strehl fitted median = %.1f%% (inner: %.1f%%)",median(sc),median(scin)),6,y0,tosys=1,height=10,color=torgb(colors(2)),justify="LA"; y0-=dy;
    plt,swrite(format="Strehl projected median = %.1f%% (inner: %.1f%%)",median(se),median(sein)),6,y0,tosys=1,height=10,color=torgb(colors(1)),justify="LA"; y0-=dy;
    plmargin; range,0;

    if (window_exists(11)==0) {
      window,11,wait=1,style="clean.gs";
      pause,500;
      system,"niri msg action focus-window --id $(niri-get-id-from-title.sh 'Yorick 11')";
      pause,100;
      system,"niri msg action consume-or-expel-window-right";
    } else window,11;
    fma;
    plot_strehl_contours,se4plc(,avg),xpos,ypos,ngrid,\
      label="Final NCPA-compensated Strehl across FoV",fovshape=fovshape,fullfield=fullfield;

    if (name) {
      window,10;
      pause,500;
      pngcrop,name+"/do_stats_hist";
      write,format="Plot saved in %s/do_stats_hist.png\n",name;
      window,11;
      pause,500;
      pngcrop,name+"/do_stats_cont";
      write,format="Plot saved in %s/do_stats_cont.png\n",name;
    }

    if (cw!=-1) window,cw;
  }

  write,format="%d rejected runs out of %d for this batch\n",rejected,nvalid;

  pause,500;

  if (disp&&(cw!=-1)) window,cw;
}

func merge_stats(marker)
{
  strehl_startv=strehl_corrv=strehl_endv=strehl_hov=rejectedv=[];
  af = findfiles("do_stats_"+marker+"_*");
  for (i=1;i<=nof(af);i++) {
    write,format="Attempting to read %s/do_stats.dat... ",af(i);
    f = openb(af(i)+"/do_stats.dat");
    restore,f,strehl_start,strehl_corr,strehl_end,strehl_ho,lambda,case,ngrid,deltafoc,flux,ron,rejected;
    close,f;
    write,format="%s\n","done.";
    grow,strehl_startv,strehl_start;
    grow,strehl_corrv,strehl_corr;
    grow,strehl_endv,strehl_end;
    grow,strehl_hov,strehl_ho;
    grow,rejectedv,rejected;
  }
  strehl_start = strehl_startv;
  strehl_corr = strehl_corrv;
  strehl_end = strehl_endv;
  strehl_ho = strehl_hov;
  rejected = rejectedv;
  dir  = "do_stats_"+marker;
  name = dir+"/do_stats.dat";
  write,format="Writing merged data to %s\n",name;
  system,"mkdir -p "+dir;
  f = createb(name);
  save,f,strehl_start,strehl_corr,strehl_end,strehl_ho,lambda,case,ngrid,deltafoc,flux,ron,rejected;
  close,f;
}
