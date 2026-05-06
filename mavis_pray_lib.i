func configuration_printout(void)
{
  write,format="\n%s\n","--------------------------------------------------------------------";
  write,format="%s\n","Pray for MAVIS (mavis_pray and pray) - System Configuration";
  write,format="Number of optics: %d, Number of extra-focal pos.: %d, number of rotation: %d\n",nopt,nof(deltafoc),dimsof(rotv)(0);
  write,format="%s","Extra focal distances: "; deltafoc_orig;
  write,format="%s","Optics conjug. altitude: "; alt;
  write,format="%s","Optics WFE [nm]: "; nm_rmsv;
  write,format="Number of modes (\033[31m%s\033[0m) per optics: ",strcase(1,usemodes); nmod;
  write,format="%s","Fit optics?: "; fit;
  perturb_from = (initphase=="screens"?"Power spectrum":"Mode coefficients");
  if (initphase=="screens") perturb_from += swrite(format=" (slope=%.3f)",ps_slope);
  write,format="Init phase perturbation from \033[31m%s\033[0m\n",perturb_from;
  for (i=1;i<=nrot;i++) { write,format="Optics rotation, config %d: ",i; rotv(,i); }
  write,format="Number of sources = %d across %.0f\" FoV, %d total, %s geometry\n",ngrid,fullfield,nof(xpos),geometry;
  write,format="source flux = %.1f, RON= %g\n",flux,ron;
  if (imask_radius_scaling!=[]) \
    write,format="Image mask radius = %.0f%%\n",imask_radius_scaling*100;
  write,format="%s\n","Current random seed stored in extern last_random_seed, use as rseed to repeat current";
  write,format="%s\n","--------------------------------------------------------------------";
}

func strehl_normalisation(&pd,&coeff,config,rotv,peak_airy)
{
  w0 = where(deltafoc==0);
  if (nof(w0)==0) \
    error,"For Strehl normalisation, you have to have one of the deltafoc=0";
  ima = [];
  for (i=1;i<=nof(w0);i++) {
    grow,ima,compute_psfs(pd,0,coeff,amp1,amp2, \
      rotv=rotv(,config.roti(w0(i))),nodisp=1,fromscreens=(initphase=="screens"));
  }
  ima = ima/ima(*,)(sum,)(-,-,); // normalisation
  strehlavg = avg(ima(*,)(max,)/peak_airy);
  // fact = sqrt(log(strehl_target)/log(strehlavg));
  fact = (log(strehl_target)/log(strehlavg))^0.5;
  if (debug) write,format="Current Strehl average = %.1f%%, scaling by %f\n",100*strehlavg,fact;
  *pd.mircube *= fact;
  *pd.truecube *= fact;
  *pd.truecoeffs *= fact;
  // *pd._def_pup *= fact; << this is definitely not necessary, as def not affected.
  // coeff *= fact; << this was creating the issue I fought with.
  // write,format="%s","\033[31mBEWARE, STREHL NORMALISATION IS BROKEN SOMEHOW (AFFECT RESULTS)\033[0m\n";
  return 0;
}

func build_bigim(data_set,xpos,ypos,variance,noeclat=)
{
  if (zoomfactor==[]) zoomfactor=3.0;
  if (variance==[]) variance = 0.0;
  dims   = dimsof(data_set);
  nim    = dims(0);
  size   = dims(2);

  nlin   = long(sqrt(nim));
  stride = long(size/zoomfactor);
  // dposp  = xpos(2)-xpos(1);
  dposp = fullfield/ngrid;

  xposp  = 1+long((xpos-min(xpos))/dposp*stride);
  yposp  = 1+long((ypos-min(ypos))/dposp*stride);

  bigim = array(0.,[2,max(xposp)+size,max(yposp)+size]);

  for (i=1;i<=nim;i++) {
    if (noeclat) {
      bigim(xposp(i):xposp(i)+size-1,yposp(i):yposp(i)+size-1) += data_set(,,i);
    } else {
      bigim(xposp(i):xposp(i)+size-1,yposp(i):yposp(i)+size-1) += eclat(data_set(,,i));
    }
  }
  // crop the image
  // xy1 = 1+(size-stride)/2;
  // xy2 = 1+long(size/2+(nlin-0.5)*stride);
  // bigim = bigim(xy1:xy2,xy1:xy2);
  bigim += random_normal(dimsof(bigim))*sqrt(variance);
  return bigim;
}

func mask_images(&pd,imask_radius_scaling)
{
  if (imask_radius_scaling==[]) return 0;
  // write,format="%s\n","-> Masking images to cut high frequencies";
  for (i=1;i<=dimsof(*pd.images)(4);i++) {
    for (j=1;j<=dimsof(*pd.images)(5);j++) {
      radius = imask_radius_scaling*(pd.size/2.+abs(deltafoc(j))*pd.size/4.);
      imask = dist(pd.size)<radius;
      (*pd.images)(,,i,j) *= imask;
    }
  }
}

func centre_image(image,pd)
{
  return image;
  cg = (image(,,-)*(*pd.xy4centring))(sum,sum,)/sum(image);
  image = float(fftshift(image,-cg(1),-cg(2)));
  // tv,image; pause,5;
  return image;
}

func get_high_order_residuals(pd,config)
{
  nopt    = nof(*pd.alt);
  nz12 = (*pd.nmod)(cum);
  nz1 = (nz12+1)(1:-1);
  nz2 = nz12(2:);
  hocube = (*pd.truecube)*0;

  for (no=1;no<=nopt;no++) {
    mask1 = (*pd.maskcube)(,,no);
    mask2 = (*pd.def)(,,nz1(no)+1) < 100;
    wpatch = where(mask2);
    d = (*pd.def)(,,nz1(no):nz2(no))(*,)(wpatch,);
    d = d(,2:); // first slot is zero
    // add TT and focus that are missing from defs:
    prepzernike,pd.size,pd.size,pd.centre,pd.centre;
    ttf = [zernike_ext(2),zernike_ext(3),zernike_ext(4)];
    ttf = ttf(*,)(wpatch,);
    d = _(ttf,d);
    ddt = d(+,)*d(+,);
    ddtm1 = LUsolve(ddt);
    dm1 = ddtm1(,+)*d(,+);
    phase = (*pd.truecube)(,,no)(*)(wpatch);
    proj = dm1(,+)*phase(+);
    rec = d(,+)*proj(+);
    rec2d = array(0.0f,[2,pd.size,pd.size]);
    rec2d(wpatch) = rec;
    window,1; fma;
    hocube(,,no) = ((*pd.truecube)(,,no)-rec2d)*mask1;
    plsys,3; pli,(*pd.truecube)(,,no)*mask1; pltitle_vp,"turbulent";
    plsys,2; pli,rec2d*mask1; pltitle_vp,"fitted";
    plsys,1; pli,hocube(,,no); pltitle_vp,"Residuals";
    // if (hitReturn()=="s") error;
    pause,100;
    // (*pd.truecube)(,,no) = rec2d;
  }
  // get strehl that correspond to the high order only:
  cubesave = *pd.truecube; // save the actual cube
  pd.truecube = &hocube;
  strehlv = init_images(pd,config,object,start_strehl);
  pd.truecube = &cubesave;
  hocube = [];
  write,format="\033[31mStrehl over FoV from high orders (fitting): avg=%.1f%%\033[0m rms=%.1f%%\n", \
    100*avg(strehlv),100*strehlv(rms);
}

func fix_diskharmonic(&dh)
{
  w = where(dh(1,1,)==0);
  pup = (abs(dh)(,,w)(,,sum)>0);
  dh = dh*pup(,,-);
  return pup;
}

func make_phase_screens(pup,lambda,nm_rms,slope,rseed=,remove_tt=,remove_foc=)
/* DOCUMENT make_phasescreens(dim,lambda,nm_rms,slope)
 * pup: pupil array (needed for computation of rms and TT removal)
 * lambda: wavelength [nm]
 * nm_rms: rms of output phase [nm]
 * slope: slope of phase power spectrum (usually around -2 for optics)
 *
 * Return the phase [2,dim,dim] in radians
 */
{
	dim = dimsof(pup)(2);
	mtf = roll(clip(dist(dim),1,)^slope);
	// random_seed,(rseed?rseed:0.3);
	pha = random(dimsof(mtf))*2*pi;
	obj = mtf*exp(1i*pha);
	phase = fft(obj,1).re;
	// phase = phase/phase(*)(rms)*nm_rms/lambda*2*pi;
	// remove piston and TT (experimental):
	w = where(pup);
	phase -= phase(w)(avg);
	if (remove_tt) {
		tt = (indices(dim)-dim/2.-0.5);
		tip = sum(phase(w)*tt(,,1)(w))/sum(tt(,,1)(w)^2);
		tilt = sum(phase(w)*tt(,,2)(w))/sum(tt(,,2)(w)^2);
		phase = phase-tip*tt(,,1)-tilt*tt(,,2);
	}
	if (remove_foc) {
		focus = dist(dim,xc=dim/2+1,yc=dim/2+1)^2;
		foc = sum(phase(w)*focus(w))/sum(focus(w)^2);
		phase = phase-foc*focus;
	}
	phase = phase/phase(w)(rms);
	// possibly we could think or re-adding the TT now that the phase
	// has been normalised on the TT less phase. Consider implementing
	phase = phase*nm_rms/lambda*2*pi;
	return phase;
}

func test_make_phase_screens(nm_rms)
{
	dim=256; lambda=550.;
	if (nm_rms==[]) nm_rms=50.;
	pup = dist(dim)<(dim/4.);
	airy = roll(abs(fft(pup,1))^2.);
	peak_airy=max(airy/sum(airy));
	phase=make_phase_screens(pup,lambda,nm_rms,-2.,remove_tt=1,remove_foc=1);
	psf = roll(abs(fft(pup*exp(1i*phase),1))^2.);
	psf = psf/sum(psf);
	write,format="Returned phase rms [nm]: %.1f\n",phase(where(pup))(rms)*lambda/2/pi;
	write,format="Strehl@%.0fnm from image, i.e. max(ima)/max(airy):  %.3f\n",lambda,max(psf)/peak_airy;
	var=phase(where(pup))(*)(rms)^2.;
	write,format="Strehl expected from Marechal approximation: %.3f\n",exp(-var); tv,psf;
	return phase;
}


func plot1(allres,nitv2,nsamp,ngrid,deltafoc,flux,ron,nsig,col=)
{
  erravg = allres(,avg);
  errrms = allres(,rms);
  w = where((abs(allres(0,)-median(allres(0,)))<nsig*allres(0,rms)));
  w2 = where((abs(allres(0,)-median(allres(0,)))>=nsig*allres(0,rms)));
  write,format="Kept %d out of %d samples, rejected %s\n",nof(w),nsamp,print(w2);
  allres = allres(,w)
  erravg = allres(,avg);
  errrms = allres(,rms);
  plg,erravg,nitv2,width=3,color=col;
  plp,erravg,nitv2,symbol="o",size=0.5,width=3,color=col;
  pleb,erravg,nitv2,dy=errrms,width=3,color=col;
  plmargin; range,0;
  xytitles,"Number of iterations","Phase error [nm]",[-0.015,0.];
  pltitle,swrite(format="Pray: %dx%d grid, PDEFD=%s, flux=%.0f, RON=%g",\
    ngrid,ngrid,print(deltafoc)(1),flux*1.,ron*1.);
  return allres;
}

func plot_failed(file)
{
  d = rdcols(file);
  rmsin = *d(1); rmsout = *d(2);
  s = sort(rmsin);
  rmsin = rmsin(s); rmsout = rmsout(s);
  nsamp = nof(rmsin);
  step = 10;
  x = y = [];
  for (i=long(min(rmsin));i<=long(max(rmsin)+step);i=i+step) {
    w = where((rmsin>=i)&(rmsin<(i+step)));
    if (nof(w)==0) continue;
    grow,x,avg(rmsin(w));
    grow,y,100.*sum(rmsout(w)>8)/nof(w);
  }
  fma; logxy,0,0; limits,square=0; limits;
  plh,y,x;
  pltitle,swrite(format="Convergence failure (%d samples)",nsamp);
  xytitles,"Input phase rms error [nm]","Failure to converge [%]",[-0.015,0];
}

func get_non_normalised_strehls(nit)
/* DOCUMENT
 * The goal of this stand alone function is to do some statistics on the
 * images generated by mavis_pray. In particular, check that the actual
 * Strehl generated from the nm_rmsv in the config are compatible with the
 * expected Strehl (from nm_rmsv), and what is the natural strehl distribution.
 * Indeed, Given we generate the optics based on a given rms per optics, these
 * optics may combine in different ways: Sometimes compensate each other
 * (possible also with a field dependency), sometimes make things worse.
 * Here, I generate nit realisations, each for N PSFs across the FoV. At the
 * end, I plot the distribution of Strehl over all images (Nxnit).
 * Turns out there is indeed a very large variation: E.g. for a 78.9nm total
 * expected rms (RSSing the individual optics rms), I get some images with as
 * low a Strehl as 10% and other with as high a Strehl as 90%. The median/average
 * is as expected though, with a value or 45% (44.4% expected).
 * This is very important to know in the case of MAVIS to manage expectations.
 * One could also think to actually use this to rotate the optics to get the
 * best overall Strehl in the FoV. But that seems very complicated to
 * implement during alignment...
 */
{
  // to init config:
  stop_after_init_images=1;
  res = mavis_pray(,8,[0,-1.5,1.5],1000,0.,,disp=1,maxiter=50);
  geometry  = "square"; // "square" or "hexagonal"
  fovshape  = "square";     // "round" if desired if not will default to square
  dmrms = 12.;
  nm_rmsv   = [10.,dmrms,dmrms,30,dmrms,47,9.05,11.43,6.9,48.25];
  nm_rms    = sqrt(sum(nm_rmsv^2));
  expected_strehl = exp(-(2*pi*sqrt(sum(nm_rmsv^2))/lambda)^2);
  strehl_normalise=0;
  // allres = array(0.,[2,2,nit]);
  allstrehl = [];
  for (n=1;n<=nit;n++) {
    write,format="\nget_non_normalised_strehls() iteration %d/%d\n\n",n,nit;
    rs = unix_time(now=1)/1e6+random(); rs -= long(rs);
    res = mavis_pray(,8,[0,-1.5,1.5],1000,0.,,disp=1,maxiter=50,noinc=1);
    grow,allstrehl,res(,-);
    // allres(,n) = res(,1);
    // cw = current_window();
    window,5,wait=1; limits,square=0;
    hy=histo2(100*allstrehl(*),hx,binsize=2.5);
    fma; plh,hy,hx;
    plg,[0.,max(hy)],[1,1]*100*expected_strehl,type=2;
    xytitles,swrite(format="Strehl@%dnm",long(lambda)),"Number in bin",[-0.015,0.];
    pth=pltitle_height; pltitle_height=10;
    // pltitle,swrite(format="Iter %d/%d, Strehl=%.1f%%+/-%.1f%% (expected=%.1f%%)",\
    //   n,nit,100*allstrehl(*)(avg),100*allstrehl(*)(rms),100*expected_strehl);
    pltitle,swrite(format="It %d, Strehl=%.1f%% (med=%.1f) +/-%.1f (expected=%.1f from Prms=%.1f)",\
      n,100*allstrehl(*)(avg),median(100*allstrehl(*)),100*allstrehl(*)(rms),100*expected_strehl,nm_rms);
    plmargin; range,0;
    window,3,wait=1;
    pltitle_height = pth;
  }
  // allres_avg = allres(,avg);
  // allres_rms = allres(,rms);
  write,format="Strehl@%dnm over FoV and %d realisations = %.1f%% (expected %.1f%%) +/- %.1f%%\n",\
    long(lambda),nit,100*avg(allstrehl(*)),100*expected_strehl,100*allstrehl(*)(rms);
  allstrehl_per_frame = allstrehl(avg,);
  write,format="%s\n","FoV-averages Strehl per realisation: ";
  stat,allstrehl_per_frame;
  // allstrehl_per_frame; allstrehl_per_frame(avg); allstrehl_per_frame(rms);
  // write,format="Average of spatial average in the FoV=%.1f%% (expected %.1f%%)\n",100*allres_avg(1),100*expected_strehl;
  // write,format="Average of spatial stdev   in the FoV=%.1f%%\n",100*allres_avg(2);
  // write,format="Stdev   of spatial average in the FoV=%.1f%%\n",100*allres_rms(1);
  // write,format="Stdev   of spatial stdev   in the FoV=%.1f%%\n",100*allres_rms(2);
  // return [allres_avg,allres_rms];
  return allstrehl;
}

func check_rms_to_strehl(nit,nmrms,ps_slope,remove_tt=,remove_foc=)
/* DOCUMENT
> for (nm=0;nm<=120;nm=nm+10) check_rms_to_strehl,500,float(nm),-2.5,remove_foc=1,remove_tt=1;
Average Strehl=100.0%+/-0.0% (set rms=0.0nm, real rms=0.0nm, expected Strehl=100.0)
Average Strehl=98.7%+/-0.0% (set rms=10.0nm, real rms=10.0nm, expected Strehl=98.7)
Average Strehl=94.9%+/-0.0% (set rms=20.0nm, real rms=20.0nm, expected Strehl=94.9)
Average Strehl=88.9%+/-0.1% (set rms=30.0nm, real rms=30.0nm, expected Strehl=88.9)
Average Strehl=81.1%+/-0.2% (set rms=40.0nm, real rms=40.0nm, expected Strehl=81.2)
Average Strehl=72.1%+/-0.4% (set rms=50.0nm, real rms=50.0nm, expected Strehl=72.2)
Average Strehl=62.4%+/-0.7% (set rms=60.0nm, real rms=60.0nm, expected Strehl=62.5)
Average Strehl=52.6%+/-1.0% (set rms=70.0nm, real rms=70.0nm, expected Strehl=52.8)
Average Strehl=43.1%+/-1.8% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
Average Strehl=35.3%+/-3.1% (set rms=90.0nm, real rms=90.0nm, expected Strehl=34.7)
Average Strehl=28.7%+/-4.2% (set rms=100.0nm, real rms=100.0nm, expected Strehl=27.1)
Average Strehl=25.5%+/-5.2% (set rms=110.0nm, real rms=110.0nm, expected Strehl=20.6)
Average Strehl=23.5%+/-5.5% (set rms=120.0nm, real rms=120.0nm, expected Strehl=15.3)

In conclusion: Strehl computed from images follow really well the expected strehl
from the Marechal approximation

Dependency on power spectrum slope:
> check_rms_to_strehl,500,80.,-1.5,remove_foc=1,remove_tt=1;
Average Strehl=43.1%+/-1.1% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
> check_rms_to_strehl,500,80.,-2.0,remove_foc=1,remove_tt=1;
Average Strehl=43.0%+/-1.4% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
> check_rms_to_strehl,500,80.,-2.5,remove_foc=1,remove_tt=1;
Average Strehl=43.3%+/-1.7% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
> check_rms_to_strehl,500,80.,-2.66,remove_foc=1,remove_tt=1;
Average Strehl=43.3%+/-1.8% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
Actual Strehl from images does not depend at first order of power spectrum slope

Dependency on TT or Focus removal:
> check_rms_to_strehl,500,80.,-2.5,remove_foc=1,remove_tt=1;
Average Strehl=43.3%+/-1.9% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
> check_rms_to_strehl,500,80.,-2.5,remove_foc=0,remove_tt=1;
Average Strehl=43.2%+/-1.8% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
> check_rms_to_strehl,500,80.,-2.5,remove_foc=0,remove_tt=0;
Average Strehl=77.6%+/-12.9% (set rms=80.0nm, real rms=80.0nm, expected Strehl=43.4)
Focus removal does not influence final strehl result (remember phase rms normalisation
is done *after* phase manipulation). But of course removing TT does.
So, that means **the marechal approximation is for a TT removed phase.**

*/
{
  dim = 128;
  lambda = 550;
  if (remove_foc==[]) remove_foc=0;
  if (remove_tt==[]) remove_tt=0;
  pup = dist(dim)<(dim/6.);
  airy = roll(abs(fft(pup,1))^2.);
  max_airy = max(airy);
  expected_strehl = exp(-(2*pi*nmrms/lambda)^2.);
  strehlv = pharms = array(0.,nit);
  wpup = where(pup);
  for (n=1;n<=nit;n++) {
    pha = make_phase_screens(pup,lambda,nmrms,ps_slope,remove_tt=remove_tt,remove_foc=remove_foc);
    pharms(n) = (pha(wpup))(rms)*lambda/(2*pi);
    psf = roll(abs(fft(pup*exp(1i*pha),1))^2.);
    tv,psf;
    strehlv(n) = max(psf)/max_airy;
  }
  write,format="Average Strehl=%.1f%%+/-%.1f%% (set rms=%.1fnm, real rms=%.1fnm, expected Strehl=%.1f)\n",
    100*strehlv(avg),100*strehlv(rms),pharms(avg),nmrms,100*expected_strehl;
}
