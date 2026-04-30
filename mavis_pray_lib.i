func configuration_printout(void)
{
  write,format="%s\n","--------------------------------------------------------------------";
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
  xy1 = 1+(size-stride)/2;
  xy2 = 1+long(size/2+(nlin-0.5)*stride);
  bigim = bigim(xy1:xy2,xy1:xy2);
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

func get_high_order_residuals(pd)
{
  nopt    = nof(*pd.alt);
  nz12 = (*pd.nmod)(cum);
  nz1 = (nz12+1)(1:-1);
  nz2 = nz12(2:);

  for (no=1;no<=nopt;no++) {
    mask1 = (*pd.maskcube)(,,no);
    mask2 = (*pd.def)(,,nz1(no)+1) < 100;
    wpatch = where(mask2);
    d = (*pd.def)(,,nz1(no):nz2(no))(*,)(wpatch,);
    d = d(,2:); // first slot is zero
    // add TT and focus that are missing from defs:
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
    plsys,3; pli,(*pd.truecube)(,,no)*mask1; pltitle_vp,"turbulent";
    plsys,2; pli,rec2d*mask1; pltitle_vp,"fitted";
    plsys,1; pli,((*pd.truecube)(,,no)-rec2d)*mask1; pltitle_vp,"Residuals";
    // if (hitReturn()=="s") error;
    pause,500;
    // (*pd.truecube)(,,no) = rec2d;
  }
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
