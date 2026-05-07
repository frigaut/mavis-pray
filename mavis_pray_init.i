func init_windows(void)
/* DOCUMENT
 * init windows
 */
{
  extern init_window_done;
  win = [4,3]

  for (i=1;i<=nof(win);i++) {
    if (!window_exists(win(i))) window,win(i),wait=1;
  }

  pause,100;
  if (!window_exists(1)) plsplit,1,3,win=1,square=1,dpi=200,margin=-0.015,style="nobox.gs",\
    vp=[0.122,0.672,0.12,1.];

  pause,100;
  nw = clip(nopt,4,13); nw = 13;
  if (!window_exists(2)) plsplit,3,nw,win=2,style="nobox.gs",square=1,dpi=190,margin=-0.02,vp=[0.277,0.656,0.14,0.95];
  pause,100;

  if (init_window_done) return 0;

  system,"./layout_windows.sh";

  init_window_done = 1;

  return 0;
}

func init_config(&config)
/* DOCUMENT
 * init config for rotations
 */
{
  config = array(config_struct,nfoc*nrot);
  config.foc = array(deltafoc,nrot)(*); // e.g. [0,-1.5,1.5,0,-1.5,1.5]
  tmp = [];
  for (i=1;i<=nrot;i++) grow,tmp,array(i,nfoc);
  config.roti = tmp; // e.g. [1,1,1,2,2,2]
  return 0;
}

func init_target_positions(geometry,diam,ngrid,gridpad,&xpos,&ypos)
/* DOCUMENT
 * computes x and y position of grid points given geometry,
 * number of points in diameter and diameter
 */
{
  if (strpart(geometry,1:3)=="squ") {
    tmp = (indices(ngrid)-ngrid/2-1);
    tmp += (odd(ngrid)?0:0.5);
    tmp *= diam/ngrid;
    xpos = tmp(,,1)(*);
    ypos = tmp(,,2)(*);
  } else if (strpart(geometry,1:3)=="hex") {
    tmp = (indices(ngrid+4)-(ngrid+4)/2-1.);
    tmp(,,1) += 0.5*(indgen(ngrid+4)%2)(-,);
    tmp(,,1) += (odd(ngrid)?0:0.5);
    tmp(,,2) *= sin(60*pi/180.);
    tmp *= diam/(ngrid-0.8);
    w = where(abs(tmp(,,1))<(diam/2.)&abs(tmp(,,2))<(diam/2.));
    xpos = tmp(,,1)(w);
    ypos = tmp(,,2)(w);
  } else error,"geometry undefined";

  if (fovshape=="round") {
    if (gridpad==[]) gridpad = 0.;
    w = where(abs(xpos,ypos)<=(diam/2.+gridpad));
    xpos = xpos(w);
    ypos = ypos(w);
  }
  return 0;
}

func init_image_centring(&pd)
{
  pd.xy4centring = &(indices(pd.size)-pd.size/2-0.5);
}

func init_defs(&pd,tiptilt=)
/* DOCUMENT init_defs
 * initialise defs and other arrays
*/
{
  require,"projection.i";

  nopt    = nof(*pd.alt);
  ntarget = nof(*pd.xpos);

  // Create mircube
  mircube = array(float,[3,size,size,nopt]);
  pd.mircube = &mircube;

  _p  = pd.pupd;
  _p1 = long(ceil(pd.centre-_p/2.));
  _p2 = long(floor(pd.centre+_p/2.));
  _p  = _p2-_p1+1;

  pd._n  = _p+4;
  pd._n1 = _p1-2; // FIXME why -2?
  pd._n2 = _p2+2;

  // Init def, mode maps for all optics
  def	= array(float,[3,pd.size,pd.size,(*pd.nmod)(sum)]);
  nz12 = (*pd.nmod)(cum);
  nz1 = (nz12+1)(1:-1);
  nz2 = nz12(2:);
  psize = pd.teldiam/pd.pupd;
  cpt = 0;

  // extern focus;
  // FIXME: there is only one array and several nopt? yes, because it is
  // way oversized and the amplitude doesn't really matter?
  // FIXED
  pd.focus = &((sqrt(3.)*dist(pd.size,xc=pd.size/2.+0.5,yc=pd.size/2.+0.5)/(pd.pupd/2.))^2.);

  // create modes per optic
  for (k=1;k<=nopt;k++) {
    patch_diam = long(ceil(pd.pupd+2.*max(abs(*pd.xpos,*pd.ypos))*4.848e-6*(abs(alt(k)))/psize));
    (*pd.patch_diam)(k) = patch_diam;
    patch_diam = long(ceil(patch_diam/2.)*2)+2;

    defdm = generate_modes(usemodes, (*pd.nmod)(k), pd.size, patch_diam, pupil);
    if (usemodes=="dh") defdm *= 5;
    if (usemodes=="kl") defdm *= 4;
    wout = where(pupil==0);
    defdm(*,)(wout,) = 1e6;
    // this returns tip,tilt,focus, etc
    // normalise:
    radeg = array(0.,(*pd.nmod)(k));
    for (i=1;i<=(*pd.nmod)(k);i++) radeg(i) = zernumero(i+1)(1); // +1 because our modes start with tip, not piston
    defdm = defdm*(1./radeg^1.5)(-,-,);
    if (k==1) {
      def(,,2:nz2(1)) = defdm(,,1:-1); // leave space for focus
    } else {
      def(,,nz1(k):nz2(k)) = defdm;
    }
    defdm = [];
  }
  def(,,1) = *pd.focus;

  pd.def = &def;

  // Init dmgsXYposcub : will be needed by get_2d_phase_from_cube
  xref = indgen(pd._n)-(pd._n+1)/2.;
  yref = indgen(pd._n)-(pd._n+1)/2.;

  dmgsxposcub = dmgsyposcub = array(float,[3,pd._n,nopt,ntarget]);

  for (n=1;n<=ntarget;n++) {
    // loop on pseudo-DMs
    for (no=1;no<=nopt;no++) {
      // offsets of the centre of beam on DM NO
      xoff = (*pd.xpos)(n)*4.848e-6*(*pd.alt)(no)/psize;
      yoff = (*pd.ypos)(n)*4.848e-6*(*pd.alt)(no)/psize;
      dmgsxposcub(,no,n) = xref + xoff;
      dmgsyposcub(,no,n) = yref + yoff;
    }
  }

  pd.dmgsxposcub = &dmgsxposcub;
  pd.dmgsyposcub = &dmgsyposcub;

  szdef = dimsof(def);
  _def_pup = array(float,[5,szdef(2),szdef(3),szdef(4),ntarget,nrot]);
  for (n=1;n<=nrot;n++) {
    for (i=1;i<=ntarget;i++) _def_pup(,,,i,n) = get_def_in_pupil_from_dir(pd,i,rotv=rotv(,n));
  }
  pd._def_pup = &_def_pup;
  return 0;
}

func init_masks(&pd)
/* DOCUMENT
 * Computes the "valid" pixels on all optics (i.e. the pixels in the cubes
 * that are seen by the PSF formation process).
 */
{
  nopt    = nof(*pd.alt);
  ntarget = nof(*pd.xpos);
  psize = pd.teldiam/pd.pupd;
  maskcube = array(0.,[3,size,size,nopt]);
  xv = yv = indgen(size);
  if (debug_masks) window,1;
  for (no=1;no<=nopt;no++) {
    for (n=1;n<=ntarget;n++) {
      xoff = (*pd.xpos)(n)*4.848e-6*(*pd.alt)(no)/psize;
      yoff = (*pd.ypos)(n)*4.848e-6*(*pd.alt)(no)/psize;
      maskcube(,,no) += bilinear(*pd.ipupil,xv+xoff,yv+yoff,grid=1);
    }
    maskcube(,,no) = maskcube(,,no)>0;
    if (debug_masks) { tv,maskcube(,,no); hitReturn; }
  }
  pd.maskcube = &maskcube;
  return 0;
}

func init_perturbation(&pd,&coeff,&cmin,&cmax)
/* DOCUMENT
 * Generate coeffs or phase screens and fill pray_data cubes
 * Also sets the cmin and cmax (parameter bounds) for vmlmb
 */
{
  nopt = nof(*pd.alt);
  coeff = cmax = cmin = [];
  nz12 = nmod(cum);
  // var = 0.;
  for (no=1;no<=nopt;no++) {
    c = weight(no)*(random(nmod(no))-0.5)/sqrt(indgen(nmod(no)));
    c(1) = 0.; // FIXME
    cmx = 20*weight(no)/sqrt(indgen(nmod(no)))*fit(no);
    cmx(1:2) = 0.; // FIXME
    // if (no>=2) cmx *= 0;
    cmn = -cmx;
    grow,coeff,c;
    grow,cmax,cmx;
    grow,cmin,cmn;
  }
  // write,var,exp(-var/2.);
  coeff = coeff(*); cmax = cmax(*); cmin = cmin(*);
  pd.truecoeffs = &coeff;
  // if (coeff_offsets!=[]) coeff -= coeff_offsets; //FIXME, will break with strehl normalisation

  // fill mircube with phases
  if (initphase=="screens") {
    for (no=1;no<=nopt;no++) {
      (*pd.mircube)(,,no) = make_phase_screens(*pd.ipupil,lambda, \
        nm_rmsv(no),ps_slope,remove_tt=1,remove_foc=1);
    }
    // in case of perturbation = screens, zero coefs:
    coeff *=0;
    *pd.truecoeffs *= 0;
    cmin = cmin * 1e3;
    cmax = cmax * 1e3;

    // delta to normalise phase screens not only in the centre-projected pupil,
    // but over all patches (beams) on the various optics:
    dim = dimsof(*pd.ipupil)(2);
    for (no=1;no<=nopt;no++) {
      // dmgsxposcub(coordinates,#optic,#target)
      // excursion of one point over all targets for this optic:
      excursion_x = (*pray_data.dmgsxposcub)(1,no,)(ptp);
      excursion_y = (*pray_data.dmgsyposcub)(1,no,)(ptp);
      npt = 5;
      stride_x = excursion_x/(npt-1);
      stride_y = excursion_y/(npt-1);
      offset_x = indgen(npt)*stride_x; offset_x -= avg(offset_x);
      offset_y = indgen(npt)*stride_y; offset_y -= avg(offset_y);
      phase = (*pd.mircube)(,,no);
      if (debug>5) write,format="phase rms: original= %f, ",phase(*)(rms);
      phase_rms = array(0.,[2,npt,npt]);
      for (i=1;i<=npt;i++) {
        for (j=1;j<=npt;j++) {
          tempphase = phase;
          npup = roll(*pd.ipupil,[offset_x(i),offset_y(j)]);
          w = where(npup);
          // tv,phase*npup;
         	// remove_tt
          tt = (indices(dim)-dim/2.-0.5);
          tip = sum(tempphase(w)*tt(,,1)(w))/sum(tt(,,1)(w)^2);
          tilt = sum(tempphase(w)*tt(,,2)(w))/sum(tt(,,2)(w)^2);
          tempphase = tempphase-tip*tt(,,1)-tilt*tt(,,2);
          // remove_foc
          // focus = dist(dim,xc=dim/2+1,yc=dim/2+1)^2;
          // foc = sum(tempphase(w)*focus(w))/sum(focus(w)^2);
          // tempphase = tempphase-foc*focus;
          // compute rms for this pupil position:
          phase_rms(i,j) = tempphase(w)(rms);
        }
      }
      phase = phase/avg(phase_rms);
      phase = phase*nm_rmsv(no)/lambda*2*pi;
      if (debug>5) write,format="adjusted after FoV average= %f\n",phase(*)(rms);
      (*pd.mircube)(,,no) = phase;
    }
  }

  // compute true phase from truecoeffs for display comparison
  pd.truecube = &array(0.,[3,size,size,nopt]);
  cpt = 0;
  for (no=1;no<=nopt;no++) {
    if (initphase=="screens") {
      (*pd.truecube)(,,no) = (*pd.mircube)(,,no);
    } else {
      zv = cpt+indgen(nmod(no));
      (*pd.truecube)(,,no) = (*pd.def)(,,zv)(,,+) * (*pd.truecoeffs)(zv)(+);
      cpt += nmod(no);
    }
    // (*pd.truecube)(,,no) *= (*pd.maskcube)(,,no);
  }
  return 0;
}

func init_images(&pd,config,&object,&start_strehl)
{
  res = [];
  write,format="%s\n","Computing initial images";
  for (i=1;i<=nfoc;i++) {
    // coeff(1) = deltafoc(i);
    // write,config.roti(i),rotv(,config.roti(i));
    grow,res,&compute_psfs(pd,deltafoc(i),coeff,amp1,amp2, \
      rotv=rotv(,config.roti(i)),nodisp=1,fromscreens=(initphase=="screens"));
  }

  // simple square object ... not used for now, but needed by pray()
  object=array(float,[2,size,size]);
  object(size/2+1,size/2+1) = flux;
  // creating images
  images = array(float,[4,size,size,ntarget,nfoc]);
  strehlv = array(0.,ntarget);
  allstv = [];
  for (n=1;n<=nfoc;n++) {
    for (i=1;i<=ntarget;i++) {
      // images are centred, object is centred res is rol (so eclat(res) is centred)
      // images(,,i,n) = fft_convolve(object,eclat((*res(n))(,,i)));
      images(,,i,n) = eclat((*res(n))(,,i))*flux;
      if (deltafoc(n)==0) {
        strehlv(i) = max(images(,,i,n)/sum(images(,,i,n)))/peak_airy;
        strehlv_at_focus = strehlv;
      }
    }
    if (deltafoc(n)==0) {
      rotvstr = strjoin(swrite(format="%.0f",rotv(,config(n).roti)),",");
      write,format="\033[31mStrehl over FoV (rot=[%s]): avg=%.1f%%\033[0m rms=%.1f%%\n", \
        rotvstr,100*avg(strehlv),100*strehlv(rms);
      grow,allstv,strehlv;
    }
  }
  if (disp) {
    // display for first defoc plane
    // disp_im = build_bigim(roll(images(,,,n)),*pd.xpos,*pd.ypos,variance);
    disp_im = build_bigim(images(,,,1),*pd.xpos,*pd.ypos,variance,noeclat=1)/flux;
    pd.original_big_image = &disp_im;
    window,1; fma; plsys,3;
    pli,disp_im; limits,square=1;
    pth = pltitle_height_vp; pltitle_height_vp = pltitle_height_vp*0+9;
    pltitle,swrite(format="%.2f-focus images - data",deltafoc(1));
    pltitle_height_vp = pth;
    pause,50;
  }

  start_strehl = [avg(allstv),allstv(rms)]; // for all rotations and deltafoc=0
  write,format="\033[31mStrehl over FoV (all rotations): avg=%.1f%%\033[0m rms=%.1f%%\n", \
          100*avg(allstv),100*allstv(rms);

  images = images + random_normal(dimsof(images))*sqrt(variance);

  pd.images = &images;

  return strehlv_at_focus;
}
