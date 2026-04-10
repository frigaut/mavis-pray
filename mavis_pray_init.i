func init_windows(win)
/* DOCUMENT
 * init windows
 */
{
  for (i=1;i<=nof(win);i++) {
    if (window_exists(win(i))) {
      window,win(i);
      if (win(i)!=5) { limits; fma; }
    }
  }
  nw = clip(nopt,4,8); nw = 8;
  if (!window_exists(2)) plsplit,2,nw,win=2,style="nobox.gs",square=1,dpi=long(dpi_target*1.5),\
    margin=-0.02,vp=[0.206+0.0115*(nw-3),0.656+0.0115*(nw-3),0.44,0.85];
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


func init_defs(&pd,tiptilt=)
/* DOCUMENT init_defs
 * initialise defs and other arrays
*/
{
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
  def	= array(float,[3,pd.size,pd.size,(*pd.nzer)(sum)]);
  nz12 = (*pd.nzer)(cum);
  psize = pd.teldiam/pd.pupd;
  cpt = 0;

  // create modes per optic
  for (k=1;k<=nopt;k++) {
    patchDiam = long(pd.pupd+2.*max(abs(*pd.xpos,*pd.ypos))*4.848e-6*(abs(alt(k)))/psize);
    prepzernike,pd.size,patchDiam,pd.centre,pd.centre;

    // extern focus;
    // FIXME: there is only one array and several nopt? yes, because it is
    // way oversized and the amplitude doesn't really matter?
    // pd.focus = &zernike_ext(4); // true focus for extra focal image calculations
    // FIXED
    pd.focus = &((sqrt(3.)*dist(pd.size,xc=pd.size/2.+0.5,yc=pd.size/2.+0.5)/(pd.pupd/2.))^2.);

    if (usemodes=="zer") {
      // we want zernike_ext but not too much (1 pixel), so we create our own mask
      mask = smooth(zernike(1))>0.5;
      for (i=4;i<=nzer(k)+3;i++) {
        cpt ++;
        def(,,cpt) = zernike_ext(i)*mask;
      }
      if (tiptilt != []) {
       write,format="%s\n","Estimation of global TipTilt";
        selz = _(4,2,3,4+indgen(nzer(k)+10));
        // get rid of spherical
        //selz(10) = max(selz)+1;
        for (i=4;i<=nzer(k)+3;i++) {
          cpt ++;
          //def(,,cpt) = zernike_ext(i)*mask;
          def(,,cpt) = zernike_ext(selz(i-3))*mask;
        }
      }
    } else if (usemodes=="kl") { // use KL
      if ((k==1)&verbose) write,format="%s\n","Using KL instead of Zernike";
      require,"yaokl.i";
      pup1 = [];
      // def *= 0;
      thispatch = ((patchDiam+2)/2)*2;
      // computes KLs, remove Tip and Tilt like
      kl = make_kl((*pd.nzer)(k)+2,thispatch,v,obas,pup1,oc=0.0,nr=128,verbose=0)(,,3:);
      // normalise so that future coef optimisation will be about the same amplitude
      kl = kl*(1./indgen((*pd.nzer)(k)+2)^0.8)(-,-,3:);
      //    kl = order_kls(kl,patchDiam,upto=20); // why is that not necessary?
      // def is supposed to have zernikes in them, with def(,,1) = true focus
      def(,,2+nz12(k):nz12(k+1)) *= 0.;
      // stick in the KL, preserving zernike focus in position 1
      def(size/2-patchDiam/2:size/2+patchDiam/2+1, \
        size/2-patchDiam/2:size/2+patchDiam/2+1,2+nz12(k):nz12(k+1)) = kl(,,2:);
      if (k==1) def(,,1+nz12(k)) = *pd.focus;
      if (k==nopt) def(,,1+nz12(k)) = *pd.focus;
    } else if (usemodes=="dh") { // use DH
      if ((k==1)&verbose) write,format="%s\n","Using DH instead of Zernike";
      require,"yaodh.i";
      sim.verbose = 0;
      pup1 = [];
      // def *= 0;
      thispatch = ((patchDiam+2)/2)*2;
      // thispatch = patchDiam; //write,thispatch;
      // computes KLs, remove Piston, Tip and Tilt
      dh = make_diskharmonic(pd.size,thispatch,(*pd.nzer)(k)+6,cobs=0.,xc=size/2+0.5,yc=size/2+0.5)(,,4:(*pd.nzer)(k)+3);
      // focus is in position 2 now. astig in 1, put astig in 2
      dh(,,2) = dh(,,1); // focus will be added below in pos 1
      // normalise so that future coef optimisation will be about the same amplitude
      dh = dh*(1./indgen((*pd.nzer)(k)+6)^0.8)(-,-,4:(*pd.nzer)(k)+3);
      // def is supposed to have zernikes in them, with def(,,1) = true focus
      def(,,1+nz12(k):nz12(k+1)) *= 0.;
      // stick in the DH, preserving zernike focus in position 1
      def(,,2+nz12(k):nz12(k+1)) = dh(,,2:);
      // if (k==2) error;
      // def = def(,,_(3,1,2,indgen(4:cpt)));}
      lpup = (def(,,1+nz12(k)) != 0);
      if (k==1) def(,,1+nz12(k)) = *pd.focus;
      if (k==nopt) def(,,1+nz12(k)) = *pd.focus;
      // if (k==1) def(,,1+nz12(k):nz12(k+1)) = def(::-1,::-1,1+nz12(k):nz12(k+1));
      // window,1; tv,def(,,nz12(k)+7); pltitle,swrite(format="%d",k); pause,1000;
    } else error,"usemodes undefined";

  }
  pd.def = &def;

  // Init dmgsXYposcub : will be needed by _get2dPhase
  xref = indgen(pd._n)-(pd._n+1)/2.;
  yref = indgen(pd._n)-(pd._n+1)/2.;

  dmgsxposcub = dmgsyposcub = array(float,[3,pd._n,nopt,ntarget]);

  for (n=1;n<=ntarget;n++) {
    // loop on pseudo-DMs
    for (no=1;no<=nopt;no++) {
      // offsets of the center of beam on DM NO
      xoff = (*pd.xpos)(n)*4.848e-6*(*pd.alt)(no)/psize;
      yoff = (*pd.ypos)(n)*4.848e-6*(*pd.alt)(no)/psize;
      dmgsxposcub(,no,n) = xref + xoff;
      dmgsyposcub(,no,n) = yref + yoff;
    }
  }

  pd.dmgsxposcub = &dmgsxposcub;
  pd.dmgsyposcub = &dmgsyposcub;

  szdef = dimsof(def);
  def_pup = array(float,[5,szdef(2),szdef(3),szdef(4),ntarget,nrot]);
  for (n=1;n<=nrot;n++) {
    for (i=1;i<=ntarget;i++) def_pup(,,,i,n) = get_def_in_pupil_from_dir(pd,i,rotv=rotv(,n));
  }
  pd.def_pup = &def_pup;
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
  nz12 = nzer(cum);
  // var = 0.;
  for (no=1;no<=nopt;no++) {
    c = weight(no)*(random(nzer(no))-0.5)/sqrt(indgen(nzer(no)));
    c(1) = 0.; // FIXME
    cmx = 10*weight(no)/sqrt(indgen(nzer(no)))*fit(no);
    cmx(1) = 0.; // FIXME
    // if (no>=2) cmx *= 0;
    cmn = -cmx;
    grow,coeff,c;
    grow,cmax,cmx;
    grow,cmin,cmn;
  }
  // write,var,exp(-var/2.);
  coeff = coeff(*); cmax = cmax(*); cmin = cmin(*);
  pd.truecoeffs = &coeff;
  if (coeff_offsets!=[]) coeff -= coeff_offsets;

  // fill mircube with phases
  if (initphase=="screens") {
    for (no=1;no<=nopt;no++) {
      (*pd.mircube)(,,no) = make_phase_screens(*pd.ipupil,lambda, \
        nm_rmsv(i),ps_slope,rseed=rseed+i*0.01,remove_tt=1,remove_foc=1);
    }
    // in case of perturbation = screens, zero coefs:
    coeff *=0;
    *pd.truecoeffs *= 0;
    cmin = cmin * 1e3;
    cmax = cmax * 1e3;
  }

  // compute true phase from truecoeffs for display comparison
  pd.truecube = &array(0.,[3,size,size,nopt]);
  cpt = 0;
  for (no=1;no<=nopt;no++) {
    if (initphase=="screens") {
      (*pd.truecube)(,,no) = (*pd.mircube)(,,no);
    } else {
      zv = cpt+indgen(nzer(no));
      (*pd.truecube)(,,no) = (*pd.def)(,,zv)(,,+) * (*pd.truecoeffs)(zv)(+);
      cpt += nzer(no);
    }
    (*pd.truecube)(,,no) *= (*pd.maskcube)(,,no);
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
      images(,,i,n) = roll((*res(n))(,,i))*flux;
      if (deltafoc(n)==0) \
        strehlv(i) = max(images(,,i,n)/sum(images(,,i,n)))/peak_airy;
    }
    if (deltafoc(n)==0) {
      rotvstr = strjoin(swrite(format="%.0f",rotv(,config(n).roti)),",");
      write,format="\033[31mStrehl over FoV (rot=[%s]): avg=%.1f%%\033[0m rms=%.1f%%\n", \
        rotvstr,100*avg(strehlv),100*strehlv(rms);
      grow,allstv,strehlv;
      if (n==1) start_strehl = [avg(strehlv),strehlv(rms)];
      // FIXME, right now the start_strehl is for 180 while end_strehl is for 0
    }

    if (disp&&(n==1)) {
      // display for first defoc plane
      disp_im = build_bigim(roll(images(,,,n)),*pd.xpos,*pd.ypos,variance);
      extern bigim2; bigim2 = disp_im;
      if (window_exists(n)) window,n;
      else window,n,wait=1,dpi=dpi_target;
      fma; pli,disp_im; limits,square=1;
      pltitle,swrite(format="%.2f-focus images - data",deltafoc(n));
      pause,50;
    }
  }
  write,format="\033[31mStrehl over FoV (all rotations): avg=%.1f%%\033[0m rms=%.1f%%\n", \
          100*avg(allstv),100*allstv(rms);
  images = images + random_normal(dimsof(images))*sqrt(variance);

  pd.images = &images;

  return 0;
}
