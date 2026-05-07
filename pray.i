/*
 * pray.i
 *
 * This file is part of the pray package
 *
 * Copyright (c) 2009, Damien Gratadour
 * Copyright (c) 2022, Francois Rigaut
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
 * Initial release D. Gratadour, Aug 2009.
 */

require,"yao.i";
require,"vops.i";
require,"pray.h";
require,"mavis_pray.h";
// require,"optimpack-mod.i";

sim = sim_struct(); // for call to make_diskharmonic()


func phase2psf(phase,pup,size,&wf,&complexwf)
/* DOCUMENT phase2psf
 *
 * psf=phase2psf(phase,pup,size,wf)
 *
 * This routine returns the PSF estimated from a phase map
 *
 *
 * KEYWORDS :
 * phase      (input)    : Phase map to be used for the PSF computation
 * imwidth    (input)    : Number of pixels on one side of the image
 * wf   (output)         : Complex amplitude in the pupil plane (size = size
 * pup       (input)     : The pupil
 *                       x size)
 *                       the PSF
 *
 * SEE ALSO phase2psfMode
 */
{
  sdim = dimsof(pup)(2);

  if (dimsof(phase)(2) != sdim) error,"phase and pup size mismatch";

  // Complex amplitude in the pupil plane
  wf = array(complex,[2,size,size]);
  (wf.re)(1:sdim,1:sdim) = cos(phase)*pup;
  (wf.im)(1:sdim,1:sdim) = sin(phase)*pup;

  // Complex amplitude in the focal plane
  complexwf = fft(wf,1);
  psf       = abs(complexwf)^2;

  // normalization
  psf      = psf/sum(psf);

  return psf;
}


func pray_j_data(&gradient_psf,object=,ft_object=,image=,ft_image=,psf=, \
                 ft_psf=,variance=)
/* DOCUMENT pray_j_data
 *
 * criterion=pray_j_data(grad_obj,grad_psf,[object=,ft_object=,image=,\
 *                       psf=,ft_psf=,variance=,guard_band=])
 *
 * This routine returns the negative log-likelihood value that measures the
 * fidelity to the data in the phase diversity process. It is the negative
 * log-likelihood of a white stationnary guaussian noise.
 *
 * The convolution is defined as :
 * (h * o)(j,k) = sum_{l,m} [h(l,m) . o(j-l,k-m)]
 * And the criterion is of the form :
 * J(o\i) = sum [1/variance(x,y) x | (h * o)(x,y) - i(x,y)|^2]
 *
 * KEYWORDS :
 * gradient_o   (output) : The gradient of the criterion with respect to the
 *                         object
 * gradient_psf (output) : The gradient of the criterion with respect to the
 *                         psf
 * object       (input)  : The object (a 2D array) for which we want the
 *                         criterion value (not needed if ft_object is passed)
 * ft_object    (input)  : The object Fourier transform (asumed centred) (can
 *                         be passed instead of object to increase the speed
 *                         of the calculation)
 * image        (input)  : The image (a 2D array) for which we want the
 *                         criterion value
 * psf          (input)  : The psf (a 2D array) for which we want the
 *                         criterion value (not needed if ft_psf is passed)
 * ft_psf       (input)  : The psf Fourier transform (asumed centred) (can
 *                         be passed instead of psf to increase the speed
 *                         of the calculation)
 * variance     (input)  : The noise variance. Can be a scalar or a 2D array
 * guard_band   (input)  : An optional guard band to estimate an object wider
 *                         than the original image. This can be useful to
 *                         avoid aliasing in the criterion and gradient
 *                         computation
 *
 * SEE ALSO: yoda, j_prior_gauss
 */
{
  // Test on the object
  if (object!=[]) ft_object=fft(object,1);
  else if (ft_object==[]) write,format="%s\n","object and ft_object missing";


  // Test on the image
  if (image!=[]) ft_image=fft(image,1);
  else if (ft_image==[]) write,format="%s\n","image and ft_image missing";

  // Test on the psf
  if (psf!=[]) ft_psf = fft(psf,1);
  else if (ft_psf==[]) write,format="%s\n","psf and ft_psf missing";

  dim2 = nof(ft_image);
  np = long(sqrt(dim2));         //image must be square

  // Test on the variance
  if ((nof(variance) != dim2) & (nof(variance) != 1)) {
    write,format="%s\n","variance should be scalar or of same size than image.";
  }

  // estimation of (h*o-i)
  ho_minus_i = 1./dim2*double(fft((ft_psf*ft_object),-1))-image;

  // gradient_psf estimation
  gradient_psf=1./dim2*double(fft(conj(ft_object)*fft(ho_minus_i/variance,1),-1));

// In the least square method :
// gradient_psf=(1./variance)*double(fft(conj(ft_object)*ho_minus_i,1));

  return 0.5 * sum((ho_minus_i)^2/variance);
}

func grad_param_psf(grad_psf,&grad_phase,modes_array,mask,ampli_pup,ampli_foc,pupd)
/* DOCUMENT gradParamPsf
 *
 * grad=gradParamPsf(grad_psf,grad_phase,modes_array,mask,ampli_pup)
 *
 * This routine returns the gradient with respect to the mode from the
 * gradient with respect to the PSF
 *
 * KEYWORDS :
 * grad_psf     (input)  : The gradient with respect to the PSF
 * grad_phase   (output) : The gradient with respect to the Phase
 *                         (pupSize x pupSize)
 * mask         (input)  : The pupil
 * modes_array  (input)  : Transformation matrix from the mode coefficients
 *                         space to the phase space (nbModes x pupsize^2)
 * ampli_pup    (input)  : The complex amplitude in the pupil plane
 *                         (pupSize x pupSize)
 *
 * SEE ALSO:
 */
{
  // modes_array contains the Zi in column ...
  // in th tomographic case, modes Array contains the modes, as viewed in the pupil for
  // the specific direction ... need to use getModesInPup before
  //pupd2 = pud^2;
  size2 = nof(mask);
  size  = sqrt(size2);

  dummy = where(mask != 0.);
  count = nof(dummy);

  //grad_phase = (2./count)*(conj(ampli_pup)*fft(grad_psf*fft(ampli_pup,-1),1)).im;

  grad_phase = (-2./size2/count)*(ampli_pup*fft(grad_psf*conj(ampli_foc),1)).im;

  // grad_phase is nil outside (1:pupd,1:pupd) so we need to recentre it for the product with the modes
  // note that modes are also nil outside the same area by definition
  grad_phase = roll(grad_phase,[size/2-pupd/2-2,size/2-pupd/2-2]);

  gradMode = modes_array(*,)(where(mask),)(+,)*grad_phase(*)(where(mask))(+);

  return gradMode;
}


func compute_psfs(&pd,xfoc,coeff,&ampli_pup,&ampli_foc,rotv=,nodisp=,fromscreens=)
/*
 */
{
  nopt     = nof(*pd.alt);
  ntarget  = dimsof(*pd.dmgsxposcub)(4);

  *pd.mircube *= 0;

  if (rotv==[]) rotv=array(0.,nopt);

  tic,4;
  cpt = 0;
  for (i=1;i<=nopt;i++) {
    zv = cpt+indgen((*pd.nmod)(i));
    // ok so it's here we have to act if we want to use predefined phases instead of coeffs:
    if (fromscreens) (*pd.mircube)(,,i) = (*pd.truecube)(,,i); \
    else (*pd.mircube)(,,i) = (*pd.def)(,,zv)(,,+)*coeff(zv)(+); // from coeffs
    // else mircube has been prefill with correct phases at each optics/altitudes
    if (rotv(i)!=0) (*pd.mircube)(,,i) = rotate2((*pd.mircube)(,,i),rotv(i),xc=pd.size/2+0.5,yc=pd.size/2+0.5);
    cpt += nmod(i);
  }

  if ((allof(rotv==0))&(nodisp!=1)) {
    window,2; fma;
    pltitle_height_vp = 5;
    isys = 1;
    for (no=1;no<=min([13,nopt]);no++) {
      // mask = smooth((*pd.truecube)(,,i)!=0)>0.9;
      mask = (*pd.maskcube)(,,no);
      plsys,isys++;
      pli,(*pd.truecube)(,,no) * mask;
      pltitle_vp,swrite(format="True phase(%d)",no),pos=-1;
      plsys,isys++;
      pli,(*pd.mircube)(,,no) * mask;
      pltitle_vp,swrite(format="Phase(%d)",no),pos=-1;
      plsys,isys++;
      pli,((*pd.truecube)(,,no)-(*pd.mircube)(,,no)) * mask;
      pltitle_vp,swrite(format="Difference(%d)",no),pos=-1;
    }
    pause,20;
  }

  // add focus to generate extra focal image
  (*pd.mircube)(,,1) += xfoc * (*pd.focus);

  bphase = array(float,[3,pd.size,pd.size,ntarget]);
  psnx = dimsof(*pd.mircube)(2);
  psny = dimsof(*pd.mircube)(3);
  skip = array(0n,nopt);

  for (k=1;k<=ntarget;k++) {
    sphase = array(float,pd._n,pd._n);
    xshifts = (*pd.dmgsxposcub)(,,k)+(pd.centre-1)(-,);
    yshifts = (*pd.dmgsyposcub)(,,k)+(pd.centre-1)(-,);

    ishifts = int(xshifts); xshifts = xshifts - ishifts;
    jshifts = int(yshifts); yshifts = yshifts - jshifts;

    err = get_2d_phase_from_cube(pd.mircube,psnx,psny,nopt,&skip,&sphase,pd._n,pd._n,\
                      &ishifts,&float(xshifts),&jshifts,&float(yshifts));

    if (err != 0) { error,"Error in getPhase2dFromDms"; }
    bphase(pd._n1:pd._n2,pd._n1:pd._n2,k) = float(sphase);
  }

  psfs = array(float,[3,pd.size,pd.size,ntarget]);
  ampli_pup = ampli_foc = array(complex,[3,size,size,ntarget]);

  for (i=1;i<=ntarget;i++) {
    psfs(,,i) = phase2psf(bphase(pd._n1:pd._n2,pd._n1:pd._n2,i),\
      (*pd.ipupil)(pd._n1:pd._n2,pd._n1:pd._n2),pd.size,amp1,amp2);
    if (centre_images) {
      ima = eclat(psfs(,,i));
      psfs(,,i) = eclat(centre_image(ima,pd));
    }
    ampli_pup(,,i) = amp1;
    ampli_foc(,,i) = amp2;
  }
  pd.ampli_pup = &ampli_pup;
  pd.ampli_foc = &ampli_foc;

  extern compute_psf_time;
  if (compute_psf_time==[]) compute_psf_time = 0.;
  compute_psf_time += tac(4);

  return psfs;
}


func get_def_in_pupil_from_dir(pd,ndir,rotv=)
/*
  return the array of modes as viewed from the pupil in a specified direction
 */
{
  tic,4;

  nopt        = nof(*pd.alt);

  psnx        = dimsof(*pd.mircube)(2);
  psny        = dimsof(*pd.mircube)(3);

  if (rotv==[]) rotv = array(0.,nopt);

  // geometry init : get the proper points coordinates : dmgsxposcub
  xshifts = (*pd.dmgsxposcub)(,,ndir)+(pd.centre-1)(-,);
  yshifts = (*pd.dmgsyposcub)(,,ndir)+(pd.centre-1)(-,);

  ishifts = int(xshifts); xshifts = xshifts - ishifts;
  jshifts = int(yshifts); yshifts = yshifts - jshifts;

  mydef = (*pd.def)*0.;
  *pd.mircube *= 0.0f;
  cpt = 0;
  skip = array(0n,nopt);
  for (no=1;no<=nopt;no++) {
    for (j=1;j<=nmod(no);j++) {
      *pd.mircube *= 0.0f;
      bphase = array(float,[2,pd.size,pd.size]);
      sphase = array(float,pd._n,pd._n);
      cpt++;
      (*pd.mircube)(,,no) = float((*pd.def)(,,cpt));
      if (rotv(no)!=0) (*pd.mircube)(,,no) = rotate2((*pd.mircube)(,,no),rotv(no),xc=pd.size/2+0.5,yc=pd.size/2+0.5);
      // this function was taking forever as we were sending the whole mircube,
      // in this version I am just sending the relevant plane as there is
      // just one set anyway - well, it's not taking much less (25% less).
      cubeslice = (*pd.mircube)(,,no)(,,-);
      err = get_2d_phase_from_cube(&cubeslice,psnx,psny,1n,&[skip(no)],&sphase,pd._n,pd._n,\
        &(ishifts(,no)),&float(xshifts(,no)),&(jshifts(,no)),&float(yshifts(,no)));
      if (err != 0) error,"Error in getPhase2dFromDms";
      bphase(pd._n1:pd._n2,pd._n1:pd._n2) = float(sphase);
      mydef(,,cpt) = bphase;
    }
  }
  extern get_def_in_pupil_from_dir_time;
  if (get_def_in_pupil_from_dir_time==[]) get_def_in_pupil_from_dir_time = 0.;
  get_def_in_pupil_from_dir_time += tac(4);

  return mydef;
}

func pray_error(coeffs,&gradient,extra=)
/* DOCUMENT pray_error
 *
 * criterion=pray_error(coefs,gradient,extra)
 *
 * This routine returns the error (global criterion) function to be
 * minimized in pray
 *
 * KEYWORDS :
 * coefs        (input)  : The parameters to be mimimized
 * gradient    (output)  : The gradient of this error function with repect
 *                         to the parameters
 * extra        (input)  : A structure containing all necessary parameters
 *                         (see the content of pday_struct() for more details).
 *
 * SEE ALSO:
 */
{
  local deltafoc;

  ncoeffs   = nof(coeffs);

  ntarget = dimsof(*extra.images)(4);
  // int crit_array
  crit_array = array(0.,nof(*extra.deltafoc)*ntarget);
  gradient_interm = array(float,ncoeffs,nof(*extra.deltafoc)*ntarget);

  //----------------------------------------------------------------
  //Extra-Focal images (we can introduce as many images as we want)
  ima = (*extra.images)(,,*);
  for (n=1;n<=nof(*extra.deltafoc);n++) {
    tmp = coeffs;
    // tmp(1) += deltafoc(n);
    // below rotv defined where? FIXME
    psfs = compute_psfs(extra,(*extra.deltafoc)(n),tmp,*pd.ampli_pup,*pd.ampli_foc,rotv=rotv(,config.roti(n)));
    for (i=1;i<=ntarget;i++) {
      ftPsf = fft(psfs(,,i),1);
      // Estimation of the criterion associated with image #i

      crit_array((n-1)*ntarget+i) = pray_j_data(gradientPsf,ft_object=*extra._ftobject,\
        image=ima(,,(n-1)*ntarget+i),ft_psf=ftPsf,variance=*extra.variance);

      ri = config(n).roti;
      tmp = grad_param_psf(gradientPsf,grad_phaseOut,(*extra._def_pup)(,,,i,ri),*extra.ipupil,\
        (*extra.ampli_pup)(,,i),(*extra.ampli_foc)(,,i),extra.pupd);

      gradient_interm(,(n-1)*ntarget+i) =  tmp;
    }
  }
  //----------------------------------------------------------------


  gradient = gradient_interm(,sum);

  //gradient*=norm;

  return crit_array(sum);
}


func pray(images,pd,deltafoc,variance,object,disp=,verbose=,threshold=,nbiter=,\
  shift_foc=)
/* DOCUMENT pray
 *
 *
 * This routine returns a set of mode coefficients that describes the
 * aberrations from a couple images (in focus and out of focus). It uses a
 * least square approach, penalized or not.
 *
 * the criterion contains only a term of fidelity to the data.
 *
 * KEYWORDS :
 * images     (input) : The couple of images (a 2D arrays) from which the
 *                      aberrant phase can be estimated
 * object     (input) : The initial object (a 2D array) which is imaged
 * guess      (input) : A guess on the restored mode coefficients (a vector
 *                      of nbModes components)
 * cobs       (input) : The size of central obscuration on the pupil
 * variance   (input) : The noise variance in the image
 * nbiter     (input) : Number max of iterations for the criterion minimization
 * disp       (input) : A flag to display (or not) the step-by-step results
 * verbose    (input) : A flag to print (or not) the details of the
 *                      deconvolution process
 * SEE ALSO:
 */
{
  extern firstdefoc;
  firstdefoc = deltafoc(1);
  local variance;

  // Init verbose and disp flags
  if (!is_set(verbose)) verbose=0;
  if (!is_set(disp)) disp=0;
  if (!shift_foc) shift_foc=0.;

  // sizes init and check
  dims = dimsof(images);
  if (dims(1) != 4) {
    error,"images format must follow : [size,size,ntarget,2]";
  }
  size = dims(2);
  ntarget = dims(4);
  if (nof(*pd.xpos) != ntarget) error,"incompatible images and xpos size";
  if (nof(*pd.ypos) != ntarget) error,"incompatible images and ypos size";


  // pray init, already done within mavis_pray()
  // status = init_defs(pd);

  // ... testing the size and content of variance
  if (variance==[]) {
    variance = sigma^2;
    sz_variance = 1;
  } else {
    sz_variance = nof(variance)
      if ((sz_variance != 1) & (sz_variance != size^2)) {
        write,format="%s\n","variance must be a scalar or of same size than image.";
        error;
      } else {
        if (min(variance) <= 0.) {
          error,swrite(format="%s\n","variance must be strictly positive");
        }
      }
    if (sz_variance == size^2) {
      dynamic_var = max(variance)/min(variance);
      if (dynamic_var > 10000.) {
        write,format="%s\n","WARNING : variance has a huge dynamics.";
        write,format="%s\n","Be prepared for minimization problems !";
      }
    }
  }

  // ... testing guess and create one if nill
  guess = array(0.0,nmod(sum));//BN (Benoit Neichel???)
  write,format="%s\n","\033[31mUsing coeff = 0 as first guess\033[0m";
  // guess = truecoeff+random_n(dimsof(truecoeff))*truecoeff(rms);
  // write,format="%s\n","\033[31mUsing true coeff+noise as first guess\033[0m";
  norm = guess;
  // param = guess*0. + 0.1;
  param = guess;

  _ftobject=fft(object,1);

  if (nof(threshold) == 0) {
    if (typeof(images) == "float") {
      conv_threshold = 1e-12;
    } else conv_threshold = 1e-16;
  } else conv_threshold = threshold;
  // conv_threshold=1e-21;
  // cmax *=10;
  // cmin *=10;

  if (verbose) write,format="%s\n","Using the VMLM-B method for minimization.";

  pd._ftobject   = &_ftobject;
  pd.variance   = &variance;
  pd.images     = &images;
  pd.norm       = &norm;
  // pd.regul      = regul;
  pd.deltafoc   = &deltafoc;
  // pd.shift_foc   = shift_foc;

  // Some initializations for the minimization process
  iter = old_fout = old_gout = cont_out = 0;
  cont_flag = 1;
  coeff_arr = array(float,(*pd.nmod)(sum));

  if (verbose) {
    write, format="%s  %s\n%s  %s\n",
      " ITER    EVAL     CPU [s]            FUNC             max(|G|)",
      " STEPLEN",
      "------  ------  ----------  -----------------------  ---------",
      "---------";
  }

  // tic;

  minim=optm_vmlmb(pray_error,param,fout,gout, \
                 lower=cmin,upper=cmax,ftol=conv_threshold,  \
                 gtol=conv_threshold,maxiter=nbiter,verb=10,observer=myobserver,extra=pd);

  param=minim;

  // tac();
  if (verbose) {
    write,format="Value of the criterion after %d iterations :",iter;
    criterion_value =  pray_error(minim,grad_fin);
    write,criterion_value;
    if (verbose>1) {
      write,format="Values of the Mode Coeffs after %d iterations :",iter;
      write,minim;
    }

  }

  return minim; //*norm;
}

func myobserver(iters,evals,rejects,t,x,f,g,gpnorm,alpha,fg,extra=)
{
  extern currentiter,objfuncvec,gradvec;
  if (iters==1) {
  	objfuncvec = [f];
    // gradvec = [g];
  }
  else {
    grow,objfuncvec,f;
    // grow,gradvec,g;
  }
  currentiter = iters;
  if (!disp) return;
  if (extra==[]) return;
  if ((iters>10) & ((iters%10)!=0)) return;
  coeff = x;
  if (window_exists(3)) window,3;
  else window,3,wait=1,dpi=long(dpi_target_small);
  fma; limits,square=0;
  plh,coeff;
  xytitles,"Coefficient number (over all optics)","Coefficient value",[-0.015,0.];
  pltitle,"Coefficients";
  if (initphase=="coefs") {
	  plh,(*extra.truecoeffs),color="red";
	  plh,(*extra.truecoeffs)-coeff,color="blue";
  }
  // coeff(1) += firstdefoc;
  psfs = compute_psfs(extra,firstdefoc,coeff);
  window,1,wait=1;
  fma;
  pth = pltitle_height_vp; pltitle_height_vp = pltitle_height_vp*0+9;
  plsys,3; pli,*extra.original_big_image; limits,square=1;
  pltitle,swrite(format="%.2f-focus images - data",firstdefoc);
  disp_im = build_bigim(psfs,*extra.xpos,*extra.ypos);
  plsys,2; pli,disp_im; limits,square=1;
  pltitle,swrite(format="%.2f-focus images - model - iter %d",firstdefoc,iters);
  plsys,1; pli,*extra.original_big_image - disp_im; limits,square=1;
  pltitle,swrite(format="%.2f-focus images - difference - iter %d",firstdefoc,iters);
  pltitle_height_vp = pth;
  if (window_exists(4)) window,4;
  else window,4,wait=1,dpi=long(dpi_target_small);
  if (iters>1) {
	  fma;
		plh,objfuncvec,indgen(iters);
		// plh,gradvec,indgen(iters);
	  logxy,0,1;
		pltitle,"Objective function";
		xytitles,"Iteration","Objective function",[-0.020,0.];
  }
  // if (iters>90) error;
}
