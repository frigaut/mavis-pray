/* projection.i
Computes projection matrices from modes at a given altitude onto modes
at another altitude.
Based on the ratio diameters, which can be readily computed from the FoV
and delta altitude for our NCPA problem.
Modes are either zernike, DH or KL.
*/
require,"yao.i";
sim = sim_struct(); // for call to make_diskharmonic()
sim.verbose = 0
if (debug==[]) debug=0;

func generate_modes(modes, nmodes, dim, patchd, &pupil)
/* DOCUMENT generate_modes(modes, nmodes, dim, patchd, &pupil)
 * Generate either zernike, KL or DH modes.
 * Each will be returned as a cube of dimensions [dim,dim,nmodes]
 * pupil is also returned.
 * We manipulate the modes to exclude piston.
 * The first modes are ordered: tip,tilt,focus, astigs
 */
{
  mod = array(0.f,[3,dim,dim,nmodes]);
  pupil = array(0.0f,[2,dim,dim]);

  if (modes=="zer") {
    prepzernike,dim,patchd,dim/2.+0.5,dim/2.+0.5;
    for (i=1;i<=nmodes;i++) mod(,,i) = zernike(i+1);
    pupil = zernike(1);
  } else if (modes=="zerext") {
    prepzernike,dim,patchd,dim/2.+0.5,dim/2.+0.5;
    for (i=1;i<=nmodes;i++) mod(,,i) = zernike_ext(i+1);
    pupil = zernike_ext(1);
  } else if (modes=="kl") {
    require,"yaokl.i";
    v = obas = pup1 = []; // goes around a make_kl/yorick bug
    kl = make_kl(nmodes,patchd,v,obas,pup1,oc=0.0,nr=128,verbose=0);
    // kl = order_kls(kl,patchd,upto=20);
    off = (dim-patchd)/2; i1 = 1+off; i2 = i1+patchd-1;
    mod(i1:i2,i1:i2,) = kl;
    pupil(i1:i2,i1:i2) = pup1;
  } else if (modes=="dh") {
    require,"yaodh.i";
    mod = make_diskharmonic(dim,patchd,nmodes+1,cobs=0.,xc=dim/2+0.5,yc=dim/2+0.5)(,,2:);
    if (nmodes>=3) mod(,,[3,4]) = mod(,,[4,3]); // switch astig45 and focus
    w = where(mod(,,min)==0);
    pupil += 1; pupil(w)=0;
    mod = mod*pupil(,,-); // fixes non zero pixels outside pupil for azimuthal order=0 modes
  }
  return mod;
}

func remove_mode(pha,pupil,mode)
{
  c = sum(pha*pupil*mode)/sum(mode*mode*pupil);
  pha = (pha-c*mode)*pupil;
  return pha;
}



// func proj_modes_from_to(modes,ratio,nmod1,nmod2)
// /* DOCUMENT proj_modes_from_to(modes,ratio,nmod1,nmod2)
//  * Compute the projection matrix from modes at plan2 to modes at plan1
//  * Plan2 is the largest patch (diameter)
//  * modes: "zer","kl" or "dh"
//  * ratio is the >1 ratio of patches = patchd2/patchd1
//  * nmod1 and nmod2 are the number of modes in plans 1 and 2.
//  * Returns p2on1 (projection of 2 on 1) array [2,nmod1,nmod2]
//  * to get coefficients on plan1 knowing coefficient on plan2 coef2:
//  * coef1 = p2on1(,+)*coef2(+)
//  * SEE ALSO: proj_modes_from_to_for_study() in projection_study.i for
//  * a full blow detailled function and analysis. This has been checked, it works.
//  */
// {
//   if (ratio<1) {
//     // error,"ratio has to be >=1";
//     if (debug>10) write,format="%s\n","Ratio < 1, inverting"
//     ratio = 1./ratio;
//   }
//   if (optim_fov_factor!=[]) ratio = 1+(ratio-1)*optim_fov_factor; // doesn't appear to work
//   sim = sim_struct(); sim.verbose = 0;
//   dim = 64;
//   patchd2 = dim;
//   patchd1 = lround(dim/ratio);
//   p2on1 = array(0.,[2,nmod1,nmod2]);
//   // can't make modes with fractional patchd. will have to up dim and patchd to
//   // do this calculation to reduce effect of rounding errors
//   // kernel diameter is patch2-patch1
//   ker = dist(dim)<=((patchd2-patchd1)/2.);
//   ker = float(ker)/sum(ker);
//   mod1 = generate_modes(modes, nmod1, dim, patchd1, pupil);
//   pupil1 = pupil;
//   mod2 = generate_modes(modes, nmod2, dim, patchd2, pupil);
//   pupil2 = pupil;
//   // modes in plan1 (smallest patch) to determine projection of phase to modes1:
//   wpup1 = where(pupil1);
//   mod1lin = mod1(*,)(wpup1,);
//   mod1lin -= avg(mod1lin);
//   mtm = mod1lin(+,)*mod1lin(+,);
//   proj = LUsolve(mtm)(+,)*mod1lin(,+);
//   for (nm2=1;nm2<=nmod2;nm2++) {
//     mode2 = mod2(,,nm2);
//     cm = fft_convolve(mod2(,,nm2),ker);
//     cm *= pupil1;
//     cm = cm(*)(wpup1);
//     cm -= avg(cm);
//     p2on1(,nm2) = proj(,+)*cm(+);
//   }
//   return p2on1;
// }


// func compute_dms_projector(modes,ratiov,nmod_opt,nmod_dm,cond=)
// /* DOCUMENT compute_dms_projector(modes,ratiov,nmodv)
//  * Computes the optimal projection matrix from one plan (plano=plan object)
//  * to several DMs.
//  * modes: "zer","kl" or "dh"
//  * ratiov: vector of ratio of DMs to plano.
//  * NOT NAYMORE: All elements have to be >= 1 and in creasing order.
//  * example: [1.5, 1.8]. See proj_modes_from_to. In this example, the ratio
//  * between plano and the DMs (2 DMs) are 1.5 and 1.8. The ratio between the proj_to_dms
//  * themselves will be cmoputed from this. numberof(ratiov) = number of DMs
//  * nmod_opt: number of modes for plano optics
//  * nmod_dm: number of modes for DMs (has to be the same)
//  */
// {
//   if (cond==[]) cond=100.;
//   pltitle_height = 14;
//   // For a 3-DM system at altitudes h_dm = [h1, h2, h3]
//   // and one optic at h_opt with nmod_opt modes
//   // each DM has nmod_dm modes (could differ; let's say all the same)

//   // === Build the joint Gramian G_tilde ===
//   // Diagonal blocks: identity if modes are orthonormal on each DM's meta-pupil
//   // Off-diagonal blocks: inter-DM cross-projection
//   ndm = numberof(ratiov);
//   G_tilde = array(0., [2,ndm*nmod_dm,ndm*nmod_dm]);
//   for (i=1;i<=ndm;i++) {
//     // diagonal block (assuming orthonormal modes on DM i's meta-pupil)
//     G_tilde(1+(i-1)*nmod_dm:i*nmod_dm, 1+(i-1)*nmod_dm:i*nmod_dm) = unit(nmod_dm);
//     for (j=i+1;j<=ndm;j++) {
//       // ratio of meta-pupil diameters between DM i and DM j
//       ratio_ij = ratiov(j)/ratiov(i);
//       P_ij = proj_modes_from_to(modes, ratio_ij, nmod_dm, nmod_dm);
//       G_tilde(1+(j-1)*nmod_dm:j*nmod_dm, 1+(i-1)*nmod_dm:i*nmod_dm) = transpose(P_ij);
//       G_tilde(1+(i-1)*nmod_dm:i*nmod_dm, 1+(j-1)*nmod_dm:j*nmod_dm) = P_ij;
//       // G_tilde(1+(j-1)*nmod_dm:j*nmod_dm, 1+(i-1)*nmod_dm:i*nmod_dm) = P_ij;
//       // G_tilde(1+(i-1)*nmod_dm:i*nmod_dm, 1+(j-1)*nmod_dm:j*nmod_dm) = transpose(P_ij);
//     }
//   }
//   if (debug>100) { tv,G_tilde; pltitle,escapechar("G_tilde"); hitReturn; }
//   // Build the RHS M_tilde
//   // For one optic, stacked over DMs
//   M_tilde = array(0.,[2,ndm*nmod_dm,nmod_opt]);
//   for (i=1;i<=ndm;i++) {
//     M_tilde(1+(i-1)*nmod_dm:i*nmod_dm, ) = proj_modes_from_to(modes, ratiov(i), nmod_dm, nmod_opt);
//   }
//   if (debug>100) { tv,M_tilde; pltitle,escapechar("M_tilde"); hitReturn; }
//   // Solve
//   ev = SVdec(G_tilde,u,vt);
//   evi = 1./ev;
//   w = where(ev<max(ev/cond));
//   if (nof(w)) evi(w) = 0.;
//   G_tilde_inv = (transpose(vt)(,+)*(diag(evi))(+,))(,+)*transpose(u)(+,);
//   // with regul
//   if (debug) write,format="G_TILDE_INV WITH REGULARIZATION=%.f RIGHT NOW!\n",1./cond;
//   G_tilde_inv = LUsolve(G_tilde+unit(ndm*nmod_dm)/cond);
//   P_joint = G_tilde_inv(,+) * M_tilde(+,);
//   if (debug>100) { tv,P_joint; pltitle,escapechar("P_joint"); hitReturn; }
//   // Apply: c_dm_stacked = P_joint(,+) * c_opt(+)
//   //
//   // CHECKS
//   if (debug>90) {
//     if (ndm!=2) {
//       write,"Check only coded for 2 DMs";
//       return P_joint
//     }
//     dim = 64;
//     mod_opt = generate_modes(modes, nmod_opt, dim, dim, pupil);
//     pupil_opt = pupil;
//     mod_dm1 = generate_modes(modes, nmod_dm, dim, lround(dim/ratiov(1)), pupil);
//     pupil_dm1 = pupil;
//     mod_dm2 = generate_modes(modes, nmod_dm, dim, lround(dim/ratiov(2)), pupil);
//     pupil_dm2 = pupil;

//     if (optmode==[]) optmode=13; else optmode=long(optmode);
//     c_opt = array(0.,nmod_opt); c_opt(optmode)=1.;
//     dmc = P_joint(,+) * c_opt(+);
//     pha_opt = mod_opt(,,+)*c_opt(+);
//     pha_dm1 = mod_dm1(,,+)*dmc(1:nmod_dm)(+);
//     pha_dm2 = mod_dm2(,,+)*dmc(1+nmod_dm:2*nmod_dm)(+);
//     tv,transpose(_(pha_opt,pha_dm1,pha_dm2));
//     limits,square=1;
//   }
//   return P_joint;
// }

func get_def(pd,noptic)
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  return (*pd.def)(,,nz1(noptic):nz2(noptic));
}

func get_coeff(pd,noptic)
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  return (*pd.coeffs)(nz1(noptic):nz2(noptic));
}

func add_to_coeff(c,pd,noptic)
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  (*pd.coeffs)(nz1(noptic):nz2(noptic)) += c;
  return pd;
}

func zero_coeff(pd,noptic)
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  (*pd.coeffs)(nz1(noptic):nz2(noptic)) *= 0;
  return pd;
}


//func project_to_dms(pd,cond=)
///* DOCUMENT project_to_dms()
// * project_to_dms(): we have mircube, for each non "active" optics OPT:
// * project_to_dms(): compute projection matrices
// * project_to_dms(): fit optics with modes -> coefs
// * project_to_dms(): projects coefs onto DMs, add phase to mircube DM plans
// * project_to_dms(): zero mircube(,,OPT)
// *
// */
//{
//  // before starting:
//  nopt = nof(*pd.alt);
//  active = *pd.active;
//  passive = 1-active;
//  wpassive = where(passive==1); npassive = nof(wpassive);
//  if (npassive==0) {
//    write,format="%s\n","No \"passive\" optics, no projection to do."
//    return mircube;
//  }
//  wactive = where(active==1); nactive = nof(wactive);
//  // write,npassive,nactive;
//  if (nactive==0) {
//    write,format="%s\n","No \"active\" optics, nothing to project on."
//    return mircube;
//  }
//  // precomputations:
//  psize = pd.teldiam/pd.pupd;
//  // active (DMs) characteristics:
//  alt_active = (*pd.alt)(wactive);
//  nmod_dm = (*pd.nmod)(wactive);
//  if (nallof(nmod_dm==nmod_dm(1))) error,"nmodes should be the same for all active optics";
//  nmod_dm = nmod_dm(1);
//
//  // Loop on passive optics. Each in turn will be projected to the active DMs
//  for (no=1;no<=npassive;no++) {
//    ipass = wpassive(no);
//    nmod_opt = (*pd.nmod)(ipass);
//    // calculation of patch ratios via altitude difference and kernel:
//    alt_pass = (*pd.alt)(ipass);
//    delta_alt = abs(alt_active-alt_pass);
//    patch_diam = pd.pupd+2.*max(abs(*pd.xpos,*pd.ypos))*4.848e-6*delta_alt/psize;
//    ratiov = patch_diam/pd.pupd;
//    patch_passive = (*pd.patch_diam)(ipass);
//    // in the calculation of the projector, we compute DM to DM projection too.
//    if (debug>10) write,format="%s\n","Computing projectors passive -> active";
//    P = compute_dms_projector(usemodes,ratiov,nmod_opt,nmod_dm,cond=cond);
//    c = get_coeff(pd,ipass);
//    cdm = P(,+)*c(+);
//    for (na=1;na<=nactive;na++) {
//      pd = add_to_coeff(cdm(1+(na-1)*nmod_dm:na*nmod_dm),pd,wactive(na));
//    }
//    pd = zero_coeff(pd,ipass);
//  }
//  // update resulting mircube:
//  window,3; fma;
//  plh,(*pd.coeffs),color="red";
//  window,2;
//  for (no=1;no<=nopt;no++) {
//    def = get_def(pd,no);
//    coef = get_coeff(pd,no);
//    // write,no,coef;
//    (*pd.mircube)(,,no) = def(,,+)*coef(+);
//    plsys,no*3-1;
//    pli,(*pd.mircube)(,,no)*(*pd.pupil)(,,no);
//    plsys,no*3; pli,array(0.,[2,64,64]);
//  }
//  redraw;
//  return pd;
//}

//func project_to_dms_simple(pd)
///* DOCUMENT project_to_dms()
// * project_to_dms(): we have mircube, for each non "active" optics OPT:
// * project_to_dms(): compute projection matrices
// * project_to_dms(): fit optics with modes -> coefs
// * project_to_dms(): projects coefs onto DMs, add phase to mircube DM plans
// * project_to_dms(): zero mircube(,,OPT)
// *
// */
//{
//  // before starting:
//  nopt = nof(*pd.alt);
//  active = *pd.active;
//  passive = 1-active;
//  wpassive = where(passive==1); npassive = nof(wpassive);
//  if (npassive==0) {
//    write,format="%s\n","No \"passive\" optics, no projection to do."
//    return mircube;
//  }
//  wactive = where(active==1); nactive = nof(wactive);
//  // write,npassive,nactive;
//  if (nactive==0) {
//    write,format="%s\n","No \"active\" optics, nothing to project on."
//    return mircube;
//  }
//  // precomputations:
//  psize = pd.teldiam/pd.pupd;
//  // active (DMs) characteristics:
//  alt_active = (*pd.alt)(wactive);
//  nmod_dm = (*pd.nmod)(wactive);
//  if (nallof(nmod_dm==nmod_dm(1))) error,"nmodes should be the same for all active optics";
//  nmod_dm = nmod_dm(1);
//
//  window,3;
//  // Loop on passive optics. Each in turn will be projected to the active DMs
//  for (no=1;no<=npassive;no++) {
//    ipass = wpassive(no);
//    if (debug>100) { write,format="\nindex passive =%d\n",ipass; }
//    // nmod_opt = (*pd.nmod)(ipass);
//    // calculation of patch ratios via altitude difference and kernel:
//    alt_pass = (*pd.alt)(ipass);
//    if (debug>100) { write,format="altitude passive =%f\n",alt_pass; }
//    // what is the closest active optics?
//    iactive = wactive(wheremin(abs(alt_active-alt_pass))(1));
//    if (debug>100) { write,format="Corresponding index active =%d\n",iactive; }
//    delta_alt = abs((*pd.alt)(iactive)-alt_pass)(1);
//    if (debug>100) { write,format="delta altitude passive -> active =%f\n",delta_alt; }
//    theta = 2*max(abs(*pd.xpos,*pd.ypos))*4.848e-6;
//    kerd = theta*delta_alt/psize; // kernel in pixel
//    if (optim_fov_fact) kerd *= optim_fov_fact;
//    if (debug>100) { write,format="Theta = %f [rd], kernel diameter = %f [pix]\n",theta,kerd; }
//    pup_pass = (*pd.pupil)(,,ipass)(,,1);
//    pup_active = (*pd.pupil)(,,iactive)(,,1);
//    // kerd = 0.1;
//    ker = dist(pd.size)<=kerd;
//    ker = float(ker)/sum(ker);
//    if (debug>100) { tv,ker; pltitle,"kernel"; if (hitReturn()=="q") error; }
//    pha2 = (*pd.mircube)(,,ipass)*pup_pass;
//    if (debug>100) { tv,pha2*pup_pass; pltitle,"phase at passive"; if (hitReturn()=="q") error; }
//    pha2conv = fft_convolve(pha2,ker);
//    if (debug>100) { tv,pha2conv*pup_pass; pltitle,"phase * kernel at passive"; if (hitReturn()=="q") error; }
//    if (debug>100) { tv,pha2conv*pup_active; pltitle,"phase * kernel at passive, active pupil"; if (hitReturn()=="q") error; }
//    if (debug>100) { tv,(*pd.mircube)(,,iactive)(,,1)*pup_active; pltitle,"phase at active"; if (hitReturn()=="q") error; }
//    (*pd.mircube)(,,iactive) += pha2conv;
//    if (debug>100) { tv,(*pd.mircube)(,,iactive)(,,1)*pup_active; pltitle,"Updated phase at active"; if (hitReturn()=="q") error; }
//    (*pd.mircube)(,,ipass) *= 0;
//  }
//  // update resulting mircube:
//  window,2;
//  for (no=1;no<=nopt;no++) {
//    (*pd.mircube)(,,no) *= (*pd.pupil)(,,no);
//    plsys,no*3-1;
//    pli,(*pd.mircube)(,,no);
//    plsys,no*3; pli,array(0.,[2,64,64]);
//  }
//  redraw;
//  return pd;
//}


func project(pd,indfrom,indto)
{
  if (debug>50) write,format="project from %d to %d\n",indfrom,indto;
  dim = pd.size;
  xy = indices(dim); //-dim/2.-0.5; // centred
  xylin = indgen(dim);
  // get the pupils
  // pupfrom = (*pd.pupil)(,,indfrom);
  // pupto = (*pd.pupil)(,,indto);
  // the "pupil" are not necessarily circular but defined by the source geometry:
  pupfrom = (*pd.maskcube)(,,indfrom);
  pupto = (*pd.maskcube)(,,indto);
  // points inside the pupils
  wfrom = where(pupfrom);
  wto = where(pupto);
  // altitudes
  altfrom = (*pd.alt)(indfrom);
  altto = (*pd.alt)(indto);
  // phases
  phafrom = (*pd.mircube)(,,indfrom);
  phato = (*pd.mircube)(,,indto);
  phaproj = pupnproj = array(0.,[2,dim,dim]);
  if (debug>50) {
    window,1; fma;
    plsys,3; pli,phafrom*pupfrom;
  }
  // proj will be used to project pixels from pup from to pup to
  // as in phase_to = proj(,+)*phase_from(+)
  // proj = array(0.,[2,nof(wto),nof(wfrom)]);
  theta = 2.*max(abs(*pd.xpos,*pd.ypos))*4.848e-6;
  psize = pd.teldiam/pd.pupd;
  patches_diam = pd.pupd+theta*abs(*pd.alt)/psize;
  ratio = (patches_diam(indto)-pd.pupd)/(patches_diam(indfrom)-pd.pupd);
  dmax = floor((dim-pd.pupd)/2.);
  pup = *pd.ipupil;
  npup = sum(pup);
  // window,3; tv,phafrom; window,1;
  npt = 0;
  // dmax = 0.5;
  for (i=-dmax;i<=dmax;i++) {
    for (j=-dmax;j<=dmax;j++) {
      // if (i>(-dmax+5)) continue;
      // if (abs(j)>2) continue;
      shiftpup = roll(pup,[i,j]);
      // just to make sure the beam is entirely contained in the "from" optics:
      if (sum(shiftpup*pupfrom)!=npup) {
        if (debug>120) write,format="%s","skip ";
        continue;
      }
      npt++;
      if (debug>120) { fma; plsys,1; pli,shiftpup; }
      // wfrom2 = where(shiftpup);
      // xfrom = xy(,,1)(wfrom2);
      // yfrom = xy(,,2)(wfrom2);
      // xto = (xfrom-dim/2-0.5)
      // shift phafrom to centre
      // phafrom2 = roll(phafrom,[-i,-j])*pupfrom;
      phafrom2 = roll(phafrom,[-i,-j]);
      if (debug>120) { plsys,2; pli,phafrom2+0.05*pup*(max(phafrom2)-min(phafrom2)); }
      // now back to "to", it has to be shifted back by:
      // xylin2 = (xylin-(dim/2+0.5))*ratio+(dim/2+0.5);
      xs = -i*ratio; ys = -j*ratio;
      // phaproj += bilinear(phafrom2,-xs+xylin,-ys+xylin,grid=1);
      phashift = bilinear(phafrom2,xs+xylin,ys+xylin,grid=1,outside=0);
      puptoproj = bilinear(pup,xs+xylin,ys+xylin,grid=1,outside=0);
      phaproj += phashift*puptoproj;
      pupnproj += puptoproj;
      if (debug>120) {
        plsys,3;
        pli,phashift+0.05*(pupto+puptoproj)*(max(phashift)-min(phashift));
        redraw;
        if (hitReturn()=="q") error;
      }
    }
  }
  // phaproj = phaproj/npt;
  phaproj = phaproj/clip(pupnproj,1e-6,);
  if (debug>50) {
    window,1;
    plsys,2;
    pli,phaproj*pupto;
    plsys,3;
    pli,phafrom*pupfrom;
    redraw;
    if (hitReturn()=="q") error;
  }
  // add phase to "to" optic:
  (*pd.mircube)(,,indto) += phaproj;
  // zero "from" optic
  (*pd.mircube)(,,indfrom) *= 0;

  return pd;
}

func simple_project(pd)
{
  // window,2;
  // for (no=1;no<=nopt;no++) {
  //   plsys,no*3-1;
  //   (*pd.mircube)(,,no) *= (*pd.pupil)(,,no);
  //   pli,(*pd.mircube)(,,no);
  //   plsys,no*3; pli,array(0.,[2,64,64]);
  // }
  // redraw;
  // pause,1000;
  // before starting:
  nopt = nof(*pd.alt);
  active = *pd.active;
  passive = 1-active;
  wpassive = where(passive==1); npassive = nof(wpassive);
  if (npassive==0) {
    write,format="%s\n","No \"passive\" optics, no projection to do."
    return mircube;
  }
  wactive = where(active==1); nactive = nof(wactive);
  // write,npassive,nactive;
  if (nactive==0) {
    write,format="%s\n","No \"active\" optics, nothing to project on."
    return mircube;
  }
  // active (DMs) characteristics:
  alt_active = (*pd.alt)(wactive);
  if (debug>100) window,3;
  // Loop on passive optics. Each in turn will be projected to the active DMs
  for (no=1;no<=npassive;no++) {
    ipass = wpassive(no);
    if (debug>100) { write,format="\nindex passive =%d\n",ipass; }
    alt_pass = (*pd.alt)(ipass);
    iactive = wactive(wheremin(abs(alt_active-alt_pass))(1));
    if (debug>100) { write,format="Corresponding index active =%d\n",iactive; }
    pd=project(pd,ipass,iactive);
  }
  if (disp) window,2;
  for (no=1;no<=nopt;no++) {
    (*pd.mircube)(,,no) *= (*pd.pupil)(,,no);
    if (disp) {
      plsys,no*3-1;
      pli,(*pd.mircube)(,,no);
      plsys,no*3;
      pli,array(0.,[2,64,64]);
    }
  }
  if (disp) redraw;

  return pd;
}
