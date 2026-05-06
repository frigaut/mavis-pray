/* projection.i
Computes projection matrices from modes at a given altitude onto modes
at another altitude.
Based on the ratio diameters, which can be readily computed from the FoV
and delta altitude for our NCPA problem.
Modes are either zernike, DH or KL.
*/
require,"yao.i";
sim = sim_struct(); // for call to make_diskharmonic()
sim.verbose = 0;
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
  } else if (modes=="kl") {
    require,"yaokl.i";
    v = obas = pup1 = []; // goes around a make_kl/yorick bug
    kl = make_kl(nmodes,patchd,v,obas,pup1,oc=0.0,nr=128,verbose=0);
    // kl = order_kls(kl,patchd,upto=20);
    off = (dim-patchd)/2; i1 = 1+off; i2 = i1+patchd-1;
    mod(i1:i2,i1:i2,) = kl;
    pupil(i1:i2,i1:i2) = pup1;
  } else if (modes=="dh") {
    require,"yaokl.i";
    mod = make_diskharmonic(dim,patchd,nmodes+1,cobs=0.,xc=dim/2+0.5,yc=dim/2+0.5)(,,2:);
    mod(,,[3,4]) = mod(,,[4,3]); // switch astig45 and focus
    w = where(mod(,,min)==0);
    pupil += 1; pupil(w)=0;
    mod = mod*pupil(,,-); // fixes non zero pixels outside pupil for azimuthal order=0 modes
  }
  return mod;
}

func proj_modes_from_to(modes,ratio,nmod1,nmod2)
/* DOCUMENT proj_modes_from_to(modes,ratio,nmod1,nmod2)
 * Compute the projection matrix from modes at plan2 to modes at plan1
 * Plan2 is the largest patch (diameter)
 * modes: "zer","kl" or "dh"
 * ratio is the >1 ratio of patches = patchd2/patchd1
 * nmod1 and nmod2 are the number of modes in plans 1 and 2.
 * Returns p2on1 (projection of 2 on 1) array [2,nmod1,nmod2]
 * to get coefficients on plan1 knowing coefficient on plan2 coef2:
 * coef1 = p2on1(,+)*coef2(+)
 */
{
  if (ratio<1) error,"ratio has to be >=1";
  dim = 64;
  patchd2 = dim;
  patchd1 = lround(dim/ratio);
  p2on1 = array(0.,[2,nmod1,nmod2]);
  // can't make modes with fractional patchd. will have to up idm and patchd to do this
  // calculation to reduce effect of rounding errors
  // kernel diameter is patch2-patch1
  ker = dist(dim)<=((patchd2-patchd1)/2.);
  ker = float(ker)/sum(ker);
  mod1 = generate_modes(modes, nmod1, dim, patchd1, pupil);
  pupil1 = pupil;
  mod2 = generate_modes(modes, nmod2, dim, patchd2, pupil);
  pupil2 = pupil;
  // tv,fft_convolve(mod2(,,13),ker)*pupil1;
  // modes in plan1 (smallest patch) to determine projection of
  // phase to modes1:
  wpup1 = where(pupil1);
  mod1lin = mod1(*,)(wpup1,);
  mtm = mod1lin(+,)*mod1lin(+,);
  proj = LUsolve(mtm)(+,)*mod1lin(,+);
  // coef = proj(,+)*mod1lin(+,4); // check, all good
  for (nm2=1;nm2<=nmod2;nm2++) {
    mode2 = mod2(,,nm2);
    cm = (fft_convolve(mod2(,,nm2),ker)*pupil1)(*)(wpup1);
    p2on1(,nm2) = proj(,+)*cm(+);
    // if (debug>100) {
    //   ar2 = mod2(,,nm2);
    //   ar2on1 = fft_convolve(mod2(,,nm2),ker)*pupil1;
    //   ar1 = p2on1(,nm2)(+)*mod1(,,+);
    //   tv,_(ar2,ar2on1,ar1); limits,square=1;
    //   if (hitReturn()=="q") return p2on1;
    // }
  }
  return p2on1;
}

func proj_to_dms(modes,ratiov,nmod_opt,nmod_dm,cond=)
/* DOCUMENT proj_to_dms(modes,ratiov,nmodv)
 * Computes the optimal projection matrix from one plan (plano=plan object)
 * to several DMs.
 * modes: "zer","kl" or "dh"
 * ratiov: vector of ratio of DMs to plano. All elements have to be >= 1
 * and in creasing order.
 * example: [1.5, 1.8]. See proj_modes_from_to. In this example, the ratio
 * between plano and the DMs (2 DMs) are 1.5 and 1.8. The ratio between the proj_to_dms
 * themselves will be cmoputed from this. numberof(ratiov) = number of DMs
 * nmod_opt: number of modes for plano optics
 * nmod_dm: number of modes for DMs (has to be the same)
 */
{
  if (cond==[]) cond=100.;
  pltitle_height = 14;
  // For a 3-DM system at altitudes h_dm = [h1, h2, h3]
  // and one optic at h_opt with nmod_opt modes
  // each DM has nmod_dm modes (could differ; let's say all the same)

  // === Build the joint Gramian G_tilde ===
  // Diagonal blocks: identity if modes are orthonormal on each DM's meta-pupil
  // Off-diagonal blocks: inter-DM cross-projection
  ndm = numberof(ratiov);
  G_tilde = array(0., [2,ndm*nmod_dm,ndm*nmod_dm]);
  for (i=1;i<=ndm;i++) {
    // diagonal block (assuming orthonormal modes on DM i's meta-pupil)
    G_tilde(1+(i-1)*nmod_dm:i*nmod_dm, 1+(i-1)*nmod_dm:i*nmod_dm) = unit(nmod_dm);
    for (j=i+1;j<=ndm;j++) {
      // ratio of meta-pupil diameters between DM i and DM j
      ratio_ij = ratiov(j)/ratiov(i);
      P_ij = proj_modes_from_to(modes, ratio_ij, nmod_dm, nmod_dm);
      // G_tilde(1+(j-1)*nmod_dm:j*nmod_dm, 1+(i-1)*nmod_dm:i*nmod_dm) = transpose(P_ij);
      // G_tilde(1+(i-1)*nmod_dm:i*nmod_dm, 1+(j-1)*nmod_dm:j*nmod_dm) = P_ij;
      G_tilde(1+(j-1)*nmod_dm:j*nmod_dm, 1+(i-1)*nmod_dm:i*nmod_dm) = P_ij;
      G_tilde(1+(i-1)*nmod_dm:i*nmod_dm, 1+(j-1)*nmod_dm:j*nmod_dm) = transpose(P_ij);
    }
  }
  if (debug>100) { tv,G_tilde; pltitle,escapechar("G_tilde"); pause,500; }
  // Build the RHS M_tilde
  // For one optic, stacked over DMs
  M_tilde = array(0.,[2,ndm*nmod_dm,nmod_opt]);
  for (i=1;i<=ndm;i++) {
    M_tilde(1+(i-1)*nmod_dm:i*nmod_dm, ) = proj_modes_from_to(modes, ratiov(i), nmod_dm, nmod_opt);
  }
  if (debug>100) { tv,M_tilde; pltitle,escapechar("M\_tilde"); pause,500; }
  // Solve
  ev = SVdec(G_tilde,u,vt);
  evi = 1./ev;
  w = where(ev<max(ev/cond));
  if (nof(w)) evi(w) = 0.;
  G_tilde_inv = (transpose(vt)(,+)*(diag(evi))(+,))(,+)*transpose(u)(+,);
  P_joint = G_tilde_inv(,+) * M_tilde(+,);
  if (debug>100) { tv,P_joint; pltitle,escapechar("P_joint"); pause,500; }
  // Apply: c_dm_stacked = P_joint(,+) * c_opt(+)
  //
  // CHECKS
  if (debug>90) {
    if (ndm!=2) error,"Check only coded for 2 DMs";
    dim = 64;
    mod_opt = generate_modes(modes, nmod_opt, dim, dim, pupil);
    pupil_opt = pupil;
    mod_dm1 = generate_modes(modes, nmod_dm, dim, lround(dim/ratiov(1)), pupil);
    pupil_dm1 = pupil;
    mod_dm2 = generate_modes(modes, nmod_dm, dim, lround(dim/ratiov(2)), pupil);
    pupil_dm2 = pupil;

    if (optmode==[]) optmode=13; else optmode=long(optmode);
    c_opt = array(0.,nmod_opt); c_opt(optmode)=1.;
    dmc = P_joint(,+) * c_opt(+);
    pha_opt = mod_opt(,,+)*c_opt(+);
    pha_dm1 = mod_dm1(,,+)*dmc(1:nmod_dm)(+);
    pha_dm2 = mod_dm2(,,+)*dmc(1+nmod_dm:2*nmod_dm)(+);
    tv,transpose(_(pha_opt,pha_dm1,pha_dm2));
    limits,square=1;
  }
  return P_joint;
}
