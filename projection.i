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
/* DOCUMENT remove_mode(pha,pupil,mode)
 * Least-squares removes the projection of "mode" from phase "pha" over
 * the area defined by "pupil", and returns the result re-masked by pupil.
 * e.g. used to remove the average focus mode from the truth phase cubes.
 */
{
  c = sum(pha*pupil*mode)/sum(mode*mode*pupil);
  pha = (pha-c*mode)*pupil;
  return pha;
}

// get_def/get_coeff/add_to_coeff/zero_coeff: small accessors into pd.def/
// pd.coeffs for a single optic "noptic", using the per-optic column ranges
// derived from pd.nmod (same indexing as the nz1/nz2 used throughout
// mavis_pray.i/pray.i). Used by optimal_project_modal()/simple_project()
// to read/update one optic's modes or coefficients without touching the
// others.
func get_def(pd,noptic)
/* DOCUMENT get_def(pd,noptic)
 * Returns the [size,size,nmod] modal basis maps for optic "noptic".
 */
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  return (*pd.def)(,,nz1(noptic):nz2(noptic));
}

func get_coeff(pd,noptic)
/* DOCUMENT get_coeff(pd,noptic)
 * Returns the current modal coefficient vector for optic "noptic".
 */
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  return (*pd.coeffs)(nz1(noptic):nz2(noptic));
}

func add_to_coeff(c,pd,noptic)
/* DOCUMENT add_to_coeff(c,pd,noptic)
 * Adds modal coefficient vector "c" to optic "noptic"'s coefficients
 * in-place (via pd.coeffs's pointer). Returns pd.
 */
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  (*pd.coeffs)(nz1(noptic):nz2(noptic)) += c;
  return pd;
}

func zero_coeff(pd,noptic)
/* DOCUMENT zero_coeff(pd,noptic)
 * Zeroes optic "noptic"'s modal coefficients in-place. Returns pd.
 */
{
  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);
  (*pd.coeffs)(nz1(noptic):nz2(noptic)) *= 0;
  return pd;
}


func project(pd,indfrom,indto)
/* DOCUMENT project(pd,indfrom,indto)
 * Real-space projection of one optic's phase onto another, used by
 * simple_project() to send a passive optic's full aberration onto its
 * single nearest active optic (no modal decomposition, no joint fit
 * across multiple DMs -- see optimal_project_modal() for that).
 *
 * Optics at different conjugate altitudes see the same field of view
 * through meta-pupils of different diameter (patches_diam, set by the
 * altitude difference and the FoV half-angle), so a uniform shift of the
 * beam footprint on optic indfrom corresponds to a *scaled* shift on
 * optic indto. This function walks the beam footprint over the maximal
 * range of FoV-induced offsets (+-dmax pixels, the full extent allowed by
 * size vs. pupd), and for each offset where the footprint stays fully
 * inside indfrom's pupil (maskcube), bilinearly resamples the phase from
 * indfrom's frame into indto's frame (scaled by "ratio") and accumulates
 * it, weighted by the resampled pupil itself (puptoproj) so that pixels
 * covered by more FoV offsets are weighted more, and the result is
 * normalized by the accumulated weight (pupnproj) at the end.
 *
 * Mutates pd.mircube in place: indto's phase gains the projected
 * contribution, indfrom's phase is zeroed. Returns pd.
 */
{
  if (debug>50) write,format="project from %d to %d\n",indfrom,indto;
  dim = pd.size;
  xylin = indgen(dim);
  // the "pupil" are not necessarily circular but defined by the source geometry:
  pupfrom = (*pd.maskcube)(,,indfrom);
  pupto = (*pd.maskcube)(,,indto);
  phafrom = (*pd.mircube)(,,indfrom);
  phaproj = pupnproj = array(0.,[2,dim,dim]);
  if (debug>50) {
    window,1; fma;
    plsys,3; pli,phafrom*pupfrom;
  }
  theta = 2.*max(abs(*pd.xpos,*pd.ypos))*4.848e-6;
  psize = pd.teldiam/pd.pupd;
  patches_diam = pd.pupd+theta*abs(*pd.alt)/psize;
  ratio = (patches_diam(indto)-pd.pupd)/(patches_diam(indfrom)-pd.pupd);
  dmax = floor((dim-pd.pupd)/2.);
  pup = *pd.ipupil;
  npup = sum(pup);
  npt = 0;
  for (i=-dmax;i<=dmax;i++) {
    for (j=-dmax;j<=dmax;j++) {
      shiftpup = roll(pup,[i,j]);
      // just to make sure the beam is entirely contained in the "from" optics:
      if (sum(shiftpup*pupfrom)!=npup) {
        if (debug>120) write,format="%s","skip ";
        continue;
      }
      npt++;
      if (debug>120) { fma; plsys,1; pli,shiftpup; }
      // shift phafrom to centre
      phafrom2 = roll(phafrom,[-i,-j]);
      if (debug>120) { plsys,2; pli,phafrom2+0.05*pup*(max(phafrom2)-min(phafrom2)); }
      // now back to "to", it has to be shifted back by:
      xs = -i*ratio; ys = -j*ratio;
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
/* DOCUMENT simple_project(pd)
 * For every passive optic, finds its single nearest active optic (DM) by
 * conjugate altitude and sends the passive optic's entire phase there via
 * project() (real-space, one-to-one). Each passive optic projects
 * independently onto whichever DM is closest in altitude -- no joint fit
 * across multiple DMs (see optimal_project_modal() for that), so a
 * passive optic between two equidistant DMs gets dumped entirely onto one
 * of them rather than split between both.
 *
 * Mutates pd.mircube in place (via repeated calls to project()), then
 * re-masks every optic's mircube by its own pupil. Returns pd.
 *
 * SEE ALSO: optimal_project_modal (joint multi-DM least-squares
 * alternative, generally better Strehl at the cost of higher DM stroke)
 */
{
  nopt = nof(*pd.alt);
  active = *pd.active;
  passive = 1-active;
  wpassive = where(passive==1); npassive = nof(wpassive);
  if (npassive==0) {
    write,format="%s\n","No \"passive\" optics, no projection to do."
    return mircube;
  }
  wactive = where(active==1); nactive = nof(wactive);
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

func optimal_project_modal(pd,cond=,tikhonov=,report=)
/* DOCUMENT optimal_project_modal(pd,cond=,tikhonov=,report=)
 * Optimal multi-DM projection of all passive optics' modal coefficients
 * onto the active (DM) optics, as a joint FoV least-squares fit. Unlike
 * simple_project(), which sends each passive optic's full aberration to
 * its single nearest active optic, this distributes the correction across
 * *all* active optics simultaneously, weighted by how well each one can
 * reproduce the passive aberration over the actual field of view.
 *
 * The fit is exact (no disk-kernel/continuous-FoV approximation): it sums
 * over the real, discrete set of target directions and rotation configs,
 * reusing pd._def_pup (already precomputed once in init_defs, and already
 * restricted to the valid pupil pixels) as the per-direction ray-tracing
 * operator R_k from modal coefficients (all optics stacked) to pupil-plane
 * phase.
 *
 * Let c_dm and c_opt be the stacked modal coefficients of the active and
 * passive optics respectively, and Rdm_k / Ropt_k the corresponding column
 * blocks of R_k for direction/rotation k. The joint least-squares problem
 *   minimize sum_k || Rdm_k c_dm - Ropt_k c_opt ||^2
 * gives the normal equations
 *   G c_dm = M c_opt,  G = sum_k Rdm_k^T Rdm_k,  M = sum_k Rdm_k^T Ropt_k
 * solved here via regularized SVD (cond= conditioning number cutoff, as
 * used for compute_dms_projector() historically) plus an optional Tikhonov
 * term (tikhonov=, in units of the mean diagonal of G) to directly trade
 * fit quality for DM stroke.
 *
 * Updates pd.coeffs: zeroes the passive optics' coefficients and adds the
 * jointly-optimal contribution to the active (DM) coefficients. Does NOT
 * update pd.mircube -- call compute_psfs() with fromscreens=0 afterwards
 * to regenerate it from the updated coefficients.
 *
 * cond=     SVD conditioning number cutoff for the joint Gramian. Defaults
 *           to dm_proj_cond if set (e.g. from mavis_pray_conf.i), else 15
 *           -- found by scanning the Strehl-vs-stroke tradeoff for the
 *           nominal MAVIS configuration; re-tune if nopt/nmod/DM altitudes
 *           change significantly.
 * tikhonov= additional ridge term added to the Gramian before inversion,
 *           as a fraction of its mean diagonal value. Defaults to
 *           dm_proj_tikhonov if set, else 1 (see cond=).
 * report=   if set, print per-DM coefficient RMS/peak after projection,
 *           for comparison against simple_project()
 *
 * SEE ALSO: simple_project (alternative, real-space nearest-DM-only projection)
 */
{
  if (cond==[]) cond = (dm_proj_cond?dm_proj_cond:15.);
  if (tikhonov==[]) tikhonov = (dm_proj_tikhonov?dm_proj_tikhonov:1.);

  active = *pd.active;
  wdm  = where(active==1);  ndm  = nof(wdm);
  wopt = where(active==0);  nopt_passive = nof(wopt);
  if (nopt_passive==0) {
    write,format="%s\n","No \"passive\" optics, no projection to do.";
    return pd;
  }
  if (ndm==0) {
    write,format="%s\n","No \"active\" optics, nothing to project on.";
    return pd;
  }

  nz12 = (*pd.nmod)(cum); nz1 = (nz12+1)(1:-1); nz2 = nz12(2:);

  dmcols = optcols = [];
  for (i=1;i<=ndm;i++) grow,dmcols,indgen(nz1(wdm(i)):nz2(wdm(i)));
  for (i=1;i<=nopt_passive;i++) grow,optcols,indgen(nz1(wopt(i)):nz2(wopt(i)));

  ndmcoef  = nof(dmcols);
  noptcoef = nof(optcols);

  szdp    = dimsof(*pd._def_pup);
  ntarget = szdp(4);
  nrotd   = szdp(5);

  G = array(0.,[2,ndmcoef,ndmcoef]);
  M = array(0.,[2,ndmcoef,noptcoef]);

  for (n=1;n<=nrotd;n++) {
    for (i=1;i<=ntarget;i++) {
      Rdm  = (*pd._def_pup)(,dmcols,i,n);
      Ropt = (*pd._def_pup)(,optcols,i,n);
      G += Rdm(+,)*Rdm(+,);
      M += Rdm(+,)*Ropt(+,);
    }
  }

  if (tikhonov) {
    meandiag = sum(G*unit(ndmcoef))/ndmcoef;
    G += (tikhonov*meandiag)*unit(ndmcoef);
  }

  ev = SVdec(G,u,vt);
  evi = 1./ev;
  w = where(ev<(max(ev)/cond));
  if (nof(w)) evi(w) = 0.;
  Ginv = (transpose(vt)(,+)*(diag(evi))(+,))(,+)*transpose(u)(+,);

  Pjoint = Ginv(,+)*M(+,); // [ndmcoef, noptcoef]

  c_opt = (*pd.coeffs)(optcols);
  (*pd.coeffs)(dmcols)  += Pjoint(,+)*c_opt(+);
  (*pd.coeffs)(optcols)  = 0.;

  if (report) {
    write,format="%s\n","Per-DM coefficient RMS/peak after optimal projection:";
    for (i=1;i<=ndm;i++) {
      c = get_coeff(pd,wdm(i));
      write,format="  optic %d: rms=%.4g peak=%.4g\n",wdm(i),c(rms),max(abs(c));
    }
  }

  return pd;
}
