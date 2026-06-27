// config is system config, including defocus values and rotation for some optics
// One could conceptually have different defoc for the various rotation, but
// in here we assume and use the same defocs (e.g. [-1,0,1]) for the various rotation
// configs (e.g. [0,0,0] and [180,0,0]). Rotations config are defined below.

struct pray_struct
{
  pointer images;                 // [size,size,ntarget,nfoc] phase-diversity images (data) fitted by pray()
  pointer original_big_image;     // mosaic (build_bigim) of in-focus images, saved for before/after display
  pointer mircube;                // [size,size,nopt] working phase cube, one map per optic, rebuilt from def(,,+)*coeff(+)
  pointer truecube;               // [size,size,nopt] "true"/reference phase per optic; truecube-mircube = residual
  pointer origcube;               // frozen copy of truecube taken before defocus-removal/projection mutates it
  pointer maskcube;               // [size,size,nopt] 0/1 per-optic illuminated-footprint mask (from ipupil, shifted)
  pointer truecoeffs;             // "true" input modal coefficients used to build truecube, for comparison to coeffs
  pointer coeffs;                 // current modal coefficient vector (length nmod(sum)), indexed via nz1/nz2 per optic
  pointer def;                    // [size,size,nmod(sum)] modal basis ("deformation") maps for all optics, concatenated
  pointer variance;               // noise-variance map/scalar weighting the pray_error data-fidelity term
  pointer norm;                   // coefficient-normalization vector (initialized from guess); effectively unused
  pointer nmod;                   // number of modes per optic; nmod(sum)/nmod(cum) give total/per-optic coeff ranges
  pointer alt;                    // conjugation altitude [m] per optic, used for beam footprint offsets and patch_diam
  pointer active;                 // 0/1 per-optic flag: active (correctable DM) vs passive, used by projection routines
  pointer patch_diam;             // per-optic meta-pupil diameter [pixels] at its conjugate altitude
  pointer ipupil;                 // [size,size] system pupil mask (make_pupil()) at ground conjugate
  pointer pupil;                  // [size,size,nopt] per-optic pupil/meta-pupil footprint, masks mircube before PSF calc
  pointer dmgsxposcub;            // [_n,nopt,ntarget] x-coord [pixels] of each target's beam footprint on each optic
  pointer dmgsyposcub;            // [_n,nopt,ntarget] y-coord, companion to dmgsxposcub
  pointer deltafoc;               // extra-focal distances (e.g. [-1.,0.,1.], replicated per rotation config)
  pointer xpos;                   // target x-positions [arcsec] in the FoV, from init_target_positions()
  pointer ypos;                   // target y-positions [arcsec], companion to xpos
  pointer focus;                  // [size,size] normalized defocus map, added to mircube for extra-focal images
  pointer ampli_pup;              // [size,size,ntarget] complex pupil-plane amplitude per target, set by compute_psfs()
  pointer ampli_foc;              // [size,size,ntarget] complex focal-plane amplitude per target, set by compute_psfs()
  pointer xy4centring;            // indices(size)-size/2-0.5, pixel grid used for image centroid computation
  pointer config;                 // array of config_struct (one per focus/rotation combo): .foc and .roti (-> rotv)
  pointer rotv;                   // [nopt,nrot] per-optic rotation angle [deg] per rotation config, applied via rotate2()
  string  fovshape;               // "square" or "round": target-grid footprint shape, also used by Strehl-contour plots
  long    disp;                   // 0/1: whether to open display windows and plot during fitting/projection
  float   fullfield;              // full FoV diameter [arcsec] spanned by the target grid
  float   avg_focus;              // average focus coefficient removed from truecube (so 0-defocus images are in-focus)
  float   shift_foc;              // intended focus-shift parameter for pray(); effectively dead/unused
  float   peak_airy;              // normalized peak of the diffraction-limited (Airy) PSF; Strehl-ratio denominator
  float   teldiam;                // telescope diameter [m]; psize = teldiam/pupd gives the pixel scale
  float   cobs;                   // central obscuration fraction, used by make_pupil() to build ipupil
  float   objective_function;     // final value of the minimized criterion (fout from optm_vmlmb) after pray()
  long    pupd;                   // pupil diameter [pixels] (illuminated disk size within the size x size array)
  long    size;                   // side length [pixels] of all [size,size,...] phase/pupil/image arrays
  long    ngrid;                  // number of targets per side of the FoV sampling grid (e.g. 4-8)
  long    centre;                 // array centre coordinate (size/2+0.5), used as xc/yc for pupil/beam indexing
  long    _n;                     // half-width-derived size of the local sub-array used to extract per-direction phase
  long    _n1;                    // lower pixel index of the _n x _n sub-array window within the full array
  long    _n2;                    // upper pixel index of that sub-array window (paired with _n1)
  pointer _def_pup;               // [nof(*_pupw),ncoeffs,ntarget,nrot] per-direction modal basis, masked to _pupw pixels
  pointer _pupw;                  // index vector of illuminated (nonzero ipupil) pixels, used to mask _def_pup etc.
  pointer _coeffs0;               // pristine pre-projection snapshot of coeffs, see simple_projection_only(reset=)
  pointer _mircube0;              // idem, snapshot of mircube
  pointer _ftobject;              // fft(object,1), precomputed once and reused in pray_error()/pray_j_data()
};

struct sim_struct { long verbose; }

struct config_struct {
  float foc;
  long roti;
}
