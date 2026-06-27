usemodes  = "zer";        // "dh" (recommended, works well), "kl" or "zer"
geometry  = "square";     // "square", "hexagonal" or "mavis"
// geometry  = "hexagonal";
geometry  = "mavis";
// fovshape  = "round";     // "round" if desired if not will default to square
fovshape  = "square";     // "round" if desired if not will default to square
// initphase = "coefs";   // "coefs" or "screens"
initphase = "screens";   // "coefs" or "screens"

if (case==[]) error,"you have to set the variable 'case' (see mavis_pray_conf.i)";

// parameters defined statically:
if (case==1) {
  alt     = [0.];
  nmod    = [200];
  nm_rmsv = [80]; // has to be defined if initphase = "screens"
  fit     = [1];
  active  = [1];
  rotv    = [[0.]];
}

if (case==2) {
  // the blow works well with
  // random_seed,0.75; res=mavis_pray(,5,[0.,-1.5,1.5,-2.5,2.5],100000,1,,disp=1,maxiter=50,modes="dh")
  // 99.1% Strehl
  alt     = [0.,4000];
  nmod    = [100,100];
  nm_rmsv = [50,50]; // has to be defined if initphase = "screens"
  fit     = [1,1];
  active  = [1,0];
  rotv    = [[0.,0.]];
}

if (case==3) {
   alt     = [0.,13500.]; //45000 seems to be the limit
   nmod    = [150,150];
   nm_rmsv = [30,30]*1.7; // has to be defined if initphase = "screens"
   fit     = [1,1];
   active  = [1,0];
   rotv    = [[0.,0],[180.,0]];
}

if (case==4) {
  alt     = [0.,6000,13500.]; //45000 seems to be the limit
  nmod    = [50,50,50]*2;
  nm_rmsv = [30,50,30]*1.2; // has to be defined if initphase = "screens"
  fit     = [1,1,1];
  active  = [1,1,1];
  rotv    = [[0.,0,0],[180,0,0]];
  // rotv = [[0.,0,0],[120,120,0],[240,240,0]];
}

if (case==5) {
  // To test projection:
  alt     = [8000.,4000.,0.,-4000];
  nmod    = [100,100,100,100];
  nm_rmsv = [50,50,50,50]; // has to be defined if initphase = "screens"
  fit     = [1,1,1,1];
  active  = [0,1,1,0];
  // rotv    = [[0.,0,0,0];
  rotv    = [[0.,0,0,0],[180.,180,0,0]];
}

if (case==6) {
  // To test projection:
  alt     = [8000.,4000.,0.,-8000];
  nmod    = [150,150,150,150];
  nm_rmsv = [50,50,50,50]; // has to be defined if initphase = "screens"
  fit     = [1,1,1,1];
  active  = [0,1,1,0];
  // rotv    = [[0.,0,0,0];
  rotv    = [[0.,0,0,0],[180.,180,0,0]];
}

if (case==9) { // updated collimator 6/4/26
  alt         = [45.5,13.6,6   ,1.2 ,0. ,-1.9,-4 ,-12.4,-23.9,-29.9]*1000; // altitude of optics, length nopt
  // nm_rmsv     = [10. ,30  ,30  ,10  ,30 ,47  ,9  ,11.0 ,6.9   ,48  ];
  nm_rmsv     = [10. ,25  ,25  ,10  ,25 ,47  ,9  ,11.0 ,6.9   ,48  ];
  nm=100; nmod = [nm  , 100 , 100 ,nm  , 100,nm  ,nm ,nm   ,nm    ,nm  ]; // number of modes per optics
  fit         = [0   ,1   ,1   ,0   ,1  ,0   ,0  ,0    ,0     ,0   ];
  active      = [0   ,1   ,1   ,0   ,1  ,0   ,0  ,0    ,0     ,0   ];
  rotv        = [[0. ,0   ,0   ,0   ,0  ,0   ,0  ,0    ,0     ,0   ],
                [180 ,180 ,180 ,180 ,180,180 ,180,90   ,0     ,0   ]];
}

if (case==10) { // updated collimator 6/4/26
  alt         = [45.5,13.6,6   ,1.2 ,0. ,-1.9,-4 ,-12.4,-23.9,-29.9]*1000; // altitude of optics, length nopt
  // nm_rmsv     = [10. ,30  ,30  ,10  ,30 ,47  ,9  ,11.0 ,6.9   ,48  ];
  nm_rmsv     = [10. ,25  ,25  ,10  ,25 ,47  ,9  ,11.0 ,6.9   ,48  ];
  nm=50; nmod = [nm  , 100 , 100 ,nm  , 100,nm  ,nm ,nm   ,nm    ,nm  ]; // number of modes per optics
  fit         = [0   ,1   ,1   ,0   ,1  ,0   ,0  ,0    ,0     ,0   ];
  active      = [0   ,1   ,1   ,0   ,1  ,0   ,0  ,0    ,0     ,0   ];
  rotv        = [[0. ,0   ,0   ,0   ,0  ,0   ,0  ,0    ,0     ,0   ],
                [180 ,180 ,180 ,180 ,180,180 ,180,90   ,0     ,0   ]];
  // rotv        = [[0. ,0   ,0   ,0   ,0  ,0   ,0  ,0    ,0     ,0   ]];
  fit       = fit*0+1;
}

if (case==11) { // updated collimator 6/4/26
  alt         = [15.5,13.6,6   ,1.2 ,0. ,-1.9,-4 ,-6.4,-9.9,-23.9]*1000; // altitude of optics, length nopt
  // nm_rmsv     = [10. ,30  ,30  ,10  ,30 ,47  ,9  ,11.0 ,6.9   ,48  ];
  nm_rmsv     = [10. ,25  ,25  ,10  ,25 ,47  ,9  ,11.0 ,6.9   ,48  ];
  nm=100; nmod = [nm  , 100 , 100 ,nm  , 100,nm  ,nm ,nm   ,nm    ,nm  ]; // number of modes per optics
  fit         = [0   ,1   ,1   ,0   ,1  ,0   ,0  ,0    ,0     ,0   ];
  active      = [0   ,1   ,1   ,0   ,1  ,0   ,0  ,0    ,0     ,0   ];
  rotv        = [[0. ,0   ,0   ,0   ,0  ,0   ,0  ,0    ,0     ,0   ],
                [180 ,180 ,180 ,180 ,180,180 ,180,90   ,0     ,0   ]];
  fit       = fit*0+1;
}


if (case==100) { // full full case
  // updated collimator 6/4/26
  //updated march 30 2026:
  alt       = [45.5,13.6,6  ,1.2,0.,-1.9,-3.3,-4.4,-8.3,-12.4,-16.5,-23.9,-29.9,-36.2]*1000; // altitude of optics, length nopt
  nm_rmsv   = [10. ,30  ,30 ,10 ,30,47  ,6.4 ,6.4 ,6.6 ,6.6  ,6.6  ,6.9   ,48   ,5];
  nmod      = [50,100,100,50,100,50,50,50,50,50,50,50,50,50]; // number of modes per optics
  fit       = [0,1,1,0,1,0,0,0,0,0,0,0,0,0];
  rotv      = [[0.,0,0,0,0,0,0,0,0,0,0,0,0,0],\
               [180,180,180,180,180,180,180,180,180,90,0,0,0,0]]; // rotation of optics, as many line as configs
  fit       = fit*0+1;

}

// just for res=get_non_normalised_strehls()
// alt       = [45.5,13.6,6   ,1.2 ,0. ,-1.9,-4 ,-12.4,-23.9 ,-29.9]*1000; // altitude of optics, length nopt
// nm_rmsv   = [10. ,30  ,30  ,30  ,30 ,47  ,9  ,11.0 ,6.9   ,48   ];
// nmod      = [2   ,2   ,2   ,2   ,2  ,2   ,2  ,2    ,2     ,2    ]; // number of modes per optics
// fit       = [0   ,1   ,1   ,0   ,1  ,0   ,0  ,0    ,0     ,0    ];
// rotv      = [[0.  ,0   ,0   ,0   ,0  ,0   ,0  ,0    ,0     ,0    ]];
// fit       = fit*0+1;

// reading optics specs from optics_data_alt_wfe.i: FULL OPTICS TRAIN
if (case==0) {
  require,"optics_data_alt_wfe.i";
  status = read_optics_data("optics_data_alt_wfe.i",opt_name,opt_alt,opt_wfe);
  nopt       = numberof(opt_name);
  wdm        = where(strmatch(opt_name,"DM"));
  wkm        = where(strmatch(opt_name,"K-Mirr2"))(1);
  alt        = opt_alt*1000; // altitude of optics, length nopt
  nm_rmsv    = opt_wfe;
  nmod = fit = active = array(0,nopt); rot1 = rot2 = array(0.,nopt);
  nmod       = nmod*0+50; nmod(wdm) = 100;
  fit(wdm)   = 1;
  rot2(1:wkm-1) = 180; rot2(wkm) = 90;
  rotv       = [rot1,rot2];
  fit        = fit*0+1;
  active(wdm) = 1;
}

// w = where(fit==0); if (nof(w)) nmod(w) = 2;
if (active==[]) active = fit*0+1;

weight             = array(2./sqrt(nof(alt)),nof(nmod)); // mode weights (static aberrations)
fullfield          = 30.; // full field in arcsec (on the side) - why 40 and not 30?
teldiam            = 8.0; // telescope diameter
lambda             = 550; // walevength [nm]
cobs               = 0.;  // central obstruction
size               = 64;  // side dimension of small PSF arrays
osampl             = 1;   // oversampling (1 or 2)
gridpad            = 1.;  // padding in arcsec to validate sources
ps_slope           = -2.5; // slope of power spectrum for init phase screens
modes_slope        = -1.5; //-1.5; // weight vs radial order for mode/defs
zoomfactor         = 3;
display            = 1;
// debug              = 1;
strehl_target      = 0.41;
strehl_target      = exp(-(2*pi*sqrt(sum(nm_rmsv^2))/lambda)^2);
strehl_normalise   = 10; // number of iterations for Strehl normalisation (recommended: 2)
centre_init_images = 0; // centre original images.
centre_pray_images = 0; // centre modelled images. DOES NOT WORK
nmoddm = nmod(where(active)(1));
projection_method  = "simple"; // "simple" or "optimal"
dm_proj_cond       = 15.0; // conditioning number fo dm projection
dm_proj_tikhonov   = 1.0; // tikhonov regularisation parameter for dm projection
proj_cond          = sqrt(nmoddm)/8.; // condition number for projection to DMs
// proj_cond          = 2.0; // condition number for projection to DMs
// Graphics parameters
dpi_target       = 160; // dpi for the "large" graphic windows
dpi_target_small = 100; // dpi for the secondary graphic windows
pltitle_height   = 14;


// CHECKS
if (fit==[]) fit = long(alt*0+1);
if (rotv==[]) rotv = alt*0.;

if ((initphase=="screens")&(nm_rmsv==[])) error,"nm_rmsv undefined while initphase=\"screens\"";
if ((initphase=="coefs")&(weight==[])) error,"nm_rmsv undefined while initphase=\"screens\"";
doa = nof(alt);
if (nof(nmod)!=doa) error,"nmod and alt do not have the same dimension";
if (nof(nm_rmsv)!=doa) error,"nm_rmsv and alt do not have the same dimension";
if (nof(fit)!=doa) error,"fit and alt do not have the same dimension";
if (dimsof(rotv)(2)!=doa) error,"rotv dimensions incompatible with alt dimension";

// For optics configuration, WFE and altitude, see optics_data_alt_wfe.i

/* FIXME
 * There is an issue. Choose osmapl=2, and a large altitude, e.g 80000.
 * with [1,30]*1.5 on the [0,87000]. spots have quite a bit of TT decorrelated, which is expected.
 * But now when the thing minimize, the spots are NOT with TT, i.e. they seem on a regular grid.
 * all the while, at the end the spots are perfecly on a grid. What's happening?
 */
