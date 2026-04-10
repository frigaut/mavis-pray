usemodes  = "zer";        // "dh" (recommended, works well), "kl" or "zer"
geometry  = "hexagonal"; // "square" or "hexagonal"
fovshape  = "round";     // "round" if desired if not will default to square
// fovshape  = "square";     // "round" if desired if not will default to square
initphase = "coefs";   // "coefs" or "screens"
// initphase = "screens";   // "coefs" or "screens"

// parameters defined statically:
alt     = [0.,6000,13500.]; //45000 seems to be the limit
nzer    = [50,50,50];
nm_rmsv = [30,50,30]*1.2; // has to be defined if initphase = "screens"
fit     = [1,1,1];
rotv    = [[0.,0,0],[180,0,0]];
// rotv = [[0.,0,0],[120,120,0],[240,240,0]];

// alt     = [0.];
// nzer    = [8];
// nm_rmsv = [80]; // has to be defined if initphase = "screens"
// fit     = [1];
// rotv    = [[0.]];

// alt     = [0.,13500.]; //45000 seems to be the limit
// nzer    = [150,150];
// nm_rmsv = [30,30]*1.7; // has to be defined if initphase = "screens"
// fit     = [1,1];
// rotv    = [[0.,0],[180.,0]];

// alt       = [-25500,-6000,0,6000,13500.,20000]; // altitude of optics, length nopt
// nzer      = [50,50,100,100,100,50]; // number of modes per optics
// fit       = [0,0,1,1,1,0];
// nm_rmsv   = [15,15,30,30,30,15]*1.5; // has to be defined if initphase = "screens"
// rotv      = [[0.,0,0,0,0,0],\
//              [180,180,0,0,0,0]]; // rotation of optics, as many line as configs

// alt       = [-25500,0,13500.,20000]; // altitude of optics, length nopt
// nzer      = [150,150,150,150]; // number of modes per optics
// fit       = [0,1,1,0];
// nm_rmsv   = [30,30,30,30]*1.3; // has to be defined if initphase = "screens"
// rotv      = [[0.,0,0,0]];

// alt     = [0.,87000.]; //45000 seems to be the limit
// nzer    = [50,50];
// nm_rmsv = [30,30]*1.5; // has to be defined if initphase = "screens"
// fit     = [1,1];
// rotv    = [[0.,0]];


//updated march 30 2026:
// alt       = [45.5,13.6,6  ,1.2,0.,-1.9,-3.3,-4.4,-8.3,-12.4,-16.5,-23.9,-29.9,-36.2]*1000; // altitude of optics, length nopt
// nm_rmsv   = [10. ,30  ,30 ,30 ,30,47  ,6.4 ,6.4 ,6.6 ,6.6  ,6.6  ,6.9   ,48   ,5];
// nzer      = [50,100,100,50,100,50,50,50,50,50,50,50,50,50]; // number of modes per optics
// fit       = [0,1,1,0,1,0,0,0,0,0,0,0,0,0];
// rotv      = [[0.,0,0,0,0,0,0,0,0,0,0,0,0,0],\
//              [180,180,180,180,180,180,180,180,180,90,0,0,0,0]]; // rotation of optics, as many line as configs
// fit       = fit*0+1;

// alt       = [45.5,13.6,6   ,1.2 ,0.,-1.9,-4 ,-12.4,-23.9,-29.9]*1000; // altitude of optics, length nopt
// nm_rmsv   = [10. ,30  ,30  ,30  ,30,47  ,9  ,11.0 ,6.9   ,48  ]*0.86;
// nzer      = [50  ,100 ,100 ,50  ,100,50 ,50 ,50   ,50    ,50]; // number of modes per optics
// fit       = [0   ,1   ,1   ,0   ,1  ,0  ,0  ,0    ,0     ,0];
// rotv      = [[0. ,0   ,0   ,0   ,0  ,0  ,0  ,0    ,0     ,0],
//              [180,180 ,180 ,180 ,180,180,180,90   ,0     ,0]];
// fit       = fit*0+1;

w = where(fit==0); if (nof(w)) nzer(w) = 2;

weight    = array(2./sqrt(nof(alt)),nof(nzer)); // mode weights (static aberrations)
// weight   *= [0.3,0.3,1,1,1,0.3];
fullfield = 30.; // full field in arcsec (on the side) - why 40 and not 30?
teldiam   = 8.0; // telescope diameter
lambda    = 550; // walevength [nm]
cobs      = 0.;  // central obstruction
size      = 64;  // side dimension of small PSF arrays
osampl    = 1;   // oversampling (1 or 2)
// gridpad   = 2.;  // padding in arcsec to ovalidate sources
gridpad   = 1.;  // padding in arcsec to ovalidate sources
ps_slope  = -2.5;
strehl_normalise = 1;
strehl_target = 0.41;
zoomfactor = 3;

// Graphics parameters
dpi_target = 160; // dpi for the "large" graphic windows
dpi_target_small = 100; // dpi for the secondary graphic windows
pltitle_height = 14;


// CHECKS
if (fit==[]) fit = long(alt*0+1);
if (rotv==[]) rotv = alt*0.;

if ((initphase=="screens")&(nm_rmsv==[])) error,"nm_rmsv undefined while initphase=\"screens\"";
if ((initphase=="coefs")&(weight==[])) error,"nm_rmsv undefined while initphase=\"screens\"";
doa = nof(alt);
if (nof(nzer)!=doa) error,"nzer and alt do not have the same dimension";
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
