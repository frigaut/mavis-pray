// config is system config, including defocus values and rotation for some optics
// One could conceptually have different defoc for the various rotation, but
// in here we assume and use the same defocs (e.g. [-1,0,1]) for the various rotation
// configs (e.g. [0,0,0] and [180,0,0]). Rotations config are defined below.

struct pray_struct
{
  pointer images;
  pointer mircube;
  pointer ftobject;
  pointer variance;
  pointer norm;
  pointer nzer;
  pointer alt;
  pointer ipupil;
  pointer def;
  pointer def_pup;
  pointer dmgsxposcub;
  pointer dmgsyposcub;
  pointer deltafoc;
  pointer xpos;
  pointer ypos;
  pointer focus;
  pointer ampli_pup;
  pointer ampli_foc;
  float   shiftFoc;
  float   teldiam;
  float   cobs;
  long    pupd;
  long    size;
  long    _n;
  long    _n1;
  long    _n2;
};

struct sim_struct { long verbose; }

struct config_struct {
  float foc;
  long roti;
}
