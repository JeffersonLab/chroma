// su3_hb_sweep.h, 2004/11/16 velytsky
// one su(3) sweep
#ifndef __su3_hb_sweep__
#define __su3_hb_sweep__
/* *************************************
 * u		link field
 * HBParams	HB parameter container
 * ************************************/
void su3_hb_sweep(multi1d<LatticeColorMatrix>& u,
                HBParams&);
#endif
