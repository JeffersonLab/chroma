// su2_hb_update.h, 2004/11/16 velytsky
// one su(2) hb update kennedy-pendleton or creutz
#ifndef __su2_hb_update__
#define __su2_hb_update__

namespace Chroma {

/* ********************************
 * u_mu 		link field
 * u_mu_staple          staple attached to the link u_mu
 * BetaMC		beta value
 * su2_index		index of su(2) subgroup of su(n)
 * sub                  subset for updating
 * NmaxHB		number of maximum hb tries
 * *******************************/
void su2_hb_update(LatticeColorMatrix& u_mu, 
			const LatticeColorMatrix& u_mu_staple, 
			Double BetaMC, 
			const int su2_index,
			const Subset& sub, 
			const int NmaxHB);

}  // end namespace Chroma

#endif
