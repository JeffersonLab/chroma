// u_staple.h, 2004/11/10 velytsky
// calculate the sum of all staples for a particular u_mu link
#ifndef __u_staple__
#define __u_staple__
/* ****************************************
 * u                    link field
 * mu                   direction of the link
 * u_mu_staple          staple attached to the link u
 * sub                  subset for updating
 * HBParams             container of HB parameters
 * ***************************************/
void u_staple(const multi1d<LatticeColorMatrix>& u,
                const int mu,
                LatticeColorMatrix& u_mu_staple,
                const OrderedSubset& sub,
                HBParams& );
#endif
