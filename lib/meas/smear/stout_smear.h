#ifndef __stout_smear_h__
#define __stout_smear_h__


namespace Chroma { 
//! Construct the "stout-smeared" links

/*!
 * Arguments:
 *
 *  \param u		gauge field (Read)
 *  \param u_smr	"stout-smeared" gauge field (Write)
 */

void stout_smear(LatticeColorMatrix& u_smear,
		 const multi1d<LatticeColorMatrix>& u,
		 int mu, 
		 const Real& sm_fact, int j_decay);


};

#endif
