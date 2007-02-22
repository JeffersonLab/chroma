//  $Id: sun_proj.cc,v 3.2 2007-02-22 21:11:50 bjoo Exp $
/*! \file
 *  \ingroup gauge
 *  \author Subsetting added by A. Hart
 *  \param[in] w            complex Nc x Nc matrix
 *  \param[out] v           the projected SU(Nc) Matrix
 *  \param[in] BlkAccu      accuracy in SU(Nc) projection
 *  \param[in] BlkMax       max number of iterations in SU(Nc) projection
 *  \param[in] mstag        an (un)ordered subset of lattice sites
 *  \brief Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 *
 *  Project a complex Nc x Nc matrix W onto SU(Nc) by maximizing Tr(VW)
 */

#include "chromabase.h"

#include "util/gauge/sun_proj.h"
#include "util/gauge/su3proj.h"
#include "util/gauge/reunit.h" 


namespace Chroma 
{ 

  template<typename S>
  inline
  void sun_proj_t(const LatticeColorMatrix& w, 
		  LatticeColorMatrix& v,
		  const Real& BlkAccu, 
		  int BlkMax,
		  const S& mstag)
  {
    START_CODE();

    Double new_tr;

    /*
     * Project back to SU(3) by maximizing tr(v w).
     * This is done by looping proj_iter times over the 3 SU(2) subgroups.
     */

    /* I need to get the number of sites in the sublattice. As far as
     * I know, mstag does not contain information about either the number
     * of sites it contains, nor of the number of subsets.
     */

    LatticeInt count;
    count = 1;
    Int numSites = sum(count,mstag);
    Double norm = Double(1)/(Nc*numSites);

    /* The initial trace */
    Double old_tr = sum(real(trace(v * w)),mstag) * norm;

    int iter = 0;
    int wrswitch = 0;			// Write out iterations?
//  Double conver = 1.0;
    Real conver = 1.0;

    while ( toBool(conver > BlkAccu)  &&  iter < BlkMax )
    {
      ++iter;

      // Loop over SU(2) subgroup index
      for(int su2_index = 0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	su3proj(v, w, su2_index, mstag);

      // Reunitarize: this is the slow bit of the code...
      reunit(v,mstag);

      // Calculate the trace
      new_tr = sum(real(trace(v * w)), mstag) * norm;

      if( wrswitch == 1 )
      {
	QDPIO::cout << "iter =     " << iter << endl;
	QDPIO::cout << "  old_tr = " << old_tr << endl;
	QDPIO::cout << "  new_tr = " << new_tr << endl;
      }

      // Normalized convergence criterion:
      conver = fabs((new_tr - old_tr) / old_tr);
      old_tr = new_tr;
    }

#if 0
    if ( wrswitch == 1 )
    {
//    push(nml,"Final_sun_proj");
//    write(nml, "iter", iter);
//    write(nml, "new_tr", new_tr);
//    pop(nml);
      QDPIO::cout << "iter = " << iter << endl;
      QDPIO::cout << "new_tr = " << new_tr << endl;
    }
#endif

    END_CODE();
  }

  void sun_proj(const LatticeColorMatrix& w, 
		LatticeColorMatrix& v,
		const Real& BlkAccu, 
		int BlkMax)
  {
    sun_proj_t(w, v, BlkAccu, BlkMax, all);
  }

  void sun_proj(const LatticeColorMatrix& w, 
		LatticeColorMatrix& v,
		const Real& BlkAccu, 
		int BlkMax,
		const Subset& mstag)
  {
    sun_proj_t(w, v, BlkAccu, BlkMax, mstag);
  }

}
