// $Id: baryon_2pt_w.cc,v 1.1 2007-05-09 17:19:44 edwards Exp $
/*! \file
 *  \brief Construct baryon 2pt correlators.
 */

#include "meas/hadron/baryon_2pt_w.h"
#include "meas/hadron/hadron_2pt_factory_w.h"

namespace Chroma 
{

  //! Anonymous namespace
  /*! @ingroup hadron */
  namespace
  {
    //! Get fermion bc
    multi1d<int>
    barSeqSourceGetBC(const multi1d<ForwardProp_t>& forward_headers)
    {
      multi1d<int> bc = getFermActBoundary(forward_headers[0].prop_header.fermact);

      for(int loop=1; loop < forward_headers.size(); ++loop)
      {
	multi1d<int> bc_b = getFermActBoundary(forward_headers[loop].prop_header.fermact);
    
	// Bummer, I wish qdp++ had a multi1d.operator!=()
	bool same = true;
	for(int i=0; i < bc.size(); ++i)
	{
	  if (bc_b[i] != bc[i]) 
	    same = false;
	}
      
	if (! same)
	{
	  QDPIO::cerr << __func__ << ": the bc in the forward props are not all equal"
		      << endl;
	  QDP_abort(1);
	}
      }

      return bc;
    }

  }



  // Default versions
  void
  BaryonSeqSourceBase::setBC(const multi1d<ForwardProp_t>& forward_headers)
  {
    getBC() = barSeqSourceGetBC(forward_headers);
  }

}  // end namespace Chroma
