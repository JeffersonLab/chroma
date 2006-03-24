// $Id: szin_write_obj_funcmap.cc,v 2.1 2006-03-24 22:16:40 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/szin_write_obj_funcmap.h"
#include "io/szin_io.h"
#include "io/writeszin.h"
#include "io/writeszinqprop_w.h"

namespace Chroma
{
 
  //! IO function map environment
  /*! \ingroup inlineio */
  namespace SZINWriteObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //! Write a propagator
      void SZINWriteLatProp(const string& buffer_id,
			    const string& file, 
			    int j_decay, int t_start, int t_end)
      {
	LatticePropagator obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
    
	Real Kappa = 0.0;
	writeSzinQprop(obj, file, Kappa);
      }

      //! Write a gauge field in floating precision
      void SZINWriteArrayLatColMat(const string& buffer_id,
				   const string& file, 
				   int j_decay, int t_start, int t_end)
      {
	SzinGauge_t szin_out;   // ignoring XML in named object

	multi1d<LatticeColorMatrix> obj = 
	  TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);

	writeSzinTrunc(szin_out, obj, j_decay,
		       t_start, t_end,
		       file);
      }


    }  // end namespace


    bool registerAll(void) 
    {
      bool success = true;
      success &= TheSZINWriteObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								     SZINWriteLatProp);
      success &= TheSZINWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								     SZINWriteArrayLatColMat);
      return success;
    }

    bool registered = registerAll();
  }

}
