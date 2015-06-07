/*! \file
 *  \brief Write object function std::map
 */

#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/szin_write_obj_funcmap.h"
#include "io/szin_io.h"
#include "io/writeszin.h"
#include "io/writeszinqprop_w.h"

namespace Chroma
{
 
  //! IO function std::map environment
  /*! \ingroup inlineio */
  namespace SZINWriteObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //! Write a propagator
      void SZINWriteLatProp(const std::string& buffer_id,
			    const std::string& file, 
			    int j_decay, int t_start, int t_end)
      {
	LatticePropagator obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
    
	Real Kappa = 0.0;
	writeSzinQprop(obj, file, Kappa);
      }

      //! Write a gauge field in floating precision
      void SZINWriteArrayLatColMat(const std::string& buffer_id,
				   const std::string& file, 
				   int j_decay, int t_start, int t_end)
      {
	SzinGauge_t szin_out;   // ignoring XML in named object

	multi1d<LatticeColorMatrix> obj = 
	  TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);

	writeSzinTrunc(szin_out, obj, j_decay,
		       t_start, t_end,
		       file);
      }

      //! Local registration flag
      bool registered = false;

    }  // end namespace



    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheSZINWriteObjFuncMap::Instance().registerFunction(std::string("LatticePropagator"), 
								       SZINWriteLatProp);
	success &= TheSZINWriteObjFuncMap::Instance().registerFunction(std::string("Multi1dLatticeColorMatrix"), 
								       SZINWriteArrayLatColMat);

	registered = true;
      }
      return success;
    }
  }

}
