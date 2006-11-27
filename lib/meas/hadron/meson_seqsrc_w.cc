// $Id: meson_seqsrc_w.cc,v 3.1 2006-11-27 04:33:35 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "meas/hadron/meson_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"

namespace Chroma 
{

  //! Anonymous namespace
  /*! @ingroup hadron */
  namespace
  {
    //! Construct hadron sequential sources
    /*! @ingroup hadron */
    template<typename T> 
    T   mesSeqSourceProject(const T& seq_src_tmp,
			    int t_sink, int j_decay)
    {
      START_CODE();

      if (j_decay < 0 || j_decay >= Nd)
      {
	QDPIO::cerr << __func__ << ": j_decay out of bounds" << endl;
	QDP_abort(1);
      }

      /*
       * Now mask out all but sink time slice
       */
      T seq_src_prop = where(Layout::latticeCoordinate(j_decay) == t_sink,
			     seq_src_tmp,
			     T(zero));
        
      END_CODE();

      return seq_src_prop;
    }
  }




}  // end namespace Chroma


  
