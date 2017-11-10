// -*- C++ -*-
/*! \file
 * \brief Generic container for producing distilled objects and hadron nodes
 */

#ifndef __unsmeared_distillation_node_h__
#define __unsmeared_distillation_node_h__

#include <hadron_dist_obj.h>
#include "util/ferm/hadron_node.h"
#include "util/ferm/distillation_soln_cache.h"
#include "util/ferm/disp_soln_cache.h"

namespace Chroma
{
  //----------------------------------------------------------------------------
  //! Generic container for colorvector and spin matrix element 
  class UnsmearedDistillationNode
  {
  public:
    //! Destructor
    virtual ~UnsmearedDistillationNode() {}

    //! Contract to form a hadron node
    virtual void constructNode(Hadron::HadronDistOperatorRep& node, 
			       const KeyHadronNode_t& hadron_node, 
			       DistillationSolnCache& soln_cache,
			       DispSolnCache& disp_soln_cache) const = 0;
  };


} // namespace Chroma

#endif
