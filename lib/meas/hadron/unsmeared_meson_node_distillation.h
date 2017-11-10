// -*- C++ -*-
/*! \file
 * \brief Meson nodes
 */

#ifndef __unsmeared_meson_node_distillation_h__
#define __unsmeared_meson_node_distillation_h__

#include "chromabase.h"
#include "meas/hadron/unsmeared_distillation_node.h"
#include "util/ferm/hadron_node.h"
#include "util/ferm/single_hadron_coeffs.h"
#include "util/ferm/spin_rep.h"
#include <hadron_dist_obj.h>

namespace Chroma
{
  //----------------------------------------------------------------------------------
  //! Holds/transmogrifies hadron operator coeffs
  class UnsmearedMesonNodeDistillation : public UnsmearedDistillationNode
  {
  public:
    //! Constructor
    UnsmearedMesonNodeDistillation(const HadronVertex_t& key_, 
				   const MapSingleHadronQuarkFlavorColorSpin_t& coeffs_,
				   const std::map<char,std::string>& flav_map_,
				   int num_vecs);

    //! Virtual destructor
    virtual ~UnsmearedMesonNodeDistillation() {}

    //! Contract to form a hadron node
    virtual void constructNode(Hadron::HadronDistOperatorRep& node, 
			       const KeyHadronNode_t& hadron_node, 
			       DistillationSolnCache& soln_cache,
			       DispSolnCache& disp_soln_cache) const;

  private:
    HadronVertex_t                        key;
    MapSingleHadronQuarkFlavorColorSpin_t coeffs;
    std::map<char,std::string>            flav_map;
    int                                   num_vecs;
    bool                                  use_derivP;
    
  private:
    // Conversion to dp
    std::vector<MatrixSpinRep_t>          diracToDrMatPlus;
    std::vector<MatrixSpinRep_t>          diracToDrMatMinus;
  };

} // namespace Chroma

#endif
