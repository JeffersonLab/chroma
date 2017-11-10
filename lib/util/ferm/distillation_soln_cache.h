// -*- C++ -*-
/*! \file
 * \brief Cache for distillation - holds solution vectors
 */

#ifndef __distillation_soln_cache_h__
#define __distillation_soln_cache_h__

#include "chromabase.h"

#include <qdp_map_obj_disk_multiple.h>
#include <qdp_map_obj_memory.h>

#include "util/ferm/hadron_node.h"
#include "util/ferm/key_prop_distillation.h"
#include "util/ferm/flavor_to_mass.h"

#include <vector>

namespace Chroma
{
  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */

  //---------------------------------------------------------------------
  //! Cache for distillation
  /*! 
   * \ingroup ferm 
   *
   * Holds unsmeared distillation solution vectors
   */
  class DistillationSolnCache
  {
  public:
    //! Default constructor
    DistillationSolnCache(int t_slice_,
			  int t_origin_,
			  QDP::MapObjectDiskMultiple<KeyPropDistillation_t, LatticeColorVector>& soln_disk_,
			  const std::map<char,std::string>& flav_map_,
			  int num_vecs_);

    //! Destructor
    virtual ~DistillationSolnCache() {} 
  
    //! Get vectors
    const LatticeColorVector& getVector(const KeyHadronNode_t::Quark_t& quark_line, 
					int colorvec_src, int spin_snk, int spin_src);

  protected:
    //! Return a mass label
    const std::string& getMass(char f) const;

    //! Return a key
    KeyPropDistillation_t createKey(const KeyHadronNode_t::Quark_t& quark, 
				    int colorvec_src, int spin_snk, int spin_src) const;

    //! Get unsmeared vectors
    LatticeColorVector& getUnsmeared(const KeyPropDistillation_t& key);

  private:
    // The active time slice
    int t_slice;

    // The origin for this config
    int t_origin;

    // Disk cache of solutions
    QDP::MapObjectDiskMultiple<KeyPropDistillation_t, LatticeColorVector>&   soln_disk;

    // The map from quark flavors to mass labels
    std::map<char,std::string>  flavor_to_mass;

    // Number of vectors to use
    int num_vecs;

    //! Unsmeared vectors
    QDP::MapObjectMemory<KeyPropDistillation_t, LatticeColorVector>  unsmeared_vectors;
  };

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
