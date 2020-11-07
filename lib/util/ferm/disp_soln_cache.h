// -*- C++ -*-
/*! \file
 * \brief Cache for displaced solution vectors
 */

#ifndef __disp_soln_cache_h__
#define __disp_soln_cache_h__


#include "chromabase.h"

#ifndef QDP_IS_QDPJIT_NO_NVPTX

#include <qdp_map_obj_disk_multiple.h>
#include <qdp_map_obj_memory.h>

//#include "util/ferm/distillation_soln_cache.h"

#include <vector>

namespace Chroma
{
  //----------------------------------------------------------------------------
  /*!
   * \ingroup ferm
   * @{
   */

  //! The key for displaced color vectors
  struct KeyDispSolnVector_t
  {
    bool                      use_derivP;
    std::vector<int>          displacement;    /*!< Orig plus/minus 1-based directional displacements */
    multi1d<int>              mom;
  };

  // Quark read
  void read(BinaryReader& bin, KeyDispSolnVector_t& param);

  // Quark write
  void write(BinaryWriter& bin, const KeyDispSolnVector_t& param);
  

  //---------------------------------------------------------------------
  //! Cache for distillation
  /*! 
   * \ingroup ferm 
   *
   * Holds unsmeared distillation solution vectors
   */
  class DispSolnCache
  {
  public:
    //! Default constructor
    DispSolnCache(const multi1d<LatticeColorMatrix>& u_smr,
		  const LatticeColorVectorSpinMatrix& soln_);

    //! Destructor
    virtual ~DispSolnCache() {} 
  
    //! Accessor
    const LatticeColorVectorSpinMatrix& getDispVector(bool use_derivP, const multi1d<int>& mom,
					const std::vector<int>& disp);

  protected:
    //! Displace an object
    const LatticeColorVectorSpinMatrix& displaceObject(const KeyDispSolnVector_t& key);
			
  private:
    //! Displacement length
    int displacement_length;
			
    //! Gauge field 
    const multi1d<LatticeColorMatrix>& u;
			
    // Disk cache of solutions
    const LatticeColorVectorSpinMatrix& soln;

    //! Unsmeared vectors
    QDP::MapObjectMemory<KeyDispSolnVector_t, LatticeColorVectorSpinMatrix>  disp_src_map;
  };

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif // QDP_IS_QDPJIT
#endif // HEADER GUARD
