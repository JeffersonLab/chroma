// -*- C++ -*-
/*! \file
 * \brief Cache for displaced solution vectors
 */

#ifndef __disp_soln_cache_h__
#define __disp_soln_cache_h__

#include "chromabase.h"

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
		  const LatticeFermion& soln_);

    //! Destructor
    virtual ~DispSolnCache() {} 
  
    //! Accessor
    const LatticeFermion& getDispVector(bool use_derivP, const multi1d<int>& mom,
					const std::vector<int>& disp);

  protected:
    //! Displace an object
    const LatticeFermion& displaceObject(const KeyDispSolnVector_t& key);
			
  private:
    //! Displacement length
    int displacement_length;
			
    //! Gauge field 
    const multi1d<LatticeColorMatrix>& u;
			
    // Disk cache of solutions
    const LatticeFermion& soln;

    //! Unsmeared vectors
    QDP::MapObjectMemory<KeyDispSolnVector_t, LatticeFermion>  disp_src_map;
  };

  /*! @} */  // end of group ferm

} // namespace Chroma

#endif
