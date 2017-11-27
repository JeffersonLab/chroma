/*! \file
 * \brief Cache for distillation - holds solution vectors
 */

#include "util/ferm/disp_soln_cache.h"
#include "meas/smear/displace.h"

namespace Chroma 
{
  // Utility functions
  namespace
  {
    //! Error output
    StandardOutputStream& operator<<(StandardOutputStream& os, const multi1d<int>& d)
    {
      if (d.size() > 0)
      {
	os << d[0];

	for(int i=1; i < d.size(); ++i)
	  os << " " << d[i];
      }

      return os;
    }
  }


  //----------------------------------------------------------------------------
  // Quark read
  void read(BinaryReader& bin, KeyDispSolnVector_t& param)
  {
    read(bin, param.use_derivP);
    read(bin, param.displacement);
    read(bin, param.mom);
  }

  // Quark write
  void write(BinaryWriter& bin, const KeyDispSolnVector_t& param)
  {
    write(bin, param.use_derivP);
    write(bin, param.displacement);
    write(bin, param.mom);
  }


  //----------------------------------------------------------------------------
  // Constructor from smeared map 
  DispSolnCache::DispSolnCache(const multi1d<LatticeColorMatrix>& u_smr,
			       const LatticeFermion& soln_)
    : displacement_length(1), u(u_smr), soln(soln_)
  {
  }


  //! Accessor
  const LatticeFermion&
  DispSolnCache::getDispVector(bool use_derivP, const multi1d<int>& mom,
			       const std::vector<int>& disp)
  {
    KeyDispSolnVector_t key;
    key.use_derivP     = use_derivP;
    key.mom            = mom;
    key.displacement   = disp;

    return displaceObject(key);
  }


  //! Accessor
  const LatticeFermion&
  DispSolnCache::displaceObject(const KeyDispSolnVector_t& key)
  {
    // If no entry, then create a displaced version of the quark
    if (! disp_src_map.exist(key))
    {
      // Only need to do more if there are displacements
      if (key.displacement.size() == 0)
      {
	// Pull out the soln vector
	//	disp_src_map.insert(key, soln);
	return soln;
      }
      else
      {
	// Have at least one displacement. Pull that one off the end of the list
	KeyDispSolnVector_t prev_key = key;

	int d = key.displacement.back();
	prev_key.displacement.pop_back();

	// Recursively get a reference to the object to be shifted 
	const LatticeFermion& disp_q = this->displaceObject(prev_key);

	// Displace or deriv the old vector
	if (d > 0)
	{
	  int disp_dir = d - 1;
	  int disp_len = displacement_length;
	  if (key.use_derivP)
	    disp_src_map.insert(key, leftRightNabla(disp_q, u, disp_dir, disp_len, key.mom[disp_dir]));
	  else
	    disp_src_map.insert(key, displace(u, disp_q, disp_len, disp_dir));
	}
	else if (d < 0)
	{
	  if (key.use_derivP)
	  {
	    QDPIO::cerr << __func__ << ": do not support (rather do not want to support) negative displacements for rightNabla\n";
	    QDP_abort(1);
	  }

	  int disp_dir = -d - 1;
	  int disp_len = -displacement_length;
	  disp_src_map.insert(key, displace(u, disp_q, disp_len, disp_dir));
	}
      }
    } // if find in map

    // The key now must exist in the map, so return the vector
    return disp_src_map[key];
  }

  /*! @} */  // end of group smear

} // namespace Chroma
