/*! \file
 * \brief Cache for distillation - holds solution vectors
 */

#include "util/ferm/distillation_soln_cache.h"

namespace Chroma 
{
  // Utility functions
  namespace
  {
    //! Error output
    StandardOutputStream& operator<<(StandardOutputStream& os, const std::list<int>& d)
    {
      bool firstP = true;

      for(std::list<int>::const_iterator qq = d.begin(); qq != d.end(); ++qq)
      {
	if (firstP)
	{
	  os << *qq;
	  firstP = false;
	}
	else 
	{
	  os << " " << *qq;
	}
      }

      return os;
    }


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

    bool isSnk(const std::string& prop_type)
    {
      return (prop_type == "SNK") ? true : false;
    }

  }


  //---------------------------------------------------------------------
  // Constructor
  DistillationSolnCache::DistillationSolnCache(int t_slice_,
					       int t_origin_,
					       QDP::MapObjectDiskMultiple<KeyPropDistillation_t, LatticeColorVector>& soln_disk_,
					       const std::map<char,std::string>& flav_map_,
					       int num_vecs_) :
    t_slice(t_slice_), t_origin(t_origin_), soln_disk(soln_disk_), flavor_to_mass(flav_map_), num_vecs(num_vecs_)
  {
  }

  //---------------------------------------------------------------------
  //! Get vectors
  const LatticeColorVector& DistillationSolnCache::getVector(const KeyHadronNode_t::Quark_t& quark_line, 
							     int colorvec_src, int spin_snk, int spin_src)
  {
    return getUnsmeared(createKey(quark_line, colorvec_src, spin_snk, spin_src));
  }


  //---------------------------------------------------------------------
  // Return a mass label
  const std::string& DistillationSolnCache::getMass(char flav) const 
  {
    std::map<char,std::string>::const_iterator f = flavor_to_mass.find(flav);
    if (f == flavor_to_mass.end())
    {
      QDPIO::cerr << __func__ << ": flavor not in map: flavor= " << flav << "\n";
      QDP_abort(1);
    }
    
    return f->second;
  }


  //---------------------------------------------------------------------
  // Return a key
  KeyPropDistillation_t DistillationSolnCache::createKey(const KeyHadronNode_t::Quark_t& quark, 
							 int colorvec_src, int spin_snk, int spin_src) const
  {
    if (quark.t_slice != t_slice)
    {
      QDPIO::cerr << __func__ << ": trying to read from a t_slice different than the active one: t_slice= " 
		  << t_slice << "\n";
      QDP_abort(1);
    }

    KeyPropDistillation_t key;

    if (isSnk(quark.prop_type))
    {
      key.t_source      = (quark.t_source + t_origin + Layout::getTimeExtent()) % Layout::getTimeExtent();
      key.t_slice       = (quark.t_slice + t_origin + Layout::getTimeExtent()) % Layout::getTimeExtent();
      key.mass          = getMass(quark.flavor);
      key.colorvec_src  = colorvec_src;
      key.spin_snk      = spin_snk;
      key.spin_src      = spin_src;
    }
    else
    {
      key.t_slice       = (quark.t_slice + t_origin + Layout::getTimeExtent()) % Layout::getTimeExtent();
      key.t_source      = key.t_slice;
      key.mass          = "";
      key.colorvec_src  = colorvec_src;
      key.spin_snk      = spin_snk;
      key.spin_src      = spin_src;
    }

    return key;
  }

  //---------------------------------------------------------------------
  //! Get unsmeared vectors
  LatticeColorVector& DistillationSolnCache::getUnsmeared(const KeyPropDistillation_t& key)
  {
    if (! unsmeared_vectors.exist(key))
    {
      // Snarf
      LatticeColorVector val;

      if (soln_disk.get(key, val) != 0)
      {
	QDPIO::cerr << __func__ << ": missing key= " << key << std::endl;

	std::vector<KeyPropDistillation_t> keys;
	soln_disk.keys(keys);

	QDPIO::cout << "Dumping keys\n";
	for(std::vector<KeyPropDistillation_t>::const_iterator k=keys.begin(); k != keys.end(); ++k)
	  QDPIO::cout << *k;
	
	QDP_abort(1);
      }

      // Succeeded. Insert into local cache.
      unsmeared_vectors.insert(key, val);
    }

    return unsmeared_vectors[key];
  }

}  // end namespace Chroma
