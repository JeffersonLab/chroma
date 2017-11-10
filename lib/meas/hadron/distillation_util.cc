/*! \file
 * \brief Distillution factory for producing keys * sources
 *
 * Distillution factory for producing keys * sources
 */

#include "meas/hadron/distillation_util.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  namespace Distillation
  {
    //----------------------------------------------------------------------------
    //! Get sink key
    KeyPropDistillation_t getSnkKey(const KeyPropElementalOperator_t& peram_key, int colorvec_src)
    {
      KeyPropDistillation_t snk_key;

      snk_key.t_source     = peram_key.t_source;
      snk_key.t_slice      = peram_key.t_slice;
      snk_key.colorvec_src = colorvec_src;
      snk_key.spin_src     = peram_key.spin_src;
      snk_key.spin_snk     = peram_key.spin_snk;
      snk_key.mass         = peram_key.mass_label;

      return snk_key;
    }


    //----------------------------------------------------------------------------
    //! Prepare a distilluted source
    LatticeColorVector getSrc(TimeSliceIOCache& eigen_source_cache, int t_source, int colorvec_src)
    {
      LatticeColorVectorF vec_srce = zero;

      vec_srce = eigen_source_cache.getVec(t_source, colorvec_src);

      return vec_srce;
    }

    //----------------------------------------------------------------------------
    //! Get active time-slices
    std::vector<bool> getActiveTSlices(int t_source, int Nt_forward, int Nt_backward)
    {
      // Initialize the active time slices
      const int Lt = Layout::getTimeExtent();

      std::vector<bool> active_t_slices(Lt);
      for(int t=0; t < Lt; ++t)
      {
	active_t_slices[t] = false;
      }

      // Forward
      for(int dt=0; dt < Nt_forward; ++dt)
      {
	int t = t_source + dt;
	active_t_slices[t % Lt] = true;
      }

      // Backward
      for(int dt=0; dt < Nt_backward; ++dt)
      {
	int t = t_source - dt;
	while (t < 0) {t += Lt;} 

	active_t_slices[t % Lt] = true;
      }

      return active_t_slices;
    }


    //----------------------------------------------------------------------------
    //! Get source keys
    KeyPropDistillation_t getSrcKey(int t_source, int colorvec_src)
    {
      KeyPropDistillation_t key;

      key.t_source     = t_source;
      key.t_slice      = t_source;
      key.colorvec_src = colorvec_src;
      key.spin_src     = -1;
      key.spin_snk     = -1;
      // key.mass         = "";   // Using an empty key

      return key;
    }

	
    //----------------------------------------------------------------------------
    //! Get sink keys
    std::list<KeyPropDistillation_t> getSnkKeys(int t_source, int colorvec_src, int Nt_forward, int Nt_backward, const std::string mass)
    {
      std::list<KeyPropDistillation_t> keys;

      std::vector<bool> active_t_slices = getActiveTSlices(t_source, Nt_forward, Nt_backward);
	
      const int Lt = Layout::getTimeExtent();

      for(int spin_source=0; spin_source < Ns; ++spin_source)
      {
	for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	{
	  for(int t=0; t < Lt; ++t)
	  {
	    if (! active_t_slices[t]) {continue;}

	    KeyPropDistillation_t key;

	    key.t_source     = t_source;
	    key.t_slice      = t;
	    key.colorvec_src = colorvec_src;
	    key.spin_src     = spin_source;
	    key.spin_snk     = spin_sink;
	    key.mass         = mass;

	    //            QDPIO::cout << key << std::flush;

	    keys.push_back(key);
	  } // for t
	} // for spin_sink
      } // for spin_source

      return keys;
    }
      
    //----------------------------------------------------------------------------
    //! Get perambulator keys
    std::list<KeyPropElementalOperator_t> getPeramKeys(int t_source, int Nt_forward, int Nt_backward, const std::string mass)
    {
      std::list<KeyPropElementalOperator_t> keys;

      std::vector<bool> active_t_slices = getActiveTSlices(t_source, Nt_forward, Nt_backward);
	
      const int Lt = Layout::getTimeExtent();

      for(int spin_source=0; spin_source < Ns; ++spin_source)
      {
	for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	{
	  for(int t=0; t < Lt; ++t)
	  {
	    if (! active_t_slices[t]) {continue;}

	    KeyPropElementalOperator_t key;

	    key.t_slice      = t;
	    key.t_source     = t_source;
	    key.spin_src     = spin_source;
	    key.spin_snk     = spin_sink;
	    key.mass_label   = mass;

	    //            QDPIO::cout << key << std::flush;

	    keys.push_back(key);
	  } // for t
	} // for spin_sink
      } // for spin_source

      return keys;
    }

	
    //----------------------------------------------------------------------------
    //! Get perambulator key time slices
    std::list<int> getTslices(int t_source, int Nt_forward, int Nt_backward)
    {
      std::list<int> keys;

      std::vector<bool> active_t_slices = getActiveTSlices(t_source, Nt_forward, Nt_backward);
	
      const int Lt = Layout::getTimeExtent();

      for(int t=0; t < Lt; ++t)
      {
	if (active_t_slices[t])
	{
	  keys.push_back(t);
	}
      } // for t

      return keys;
    }

	
  } // namespace Connected

} // namespace Chroma
