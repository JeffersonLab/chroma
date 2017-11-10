// -*- C++ -*-
/*! \file
 * \brief Distillution factory for producing keys * sources
 *
 * Distillution factory for producing keys * sources
 */

#ifndef __distillation_util_w_h__
#define __distillation_util_w_h__

#include "chromabase.h"
#include "util/ferm/timeslice_io_cache.h"
#include "util/ferm/key_prop_distillation.h"
#include "util/ferm/key_prop_matelem.h"

namespace Chroma 
{ 
  //----------------------------------------------------------------------------
  namespace Distillation
  {
    //----------------------------------------------------------------------------
    //! Get sink key
    KeyPropDistillation_t getSnkKey(const KeyPropElementalOperator_t& peram_key, int colorvec_src);

    //----------------------------------------------------------------------------
    //! Prepare a distilluted source
    LatticeColorVector getSrc(TimeSliceIOCache& eigen_source_cache, int t_source, int colorvec_src);

    //----------------------------------------------------------------------------
    //! Get active time-slices
    std::vector<bool> getActiveTSlices(int t_source, int Nt_forward, int Nt_backward);

    //----------------------------------------------------------------------------
    //! Get source key
    KeyPropDistillation_t getSrcKey(int t_source, int colorvec_src);

    //----------------------------------------------------------------------------
    //! Get sink keys
    std::list<KeyPropDistillation_t> getSnkKeys(int t_source, int colorvec_src, int Nt_forward, int Nt_backward, const std::string mass);

    //----------------------------------------------------------------------------
    //! Get perambulator keys
    std::list<KeyPropElementalOperator_t> getPeramKeys(int t_source, int Nt_forward, int Nt_backward, const std::string mass);

    //----------------------------------------------------------------------------
    //! Get perambulator key time slices
    std::list<int> getTslices(int t_source, int Nt_forward, int Nt_backward);

  }

}

#endif
