// -*- C++ -*-
// $Id: extfield_aggregate_w.h,v 1.1 2006-09-19 17:53:36 edwards Exp $
/*! \file
 *  \brief External field functions
 */

#ifndef __extfield_functions_w_h__
#define __extfield_functions_w_h__

#include "actions/ferm/fermstates/extfield.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup fermstates */
  namespace ExternalFieldEnv
  {
    extern const bool registered;

  
    //! Construct zero field
    /*!
     * \ingroup fermstates
     */
    class ZeroExternalField : public ExternalField
    {
    public:
      //! Full constructor
      ZeroExternalField() {}

      //! Return the field
      LatticeComplex operator()(int dir) const;
    };

  }  // end namespace


#if 0
  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DerivQuarkDisplacementEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DerivQuarkDisplacementEnv::Params& param);
#endif

}  // end namespace Chroma

#endif
