// -*- C++ -*-
// $Id: extfield_aggregate_w.h,v 3.1 2006-07-20 20:06:52 edwards Exp $
/*! \file
 *  \brief External field functions
 */

#ifndef __extfield_functions_w_h__
#define __extfield_functions_w_h__

#include "actions/ferm/fermacts/extfield.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup fermacts */
  namespace ExternalFieldEnv
  {
    extern const bool registered;

  
    //! Construct zero field
    /*!
     * \ingroup fermacts
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
