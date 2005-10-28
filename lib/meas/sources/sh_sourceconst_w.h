// -*- C++ -*-
// $Id: sh_sourceconst_w.h,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Shell source construction
 */

#ifndef __sh_sourceconst_w_h__
#define __sh_sourceconst_w_h__

#include "meas/sources/source_construction.h"
#include "meas/sources/shsrc_params.h"

namespace Chroma
{

  //! Name and registration
  namespace PropShellSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;
    //! Name to be used
  }
  

  //! Shell source construction
  /*! @ingroup sources
   *
   * Create a point propagator source
   */
  class PropShellSourceConst : public SourceConstruction<LatticePropagator>
  {
  public:
    //! Full constructor
    PropShellSourceConst(const ShellSourceConstParams& p) : params(p) {}

    //! Construct the source
    LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    PropShellSourceConst() {}

  private:
    ShellSourceConstParams  params;   /*!< source params */
  };

}  // end namespace Chroma


#endif
