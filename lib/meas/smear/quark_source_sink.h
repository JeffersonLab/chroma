// -*- C++ -*-
// $Id: quark_source_sink.h,v 2.3 2005-11-22 22:00:27 edwards Exp $

/*! @file
 * @brief Quark source or sink smearing
 */

#ifndef __quark_source_sink_h__
#define __quark_source_sink_h__

#include "chromabase.h"
#include "meas/smear/link_smearing_factory.h"

namespace Chroma
{
  //! Base class for quark source and sink smearing
  /*! @ingroup smear
   *
   * Supports creation and application of smearing (with link smearing)
   * on quarks, and potentially displacements. Basically the construction
   * of a "source" or "sink" state on a pre-existing quark object
   */
  template<typename T>
  class QuarkSourceSink
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~QuarkSourceSink() {}

    //! Smear the quark
    /*!
     * \param obj      Object to source or sink smear ( Modify )
     */
    virtual void operator()(T& obj) const = 0;


  protected:
    //! Potentially smear the gauge field
    /*!
     * \param u                   Gauge field to smear ( Modify )
     * \param link_smearing       XML of link smearing ( Read )
     * \param link_smearing_type  link smearing type ( Read )
     */
    virtual void create(multi1d<LatticeColorMatrix>& u,
			std::string link_smearing,
			std::string link_smearing_type)
    {
      linkSmear(u, std::string("/LinkSmearing"), link_smearing, link_smearing_type);
    }


  };

}


#endif
