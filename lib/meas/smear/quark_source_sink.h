// -*- C++ -*-
// $Id: quark_source_sink.h,v 3.1 2006-06-12 02:13:47 edwards Exp $

/*! @file
 * @brief Quark source or sink smearing
 */

#ifndef __quark_source_sink_h__
#define __quark_source_sink_h__

#include "chromabase.h"
#include "handle.h"
#include "meas/smear/link_smearing_factory.h"
#include "io/xml_group_reader.h"

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
     * \param link_smearing       group holding XML of link smearing ( Read )
     */
    virtual void create(multi1d<LatticeColorMatrix>& u,
			const GroupXML_t& link_smearing)
    {
      //
      // Smear the gauge field if needed
      //
      std::istringstream  xml_l(link_smearing.xml);
      XMLReader  linktop(xml_l);
      const string link_path = "/LinkSmearing";
	
      Handle< LinkSmearing >
	linkSmearing(TheLinkSmearingFactory::Instance().createObject(link_smearing.id,
								     linktop,
								     link_path));
      (*linkSmearing)(u);
    }


  };

}


#endif
