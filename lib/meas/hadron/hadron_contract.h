// -*- C++ -*-
// $Id: hadron_contract.h,v 3.3 2007-06-12 16:09:37 edwards Exp $
/*! \file
 *  \brief Construct hadron correlators
 */

#ifndef __hadron_contract_h__
#define __hadron_contract_h__

#include "chromabase.h"
#include "handle.h"
#include "io/qprop_io.h"
#include <list>

namespace Chroma
{
  //! The result of hadron contractions
  /*! @ingroup hadron */
  struct HadronContractResult_t
  {
    XMLBufferWriter     xml;    /*!< XML about each corr group - used to drive the stripper */
    BinaryBufferWriter  bin;    /*!< Holds momentum projected correlators */

    XMLBufferWriter     xml_regres;  /*!< Sample XML used for regression checking */
  };
  

  //! Construct hadron correlators
  /*! @ingroup hadron
   *
   * Supports creation of hadron correlators
   */
  class HadronContract
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~HadronContract() {}

    //! Construct the correlators
    virtual std::list< Handle<HadronContractResult_t> > operator()(const multi1d<LatticeColorMatrix>& u,
								   const std::string& xml_group,
								   const std::string& id_tag) = 0;
    
  protected:
    //! Convenience function to read propagator
    virtual ForwardProp_t readForwardPropHeader(const std::string& prop_id) const;

    //! Convenience function to get t_srce from headers
    virtual multi1d<int> getTSrce(const multi1d<ForwardProp_t>& forward_headers) const;

    //! Convenience function to get decay_dir from headers
    virtual int getDecayDir(const multi1d<ForwardProp_t>& forward_headers) const;
  };

}


#endif
