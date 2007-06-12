// -*- C++ -*-
// $Id: hadron_2pt.h,v 1.4 2007-06-12 16:09:37 edwards Exp $
/*! \file
 *  \brief Construct hadron 2pt correlators
 */

#ifndef __hadron_2pt_h__
#define __hadron_2pt_h__

#include "meas/hadron/hadron_contract.h"
#include "util/ft/sftmom.h"

namespace Chroma
{
  //! The result of hadron 2pt correlators
  /*! @ingroup hadron */
  struct Hadron2PtContract_t
  {
    XMLBufferWriter    xml;   /*!< XML about each correlator group. Used to drive the stripper */
    LatticeComplex     corr;  /*!< Holds correlator after contraction, but before mom. projection */
  };


  //! The result of hadron 2pt correlators
  /*! @ingroup hadron */
  struct Hadron2PtCorrs_t
  {
    Handle<HadronContractResult_t> serialize() const;  /*!< Serialization function */

    //! Momentum projected correlator
    struct Mom_t
    {
      multi1d<int>       mom;    /*!< D-1 momentum of this correlator*/
      multi1d<DComplex>  corr;   /*!< Momentum projected correlator */
    };

    XMLBufferWriter     xml;    /*!< XML about each corr group - used to drive the stripper */
    std::list<Mom_t>    corrs;  /*!< Holds momentum projected correlators */

    XMLBufferWriter     xml_regres;  /*!< Sample XML used for regression checking */
  };
  

  //! Construct hadron 2pt correlators
  /*! @ingroup hadron
   *
   * Supports creation of hadron 2pt correlators
   *
   */
  class Hadron2PtCorr : public HadronContract
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~Hadron2PtCorr() {}

    //! Construct the correlators
    virtual std::list< Handle<HadronContractResult_t> > operator()(
      const multi1d<LatticeColorMatrix>& u,
      const std::string& xml_group,
      const std::string& id_tag) = 0;
    
  protected:
    //! Convenience function to project onto fixed momenta
    virtual std::list< Handle<HadronContractResult_t> > project(
      const std::list< Handle<Hadron2PtContract_t> >& had_list,
      const SftMomParams_t& p) const;
  };

}


#endif
