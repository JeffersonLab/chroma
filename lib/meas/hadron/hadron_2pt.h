// -*- C++ -*-
// $Id: hadron_2pt.h,v 1.3 2007-06-10 14:49:06 edwards Exp $
/*! \file
 *  \brief Construct hadron 2pt correlators
 */

#ifndef __hadron_2pt_h__
#define __hadron_2pt_h__

#include "meas/hadron/hadron_contract.h"

namespace Chroma
{
  //! The result of hadron 2pt correlators
  /*! @ingroup hadron */
  struct Hadron2PtContract_t
  {
    std::string        xml;   /*!< XML about each correlator group. Used to drive the stripper */
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

    std::string         xml;    /*!< XML about each corr group - used to drive the stripper */
    std::list<Mom_t>    corrs;  /*!< Holds momentum projected correlators */
  };
  

  //! The result of hadron 2pt correlators
  /*! @ingroup hadron */
  struct Hadron2PtCorrParams_t
  {
    int           mom2_max;           /*!< (mom - mom_origin)^2 <= mom2_max */
    multi1d<int>  mom_origin;         /*!< Origin for the momentum */
    bool          avg_equiv_mom;      /*!< average over equivalent momenta */
    multi1d<int>  t_srce;             /*<! Origin of the prop. Here used to offset momentum phases */
    int           decay_dir;          /*!< Decay direction */
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
      const std::list< Hadron2PtContract_t >& had_list,
      const Hadron2PtCorrParams_t& p) const;
  };

}


#endif
