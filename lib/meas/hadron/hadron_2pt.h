// -*- C++ -*-
// $Id: hadron_2pt.h,v 1.1 2007-05-09 17:19:44 edwards Exp $
/*! \file
 *  \brief Construct hadron correlators
 */

#ifndef __hadron_2pt_h__
#define __hadron_2pt_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "io/qprop_io.h"
#include <list>

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
    //! Momentum projected correlator
    struct Mom_t
    {
      multi1d<int>       mom;    /*!< D-1 momentum of this correlator*/
      multi1d<DComplex>  corr;   /*!< Momentum projected correlator */
    };

    std::string         xml;    /*!< XML about each correlator group. Used to drive the stripper */
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
  class HadronCorrelator
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~HadronCorrelator() {}

    //! Construct the correlators
    virtual std::list<Hadron2PtCorrs_t> operator()(const multi1d<LatticeColorMatrix>& u,
						   const std::string& xml_group,
						   const std::string& id_tag) = 0;

  protected:
    //! Project onto fixed momenta
    virtual std::list<Hadron2PtCorrs_t> project(const std::list<Hadron2PtContract_t>& had_list,
						const Hadron2PtCorrParams_t& p) const;
    
    //! Convenience function to read propagator
    virtual ForwardProp_t readPropHeader(const std::string& prop_id) const;

    //! Convenience function to get t_srce from headers
    virtual multi1d<int> getTSrce(const multi1d<ForwardProp_t>& forward_headers) const;

    //! Convenience function to get decay_dir from headers
    virtual int getDecayDir(const multi1d<ForwardProp_t>& forward_headers) const;
  };

}


#endif
