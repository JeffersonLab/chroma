// -*- C++ -*-
// $Id: mesonseqsrc_w.h,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#ifndef __mesonseqsrc_w_h__
#define __mesonseqsrc_w_h__

#include "meas/hadron/hadron_seqsource.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleMesonSeqSourceEnv
  {
    extern const bool registered;

  
    //! Simple meson sequential source parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      multi1d<int>     sink_mom;        /*!< sink momentum */
      int              t_sink;          /*!< time slice of sink */
      int              j_decay;         /*!< Decay direction */
    };


    //! Simple meson sequential source construction
    /*! @ingroup hadron
     *
     * Create a simple meson sequential propagator source
     */
    class SimpleMesonSeqSource : public HadronSeqSource<LatticePropagator>
    {
    public:
      //! Full constructor
      SimpleMesonSeqSource(const Params& p, int gamma) : params(p), gamma_sink(gamma) {}

      //! Default destructor
      ~SimpleMesonSeqSource() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props) const;

    private:
      //! Hide partial constructor
      SimpleMesonSeqSource() {}

    private:
      Params  params;   /*!< Seqsource params */
      int gamma_sink;   /*!< The gamma matrix at the sink */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleMesonSeqSourceEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const SimpleMesonSeqSourceEnv::Params& param);


}  // end namespace Chroma

#endif
