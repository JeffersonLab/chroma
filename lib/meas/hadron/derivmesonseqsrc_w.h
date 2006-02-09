// -*- C++ -*-
// $Id: derivmesonseqsrc_w.h,v 2.1 2006-02-09 02:25:24 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#ifndef __mesonseqsrc_w_h__
#define __mesonseqsrc_w_h__

namespace Chroma 
{

  //! Meson sequential sources
  /*! \ingroup hadron */
  namespace MesonSeqSourceCallMapEnv
  { 
    extern bool registered;   // forward decl
  }

}  // end namespace Chroma


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
      SimpleMesonSeqSource(int gamma_sink_) : gamma_sink(gamma_sink_) {}

      //! Default destructor
      ~SimpleMesonSeqSource() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& forward_props,
				   const multi1d<int>& sink_mom, 
				   int t_sink, int j_decay) const;

    private:
      //! Hide partial constructor
      SimpleMesonSeqSource() {}

    private:
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
