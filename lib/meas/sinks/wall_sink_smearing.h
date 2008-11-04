// -*- C++ -*-
// $Id: wall_sink_smearing.h,v 1.3 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Wall sink smearing
 */

#ifndef __wall_sink_smearing_h__
#define __wall_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sinks */
  namespace WallQuarkSinkSmearingEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();


    //! Wall sink parameters
    /*! @ingroup sinks */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      int j_decay;    /*! Direction of decay. The wall projection is orthog to this dir */
    };


    //! Wall sink smearing
    /*! @ingroup sinks
     *
     * Make a wall propagator sink
     */
    template<typename T>
    class SinkSmear : public QuarkSourceSink<T>
    {
    public:
      //! Full constructor
      /*! The gauge field is not used, so discard */
      SinkSmear(const Params& p, const multi1d<LatticeColorMatrix>& u) : params(p) {}

      //! Smear the sink
      void operator()(T& obj) const;

    private:
      //! Hide partial constructor
      SinkSmear() {}

    private:
      Params  params;                         /*!< sink params */
    };

  } // end namespace


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, WallQuarkSinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const WallQuarkSinkSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
