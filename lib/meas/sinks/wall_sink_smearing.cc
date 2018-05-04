/*! \file
 *  \brief Wall sink smearing
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/wall_sink_smearing.h"
#include "util/ft/sftmom.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const std::string& path, WallQuarkSinkSmearingEnv::Params& param)
  {
    WallQuarkSinkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const WallQuarkSinkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace WallQuarkSinkSmearingEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceSink<LatticePropagator>* createProp(XMLReader& xml_in,
						     const std::string& path,
						   const multi1d<LatticeColorMatrix>& u)
      {
	return new SinkSmear<LatticePropagator>(Params(xml_in, path), u);
      }

      //! Callback function
      QuarkSourceSink<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
								  const std::string& path,
								  const multi1d<LatticeColorMatrix>& u)
      {
	return new SinkSmear<LatticeStaggeredPropagator>(Params(xml_in, path), u);
      }

      //! Callback function
      QuarkSourceSink<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  const multi1d<LatticeColorMatrix>& u)
      {
	return new SinkSmear<LatticeFermion>(Params(xml_in, path), u);
      }

      //! Name to be used
      const std::string name("WALL_SINK");

      //! Local registration flag
      bool registered = false;
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropSinkSmearingFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheStagPropSinkSmearingFactory::Instance().registerObject(name, createStagProp);
	success &= Chroma::TheFermSinkSmearingFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }


    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }

      read(paramtop, "j_decay", j_decay);

      // Sanity check
      if (j_decay < 0 || j_decay >= Nd)
      {
	QDPIO::cerr << name << ": invalid params.j_decay=" << j_decay << std::endl;
	QDP_abort(1);
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      int version = 1;
      push(xml, path);
      write(xml, "version", version);
      write(xml, "SinkType", WallQuarkSinkSmearingEnv::getName());
      write(xml, "j_decay", j_decay);
      pop(xml);
    }


    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticePropagator>::operator()(LatticePropagator& quark_sink) const
    {
      QDPIO::cout << "Wall sink" << std::endl;
 
      // Project onto zero mom. at each time slice, the put back into
      // original field
      SftMom phases(0, true, params.j_decay);
      multi1d<DPropagator> slice_prop(sumMulti(quark_sink, phases.getSet()));

      for(int t=0; t < phases.numSubsets(); ++t)
	quark_sink[phases.getSet()[t]] = slice_prop[t];

      return;
    }



    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark_sink) const
    {
      QDPIO::cout << "Wall sink" << std::endl;
 
      // Project onto zero mom. at each time slice, the put back into
      // original field
      SftMom phases(0, true, params.j_decay);
      multi1d<DStaggeredPropagator> slice_prop(sumMulti(quark_sink, phases.getSet()));

      for(int t=0; t < phases.numSubsets(); ++t)
	quark_sink[phases.getSet()[t]] = slice_prop[t];

      return;
    }



    //! Construct the sink smearing
    template<>
    void
    SinkSmear<LatticeFermion>::operator()(LatticeFermion& quark_sink) const
    {
      QDPIO::cout << "Wall sink" << std::endl;
 
      // Project onto zero mom. at each time slice, the put back into
      // original field
      SftMom phases(0, true, params.j_decay);
      multi1d<DFermion> slice_prop(sumMulti(quark_sink, phases.getSet()));

      for(int t=0; t < phases.numSubsets(); ++t)
	quark_sink[phases.getSet()[t]] = slice_prop[t];

      return;
    }
  }

}
