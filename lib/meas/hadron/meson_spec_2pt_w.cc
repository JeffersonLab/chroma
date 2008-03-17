// $Id: meson_spec_2pt_w.cc,v 1.1 2008-03-17 15:23:58 edwards Exp $
/*! \file
 *  \brief Construct meson 2pt correlators leaving all spin indices open
 */

#error "STILL WORKING ON THIS"

#include "meas/hadron/meson_spec_2pt_w.h"
#include "meas/hadron/hadron_contract_factory.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, MesonSpec2PtEnv::Params& param)
  {
    MesonSpec2PtEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const MesonSpec2PtEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Meson correlators
  /*! 
   * \ingroup hadron 
   *
   * @{
   */
  namespace MesonSpec2PtEnv
  { 
    //! Anonymous namespace
    namespace
    {
      //! Construct the correlator
      LatticeComplex mesXCorr(const LatticePropagator& quark_prop_1,
			      const LatticePropagator& quark_prop_2, 
			      int gamma_value)
      {
	// Construct the anti-quark propagator from quark_prop_2
	int G5 = Ns*Ns-1;
	LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

	return LatticeComplex(trace(adj(anti_quark_prop) * (Gamma(gamma_value) *
							    quark_prop_1 * Gamma(gamma_value))));
      }


      //-------------------- callback functions ---------------------------------------

      //! Construct pion correlator
      HadronContract* mesDiagGammaCorrs(XMLReader& xml_in,
					const std::string& path)
      {
	return new MesonSpecCorrs(Params(xml_in, path));   // all gammas
      }


      //! Local registration flag
      bool registered = false;

    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
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
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", mom2_max);
      read(paramtop, "avg_equiv_mom", avg_equiv_mom);
      read(paramtop, "mom_origin", mom_origin);
      read(paramtop, "first_id", first_id);
      read(paramtop, "second_id", second_id);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "mom2_max", mom2_max);
      write(xml, "avg_equiv_mom", avg_equiv_mom);
      write(xml, "mom_origin", mom_origin);

      write(xml, "first_id", first_id);
      write(xml, "second_id", second_id);

      pop(xml);
    }


    //! Do a SFT and serialize the output
    void serializeSFT(BinaryWriter& bin, const SftMom& phases, const LatticeComplex& prop)
    {
      multi2d<DComplex> hsum(phases.sft(prop));
      int length  = phases.numSubsets();

      write(bin, phases.numMom());

      for(int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
      {
	multi1d<DComplex>  corr(length);   /*!< Momentum projected correlator */
	for (int t=0; t < length; ++t) 
	{
//        int t_eff = (t - t0 + length) % length;  // NOTE: no longer shifting source
	  corr[t] = hsum[sink_mom_num][t];
	}

	write(bin, phases.numToMom(sink_mom_num));
	write(bin, corr);
      }
    }


    //! Meson data
    struct MesonSpecData_t
    {
      LatticePropagator  quark_prop1;
      LatticePropagator  quark_prop2;
      Handle<SftMom>     phases;
    };

    //! Do some initialization
    void init(MesonSpecData_t& data,
	      XMLWriter& xml, const string& path, const string& id_tag, 
	      const Params& params)
    {
      MesonSpecData_t data;

      multi1d<ForwardProp_t> forward_headers(2);
      forward_headers[0] = readForwardPropHeader(params.first_id);
      forward_headers[1] = readForwardPropHeader(params.second_id);
      
      push(xml, path);
      write(xml, id_tag, "meson_spec");
      write(xml, "PropHeaders", forward_headers);
      pop(xml);

      multi1d<int> t_srce = getTSrce(forward_headers);
      int decay_dir       = getDecayDir(forward_headers);

      // Get references to the props
      data.quark_prop1 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.first_id);
      data.quark_prop2 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.second_id);

      // Parameters needed for the momentum projection
      SftMomParams_t sft_params;
      sft_params.mom2_max      = params.mom2_max;
      sft_params.origin_offset = t_srce;
      sft_params.mom_offset    = params.mom_origin;
      sft_params.avg_equiv_mom = params.avg_equiv_mom;
      sft_params.decay_dir     = decay_dir;

      data.phases(new SftMom(sft_params));
    }


    // Construct all the correlators
    std::list< Handle<HadronContractResult_t> >
    MesonSpecCorrs::operator()(const multi1d<LatticeColorMatrix>& u,
			       const std::string& xml_group,
			       const std::string& id_tag)
    {
      START_CODE();

      QDPIO::cout << "Hadron2Pt: diagonal_gamma_mesons" << endl;

      HadronContractResult_t corr;
      MesonSpecData_t data;
      init(data, corr->bin, xml_group, id_tag, params);

      // Length of lattice in decay direction
      int length  = data.phases->numSubsets();
      int G5      = Ns*Ns-1;

      LatticePropagator antiquark_1 = adj(Gamma(G5) * data.quark_prop1 * Gamma(G5));
      LatticeComplex m_prop;

      for(int sf_2=0; sf_2 < Ns; ++sf_2)           // sf_2
	for(int sf_1=0; sf_1 < Ns; ++sf_1)         // sf_1
	  for(int si_2=0; si_2 < Ns; ++si_2)       // si_2
	    for(int si_1=0; si_1 < Ns; ++si_1)     // si_1
	    {
	      write(bin, si_1);
	      write(bin, sf_1);
	      write(bin, si_2);
	      write(bin, sf_2);

	      // Contract over color indices with antisym tensors
	      m_prop = 
		traceColor(peekSpin(antiquark_1,
				    sf_1,si_1)   // (sf_1,si_1)
			   * peekSpin(quark_prop2,
				      sf_2,si_2)); // (sf_2,si_2)

	      serializeSFT(corr->bin, *(data.phases), m_prop);
	    }

      std::list< Handle<HadronContractResult_t> > corrs;
      corrs.push_back(corr);

      END_CODE();

      return corrs;
    }


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Chroma::TheHadronContractFactory::Instance().registerObject(string("meson_spec"),
									       mesDiagGammaCorrs);

	registered = true;
      }
      return success;
    }

  }  // end namespace MesonSpec2PtEnv

  /*! @} */   // end of group io

}  // end namespace Chroma


  
