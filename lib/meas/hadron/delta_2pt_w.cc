// $Id: delta_2pt_w.cc,v 3.1 2007-12-04 04:08:53 kostas Exp $
/*! \file
 *  \brief Construct meson 2pt correlators.
 */

#include "meas/hadron/delta_2pt_w.h"
#include "meas/hadron/hadron_contract_factory.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, Delta2PtEnv::Params& param)
  {
    Delta2PtEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const Delta2PtEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Meson correlators
  /*! 
   * \ingroup hadron 
   *
   * @{
   */
  namespace Delta2PtEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------

      //! Construct pion correlator
      HadronContract* mesDeltaCorrs(XMLReader& xml_in,
					const std::string& path)
      {
	return new DeltaCorrs(Params(xml_in, path));   // all gammas
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


    // Construct all the correlators
    std::list< Handle<HadronContractResult_t> >
    DeltaCorrs::operator()(const multi1d<LatticeColorMatrix>& u,
				    const std::string& xml_group,
				    const std::string& id_tag)
    {
      START_CODE();

      QDPIO::cout << "Hadron2Pt: Delta" << endl;

      multi1d<ForwardProp_t> forward_headers(2);
      forward_headers[0] = readForwardPropHeader(params.first_id);
      forward_headers[1] = readForwardPropHeader(params.second_id);
      
      multi1d<int> t_srce = getTSrce(forward_headers);
      int decay_dir       = getDecayDir(forward_headers);

      // Get references to the props
      const LatticePropagator& quark_prop1 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.first_id);
      const LatticePropagator& quark_prop2 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.second_id);

      // Parameters needed for the momentum projection
      SftMomParams_t sft_params;
      sft_params.mom2_max      = params.mom2_max;
      sft_params.origin_offset = t_srce;
      sft_params.mom_offset    = params.mom_origin;
      sft_params.avg_equiv_mom = params.avg_equiv_mom;
      sft_params.decay_dir     = decay_dir;

      std::list< Handle<Hadron2PtContract_t> > hadron;   // holds the contract lattice correlator

      for(int gamma_value=0; gamma_value < Ns*Ns; ++gamma_value)
      {
	Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);

	push(had->xml, xml_group);
	write(had->xml, id_tag, "delta");
	write(had->xml, "gamma_value", gamma_value);
	write(had->xml, "PropHeaders", forward_headers);
	pop(had->xml);

	//here is where the contraction goes
	//had->corr = mesXCorr(quark_prop1, quark_prop2, gamma_value);

	hadron.push_back(had);  // push onto end of list
      }

      END_CODE();

      return this->project(hadron, sft_params);
    }


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Chroma::TheHadronContractFactory::Instance().registerObject(string("Delta"), mesDeltaCorrs);

	registered = true;
      }
      return success;
    }

  }  // end namespace Delta2PtEnv

  /*! @} */   // end of group io

}  // end namespace Chroma


  
