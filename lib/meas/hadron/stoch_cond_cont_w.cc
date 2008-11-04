// $Id: stoch_cond_cont_w.cc,v 3.2 2008-11-04 18:43:57 edwards Exp $
/*! \file
 * \brief Stoch quark condensates
 */

#include "meas/hadron/stoch_cond_cont_w.h"
#include "meas/hadron/hadron_contract_factory.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"

#include "meas/inline/io/named_objmap.h"


namespace Chroma 
{ 

  // Read parameters
  void read(XMLReader& xml, const string& path, StochCondContEnv::Params& param)
  {
    StochCondContEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const StochCondContEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  namespace StochCondContEnv 
  { 
    namespace
    {
      HadronContract* mesStochCondCont(XMLReader& xml_in, 
				       const std::string& path) 
      {
	return new StochCondCont(Params(xml_in, path));
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
	/**************************************************************************/
	break;

      default :
	/**************************************************************************/

	QDPIO::cerr << "StochCond: Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", mom2_max);
      read(paramtop, "avg_equiv_mom", avg_equiv_mom);
      read(paramtop, "mom_origin", mom_origin);
      read(paramtop, "soln_files", soln_files);
    }


    // Reader for input parameters
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", mom2_max);
      write(xml, "avg_equiv_mom", avg_equiv_mom);
      write(xml, "mom_origin", mom_origin);
      write(xml, "soln_files", soln_files);

      pop(xml);
    }


    //--------------------------------------------------------------

    //! Structure holding a source and its solutions
    struct QuarkSourceSolutions_t
    {
      //! Structure holding solutions
      struct QuarkSolution_t
      {
	LatticeFermion     source;
	LatticeFermion     soln;
	PropSourceConst_t  source_header;
	ChromaProp_t       prop_header;
      };

      int   j_decay;
      Seed  seed;
      multi1d<QuarkSolution_t>  dilutions;
    };


    //--------------------------------------------------------------
    // Construct some condensates
    std::list< Handle<HadronContractResult_t> >
    StochCondCont::operator()(const multi1d<LatticeColorMatrix>& u,
			      const std::string& xml_group,
			      const std::string& id_tag)
    {
      START_CODE();

      QDPIO::cout << "Stochastic Condensates" << endl;

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      StopWatch swatch;

      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);

      //
      // Read the source and solutions
      //
      swatch.reset();
      swatch.start();
      QuarkSourceSolutions_t  quark;

      try
      {
	QDPIO::cout << "Attempt to read solutions" << endl;
	quark.dilutions.resize(params.soln_files.size());

	QDPIO::cout << "dilutions.size= " << quark.dilutions.size() << endl;
	for(int i=0; i < quark.dilutions.size(); ++i)
	{
	  XMLReader file_xml, record_xml;

	  QDPIO::cout << "reading file= " << params.soln_files[i] << endl;
	  QDPFileReader from(file_xml, params.soln_files[i], QDPIO_SERIAL);
	  read(from, record_xml, quark.dilutions[i].soln);
	  close(from);
	
	  read(record_xml, "/Propagator/PropSource", quark.dilutions[i].source_header);
	  read(record_xml, "/Propagator/ForwardProp", quark.dilutions[i].prop_header);
	}
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting headers: " << e << endl;
	QDP_abort(1);
      }
      swatch.stop();

      QDPIO::cout << "Sources and solutions successfully read: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;



      //
      // Check for each quark source that the solutions have their diluted
      // on every site only once
      //
      swatch.reset();
      swatch.start();

      try
      {
	bool first = true;
	int  N;
	LatticeFermion quark_noise;      // noisy source on entire lattice

	for(int i=0; i < quark.dilutions.size(); ++i)
	{
	  std::istringstream  xml_s(quark.dilutions[i].source_header.source.xml);
	  XMLReader  sourcetop(xml_s);
//	QDPIO::cout << "Source = " << quark.dilutions[i].source_header.source.id << endl;

	  if (quark.dilutions[i].source_header.source.id != DiluteZNQuarkSourceConstEnv::getName())
	  {
	    QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::getName() << endl;
	    QDP_abort(1);
	  }

	  QDPIO::cout << "Dilution num= " << i << endl;

	  // Manually create the params so I can peek into them and use the source constructor
	  DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
							 quark.dilutions[i].source_header.source.path);
	  DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
	  if (first) 
	  {
	    first = false;

	    quark.j_decay = srcParams.j_decay;

	    // Grab N
	    N = srcParams.N;

	    // Set the seed to desired value
	    quark.seed = srcParams.ran_seed;
	    QDP::RNG::setrn(quark.seed);

	    // Create the noisy quark source on the entire lattice
	    zN_src(quark_noise, N);
	  }

	  // The seeds must always agree - here the seed is the unique id of the source
	  if ( toBool(srcParams.ran_seed != quark.seed) )
	  {
	    QDPIO::cerr << "dilution=" << i << " seed does not match" << endl;
	    QDP_abort(1);
	  }

	  // The N's must always agree
	  if ( toBool(srcParams.N != N) )
	  {
	    QDPIO::cerr << "dilution=" << i << " N does not match" << endl;
	    QDP_abort(1);
	  }

	  // Use a trick here, create the source and subtract it from the global noisy
	  // Check at the end that the global noisy is zero everywhere.
	  // NOTE: the seed will be set every call
	  quark.dilutions[i].source = srcConst(u);
	  quark_noise -= quark.dilutions[i].source;

#if 0
	  // Diagnostic
	  {
	    // Keep a copy of the phases with NO momenta
	    SftMom phases_nomom(0, true, quark.dilutions[i].source_header.j_decay);

	    multi1d<Double> source_corr = sumMulti(localNorm2(quark.dilutions[i].source), 
						   phases_nomom.getSet());
	      
	    multi1d<Double> soln_corr = sumMulti(localNorm2(quark.dilutions[i].soln), 
						 phases_nomom.getSet());

	  }
#endif
	} // end for i

	Double dcnt = norm2(quark_noise);
	if (toDouble(dcnt) != 0.0)  // problematic - seems to work with unnormalized sources 
	{
	  QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	  QDP_abort(1);
	}
      } // end try
      catch(const std::string& e) 
      {
	QDPIO::cerr << ": Caught Exception creating source: " << e << endl;
	QDP_abort(1);
      }

      swatch.stop();

      QDPIO::cout << "Sources saturated: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;

      //
      // Condensates
      //
      // Parameters needed for the momentum projection
      SftMomParams_t sft_params;

      sft_params.origin_offset.resize(Nd);
      sft_params.origin_offset = 0;
      sft_params.mom2_max      = params.mom2_max;
      sft_params.mom_offset    = params.mom_origin;
      sft_params.avg_equiv_mom = params.avg_equiv_mom;
      sft_params.decay_dir     = quark.j_decay;

      // Start operator contractions
      swatch.reset();
      swatch.start();

      std::list< Handle<Hadron2PtContract_t> > hadron;   // holds the contract lattice correlator

      for(int gamma_value=0; gamma_value < Ns*Ns; ++gamma_value)
      {
	Handle<Hadron2PtContract_t> had(new Hadron2PtContract_t);

	push(had->xml, xml_group);
	write(had->xml, id_tag, "stoch_diagonal_gamma_condensates");
	write(had->xml, "gamma_value", gamma_value);
	write(had->xml, "mom2_max", sft_params.mom2_max);
	write(had->xml, "avg_equiv_mom", sft_params.avg_equiv_mom);
	write(had->xml, "decay_dir", sft_params.decay_dir);
	pop(had->xml);

	had->corr = zero;
	for(int i=0; i < quark.dilutions.size(); ++i)
	{
	  had->corr += 
	    localInnerProduct(quark.dilutions[i].source, Gamma(gamma_value) * quark.dilutions[i].soln);
	} // end for i
	
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
	success &= Chroma::TheHadronContractFactory::Instance().registerObject(string("stoch_diagonal_gamma_condensates"),
									       mesStochCondCont);

	registered = true;
      }
      return success;
    }

  } // namespace StochCondContEnv

} // namespace Chroma
