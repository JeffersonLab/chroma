// $Id: simple_baryon_operator_w.cc,v 1.2 2006-05-21 04:40:21 edwards Exp $
/*! \file
 *  \brief Construct simple baryon operators
 */

#include "meas/hadron/simple_baryon_operator_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/hadron/barspinmat_w.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SimpleBaryonOperatorEnv::Params& param)
  {
    SimpleBaryonOperatorEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const SimpleBaryonOperatorEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace SimpleBaryonOperatorEnv
  { 

    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      try
      {
	XMLReader paramtop(xml, path);

	int version;
	read(paramtop, "version", version);

	switch (version) 
	{
	case 2:
	  break;
	  
	default:
	  QDPIO::cerr << name << ": parameter version " << version 
		      << " unsupported." << endl;
	  QDP_abort(1);
	}

	{
	  XMLReader xml_tmp(paramtop, "SourceQuarkSmearing");
	  std::ostringstream os;
	  xml_tmp.print(os);
	  read(xml_tmp, "wvf_kind", source_quark_smearing_type);
	  source_quark_smearing = os.str();
	}

	{
	  XMLReader xml_tmp(paramtop, "SinkQuarkSmearing");
	  std::ostringstream os;
	  xml_tmp.print(os);
	  read(xml_tmp, "wvf_kind", sink_quark_smearing_type);
	  sink_quark_smearing = os.str();
	}

	{
	  XMLReader xml_tmp(paramtop, "LinkSmearing");
	  std::ostringstream os;
	  xml_tmp.print(os);
	  read(xml_tmp, "LinkSmearingType", link_smearing_type);
	  link_smearing = os.str();
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception in params: " << e << endl;
	QDP_abort(1);
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "BaryonOperatorType", SimpleBaryonOperatorEnv::name);
      xml << source_quark_smearing;
      xml << sink_quark_smearing;
      xml << link_smearing;

      pop(xml);
    }



    //! Full constructor
    BarNuclCg5::BarNuclCg5(const Params& p, const multi1d<LatticeColorMatrix>& u_) : 
      params(p), u_smr(u_)
    {
      // Factory constructions
      try
      {
	// Smear the gauge field if needed
        QDPIO::cout << "Link smearing type = " << params.link_smearing_type << endl;
	linkSmear(u_smr, std::string("/LinkSmearing"), 
		  params.link_smearing, params.link_smearing_type);


	// Create the source quark smearing object
	{
	  std::istringstream  xml_s(params.source_quark_smearing);
	  XMLReader  smeartop(xml_s);
	  const string smear_path = "/SourceQuarkSmearing";
	
	  sourceQuarkSmearing = 
	    TheFermSmearingFactory::Instance().createObject(params.source_quark_smearing_type,
							    smeartop,
							    smear_path);
	}

	// Create the sink quark smearing object
	{
	  std::istringstream  xml_s(params.sink_quark_smearing);
	  XMLReader  smeartop(xml_s);
	  const string smear_path = "/SinkQuarkSmearing";
	
	  sinkQuarkSmearing = 
	    TheFermSmearingFactory::Instance().createObject(params.sink_quark_smearing_type,
							    smeartop,
							    smear_path);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception smearing: " << e << endl;
	QDP_abort(1);
      }
    }


    //! Manipulate the quark fields
    void
    BarNuclCg5::quarkManip(multi1d<LatticeFermion>& q,
			   const LatticeFermion& q1, 
			   const LatticeFermion& q2, 
			   const LatticeFermion& q3,
			   enum PlusMinus isign) const
    {
      START_CODE();

      q.resize(3);
      q[0] = q1;
      q[1] = q2;
      q[2] = q3;

      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      switch (isign)
      {
      case PLUS:
	// Sink smear the quarks
	for(int i=0; i < q.size(); ++i)
	  (*sinkQuarkSmearing)(q[i], u_smr);
	break;

      case MINUS:
	// Source smear the quarks
	for(int i=0; i < q.size(); ++i)
	  (*sourceQuarkSmearing)(q[i], u_smr);
	break;

      default:
	QDPIO::cerr << name << ": illegal isign" << endl;
	QDP_abort(1);
      }
  
      END_CODE();
    }


    //! Compute the operator
    multi1d<LatticeComplex> 
    BarNuclCg5::operator()(const LatticeFermion& q1, 
			   const LatticeFermion& q2, 
			   const LatticeFermion& q3,
			   enum PlusMinus isign) const
    {
      START_CODE();

      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      multi1d<LatticeFermion> q;
      quarkManip(q, q1, q2, q3, isign);
  
      // The return
      multi1d<LatticeComplex> d(Ns);
      d = zero;

      // C gamma_5 = Gamma(5)
      SpinMatrix Cg5 = BaryonSpinMats::Cg5();

      for(int k=0; k < Ns; ++k)
      {
	LatticeSpinMatrix di_quark = zero;

	for(int j=0; j < Ns; ++j)
	{
	  for(int i=0; i < Ns; ++i)
	  {
	    // Contract over color indices with antisym tensors
	    LatticeComplex b_oper = colorContract(peekSpin(q[0], i),
						  peekSpin(q[1], j),
						  peekSpin(q[2], k));

	    pokeSpin(di_quark, b_oper, j, i);
	  }
	}

	d[k] += traceSpin(Cg5 * di_quark);
      }

      END_CODE();

      return d;
    }


    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! Nucleon = (u C gamma_5 d) u
      /*!
       * \ingroup hadron
       *
       * C gamma_5 = Gamma(5) = - (C gamma_5)^T
       */
      BaryonOperator<LatticeFermion>* barNuclCg5(XMLReader& xml_in,
						 const std::string& path,
						 const multi1d<LatticeColorMatrix>& u)
      {
	return new BarNuclCg5(Params(xml_in, path), u);
      }

    }  // end anonymous namespace


    //! Baryon operators
    /*! \ingroup hadron */
    bool registerAll(void) 
    {
      bool success = true;

      //! Register all the factories
      success &= Chroma::TheWilsonBaryonOperatorFactory::Instance().registerObject(name,
										   barNuclCg5);

      return success;
    }

    const std::string name = "NUCLEON";

    const bool registered = registerAll();

  } // namespace BaryonOperatorCallMapEnv


}  // end namespace Chroma
