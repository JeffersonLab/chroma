// $Id: group_baryon_operator_w.cc,v 1.5 2006-05-21 04:49:51 edwards Exp $
/*! \file
 *  \brief Construct group baryon operators
 */

#include "meas/hadron/group_baryon_operator_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/displacement.h"
#include "util/ferm/diractodr.h"

#include <map>

using std::map;

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, GroupBaryonOperatorEnv::Params& param)
  {
    GroupBaryonOperatorEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const GroupBaryonOperatorEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace GroupBaryonOperatorEnv
  { 

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

      read(paramtop, "operator_coeff_file", operator_coeff_file);
      read(paramtop, "displacement_length", displacement_length);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "BaryonOperatorType", GroupBaryonOperatorEnv::name);
      write(xml, "operator_coeff_file", operator_coeff_file);
      write(xml, "displacement_length", displacement_length);
      xml << source_quark_smearing;
      xml << sink_quark_smearing;
      xml << link_smearing;

      pop(xml);
    }



    //! Full constructor
    GroupBaryon::GroupBaryon(const Params& p, const multi1d<LatticeColorMatrix>& u_) : 
      params(p), u_smr(u_)
    {
      readCoeffs(coeffs);

      // The spin basis matrix to goto Dirac
      rotate_mat = adj(DiracToDRMat());

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


    //! Reader
    void GroupBaryon::readCoeffs(multi1d< multi1d< CoeffTerm_t > >& coef)
    {
      TextReader reader(params.operator_coeff_file);

      int num_rows;
      reader >> num_rows;
      coeffs.resize(num_rows);

      for(int n=0; n < coeffs.size(); ++n)
      {
	int num_terms;
	reader >> num_terms;

	coeffs[n].resize(num_terms);
	
	for(int l=0; l < coeffs[n].size(); ++l)
	{
	  CoeffTerm_t& term = coeffs[n][l];
	  term.quark.resize(3);

	  // Make spin index 0 based
	  {
	    multi1d<int> spin(3);
	    reader >> spin[0] >> spin[1] >> spin[2];

	    for(int i=0; i < spin.size(); ++i)
	      term.quark[i].spin = spin[i] - 1;
	  }

	  // Convert displacements to  disp_dir, disp_len
	  {
	    multi1d<int> displacement(3);
	    reader >> displacement[0] >> displacement[1] >> displacement[2];

	    for(int i=0; i < displacement.size(); ++i)
	    {
	      term.quark[i].displacement = displacement[i];

	      if (displacement[i] == 0)
	      {
		term.quark[i].disp_dir = 0;
		term.quark[i].disp_len = 0;
	      }
	      else if (displacement[i] > 0)
	      {
		term.quark[i].disp_dir = term.quark[i].displacement - 1;
		term.quark[i].disp_len = params.displacement_length;
	      }
	      else
	      {
		term.quark[i].disp_dir = -term.quark[i].displacement - 1;
		term.quark[i].disp_len = -params.displacement_length;
	      }
	    }
	  }

	  // Read the garbage around a complex
	  {
	    Real re, im;
	    char lparen, comma, rparen;

	    reader >> lparen >> re >> comma >> im >> rparen;
	    term.coeff = cmplx(re,im);
	  }
	}
      }
	
      reader.close();
    }


    //! Construct array of maps of displacements
    void
    GroupBaryon::displaceQuarks(multi1d< map<int,LatticeFermion> >& disp_quarks,
				const multi1d<LatticeFermion>& q,
				enum PlusMinus isign) const
    {
      START_CODE();

//      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      disp_quarks.resize(3);

      for(int n=0; n < coeffs.size(); ++n)
      {
	for(int l=0; l < coeffs[n].size(); ++l)
	{
	  const CoeffTerm_t& term = coeffs[n][l];
	    
	  for(int i=0; i < disp_quarks.size(); ++i)
	  {
	    // Make some shorthands to ease my brain
	    map<int,LatticeFermion>& disp_q = disp_quarks[i];
	    const CoeffTerm_t::QuarkTerm_t& term_q = term.quark[i];

	    // If no entry, then create a displaced version of the quark
	    if (disp_q.find(term_q.displacement) == disp_q.end())
	    {
//	      cout << __func__ 
//		   << ": n=" << n
//		   << " l=" << l
//		   << " i=" << i 
//		   << " disp=" << term.quark[i].displacement
//		   << " len=" << term.quark[i].disp_len
//		   << " dir=" << term.quark[i].disp_dir
//		   << endl;

	      LatticeFermion qq = q[i];

	      switch (isign)
	      {
	      case PLUS:
		displacement(u_smr, qq, term_q.disp_len, term_q.disp_dir);
		break;

	      case MINUS:
		displacement(u_smr, qq, -term_q.disp_len, term_q.disp_dir);
		break;
	      }

	      // Insert
	      disp_q.insert(std::make_pair(term_q.displacement, qq));
	    }
	  } // for i
	} // for l
      }  // for n

//      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    }
    

    //! First displace then smear the quarks
    void
    GroupBaryon::displaceSmearQuarks(multi1d< map<int,LatticeFermion> >& disp_quarks,
				     const LatticeFermion& q1, 
				     const LatticeFermion& q2, 
				     const LatticeFermion& q3,
				     enum PlusMinus isign) const
    {
      START_CODE();

//      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      multi1d<LatticeFermion> q(3);

      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[0] = rotateMat() * q1;
      q[1] = rotateMat() * q2;
      q[2] = rotateMat() * q3;

      // Displace
      displaceQuarks(disp_quarks, q, isign);


      // Source smear after the displacements
      for(int i=0; i < disp_quarks.size(); ++i)
      {
	// Make some shorthands to ease my brain
	map<int,LatticeFermion>& disp_q = disp_quarks[i];

	// Loop over all keys
	for(std::map<int,LatticeFermion>::const_iterator mm = disp_q.begin();
	    mm != disp_q.end();
	    ++mm)
	{
	  (*sourceQuarkSmearing)(disp_q[mm->first], u_smr);
	}
      }

//      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    }


    //! First smear then displace the quarks
    void
    GroupBaryon::smearDisplaceQuarks(multi1d< map<int,LatticeFermion> >& disp_quarks,
				     const LatticeFermion& q1, 
				     const LatticeFermion& q2, 
				     const LatticeFermion& q3,
				     enum PlusMinus isign) const
    {
      START_CODE();

//      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      multi1d<LatticeFermion> q(3);

      // Rotate the quarks up-front.
      // NOTE: this step assumes the spin rotation commutes with the
      // quark smearing and the displacement
      // However, we are careful about the ordering of the smearing and
      // the displacement
      q[0] = rotateMat() * q1;
      q[1] = rotateMat() * q2;
      q[2] = rotateMat() * q3;

      // Sink smear the quarks
      for(int i=0; i < q.size(); ++i)
	(*sinkQuarkSmearing)(q[i], u_smr);

      // Displace after the smearing
      displaceQuarks(disp_quarks, q, isign);

//      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    }


    //! Manipulate the quark fields
    void
    GroupBaryon::quarkManip(multi1d< map<int,LatticeFermion> >& disp_quarks,
			    const LatticeFermion& q1, 
			    const LatticeFermion& q2, 
			    const LatticeFermion& q3,
			    enum PlusMinus isign) const
    {
      START_CODE();

//      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      switch (isign)
      {
      case PLUS:
	// Sink
	smearDisplaceQuarks(disp_quarks, q1, q2, q3, isign);
	break;

      case MINUS:
	// Source
	displaceSmearQuarks(disp_quarks, q1, q2, q3, isign);
	break;

      default:
	QDPIO::cerr << name << ": illegal isign" << endl;
	QDP_abort(1);
      }
  
//      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      END_CODE();
    }


    //! Compute the operator
    multi1d<LatticeComplex> 
    GroupBaryon::operator()(const LatticeFermion& q1, 
			    const LatticeFermion& q2, 
			    const LatticeFermion& q3,
			    enum PlusMinus isign) const
    {
      START_CODE();

//      cout << __PRETTY_FUNCTION__ << ": entering" << endl;

      // The result of displace and smearing (in some unspecified order here)
      multi1d< map<int,LatticeFermion> > disp_quarks;

      // Depending on whether this is the sink or source, do the appropriate
      // combination of smearing and displacing
      quarkManip(disp_quarks, q1, q2, q3, isign);
  
      // The return
      multi1d<LatticeComplex> d(coeffs.size());

      for(int n=0; n < coeffs.size(); ++n)
      {
	d[n] = zero;

	for(int l=0; l < coeffs[n].size(); ++l)
	{
	  const CoeffTerm_t& term = coeffs[n][l];
	    
	  // Contract over color indices with antisym tensors
	  LatticeComplex b_oper = 
	    colorContract(peekSpin(disp_quarks[0].find(term.quark[0].displacement)->second,
				   term.quark[0].spin),
			  peekSpin(disp_quarks[1].find(term.quark[1].displacement)->second,
				   term.quark[1].spin),
			  peekSpin(disp_quarks[2].find(term.quark[2].displacement)->second,
				   term.quark[2].spin));
	  
	  d[n] += term.coeff * b_oper;
	}
      }

      END_CODE();

//      cout << __PRETTY_FUNCTION__ << ": exiting" << endl;

      return d;
    }



    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! Call-back function
      BaryonOperator<LatticeFermion>* groupBaryon(XMLReader& xml_in,
						  const std::string& path,
						  const multi1d<LatticeColorMatrix>& u)
      {
	return new GroupBaryon(Params(xml_in, path), u);
      }

    }  // end anonymous namespace


    //! Baryon operators
    /*! \ingroup hadron */
    bool registerAll(void) 
    {
      bool success = true;

      // Required stuff
      success &= LinkSmearingEnv::registered;
      success &= QuarkSmearingEnv::registered;

      //! Register all the factories
      success &= Chroma::TheWilsonBaryonOperatorFactory::Instance().registerObject(name, 
										   groupBaryon);

      return success;
    }

    const std::string name = "GROUP_BARYON";

    const bool registered = registerAll();

  } // namespace BaryonOperatorCallMapEnv


}  // end namespace Chroma
