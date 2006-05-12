// $Id: group_baryon_operator_w.cc,v 1.1 2006-05-12 03:38:01 edwards Exp $
/*! \file
 *  \brief Construct group baryon operators
 */

#include "meas/hadron/group_baryon_operator_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/hadron/barspinmat_w.h"
#include "meas/smear/displacement.h"

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
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
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
      write(xml, "displacement_length", displacement_length);
      pop(xml);
    }



    //! Full constructor
    GroupBaryon::GroupBaryon(const Params& p, const multi1d<LatticeColorMatrix>& u_) : params(p), u(u_)
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

	  // Make spin index 0 based
	  {
	    multi1d<int> spin(3);
	    reader >> spin[0] >> spin[1] >> spin[2];
	    term.spin = spin;

	    for(int i=0; i < spin.size(); ++i)
	      term.spin -= 1;
	  }

	  // Convert displacements to  disp_dir, disp_len
	  {
	    multi1d<int> displacement(3);
	    reader >> displacement[0] >> displacement[1] >> displacement[2];

	    term.disp_dir.resize(displacement.size());
	    term.disp_len.resize(displacement.size());

	    for(int i=0; i < displacement.size(); ++i)
	    {
	      if (displacement[i] == 0)
	      {
		term.disp_dir[i] = 0;
		term.disp_len[i] = 0;
	      }
	      else if (displacement[i] > 0)
	      {
		term.disp_dir[i] = displacement[i] - 1;
		term.disp_len[i] = params.displacement_length;
	      }
	      else
	      {
		term.disp_dir[i] = -displacement[i] - 1;
		term.disp_len[i] = -params.displacement_length;
	      }
	    }
	  }

	  // Complex coeff;  /* Oh bleeping darn yuk */
	  // std::string complex_coeff;

	  // Complication - hack around the parentheses in the operator file
	  // Assume it is been preprocessed to space separated floats
	  {
	    Real re, im;

	    reader >> re >> im;
	    term.coeff = cmplx(re,im);
	  }
	}
      }
	
      reader.close();
    }


    //! Compute the operator
    multi1d<LatticeComplex> 
    GroupBaryon::operator()(const LatticeFermion& q1, 
			    const LatticeFermion& q2, 
			    const LatticeFermion& q3,
			    enum PlusMinus isign) const
    {
      START_CODE();

      multi1d<LatticeComplex> d(coeffs.size());

      for(int n=0; n < coeffs.size(); ++n)
      {
	d[n] = zero;

	for(int l=0; l < coeffs[n].size(); ++l)
	{
	  const CoeffTerm_t& term = coeffs[n][l];
	    
	  multi1d<LatticeFermion> q(3);
	  q[0] = q1;
	  q[1] = q2;
	  q[2] = q3;

	  for(int i=0; i < q.size(); ++i)
	  {
	    switch (isign)
	    {
	    case PLUS:
	      displacement(u, q[i], term.disp_len[i], term.disp_dir[i]);
	      break;

	    case MINUS:
	      displacement(u, q[i], -term.disp_len[i], term.disp_dir[i]);
	      break;

	    default:
	      QDP_error_exit("illegal isign in GroupBaryon");
	    }
	  }

	  // Contract over color indices with antisym tensors
	  LatticeComplex b_oper = colorContract(peekSpin(q[0], term.spin[0]),
						peekSpin(q[1], term.spin[1]),
						peekSpin(q[2], term.spin[2]));

	  d[n] += term.coeff * b_oper;
	}
      }

      END_CODE();

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

      //! Register all the factories
      success &= Chroma::TheWilsonBaryonOperatorFactory::Instance().registerObject(string("GROUP_BARYON"), 
										   groupBaryon);

      return success;
    }

    const bool registered = registerAll();

  } // namespace BaryonOperatorCallMapEnv


}  // end namespace Chroma
