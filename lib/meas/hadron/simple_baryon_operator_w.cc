// $Id: simple_baryon_operator_w.cc,v 1.1 2006-05-12 03:38:01 edwards Exp $
/*! \file
 *  \brief Construct simple baryon operators
 */

#include "meas/hadron/simple_baryon_operator_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/hadron/barspinmat_w.h"

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
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
      pop(xml);
    }



    //! Full constructor
    BarNuclCg5::BarNuclCg5(const Params& p, const multi1d<LatticeColorMatrix>& u_) : params(p), u(u_)
    {
    }


    //! Compute the operator
    multi1d<LatticeComplex> 
    BarNuclCg5::operator()(const LatticeFermion& q1, 
			   const LatticeFermion& q2, 
			   const LatticeFermion& q3,
			   enum PlusMinus isign) const
    {
      START_CODE();

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
	    LatticeComplex b_oper = colorContract(peekSpin(q1, i),
						  peekSpin(q2, j),
						  peekSpin(q3, k));

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
      success &= Chroma::TheWilsonBaryonOperatorFactory::Instance().registerObject(string("NUCLEON"), 
										   barNuclCg5);

      return success;
    }

    const bool registered = registerAll();

  } // namespace BaryonOperatorCallMapEnv


}  // end namespace Chroma
