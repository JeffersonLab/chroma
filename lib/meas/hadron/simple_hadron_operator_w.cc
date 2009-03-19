// $Id: simple_hadron_operator_w.cc,v 1.2 2009-03-19 17:17:20 mcneile Exp $
/*! \file
 *  \brief Construct simple baryon operators
 */

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3



#include "meas/hadron/simple_hadron_operator_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/hadron/barspinmat_w.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

using namespace std;
namespace Chroma 
{

  //! Hadron sequential sources
  /*! \ingroup hadron */
  namespace SimpleHadronOperatorEnv
  { 
    map<string, HadronOperator<LatticeFermion>* (*)(const GroupXML_t&)> Ops ;
    
    Baryon::Baryon(const GroupXML_t& p):HadronOperator<LatticeFermion>(p){
      
      std::istringstream  xml_l(p.xml);
      XMLReader  xmltop(xml_l);
      QDPIO::cout << "Baryon state is = " <<p.id ; 
      QDPIO::cout << endl;
      string DiqGammaStr ;
      read(xmltop,"DiqGamma",DiqGammaStr); 
      if(DiqGammaStr == "CG5"){
	DiqGamma =  BaryonSpinMats::Cg5();
      }
      else if(DiqGammaStr == "CGmu"){
	int mu ;
	read(xmltop,"mu",mu);
	DiqGamma =  BaryonSpinMats::Cgmu(mu); // mu = 1 2 3 4 . 4 is time
      }
      else{
	throw "Unknown diquark operator "+DiqGammaStr ;
      }
    }


    //! Compute the operator
    multi1d<LatticeComplex> 
    Baryon::operator()(const multi1d<LatticeFermion>& q, 
		       enum PlusMinus isign) const
    {
      START_CODE();

      // The return
      multi1d<LatticeComplex> d(Ns);
      d = zero;

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

	d[k] += traceSpin(DiqGamma * di_quark);
      }

      END_CODE();

      return d;
    }


    //! Anonymous namespace
    namespace
    {

      
      //-------------------- callback functions ------------------------------

      //! Baryon operator
      /*!
       * \ingroup hadron
       *
       * simple baryon operator
       */
      HadronOperator<LatticeFermion>* baryon(const GroupXML_t& gxml){
	return new Baryon(gxml);
      }


      //! Local registration flag
      bool registered = false;

    }  // end anonymous namespace


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	Ops["BARYON"] =  &baryon ;

	registered = true;
      }
      return success;
    }

  } // namespace HadronOperatorEnv


}  // end namespace Chroma

#endif
#endif
#endif
