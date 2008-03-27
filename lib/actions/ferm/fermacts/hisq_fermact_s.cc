// $Id: hisq_fermact_s.cc,v 1.8 2008-03-27 10:17:34 mcneile Exp $
/*! \file
 *  \brief Hisq staggered fermion action
 */

/*

Highly improved staggered quarks on the lattice, with applications to char\m physics.
By HPQCD Collaboration and UKQCD Collaboration (E. Follana et al.). Oct 20\06. 21pp.
Published in Phys.Rev.D75:054502,2007.
e-Print: hep-lat/0610092

Note
-----

My current understanding is that Hisq and Asqtad then
differ by the way the fat links are constructed.

This is meant to be the minimum modification to add
HISQ to chroma. I want to try not do too much cut and pasting
of Asqtad code, so that the chroma class police
are kept happy. Reuse and all that!


The correction to the Naik term, known as epsilon
(equation 24 in the above paper) is included now as an external argument.
This is the way that Christine wanted it done. I assume that 
epsilon will be computed to one loop order in perturbation
theory (one day).
*/

#include "chromabase.h"

#include "actions/ferm/fermacts/fermact_factory_s.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_s.h"

#include "actions/ferm/linop/asqtad_mdagm_s.h"
#include "actions/ferm/linop/asqtad_linop_s.h"
#include "actions/ferm/fermacts/hisq_fermact_s.h"
#include "util/gauge/stag_phases_s.h"
#include "util/gauge/reunit.h"
#include "util/gauge/sun_proj.h"
#include "meas/gfix/polar_dec.h"

// DEBUG
#include "util/gauge/unit_check.h"
#include "meas/glue/mesplq.h"

namespace Chroma 
{ 

  //! Hooks to register the class with the fermact factory
  namespace HisqFermActEnv
  {
    //! Callback function
    StaggeredTypeFermAct<LatticeStaggeredFermion,
			 multi1d<LatticeColorMatrix>,
			 multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
									const std::string& path)
    {
      return new HisqFermAct(StaggeredCreateFermStateEnv::reader(xml_in, path), 
			       HisqFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeStaggeredFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "HISQ";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheStagFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheStagTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
  }


  //! Produce a linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the entire lattice
   *
   * \param u_fat, u_triple 	 fat7 and triple links    (Read)
   * \u has already had KS phases multiplied in.
   */
  EvenOddLinearOperator<LatticeStaggeredFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* 
  HisqFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {

    // Why in fact are we casting to the base class on both sides of
    // this assignment ? The answer is so that we can use a proxy.
    // Both the Proxy and the ConnectState inherit from the BaseClass
    // and can be cast to and from the base class. However the Proxy
    // and the connect state cannot be directly cast into each other.
    // Which is why we have a virtual base class in the first place.
    //
    // So We cast the ConnectState to an AsqtadConnectStateBase
    // This we can do at our leisure from either AsqtadConnectState
    // OR from the Proxy. We then get access to all the virtual methods
    // in the AsqtadConnectState. Only Restriction: We have to use the
    // get() methods as they are all the base class provides.
    return new AsqtadLinOp(state.cast<AsqtadConnectStateBase>(), param.Mass);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the checkerboarded lattice
   *
   * \param u_fat, u_triple 	 fat7 and triple links	       (Read)
   */
  DiffLinearOperator<LatticeStaggeredFermion, 
		     multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> >* 
  HisqFermAct::lMdagM(Handle< FermState<T,P,Q> > state) const
  {
    return new AsqtadMdagM(state.cast<AsqtadConnectStateBase>(), param.Mass);
  }


  //! Create a state -- this multiplies in the K-S phases computes the fat and triple links etc
  AsqtadConnectStateBase*
  HisqFermAct::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    multi1d<LatticeColorMatrix> u_with_phases(Nd);
    multi1d<LatticeColorMatrix> u_fat(Nd);
    multi1d<LatticeColorMatrix> u_fat_I(Nd);
    multi1d<LatticeColorMatrix> u_triple(Nd);

     QDPIO::cout << "HISQ setting up the fat links\n"  ;

    // First put in the BC
    u_with_phases = u_;
    getFermBC().modify(u_with_phases);

#if 0 
    // DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
    fat7_param pp ; 
     QDPIO::cout << "HISQ hacked to do ASQTAD" << endl ; 
    pp.c_1l = (Real)(5) / (Real)(8);
    pp.c_3l = (Real)(-1) / ((Real)(16));
    pp.c_5l = - pp.c_3l / ((Real)(4));
    pp.c_7l = - pp.c_5l / ((Real)(6));
    pp.c_Lepage = pp.c_3l ; 
    Fat7_Links(u_with_phases, u_fat, pp);
    Real one = (Real) 1.0 ; 
    Triple_Links(u_with_phases, u_triple, one);
    // DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
#endif
    Real ep = param.epsilon ; 
#if 0
    QDPIO::cout << "HISQ epsilon = " << ep << "\n"; ;
    QDPIO::cout << "HISQ Mass = "    <<  param.Mass << "\n"; ;
    QDPIO::cout << "HISQ u0 = "      <<  param.u0 << "\n"; ;
#endif
    QDPIO::cout << "Setting up HISQ inverter\n"; ;

    // Create Fat7 links. This uses the same
    // coefficients as Asqtad, but with zero
    // Lepage term
    fat7_param pp ; 


    pp.c_1l = (Real)(1)/(Real)(8) ;   // lepage contributes here
    pp.c_3l = (Real)(1) / ((Real)(16));
    pp.c_5l = pp.c_3l / ((Real)(4));    // 1/64
    pp.c_7l = pp.c_5l / ((Real)(6));    // -1/384
    pp.c_Lepage =  0.0 ; 


    Fat7_Links(u_with_phases, u_fat_I, pp);

    // reunitarise (using polar method)
   LatticeReal  alpha ; // complex phase (not needed here)
   Real  JacAccu = 0.00000000001 ;
   int JacMax = 100 ; 
   LatticeColorMatrix  w ;
   QDPIO::cout << "SU3 polar projection Accuracy " << JacAccu << " max iters = " <<  JacMax << endl;

   for(int i = 0; i < Nd; i++) 
     {
       w = u_fat_I[i] ;
       polar_dec(w, u_fat_I[i],alpha, JacAccu, JacMax) ;
     }

   // with HISQ the three links are fat
   Real UU0 = (Real) 1.0 ;  // tadpole 
   Real c_3 ; 
   c_3 = (Real)(-1 - ep) / (Real)(24);
   Triple_Links(u_fat_I, u_triple, UU0, c_3);
   //   Triple_Links(u_fat_I, u_triple, UU0);

    // fatten again with different coefficient of fat term
   pp.c_1l = (Real)(8 + ep) / (Real)(8)  ;
   pp.c_3l = (Real)(1) / ((Real)(16)); //  -0.0625 or -1/16
   pp.c_5l = pp.c_3l / ((Real)(4));   //   0.01562500 or 1/64
   pp.c_7l = pp.c_5l / ((Real)(6));   // .00260416666 or 1/384
   pp.c_Lepage = -2.0 * pp.c_3l ;  // double Lepage term  -1/8

    Fat7_Links(u_fat_I, u_fat, pp);

   for(int i = 0; i < Nd; i++) 
     {
        u_fat[i]     *= StagPhases::alpha(i);
        u_triple[i]  *= StagPhases::alpha(i);
     }

    // ---------------------------------------------------

    return new AsqtadConnectState(cfs->getFermBC(), u_with_phases, u_fat, u_triple);
  }

} // End Namespace Chroma

