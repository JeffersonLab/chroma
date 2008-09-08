// $Id: t_solver_accum.cc,v 1.2 2008-09-08 18:41:16 bjoo Exp $
/*! \file
 *  \brief Test 4d fermion actions
 */

#include <iostream>
#include <cstdio>

#include "chroma.h"


#include "update/molecdyn/monomial/read_rat_approx.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_aggregate.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_aggregate.h"

using namespace Chroma;

//! To insure linking of code, place the registered code flags here
/*! This is the bit of code that dictates what fermacts are in use */
bool linkage_hack()
{
  bool foo = true;
  // All actions
  foo &= WilsonTypeFermActsEnv::registerAll();
  foo &= WilsonTypeFermActs5DEnv::registerAll();
  foo &= MdagMMultiSysSolverEnv::registerAll();
  foo &= MdagMMultiSysSolverArrayEnv::registerAll();
  foo &= MdagMMultiSysSolverAccumulateEnv::registerAll();
  foo &= MdagMMultiSysSolverAccumulateArrayEnv::registerAll();
  
  return foo;
}



int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  QDPIO::cout << "linkage=" << linkage_hack() << endl;


  multi1d<int> nrow(Nd);
  for(int i=0; i < Nd-1; i++) { 
    nrow[i] = 4;
  }
  nrow[Nd-1] = 16;

  // Specify lattice size, shape, etc.
  Layout::setLattSize(nrow);
  Layout::create();

  QDPIO::cout << "t_temp_prec" << endl;
  
  // Start up a weak field
  struct Cfg_t config = { CFG_TYPE_WEAK_FIELD, "dummy" };
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, config);
  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);
  
  // Instantiate XML writer for XMLDAT
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "t_temp_prec");
  proginfo(xml_out);    // Print out basic program info


  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();  
  pop(xml_out);

  XMLReader xml_in("./t_solver_accum.ini.xml");
  ReadRatApproxEnv::Params rat_approx_param;
  GroupXML_t fermact_group, fermact_5d_group,
    inv_param_group;

  try {
    XMLReader paramtop(xml_in, "/Params");
    fermact_group=readXMLGroup(paramtop, 
					  "./FermionAction",
					  "FermAct");

    fermact_5d_group = readXMLGroup(paramtop, 
					  "./FermionAction5D",
					  "FermAct");

    QDPIO::cout << "fermact group read" << endl;
     inv_param_group=readXMLGroup(paramtop, 
					    "./InvertParam", 
					    "invType");

     QDPIO::cout << "inv_param_group read" << endl;

    
    read(paramtop, "./Approximation", rat_approx_param);
    QDPIO::cout << "Approx Read" << endl;
  }
  catch(const std::string& e) { 
    QDPIO::cerr << "XML Error returned : " << e << endl;
    QDP_abort(1);
  }

  // Typedefs to save typing
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;
  
  try {
    std::istringstream fermact_xml_is( fermact_group.xml);
    XMLReader fermact_reader(fermact_xml_is);
    // Create unprec and Prec actions
    Handle< FermAct4D<T,P,Q> >
      S(  TheWilsonTypeFermActFactory::Instance().createObject(fermact_group.id,
							       fermact_reader,
							       fermact_group.path)     );
    
  
   // Now create a FermState
    Handle< FermState<T, P, Q> >  fs( S->createState(u) );
    Handle< MdagMMultiSystemSolver<T> > multi_solver( S->mInvMdagM(fs, inv_param_group) );
    Handle< MdagMMultiSystemSolverAccumulate<T> > multi_solver_acc( S->mInvMdagMAcc(fs, inv_param_group) );

    Handle< LinearOperator<T> > M(S->linOp(fs));


    const Subset& subset = M->subset();
    
    LatticeFermion rhs;
    gaussian(rhs, subset);           // Use this as a source
    
    
    LatticeFermion psi1, psi2; // The two answers
    
    RemezCoeff_t pfe, ipfe;
    ReadRatApprox approx_reader(rat_approx_param);
    approx_reader(pfe, ipfe);
    
    // Work with pfe
    QDPIO::cout << "Working With PFE:  A = " << pfe.norm << endl;
    for(int i=0; i < pfe.pole.size(); i++) { 
      QDPIO::cout << "  res["<<i<<"]=" << pfe.res[i] << " pole[" << i 
		  << "]=" << pfe.pole[i] << endl;
    }

    {
      {
	
	multi1d<LatticeFermion> temp_ferms(pfe.pole.size());
	for(int i=0; i < temp_ferms.size(); i++) { 
	  temp_ferms[i][subset] = zero;
	}
	(*multi_solver)(temp_ferms, pfe.pole,rhs);
	

	psi1[subset] = pfe.norm*rhs;

	for(int i=0; i < temp_ferms.size(); i++) { 
	  psi1[subset] += pfe.res[i] * temp_ferms[i];
	}
      }
    
    /* Do the accumulated solver */
    (*multi_solver_acc)(psi2, 
			pfe.norm, 
			pfe.res,
			pfe.pole,
			rhs);
    
    
    // Compare psi1 and psi2;
    LatticeFermion diff;
    diff[subset] = psi1-psi2;
    Double n1 = sqrt(norm2(diff, subset));
    Double n2 = sqrt(norm2(psi1, subset));
    
    QDPIO::cout << "|| psi1 - psi2 ||=" << n1 << endl;
    QDPIO::cout << "|| psi1 - psi2 ||/||psi1||= " << n1/n2 << endl;

    }



    // Now Array versions
    std::istringstream fermact_5d_xml_is( fermact_5d_group.xml);
    XMLReader fermact_5d_reader(fermact_5d_xml_is);

    Handle< FermAct5D<T,P,Q> >
      S5(  TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_5d_group.id,
							       fermact_5d_reader,
							       fermact_5d_group.path)     );
    
  
   // Now create a FermState
    Handle< FermState<T, P, Q> >  fs5( S5->createState(u) );
    Handle< MdagMMultiSystemSolverArray<T> > multi_solver_5d( S5->mInvMdagM(fs, inv_param_group) );
    Handle< MdagMMultiSystemSolverAccumulateArray<T> > multi_solver_acc_5d( S5->mInvMdagMAcc(fs, inv_param_group) );

    Handle< LinearOperatorArray<T> > M5(S5->linOp(fs));
  
    
    const Subset& sub5 = M5->subset();
    const int N5 = M5->size();

    multi1d<LatticeFermion> psi5(N5);
    multi1d<LatticeFermion> psi5_2(N5);
    multi1d<LatticeFermion> rhs5(N5);
    for(int n=0; n < N5; n++) { 
      psi5[n][sub5] = zero;
      psi5_2[n][sub5] = zero;
      gaussian(rhs5[n],sub5);
    }

    {
      multi1d<multi1d<LatticeFermion> > tmp5(pfe.pole.size());
      for(int i=0; i< pfe.pole.size(); i++) { 
	tmp5[i].resize(N5);
	for(int n=0; n < N5; n++) { 
	  tmp5[i][n][sub5]=zero;
	}
      }

      (*multi_solver_5d)(tmp5, pfe.pole, rhs5);

      for(int n=0; n < N5; n++) {
	psi5[n][sub5] = pfe.norm*rhs5[n];
      }
      for(int i=0; i < pfe.pole.size(); i++) { 
	for(int n=0; n < N5; n++) { 
	  psi5[n][sub5] += pfe.res[i] * tmp5[i][n];
	}
      }
    }

    (*multi_solver_acc_5d)(psi5_2, pfe.norm, pfe.res, pfe.pole, rhs5);

    multi1d<LatticeFermion> diff5(N5);
    for(int n=0; n < N5; n++) { 
      diff5[n][sub5] = psi5_2[n] - psi5[n];
    }
    Double norm5_diff = norm2(diff5, sub5);
    Double norm5_psi = norm2(psi5, sub5);
    QDPIO::cout << "|| psi5_1 - psi5_2 ||=" << norm5_diff << endl;
    QDPIO::cout << "|| psi5_1 - psi5_2 ||/||psi5_1||= " << norm5_diff/norm5_psi << endl;

  } 
  catch(const std::string& e) { 
    QDPIO::cerr << "Caught exception: " << e << endl;
  }

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
