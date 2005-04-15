// -*- C++ -*-
// $Id: inline_eigbnds.cc,v 1.4 2005-04-15 11:23:25 edwards Exp $
/*! \file
 * \brief Inline measurements for eigenvalue bounds
 *
 * Measure eigenvalue bounds of M^dag*M
 */

#include "meas/inline/eig/inline_eigbnds.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/eig/eig_spec.h"
#include "meas/eig/eig_spec_array.h"

#include "actions/ferm/linop/lmdagm.h"
#include "actions/ferm/linop/lopscl.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

namespace Chroma { 

  namespace InlineEigBndsMdagMEnv { 

    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineEigBndsMdagM(InlineEigBndsMdagMParams(xml_in, path));
    }

    const std::string name = "EIGBNDSMDAGM";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Ritz input
  void read(XMLReader& xml, const string& path, 
	    InlineEigBndsMdagMParams::RitzParams_t& param)
  {
    XMLReader paramtop(xml, path);
    
    read(paramtop, "RsdR", param.RsdR);
    read(paramtop, "RsdA", param.RsdA);
    read(paramtop, "RsdZero", param.RsdZero);
    read(paramtop, "ProjApsiP", param.ProjApsiP);
    read(paramtop, "Nmin", param.Nmin);
    read(paramtop, "MaxCG", param.MaxCG);
    read(paramtop, "Nrenorm", param.Nrenorm);
  }


  // Param stuff
  InlineEigBndsMdagMParams::InlineEigBndsMdagMParams() { frequency = 0; }

  InlineEigBndsMdagMParams::InlineEigBndsMdagMParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      read(paramtop, "Frequency", frequency);
      read(paramtop, "Ritz", ritz);

      if (paramtop.count("usePV") == 1)
	read(paramtop, "usePV", usePV);
      else
	usePV = false;

      // Generic Wilson-Type stuff
      string fa;
      read(paramtop, "./FermionAction/FermAct", fa);
      fermact = TheFermionActionFactory::Instance().createObject(fa,
								 paramtop,
								 string("./FermionAction"));
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  // "Do" helper on a 4D action
  void 
  InlineEigBndsMdagM::do4d(Handle< const LinearOperator<LatticeFermion> > MM,
			   unsigned long update_no,
			   XMLWriter& xml_out) 
  {
    QDPIO::cout << "4D eig bnds" << endl;

    push(xml_out, "EigBndsMdagM");
    write(xml_out, "update_no", update_no);

    int n_eig = 1;
    multi1d<Real> lambda(n_eig);
    multi1d<Real> check_norm(n_eig);
    multi1d<LatticeFermion> psi(n_eig);
    int n_CG_count;

    for(int i =0; i < n_eig; i++)
      gaussian(psi[i]);

    QDPIO::cout << "Look for lowest ev" << endl;

    push(xml_out,"LowestEv");
    EigSpecRitzCG(*MM, 
		  lambda, 
		  psi, 
		  n_eig,
		  params.ritz.Nrenorm, 
		  params.ritz.Nmin, 
		  params.ritz.MaxCG,
		  params.ritz.RsdR,
		  params.ritz.RsdA,
		  params.ritz.RsdZero,
		  params.ritz.ProjApsiP,
		  n_CG_count,
		  xml_out);
    pop(xml_out); // LowestEv

    QDPIO::cout << "Look for highest ev" << endl;
    Handle<const LinearOperator<LatticeFermion> > MinusMM(new lopscl<LatticeFermion, Real>(MM, Real(-1.0)));
  
    // Look for highest ev
    for(int i =0; i < n_eig; i++)
      gaussian(psi[i]);
    
    push(xml_out,"HighestEv");
    EigSpecRitzCG(*MinusMM,
		  lambda,
		  psi,
		  n_eig,
		  params.ritz.Nrenorm,
		  params.ritz.Nmin,
		  params.ritz.MaxCG,
		  params.ritz.RsdR,
		  params.ritz.RsdA,
		  params.ritz.RsdZero,
		  params.ritz.ProjApsiP,
		  n_CG_count,
		  xml_out);
    pop(xml_out); // HighestEv

    pop(xml_out); // pop("EigBndsMdagM");
  } 



  // "Do" helper on a 5D action
  void 
  InlineEigBndsMdagM::do5d(Handle< const LinearOperator< multi1d<LatticeFermion> > > MM,
			   unsigned long update_no,
			   XMLWriter& xml_out) 
  {
    QDPIO::cout << "5D eig bnds" << endl;

    push(xml_out, "EigBndsMdagM");
    write(xml_out, "update_no", update_no);

    int n_eig = 1;
    const int N5 = MM->size();
    multi1d<Real> lambda(n_eig);
    multi1d<Real> check_norm(n_eig);
    multi2d<LatticeFermion> psi(n_eig, N5);
    int n_CG_count;

    for(int i=0; i < n_eig; i++)
      for(int n=0; n < N5; n++)
	gaussian(psi[i][n]);

    QDPIO::cout << "Look for lowest ev" << endl;

    push(xml_out,"LowestEv");
    EigSpecRitzCG(*MM, 
		  lambda, 
		  psi, 
		  n_eig,
		  params.ritz.Nrenorm, 
		  params.ritz.Nmin, 
		  params.ritz.MaxCG,
		  params.ritz.RsdR,
		  params.ritz.RsdA,
		  params.ritz.RsdZero,
		  params.ritz.ProjApsiP,
		  n_CG_count,
		  xml_out);
    pop(xml_out); // LowestEv

    QDPIO::cout << "Look for highest ev" << endl;
    Handle<const LinearOperator< multi1d<LatticeFermion> > > MinusMM(new lopscl< multi1d<LatticeFermion>, Real>(MM, Real(-1.0)));
  
    // Look for highest ev
    for(int i=0; i < n_eig; i++)
      for(int n=0; n < N5; n++)
	gaussian(psi[i][n]);

    push(xml_out,"HighestEv");
    EigSpecRitzCG(*MinusMM,
		  lambda,
		  psi,
		  n_eig,
		  params.ritz.Nrenorm,
		  params.ritz.Nmin,
		  params.ritz.MaxCG,
		  params.ritz.RsdR,
		  params.ritz.RsdA,
		  params.ritz.RsdZero,
		  params.ritz.ProjApsiP,
		  n_CG_count,
		  xml_out);
    pop(xml_out); // HighestEv

    pop(xml_out); // pop("EigBndsMdagM");
  } 




  // The "do" function
  void 
  InlineEigBndsMdagM::operator()(const multi1d<LatticeColorMatrix>& u,
				 XMLBufferWriter& gauge_xml,
				 unsigned long update_no,
				 XMLWriter& xml_out) 
  {
    // Check success of the downcast 
    // Possible actions
    const FermAct4D<LatticeFermion>* S_4d = 
      dynamic_cast<const FermAct4D<LatticeFermion>*>(params.fermact.operator->()); // get pointer

    // Possible actions
    const FermAct5D<LatticeFermion>* S_5d = 
      dynamic_cast<const FermAct5D<LatticeFermion>*>(params.fermact.operator->()); // get pointer

    Handle< const ConnectState > connect_state(params.fermact->createState(u));

    if (S_4d != 0)
    {
      Handle< const LinearOperator<LatticeFermion> > MM(S_4d->lMdagM(connect_state));
      this->do4d(MM, update_no, xml_out);
    }
    else if (S_5d != 0)
    {
      if (! params.usePV)
      {
	//! Find evs of base operator
	Handle< const LinearOperator< multi1d<LatticeFermion> > > MM(S_5d->lMdagM(connect_state));
	this->do5d(MM, update_no, xml_out);
      }
      else
      {
	//! Find evs of PV operator
	Handle< const LinearOperator< multi1d<LatticeFermion> > > 
	  MM(new lmdagm< multi1d<LatticeFermion> >(S_5d->linOpPV(connect_state)));

	this->do5d(MM, update_no, xml_out);
      }
    }
    else
    {
      throw string("InlineEigBndsMdagM: no suitable cast found");
    }
  } 

};
