// -*- C++ -*-
// $Id: inline_eigbnds.cc,v 3.1 2006-04-11 04:18:23 edwards Exp $
/*! \file
 * \brief Inline measurements for eigenvalue bounds
 *
 * Measure eigenvalue bounds of M^dag*M
 */

#include "meas/inline/eig/inline_eigbnds.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/eig/eig_spec.h"
#include "meas/eig/eig_spec_array.h"
#include "meas/inline/io/named_objmap.h"

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
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, 
										   createMeasurement);
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

  //! Ritz output
  void write(XMLWriter& xml, const string& path, 
	     InlineEigBndsMdagMParams::RitzParams_t& param)
  {
    push(xml, path);
    
    write(xml, "RsdR", param.RsdR);
    write(xml, "RsdA", param.RsdA);
    write(xml, "RsdZero", param.RsdZero);
    write(xml, "ProjApsiP", param.ProjApsiP);
    write(xml, "Nmin", param.Nmin);
    write(xml, "MaxCG", param.MaxCG);
    write(xml, "Nrenorm", param.Nrenorm);

    pop(xml);
  }


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineEigBndsMdagMParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineEigBndsMdagMParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }


  // Param stuff
  InlineEigBndsMdagMParams::InlineEigBndsMdagMParams()
  { 
    frequency = 0; 
  }

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

      XMLReader xml_tmp(paramtop, "./FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      ferm_act = os.str();
 
      // Ids
      read(paramtop, "NamedObject", named_obj);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  // Writer
  void
  InlineEigBndsMdagMParams::write(XMLWriter& xml, const std::string& path) 
  {
    push(xml, path);

    QDP::write(xml, "Frequency", frequency);
    QDP::write(xml, "usePV", usePV);
    xml << ferm_act;
    Chroma::write(xml, "Ritz", ritz);

    // Ids
    Chroma::write(xml, "NamedObject", named_obj);

    pop(xml);
  }


  // Constructor
  InlineEigBndsMdagM::InlineEigBndsMdagM(const InlineEigBndsMdagMParams& p) : params(p) 
  {
    std::istringstream is(params.ferm_act);
    XMLReader fermact_reader(is);

    // Get the name of the ferm act
    string fa;
    try 
    { 
      read(fermact_reader, "/FermionAction/FermAct", fa);
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    // Generic Wilson-Type stuff
    fermact = TheFermionActionFactory::Instance().createObject(fa,
							       fermact_reader,
							       string("./FermionAction"));
  }

  // "Do" helper on a 4D action
  void 
  InlineEigBndsMdagM::do4d(Handle< LinearOperator<LatticeFermion> > MM,
			   unsigned long update_no,
			   XMLWriter& xml_out) 
  {
    QDPIO::cout << "4D eig bnds" << endl;

    push(xml_out, "EigBndsMdagM");
    write(xml_out, "update_no", update_no);
    params.write(xml_out, "Input");

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
    Handle< LinearOperator<LatticeFermion> > MinusMM(new lopscl<LatticeFermion, Real>(MM, Real(-1.0)));
  
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
  InlineEigBndsMdagM::do5d(Handle< LinearOperatorArray<LatticeFermion> > MM,
			   unsigned long update_no,
			   XMLWriter& xml_out) 
  {
    QDPIO::cout << "5D eig bnds" << endl;

    push(xml_out, "EigBndsMdagM");
    write(xml_out, "update_no", update_no);
    params.write(xml_out, "Input");

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
    Handle< LinearOperatorArray<LatticeFermion> > MinusMM(new lopsclArray<LatticeFermion, Real>(MM, Real(-1.0)));
  
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
  InlineEigBndsMdagM::operator()(unsigned long update_no,
				 XMLWriter& xml_out) 
  {
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // Check success of the downcast 
    // Possible actions
    FermAct4D<T,P,Q>* S_4d = 
      dynamic_cast< FermAct4D<T,P,Q>*>(fermact.operator->()); // get pointer

    // Possible actions
    FermAct5D<T,P,Q>* S_5d = 
      dynamic_cast< FermAct5D<T,P,Q>*>(fermact.operator->()); // get pointer


    multi1d<LatticeColorMatrix> u;
    try
    {
      // Grab the object
      u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
    }
    catch (std::bad_cast) 
    {
      QDPIO::cerr << InlineEigBndsMdagMEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineEigBndsMdagMEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    

    Handle< FermState<T,P,Q> > connect_state(fermact->createState(u));

    if (S_4d != 0)
    {
      Handle< LinearOperator<LatticeFermion> > MM(S_4d->lMdagM(connect_state));
      this->do4d(MM, update_no, xml_out);
    }
    else if (S_5d != 0)
    {
      if (! params.usePV)
      {
	//! Find evs of base operator
	Handle< LinearOperatorArray<LatticeFermion> > MM(S_5d->lMdagM(connect_state));
	this->do5d(MM, update_no, xml_out);
      }
      else
      {
	//! Find evs of PV operator
	Handle< LinearOperatorArray<LatticeFermion> > 
	  MM(new MdagMLinOpArray<LatticeFermion>(S_5d->linOpPV(connect_state)));

	this->do5d(MM, update_no, xml_out);
      }
    }
    else
    {
      QDPIO::cerr << InlineEigBndsMdagMEnv::name << ": no suitable cast found" << endl;
      QDP_abort(1);
    }
  } 

}
