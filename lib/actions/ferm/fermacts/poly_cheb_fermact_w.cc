// $Id: poly_cheb_fermact_w.cc,v 2.1 2006-02-09 22:27:01 edwards Exp $
/*! \file
 *  \brief Chebyshev polynomial fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/poly_cheb_fermact_w.h"
#include "actions/ferm/linop/polprec_op.h"
#include "actions/ferm/linop/polynomial_op.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace PolyChebFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
										    const std::string& path)
    {
      return new PolyChebFermAct(WilsonTypeFermBCEnv::reader(xml_in, path), 
				 PolyChebFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "POLYNOMIAL_CHEBYSHEV";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	   & Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Read parameters
  PolyChebFermActParams::PolyChebFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the chebyshev params
    {
      XMLReader chebtop(xml, "PolyParams");
      read(chebtop, "degree", polyParams.degree);
      read(chebtop, "UpperBound", polyParams.UpperBound);
      read(chebtop, "LowerBound", polyParams.LowerBound);
      read(chebtop, "order", polyParams.order);
    }

    // Read the underlying fermion action
    try 
    { 
      if(in.count("AuxFermAct") == 1 )
      { 
	XMLReader xml_tmp(in, "AuxFermAct");
	std::ostringstream os;
	xml_tmp.print(os);
	AuxFermAct = os.str();
      }
      else
      {
	throw std::string("No auxilliary action");
      }
    }
    catch(const string& e) {
      QDPIO::cerr << "Caught Exception reading Chebyshev Fermact params: " << e << endl;
      QDP_abort(1);
    }

  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, PolyChebFermActParams& param)
  {
    PolyChebFermActParams tmp(xml, path);
    param = tmp;
  }


  //! Initializer
  void
  PolyChebFermAct::init()
  {
    QDPIO::cout << "Constructing PolyChebFermAct from params" << endl;
    std::istringstream  xml_s(param.AuxFermAct);
    XMLReader  fermacttop(xml_s);
    const string fermact_path = "/AuxFermAct";

    struct UnprecCastFailure {
      UnprecCastFailure(std::string e) : auxfermact(e) {};
      const string auxfermact;
    };

    try
    {
      string auxfermact;
      read(fermacttop, fermact_path + "/FermAct", auxfermact);
      QDPIO::cout << "AuxFermAct: " << auxfermact << endl;

      // Generic Wilson-Type stuff
      FermionAction<LatticeFermion>* S_f =
	TheFermionActionFactory::Instance().createObject(auxfermact,
							 fermacttop,
							 fermact_path);

      UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* S_aux; 
      S_aux = dynamic_cast<UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >*>(S_f);

      // Dumbass User specifies something that is not UnpreWilsonTypeFermAct
      // dynamic_cast MUST be checked for 0
      if( S_aux == 0 ) throw UnprecCastFailure(auxfermact);
     

      // Drop AuxFermAct into a Handle immediately.
      // This should free things up at the end
      Handle<UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > >  S_w(S_aux);
      fermact = S_w;
    }
    catch( const UnprecCastFailure& e) 
    {
      // Breakage Scenario
      QDPIO::cerr << "Unable to upcast auxiliary fermion action to "
		  << "UnprecWilsonTypeFermAct " << endl;
//      QDPIO::cerr << OvlapPartFrac4DFermActEnv::name << " does not support even-odd preconditioned "
//		  << "auxiliary FermActs" << endl;
      QDPIO::cerr << "You passed : " << endl;
      QDPIO::cerr << e.auxfermact << endl;
      QDP_abort(1);
    }
    catch (const std::exception& e) {
      // General breakage Scenario
      QDPIO::cerr << "Error reading data: " << e.what() << endl;
      throw;
    }

  }

    //! Produce a linear operator for this action
  const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >*
  PolyChebFermAct::linOp(Handle<const ConnectState> state) const
  {
    return fermact->linOp(state);
  }

  //! Produce a linear operator M^dag.M for this action
  const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >*
  PolyChebFermAct::lMdagM(Handle<const ConnectState> state) const
  {
    return fermact->lMdagM(state);
  }

  //! Produce the gamma_5 hermitian operator H_w
  const LinearOperator<LatticeFermion>*
  PolyChebFermAct::hermitianLinOp(Handle< const ConnectState> state) const
  {
    return fermact->hermitianLinOp(state);
  }

  //! Produce a linear operator for this action
  const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >*
  PolyChebFermAct::polyLinOp(Handle<const ConnectState> state) const
  {
    Handle< const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > 
      MdagM(fermact->lMdagM(state));

    return new lpoly< LatticeFermion, multi1d<LatticeColorMatrix> >(MdagM,
								    param.degree, 
								    param.LowerBound, 
								    param.UpperBound, 
								    param.order);
  }

  //! Produce a linear operator for this action
  const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >*
  PolyChebFermAct::polyPrecLinOp(Handle<const ConnectState> state) const
  {
    Handle< const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > 
      M(fermact->linOp(state));

    Handle< const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > 
      Pol(polyLinOp(state));

    return new lpoly< LatticeFermion, multi1d<LatticeColorMatrix> >(M, Pol);
  }

}
