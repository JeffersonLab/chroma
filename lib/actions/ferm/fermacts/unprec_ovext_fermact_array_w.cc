// $Id: unprec_ovext_fermact_array_w.cc,v 1.14 2005-01-02 05:21:10 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovext_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/invert/invcg2_array.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecOvExtFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
										      const std::string& path)
    {
      return new UnprecOvExtFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					 UnprecOvExtFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_OVEXT";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	   & Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Read parameters
  UnprecOvExtFermActArrayParams::UnprecOvExtFermActArrayParams(XMLReader& xml, 
							       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "OverMass", OverMass);
    read(paramtop, "Mass", Mass);
    read(paramtop, "N5", N5);

    if (paramtop.count("a5") != 0) 
      read(paramtop, "a5", a5);
    else
      a5 = 1.0;
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecOvExtFermActArrayParams& param)
  {
    UnprecOvExtFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* 
  UnprecOvExtFermActArray::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecOvExtLinOpArray(state->getLinks(),OverMass,Mass,N5);
  }


  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator< multi1d<LatticeFermion> >* 
  UnprecOvExtFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }


  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*! \ingroup qprop
   *
   * Propagator solver for Extended overlap fermions
   */
  template<typename T>
  class OvExt5DQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param Mass_      quark mass ( Read )
     */
    OvExt5DQprop(Handle< const LinearOperator< multi1d<T> > > A_,
		 const Real& Mass_,
		 const InvertParam_t& invParam_) : 
      A(A_), Mass(Mass_), invParam(invParam_) {}

    //! Destructor is automatic
    ~OvExt5DQprop() {}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    int operator() (T& psi, const T& chi) const
    {
      START_CODE();

      const int  N5 = A->size();   // array size better match
      int n_count;
  
      int G5 = Ns*Ns - 1;

      // Initialize the 5D fields
      multi1d<LatticeFermion> chi5(N5);
      multi1d<LatticeFermion> psi5(N5);
      psi5 = zero;
      chi5 = zero;

      psi5[0] = psi;
      chi5[0] = Gamma(G5) * chi;

      switch(invParam.invType)
      {
      case CG_INVERTER: 
	// psi5 = (H_o)^(-2) chi5
	InvCG2(*A, chi5, psi5, invParam.RsdCG, invParam.MaxCG, n_count);

	// chi5 = H_o * (H_o)^(-2) * gamma_5 * chi
	(*A)(chi5, psi5, MINUS);
	break;
  
      case MR_INVERTER:
      case BICG_INVERTER:
	QDP_error_exit("Unsupported inverter type", invParam.invType);
	break;
  
      default:
	QDP_error_exit("Unknown inverter type", invParam.invType);
      }
      
      if ( n_count == invParam.MaxCG )
	QDP_error_exit("no convergence in the inverter", n_count);
  
      // Overall normalization
      Real ftmp1 = Real(1) / Real(1 - Mass);

      // Normalize and remove contact term
      psi = ftmp1*(chi5[0] - chi);

      END_CODE();

      return n_count;
    }

  private:
    // Hide default constructor
    OvExt5DQprop() {}

    Handle< const LinearOperator< multi1d<T> > > A;
    const Real Mass;
    const InvertParam_t invParam;
  };

 
  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  const SystemSolver<LatticeFermion>* 
  UnprecOvExtFermActArray::qprop(Handle<const ConnectState> state,
				 const InvertParam_t& invParam) const
  {
    return new OvExt5DQprop<LatticeFermion>(Handle< const LinearOperator< multi1d<LatticeFermion> > >(linOp(state)),
					    getQuarkMass(),
					    invParam);
  }
  


}
