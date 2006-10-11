// $Id: nef_quarkprop4_w.cc,v 3.3 2006-10-11 15:42:26 edwards Exp $
// $Log: nef_quarkprop4_w.cc,v $
// Revision 3.3  2006-10-11 15:42:26  edwards
// Removed use of numRetries !! This was a misguided attempt to get restarts
// into the code. Will try another way via the system solvers and the
// InvertParam xml group.
//
// Revision 3.2  2006/07/03 15:26:09  edwards
// Changed FermionAction API to produce system solver classes. No longer have
// a fixed InvertParam set. Removed old inverter enum and support. Have default
// creation for inverters to solve M*psi=chi, MdagM*psi=chi, and multi-shift
// support. Some older files still have explicit calls to CG solver.
//
// Revision 3.1  2006/06/11 06:30:32  edwards
// Change in interface. The quarkProp routines now all take a "numRetries"
// variable for the number of times to call the qprop routine. The propagator
// regression tests have all been up to version 10 to read this new variable.
// The SystemSolver routines now all return a  SystemSolverResults_t  struct
// instead of an int of "n_count" . A residual of the unpreconditioned
// system of equations for the qprop and qpropT is computed.
//
// Revision 3.0  2006/04/03 04:58:52  edwards
// Major overhaul of fermion and gauge action interface. Basically,
// all fermacts and gaugeacts now carry out  <T,P,Q>  template parameters. These are
// the fermion type, the "P" - conjugate momenta, and "Q" - generalized coordinates
// in the sense of Hamilton's equations. The fermbc's have been rationalized to never
// be over multi1d<T>. The "createState" within the FermionAction is now fixed meaning
// the "u" fields are now from the coordinate type. There are now "ConnectState" that
// derive into FermState<T,P,Q> and GaugeState<P,Q>.
//
// Revision 2.0  2005/09/25 21:04:30  edwards
// Moved to version 2.0
//
// Revision 1.14  2005/01/14 20:13:06  edwards
// Removed all using namespace QDP/Chroma from lib files. The library
// should now be 100% in the Chroma namespace. All mainprogs need a
// using namespace Chroma.
//
// Revision 1.13  2005/01/02 05:21:10  edwards
// Rearranged top-level fermion actions and linear operators.
// The deriv is pushed up much higher. This has the effect that
// now the "P" type (conjugate momentum type) is carried around
// in DiffFermAct4D<T,P> and below carry this "P" template param.
// All fermacts now have at least default deriv support (possibly
// Not Implemented error). Made qprop and qpropT now go through
// factory functions.
//
// Revision 1.12  2004/12/29 22:13:41  edwards
// Rearranged the top-level FermionAction structure. Moved up
// 5d-style actions to be a  FermAct5D and there is a FermAct4D.
// Consolidated quarkprop routines and fermact factories.
//
// Revision 1.11  2004/11/24 04:16:56  kostas
// code works now
//
// Revision 1.10  2004/11/20 22:28:45  edwards
// Narrowed include dependencies.
//
// Revision 1.9  2004/11/20 21:16:41  edwards
// Tried to simplify number of include lines.
//
// Revision 1.8  2004/10/29 13:36:13  bjoo
// Added generalised preconditioned nef operator, and a zolotarev facade to it. Can now do Chiu! Which is just as well because I am just about to go down with the flu
//
// Revision 1.7  2004/10/21 16:43:20  bjoo
// UNPRECONDITIONED_ZOLO_NEF
//
// Revision 1.6  2004/10/03 01:21:19  edwards
// Removed Dminus on array. Changed Dminus to also require N5 slice.
// Changed NEF linop to require an array of b5/c5. Changed NEF fermact to
// possibly except an array or scalar for b5/c5.
//
// Revision 1.5  2004/09/19 02:41:02  edwards
// Only indented. NOTE: this routine could soon (?) be merged
// as the general dwf_quarkProp4 routine. Needs, currents, etc.
//
// Revision 1.4  2004/09/16 16:40:53  kostas
// nef_quarkprop4 now compiles
//
// Revision 1.3  2004/09/16 15:20:55  kostas
// fixed up mres for NEF
//
// Revision 1.2  2004/09/03 14:24:36  kostas
// added mres measurement for NEF fermions
//
/*! \file
 * \brief Full quark propagator solver for domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */


#include "chromabase.h"
#include "actions/ferm/qprop/nef_quarkprop4_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"


namespace Chroma
{
  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to Domain Wall fermions (Array) fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param t_src    time slice of source ( Read )
   * \param j_decay  direction of decay ( Read )
   * \param invParam inverter params ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  template<typename T, typename P, typename Q, template<class,class,class> class C>
  void nef_quarkProp_a(LatticePropagator& q_sol, 
		       XMLWriter& xml_out,
		       const LatticePropagator& q_src,
		       int t_src, int j_decay,
		       const C<T,P,Q>& S_f,
		       Handle< FermState<T,P,Q> > state,
		       const GroupXML_t& invParam,
		       int& ncg_had)
  {
    START_CODE();

    push(xml_out, "DWF_QuarkProp4");

    ncg_had = 0;

    // Setup solver
    Handle< SystemSolverArray<LatticeFermion> > qpropT(S_f.qpropT(state,invParam));

    multi1d<LatticePropagator> prop5d(S_f.size()) ;
    LatticePropagator q_mp  ;

    multi1d<LatticeFermion> psi(S_f.size()) ;
    
    // This version loops over all color and spin indices
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	QDPIO::cout<<"nef_quarkProp:: doing color  : "<< color_source;
	QDPIO::cout<<" and spin : "<< spin_source<<endl  ;

	psi = zero ;  // note this is ``zero'' and not 0
	LatticeFermion tmp,tt ;
	tmp = zero ;
	PropToFerm(q_src, tmp, color_source, spin_source);
	   
	/* 
	 * Normalize the source in case it is really huge or small - 
	 * a trick to avoid overflows or underflows
	 */

	Real fact = 1.0;
	Real nrm = sqrt(norm2(tmp));
	if (toFloat(nrm) != 0.0)
	  fact /= nrm;

	// Rescale
	tmp *= fact;

	QDPIO::cout<<"Normalization Factor: "<< fact<<endl ;

	int N5(S_f.size());
	QDPIO::cout << "N5=" << N5 << endl;

	multi1d<LatticeFermion> chi(N5) ;
	chi = zero ;
	// Split the source to oposite walls according to chirality
	// and apply Dminus
	tt = chiralProjectPlus(tmp) ;
	S_f.Dminus(chi[0   ],tt,state,PLUS,0);
	tt = chiralProjectMinus(tmp) ; 
	S_f.Dminus(chi[N5-1],tt,state,PLUS,N5-1);

      

	// now we are ready invert
	// Compute the propagator for given source color/spin.	   
	{
	  SystemSolverResults_t result = (*qpropT)(psi,chi);
	  ncg_had += result.n_count;

	  push(xml_out,"Qprop");
	  write(xml_out, "color_source", color_source);
	  write(xml_out, "spin_source", spin_source);
	  write(xml_out, "n_count", result.n_count);
	  write(xml_out, "resid", result.resid);
	  pop(xml_out);
	}

	// Unnormalize the source following the inverse 
	// of the normalization above
	fact = Real(1) / fact;
	for(int i = 0; i < N5; ++i)
	  psi[i] *= fact; 

	/*
	 * Move the solution to the appropriate components
	 * of quark propagator.
	 */
	   
	//First the 5D quark propagator
	for(int s(0);s<S_f.size();s++)
	  FermToProp(psi[s], prop5d[s], color_source, spin_source);
	// Now get the 4D propagator too

	tmp = chiralProjectMinus(psi[0]) + chiralProjectPlus(psi[N5-1]) ;
	// move solution to the appropriate components of the 4d
	// quark propagator
	FermToProp(tmp, q_sol, color_source, spin_source);

	// move solution to the appropriate components of the 4d
	// midpoint quark propagator 
	tmp = chiralProjectPlus(psi[N5/2 - 1]) + chiralProjectMinus(psi[N5/2]) ;
	FermToProp(tmp, q_mp, color_source, spin_source);

	   
      }	/* end loop over spin_source */
    } /* end loop over color_source */

    LatticeComplex cfield ;

    /**
     // constuct the conserved axial current correlator
     nef_conserved_axial_ps_corr(cfield,state->getLinks(),prop5d,j_decay);
    **/
		       
    multi1d<DComplex> corr ;  
   
    SftMom trick(0,false,j_decay) ;
   
    corr = sumMulti(cfield, trick.getSet());
    // Length of lattice in time direction
    int length = trick.numSubsets();
    multi1d<Real> mesprop(length);
    for(int t(0);t<length; t++){
      int t_eff( (t - t_src + length) % length ) ;
      mesprop[t_eff] = real(corr[t]) ; 
    }

    push(xml_out, "time_direction");
    write(xml_out, "t_dir",j_decay);
    pop(xml_out);

    /**
       push(xml_out, "DWF_ConservedAxial");
       write(xml_out, "mesprop", mesprop); 
       pop(xml_out);

       // The local axial corruent pseudoscalar correlator
       int d(1<<j_decay);
       cfield = trace( adj(q_sol)*Gamma(d)*q_sol ) ;
       corr = sumMulti(cfield, trick.getSet()) ;
       for(int t(0);t<length; t++){
       int t_eff( (t - t_src + length) % length ) ;
       mesprop[t_eff] = -real(corr[t]) ; // sign fix
       }
       push(xml_out, "DWF_LocalAxial");
       write(xml_out, "mesprop", mesprop); 
       pop(xml_out);
    **/

    //Now the midpoint Pseudoscalar correlator
    multi1d<Double> tmp(length);
    tmp = sumMulti(localNorm2(q_mp), trick.getSet());
    for(int t(0);t<length; t++){
      int t_eff( (t - t_src + length) % length ) ;
      mesprop[t_eff] = tmp[t] ; 
    }
  
    push(xml_out, "DWF_MidPoint_Pseudo");
    write(xml_out, "mesprop", mesprop);
    pop(xml_out);


    tmp = sumMulti(localNorm2(q_sol), trick.getSet());
    for(int t(0);t<length; t++){
      int t_eff( (t - t_src + length) % length ) ;
      mesprop[t_eff] = tmp[t] ; // only need the zero momentum
    }
    push(xml_out, "DWF_Psuedo_Pseudo");
    write(xml_out, "mesprop", mesprop);
    pop(xml_out);

    pop(xml_out);   // DWF_QuarkProp

    /**
       check_nef_ward_identity(state->getLinks(),prop5d,q_src,
       q_sol,q_mp,S_f.quark_mass(),
       j_decay);
    **/

    END_CODE();
  }


  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to Domain Wall fermions (Array) fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param t_src    time slice of source ( Read )
   * \param j_decay  direction of decay ( Read )
   * \param invParam inverter params ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */

  void nef_quarkProp4(LatticePropagator& q_sol, 
		      XMLWriter& xml_out,
		      const LatticePropagator& q_src,
		      int t_src, int j_decay,
		      const UnprecDWFermActBaseArray<LatticeFermion,
		      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& S_f,
		      Handle< FermState<LatticeFermion, 
		      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
		      const GroupXML_t& invParam,
		      int& ncg_had)
  {
    nef_quarkProp_a<LatticeFermion,multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix>,
      UnprecDWFermActBaseArray>(
	q_sol, 
	xml_out, 
	q_src, 
	t_src, 
	j_decay, 
	S_f, 
	state,
	invParam,
	ncg_had);
  }

  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to Domain Wall fermions (Array) fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param t_src    time slice of source ( Read )
   * \param j_decay  direction of decay ( Read )
   * \param invParam inverter params ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */

  void nef_quarkProp4(LatticePropagator& q_sol, 
		      XMLWriter& xml_out,
		      const LatticePropagator& q_src,
		      int t_src, int j_decay,
		      const EvenOddPrecDWFermActBaseArray<LatticeFermion,
		      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& S_f,
		      Handle< FermState<LatticeFermion,
		      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
		      const GroupXML_t& invParam,
		      int& ncg_had)
  {
    nef_quarkProp_a<LatticeFermion,multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix>,
      EvenOddPrecDWFermActBaseArray>(
	q_sol, 
	xml_out, 
	q_src, 
	t_src, 
	j_decay, 
	S_f, 
	state,
	invParam,
	ncg_had);
  }

}
