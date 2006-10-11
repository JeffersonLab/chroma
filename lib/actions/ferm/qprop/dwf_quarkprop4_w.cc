// $Id: dwf_quarkprop4_w.cc,v 3.3 2006-10-11 15:42:26 edwards Exp $
// $Log: dwf_quarkprop4_w.cc,v $
// Revision 3.3  2006-10-11 15:42:26  edwards
// Removed use of numRetries !! This was a misguided attempt to get restarts
// into the code. Will try another way via the system solvers and the
// InvertParam xml group.
//
// Revision 3.2  2006/07/20 20:05:50  edwards
// Added norm.
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
// Revision 1.25  2005/02/21 19:28:59  edwards
// Changed initHeader's to be instead the appropriate default constructor
// for the ChiralParam, AnisoParam and QQQ structs. Removed all
// calls to initHeader. Changed quarkProp routines to instead take
// the appropriate SystemSolver instead of fermact. Added AnisoParam
// support for DWF and SSE DWF.
//
// Revision 1.24  2005/01/14 20:13:06  edwards
// Removed all using namespace QDP/Chroma from lib files. The library
// should now be 100% in the Chroma namespace. All mainprogs need a
// using namespace Chroma.
//
// Revision 1.23  2005/01/07 04:53:53  edwards
// Some debugging output added.
//
// Revision 1.22  2005/01/02 05:21:10  edwards
// Rearranged top-level fermion actions and linear operators.
// The deriv is pushed up much higher. This has the effect that
// now the "P" type (conjugate momentum type) is carried around
// in DiffFermAct4D<T,P> and below carry this "P" template param.
// All fermacts now have at least default deriv support (possibly
// Not Implemented error). Made qprop and qpropT now go through
// factory functions.
//
// Revision 1.21  2004/12/29 22:13:41  edwards
// Rearranged the top-level FermionAction structure. Moved up
// 5d-style actions to be a  FermAct5D and there is a FermAct4D.
// Consolidated quarkprop routines and fermact factories.
//
// Revision 1.20  2004/10/19 03:25:03  edwards
// Added SSE even/odd CG inverter hooks for now.
//
// Revision 1.19  2004/09/08 02:48:26  edwards
// Big switch-over of fermact IO. New fermact startup mechanism -
// now using Singleton Factory object. Moved  quarkprop4 to be
// a virtual func with top level FermionAction. Disconnected
// zolo4d and 5d fermacts temporarily. Removing usage of old
// fermact and subsequently md,gaugeact  IO mechanisms.
//
// Revision 1.18.2.2  2004/09/02 16:00:07  bjoo
// Trolled beneath the fermact mountains and changed the invocations
// of the inverters to conform to the new inverter interface (invParam -- vs RsdCG, MaxCG and friends) - also in the qprop_files too.
//
// I HAVE REMOVED ZOLOTAREV4D and ZOLOTAREV5D from the build as these
// will need extensive reworking to keep those params inside them
// (simply too many).  I will get them to conform later. Use the
// production branch if you need access to those
//
// I have added the various LOKI based files (Alexandrescu) to the
// Makefile.am so that things install correctly.
//
// I made the propagator code go through the factory invokation to
// do a propagator calculation with prec wilson fermions.
//
// Interestingly the linkage_hack appeared to be optimised away in its
// old form. I turned it into a function that returns the "foo" within,
// so it cannot be compiled away. Interestingly, it still doesn't actually
// need to be CALLED.
//
// Main achievements:
// 	Library compiles
// 	Propagator runs with suitable fermacts (I should check DWF etc)
//
// Revision 1.18.2.1  2004/09/01 15:13:10  edwards
// Start of major changes to support maps.
//
// Revision 1.18  2004/07/28 02:38:02  edwards
// Changed {START,END}_CODE("foo") to {START,END}_CODE().
//
// Revision 1.17  2004/02/23 03:05:11  edwards
// Pass in j_decay.
//
// Revision 1.16  2004/02/11 12:51:33  bjoo
// Stripped out Read() and Write()
//
// Revision 1.15  2004/02/10 22:59:46  kostas
// fixed a comment
//
// Revision 1.14  2004/02/06 16:47:41  kostas
// Ward identity check works. Fixed the local current sign
//
// Revision 1.13  2004/02/05 21:15:05  kostas
// fixed all bugs. Needs some clean up still
//
// Revision 1.12  2004/02/05 20:10:38  kostas
// Changed remaining getSubset to getSet
//
// Revision 1.11  2004/02/05 19:19:26  kostas
// few bugs fixed
//
// Revision 1.10  2004/02/03 20:04:53  edwards
// Changed phases.getSubset() to phases.getSet(). Removed passing j_decay
// into curcor2
//
// Revision 1.9  2004/01/30 21:35:49  kostas
// added uprec_dwf support
//
// Revision 1.8  2004/01/29 19:48:26  kostas
// added the t_src as input parameter
//
// Revision 1.7  2004/01/29 17:37:33  edwards
// Implemented some flop optimizations.
//
// Revision 1.6  2004/01/29 16:18:52  kostas
// Removed the if 1 directive
//
// Revision 1.5  2004/01/29 16:17:50  kostas
// Fixed an ax-ps correlation function bug
//
// Revision 1.4  2004/01/29 14:59:46  kostas
// Fixed the bugs. Code compiles now
//
/*! \file
 * \brief Full quark propagator solver for domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */


#include "chromabase.h"
#include "actions/ferm/qprop/dwf_quarkprop4_w.h"
#include "actions/ferm/linop/dwffld_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

namespace Chroma
{

  void dwf_conserved_axial_ps_corr(LatticeComplex& corr,
				   const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& p5d, 
				   const int mu) ;

  void check_dwf_ward_identity(const multi1d<LatticeColorMatrix>& u,
			       const multi1d<LatticePropagator>& p5d,
			       const LatticePropagator& src,
			       const LatticePropagator& q_q,
			       const LatticePropagator& q_mp_q,
			       const Real& m_q,
			       int j_decay)
  {
    QDPIO::cout<<"check_dwf_ward_identity: Checking the chiral Ward Identity...";
    QDPIO::cout<<endl ;

    LatticeComplex divA = zero;
    for(int mu(0);mu<Nd;mu++){
      LatticeComplex tt ;
      dwf_conserved_axial_ps_corr(tt,u,p5d,mu);
      divA += tt - shift(tt,BACKWARD,mu) ; 
    }

    LatticeComplex mpps_ps = localNorm2(q_mp_q) ;
    LatticeComplex ps_ps = localNorm2(q_q);
    LatticeComplex q_bar_q = localInnerProduct(src,q_q);
    LatticeComplex diff = divA - 2.0 * m_q * ps_ps - 2.0*mpps_ps + 2.0*q_bar_q;

    SftMom trick(0,false,j_decay) ;
    multi1d<Double> corr = sumMulti(localNorm2(diff), trick.getSet());
    QDPIO::cout<<"check_dwf_ward_identity: ";
    QDPIO::cout<<"Ward Identity violation per timeslice: "<<endl;
    for(int t(0);t<corr.size(); t++){
      QDPIO::cout<<"        "<<t<<"                "<< sqrt(corr[t])<<endl ;
    }

    QDPIO::cout<<"check_dwf_ward_identity: Ward Identity violation: ";
    QDPIO::cout<<sqrt(norm2(diff))<<endl ;

    QDPIO::cout<<"check_dwf_ward_identity: |divA|^2    : "<<norm2(divA)<<endl;
    QDPIO::cout<<"check_dwf_ward_identity: |ps_ps|^2   : "<<norm2(ps_ps)<<endl;
    QDPIO::cout<<"check_dwf_ward_identity: |mpps_ps|^2 : "<<norm2(mpps_ps)<<endl;
    QDPIO::cout<<"check_dwf_ward_identity: |q_bar_q|^2 : "<<norm2(q_bar_q)<<endl;
    Double gmor( sqrt(norm2(sum(m_q*ps_ps + mpps_ps - q_bar_q))) );
    QDPIO::cout<<"check_dwf_ward_identity: GMOR        : "<<gmor<<endl;
  
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
   * \param qpropT   5D inverter ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  void dwf_quarkProp4(LatticePropagator& q_sol, 
		      XMLWriter& xml_out,
		      const LatticePropagator& q_src,
		      int t_src, int j_decay,
		      Handle< SystemSolverArray<LatticeFermion> > qpropT,
		      Handle< FermState<LatticeFermion,
		      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
		      const Real& m_q,
		      int& ncg_had)
  {
    START_CODE();

    QDPIO::cout << "entering DWF_QuarkProp4" << endl;

    push(xml_out, "DWF_QuarkProp4");

    ncg_had = 0;

    int N5 = qpropT->size();
    multi1d<LatticePropagator> prop5d(N5);
    LatticePropagator q_mp;

    multi1d<LatticeFermion> psi(N5);
    
    // This version loops over all color and spin indices
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = 0; spin_source < Ns; ++spin_source)
      {
	QDPIO::cout<<"dwf_quarkProp:: doing color  : "<< color_source;
	QDPIO::cout<<" and spin : "<< spin_source<<endl  ;

	psi = zero ;  // note this is ``zero'' and not 0
	LatticeFermion tmp = zero;
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
	  

	//QDPIO::cout<<"Normalization Factor: "<< fact<<endl ;

	multi1d<LatticeFermion> chi(N5) ;
	chi = zero ;
	// Split the source to oposite walls according to chirality
	chi[0   ] = chiralProjectPlus(tmp) ;
	chi[N5-1] = chiralProjectMinus(tmp) ; 

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
	for(int s(0);s<N5;s++)
	  FermToProp(psi[s], prop5d[s], color_source, spin_source);
	// Now get the 4D propagator too

	tmp = chiralProjectMinus(psi[0]) + chiralProjectPlus(psi[N5-1]);

	// move solution to the appropriate components of the 4d
	// quark propagator
	FermToProp(tmp, q_sol, color_source, spin_source);

	// move solution to the appropriate components of the 4d
	// midpoint quark propagator 
	tmp = chiralProjectPlus(psi[N5/2 - 1]) + chiralProjectMinus(psi[N5/2]) ;
	FermToProp(tmp, q_mp, color_source, spin_source);
	
      }	/* end loop over spin_source */
    } /* end loop over color_source */

    // Construct the norm of the chiral field throughout the bulk
    {
      multi1d<Double> bulk_norm2(N5/2);
      bulk_norm2 = zero;

      for(int s=0; s < bulk_norm2.size(); ++s)
      {
	bulk_norm2[s] = norm2(chiralProjectMinus(prop5d[s]) + chiralProjectPlus(prop5d[N5-1-s]));
      }

      write(xml_out, "bulk_norm2", bulk_norm2);
    }

    // constuct the conserved axial current correlator
    LatticeComplex cfield ;
    dwf_conserved_axial_ps_corr(cfield,state->getLinks(),prop5d,j_decay);
			       
	
    SftMom trick(0,false,j_decay) ;
    int length = trick.numSubsets();   // Length of lattice in time direction
   
    multi1d<Real> mesprop(length);
    {
      multi1d<DComplex> corr = sumMulti(cfield, trick.getSet());
      for(int t(0);t<length; t++){
	int t_eff( (t - t_src + length) % length ) ;
	mesprop[t_eff] = real(corr[t]) ; 
      }
    }

    push(xml_out, "time_direction");
    write(xml_out, "t_dir",j_decay);
    pop(xml_out);
    push(xml_out, "DWF_ConservedAxial");
    write(xml_out, "mesprop", mesprop); 
    pop(xml_out);

    // The local axial corruent pseudoscalar correlator
    int d(1<<j_decay);
    cfield = localInnerProduct(q_sol,Gamma(d)*q_sol);
    {
      multi1d<DComplex> corr = sumMulti(cfield, trick.getSet()) ;
      for(int t(0);t<length; t++){
	int t_eff( (t - t_src + length) % length ) ;
	mesprop[t_eff] = -real(corr[t]) ; // sign fix
      }
    }
    push(xml_out, "DWF_LocalAxial");
    write(xml_out, "mesprop", mesprop); 
    pop(xml_out);

    //Now the midpoint Pseudoscalar correlator
    multi1d<Double> tmp = sumMulti(localNorm2(q_mp), trick.getSet());
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

    check_dwf_ward_identity(state->getLinks(),prop5d,q_src,
			    q_sol,q_mp,m_q,
			    j_decay);


    QDPIO::cout << "exiting DWF_QuarkProp4" << endl;

    END_CODE();
  }


  /*!
    Corelation function:
    \f[
    C(t) = \sum_x \sum_s \left[\bar{\Psi}(x+\hat{\mu},t,s) U^{\dagger}_{\mu}(x)
    \frac{1+\gamma_{\mu}}{2} 
    \Psi (x,t,s) - 
    \bar{\Psi}(x,t,s) U_{\mu}(x)
    \frac{1-\gamma_{\mu}}{2} 
    \Psi (x+\hat{\mu},t,s)\right] 
    \bar{q}(0,0)\gamma_5 q(0,0) sign(s - \frac{L_s - 1}{2})
    \f]
  **/
  void dwf_conserved_axial_ps_corr(LatticeComplex& corr,
				   const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& p5d, 
				   const int mu)
  {
    // gamma mapping G_0 --> Gamma(1)
    //               G_1 --> Gamma(2)
    //               G_2 --> Gamma(4)
    //               G_3 --> Gamma(8)
    int d(1<<mu);
    int g5(Ns*Ns - 1) ;
    int N5(p5d.size());

    corr = zero;

    multi1d<LatticePropagator> us_p5d(N5) ;
    for(int s(0); s<N5;s++)
      us_p5d[s] = u[mu]*shift(p5d[s],FORWARD,mu) ;
  
    for(int s(0); s<N5;s++){
      // first the 1-gamma_mu term 
      LatticeComplex  C;
      C=0.5*(  trace(adj(   p5d[N5-1-s])*Gamma(g5)*Gamma(d)*us_p5d[s]) -
	       trace(adj(   p5d[N5-1-s])*Gamma(g5)         *us_p5d[s]) +
	       //now the 1+gamma_mu term
	       trace(adj(us_p5d[N5-1-s])*Gamma(g5)*Gamma(d)*   p5d[s]) +
	       trace(adj(us_p5d[N5-1-s])*Gamma(g5)         *   p5d[s])   );

      if(s<N5/2) 
	corr -= C ;
      else
	corr += C ;
    }  
  
  }

}
