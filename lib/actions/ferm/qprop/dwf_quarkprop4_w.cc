// $Id: dwf_quarkprop4_w.cc,v 1.16 2004-02-11 12:51:33 bjoo Exp $
// $Log: dwf_quarkprop4_w.cc,v $
// Revision 1.16  2004-02-11 12:51:33  bjoo
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
#include "fermact.h"
#include "actions/ferm/qprop/dwf_quarkprop4_w.h"
#include "actions/ferm/linop/dwffld_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

using namespace QDP;

void dwf_conserved_axial_ps_corr(LatticeComplex& corr,
				 const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& p5d, 
				 const int mu) ;

void check_dwf_ward_identity(const multi1d<LatticeColorMatrix>& u,
			     const multi1d<LatticePropagator>& p5d,
			     const LatticePropagator& src,
			     const LatticePropagator& q_q,
			     const LatticePropagator& q_mp_q,
			     const Real& m_q)
{
  QDPIO::cout<<"check_dwf_ward_identity: Checking the chiral Ward Identity...";
  QDPIO::cout<<endl ;

  LatticeComplex divA ;
  divA = 0.0 ;
  for(int mu(0);mu<Nd;mu++){
    LatticeComplex tt ;
    dwf_conserved_axial_ps_corr(tt,u,p5d,mu);
    divA += tt - shift(tt,BACKWARD,mu) ; 
  }

  LatticeComplex mpps_ps, ps_ps, q_bar_q ;
  mpps_ps = localNorm2(q_mp_q) ;
  ps_ps = localNorm2(q_q) ;

  q_bar_q = trace(adj(src)*q_q) ;

  LatticeComplex diff ;

  diff = divA - 2.0 * m_q * ps_ps - 2.0*mpps_ps + 2.0*q_bar_q;

  multi1d<Double> corr ;     
  SftMom trick(0,false,3) ;
  corr = sumMulti(localNorm2(diff), trick.getSet());
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
  diff = m_q*ps_ps + mpps_ps - q_bar_q ;
  Double gmor( sqrt(norm2(sum(diff))) ) ;
  QDPIO::cout<<"check_dwf_ward_identity: GMOR        : "<<gmor<<endl;
  
}


//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */
template<typename T, template<class> class C>
static 
void dwf_quarkProp4_a(LatticePropagator& q_sol, 
		      XMLWriter& xml_out,
		      const LatticePropagator& q_src,
		      const int t_src,
		      const C<T>& S_f,
		      Handle<const ConnectState> state,
		      enum InvType invType,
		      const Real& RsdCG, 
		      int MaxCG, int& ncg_had)
{
  START_CODE("dwf_quarkProp4");

  push(xml_out, "DWF_QuarkProp4");

  int t_dir = Nd-1 ; // Hard code Nd-1 to be the time direction 

  ncg_had = 0;

  multi1d<LatticePropagator> prop5d(S_f.size()) ;
  LatticePropagator q_mp  ;

  multi1d<LatticeFermion> psi(S_f.size()) ;
    
  // This version loops over all color and spin indices
  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    for(int spin_source = 0; spin_source < Ns; ++spin_source)
    {
      QDPIO::cout<<"dwf_quarkProp4:: doing color  : "<< color_source;
      QDPIO::cout<<" and spin : "<< spin_source<<endl  ;

      psi = zero ;  // note this is ``zero'' and not 0
      LatticeFermion tmp ;
      tmp = zero ;
      PropToFerm(q_src, tmp, color_source, spin_source);
	   
      /* 
       * Normalize the source in case it is really huge or small - 
       * a trick to avoid overflows or underflows
       */


      Real fact = Real(1) / sqrt(norm2(tmp));
      tmp *= fact ;

      //QDPIO::cout<<"Normalization Factor: "<< fact<<endl ;

      int N5(S_f.size());
      multi1d<LatticeFermion> chi(N5) ;
      chi = zero ;
      // Split the source to oposite walls according to chirality
      chi[0   ] = chiralProjectPlus(tmp) ;
      chi[N5-1] = chiralProjectMinus(tmp) ; 


      // now we are ready invert
	   

      int n_count;
      // Compute the propagator for given source color/spin.	   
      S_f.qpropT(psi, state, chi, invType, RsdCG, MaxCG, n_count);
      ncg_had += n_count;

      push(xml_out,"Qprop");
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", n_count);
      pop(xml_out);

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

  // constuct the conserved axial current correlator
  LatticeComplex cfield ;
  dwf_conserved_axial_ps_corr(cfield,state->getLinks(),prop5d,t_dir);
			       
  multi1d<DComplex> corr ;  
   
  SftMom trick(0,false,t_dir) ;
   
  corr = sumMulti(cfield, trick.getSet());
  // Length of lattice in time direction
  int length = trick.numSubsets();
  multi1d<Real> mesprop(length);
  for(int t(0);t<length; t++){
    int t_eff( (t - t_src + length) % length ) ;
    mesprop[t_eff] = real(corr[t]) ; 
  }

  push(xml_out, "time_direction");
  write(xml_out, "t_dir",t_dir);
  pop(xml_out);
  push(xml_out, "DWF_ConservedAxial");
  write(xml_out, "mesprop", mesprop); 
  pop(xml_out);

  // The local axial corruent pseudoscalar correlator
  int d(1<<t_dir);
  cfield = trace( adj(q_sol)*Gamma(d)*q_sol ) ;
  corr = sumMulti(cfield, trick.getSet()) ;
  for(int t(0);t<length; t++){
    int t_eff( (t - t_src + length) % length ) ;
    mesprop[t_eff] = -real(corr[t]) ; // sign fix
  }
  push(xml_out, "DWF_LocalAxial");
  write(xml_out, "mesprop", mesprop); 
  pop(xml_out);

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

  pop(xml_out);   // DWF_QuarkProp4

  check_dwf_ward_identity(state->getLinks(),prop5d,q_src,
			  q_sol,q_mp,S_f.quark_mass());


  END_CODE("dwf_quarkProp4");
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

  LatticeComplex  C ;
  corr = 0.0 ;

  int N5(p5d.size());

  multi1d<LatticePropagator> us_p5d(N5) ;
  for(int s(0); s<N5;s++)
    us_p5d[s] = u[mu]*shift(p5d[s],FORWARD,mu) ;
  
  for(int s(0); s<N5;s++){
    // first the 1-gamma_mu term 
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


//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    const int t_src,
		    const EvenOddPrecDWFermActBaseArray<LatticeFermion>& S_f,
		    Handle<const ConnectState> state,
		    enum InvType invType,
		    const Real& RsdCG, 
		    int MaxCG, int& ncg_had)
{
  dwf_quarkProp4_a(q_sol, xml_out, q_src, t_src, S_f, state, invType, RsdCG, MaxCG, ncg_had);
}

//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void dwf_quarkProp4(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    const int t_src,
		    const UnprecDWFermActBaseArray<LatticeFermion>& S_f,
		    Handle<const ConnectState> state,
		    enum InvType invType,
		    const Real& RsdCG, 
		    int MaxCG, int& ncg_had)
{
  dwf_quarkProp4_a(q_sol, xml_out, q_src, t_src, S_f, state, invType, RsdCG, MaxCG, ncg_had);
}

