// $Id: dwf_quarkprop4_w.cc,v 1.11 2004-02-05 19:19:26 kostas Exp $
// $Log: dwf_quarkprop4_w.cc,v $
// Revision 1.11  2004-02-05 19:19:26  kostas
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

  LatticeComplex mpps_ps, ps_ps, contact ;
  mpps_ps = localNorm2(q_mp_q) ;
  ps_ps = localNorm2(q_q) ;

  //contact = 

  LatticeComplex diff ;

  diff = divA - 2.0 * m_q * ps_ps - 2.0*mpps_ps ;

  QDPIO::cout<<"check_dwf_ward_identity: Ward Identity violation: ";
  QDPIO::cout<<norm2(diff)<<endl ;

  QDPIO::cout<<"check_dwf_ward_identity: |divA|^2    : "<< norm2(divA)<<endl;
  QDPIO::cout<<"check_dwf_ward_identity: |ps_ps|^2   : "<< norm2(ps_ps)<<endl;
  QDPIO::cout<<"check_dwf_ward_identity: |mpps_ps|^2 : "<< norm2(mpps_ps)<<endl;

  QDPIO::cout<<"check_dwf_ward_identity: Sanity check: "<<endl;
  
  for(int s(0);s<p5d.size();s++){
    QDPIO::cout<<"                                     "<<norm2(p5d[s])<<endl ;
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

  ncg_had = 0;

  multi1d<LatticePropagator> prop5d(S_f.size()) ;
  LatticePropagator q_mp  ;

  multi1d<LatticeFermion> psi(S_f.size()) ;
 
  psi = zero;  // note this is ``zero'' and not 0
   
  // This version loops over all color and spin indices
  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    for(int spin_source = 0; spin_source < Ns; ++spin_source)
    {
      psi = zero ;
      LatticeFermion tmp ;
      tmp = zero ;
      PropToFerm(q_src, tmp, color_source, spin_source);

      // Use the last initial guess as the current initial guess
	   
      /* 
       * Normalize the source in case it is really huge or small - 
       * a trick to avoid overflows or underflows
       */


      Real fact = Real(1) / sqrt(norm2(tmp));
      tmp *= fact ;

QDPIO::cout<<"Normalization Factor: "<< fact<<endl ;

      int N5(S_f.size());
      multi1d<LatticeFermion> chi(N5) ;
      chi = zero ;
      // Split the source to oposite walls according to chirality
      chi[0   ] = chiralProjectPlus(tmp) ;
      chi[N5-1] = chiralProjectMinus(tmp) ; 

QDPIO::cout<<"|chi|^2 : "<< norm2(chi)<<endl  ;

      // DwfFld(chi, tmp, PLUS); // It does not do what I thought it does

      // now we are ready invert
	   
      // Compute the propagator for given source color/spin.
      int n_count;
	   
      S_f.qpropT(psi, state, chi, invType, RsdCG, MaxCG, n_count);
      ncg_had += n_count;
QDPIO::cout<<"psi.size() : "<< psi.size()<<endl  ;
QDPIO::cout<<"color  : "<< color_source;
QDPIO::cout<<" spin : "<< spin_source<<endl  ;
 for(int s(0);s<N5;s++)
   QDPIO::cout<<"|psi["<<s<<"]|^2 : "<< norm2(psi[s])<<endl  ;
      push(xml_out,"Qprop");
      write(xml_out, "RsdCG", RsdCG);
      Write(xml_out, n_count);
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
  dwf_conserved_axial_ps_corr(cfield,state->getLinks(),prop5d,3);
			       
  multi1d<DComplex> corr ;  
   
  SftMom trick(0,false,3) ;
   
  corr = sumMulti(cfield, trick.getSet());
  // Length of lattice in decay direction
  int length = trick.numSubsets();
  multi1d<Real> mesprop(length);
  //multi1d<Complex> axps(length);
  for(int t(0);t<length; t++){
    int t_eff( (t - t_src + length) % length ) ;
    mesprop[t_eff] = real(corr[t]) ; 
  }

  //push(xml_out, "t_dir");
  //Write(xml_out, 3);
  //pop(xml_out);
  push(xml_out, "DWF_ConservedAxial");
  //write(xml_out, "corr", corr);
  Write(xml_out, mesprop); 
  pop(xml_out);

  // The local axial corruent pseudoscalar correlator
  cfield = trace( adj(q_sol)*Gamma(8)*q_sol ) ;
  corr = sumMulti(cfield, trick.getSubset()) ;
  for(int t(0);t<length; t++){
    int t_eff( (t - t_src + length) % length ) ;
    mesprop[t_eff] = real(corr[t]) ; 
  }
  push(xml_out, "DWF_LocalAxial");
  Write(xml_out, mesprop); 
  pop(xml_out);

  //Now the midpoint Pseudoscalar correlator
  multi1d<Double> tmp(length);
  tmp = sumMulti(localNorm2(q_mp), trick.getSet());
  for(int t(0);t<length; t++){
    int t_eff( (t - t_src + length) % length ) ;
    mesprop[t_eff] = tmp[t] ; // only need the zero momentum
  }
  
  push(xml_out, "DWF_MidPoint_Pseudo");
  Write(xml_out, mesprop);
  pop(xml_out);


  tmp = sumMulti(localNorm2(q_sol), trick.getSubset());
  for(int t(0);t<length; t++){
    int t_eff( (t - t_src + length) % length ) ;
    mesprop[t_eff] = tmp[t] ; // only need the zero momentum
  }
  push(xml_out, "DWF_Psuedo_Pseudo");
  Write(xml_out, mesprop);
  pop(xml_out);

  pop(xml_out);   // DWF_QuarkProp4

  check_dwf_ward_identity(state->getLinks(),prop5d,q_sol,q_mp,S_f.quark_mass());


  END_CODE("dwf_quarkProp4");
}


/*!
  Corelation function:
\f[
  C(t) = \sum_x \sum_s \left[\bar{\Psi}(x+\hat{\mu},t,s)\frac{1+\gamma_{\mu}}{2} 
                                  \Psi (x,t,s) - 
                             \bar{\Psi}(x,t,s)\frac{1-\gamma_{\mu}}{2} 
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
  int d(1);
  for(int c(0);c<mu;c++)
    d *= 2 ;

  int g5(Ns*Ns - 1) ;

  LatticeComplex  C(0.0) ;
 
  int N5(p5d.size());
  for(int s(0); s<N5;s++){
    // first the 1-gamma_mu term 
    C=.5*(trace(adj(p5d[N5-1-s])*Gamma(g5)*Gamma(d)*u[mu]*shift(p5d[s],FORWARD,mu))-
	  trace(adj(p5d[N5-1-s])*Gamma(g5)         *u[mu]*shift(p5d[s],FORWARD,mu))+
	  //now the 1+gamma_mu term
	  trace(shift(adj(p5d[N5-1-s]),FORWARD,mu)*Gamma(g5)*Gamma(d)*adj(u[mu])*p5d[s])+
	  trace(shift(adj(p5d[N5-1-s]),FORWARD,mu)*Gamma(g5)         *adj(u[mu])*p5d[s])
	  );
	  
      
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

