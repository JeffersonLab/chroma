// $Id: nef_quarkprop4_w.cc,v 1.7 2004-10-21 16:43:20 bjoo Exp $
// $Log: nef_quarkprop4_w.cc,v $
// Revision 1.7  2004-10-21 16:43:20  bjoo
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
#include "fermact.h"
#include "actions/ferm/qprop/dwf_quarkprop4_w.h"
#include "actions/ferm/linop/dwffld_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

using namespace QDP;

//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param t_src    time slice of source ( Read )
 * \param j_decay  direction of decay ( Read )
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */
template<class C>
static 
void nef_quarkProp4_a(LatticePropagator& q_sol, 
		      XMLWriter& xml_out,
		      const LatticePropagator& q_src,
		      int t_src, int j_decay,
		      const C& S_f,
		      Handle<const ConnectState> state,
		      const InvertParam_t& invParam,
		      int& ncg_had)
{
  START_CODE();

  push(xml_out, "DWF_QuarkProp4");

  ncg_had = 0;

  multi1d<LatticePropagator> prop5d(S_f.size()) ;
  LatticePropagator q_mp  ;

  multi1d<LatticeFermion> psi(S_f.size()) ;
    
  // This version loops over all color and spin indices
  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    for(int spin_source = 0; spin_source < Ns; ++spin_source)
    {
      QDPIO::cout<<"nef_quarkProp4:: doing color  : "<< color_source;
      QDPIO::cout<<" and spin : "<< spin_source<<endl  ;

      psi = zero ;  // note this is ``zero'' and not 0
      LatticeFermion tmp,tt ;
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
      // and apply Dminus
      tt = chiralProjectPlus(tmp) ;
      S_f.Dminus(chi[0   ],tt,state,PLUS,0);
      tt = chiralProjectMinus(tmp) ; 
      S_f.Dminus(chi[N5-1],tt,state,PLUS,N5-1);

      

      // now we are ready invert
	   

      int n_count;
      // Compute the propagator for given source color/spin.	   
      S_f.qpropT(psi, state, chi, invParam, n_count);
      ncg_had += n_count;

      push(xml_out,"Qprop");
      write(xml_out, "RsdCG", invParam.RsdCG);
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

  pop(xml_out);   // DWF_QuarkProp4

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
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void UnprecNEFFermActArray::dwf_quarkProp4(LatticePropagator& q_sol, 
					   XMLWriter& xml_out,
					   const LatticePropagator& q_src,
					   int t_src, int j_decay,
					   Handle<const ConnectState> state,
					   const InvertParam_t& invParam,
					   int& ncg_had)
{
  nef_quarkProp4_a<UnprecNEFFermActArray>(q_sol, 
					  xml_out, 
					  q_src, 
					  t_src, 
					  j_decay, 
					  *this, 
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
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void UnprecZoloNEFFermActArray::dwf_quarkProp4(LatticePropagator& q_sol, 
					   XMLWriter& xml_out,
					   const LatticePropagator& q_src,
					   int t_src, int j_decay,
					   Handle<const ConnectState> state,
					   const InvertParam_t& invParam,
					   int& ncg_had)
{
  nef_quarkProp4_a<UnprecZoloNEFFermActArray>(q_sol, 
					      xml_out, 
					      q_src, 
					      t_src, 
					      j_decay, 
					      *this, 
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
 * \param invType  inverter type ( Read )
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */
void 
EvenOddPrecNEFFermActArray::dwf_quarkProp4(LatticePropagator& q_sol, 
					   XMLWriter& xml_out,
					   const LatticePropagator& q_src,
					   int t_src, int j_decay,
					   Handle<const ConnectState> state,
					   const InvertParam_t& invParam,
					   int& ncg_had)
{
  nef_quarkProp4_a<EvenOddPrecNEFFermActArray>(q_sol, 
					       xml_out, 
					       q_src, 
					       t_src, 
					       j_decay, 
					       *this, 
					       state, 
					       invParam, 
					       ncg_had);
}

