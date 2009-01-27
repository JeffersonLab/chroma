// $Id: Ql_3pt_w.cc,v 1.3 2009-01-27 15:39:32 caubin Exp $ 
/*! \file
 *  \brief Static-Light 3pt function
 */

#include "Ql_3pt_w.h"
#include "barQll_w.h"

namespace Chroma {

//! Heavy-light meson 2-pt function
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * The heavy quark is inserted in the infinitely heavy quark limit
 * by a Wilson-Line without spin indices.
 * We are effectively propagating a spin-1/2 light quark
 * This generates with two quark propagators all three-point
 * functions, with all 16 gamma matrix insertions, of the form:
 *
 * Mat1 = Sc(z-y) [Gamma] Su(y-z) gamma5
 * Mat2 = Sb(x-z) gamma5 Ss(z-x) [Gamma]
 * 
 * And we output the 16 gamma's.
 *                                                                                                       
 * \param u                  gauge field (Read)                                                          
 * \param quark_propagator1   quark propagator1 ( Read )                                                 
 * \param quark_propagator2   quark propagator2 ( Read )                                                 
 * \param heavy_quark_propagator1   heavy quark propagator1 ( Read; Opt, if not static  )                                                 
 * \param heavy_quark_propagator2   heavy quark propagator2 ( Read; Opt, if not static )                                                 
 * \param src_coord          cartesian coordinates of the source ( Read )                                
 * \param snk_coord          cartesian coordinates of the sink ( Read )                       

 * \param heavy_src          cartesian coordinates of the heavy quark if not static ( Read )                                  
 *                           Not needed if both are static or both are not static.
 * \param phases             object holds list of momenta and Fourier phases ( Read )                    
 * \param xml                xml file object ( Read )                                                    
 * \param xml_group          group name for xml data ( Read )                                            
 *                                                                                                       
 */
  
void QlQl(const multi1d<LatticeColorMatrix>& u,
	    const LatticePropagator& quark_propagator1,
	    const LatticePropagator& quark_propagator2,
	    const multi1d<int>& src_coord,
	    const multi1d<int>& snk_coord,
	    const SftMom& phases,
	    XMLWriter& xml,
	    const string& xml_group)
{
    START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  LatticeColorMatrix Qprop1; // b quark
  LatticeColorMatrix Qprop2; // c quark

  HeavyQuarkProp(Qprop1,u,src_coord,length);
  HeavyQuarkPropBack(Qprop2,u,snk_coord,length);

  // S_proj_unpol = (1/2)(1 + gamma_4)
  // This is the heavy quark dirac structure, but I am using the notation
  // that is currently used in Chroma!!!!!!!!

  SpinMatrix g_one = 1.0;
  SpinMatrix S_proj_unpol = 0.5 * (g_one + (g_one * Gamma(8)));
  multi1d<DComplex> Hq;
  
  for(int gind=0; gind<16;++gind){// Loop over gamma matrices...
    LatticePropagator Mat1 , Mat2;
    LatticeComplex Hq_prop , Hq_prop_m ;

    Mat1 = quark_propagator1 * adj(Qprop1*S_proj_unpol) * Gamma(15) * Gamma(gind) ;
    Mat2 = Qprop2*S_proj_unpol * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
    
    Hq_prop = trace(Mat1)*trace(Mat2);
    Hq_prop_m = traceColor(traceSpin(Mat1)*traceSpin(Mat2));
    
    // Project onto zero momentum
    multi2d<DComplex> hsumHq, hsumHq_m ;
    hsumHq = phases.sft(Hq_prop);
    hsumHq_m = phases.sft(Hq_prop_m);
    multi2d<DComplex> HQprop(num_mom,length) ;
    multi2d<DComplex> HQprop_m(num_mom,length) ;
    
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - src_coord[Nd-1] + length) % length;
	  HQprop[sink_mom_num][t_eff] = hsumHq[sink_mom_num][t];
	  HQprop_m[sink_mom_num][t_eff] = hsumHq_m[sink_mom_num][t];
	}
    
    //XMLWriter xml_bar(xml);
    push(xml, xml_group);
    write(xml, "gamma_value", gind);
    write(xml, "ThreePointColorUnmixed", HQprop[0]);
    //    pop(xml);
    //    push(xml, xml_group);
    write(xml, "ThreePointColorMixed", HQprop_m[0]);
    pop(xml);
    
    END_CODE();
  }// gind loop
}// QlQl routine
  
  
void QlQl(const multi1d<LatticeColorMatrix>& u,
	    const LatticePropagator& quark_propagator1,
	    const LatticePropagator& quark_propagator2,
		const LatticePropagator& heavy_quark_propagator,
	    const multi1d<int>& src_coord,
	    const multi1d<int>& snk_coord,
	    const multi1d<int>& heavy_src, 
	    const SftMom& phases,
	    XMLWriter& xml,
	    const string& xml_group)
{
    START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  LatticeColorMatrix Qprop; // b quark

  if (heavy_src==src_coord)
	HeavyQuarkProp(Qprop,u,snk_coord,length);
  else 
	HeavyQuarkProp(Qprop,u,src_coord,length);

  // S_proj_unpol = (1/2)(1 + gamma_4)
  // This is the heavy quark dirac structure, but I am using the notation
  // that is currently used in Chroma!!!!!!!!

  SpinMatrix g_one = 1.0;
  SpinMatrix S_proj_unpol = 0.5 * (g_one + (g_one * Gamma(8)));
  multi1d<DComplex> Hq;
  
  for(int gind=0; gind<16;++gind){// Loop over gamma matrices...
    LatticePropagator Mat1 , Mat2;
    LatticeComplex Hq_prop , Hq_prop_m ;

	if (heavy_src==src_coord){
	  Mat1 = quark_propagator1 * adj(heavy_quark_propagator) * Gamma(15) * Gamma(gind) ;
	  Mat2 = Qprop*S_proj_unpol * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
	}
	else{
	  Mat1 = quark_propagator1 * adj(Qprop*S_proj_unpol) * Gamma(15) * Gamma(gind) ;
	  Mat2 = heavy_quark_propagator * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
    }

    Hq_prop = trace(Mat1)*trace(Mat2);
    Hq_prop_m = traceColor(traceSpin(Mat1)*traceSpin(Mat2));
    
    // Project onto zero momentum
    multi2d<DComplex> hsumHq, hsumHq_m ;
    hsumHq = phases.sft(Hq_prop);
    hsumHq_m = phases.sft(Hq_prop_m);
    multi2d<DComplex> HQprop(num_mom,length) ;
    multi2d<DComplex> HQprop_m(num_mom,length) ;
    
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - src_coord[Nd-1] + length) % length;
	  HQprop[sink_mom_num][t_eff]  = hsumHq[sink_mom_num][t];
	  HQprop_m[sink_mom_num][t_eff]  = hsumHq_m[sink_mom_num][t];
	}
    
    push(xml, xml_group);
    write(xml, "gamma_value", gind);
    write(xml, "ThreePointColorUnmixed", HQprop[0]);
    write(xml, "ThreePointColorMixed", HQprop_m[0]);
    pop(xml);
    
    END_CODE();
  }// gind loop
}// QlQl routine
  
  
void QlQl(const multi1d<LatticeColorMatrix>& u,
	    const LatticePropagator& quark_propagator1,
	    const LatticePropagator& quark_propagator2,
	    const LatticePropagator& heavy_quark_propagator1,
	    const LatticePropagator& heavy_quark_propagator2,
	    const multi1d<int>& src_coord,
	    const multi1d<int>& snk_coord,
	    const SftMom& phases,
	    XMLWriter& xml,
	    const string& xml_group)
{
    START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  // S_proj_unpol = (1/2)(1 + gamma_4)
  // This is the heavy quark dirac structure, but I am using the notation
  // that is currently used in Chroma!!!!!!!!

  multi1d<DComplex> Hq;
  
  for(int gind=0; gind<16;++gind){// Loop over gamma matrices...
    LatticePropagator Mat1 , Mat2;
    LatticeComplex Hq_prop , Hq_prop_m ;

    Mat1 = quark_propagator1 * adj(heavy_quark_propagator1) * Gamma(15) * Gamma(gind) ;
    Mat2 = heavy_quark_propagator2 * adj(quark_propagator2) * Gamma(15) * Gamma(gind);
    
    Hq_prop = trace(Mat1)*trace(Mat2);
    Hq_prop_m = traceColor(traceSpin(Mat1)*traceSpin(Mat2));
    
    // Project onto zero momentum
    multi2d<DComplex> hsumHq, hsumHq_m ;
    hsumHq = phases.sft(Hq_prop);
    hsumHq_m = phases.sft(Hq_prop_m);
    multi2d<DComplex> HQprop(num_mom,length) ;
    multi2d<DComplex> HQprop_m(num_mom,length) ;
    
    for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
      for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - src_coord[Nd-1] + length) % length;
	  HQprop[sink_mom_num][t_eff]  = hsumHq[sink_mom_num][t];
	  HQprop_m[sink_mom_num][t_eff]  = hsumHq_m[sink_mom_num][t];
	}
    
    //XMLWriter xml_bar(xml);
    push(xml, xml_group);
    write(xml, "gamma_value", gind);
    write(xml, "ThreePointColorUnmixed", HQprop[0]);
    //    pop(xml);
    //    push(xml, xml_group);
    write(xml, "ThreePointColorMixed", HQprop_m[0]);
    pop(xml);
    
    END_CODE();
  }// gind loop
}// QlQl routine
  
  
}// Namespace Chroma
