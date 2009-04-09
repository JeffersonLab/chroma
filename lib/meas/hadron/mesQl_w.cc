// $Id: mesQl_w.cc,v 1.4 2009-04-09 22:57:38 caubin Exp $ 
/*! \file
 *  \brief Heavy Meson (Qlbar)  2-pt function : Orginos and Savage
 */

#include "mesQl_w.h"
#include "barQll_w.h"

namespace Chroma {

//! Heavy-light meson 2-pt function
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * Construct propagators for a heavy-light pseudoscalar meson.
 * In the heavy quark limit the D and the D* are degenerate.  
 * The heavy quark is inserted in the infinitely heavy quark limit
 * by a Wilson-Line without spin indices.
 * We are effectively propagating a spin-1/2 light quark
 *
 * \param u                  gauge field (Read) 
 * \param quark_propagator   quark propagator ( Read )
 * \param src_coord          cartesian coordinates of the source ( Read )
 * \param phases             object holds list of momenta and Fourier phases ( Read )
 * \param xml                xml file object ( Read )
 * \param xml_group          group name for xml data ( Read )
 * \param bc                 boundary condition (default = 0 --> Dirichlet)
 *
 */

void Qlbar(const multi1d<LatticeColorMatrix>& u, 
	   const LatticePropagator& quark_propagator,
	   const multi1d<int>& src_coord, 
	   const SftMom& phases,
	   XMLWriter& xml,
	   const string& xml_group,
	   const int bc)
{
  START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  LatticeColorMatrix Qprop;

  HeavyQuarkProp(Qprop,u,src_coord,length, bc);

  multi1d<DComplex> Hq;
  LatticeComplex Hq_prop ;
 
  // S_proj_unpol = (1/2)(1 + gamma_4)
  // This is the heavy quark dirac structure, but I am using the notation
  // that is currently used in Chroma!!!!!!!!

  SpinMatrix g_one = 1.0;
  SpinMatrix S_proj_unpol = 0.5 * (g_one + (g_one * Gamma(8)));

  // 

  Hq_prop = trace(Qprop * S_proj_unpol *adj(quark_propagator) );
 
 

  // Project onto zero momentum
  multi2d<DComplex> hsumHq ;
  hsumHq = phases.sft(Hq_prop);
  multi2d<DComplex> HQprop(num_mom,length) ;

  for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
    for(int t = 0; t < length; ++t)
      {
	int t_eff = (t - src_coord[Nd-1] + length) % length;
	HQprop[sink_mom_num][t_eff]  = hsumHq[sink_mom_num][t];
      }


  //XMLWriter xml_bar(xml);
  push(xml, xml_group);
  write(xml, "HeavyLight", HQprop[0]);
  pop(xml);

  END_CODE();
}

//! Heavy-light meson 2-pt function with backwards moving static quark
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions! 
 *
 * Construct propagators for a heavy-light pseudoscalar meson.
 * In the heavy quark limit the D and the D* are degenerate.  
 * The heavy quark is inserted in the infinitely heavy quark limit
 * by a Wilson-Line without spin indices.
 * We are effectively propagating a spin-1/2 light quark
 *
 * \param u                  gauge field (Read) 
 * \param quark_propagator   quark propagator ( Read )
 * \param src_coord          cartesian coordinates of the source ( Read )
 * \param phases             object holds list of momenta and Fourier phases ( Read )
 * \param xml                xml file object ( Read )
 * \param xml_group          group name for xml data ( Read )
 * \param bc                 boundary condition (default = 0 --> Dirichlet)
 *
 */

void QlbarBACK(const multi1d<LatticeColorMatrix>& u, 
	   const LatticePropagator& quark_propagator,
	   const multi1d<int>& src_coord, 
	   const SftMom& phases,
	   XMLWriter& xml,
	   const string& xml_group,
	   const int bc)
{
  START_CODE();
  
  if ( Ns != 4 )		/* Code is specific to Ns=4 */
    return;

  int length = phases.numSubsets() ;
  int num_mom = phases.numMom();
  
  LatticeColorMatrix Qprop;

  HeavyQuarkPropBack(Qprop,u,src_coord,length,bc);

  multi1d<DComplex> Hq;
  LatticeComplex Hq_prop ;
 
  // S_proj_unpol = (1/2)(1 + gamma_4)
  // This is the heavy quark dirac structure, but I am using the notation
  // that is currently used in Chroma!!!!!!!!

  SpinMatrix g_one = 1.0;
  SpinMatrix S_proj_unpol = 0.5 * (g_one + (g_one * Gamma(8)));

  Hq_prop = trace(Qprop * S_proj_unpol *adj(quark_propagator) );

  // Project onto zero momentum
  multi2d<DComplex> hsumHq ;
  hsumHq = phases.sft(Hq_prop);
  multi2d<DComplex> HQprop(num_mom,length) ;

  for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
    for(int t = 0; t < length; ++t)
      {
	int t_eff = (t - src_coord[Nd-1] + length) % length;
	HQprop[sink_mom_num][t_eff]  = hsumHq[sink_mom_num][t];
      }

  //XMLWriter xml_bar(xml);
  push(xml, xml_group);
  write(xml, "HeavyLight", HQprop[0]);
  pop(xml);

  END_CODE();
}

}
