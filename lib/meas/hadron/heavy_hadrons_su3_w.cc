// $Id: heavy_hadrons_su3_w.cc,v 1.3 2008-12-21 21:22:36 edwards Exp $ 
/*! \file
 *  \brief Heavy hadrons in su3 : Detmold
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "barQll_w.h"
#include "heavy_hadron_potentials_w.h"
#include "heavy_hadrons_su3_w.h"

namespace Chroma
{

  //! Heavy hadron spectrum for SU(3) isospin limit
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * \param u                  gauge field (Read) 
   * \param quark1             light quark propagator ( Read )
   * \param quark2             strange quark propagator ( Read )
   * \param src                cartesian coordinates of one source "0"( Read )
   * \param phases             object holds list of momenta and Fourier phases ( Read )
   * \param xml                xml file object ( Read )
   * \param xml_group          group name for xml data ( Read )
   *
   */

  void static_light_su3(const multi1d<LatticeColorMatrix>& u, 
			const LatticePropagator& quark1,
			const LatticePropagator& quark2,
			const multi1d<int>& src,
			const SftMom& phases,
			XMLWriter& xml,
			const string& xml_group)
  {
    START_CODE();
  
    if ( Ns != 4 )		/* Code is specific to Ns=4 */
      return;

    int length = phases.numSubsets() ;
    int Nt = length;
    int num_mom = phases.numMom();

    //Spin matrices
    SpinMatrix g_one = 1.0;
    SpinMatrix Cg5 = g_one * Gamma(5);
    SpinMatrix Cg3 = - g_one * Gamma(14);
    SpinMatrix Cgminus = g_one * Gamma(11) + timesI(g_one * Gamma(8)); 
    SpinMatrix Cgplus  = g_one * Gamma(11) + timesMinusI(g_one * Gamma(8)); 
    SpinMatrix G5 = g_one * Gamma(15);
    SpinMatrix OnePlusSigma3 = g_one + timesMinusI(g_one*Gamma(3));
    SpinMatrix OneMinusSigma3 = g_one + timesI(g_one*Gamma(3));
    SpinMatrix Gup = g_one * Gamma(15) + timesMinusI((g_one*Gamma(3))*Gamma(15));
    SpinMatrix Gdown = g_one * Gamma(15) + timesI((g_one*Gamma(3))*Gamma(15));
    SpinMatrix PosEnergyProj =  0.5*((g_one + Gamma(8) * g_one));
    SpinMatrix NegEnergyProj =  0.5*((g_one - Gamma(8) * g_one));
    SpinMatrix OnePlusSigma3PEP = OnePlusSigma3 * PosEnergyProj;
    SpinMatrix OneMinusSigma3PEP = OneMinusSigma3 * PosEnergyProj;
    SpinMatrix G5PEP = G5 * PosEnergyProj;


    // Make heavy quark propagators
    LatticeColorMatrix Qprop1;
    HeavyQuarkProp(Qprop1,u,src,length);

    //peek propagators: reduce lattice problem to two site problem
    multi1d<DPropagator> U1,D1,S1;
    multi1d<ColorMatrix> Q1, antiQ1;
    Q1.resize(length);
    antiQ1.resize(length);
    U1.resize(length);
    S1.resize(length);
    U1=0; S1=0;
    multi1d<int> currsrc=src;
    for (int t=0; t<length; t++)
    {
      currsrc[Nd-1]=t;
      Q1[t] = peekSite(Qprop1,currsrc);  // HQ prop from src1 == "0"
      antiQ1[t] = adj(Q1[t]);  // antiHQ prop from src1 == "0"

      U1[t] = peekSite(quark1,currsrc);  // light propagator 
      S1[t] = peekSite(quark2,currsrc);  // strange propagator

    }

    D1=U1;


    // Pseudoscalar/vector meson blocks (these are degenerate)

    QDPIO::cout<<"Making B^+ meson blocks\n";
    HeavyMesonBlock H_U1_G5_0(Nt,U1,G5,antiQ1,NegEnergyProj);

    // B meson
    QDPIO::cout<<"  Contracting B meson \n";
    multi1d<DComplex> bmes;
    bmes.resize(length);  
    bmes = bcontract(H_U1_G5_0,G5);

    QDPIO::cout<<"Making B_s meson blocks\n";
    HeavyMesonBlock H_S1_G5_0(Nt,S1,G5,antiQ1,NegEnergyProj);

    // B_s meson
    QDPIO::cout<<"  Contracting B_s meson \n";
    multi1d<DComplex> bsmes;
    bsmes.resize(length);  
    bsmes = bcontract(H_S1_G5_0,G5);


    // BARYON BLOCKS 

    // Add epsilon tensors onto ends of HQ propagators
    multiNd<DComplex> HQB0;
    multiNd<DComplex> antiHQB0;
    multi1d<int> HQBarray;
    HQBarray.resize(4);
    HQBarray= Nc; HQBarray[0]=length;
    HQB0.resize(HQBarray);
    HQB0 = HBQfunc(Q1);
    antiHQB0.resize(HQBarray);
    antiHQB0 = HBQfunc(antiQ1);

    QDPIO::cout<<"Making lambdaB blocks\n";
    // Lambda_b
    // Note we will use spin to differentiate the Lambda_b from the \Sigma_B^0
    // and simply take Bud to be the flavour structure
    QllBlock B_U1_D1_0_Cg5(Nt,U1,D1,Cg5,HQB0);
  
    // Lambda_b 2 pt
    QDPIO::cout<<"  Contracting Lambda_b \n";
    multi1d<DComplex> lambdab;
    lambdab.resize(length); 
    lambdab = lambdabcontract(B_U1_D1_0_Cg5, Cg5);

    QDPIO::cout<<"Making xiB0 blocks\n";
    // Lambda_b
    // Note we will use spin to differentiate the Lambda_b from the \Sigma_B^0
    // and simply take Bud to be the flavour structure
    QllBlock B_U1_S1_0_Cg5(Nt,U1,S1,Cg5,HQB0);
  
    // Xi_b^0 2 pt
    QDPIO::cout<<"  Contracting xi_b^0 \n";
    multi1d<DComplex> xibzero;
    xibzero.resize(length); 
    xibzero = lambdabcontract(B_U1_S1_0_Cg5, Cg5);

    QDPIO::cout<<"Making sigmaB^+ blocks\n";
    QllBlock B_U1_U1_0_Cg3(Nt,U1,U1,Cg3,HQB0);
    QllBlock B_U1_U1_0_Cgplus(Nt,U1,U1,Cgplus,HQB0);
    QllBlock B_U1_U1_0_Cgminus(Nt,U1,U1,Cgminus,HQB0);

    // Sigma_b^+ and (=) \Sigma_b^- baryon 2pt
    QDPIO::cout<<"  Contracting Sigma_b^+ \n";
    multi1d<DComplex> sigmabplusJ1m1;
    multi1d<DComplex> sigmabplusJ1m0;
    multi1d<DComplex> sigmabplusJ1mneg1;
    sigmabplusJ1m1.resize(length);
    sigmabplusJ1m0.resize(length);
    sigmabplusJ1mneg1.resize(length);

    sigmabplusJ1m1 = sigmabpluscontract(B_U1_U1_0_Cgplus, Cgplus);
    sigmabplusJ1m0 = sigmabpluscontract(B_U1_U1_0_Cg3, Cg3);
    sigmabplusJ1mneg1 = sigmabpluscontract(B_U1_U1_0_Cgminus, Cgminus);

    QDPIO::cout<<"Making xiB^prime0 blocks\n";
    QllBlock B_U1_S1_0_Cg3(Nt,U1,S1,Cg3,HQB0);
    QllBlock B_U1_S1_0_Cgplus(Nt,U1,S1,Cgplus,HQB0);
    QllBlock B_U1_S1_0_Cgminus(Nt,U1,S1,Cgminus,HQB0);

    // xiB^prime0 2pt
    QDPIO::cout<<"  Contracting xiB^prime0 \n";
    multi1d<DComplex> xibprime0J1m1;
    multi1d<DComplex> xibprime0J1m0;
    multi1d<DComplex> xibprime0J1mneg1;
    xibprime0J1m1.resize(length);
    xibprime0J1m0.resize(length);
    xibprime0J1mneg1.resize(length);

    xibprime0J1m1 = sigmabpluscontract(B_U1_S1_0_Cgplus, Cgplus);
    xibprime0J1m0 = sigmabpluscontract(B_U1_S1_0_Cg3, Cg3);
    xibprime0J1mneg1 = sigmabpluscontract(B_U1_S1_0_Cgminus, Cgminus);

    QDPIO::cout<<"Making omegaB^- blocks\n";
    QllBlock B_S1_S1_0_Cg3(Nt,S1,S1,Cg3,HQB0);
    QllBlock B_S1_S1_0_Cgplus(Nt,S1,S1,Cgplus,HQB0);
    QllBlock B_S1_S1_0_Cgminus(Nt,S1,S1,Cgminus,HQB0);

    // omegaB^- 2pt
    QDPIO::cout<<"  Contracting omegaB^- \n";
    multi1d<DComplex> omegabminusJ1m1;
    multi1d<DComplex> omegabminusJ1m0;
    multi1d<DComplex> omegabminusJ1mneg1;
    omegabminusJ1m1.resize(length);
    omegabminusJ1m0.resize(length);
    omegabminusJ1mneg1.resize(length);

    omegabminusJ1m1 = sigmabpluscontract(B_S1_S1_0_Cgplus, Cgplus);
    omegabminusJ1m0 = sigmabpluscontract(B_S1_S1_0_Cg3, Cg3);
    omegabminusJ1mneg1 = sigmabpluscontract(B_S1_S1_0_Cgminus, Cgminus);



    push(xml, xml_group);
    // Hadrons
    write(xml, "Bu",bmes);
    write(xml, "Bs",bsmes);
    write(xml,"lambdab", lambdab);
    write(xml,"xibzero", xibzero);
    write(xml,"sigmabplusJ1m1", sigmabplusJ1m1);
    write(xml,"sigmabplusJ1m0", sigmabplusJ1m0);
    write(xml,"sigmabplusJ1mneg1", sigmabplusJ1mneg1);
    write(xml,"xibprime0J1m1", xibprime0J1m1);
    write(xml,"xibprime0J1m0", xibprime0J1m0);
    write(xml,"xibprime0J1mneg1", xibprime0J1mneg1);
    write(xml,"omegabminusJ1m1", omegabminusJ1m1);
    write(xml,"omegabminusJ1m0", omegabminusJ1m0);
    write(xml,"omegabminusJ1mneg1", omegabminusJ1mneg1);
    pop(xml);

    END_CODE();
  }



}
