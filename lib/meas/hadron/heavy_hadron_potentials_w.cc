// $Id: heavy_hadron_potentials_w.cc,v 1.3 2008-12-21 21:45:47 edwards Exp $ 
/*! \file
 *  \brief Potential between 2 heavy hadrons : Detmold
 *  Correlators checked independentely by Savage
 *

 * Includes Lambda_b's etc

 * Soon also will have the BBB "potential"

 */

#include "util/ferm/antisymtensor.h"
#include "barQll_w.h"
#include "mesQl_w.h"
#include "heavy_hadron_potentials_w.h"

namespace Chroma
{
  //! Heavy hadron potentials for SU(2) isospin limit
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *

   * \param u                  gauge field (Read) 
   * \param quark1             quark propagator, src at 0 ( Read )
   * \param quark2             quark propagator, src at R ( Read )
   * \param src1               cartesian coordinates of one source "0"( Read )
   * \param src2               cartesian coordinates of the other source "R" ( Read )
   * \param phases             object holds list of momenta and Fourier phases ( Read )
   * \param xml                xml file object ( Read )
   * \param xml_group          group name for xml data ( Read )
   *
   */


  void QllQllPOT(const multi1d<LatticeColorMatrix>& u, 
		 const LatticePropagator& quark1,
		 const LatticePropagator& quark2,
		 const multi1d<int>& src1, 
		 const multi1d<int>& src2, 
		 const SftMom& phases,
		 XMLWriter& xml,
		 const string& xml_group)
  {
    START_CODE();
  
    if ( Ns != 4 )		/* Code is specific to Ns=4 */
      return;

    if (src1[Nd-1] != src2[Nd-1] ){
      QDPIO::cerr<<"Sources must be at the same time slice\n";
      return ;
    }

    QDPIO::cout << "Sources at : " << src1[0] << " " << src1[1] 
		<< " " << src1[2] << " " << src1[3] << endl
		<< "       and : " << src2[0] << " " << src2[1] 
		<< " " << src2[2] << " " << src2[3] << endl; 

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
    LatticeColorMatrix Qprop2;
    HeavyQuarkProp(Qprop1,u,src1,length);
    HeavyQuarkProp(Qprop2,u,src2,length);

    //peek propagators: reduce lattice problem to two site problem
    multi1d<DPropagator> U1,U2,U3,U4;
    multi1d<DPropagator> D1,D2,D3,D4;
    multi1d<DPropagator> antiU1,antiU2,antiU3,antiU4;
    multi1d<DPropagator> antiD1,antiD2,antiD3,antiD4;
    multi1d<ColorMatrix> Q1, Q2, antiQ1, antiQ2;
    Q1.resize(length);
    Q2.resize(length);
    antiQ1.resize(length);
    antiQ2.resize(length);
    U1.resize(length); U2.resize(length); U3.resize(length);U4.resize(length);
    D1.resize(length); D2.resize(length); D3.resize(length);D4.resize(length);
    antiU1.resize(length); antiU2.resize(length); antiU3.resize(length);antiU4.resize(length);
    antiD1.resize(length); antiD2.resize(length); antiD3.resize(length);antiD4.resize(length);
    U1=0;U2=0;U3=0;U4=0;
    multi1d<int> currsrc1=src1;
    multi1d<int> currsrc2=src2;
    for (int t=0; t<length; t++)
    {
      currsrc1[Nd-1]=t;
      currsrc2[Nd-1]=t;
      Q1[t] = peekSite(Qprop1,currsrc1);  // HQ prop from src1 == "0"
      Q2[t] = peekSite(Qprop2,currsrc2);  // HQ prop from src2 == "R"
      antiQ1[t] = adj(Q1[t]);  // antiHQ prop from src1 == "0"
      antiQ2[t] = adj(Q2[t]);  // antiHQ prop from src2 == "R"

      U1[t] = peekSite(quark1,currsrc1);  // propagator from 0 to 0
      U2[t] = peekSite(quark2,currsrc2);  // propagator from R to R
      U3[t] = peekSite(quark2,currsrc1);  // propagator from R to 0
      U4[t] = peekSite(quark1,currsrc2);  // propagator from 0 to R

      antiU1[t]= Gamma(15) * adj(U1[t]) * Gamma(15);
      antiU2[t]= Gamma(15) * adj(U2[t]) * Gamma(15);
      antiU3[t]= Gamma(15) * adj(U3[t]) * Gamma(15);
      antiU4[t]= Gamma(15) * adj(U4[t]) * Gamma(15);

    }

    D1 = U1; D2 = U2; D3 = U3; D4 = U4; // D quarks are the same as U quarks
    antiD1 = antiU1; antiD2 = antiU2; antiD3 = antiU3; antiD4 = antiU4;

    // MESON BLOCKS 

    QDPIO::cout<<"Making B^+ (bbar u) meson blocks\n";

    HeavyMesonBlock H_U4_G5_R(Nt,U4,G5,antiQ2,NegEnergyProj);
    HeavyMesonBlock H_U2_G5_R(Nt,U2,G5,antiQ2,NegEnergyProj);
    HeavyMesonBlock H_U3_G5_0(Nt,U3,G5,antiQ1,NegEnergyProj);
    HeavyMesonBlock H_U1_G5_0(Nt,U1,G5,antiQ1,NegEnergyProj);

    HeavyMesonBlock H_U4_Gup_R(Nt,U4,Gup,antiQ2,NegEnergyProj);
    HeavyMesonBlock H_U2_Gup_R(Nt,U2,Gup,antiQ2,NegEnergyProj);
    HeavyMesonBlock H_U3_Gup_0(Nt,U3,Gup,antiQ1,NegEnergyProj);
    HeavyMesonBlock H_U1_Gup_0(Nt,U1,Gup,antiQ1,NegEnergyProj);

    HeavyMesonBlock H_U4_Gdown_R(Nt,U4,Gdown,antiQ2,NegEnergyProj);
    HeavyMesonBlock H_U2_Gdown_R(Nt,U2,Gdown,antiQ2,NegEnergyProj);
    HeavyMesonBlock H_U3_Gdown_0(Nt,U3,Gdown,antiQ1,NegEnergyProj);
    HeavyMesonBlock H_U1_Gdown_0(Nt,U1,Gdown,antiQ1,NegEnergyProj);

    QDPIO::cout<<"Making B^0 (bbar d) meson blocks\n";
    HeavyMesonBlock H_D4_G5_R = H_U4_G5_R;
    HeavyMesonBlock H_D2_G5_R = H_U2_G5_R;
    HeavyMesonBlock H_D3_G5_0 = H_U3_G5_0;
    HeavyMesonBlock H_D1_G5_0 = H_U1_G5_0;

    HeavyMesonBlock H_D4_Gup_R = H_U4_Gup_R;
    HeavyMesonBlock H_D2_Gup_R = H_U2_Gup_R;
    HeavyMesonBlock H_D3_Gup_0 = H_U3_Gup_0;
    HeavyMesonBlock H_D1_Gup_0 = H_U1_Gup_0;

    HeavyMesonBlock H_D4_Gdown_R = H_U4_Gdown_R;
    HeavyMesonBlock H_D2_Gdown_R = H_U2_Gdown_R;
    HeavyMesonBlock H_D3_Gdown_0 = H_U3_Gdown_0;
    HeavyMesonBlock H_D1_Gdown_0 = H_U1_Gdown_0;

    QDPIO::cout<<"Making Bbar^0 (b dbar) meson blocks\n";
    HeavyMesonBlock Hbar_D4_G5_R(Nt,antiD4,G5,Q2,PosEnergyProj);
    HeavyMesonBlock Hbar_D2_G5_R(Nt,antiD2,G5,Q2,PosEnergyProj);
    HeavyMesonBlock Hbar_D3_G5_0(Nt,antiD3,G5,Q1,PosEnergyProj);
    HeavyMesonBlock Hbar_D1_G5_0(Nt,antiD1,G5,Q1,PosEnergyProj);

    HeavyMesonBlock Hbar_D4_Gup_R(Nt,antiD4,Gup,Q2,PosEnergyProj);
    HeavyMesonBlock Hbar_D2_Gup_R(Nt,antiD2,Gup,Q2,PosEnergyProj);
    HeavyMesonBlock Hbar_D3_Gup_0(Nt,antiD3,Gup,Q1,PosEnergyProj);
    HeavyMesonBlock Hbar_D1_Gup_0(Nt,antiD1,Gup,Q1,PosEnergyProj);

    HeavyMesonBlock Hbar_D4_Gdown_R(Nt,antiD4,Gdown,Q2,PosEnergyProj);
    HeavyMesonBlock Hbar_D2_Gdown_R(Nt,antiD2,Gdown,Q2,PosEnergyProj);
    HeavyMesonBlock Hbar_D3_Gdown_0(Nt,antiD3,Gdown,Q1,PosEnergyProj);
    HeavyMesonBlock Hbar_D1_Gdown_0(Nt,antiD1,Gdown,Q1,PosEnergyProj);


    // B meson
    QDPIO::cout<<"  Contracting B meson \n";
    multi1d<DComplex> bmes, bmesup, bmesdown, bmes0, bmes0up, bmes0down, bmesR, bmesRup, bmesRdown;
    bmes.resize(length);  bmesup.resize(length);  bmesdown.resize(length);  
    bmes0.resize(length);   bmes0up.resize(length);  bmes0down.resize(length);  
    bmesR.resize(length); bmesRup.resize(length);  bmesRdown.resize(length);  
    bmes0 = bcontract(H_U1_G5_0,G5);
    bmesR = bcontract(H_U2_G5_R,G5);
    bmes0up = bcontract(H_U1_Gup_0,Gup);
    bmes0down = bcontract(H_U1_Gdown_0,Gdown);
    bmesRup = bcontract(H_U2_Gup_R,Gup);
    bmesRdown = bcontract(H_U2_Gdown_R,Gdown);
    for (int t=0; t<length;t++){
      bmes[t] = 0.5 * (bmes0[t] + bmesR[t]);
      bmesup[t] = 0.5 * (bmes0up[t] + bmesRup[t]);
      bmesdown[t] = 0.5 * (bmes0down[t] + bmesRdown[t]);
    };

    // B B I=I_3=1: Bu Bu
    QDPIO::cout<<"  Contracting Bu--Bu \n";
    multi1d<DComplex> BuBuJ0m0,BuBuJ1m1,BuBuJ1m0,BuBuJ1mneg1;
    // J=0
    BuBuJ0m0.resize(length);
    BuBuJ0m0 = m1contract(H_U1_Gup_0, H_U2_Gdown_R,
			  H_U3_Gup_0, H_U4_Gdown_R,
			  Gup, Gdown);
    BuBuJ0m0 += m1contract(H_U1_Gdown_0, H_U2_Gup_R,
			   H_U3_Gdown_0, H_U4_Gup_R,
			   Gdown, Gup);
    BuBuJ0m0 -= m1contract(H_U1_Gup_0, H_U2_Gdown_R,
			   H_U3_Gup_0, H_U4_Gdown_R,
			   Gdown, Gup);
    BuBuJ0m0 -= m1contract(H_U1_Gdown_0, H_U2_Gup_R,
			   H_U3_Gdown_0, H_U4_Gup_R,
			   Gup, Gdown);
    BuBuJ0m0 *= 0.5;

    //J=1
    BuBuJ1m1.resize(length);
    BuBuJ1m0.resize(length);
    BuBuJ1mneg1.resize(length);

    //J=1, m=1
    BuBuJ1m1 = m1contract(H_U1_Gup_0, H_U2_Gup_R,
			  H_U3_Gup_0, H_U4_Gup_R,
			  Gup, Gup);

    //J=1, m=0
    BuBuJ1m0 = m1contract(H_U1_Gup_0, H_U2_Gdown_R,
			  H_U3_Gup_0, H_U4_Gdown_R,
			  Gup, Gdown);
    BuBuJ1m0 += m1contract(H_U1_Gdown_0, H_U2_Gup_R,
			   H_U3_Gdown_0, H_U4_Gup_R,
			   Gdown, Gup);
    BuBuJ1m0 += m1contract(H_U1_Gup_0, H_U2_Gdown_R,
			   H_U3_Gup_0, H_U4_Gdown_R,
			   Gdown, Gup);
    BuBuJ1m0 += m1contract(H_U1_Gdown_0, H_U2_Gup_R,
			   H_U3_Gdown_0, H_U4_Gup_R,
			   Gup, Gdown);
    BuBuJ1m0 *= 0.5;

    //J=1, m=-1
    BuBuJ1mneg1 = m1contract(H_U1_Gdown_0, H_U2_Gdown_R,
			     H_U3_Gdown_0, H_U4_Gdown_R,
			     Gdown, Gdown);

    // B B I=0,1 with I_3=0: Bu Bd
    // since there are two different particles, we might symmetrise over 
    // Bu(r) Bd(0) -->  Bu(r) Bd(0)  and Bd(r) Bu(0) -->  Bd(r) Bu(0)  
    // but since the u and d quark propagators are the same this is
    // pointless in this case

    QDPIO::cout<<"  Contracting Bu--Bd \n";
    multi1d<DComplex> BuBdJ0m0,BuBdJ1m1,BuBdJ1m0,BuBdJ1mneg1;
    BuBdJ0m0.resize(length);
    BuBdJ0m0 = m2contract(H_U1_Gup_0, H_D2_Gdown_R,
			  Gup, Gdown);
    BuBdJ0m0 -= m2contract(H_U1_Gdown_0, H_D2_Gup_R,
			   Gup, Gdown);
    BuBdJ0m0 += m2contract(H_U1_Gdown_0, H_D2_Gup_R,
			   Gdown, Gup);
    BuBdJ0m0 -= m2contract(H_U1_Gup_0, H_D2_Gdown_R,
			   Gdown, Gup);
    BuBdJ0m0 *= 0.5;


    BuBdJ1m1.resize(length);
    BuBdJ1m0.resize(length);
    BuBdJ1mneg1.resize(length);

    SpinMatrix G5Gup = Gamma(15)* Gup;

    BuBdJ1m1 = m2contract(H_U1_Gup_0, H_D2_Gup_R,
			  Gup, Gup);
  
    BuBdJ1m0 = m2contract(H_U1_Gup_0, H_D2_Gdown_R,
			  Gup, Gdown);
    BuBdJ1m0 += m2contract(H_U1_Gdown_0, H_D2_Gup_R,
			   Gup, Gdown);
    BuBdJ1m0 += m2contract(H_U1_Gdown_0, H_D2_Gup_R,
			   Gdown, Gup);
    BuBdJ1m0 += m2contract(H_U1_Gup_0, H_D2_Gdown_R,
			   Gdown, Gup);
    BuBdJ1m0 *= 0.5;


    BuBdJ1mneg1 = m2contract(H_U1_Gdown_0, H_D2_Gdown_R,
			     Gdown, Gdown);

    // B Bbar I=0,1 with I_3=0: Bu Bdbar
    // since there are two different particles, we might symmetrise over 
    // Bu(r) Bd(0) -->  Bu(r) Bd(0)  and Bd(r) Bu(0) -->  Bd(r) Bu(0)  
    // but since the u and d quark propagators are the same this is
    // pointless in this case

    QDPIO::cout<<"  Contracting Bu--Bdbar \n";
    multi1d<DComplex> BuBdbarJ0m0,BuBdbarJ1m1,BuBdbarJ1m0,BuBdbarJ1mneg1;
    BuBdbarJ0m0.resize(length);
    BuBdbarJ0m0 = m2contract(H_U1_Gup_0, Hbar_D2_Gdown_R,
			     Gup, Gdown);
    BuBdbarJ0m0 -= m2contract(H_U1_Gdown_0, Hbar_D2_Gup_R,
			      Gup, Gdown);
    BuBdbarJ0m0 += m2contract(H_U1_Gdown_0, Hbar_D2_Gup_R,
			      Gdown, Gup);
    BuBdbarJ0m0 -= m2contract(H_U1_Gup_0, Hbar_D2_Gdown_R,
			      Gdown, Gup);
    BuBdbarJ0m0 *= 0.5;


    BuBdbarJ1m1.resize(length);
    BuBdbarJ1m0.resize(length);
    BuBdbarJ1mneg1.resize(length);

    BuBdbarJ1m1 = m2contract(H_U1_Gup_0, Hbar_D2_Gup_R,
			     Gup, Gup);


    BuBdbarJ1m0 = m2contract(H_U1_Gup_0, Hbar_D2_Gdown_R,
			     Gup, Gdown);
    BuBdbarJ1m0 += m2contract(H_U1_Gdown_0, Hbar_D2_Gup_R,
			      Gup, Gdown);
    BuBdbarJ1m0 += m2contract(H_U1_Gdown_0, Hbar_D2_Gup_R,
			      Gdown, Gup);
    BuBdbarJ1m0 += m2contract(H_U1_Gup_0, Hbar_D2_Gdown_R,
			      Gdown, Gup);
    BuBdbarJ1m0 *= 0.5;


    BuBdbarJ1mneg1 = m2contract(H_U1_Gdown_0, Hbar_D2_Gdown_R,
				Gdown, Gdown);


    // BARYON BLOCKS 

    // Add epsilon tensors onto ends of HQ propagators
    multiNd<DComplex> HQBR, HQB0;
    multiNd<DComplex> antiHQBR, antiHQB0;
    multi1d<int> HQBarray;
    HQBarray.resize(4);
    HQBarray= Nc; HQBarray[0]=length;
    HQBR.resize(HQBarray);
    HQB0.resize(HQBarray);
    HQB0 = HBQfunc(Q1);
    HQBR = HBQfunc(Q2);
    antiHQBR.resize(HQBarray);
    antiHQB0.resize(HQBarray);
    antiHQB0 = HBQfunc(antiQ1);
    antiHQBR = HBQfunc(antiQ2);

    QDPIO::cout<<"Making lambdaB blocks\n";
    // Lambda_b
    // Note we will use spin to differentiate the Lambda_b from the \Sigma_B^0
    // and simply take Bud to be the flavour structure
    QllBlock B_U4_D2_R_Cg5(Nt,U4,D2,Cg5,HQBR);
    QllBlock B_U3_D1_0_Cg5(Nt,U3,D1,Cg5,HQB0);
    QllBlock B_U2_D4_R_Cg5(Nt,U2,D4,Cg5,HQBR);
    QllBlock B_U1_D3_0_Cg5(Nt,U1,D3,Cg5,HQB0);
    QllBlock B_U4_D4_R_Cg5(Nt,U4,D4,Cg5,HQBR);
    QllBlock B_U3_D3_0_Cg5(Nt,U3,D3,Cg5,HQB0);
    QllBlock B_U2_D2_R_Cg5(Nt,U2,D2,Cg5,HQBR);
    QllBlock B_U1_D1_0_Cg5(Nt,U1,D1,Cg5,HQB0);
  
    // Lambda_b 2 pt
    QDPIO::cout<<"  Contracting Lambda_b \n";
    multi1d<DComplex> lambdab,lambdaba,lambdabb;
    lambdab.resize(length); lambdaba.resize(length); lambdabb.resize(length);
    lambdaba = lambdabcontract(B_U1_D1_0_Cg5, Cg5);
    lambdabb = lambdabcontract(B_U2_D2_R_Cg5, Cg5);
    for (int t=0; t<length;t++){
      lambdab[t]= 0.5 * ( lambdaba[t]+lambdabb[t]);
    };
  
    // Lambda_b Lambda_b contractions
    QDPIO::cout<<"  Contracting Lambda_b--Lambda_b \n";
    multi1d<DComplex> lambdablambdab;
    lambdablambdab.resize(length);
    lambdablambdab = c1contract(B_U1_D1_0_Cg5,  B_U1_D3_0_Cg5, 
				B_U2_D2_R_Cg5,  B_U2_D4_R_Cg5, 
				B_U3_D1_0_Cg5,  B_U3_D3_0_Cg5, 
				B_U4_D2_R_Cg5,  B_U4_D4_R_Cg5,
				Cg5, Cg5);

    QDPIO::cout<<"Making sigmaB^+ blocks\n";
    QllBlock B_U4_U2_R_Cg3(Nt,U4,U2,Cg3,HQBR);
    QllBlock B_U3_U1_0_Cg3(Nt,U3,U1,Cg3,HQB0);
    QllBlock B_U2_U4_R_Cg3(Nt,U2,U4,Cg3,HQBR);
    QllBlock B_U1_U3_0_Cg3(Nt,U1,U3,Cg3,HQB0);
    QllBlock B_U4_U4_R_Cg3(Nt,U4,U4,Cg3,HQBR);
    QllBlock B_U3_U3_0_Cg3(Nt,U3,U3,Cg3,HQB0);
    QllBlock B_U2_U2_R_Cg3(Nt,U2,U2,Cg3,HQBR);
    QllBlock B_U1_U1_0_Cg3(Nt,U1,U1,Cg3,HQB0);

    QllBlock B_U4_U2_R_Cgplus(Nt,U4,U2,Cgplus,HQBR);
    QllBlock B_U3_U1_0_Cgplus(Nt,U3,U1,Cgplus,HQB0);
    QllBlock B_U2_U4_R_Cgplus(Nt,U2,U4,Cgplus,HQBR);
    QllBlock B_U1_U3_0_Cgplus(Nt,U1,U3,Cgplus,HQB0);
    QllBlock B_U4_U4_R_Cgplus(Nt,U4,U4,Cgplus,HQBR);
    QllBlock B_U3_U3_0_Cgplus(Nt,U3,U3,Cgplus,HQB0);
    QllBlock B_U2_U2_R_Cgplus(Nt,U2,U2,Cgplus,HQBR);
    QllBlock B_U1_U1_0_Cgplus(Nt,U1,U1,Cgplus,HQB0);

    QllBlock B_U4_U2_R_Cgminus(Nt,U4,U2,Cgminus,HQBR);
    QllBlock B_U3_U1_0_Cgminus(Nt,U3,U1,Cgminus,HQB0);
    QllBlock B_U2_U4_R_Cgminus(Nt,U2,U4,Cgminus,HQBR);
    QllBlock B_U1_U3_0_Cgminus(Nt,U1,U3,Cgminus,HQB0);
    QllBlock B_U4_U4_R_Cgminus(Nt,U4,U4,Cgminus,HQBR);
    QllBlock B_U3_U3_0_Cgminus(Nt,U3,U3,Cgminus,HQB0);
    QllBlock B_U2_U2_R_Cgminus(Nt,U2,U2,Cgminus,HQBR);
    QllBlock B_U1_U1_0_Cgminus(Nt,U1,U1,Cgminus,HQB0);

    QDPIO::cout<<"Making sigmaB^0 blocks\n";
    QllBlock B_U4_D2_R_Cg3(Nt,U4,D2,Cg3,HQBR);
    QllBlock B_U3_D1_0_Cg3(Nt,U3,D1,Cg3,HQB0);
    QllBlock B_U2_D4_R_Cg3(Nt,U2,D4,Cg3,HQBR);
    QllBlock B_U1_D3_0_Cg3(Nt,U1,D3,Cg3,HQB0);
    QllBlock B_U4_D4_R_Cg3(Nt,U4,D4,Cg3,HQBR);
    QllBlock B_U3_D3_0_Cg3(Nt,U3,D3,Cg3,HQB0);
    QllBlock B_U2_D2_R_Cg3(Nt,U2,D2,Cg3,HQBR);
    QllBlock B_U1_D1_0_Cg3(Nt,U1,D1,Cg3,HQB0);

    QllBlock B_U4_D2_R_Cgplus(Nt,U4,D2,Cgplus,HQBR);
    QllBlock B_U3_D1_0_Cgplus(Nt,U3,D1,Cgplus,HQB0);
    QllBlock B_U2_D4_R_Cgplus(Nt,U2,D4,Cgplus,HQBR);
    QllBlock B_U1_D3_0_Cgplus(Nt,U1,D3,Cgplus,HQB0);
    QllBlock B_U4_D4_R_Cgplus(Nt,U4,D4,Cgplus,HQBR);
    QllBlock B_U3_D3_0_Cgplus(Nt,U3,D3,Cgplus,HQB0);
    QllBlock B_U2_D2_R_Cgplus(Nt,U2,D2,Cgplus,HQBR);
    QllBlock B_U1_D1_0_Cgplus(Nt,U1,D1,Cgplus,HQB0);

    QllBlock B_U4_D2_R_Cgminus(Nt,U4,D2,Cgminus,HQBR);
    QllBlock B_U3_D1_0_Cgminus(Nt,U3,D1,Cgminus,HQB0);
    QllBlock B_U2_D4_R_Cgminus(Nt,U2,D4,Cgminus,HQBR);
    QllBlock B_U1_D3_0_Cgminus(Nt,U1,D3,Cgminus,HQB0);
    QllBlock B_U4_D4_R_Cgminus(Nt,U4,D4,Cgminus,HQBR);
    QllBlock B_U3_D3_0_Cgminus(Nt,U3,D3,Cgminus,HQB0);
    QllBlock B_U2_D2_R_Cgminus(Nt,U2,D2,Cgminus,HQBR);
    QllBlock B_U1_D1_0_Cgminus(Nt,U1,D1,Cgminus,HQB0);


    QDPIO::cout<<"Making sigmaB^- blocks\n";
    QllBlock B_D4_D2_R_Cg3(Nt,D4,D2,Cg3,HQBR);
    QllBlock B_D3_D1_0_Cg3(Nt,D3,D1,Cg3,HQB0);
    QllBlock B_D2_D4_R_Cg3(Nt,D2,D4,Cg3,HQBR);
    QllBlock B_D1_D3_0_Cg3(Nt,D1,D3,Cg3,HQB0);
    QllBlock B_D4_D4_R_Cg3(Nt,D4,D4,Cg3,HQBR);
    QllBlock B_D3_D3_0_Cg3(Nt,D3,D3,Cg3,HQB0);
    QllBlock B_D2_D2_R_Cg3(Nt,D2,D2,Cg3,HQBR);
    QllBlock B_D1_D1_0_Cg3(Nt,D1,D1,Cg3,HQB0);

    QllBlock B_D4_D2_R_Cgplus(Nt,D4,D2,Cgplus,HQBR);
    QllBlock B_D3_D1_0_Cgplus(Nt,D3,D1,Cgplus,HQB0);
    QllBlock B_D2_D4_R_Cgplus(Nt,D2,D4,Cgplus,HQBR);
    QllBlock B_D1_D3_0_Cgplus(Nt,D1,D3,Cgplus,HQB0);
    QllBlock B_D4_D4_R_Cgplus(Nt,D4,D4,Cgplus,HQBR);
    QllBlock B_D3_D3_0_Cgplus(Nt,D3,D3,Cgplus,HQB0);
    QllBlock B_D2_D2_R_Cgplus(Nt,D2,D2,Cgplus,HQBR);
    QllBlock B_D1_D1_0_Cgplus(Nt,D1,D1,Cgplus,HQB0);

    QllBlock B_D4_D2_R_Cgminus(Nt,D4,D2,Cgminus,HQBR);
    QllBlock B_D3_D1_0_Cgminus(Nt,D3,D1,Cgminus,HQB0);
    QllBlock B_D2_D4_R_Cgminus(Nt,D2,D4,Cgminus,HQBR);
    QllBlock B_D1_D3_0_Cgminus(Nt,D1,D3,Cgminus,HQB0);
    QllBlock B_D4_D4_R_Cgminus(Nt,D4,D4,Cgminus,HQBR);
    QllBlock B_D3_D3_0_Cgminus(Nt,D3,D3,Cgminus,HQB0);
    QllBlock B_D2_D2_R_Cgminus(Nt,D2,D2,Cgminus,HQBR);
    QllBlock B_D1_D1_0_Cgminus(Nt,D1,D1,Cgminus,HQB0);


    QDPIO::cout<<"Making anti-sigmaB^- blocks\n";
    QllBlock Bbar_D4_D2_R_Cg3(Nt,antiD4,antiD2,Cg3,antiHQBR);
    QllBlock Bbar_D3_D1_0_Cg3(Nt,antiD3,antiD1,Cg3,antiHQB0);
    QllBlock Bbar_D2_D4_R_Cg3(Nt,antiD2,antiD4,Cg3,antiHQBR);
    QllBlock Bbar_D1_D3_0_Cg3(Nt,antiD1,antiD3,Cg3,antiHQB0);
    QllBlock Bbar_D4_D4_R_Cg3(Nt,antiD4,antiD4,Cg3,antiHQBR);
    QllBlock Bbar_D3_D3_0_Cg3(Nt,antiD3,antiD3,Cg3,antiHQB0);
    QllBlock Bbar_D2_D2_R_Cg3(Nt,antiD2,antiD2,Cg3,antiHQBR);
    QllBlock Bbar_D1_D1_0_Cg3(Nt,antiD1,antiD1,Cg3,antiHQB0);

    QllBlock Bbar_D4_D2_R_Cgplus(Nt,antiD4,antiD2,Cgplus,antiHQBR);
    QllBlock Bbar_D3_D1_0_Cgplus(Nt,antiD3,antiD1,Cgplus,antiHQB0);
    QllBlock Bbar_D2_D4_R_Cgplus(Nt,antiD2,antiD4,Cgplus,antiHQBR);
    QllBlock Bbar_D1_D3_0_Cgplus(Nt,antiD1,antiD3,Cgplus,antiHQB0);
    QllBlock Bbar_D4_D4_R_Cgplus(Nt,antiD4,antiD4,Cgplus,antiHQBR);
    QllBlock Bbar_D3_D3_0_Cgplus(Nt,antiD3,antiD3,Cgplus,antiHQB0);
    QllBlock Bbar_D2_D2_R_Cgplus(Nt,antiD2,antiD2,Cgplus,antiHQBR);
    QllBlock Bbar_D1_D1_0_Cgplus(Nt,antiD1,antiD1,Cgplus,antiHQB0);

    QllBlock Bbar_D4_D2_R_Cgminus(Nt,antiD4,antiD2,Cgminus,antiHQBR);
    QllBlock Bbar_D3_D1_0_Cgminus(Nt,antiD3,antiD1,Cgminus,antiHQB0);
    QllBlock Bbar_D2_D4_R_Cgminus(Nt,antiD2,antiD4,Cgminus,antiHQBR);
    QllBlock Bbar_D1_D3_0_Cgminus(Nt,antiD1,antiD3,Cgminus,antiHQB0);
    QllBlock Bbar_D4_D4_R_Cgminus(Nt,antiD4,antiD4,Cgminus,antiHQBR);
    QllBlock Bbar_D3_D3_0_Cgminus(Nt,antiD3,antiD3,Cgminus,antiHQB0);
    QllBlock Bbar_D2_D2_R_Cgminus(Nt,antiD2,antiD2,Cgminus,antiHQBR);
    QllBlock Bbar_D1_D1_0_Cgminus(Nt,antiD1,antiD1,Cgminus,antiHQB0);

    // Sigma_b^+ and (=) \Sigma_b^- baryon 2pt
    QDPIO::cout<<"  Contracting Sigma_b^+ \n";
    multi1d<DComplex> sigmabplusJ1m1, sigmabplusJ1m1a, sigmabplusJ1m1b;
    multi1d<DComplex> sigmabplusJ1m0,sigmabplusJ1m0a,sigmabplusJ1m0b;
    multi1d<DComplex> sigmabplusJ1mneg1,sigmabplusJ1mneg1a,sigmabplusJ1mneg1b;

    sigmabplusJ1m1.resize(length); sigmabplusJ1m1a.resize(length); sigmabplusJ1m1b.resize(length);
    sigmabplusJ1m0.resize(length);sigmabplusJ1m0a.resize(length);sigmabplusJ1m0b.resize(length);
    sigmabplusJ1mneg1.resize(length);sigmabplusJ1mneg1a.resize(length);sigmabplusJ1mneg1b.resize(length);

    sigmabplusJ1m1a = sigmabpluscontract(B_U1_U1_0_Cgplus, Cgplus);
    sigmabplusJ1m1b = sigmabpluscontract(B_U2_U2_R_Cgplus, Cgplus);
    sigmabplusJ1m0a = sigmabpluscontract(B_U1_U1_0_Cg3, Cg3);
    sigmabplusJ1m0b = sigmabpluscontract(B_U2_U2_R_Cg3, Cg3);
    sigmabplusJ1mneg1a = sigmabpluscontract(B_U1_U1_0_Cgminus, Cgminus);
    sigmabplusJ1mneg1b = sigmabpluscontract(B_U2_U2_R_Cgminus, Cgminus);
    for (int t=0; t<length;t++){
      sigmabplusJ1m1[t] = 0.5 * (sigmabplusJ1m1a[t]+sigmabplusJ1m1b[t]);
      sigmabplusJ1m0[t] = 0.5 * (sigmabplusJ1m0a[t]+sigmabplusJ1m0b[t]);
      sigmabplusJ1mneg1[t] = 0.5 * (sigmabplusJ1mneg1a[t]+sigmabplusJ1mneg1b[t]);
    };


    // Sigma_b^0 2pt : average over two positions
    QDPIO::cout<<"  Contracting Sigma_b^0 \n";
    multi1d<DComplex> sigmabzeroJ1m1, sigmabzeroJ1m0, sigmabzeroJ1mneg1;

    sigmabzeroJ1m1 = lambdabcontract(B_U1_D1_0_Cgplus, Cgplus);
    sigmabzeroJ1m1 += lambdabcontract(B_U2_D2_R_Cgplus, Cgplus);
    sigmabzeroJ1m1 *= 0.5;

    sigmabzeroJ1m0 = lambdabcontract(B_U1_D1_0_Cg3, Cg3);
    sigmabzeroJ1m0 += lambdabcontract(B_U2_D2_R_Cg3, Cg3);
    sigmabzeroJ1m0 *= 0.5;

    sigmabzeroJ1mneg1 = lambdabcontract(B_U1_D1_0_Cgminus, Cgminus);
    sigmabzeroJ1mneg1 += lambdabcontract(B_U2_D2_R_Cgminus, Cgminus);
    sigmabzeroJ1mneg1 *= 0.5;


    // Sigma_b^+ \Sigma_b^+ contractions: I=I_3=2,  J=0,1,2
    QDPIO::cout<<"  Contracting Sigma_b^+--Sigma_b^+ \n";
    multi2d<DComplex> sigmabplussigmabplus;
    sigmabplussigmabplus.resize(9,length); // Store all the spin states  
    sigmabplussigmabplus = c5J2corr(// Spin up blocks
      B_U1_U1_0_Cgplus,  B_U1_U3_0_Cgplus,
      B_U2_U2_R_Cgplus,  B_U2_U4_R_Cgplus,
      B_U3_U1_0_Cgplus,  B_U3_U3_0_Cgplus,
      B_U4_U2_R_Cgplus,  B_U4_U4_R_Cgplus,
      // Spin zero blocks
      B_U1_U1_0_Cg3,  B_U1_U3_0_Cg3,
      B_U2_U2_R_Cg3,  B_U2_U4_R_Cg3,
      B_U3_U1_0_Cg3,  B_U3_U3_0_Cg3,
      B_U4_U2_R_Cg3,  B_U4_U4_R_Cg3,
      // Spin down blocks
      B_U1_U1_0_Cgminus,  B_U1_U3_0_Cgminus,
      B_U2_U2_R_Cgminus,  B_U2_U4_R_Cgminus,
      B_U3_U1_0_Cgminus,  B_U3_U3_0_Cgminus,
      B_U4_U2_R_Cgminus,  B_U4_U4_R_Cgminus);

    // Sigma_b^+ \Sigma_b^0 contractions: I_3=1 (I=1,2),  J=0,1,2
    QDPIO::cout<<"  Contracting Sigma_b^+--Sigma_b^0 \n";
    multi2d<DComplex> sigmabplussigmabzero;
    sigmabplussigmabzero.resize(9,length); // Store all the spin states  
    sigmabplussigmabzero = c4J2corr(// Spin up blocks
      B_U1_D1_0_Cgplus,  B_U2_U2_R_Cgplus,
      B_U2_U4_R_Cgplus,  B_U3_D1_0_Cgplus, 
      B_U4_U2_R_Cgplus,
      // Spin zero blocks
      B_U1_D1_0_Cg3,  B_U2_U2_R_Cg3,
      B_U2_U4_R_Cg3,  B_U3_D1_0_Cg3, 
      B_U4_U2_R_Cg3,
      // Spin down blocks
      B_U1_D1_0_Cgminus,  B_U2_U2_R_Cgminus,
      B_U2_U4_R_Cgminus,  B_U3_D1_0_Cgminus, 
      B_U4_U2_R_Cgminus);

  
    // Sigma_b^+ \Sigma_b^- contractions: I_3=0 (I=0,1,2),  J=0,1,2
    QDPIO::cout<<"  Contracting Sigma_b^+--Sigma_b^- \n";
    multi2d<DComplex> sigmabplussigmabminus;
    sigmabplussigmabminus.resize(9,length); // Store all the spin states  
    sigmabplussigmabminus = c6J2corr(B_D1_D1_0_Cgplus, B_U2_U2_R_Cgplus,
				     B_D1_D1_0_Cg3, B_U2_U2_R_Cg3,
				     B_D1_D1_0_Cgminus, B_U2_U2_R_Cgminus);

    // Sigma_b^+ \Lambda_b contractions: I=I_3=1, J=1
    QDPIO::cout<<"  Contracting Sigma_b^+--Lambda_b \n";
    multi1d<DComplex> sigmabpluslambdabJ1m1, sigmabpluslambdabJ1m0, sigmabpluslambdabJ1mneg1;
    sigmabpluslambdabJ1m1.resize(length);
    sigmabpluslambdabJ1m0.resize(length);
    sigmabpluslambdabJ1mneg1.resize(length);

    sigmabpluslambdabJ1m1 = c4contract(B_U1_D1_0_Cg5, B_U2_U2_R_Cgplus,
				       B_U2_U4_R_Cgplus, B_U3_D1_0_Cg5, 
				       B_U4_U2_R_Cgplus,
				       Cgplus, Cg5);
    sigmabpluslambdabJ1m0 = c4contract(B_U1_D1_0_Cg5, B_U2_U2_R_Cg3,
				       B_U2_U4_R_Cg3, B_U3_D1_0_Cg5, 
				       B_U4_U2_R_Cg3,
				       Cg3, Cg5);
    sigmabpluslambdabJ1mneg1 = c4contract(B_U1_D1_0_Cg5, B_U2_U2_R_Cgminus,
					  B_U2_U4_R_Cgminus, B_U3_D1_0_Cg5, 
					  B_U4_U2_R_Cgminus,
					  Cgminus, Cg5);
  
    // Sigma_b^+ \Lambda_b transition contractions: I=I_3=1, J=1
    // Since there are two different particles, we  symmetrise over 
    // A(r) B(0) -->  A(r) B(0)  and B(r) A(0) -->  B(r) A(0)  
    QDPIO::cout<<"  Contracting Sigma_b^+--Lambda_b transition \n";
    multi1d<DComplex> sigmabpluslambdabtransitionJ1m1, sigmabpluslambdabtransitionJ1m0, sigmabpluslambdabtransitionJ1mneg1;
    sigmabpluslambdabtransitionJ1m1.resize(length);
    sigmabpluslambdabtransitionJ1m0.resize(length);
    sigmabpluslambdabtransitionJ1mneg1.resize(length);

    sigmabpluslambdabtransitionJ1m1 = c7contract(B_U2_D4_R_Cg5, B_U1_U3_0_Cgplus, B_U4_D4_R_Cg5, 
						 B_U3_U1_0_Cgplus, B_U3_U3_0_Cgplus,
						 Cg5, Cgplus);
    sigmabpluslambdabtransitionJ1m0 = c7contract(B_U2_D4_R_Cg5, B_U1_U3_0_Cg3, B_U4_D4_R_Cg5, 
						 B_U3_U1_0_Cg3, B_U3_U3_0_Cg3,
						 Cg5, Cg3);
    sigmabpluslambdabtransitionJ1mneg1 = c7contract(B_U2_D4_R_Cg5, B_U1_U3_0_Cgminus, B_U4_D4_R_Cg5, 
						    B_U3_U1_0_Cgminus, B_U3_U3_0_Cgminus,
						    Cg5, Cgminus);
//    sigmabpluslambdabtransitionJ1m1 = c7contract(B_U1_D3_0_Cg5, B_U2_U4_R_Cgplus, B_U3_D3_0_Cg5, 
// 					       B_U4_U2_R_Cgplus, B_U4_U4_R_Cgplus,
// 					       Cg5, Cgplus);
//    sigmabpluslambdabtransitionJ1m0 = c7contract(B_U1_D3_0_Cg5, B_U2_U4_R_Cg3, B_U3_D3_0_Cg5, 
// 					       B_U4_U2_R_Cg3, B_U4_U4_R_Cg3,
// 					       Cg5, Cg3);
//    sigmabpluslambdabtransitionJ1mneg1 = c7contract(B_U1_D3_0_Cg5, B_U2_U4_R_Cgminus, B_U3_D3_0_Cg5, 
// 						  B_U4_U2_R_Cgminus, B_U4_U4_R_Cgminus,
// 						  Cg5, Cgminus);

    QDPIO::cout<<"  Contracting Sigma_b^+--antiSigma_b^- \n";
    // \overline{Sigma}_b^+ \Sigma_b^- contractions: I=I_3=1, J=0,1,2
    // Since there are two different particles, we  symmetrise over 
    // A(r) B(0) -->  A(r) B(0)  and B(r) A(0) -->  B(r) A(0)  
    multi2d<DComplex> sigmabplussigmabarbminusa,sigmabplussigmabarbminusb,sigmabplussigmabarbminus;
    sigmabplussigmabarbminusa.resize(9,length); // Store all the spin states  
    sigmabplussigmabarbminusb.resize(9,length); // Store all the spin states  
    sigmabplussigmabarbminus.resize(9,length); // Store all the spin states  
    sigmabplussigmabarbminusa = c6J2corr(Bbar_D1_D1_0_Cgplus, B_U2_U2_R_Cgplus,
					 Bbar_D1_D1_0_Cg3, B_U2_U2_R_Cg3,
					 Bbar_D1_D1_0_Cgminus, B_U2_U2_R_Cgminus);
    sigmabplussigmabarbminusb = c6J2corr(B_U1_U1_0_Cgplus, Bbar_D2_D2_R_Cgplus,
					 B_U1_U1_0_Cg3, Bbar_D2_D2_R_Cg3,
					 B_U1_U1_0_Cgminus, Bbar_D2_D2_R_Cgminus);
    for (int i=0;i<9;i++){   for (int t=0; t<length;t++){
	sigmabplussigmabarbminus[i][t] = 0.5 * (sigmabplussigmabarbminusa[i][t] +  sigmabplussigmabarbminusb[i][t]);
      }}

    // MESON-BARYON CONTRACTIONS

    QDPIO::cout<<"  Contracting Bu--Lambda_b \n";
    // Bu \Lambda_b contractions: I=I_3=1/2, J=1/2, m=1/2
    // Since there are two different particles, we  symmetrise over 
    // A(r) B(0) -->  A(r) B(0)  and B(r) A(0) -->  B(r) A(0)  
    multi1d<DComplex> bulambdabJhalfmhalf;
    bulambdabJhalfmhalf.resize(length);

    bulambdabJhalfmhalf = d1contract(B_U1_D1_0_Cg5, B_U3_D1_0_Cg5, 
				     H_U2_Gup_R, H_U4_Gup_R,
				     Gup, Cg5);
    bulambdabJhalfmhalf += d1contract(B_U2_D2_R_Cg5, B_U4_D2_R_Cg5, 
				      H_U1_Gup_0, H_U3_Gup_0,
				      Gup, Cg5);
    bulambdabJhalfmhalf *= 0.5;

    multi1d<DComplex> bulambdabJhalfmneghalf;
    bulambdabJhalfmneghalf.resize(length);
    bulambdabJhalfmneghalf = d1contract(B_U1_D1_0_Cg5, B_U3_D1_0_Cg5, 
					H_U2_Gdown_R, H_U4_Gdown_R,
					Gdown, Cg5);
    bulambdabJhalfmneghalf += d1contract(B_U2_D2_R_Cg5, B_U4_D2_R_Cg5, 
					 H_U1_Gdown_0, H_U3_Gdown_0,
					 Gdown, Cg5);
    bulambdabJhalfmneghalf *= 0.5;

    QDPIO::cout<<"  Contracting Bu--Sigma_b^+ \n";
    // Bu \Sigma_b^+ contractions: I=I_3=3/2, J=1/2, 3/2
    // Since there are two different particles, we  symmetrise over 
    // A(r) B(0) -->  A(r) B(0)  and B(r) A(0) -->  B(r) A(0)  
    multi2d<DComplex> busigmabplus,busigmabplusa,busigmabplusb;
    busigmabplus.resize(6,length);
    busigmabplusa.resize(6,length);
    busigmabplusb.resize(6,length);
    busigmabplusa = d2J32corr(// Spin up blocks
      B_U1_U1_0_Cgplus, B_U3_U1_0_Cgplus, B_U1_U3_0_Cgplus, 
      H_U2_Gup_R, H_U4_Gup_R,
      // Spin zero blocks
      B_U1_U1_0_Cg3, B_U3_U1_0_Cg3, B_U1_U3_0_Cg3, 
      // Spin down blocks
      B_U1_U1_0_Cgminus, B_U3_U1_0_Cgminus, B_U1_U3_0_Cgminus, 
      H_U2_Gdown_R, H_U4_Gdown_R);
    busigmabplusb = d2J32corr(// Spin up blocks
      B_U2_U2_R_Cgplus, B_U4_U2_R_Cgplus, B_U2_U4_R_Cgplus, 
      H_U1_Gup_0, H_U3_Gup_0,
      // Spin zero blocks
      B_U2_U2_R_Cg3, B_U4_U2_R_Cg3, B_U2_U4_R_Cg3, 
      // Spin down blocks
      B_U2_U2_R_Cgminus, B_U4_U2_R_Cgminus, B_U2_U4_R_Cgminus, 
      H_U1_Gdown_0, H_U3_Gdown_0);
    for (int i=0;i<6;i++){   for (int t=0; t<length;t++){
	busigmabplus[i][t] = 0.5 *(busigmabplusa[i][t] + busigmabplusb[i][t]);
      }}

    QDPIO::cout<<"  Contracting Bd--Sigma_b^+ \n";
    // Bd \Sigma_b^+ contractions: I_3=1/2 (I=1/2, 3/2), J=1/2, 3/2
    // Since there are two different particles, we  symmetrise over 
    // A(r) B(0) -->  A(r) B(0)  and B(r) A(0) -->  B(r) A(0)  
    multi2d<DComplex> bdsigmabplus,bdsigmabplusa,bdsigmabplusb;
    bdsigmabplus.resize(6,length);
    bdsigmabplusa.resize(6,length);
    bdsigmabplusb.resize(6,length);
    bdsigmabplusa = d3J32corr(// Spin up blocks
      B_U1_U1_0_Cgplus,
      H_D2_Gup_R,
      // Spin zero blocks
      B_U1_U1_0_Cg3,
      // Spin down blocks
      B_U1_U1_0_Cgminus,
      H_D2_Gdown_R);
    bdsigmabplusb = d3J32corr(// Spin up blocks
      B_U2_U2_R_Cgplus,
      H_D1_Gup_0,
      // Spin zero blocks
      B_U2_U2_R_Cg3,
      // Spin down blocks
      B_U2_U2_R_Cgminus,
      H_D1_Gdown_0);
    for (int i=0;i<6;i++){   for (int t=0; t<length;t++){
	bdsigmabplus[i][t] = 0.5 *(bdsigmabplusa[i][t] + bdsigmabplusb[i][t]);
      }}

    QDPIO::cout<<"  Contracting Bbard--Sigma_b^+ \n";
    // Bbard \Sigma_b^+ (Q dbar Quu) contractions: I_3=1/2 (I=1/2, 3/2), J=1/2, 3/2
    // Since there are two different particles, we  symmetrise over 
    // A(r) B(0) -->  A(r) B(0)  and B(r) A(0) -->  B(r) A(0)  
    multi2d<DComplex> bbardsigmabplus,bbardsigmabplusa,bbardsigmabplusb;
    bbardsigmabplus.resize(6,length);
    bbardsigmabplusa.resize(6,length);
    bbardsigmabplusb.resize(6,length);
    bbardsigmabplusa = d3J32corr(// Spin up blocks
      B_U1_U1_0_Cgplus,
      Hbar_D2_Gup_R,
      // Spin zero blocks
      B_U1_U1_0_Cg3,
      // Spin down blocks
      B_U1_U1_0_Cgminus,
      Hbar_D2_Gdown_R);
    bbardsigmabplusb = d3J32corr(// Spin up blocks
      B_U2_U2_R_Cgplus,
      Hbar_D1_Gup_0,
      // Spin zero blocks
      B_U2_U2_R_Cg3,
      // Spin down blocks
      B_U2_U2_R_Cgminus,
      Hbar_D1_Gdown_0);
    for (int i=0;i<6;i++){   for (int t=0; t<length;t++){
	bbardsigmabplus[i][t] = 0.5 *(bbardsigmabplusa[i][t] + bbardsigmabplusb[i][t]);
      }}



    push(xml, xml_group);
    // Hadrons
    write(xml, "Bmes",bmes);
    write(xml, "Bmes0",bmes0);
    write(xml, "BmesR", bmesR);
    write(xml, "BmesJ0.5m0.5",bmesup);
    write(xml, "Bmes0J0.5m0.5",bmes0up);
    write(xml, "BmesRJ0.5m0.5", bmesRup);
    write(xml, "BmesJ0.5mneg0.5",bmesdown);
    write(xml, "Bmes0J0.5mneg0.5",bmes0down);
    write(xml, "BmesRJ0.5mneg0.5", bmesRdown);
    write(xml,"lambdab", lambdab);
    write(xml,"lambdab0", lambdaba);
    write(xml,"lambdabR", lambdabb);
    write(xml,"sigmabplusJ1m1", sigmabplusJ1m1);
    write(xml,"sigmabplusJ1m0", sigmabplusJ1m0);
    write(xml,"sigmabplusJ1mneg1", sigmabplusJ1mneg1);
    write(xml,"sigmabplusJ1m1_0", sigmabplusJ1m1a);
    write(xml,"sigmabplusJ1m0_0", sigmabplusJ1m0a);
    write(xml,"sigmabplusJ1mneg1_0", sigmabplusJ1mneg1a);
    write(xml,"sigmabplusJ1m1_R", sigmabplusJ1m1b);
    write(xml,"sigmabplusJ1m0_R", sigmabplusJ1m0b);
    write(xml,"sigmabplusJ1mneg1_R", sigmabplusJ1mneg1b);
    write(xml,"sigmabzeroJ1m1", sigmabzeroJ1m1);
    write(xml,"sigmabzeroJ1m0", sigmabzeroJ1m0);
    write(xml,"sigmabzeroJ1mneg1", sigmabzeroJ1mneg1);
    // Meson-meson  
    write(xml, "Bu_Bu_J0m0", BuBuJ0m0);
    write(xml, "Bu_Bu_J1m1", BuBuJ1m1);
    write(xml, "Bu_Bu_J1m0",BuBuJ1m0);
    write(xml, "Bu_Bu_J1mneg1", BuBuJ1mneg1);
    write(xml, "Bu_Bd_J0m0", BuBdJ0m0);
    write(xml, "Bu_Bd_J1m1", BuBdJ1m1);
    write(xml, "Bu_Bd_J1m0", BuBdJ1m0);
    write(xml, "Bu_Bd_J1mneg1", BuBdJ1mneg1);
    write(xml, "Bu_Bdbar_J0m0", BuBdbarJ0m0);
    write(xml, "Bu_Bdbar_J1m1", BuBdbarJ1m1);
    write(xml, "Bu_Bdbar_J1m0", BuBdbarJ1m0);
    write(xml, "Bu_Bdbar_J1mneg1", BuBdbarJ1mneg1);
    // Baryon-baryon
    write(xml, "lambdaB_lambdaB_J0m0", lambdablambdab);
    write(xml, "sigmabplus_sigmabplus_J2m2", sigmabplussigmabplus[0]);
    write(xml, "sigmabplus_sigmabplus_J2m1", sigmabplussigmabplus[1]);
    write(xml, "sigmabplus_sigmabplus_J2m0", sigmabplussigmabplus[2]);
    write(xml, "sigmabplus_sigmabplus_J2mneg1", sigmabplussigmabplus[3]);
    write(xml, "sigmabplus_sigmabplus_J2mneg2", sigmabplussigmabplus[4]);
    write(xml, "sigmabplus_sigmabplus_J1m1", sigmabplussigmabplus[5]);
    write(xml, "sigmabplus_sigmabplus_J1m0", sigmabplussigmabplus[6]);
    write(xml, "sigmabplus_sigmabplus_J1mneg1", sigmabplussigmabplus[7]);
    write(xml, "sigmabplus_sigmabplus_J0m0", sigmabplussigmabplus[8]);

    write(xml, "sigmabplus_sigmabzero_J2m2", sigmabplussigmabzero[0]);
    write(xml, "sigmabplus_sigmabzero_J2m1", sigmabplussigmabzero[1]);
    write(xml, "sigmabplus_sigmabzero_J2m0", sigmabplussigmabzero[2]);
    write(xml, "sigmabplus_sigmabzero_J2mneg1", sigmabplussigmabzero[3]);
    write(xml, "sigmabplus_sigmabzero_J2mneg2", sigmabplussigmabzero[4]);
    write(xml, "sigmabplus_sigmabzero_J1m1", sigmabplussigmabzero[5]);
    write(xml, "sigmabplus_sigmabzero_J1m0", sigmabplussigmabzero[6]);
    write(xml, "sigmabplus_sigmabzero_J1mneg1", sigmabplussigmabzero[7]);
    write(xml, "sigmabplus_sigmabzero_J0m0", sigmabplussigmabzero[8]);

    write(xml, "sigmabplus_sigmabminus_J2m2", sigmabplussigmabminus[0]);
    write(xml, "sigmabplus_sigmabminus_J2m1", sigmabplussigmabminus[1]);
    write(xml, "sigmabplus_sigmabminus_J2m0", sigmabplussigmabminus[2]);
    write(xml, "sigmabplus_sigmabminus_J2mneg1", sigmabplussigmabminus[3]);
    write(xml, "sigmabplus_sigmabminus_J2mneg2", sigmabplussigmabminus[4]);
    write(xml, "sigmabplus_sigmabminus_J1m1", sigmabplussigmabminus[5]);
    write(xml, "sigmabplus_sigmabminus_J1m0", sigmabplussigmabminus[6]);
    write(xml, "sigmabplus_sigmabminus_J1mneg1", sigmabplussigmabminus[7]);
    write(xml, "sigmabplus_sigmabminus_J0m0", sigmabplussigmabminus[8]);

    write(xml, "sigmabplus_antisigmabminus_J2m2", sigmabplussigmabarbminus[0]);
    write(xml, "sigmabplus_antisigmabminus_J2m1", sigmabplussigmabarbminus[1]);
    write(xml, "sigmabplus_antisigmabminus_J2m0", sigmabplussigmabarbminus[2]);
    write(xml, "sigmabplus_antisigmabminus_J2mneg1", sigmabplussigmabarbminus[3]);
    write(xml, "sigmabplus_antisigmabminus_J2mneg2", sigmabplussigmabarbminus[4]);
    write(xml, "sigmabplus_antisigmabminus_J1m1", sigmabplussigmabarbminus[5]);
    write(xml, "sigmabplus_antisigmabminus_J1m0", sigmabplussigmabarbminus[6]);
    write(xml, "sigmabplus_antisigmabminus_J1mneg1", sigmabplussigmabarbminus[7]);
    write(xml, "sigmabplus_antisigmabminus_J0m0", sigmabplussigmabarbminus[8]);

    write(xml, "sigmabplus_lambda_J1m1", sigmabpluslambdabJ1m1);
    write(xml, "sigmabplus_lambda_J1m0", sigmabpluslambdabJ1m0);
    write(xml, "sigmabplus_lambda_J1mneg1", sigmabpluslambdabJ1mneg1);

    write(xml, "sigmabplus_lambda_transition_J1m1", sigmabpluslambdabtransitionJ1m1);
    write(xml, "sigmabplus_lambda_transition_J1m0", sigmabpluslambdabtransitionJ1m0);
    write(xml, "sigmabplus_lambda_transition_J1mneg1", sigmabpluslambdabtransitionJ1mneg1);

    // Meson-baryon
    write(xml, "bu_lambdab_J0.5m0.5", bulambdabJhalfmhalf);
    write(xml, "bu_lambdab_J0.5mneg0.5", bulambdabJhalfmneghalf);

    write(xml, "bu_sigmabplus_J1.5m1.5", busigmabplus[0]);
    write(xml, "bu_sigmabplus_J1.5m0.5", busigmabplus[1]);
    write(xml, "bu_sigmabplus_J1.5mneg0.5", busigmabplus[2]);
    write(xml, "bu_sigmabplus_J1.5mneg1.5", busigmabplus[3]);
    write(xml, "bu_sigmabplus_J0.5m0.5", busigmabplus[4]);
    write(xml, "bu_sigmabplus_J0.5mneg0.5", busigmabplus[5]);

    write(xml, "bd_sigmabplus_J1.5m1.5", bdsigmabplus[0]);
    write(xml, "bd_sigmabplus_J1.5m0.5", bdsigmabplus[1]);
    write(xml, "bd_sigmabplus_J1.5mneg0.5", bdsigmabplus[2]);
    write(xml, "bd_sigmabplus_J1.5mneg1.5", bdsigmabplus[3]);
    write(xml, "bd_sigmabplus_J0.5m0.5", bdsigmabplus[4]);
    write(xml, "bd_sigmabplus_J0.5mneg0.5", bdsigmabplus[5]);

    write(xml, "bbard_sigmabplus_J1.5m1.5", bbardsigmabplus[0]);
    write(xml, "bbard_sigmabplus_J1.5m0.5", bbardsigmabplus[1]);
    write(xml, "bbard_sigmabplus_J1.5mneg0.5", bbardsigmabplus[2]);
    write(xml, "bbard_sigmabplus_J1.5mneg1.5", bbardsigmabplus[3]);
    write(xml, "bbard_sigmabplus_J0.5m0.5", bbardsigmabplus[4]);
    write(xml, "bbard_sigmabplus_J0.5mneg0.5", bbardsigmabplus[5]);

    pop(xml);

    END_CODE();
  }






  multiNd<DComplex> HBQfunc(const multi1d<ColorMatrix>&  HQ)
  {
    int length = HQ.size();
    multiNd<DComplex> result;
    multi1d<int> Barray;
    Barray.resize(4); Barray= Nc; Barray[0]=length; 
    result.resize(Barray);
    result = 0;
    DComplex thisres;
    for (int t=0; t<length; t++){
      Barray[0]=t;
      for (int aa=0; aa<Nc; aa++){
	Barray[1]=aa;
	for (int bb=0; bb<Nc; bb++){
	  Barray[2]=bb;
	  for (int c=0; c<Nc; c++){
	    Barray[3]=c;
	    thisres = 0;
	    for (int cc=0;cc<Nc; cc++) // summed colour index
	      thisres += antiSymTensor3d(aa,bb,cc) * peekColor(HQ[t],cc,c);
	    result[Barray] = thisres;
	  }
	}
      } 
    }
    return result;
  }

  multiNd<DComplex> antiHBQfunc(const multi1d<ColorMatrix>& HQ)
  {
    int length = HQ.size();
    multiNd<DComplex> result;
    multi1d<int> Barray;
    Barray.resize(4); Barray= Nc; Barray[0]=length; 
    result.resize(Barray);
    result = 0;
    DComplex thisres;
    for (int t=0; t<length; t++){
      Barray[0]=t;
      for (int aa=0; aa<Nc; aa++){
	Barray[1]=aa;
	for (int bb=0; bb<Nc; bb++){
	  Barray[2]=bb;
	  for (int c=0; c<Nc; c++){
	    Barray[3]=c;
	    thisres = 0;
	    for (int cc=0;cc<Nc; cc++) // summed colour index
	      thisres += antiSymTensor3d(aa,bb,cc) * conj(peekColor(HQ[t],c,cc)); // conjugate transpose
	    result[Barray] = thisres;
	  }
	}
      } 
    }
    return result;
  }



  multi1d<DComplex> c1contract(const QllBlock& BzU1zD1z0zCg5, const QllBlock& BzU1zD3z0zCg5, 
			       const QllBlock& BzU2zD2zRzCg5, const QllBlock& BzU2zD4zRzCg5, 
			       const QllBlock& BzU3zD1z0zCg5, const QllBlock& BzU3zD3z0zCg5, 
			       const QllBlock& BzU4zD2zRzCg5, const QllBlock& BzU4zD4zRzCg5,
			       const SpinMatrix& S1, const SpinMatrix& S2)
  // Contractions for Lambda_b (R) \Lambda_b (0) 
  // spin matrix S1 is at R and S2 is at 0
  // 

  {
    int length=BzU1zD1z0zCg5.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde, S2tilde; // \gamma_0 S_i^\dagger\gamma_0

    S1tilde = Gamma(8) * adj(S1) * Gamma(8);
    S2tilde = Gamma(8) * adj(S2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res, c3res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,is,js) * peekSpin(S2tilde,ks,ls);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			for (int mc=0; mc<Nc; mc++){  // summed colour index
			  for (int nc=0; nc<Nc; nc++){  // summed colour index
			    c3res =0;
			    if (jc != mc && jc != nc && mc != nc) 
			    {
			      c3res = -BzU1zD1z0zCg5(t,jc,mc,nc,ks,ls)*BzU2zD2zRzCg5(t,ic,kc,lc,is,js) + 
				BzU1zD3z0zCg5(t,jc,kc,nc,ks,js)*BzU2zD4zRzCg5(t,ic,mc,lc,is,ls) + 
				BzU3zD1z0zCg5(t,ic,mc,nc,is,ls)*BzU4zD2zRzCg5(t,jc,kc,lc,ks,js) - 
				BzU3zD3z0zCg5(t,ic,kc,nc,is,js)*BzU4zD4zRzCg5(t,jc,mc,lc,ks,ls);


			      c3res *= antiSymTensor3d(jc,mc,nc);
			    }
			    c2res += c3res;
			  } // colour nc
			} //colour mc
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }



  multi1d<DComplex> c4contract(const QllBlock& BzU1zD1z0zCgjj, const QllBlock& BzU2zU2zRzCgii,
			       const QllBlock& BzU2zU4zRzCgii, const QllBlock& BzU3zD1z0zCgjj, 
			       const QllBlock& BzU4zU2zRzCgii,
			       const SpinMatrix& S1, const SpinMatrix& S2)
  // Contractions for Sigma_b^+ (R) \Sigma_b^0 (0) 
  // and  Sigma_b^+ (R) \Lambda_b (0)
  // spin matrix S1 is at R and S2 is at 0
  // 

  {
    int length=BzU1zD1z0zCgjj.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde, S2tilde; // \gamma_0 S_i^\dagger\gamma_0

    S1tilde = Gamma(8) * adj(S1) * Gamma(8);
    S2tilde = Gamma(8) * adj(S2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res, c3res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,is,js) * peekSpin(S2tilde,ks,ls);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			for (int mc=0; mc<Nc; mc++){  // summed colour index
			  for (int nc=0; nc<Nc; nc++){  // summed colour index
			    c3res =0;
			    if (jc != mc && jc != nc && mc != nc) 
			    {
			      c3res = -BzU1zD1z0zCgjj(t,jc,mc,nc,ks,ls)*BzU2zU2zRzCgii(t,ic,kc,lc,is,js) + 
				BzU1zD1z0zCgjj(t,jc,mc,nc,ks,ls)*BzU2zU2zRzCgii(t,kc,ic,lc,js,is) - 
				BzU2zU4zRzCgii(t,kc,jc,lc,js,ks)*BzU3zD1z0zCgjj(t,ic,mc,nc,is,ls) + 
				BzU2zU4zRzCgii(t,ic,jc,lc,is,ks)*BzU3zD1z0zCgjj(t,kc,mc,nc,js,ls) - 
				BzU3zD1z0zCgjj(t,kc,mc,nc,js,ls)*BzU4zU2zRzCgii(t,jc,ic,lc,ks,is) + 
				BzU3zD1z0zCgjj(t,ic,mc,nc,is,ls)*BzU4zU2zRzCgii(t,jc,kc,lc,ks,js);
				
			      c3res *= antiSymTensor3d(jc,mc,nc);
			    }
			    c2res += c3res;
			  } // colour nc
			} //colour mc
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }


  multi1d<DComplex> c5contract(const QllBlock& BzU1zU1z0zCgjj, const QllBlock& BzU1zU3z0zCgjj,
			       const QllBlock& BzU2zU2zRzCgii, const QllBlock& BzU2zU4zRzCgii,
			       const QllBlock& BzU3zU1z0zCgjj, const QllBlock& BzU3zU3z0zCgjj,
			       const QllBlock& BzU4zU2zRzCgii, const QllBlock& BzU4zU4zRzCgii,
			       const SpinMatrix& S1, const SpinMatrix& S2)
  // Contractions for Sigma_b^+ \Sigma_b^+
  // spin matrix S1 is at R and S2 is at 0
  // 

  {
    int length=BzU1zU1z0zCgjj.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde, S2tilde; // \gamma_0 S_i^\dagger\gamma_0

    S1tilde = Gamma(8) * adj(S1) * Gamma(8);
    S2tilde = Gamma(8) * adj(S2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res, c3res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,is,js) * peekSpin(S2tilde,ks,ls);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			for (int mc=0; mc<Nc; mc++){  // summed colour index
			  for (int nc=0; nc<Nc; nc++){  // summed colour index
			    c3res =0;
			    if (jc != mc && jc != nc && mc != nc) 
			    {
			      c3res = - BzU1zU1z0zCgjj(t,jc,mc,nc,ks,ls)*BzU2zU2zRzCgii(t,ic,kc,lc,is,js) + 
				BzU1zU1z0zCgjj(t,mc,jc,nc,ls,ks)*BzU2zU2zRzCgii(t,ic,kc,lc,is,js) + 
				BzU1zU1z0zCgjj(t,jc,mc,nc,ks,ls)*BzU2zU2zRzCgii(t,kc,ic,lc,js,is) - 
				BzU1zU1z0zCgjj(t,mc,jc,nc,ls,ks)*BzU2zU2zRzCgii(t,kc,ic,lc,js,is) - 
				BzU1zU3z0zCgjj(t,mc,kc,nc,ls,js)*BzU2zU4zRzCgii(t,ic,jc,lc,is,ks) + 
				BzU1zU3z0zCgjj(t,jc,kc,nc,ks,js)*BzU2zU4zRzCgii(t,ic,mc,lc,is,ls) + 
				BzU1zU3z0zCgjj(t,mc,ic,nc,ls,is)*BzU2zU4zRzCgii(t,kc,jc,lc,js,ks) - 
				BzU1zU3z0zCgjj(t,jc,ic,nc,ks,is)*BzU2zU4zRzCgii(t,kc,mc,lc,js,ls) + 
				BzU2zU4zRzCgii(t,kc,mc,lc,js,ls)*BzU3zU1z0zCgjj(t,ic,jc,nc,is,ks) - 
				BzU2zU4zRzCgii(t,kc,jc,lc,js,ks)*BzU3zU1z0zCgjj(t,ic,mc,nc,is,ls) - 
				BzU2zU4zRzCgii(t,ic,mc,lc,is,ls)*BzU3zU1z0zCgjj(t,kc,jc,nc,js,ks) + 
				BzU2zU4zRzCgii(t,ic,jc,lc,is,ks)*BzU3zU1z0zCgjj(t,kc,mc,nc,js,ls) + 
				BzU1zU3z0zCgjj(t,mc,kc,nc,ls,js)*BzU4zU2zRzCgii(t,jc,ic,lc,ks,is) - 
				BzU3zU1z0zCgjj(t,kc,mc,nc,js,ls)*BzU4zU2zRzCgii(t,jc,ic,lc,ks,is) - 
				BzU1zU3z0zCgjj(t,mc,ic,nc,ls,is)*BzU4zU2zRzCgii(t,jc,kc,lc,ks,js) + 
				BzU3zU1z0zCgjj(t,ic,mc,nc,is,ls)*BzU4zU2zRzCgii(t,jc,kc,lc,ks,js) - 
				BzU1zU3z0zCgjj(t,jc,kc,nc,ks,js)*BzU4zU2zRzCgii(t,mc,ic,lc,ls,is) + 
				BzU3zU1z0zCgjj(t,kc,jc,nc,js,ks)*BzU4zU2zRzCgii(t,mc,ic,lc,ls,is) + 
				BzU1zU3z0zCgjj(t,jc,ic,nc,ks,is)*BzU4zU2zRzCgii(t,mc,kc,lc,ls,js) - 
				BzU3zU1z0zCgjj(t,ic,jc,nc,is,ks)*BzU4zU2zRzCgii(t,mc,kc,lc,ls,js) - 
				BzU3zU3z0zCgjj(t,ic,kc,nc,is,js)*BzU4zU4zRzCgii(t,jc,mc,lc,ks,ls) + 
				BzU3zU3z0zCgjj(t,kc,ic,nc,js,is)*BzU4zU4zRzCgii(t,jc,mc,lc,ks,ls) + 
				BzU3zU3z0zCgjj(t,ic,kc,nc,is,js)*BzU4zU4zRzCgii(t,mc,jc,lc,ls,ks) - 
				BzU3zU3z0zCgjj(t,kc,ic,nc,js,is)*BzU4zU4zRzCgii(t,mc,jc,lc,ls,ks);
				
				
			      c3res *= antiSymTensor3d(jc,mc,nc);
			    }
			    c2res += c3res;
			  } // colour nc
			} //colour mc
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }



  multi1d<DComplex> c6contract(const QllBlock& BzD1zD1z0zCgjj, const QllBlock& BzU2zU2zRzCgii,
			       const SpinMatrix& S1, const SpinMatrix& S2)
  // Contractions for Sigma_b^+ \Sigma_b^- and  Sigma_b^+ \overline{\Sigma}_b^
  // spin matrix S1 is at R and S2 is at 0
  // 

  {
    int length=BzD1zD1z0zCgjj.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde, S2tilde; // \gamma_0 S_i^\dagger\gamma_0

    S1tilde = Gamma(8) * adj(S1) * Gamma(8);
    S2tilde = Gamma(8) * adj(S2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res, c3res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,is,js) * peekSpin(S2tilde,ks,ls);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			for (int mc=0; mc<Nc; mc++){  // summed colour index
			  for (int nc=0; nc<Nc; nc++){  // summed colour index
			    c3res =0;
			    if (jc != mc && jc != nc && mc != nc) 
			    {
			      c3res = -(BzD1zD1z0zCgjj(t,jc,mc,nc,ks,ls)*BzU2zU2zRzCgii(t,ic,kc,lc,is,js)) + 
				BzD1zD1z0zCgjj(t,mc,jc,nc,ls,ks)*BzU2zU2zRzCgii(t,ic,kc,lc,is,js) + 
				BzD1zD1z0zCgjj(t,jc,mc,nc,ks,ls)*BzU2zU2zRzCgii(t,kc,ic,lc,js,is) - 
				BzD1zD1z0zCgjj(t,mc,jc,nc,ls,ks)*BzU2zU2zRzCgii(t,kc,ic,lc,js,is);
				
			      c3res *= antiSymTensor3d(jc,mc,nc);
			    }
			    c2res += c3res;
			  } // colour nc
			} //colour mc
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }


  multi1d<DComplex> c7contract(const QllBlock& BzU1zD3z0zCG5, const QllBlock& BzU2zU4zRzCGii, const QllBlock& BzU3zD3z0zCG5, 
			       const QllBlock& BzU4zU2zRzCGii, const QllBlock& BzU4zU4zRzCGii,
			       const SpinMatrix& S1, const SpinMatrix& S2)
  // Contractions for Sigma_b^+ \Lambda_b transition
  // spin matrix S1 is CG5 and S2 is at CGi
  // 

  {
    int length=BzU1zD3z0zCG5.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde, S2tilde; // \gamma_0 S_i^\dagger\gamma_0

    S1tilde = Gamma(8) * adj(S1) * Gamma(8);
    S2tilde = Gamma(8) * adj(S2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res, c3res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,is,js) * peekSpin(S2tilde,ks,ls);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			for (int mc=0; mc<Nc; mc++){  // summed colour index
			  for (int nc=0; nc<Nc; nc++){  // summed colour index
			    c3res =0;
			    if (jc != mc && jc != nc && mc != nc) 
			    {
			      c3res = 
				BzU1zD3z0zCG5(t,mc,kc,nc,ls,js)*BzU2zU4zRzCGii(t,ic,jc,lc,is,ks) - 
				BzU1zD3z0zCG5(t,jc,kc,nc,ks,js)*BzU2zU4zRzCGii(t,ic,mc,lc,is,ls) - 
				BzU1zD3z0zCG5(t,mc,kc,nc,ls,js)*BzU4zU2zRzCGii(t,jc,ic,lc,ks,is) + 
				BzU1zD3z0zCG5(t,jc,kc,nc,ks,js)*BzU4zU2zRzCGii(t,mc,ic,lc,ls,is) + 
				BzU3zD3z0zCG5(t,ic,kc,nc,is,js)*BzU4zU4zRzCGii(t,jc,mc,lc,ks,ls) - 
				BzU3zD3z0zCG5(t,ic,kc,nc,is,js)*BzU4zU4zRzCGii(t,mc,jc,lc,ls,ks);
				
			      c3res *= antiSymTensor3d(jc,mc,nc);
			    }
			    c2res += c3res;
			  } // colour nc
			} //colour mc
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }


  multi1d<DComplex> d1contract(const QllBlock& BzU1zD1z0zCG5, const QllBlock& BzU3zD1z0zCG5, 
			       const HeavyMesonBlock& HzU2zRzG5, const HeavyMesonBlock& HzU4zRzG5,
			       const SpinMatrix& mesonS1, const SpinMatrix& baryonS2)
  // Contractions for Bu (R)  Lambda_b (0)
  // spin matrix S1 is at R and S2 is at 0
  // 

  {
    int length=BzU1zD1z0zCG5.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde,S2tilde; // \gamma_0 S_i^\dagger\gamma_0
    S1tilde = mesonS1;
    //  S1tilde = adj(mesonS1);
    S2tilde = Gamma(8) * adj(baryonS2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,ks,ls) * peekSpin(S2tilde,is,js);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			c2res += -BzU1zD1z0zCG5(t,ic,kc,lc,is,js)*HzU2zRzG5(t,jc,jc,ls,ks) +
			  BzU3zD1z0zCG5(t,jc,kc,lc,ks,js)*HzU4zRzG5(t,jc,ic,ls,is);
			  			
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }



  multi1d<DComplex> d2contract(const QllBlock& BzU1zU1z0zCGi, const QllBlock& BzU3zU1z0zCGi, const QllBlock& BzU1zU3z0zCGi, 
			       const HeavyMesonBlock& HzU2zRzG5, const HeavyMesonBlock& HzU4zRzG5,
			       const SpinMatrix& mesonS1, const SpinMatrix& baryonS2)
  // Contractions for Bu (R)  Sigma_b ^+(0)
  // spin matrix S1 is at R and S2 is at 0
  // 

  {
    int length=BzU1zU1z0zCGi.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde,S2tilde; // \gamma_0 S_i^\dagger\gamma_0
    S1tilde = mesonS1;
    S2tilde = Gamma(8) * adj(baryonS2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,ks,ls) * peekSpin(S2tilde,is,js);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			c2res += -BzU1zU1z0zCGi(t,ic,kc,lc,is,js)*HzU2zRzG5(t,jc,jc,ls,ks) + 
			  BzU1zU1z0zCGi(t,kc,ic,lc,js,is)*HzU2zRzG5(t,jc,jc,ls,ks) - 
			  BzU1zU3z0zCGi(t,kc,jc,lc,js,ks)*HzU4zRzG5(t,jc,ic,ls,is) + 
			  BzU3zU1z0zCGi(t,jc,kc,lc,ks,js)*HzU4zRzG5(t,jc,ic,ls,is) + 
			  BzU1zU3z0zCGi(t,ic,jc,lc,is,ks)*HzU4zRzG5(t,jc,kc,ls,js) - 
			  BzU3zU1z0zCGi(t,jc,ic,lc,ks,is)*HzU4zRzG5(t,jc,kc,ls,js);
			
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }


  multi1d<DComplex> d3contract(const QllBlock& BzU1zU1z0zCGi, const HeavyMesonBlock& HzD2zRzG5, 
			       const SpinMatrix& mesonS1, const SpinMatrix& baryonS2)
  // Contractions for Bd (R)  Sigma_b ^+(0)
  // spin matrix S1 is at R and S2 is at 0
  // 

  {
    int length=BzU1zU1z0zCGi.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde,S2tilde; // \gamma_0 S_i^\dagger\gamma_0
    S1tilde = mesonS1;
    S2tilde = Gamma(8) * adj(baryonS2) * Gamma(8);

    DComplex tmpSpin, c2res, c1res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){ // summed spin index
	for (int js=0; js<Nd; js++){  // summed spin index
	  for (int ks=0; ks<Nd; ks++){  // summed spin index
	    for (int ls=0; ls<Nd; ls++){  // summed spin index
	      tmpSpin = peekSpin(S1tilde,ks,ls) * peekSpin(S2tilde,is,js);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){ // summed colour index
		for (int kc=0; kc<Nc; kc++){  // summed colour index
		  for (int lc=0; lc<Nc; lc++){  // summed colour index
		    if (ic != kc && ic != lc && kc != lc) 
		    {
		      c2res =0;
		      for (int jc=0; jc<Nc; jc++){  // summed colour index
			c2res += (BzU1zU1z0zCGi(t,kc,ic,lc,js,is)-BzU1zU1z0zCGi(t,ic,kc,lc,is,js))*HzD2zRzG5(t,jc,jc,ls,ks);
		      } // colour jc
		      c1res += c2res * antiSymTensor3d(ic,kc,lc);
		    } // if epsilon
		  } // colour lc
		} // colour kc
	      } // colour ic
	      result[t] += tmpSpin * c1res;
	    }  // spin
	  } // spin
	} // spin
      } // spin
    } // t
  
    return result;
  }




  multi1d<DComplex> m1contract( const HeavyMesonBlock& HzU1z0zG5,  const HeavyMesonBlock& HzU2zRzG5,
				const HeavyMesonBlock& HzU3z0zG5,  const HeavyMesonBlock& HzU4zRzG5,
				const SpinMatrix& S1, const SpinMatrix& S2)
  // Contractions for B B potential
  // spin matrix S1 is at 0 and S2 is at R
  // blocks 1,3 are at R and the others at 0

  {
    int length=HzU3z0zG5.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    DComplex tmpSpin, c1res;
  
    result = 0;

    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){
	for (int js=0; js<Nd; js++){
	  for (int ks=0; ks<Nd; ks++){
	    for (int ls=0; ls<Nd; ls++){
	      tmpSpin = peekSpin(S1,is,js) * peekSpin(S2,ks,ls);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){
		for (int jc=0; jc<Nc; jc++){
		  c1res += HzU1z0zG5(t,ic,ic,js,is)*HzU2zRzG5(t,jc,jc,ls,ks);
		  c1res -= HzU3z0zG5(t,ic,jc,js,ks)*HzU4zRzG5(t,jc,ic,ls,is);
		}
	      }	  
	      result[t] += tmpSpin * c1res;
	    }
	  }
	}
      }
    }
    return result;
  }


  multi1d<DComplex> m2contract( const HeavyMesonBlock& HzU1z0zG5,  const HeavyMesonBlock& HzD2zRzG5,
				const SpinMatrix& S1, const SpinMatrix& S2)
  // Contractions for B Bs ( or Bu Bd in I=0) potential
  // spin matrix S1 is at 0 and S2 is at R
  // blocks 1 are at R and block 2 at 0
  {
    int length=HzD2zRzG5.length();
    multi1d<DComplex> result;
    result.resize(length);

    DComplex tmpSpin, c1res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){
	for (int js=0; js<Nd; js++){
	  for (int ks=0; ks<Nd; ks++){
	    for (int ls=0; ls<Nd; ls++){
	      tmpSpin = peekSpin(S1,is,js) * peekSpin(S2,ks,ls);
	      c1res=0;
	      for (int ic=0; ic<Nc; ic++){
		for (int jc=0; jc<Nc; jc++){		
		  c1res += HzD2zRzG5(t,jc,jc,ls,ks)*HzU1z0zG5(t,ic,ic,js,is);
		}
	      }	  
	      result[t] += tmpSpin * c1res;
	    }
	  }
	}
      }
    }
    return result;
  }





  multi1d<DComplex> lambdabcontract( const QllBlock& BzU1zD1z0zCG5,
				     const SpinMatrix& S1)
  // Contractions for lambda_B baryon
  {
    int length=BzU1zD1z0zCG5.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde; // \gamma_0 S_i^\dagger\gamma_0

    S1tilde = Gamma(8) * adj(S1) * Gamma(8);

    DComplex tmpSpin,  c1res;
  
    result = 0;
    for (int t=0; t<length; t++){
      result[t] =0;
      for (int is=0; is<Nd; is++){
	for (int js=0; js<Nd; js++){
	  tmpSpin = peekSpin(S1tilde,is,js);
	  c1res=0;
	  for (int ic=0; ic<Nc; ic++){
	    for (int kc=0; kc<Nc; kc++){
	      for (int lc=0; lc<Nc; lc++){
		if (ic != kc && ic != lc && kc != lc) 
		{
		  c1res += BzU1zD1z0zCG5(t, ic, kc, lc, is, js) * antiSymTensor3d(ic,kc,lc);    
		} // if epsilon
	      } // colour lc
	    } // colour kc
	  } // colour ic
	  result[t] += tmpSpin * c1res;
	} // spin js
      } // spin is
    } //t
    return result;
  }

  multi1d<DComplex> sigmabpluscontract( const QllBlock& BzU1zU1z0zCGi,
					const SpinMatrix& S1)
  // Contractions for lambda_B baryon
  {
    int length=BzU1zU1z0zCGi.length();
    multi1d<DComplex> result;
    result.resize(length);
  
    SpinMatrix S1tilde; // \gamma_0 S_i^\dagger\gamma_0

    S1tilde = Gamma(8) * adj(S1) * Gamma(8);

    DComplex tmpSpin,  c1res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){
	for (int js=0; js<Nd; js++){
	  tmpSpin = peekSpin(S1tilde,is,js);
	  c1res=0;
	  for (int ic=0; ic<Nc; ic++){
	    for (int kc=0; kc<Nc; kc++){
	      for (int lc=0; lc<Nc; lc++){
		if (ic != kc && ic != lc && kc != lc) 
		{
		  c1res += (-BzU1zU1z0zCGi(t,ic,kc,lc,is,js) + BzU1zU1z0zCGi(t,kc,ic,lc,js,is))
		    * antiSymTensor3d(ic,kc,lc);    
		} // if epsilon
	      } // colour lc
	    } // colour kc
	  } // colour ic
	  result[t] += tmpSpin * c1res;
	} // spin js
      } // spin is
    } //t
    return result;
  }

  multi1d<DComplex> bcontract( const HeavyMesonBlock& H1,
			       const SpinMatrix& S1)
  // Contractions for B meson
  {
    int length=H1.length();
    int size= H1.size();
    multi1d<DComplex> result;
    result.resize(length);
  
    DComplex tmpSpin,  c1res;
  
    result = 0;
    for (int t=0; t<length; t++){
      for (int is=0; is<Nd; is++){
	for (int js=0; js<Nd; js++){
	  tmpSpin = peekSpin(S1,is,js);
	  c1res=0;
	  for (int ic=0; ic<Nc; ic++){
	    c1res += H1(t,ic,ic,js,is) ;     
	  }
	  result[t] += tmpSpin * c1res;
	}
      }
    }
    return result;
  }



  multi2d<DComplex> c4J2corr(
    // Spin up blocks
    const QllBlock& BzU1zD1z0zCgplus, const QllBlock& BzU2zU2zRzCgplus,
    const QllBlock& BzU2zU4zRzCgplus, const QllBlock& BzU3zD1z0zCgplus, 
    const QllBlock& BzU4zU2zRzCgplus,
    // Spin zero blocks
    const QllBlock& BzU1zD1z0zCg3, const QllBlock& BzU2zU2zRzCg3,
    const QllBlock& BzU2zU4zRzCg3, const QllBlock& BzU3zD1z0zCg3, 
    const QllBlock& BzU4zU2zRzCg3,
    // Spin down blocks
    const QllBlock& BzU1zD1z0zCgminus, const QllBlock& BzU2zU2zRzCgminus,
    const QllBlock& BzU2zU4zRzCgminus, const QllBlock& BzU3zD1z0zCgminus, 
    const QllBlock& BzU4zU2zRzCgminus)
  {
    // Calculate all 9 spin states of the J=1 \otimes J=1 baryon system

    SpinMatrix g_one = 1.0;
    SpinMatrix Cg5 = g_one * Gamma(5);
    SpinMatrix Cg3 = - g_one * Gamma(14);
    SpinMatrix Cgminus = g_one * Gamma(11) + timesI(g_one * Gamma(8)); 
    SpinMatrix Cgplus  = g_one * Gamma(11) + timesMinusI(g_one * Gamma(8)); 

    int length = BzU1zD1z0zCgplus.length();
    multi2d<DComplex> c4J2result;
    c4J2result.resize(9,length);

    multi1d<DComplex> cJ2m2(length), cJ2m1(length), cJ2m0(length), cJ2mneg1(length), cJ2mneg2(length);
    multi1d<DComplex> cJ1m1(length), cJ1m0(length), cJ1mneg1(length), cJ0m0(length);

    // Below cIJ_KL correspons to spin 
    //  I at 0,0
    //  J at r,0
    //  K at x,0
    //  L at x+r,0


  
    // m=2
    multi1d<DComplex> c11_11(length);
    c11_11 = c4contract(BzU1zD1z0zCgplus, BzU2zU2zRzCgplus,
			BzU2zU4zRzCgplus, BzU3zD1z0zCgplus, 
			BzU4zU2zRzCgplus,
			Cgplus, Cgplus);


    //  m=1
    multi1d<DComplex> c10_10(length), c10_01(length), c01_10(length), c01_01(length);
    c10_10 = c4contract(BzU1zD1z0zCgplus, BzU2zU2zRzCg3,
			BzU2zU4zRzCg3, BzU3zD1z0zCgplus, 
			BzU4zU2zRzCg3,
			Cg3, Cgplus);
    c10_01 = c4contract(BzU1zD1z0zCgplus, BzU2zU2zRzCg3,
			BzU2zU4zRzCg3, BzU3zD1z0zCgplus, 
			BzU4zU2zRzCg3,
			Cgplus, Cg3);
    c01_10 = c4contract(BzU1zD1z0zCg3, BzU2zU2zRzCgplus,
			BzU2zU4zRzCgplus, BzU3zD1z0zCg3, 
			BzU4zU2zRzCgplus,
			Cg3, Cgplus);
    c01_01 = c4contract(BzU1zD1z0zCg3, BzU2zU2zRzCgplus,
			BzU2zU4zRzCgplus, BzU3zD1z0zCg3, 
			BzU4zU2zRzCgplus,
			Cgplus, Cg3);


    // m=0
    multi1d<DComplex> c1m1_1m1(length), c1m1_00(length), c1m1_m11(length);
    multi1d<DComplex> c00_1m1(length), c00_00(length), c00_m11(length);
    multi1d<DComplex> cm11_1m1(length), cm11_00(length), cm11_m11(length);
    c1m1_1m1 = c4contract(BzU1zD1z0zCgplus, BzU2zU2zRzCgminus,
			  BzU2zU4zRzCgminus, BzU3zD1z0zCgplus, 
			  BzU4zU2zRzCgminus,
			  Cgminus, Cgplus);
    c1m1_00 = c4contract(BzU1zD1z0zCgplus, BzU2zU2zRzCgminus,
			 BzU2zU4zRzCgminus, BzU3zD1z0zCgplus, 
			 BzU4zU2zRzCgminus,
			 Cg3, Cg3);
    c1m1_m11 = c4contract(BzU1zD1z0zCgplus, BzU2zU2zRzCgminus,
			  BzU2zU4zRzCgminus, BzU3zD1z0zCgplus, 
			  BzU4zU2zRzCgminus,
			  Cgplus, Cgminus);
    c00_1m1 = c4contract(BzU1zD1z0zCg3, BzU2zU2zRzCg3,
			 BzU2zU4zRzCg3, BzU3zD1z0zCg3, 
			 BzU4zU2zRzCg3,
			 Cgminus, Cgplus);
    c00_00 = c4contract(BzU1zD1z0zCg3, BzU2zU2zRzCg3,
			BzU2zU4zRzCg3, BzU3zD1z0zCg3, 
			BzU4zU2zRzCg3,
			Cg3, Cg3);
    c00_m11 = c4contract(BzU1zD1z0zCg3, BzU2zU2zRzCg3,
			 BzU2zU4zRzCg3, BzU3zD1z0zCg3, 
			 BzU4zU2zRzCg3,
			 Cgplus, Cgminus);
    cm11_1m1 = c4contract(BzU1zD1z0zCgminus, BzU2zU2zRzCgplus,
			  BzU2zU4zRzCgplus, BzU3zD1z0zCgminus, 
			  BzU4zU2zRzCgplus,
			  Cgminus, Cgplus);
    cm11_00 = c4contract(BzU1zD1z0zCgminus, BzU2zU2zRzCgplus,
			 BzU2zU4zRzCgplus, BzU3zD1z0zCgminus, 
			 BzU4zU2zRzCgplus,
			 Cg3, Cg3);
    cm11_m11 = c4contract(BzU1zD1z0zCgminus, BzU2zU2zRzCgplus,
			  BzU2zU4zRzCgplus, BzU3zD1z0zCgminus, 
			  BzU4zU2zRzCgplus,
			  Cgplus, Cgminus);



    //  m=-1
    multi1d<DComplex> cm10_m10(length), cm10_0m1(length), c0m1_m10(length), c0m1_0m1(length);
    cm10_m10 = c4contract(BzU1zD1z0zCgminus, BzU2zU2zRzCg3,
			  BzU2zU4zRzCg3, BzU3zD1z0zCgminus, 
			  BzU4zU2zRzCg3,
			  Cg3, Cgminus);
    cm10_0m1 = c4contract(BzU1zD1z0zCgminus, BzU2zU2zRzCg3,
			  BzU2zU4zRzCg3, BzU3zD1z0zCgminus, 
			  BzU4zU2zRzCg3,
			  Cgminus, Cg3);
    c0m1_m10 = c4contract(BzU1zD1z0zCg3, BzU2zU2zRzCgminus,
			  BzU2zU4zRzCgminus, BzU3zD1z0zCg3, 
			  BzU4zU2zRzCgminus,
			  Cg3, Cgminus);
    c0m1_0m1 = c4contract(BzU1zD1z0zCg3, BzU2zU2zRzCgminus,
			  BzU2zU4zRzCgminus, BzU3zD1z0zCg3, 
			  BzU4zU2zRzCgminus,
			  Cgminus, Cg3);


    // m=-2
    multi1d<DComplex> cm1m1_m1m1(length);
    cm1m1_m1m1 = c4contract(BzU1zD1z0zCgminus, BzU2zU2zRzCgminus,
			    BzU2zU4zRzCgminus, BzU3zD1z0zCgminus, 
			    BzU4zU2zRzCgminus,
			    Cgminus, Cgminus);






    for (int t=0; t < length; t++) {
      cJ2m2[t] = c11_11[t];
      cJ2m1[t] = (c01_01[t] + c01_10[t] + c10_01[t] + c10_10[t]) /2.;
      cJ1m1[t] = (c01_01[t] - c01_10[t] - c10_01[t] + c10_10[t]) /2.;
      cJ2m0[t] = (cm11_m11[t] + 2 * cm11_00[t] + cm11_1m1[t] 
		  + 2 * c00_m11[t] + 4 * c00_00[t] + 2 * c00_1m1[t]
		  + c1m1_m11[t] + 2 * c1m1_00[t] + c1m1_1m1[t])/6.;
      cJ1m0[t] = (cm11_m11[t] - cm11_1m1[t] - c1m1_m11[t] +  c1m1_1m1[t])/2.;
      cJ0m0[t] = (cm11_m11[t] - cm11_00[t] + cm11_1m1[t] 
		  - c00_m11[t] + c00_00[t] - c00_1m1[t]
		  + c1m1_m11[t] - c1m1_00[t] + c1m1_1m1[t])/3.;
      cJ2mneg1[t] = (c0m1_0m1[t] + c0m1_m10[t] + cm10_0m1[t] + cm10_m10[t])/2.;
      cJ1mneg1[t] = (c0m1_0m1[t] - c0m1_m10[t] - cm10_0m1[t] + cm10_m10[t])/2.;
      cJ2mneg2[t] = cm1m1_m1m1[t];

      c4J2result[0][t] = cJ2m2[t];    // J=2, m=2 state
      c4J2result[1][t] = cJ2m1[t];    // J=2, m=1 state
      c4J2result[2][t] = cJ2m0[t];    // J=2, m=0 state
      c4J2result[3][t] = cJ2mneg1[t]; // J=2, m=-1 state
      c4J2result[4][t] = cJ2mneg2[t]; // J=2, m=-2 state
      c4J2result[5][t] = cJ1m1[t];;   // J=1, m=1 state
      c4J2result[6][t] = cJ1m0[t]; // J=1, m=0 state
      c4J2result[7][t] = cJ1mneg1[t]; // J=1, m=-1 state
      c4J2result[8][t] = cJ0m0[t];    // J=0, m=0 state
    }
  
    return c4J2result;
  }






  multi2d<DComplex> c5J2corr(
    // Spin up blocks
    const QllBlock& BzU1zU1z0zCgplus, const QllBlock& BzU1zU3z0zCgplus,
    const QllBlock& BzU2zU2zRzCgplus, const QllBlock& BzU2zU4zRzCgplus,
    const QllBlock& BzU3zU1z0zCgplus, const QllBlock& BzU3zU3z0zCgplus,
    const QllBlock& BzU4zU2zRzCgplus, const QllBlock& BzU4zU4zRzCgplus,
    // Spin zero blocks
    const QllBlock& BzU1zU1z0zCg3, const QllBlock& BzU1zU3z0zCg3,
    const QllBlock& BzU2zU2zRzCg3, const QllBlock& BzU2zU4zRzCg3,
    const QllBlock& BzU3zU1z0zCg3, const QllBlock& BzU3zU3z0zCg3,
    const QllBlock& BzU4zU2zRzCg3, const QllBlock& BzU4zU4zRzCg3,
    // Spin down blocks
    const QllBlock& BzU1zU1z0zCgminus, const QllBlock& BzU1zU3z0zCgminus,
    const QllBlock& BzU2zU2zRzCgminus, const QllBlock& BzU2zU4zRzCgminus,
    const QllBlock& BzU3zU1z0zCgminus, const QllBlock& BzU3zU3z0zCgminus,
    const QllBlock& BzU4zU2zRzCgminus, const QllBlock& BzU4zU4zRzCgminus)
  {
    // Calculate all 9 spin states of the J=1 \otimes J=1 baryon system

    SpinMatrix g_one = 1.0;
    SpinMatrix Cg5 = g_one * Gamma(5);
    SpinMatrix Cg3 = - g_one * Gamma(14);
    SpinMatrix Cgminus = g_one * Gamma(11) + timesI(g_one * Gamma(8)); 
    SpinMatrix Cgplus  = g_one * Gamma(11) + timesMinusI(g_one * Gamma(8)); 

    int length = BzU1zU1z0zCgplus.length();
    multi2d<DComplex> c5J2result;
    c5J2result.resize(9,length);

    multi1d<DComplex> cJ2m2(length), cJ2m1(length), cJ2m0(length), cJ2mneg1(length), cJ2mneg2(length);
    multi1d<DComplex> cJ1m1(length), cJ1m0(length), cJ1mneg1(length), cJ0m0(length);

    // Below cIJ_KL correspons to spin 
    //  I at 0,0
    //  J at r,0
    //  K at x,0
    //  L at x+r,0
  
    // m=2
    multi1d<DComplex> c11_11(length);
    c11_11 = c5contract( BzU1zU1z0zCgplus,  BzU1zU3z0zCgplus,
			 BzU2zU2zRzCgplus,  BzU2zU4zRzCgplus,
			 BzU3zU1z0zCgplus,  BzU3zU3z0zCgplus,
			 BzU4zU2zRzCgplus,  BzU4zU4zRzCgplus,
			 Cgplus, Cgplus);


    //  m=1
    multi1d<DComplex> c10_10(length), c10_01(length), c01_10(length), c01_01(length);
    c10_10 = c5contract( BzU1zU1z0zCgplus,  BzU1zU3z0zCgplus,
			 BzU2zU2zRzCg3,  BzU2zU4zRzCg3,
			 BzU3zU1z0zCgplus,  BzU3zU3z0zCgplus,
			 BzU4zU2zRzCg3,  BzU4zU4zRzCg3,
			 Cg3, Cgplus);
    c10_01 = c5contract( BzU1zU1z0zCgplus,  BzU1zU3z0zCgplus,
			 BzU2zU2zRzCg3,  BzU2zU4zRzCg3,
			 BzU3zU1z0zCgplus,  BzU3zU3z0zCgplus,
			 BzU4zU2zRzCg3,  BzU4zU4zRzCg3,
			 Cgplus, Cg3);
    c01_10 = c5contract( BzU1zU1z0zCg3,  BzU1zU3z0zCg3,
			 BzU2zU2zRzCgplus,  BzU2zU4zRzCgplus,
			 BzU3zU1z0zCg3,  BzU3zU3z0zCg3,
			 BzU4zU2zRzCgplus,  BzU4zU4zRzCgplus,
			 Cg3, Cgplus);
    c01_01 = c5contract(  BzU1zU1z0zCg3,  BzU1zU3z0zCg3,
			  BzU2zU2zRzCgplus,  BzU2zU4zRzCgplus,
			  BzU3zU1z0zCg3,  BzU3zU3z0zCg3,
			  BzU4zU2zRzCgplus,  BzU4zU4zRzCgplus,
			  Cgplus, Cg3);


    // m=0
    multi1d<DComplex> c1m1_1m1(length), c1m1_00(length), c1m1_m11(length);
    multi1d<DComplex> c00_1m1(length), c00_00(length), c00_m11(length);
    multi1d<DComplex> cm11_1m1(length), cm11_00(length), cm11_m11(length);
    c1m1_1m1 = c5contract( BzU1zU1z0zCgplus,  BzU1zU3z0zCgplus,
			   BzU2zU2zRzCgminus,  BzU2zU4zRzCgminus,
			   BzU3zU1z0zCgplus,  BzU3zU3z0zCgplus,
			   BzU4zU2zRzCgminus,  BzU4zU4zRzCgminus,
			   Cgminus, Cgplus);
    c1m1_00 = c5contract( BzU1zU1z0zCgplus,  BzU1zU3z0zCgplus,
			  BzU2zU2zRzCgminus,  BzU2zU4zRzCgminus,
			  BzU3zU1z0zCgplus,  BzU3zU3z0zCgplus,
			  BzU4zU2zRzCgminus,  BzU4zU4zRzCgminus,
			  Cg3, Cg3);
    c1m1_m11 = c5contract( BzU1zU1z0zCgplus,  BzU1zU3z0zCgplus,
			   BzU2zU2zRzCgminus,  BzU2zU4zRzCgminus,
			   BzU3zU1z0zCgplus,  BzU3zU3z0zCgplus,
			   BzU4zU2zRzCgminus,  BzU4zU4zRzCgminus,
			   Cgplus, Cgminus);
    c00_1m1 = c5contract( BzU1zU1z0zCg3,  BzU1zU3z0zCg3,
			  BzU2zU2zRzCg3,  BzU2zU4zRzCg3,
			  BzU3zU1z0zCg3,  BzU3zU3z0zCg3,
			  BzU4zU2zRzCg3,  BzU4zU4zRzCg3,
			  Cgminus, Cgplus);
    c00_00 = c5contract( BzU1zU1z0zCg3,  BzU1zU3z0zCg3,
			 BzU2zU2zRzCg3,  BzU2zU4zRzCg3,
			 BzU3zU1z0zCg3,  BzU3zU3z0zCg3,
			 BzU4zU2zRzCg3,  BzU4zU4zRzCg3,
			 Cg3, Cg3);
    c00_m11 = c5contract( BzU1zU1z0zCg3,  BzU1zU3z0zCg3,
			  BzU2zU2zRzCg3,  BzU2zU4zRzCg3,
			  BzU3zU1z0zCg3,  BzU3zU3z0zCg3,
			  BzU4zU2zRzCg3,  BzU4zU4zRzCg3,
			  Cgplus, Cgminus);
    cm11_1m1 = c5contract( BzU1zU1z0zCgminus,  BzU1zU3z0zCgminus,
			   BzU2zU2zRzCgplus,  BzU2zU4zRzCgplus,
			   BzU3zU1z0zCgminus,  BzU3zU3z0zCgminus,
			   BzU4zU2zRzCgplus,  BzU4zU4zRzCgplus,
			   Cgminus, Cgplus);
    cm11_00 = c5contract( BzU1zU1z0zCgminus,  BzU1zU3z0zCgminus,
			  BzU2zU2zRzCgplus,  BzU2zU4zRzCgplus,
			  BzU3zU1z0zCgminus,  BzU3zU3z0zCgminus,
			  BzU4zU2zRzCgplus,  BzU4zU4zRzCgplus,
			  Cg3, Cg3);
    cm11_m11 = c5contract( BzU1zU1z0zCgminus,  BzU1zU3z0zCgminus,
			   BzU2zU2zRzCgplus,  BzU2zU4zRzCgplus,
			   BzU3zU1z0zCgminus,  BzU3zU3z0zCgminus,
			   BzU4zU2zRzCgplus,  BzU4zU4zRzCgplus,
			   Cgplus, Cgminus);



    //  m=-1
    multi1d<DComplex> cm10_m10(length), cm10_0m1(length), c0m1_m10(length), c0m1_0m1(length);
    cm10_m10 = c5contract( BzU1zU1z0zCgminus,  BzU1zU3z0zCgminus,
			   BzU2zU2zRzCg3,  BzU2zU4zRzCg3,
			   BzU3zU1z0zCgminus,  BzU3zU3z0zCgminus,
			   BzU4zU2zRzCg3,  BzU4zU4zRzCg3,
			   Cg3, Cgminus);
    cm10_0m1 = c5contract( BzU1zU1z0zCgminus,  BzU1zU3z0zCgminus,
			   BzU2zU2zRzCg3,  BzU2zU4zRzCg3,
			   BzU3zU1z0zCgminus,  BzU3zU3z0zCgminus,
			   BzU4zU2zRzCg3,  BzU4zU4zRzCg3,
			   Cgminus, Cg3);
    c0m1_m10 = c5contract( BzU1zU1z0zCg3,  BzU1zU3z0zCg3,
			   BzU2zU2zRzCgminus,  BzU2zU4zRzCgminus,
			   BzU3zU1z0zCg3,  BzU3zU3z0zCg3,
			   BzU4zU2zRzCgminus,  BzU4zU4zRzCgminus,
			   Cg3, Cgminus);
    c0m1_0m1 = c5contract(  BzU1zU1z0zCg3,  BzU1zU3z0zCg3,
			    BzU2zU2zRzCgminus,  BzU2zU4zRzCgminus,
			    BzU3zU1z0zCg3,  BzU3zU3z0zCg3,
			    BzU4zU2zRzCgminus,  BzU4zU4zRzCgminus,
			    Cgminus, Cg3);


    // m=-2
    multi1d<DComplex> cm1m1_m1m1(length);
    cm1m1_m1m1 = c5contract( BzU1zU1z0zCgminus,  BzU1zU3z0zCgminus,
			     BzU2zU2zRzCgminus,  BzU2zU4zRzCgminus,
			     BzU3zU1z0zCgminus,  BzU3zU3z0zCgminus,
			     BzU4zU2zRzCgminus,  BzU4zU4zRzCgminus,
			     Cgminus, Cgminus);






    for (int t=0; t < length; t++) {
      cJ2m2[t] = c11_11[t];
      cJ2m1[t] = (c01_01[t] + c01_10[t] + c10_01[t] + c10_10[t]) /2.;
      cJ1m1[t] = (c01_01[t] - c01_10[t] - c10_01[t] + c10_10[t]) /2.;
      cJ2m0[t] = (cm11_m11[t] + 2 * cm11_00[t] + cm11_1m1[t] 
		  + 2 * c00_m11[t] + 4 * c00_00[t] + 2 * c00_1m1[t]
		  + c1m1_m11[t] + 2 * c1m1_00[t] + c1m1_1m1[t])/6.;
      cJ1m0[t] = (cm11_m11[t] - cm11_1m1[t] - c1m1_m11[t] +  c1m1_1m1[t])/2.;
      cJ0m0[t] = (cm11_m11[t] - cm11_00[t] + cm11_1m1[t] 
		  - c00_m11[t] + c00_00[t] - c00_1m1[t]
		  + c1m1_m11[t] - c1m1_00[t] + c1m1_1m1[t])/3.;
      cJ2mneg1[t] = (c0m1_0m1[t] + c0m1_m10[t] + cm10_0m1[t] + cm10_m10[t])/2.;
      cJ1mneg1[t] = (c0m1_0m1[t] - c0m1_m10[t] - cm10_0m1[t] + cm10_m10[t])/2.;
      cJ2mneg2[t] = cm1m1_m1m1[t];

      c5J2result[0][t] = cJ2m2[t];    // J=2, m=2 state
      c5J2result[1][t] = cJ2m1[t];    // J=2, m=1 state
      c5J2result[2][t] = cJ2m0[t];    // J=2, m=0 state
      c5J2result[3][t] = cJ2mneg1[t]; // J=2, m=-1 state
      c5J2result[4][t] = cJ2mneg2[t]; // J=2, m=-2 state
      c5J2result[5][t] = cJ1m1[t];;   // J=1, m=1 state
      c5J2result[6][t] = cJ1m0[t]; // J=1, m=0 state
      c5J2result[7][t] = cJ1mneg1[t]; // J=1, m=-1 state
      c5J2result[8][t] = cJ0m0[t];    // J=0, m=0 state
    }
  
    return c5J2result;
  }




  multi2d<DComplex> c6J2corr(
    // Spin up blocks
    const QllBlock& BzD1zD1z0zCgplus, const QllBlock& BzU2zU2zRzCgplus,
    // Spin zero blocks
    const QllBlock& BzD1zD1z0zCg3, const QllBlock& BzU2zU2zRzCg3,
    // Spin down blocks
    const QllBlock& BzD1zD1z0zCgminus, const QllBlock& BzU2zU2zRzCgminus)
  {
    // Calculate all 9 spin states of the J=1 \otimes J=1 baryon system

    SpinMatrix g_one = 1.0;
    SpinMatrix Cg5 = g_one * Gamma(5);
    SpinMatrix Cg3 = - g_one * Gamma(14);
    SpinMatrix Cgminus = g_one * Gamma(11) + timesI(g_one * Gamma(8)); 
    SpinMatrix Cgplus  = g_one * Gamma(11) + timesMinusI(g_one * Gamma(8)); 

    int length = BzD1zD1z0zCgplus.length();
    multi2d<DComplex> c6J2result;
    c6J2result.resize(9,length);

    multi1d<DComplex> cJ2m2(length), cJ2m1(length), cJ2m0(length), cJ2mneg1(length), cJ2mneg2(length);
    multi1d<DComplex> cJ1m1(length), cJ1m0(length), cJ1mneg1(length), cJ0m0(length);

    // Below cIJ_KL correspons to spin 
    //  I at 0,0
    //  J at r,0
    //  K at x,0
    //  L at x+r,0
  
    // m=2
    multi1d<DComplex> c11_11(length);
    c11_11 = c6contract( BzD1zD1z0zCgplus,  BzU2zU2zRzCgplus,
			 Cgplus, Cgplus);

    //  m=1
    multi1d<DComplex> c10_10(length), c10_01(length), c01_10(length), c01_01(length);
    c10_10 = c6contract( BzD1zD1z0zCgplus,  BzU2zU2zRzCg3,
			 Cg3, Cgplus);
    c10_01 = c6contract( BzD1zD1z0zCgplus,  BzU2zU2zRzCg3,
			 Cgplus, Cg3);
    c01_10 = c6contract( BzD1zD1z0zCg3,  BzU2zU2zRzCgplus,
			 Cg3, Cgplus);
    c01_01 = c6contract(  BzD1zD1z0zCg3,  BzU2zU2zRzCgplus,
			  Cgplus, Cg3);

    // m=0
    multi1d<DComplex> c1m1_1m1(length), c1m1_00(length), c1m1_m11(length);
    multi1d<DComplex> c00_1m1(length), c00_00(length), c00_m11(length);
    multi1d<DComplex> cm11_1m1(length), cm11_00(length), cm11_m11(length);
    c1m1_1m1 = c6contract( BzD1zD1z0zCgplus,  BzU2zU2zRzCgminus,
			   Cgminus, Cgplus);
    c1m1_00 = c6contract( BzD1zD1z0zCgplus,  BzU2zU2zRzCgminus,
			  Cg3, Cg3);
    c1m1_m11 = c6contract( BzD1zD1z0zCgplus,  BzU2zU2zRzCgminus,
			   Cgplus, Cgminus);
    c00_1m1 = c6contract( BzD1zD1z0zCg3,  BzU2zU2zRzCg3,
			  Cgminus, Cgplus);
    c00_00 = c6contract( BzD1zD1z0zCg3,  BzU2zU2zRzCg3,
			 Cg3, Cg3);
    c00_m11 = c6contract( BzD1zD1z0zCg3,  BzU2zU2zRzCg3,
			  Cgplus, Cgminus);
    cm11_1m1 = c6contract( BzD1zD1z0zCgminus,  BzU2zU2zRzCgplus,
			   Cgminus, Cgplus);
    cm11_00 = c6contract( BzD1zD1z0zCgminus,  BzU2zU2zRzCgplus,
			  Cg3, Cg3);
    cm11_m11 = c6contract( BzD1zD1z0zCgminus,  BzU2zU2zRzCgplus,
			   Cgplus, Cgminus);

    //  m=-1
    multi1d<DComplex> cm10_m10(length), cm10_0m1(length), c0m1_m10(length), c0m1_0m1(length);
    cm10_m10 = c6contract( BzD1zD1z0zCgminus,  BzU2zU2zRzCg3,
			   Cg3, Cgminus);
    cm10_0m1 = c6contract( BzD1zD1z0zCgminus,  BzU2zU2zRzCg3,
			   Cgminus, Cg3);
    c0m1_m10 = c6contract( BzD1zD1z0zCg3,  BzU2zU2zRzCgminus,
			   Cg3, Cgminus);
    c0m1_0m1 = c6contract(  BzD1zD1z0zCg3,  BzU2zU2zRzCgminus,
			    Cgminus, Cg3);

    // m=-2
    multi1d<DComplex> cm1m1_m1m1(length);
    cm1m1_m1m1 = c6contract( BzD1zD1z0zCgminus,  BzU2zU2zRzCgminus,
			     Cgminus, Cgminus);






    for (int t=0; t < length; t++) {
      cJ2m2[t] = c11_11[t];
      cJ2m1[t] = (c01_01[t] + c01_10[t] + c10_01[t] + c10_10[t]) /2.;
      cJ1m1[t] = (c01_01[t] - c01_10[t] - c10_01[t] + c10_10[t]) /2.;
      cJ2m0[t] = (cm11_m11[t] + 2 * cm11_00[t] + cm11_1m1[t] 
		  + 2 * c00_m11[t] + 4 * c00_00[t] + 2 * c00_1m1[t]
		  + c1m1_m11[t] + 2 * c1m1_00[t] + c1m1_1m1[t])/6.;
      cJ1m0[t] = (cm11_m11[t] - cm11_1m1[t] - c1m1_m11[t] +  c1m1_1m1[t])/2.;
      cJ0m0[t] = (cm11_m11[t] - cm11_00[t] + cm11_1m1[t] 
		  - c00_m11[t] + c00_00[t] - c00_1m1[t]
		  + c1m1_m11[t] - c1m1_00[t] + c1m1_1m1[t])/3.;
      cJ2mneg1[t] = (c0m1_0m1[t] + c0m1_m10[t] + cm10_0m1[t] + cm10_m10[t])/2.;
      cJ1mneg1[t] = (c0m1_0m1[t] - c0m1_m10[t] - cm10_0m1[t] + cm10_m10[t])/2.;
      cJ2mneg2[t] = cm1m1_m1m1[t];

      c6J2result[0][t] = cJ2m2[t];    // J=2, m=2 state
      c6J2result[1][t] = cJ2m1[t];    // J=2, m=1 state
      c6J2result[2][t] = cJ2m0[t];    // J=2, m=0 state
      c6J2result[3][t] = cJ2mneg1[t]; // J=2, m=-1 state
      c6J2result[4][t] = cJ2mneg2[t]; // J=2, m=-2 state
      c6J2result[5][t] = cJ1m1[t];;   // J=1, m=1 state
      c6J2result[6][t] = cJ1m0[t]; // J=1, m=0 state
      c6J2result[7][t] = cJ1mneg1[t]; // J=1, m=-1 state
      c6J2result[8][t] = cJ0m0[t];    // J=0, m=0 state
    }
  
    return c6J2result;
  }


  multi2d<DComplex> d2J32corr(
    // Spin up blocks
    const QllBlock& BzU1zU1z0zCgplus, const QllBlock& BzU3zU1z0zCgplus, const QllBlock& BzU1zU3z0zCgplus, 
    const HeavyMesonBlock& HzU2zRzGup, const HeavyMesonBlock& HzU4zRzGup,
    // Spin zero blocks
    const QllBlock& BzU1zU1z0zCg3, const QllBlock& BzU3zU1z0zCg3, const QllBlock& BzU1zU3z0zCg3, 
    // Spin down blocks
    const QllBlock& BzU1zU1z0zCgminus, const QllBlock& BzU3zU1z0zCgminus, const QllBlock& BzU1zU3z0zCgminus, 
    const HeavyMesonBlock& HzU2zRzGdown, const HeavyMesonBlock& HzU4zRzGdown)
  {
    // Calculate all 6 spin states of the J=1 \otimes J=1/2 baryon - meson system

    SpinMatrix g_one = 1.0;
    SpinMatrix Cg5 = g_one * Gamma(5);
    SpinMatrix Cg3 = - g_one * Gamma(14);
    SpinMatrix Cgminus = g_one * Gamma(11) + timesI(g_one * Gamma(8)); 
    SpinMatrix Cgplus  = g_one * Gamma(11) + timesMinusI(g_one * Gamma(8)); 
    SpinMatrix Gup = g_one * Gamma(15) + timesMinusI((g_one*Gamma(3))*Gamma(15));
    SpinMatrix Gdown = g_one * Gamma(15) + timesI((g_one*Gamma(3))*Gamma(15));


    int length = BzU1zU1z0zCgplus.length();
    multi2d<DComplex> d2J32result;
    d2J32result.resize(6,length);

    multi1d<DComplex> cJ32m32(length), cJ32m12(length), cJ32mneg12(length), cJ32mneg32(length);
    multi1d<DComplex> cJ12m12(length), cJ12mneg12(length);

    // Below cIJ_KL correspons to spin 
    //  I at 0,0
    //  J at r,0
    //  K at x,0
    //  L at x+r,0
  
    // notation is spin of spin 1 first then of spin1/2 and I don't keep the /2
  

    // m=3/2
    multi1d<DComplex> c11_11(length);
    c11_11 = d2contract(BzU1zU1z0zCgplus, BzU3zU1z0zCgplus, BzU1zU3z0zCgplus, 
			HzU2zRzGup, HzU4zRzGup, 
			Gup, Cgplus);

    //  m=1/2
    multi1d<DComplex> c01_01(length), c1m1_1m1(length), c1m1_01(length), c01_1m1(length);
    c01_01 = d2contract(BzU1zU1z0zCg3, BzU3zU1z0zCg3, BzU1zU3z0zCg3, 
			HzU2zRzGup, HzU4zRzGup, 
			Gup, Cg3);
    c01_1m1 = d2contract(BzU1zU1z0zCg3, BzU3zU1z0zCg3, BzU1zU3z0zCg3, 
			 HzU2zRzGup, HzU4zRzGup, 
			 Gdown, Cgplus);
    c1m1_01 = d2contract(BzU1zU1z0zCgplus, BzU3zU1z0zCgplus, BzU1zU3z0zCgplus, 
			 HzU2zRzGdown, HzU4zRzGdown, 
			 Gup, Cg3);
    c1m1_1m1 = d2contract(BzU1zU1z0zCgplus, BzU3zU1z0zCgplus, BzU1zU3z0zCgplus, 
			  HzU2zRzGdown, HzU4zRzGdown, 
			  Gdown, Cgplus);

    // m=-1/2
    multi1d<DComplex> c0m1_0m1(length), cm11_m11(length), cm11_0m1(length), c0m1_m11(length);
    c0m1_0m1 = d2contract(BzU1zU1z0zCg3, BzU3zU1z0zCg3, BzU1zU3z0zCg3, 
			  HzU2zRzGdown, HzU4zRzGdown, 
			  Gdown, Cg3);
    c0m1_m11 = d2contract(BzU1zU1z0zCg3, BzU3zU1z0zCg3, BzU1zU3z0zCg3, 
			  HzU2zRzGdown, HzU4zRzGdown, 
			  Gup, Cgminus);
    cm11_0m1 = d2contract(BzU1zU1z0zCgminus, BzU3zU1z0zCgminus, BzU1zU3z0zCgminus, 
			  HzU2zRzGup, HzU4zRzGup,
			  Gdown, Cg3);
    cm11_m11 = d2contract(BzU1zU1z0zCgminus, BzU3zU1z0zCgminus, BzU1zU3z0zCgminus, 
			  HzU2zRzGup, HzU4zRzGup,
			  Gup, Cgminus);

    // m=-3/2
    multi1d<DComplex> cm1m1_m1m1(length);
    cm1m1_m1m1 = d2contract(BzU1zU1z0zCgminus, BzU3zU1z0zCgminus, BzU1zU3z0zCgminus, 
			    HzU2zRzGdown, HzU4zRzGdown,
			    Gdown, Cgminus);


    for (int t=0; t < length; t++) {
      // notation is spin of spin 1 first then of spin1/2 and I don't keep the /2
      cJ32m32[t] = c11_11[t];
      cJ32m12[t] = (c1m1_1m1[t] + sqrt(2.) * c1m1_01[t] + sqrt(2.) * c01_1m1[t] + 2.* c01_01[t]) /3.;
      cJ32mneg12[t] = (2.* c0m1_0m1[t] + sqrt(2.) * c0m1_m11[t] + sqrt(2.) *  cm11_0m1[t] + cm11_m11[t]) /3.;
      cJ32mneg32[t] = cm1m1_m1m1[t];
      cJ12m12[t] = (2. * c1m1_1m1[t] - sqrt(2.) * c1m1_01[t] - sqrt(2.) * c01_1m1[t] +  c01_01[t]) /3.;
      cJ12mneg12[t] =  (c0m1_0m1[t] - sqrt(2.) * c0m1_m11[t] - sqrt(2.) *  cm11_0m1[t] + 2. * cm11_m11[t]) /3.;

      d2J32result[0][t] = cJ32m32[t];    // J=3/2, m=3/2 state
      d2J32result[1][t] = cJ32m12[t];    // J=3/2, m=1/2 state
      d2J32result[2][t] = cJ32mneg12[t];    // J=3/2, m=-1/2 state
      d2J32result[3][t] = cJ32mneg32[t]; // J=3/2, m=-3/2 state
      d2J32result[4][t] = cJ12m12[t]; // J=1/2, m=1/2 state
      d2J32result[5][t] = cJ12mneg12[t];;   // J=1/2, m=-1/2 state
    }
  
    return d2J32result;
  }



  multi2d<DComplex> d3J32corr(
    // Spin up blocks
    const QllBlock& BzU1zU1z0zCgplus, 
    const HeavyMesonBlock& HzU2zRzGup, 
    // Spin zero blocks
    const QllBlock& BzU1zU1z0zCg3,
    // Spin down blocks
    const QllBlock& BzU1zU1z0zCgminus,
    const HeavyMesonBlock& HzU2zRzGdown)
  {
    // Calculate all 6 spin states of the J=1 \otimes J=1/2 baryon - meson system

    SpinMatrix g_one = 1.0;
    SpinMatrix Cg5 = g_one * Gamma(5);
    SpinMatrix Cg3 = - g_one * Gamma(14);
    SpinMatrix Cgminus = g_one * Gamma(11) + timesI(g_one * Gamma(8)); 
    SpinMatrix Cgplus  = g_one * Gamma(11) + timesMinusI(g_one * Gamma(8)); 
    SpinMatrix Gup = g_one * Gamma(15) + timesMinusI((g_one*Gamma(3))*Gamma(15));
    SpinMatrix Gdown = g_one * Gamma(15) + timesI((g_one*Gamma(3))*Gamma(15));


    int length = BzU1zU1z0zCgplus.length();
    multi2d<DComplex> d3J32result;
    d3J32result.resize(6,length);

    multi1d<DComplex> cJ32m32(length), cJ32m12(length), cJ32mneg12(length), cJ32mneg32(length);
    multi1d<DComplex> cJ12m12(length), cJ12mneg12(length);

    // Below cIJ_KL correspons to spin 
    //  I at 0,0
    //  J at r,0
    //  K at x,0
    //  L at x+r,0
  
    // notation is spin of spin 1 first then of spin1/2 and I don't keep the /2
  

    // m=3/2
    multi1d<DComplex> c11_11(length);
    c11_11 = d3contract(BzU1zU1z0zCgplus, 
			HzU2zRzGup, 
			Gup, Cgplus);

    //  m=1/2
    multi1d<DComplex> c01_01(length), c1m1_1m1(length), c1m1_01(length), c01_1m1(length);
    c01_01 = d3contract(BzU1zU1z0zCg3,
			HzU2zRzGup, 
			Gup, Cg3);
    c01_1m1 = d3contract(BzU1zU1z0zCg3, 
			 HzU2zRzGup, 
			 Gdown, Cgplus);
    c1m1_01 = d3contract(BzU1zU1z0zCgplus,
			 HzU2zRzGdown,
			 Gup, Cg3);
    c1m1_1m1 = d3contract(BzU1zU1z0zCgplus,
			  HzU2zRzGdown,
			  Gdown, Cgplus);

    // m=-1/2
    multi1d<DComplex> c0m1_0m1(length), cm11_m11(length), cm11_0m1(length), c0m1_m11(length);
    c0m1_0m1 = d3contract(BzU1zU1z0zCg3,
			  HzU2zRzGdown,
			  Gdown, Cg3);
    c0m1_m11 = d3contract(BzU1zU1z0zCg3,
			  HzU2zRzGdown,
			  Gup, Cgminus);
    cm11_0m1 = d3contract(BzU1zU1z0zCgminus,
			  HzU2zRzGup,
			  Gdown, Cg3);
    cm11_m11 = d3contract(BzU1zU1z0zCgminus,
			  HzU2zRzGup,
			  Gup, Cgminus);

    // m=-3/2
    multi1d<DComplex> cm1m1_m1m1(length);
    cm1m1_m1m1 = d3contract(BzU1zU1z0zCgminus,
			    HzU2zRzGdown,
			    Gdown, Cgminus);


    for (int t=0; t < length; t++) {
      // notation is spin of spin 1 first then of spin1/2 and I don't keep the /2
      cJ32m32[t] = c11_11[t];
      cJ32m12[t] = (c1m1_1m1[t] + sqrt(2.) * c1m1_01[t] + sqrt(2.) * c01_1m1[t] + 2.* c01_01[t]) /3.;
      cJ32mneg12[t] = (2.* c0m1_0m1[t] + sqrt(2.) * c0m1_m11[t] + sqrt(2.) *  cm11_0m1[t] + cm11_m11[t]) /3.;
      cJ32mneg32[t] = cm1m1_m1m1[t];
      cJ12m12[t] = (2. * c1m1_1m1[t] - sqrt(2.) * c1m1_01[t] - sqrt(2.) * c01_1m1[t] +  c01_01[t]) /3.;
      cJ12mneg12[t] =  (c0m1_0m1[t] - sqrt(2.) * c0m1_m11[t] - sqrt(2.) *  cm11_0m1[t] + 2. * cm11_m11[t]) /3.;

      d3J32result[0][t] = cJ32m32[t];    // J=3/2, m=3/2 state
      d3J32result[1][t] = cJ32m12[t];    // J=3/2, m=1/2 state
      d3J32result[2][t] = cJ32mneg12[t];    // J=3/2, m=-1/2 state
      d3J32result[3][t] = cJ32mneg32[t]; // J=3/2, m=-3/2 state
      d3J32result[4][t] = cJ12m12[t]; // J=1/2, m=1/2 state
      d3J32result[5][t] = cJ12mneg12[t];;   // J=1/2, m=-1/2 state
    }
  
    return d3J32result;
  }


}
