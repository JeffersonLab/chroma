// $Id: wloop.cc,v 1.1 2008-06-26 14:58:35 mcneile Exp $
/*! \file
 *  \brief Compute simple Wilson loops for use in
 *         alpha_s calculation. These are the Wilson loops used by
 *         the HPQCD collaboration.
 *
 * Accurate determinations of alpha(s) from realistic lattice QCD.
 * By HPQCD Collaboration and UKQCD Collaboration 
 * Published in Phys.Rev.Lett.95:052002,2005.
 * e-Print: hep-lat/0503005
 *
 */

#include "chromabase.h"
#include "meas/glue/wloop.h"
#include "meas/glue/polylp.h"

namespace Chroma 
{

  // Primitive way for now to indicate the time direction
  static int tDir() {return Nd-1;}

  //! Return the value of the 11, 12, 13, 14 Wilson loops
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param plane_plaq  plane plaquette average (Write)
   * \param link        space-time average link (Write)
   */

  void Wloop(const multi1d<LatticeColorMatrix>& u, 
	     multi2d<Double>& plane_plaq_11,
	     multi2d<Double>& plane_plaq_12, 
	     multi2d<Double>& plane_plaq_13, 
	     multi2d<Double>& plane_plaq_14 )
  {
    START_CODE();

    plane_plaq_11.resize(Nd,Nd);
    plane_plaq_12.resize(Nd,Nd);
    plane_plaq_13.resize(Nd,Nd);
    plane_plaq_14.resize(Nd,Nd);

    // Compute the average Wilson loops
    for(int mu=1; mu < Nd; ++mu)
    {
      for(int nu=0; nu < mu; ++nu)
      {

	/** compute the 1x1 Wilson loop **/
	/* staple_0 = u(x,mu) * u(x+mu,nu)*u_dag(x+nu,mu) */
	LatticeColorMatrix tmp0 = adj(shift(u[mu],FORWARD,nu)) ;
	LatticeColorMatrix staple_0 = u[mu]*shift(u[nu],FORWARD,mu) * tmp0 ;

	/*      = sum(tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu))) */
	plane_plaq_11[mu][nu] = sum(real(trace(staple_0 * adj(u[nu]))));

	/** compute the 1x2 Wilson loop **/
	LatticeColorMatrix staple_1 = u[mu]*shift(staple_0 ,FORWARD,mu) * tmp0 ;
	plane_plaq_12[mu][nu] = sum(real(trace(staple_1 * adj(u[nu]))));

	/** compute the 1x3 Wilson loop **/
	LatticeColorMatrix staple_2 = u[mu]*shift(staple_1 ,FORWARD,mu) * tmp0 ;
	plane_plaq_13[mu][nu] = sum(real(trace(staple_2 * adj(u[nu]))));

	/** compute the 1x4 Wilson loop **/
	LatticeColorMatrix staple_3 = u[mu]*shift(staple_2 ,FORWARD,mu) * tmp0 ;
	plane_plaq_14[mu][nu] = sum(real(trace(staple_3 * adj(u[nu]))));

      }
    }

    // Normalize the planes
    for(int mu=1; mu < Nd; ++mu)
      for(int nu=0; nu < mu; ++nu)
      {
	plane_plaq_11[mu][nu] /= Double(Layout::vol()*Nc);
	plane_plaq_11[nu][mu] = plane_plaq_11[mu][nu];

	plane_plaq_12[mu][nu] /= Double(Layout::vol()*Nc);
	plane_plaq_12[nu][mu] = plane_plaq_12[mu][nu];

	plane_plaq_13[mu][nu] /= Double(Layout::vol()*Nc);
	plane_plaq_13[nu][mu] = plane_plaq_13[mu][nu];

	plane_plaq_14[mu][nu] /= Double(Layout::vol()*Nc);
	plane_plaq_14[nu][mu] = plane_plaq_14[mu][nu];

      }

    END_CODE();
}



  //! Return the value of the bent 12 Wilson loops
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param mu          direction for plane 0 (Read)
   * \param nu          direction for plane 0 (Read)
   * \param eta         direction for plane 1 (Read)
   * \param sigma       direction for plane 1 (Read)
   * \param ans         bent rectangle  (Write)
   */

  void Wloop_bent(const multi1d<LatticeColorMatrix>& u, 
		  int mu, int nu, int eta,
		  Double & ans)
  {
    START_CODE();

    ans = zero ; 

    LatticeColorMatrix staple_0 = u[mu]*shift(u[nu],FORWARD,mu) * 
      adj(shift(u[mu],FORWARD,nu)) ;

    LatticeColorMatrix staple_1 = u[eta]*shift(u[nu],FORWARD,eta) * 
      adj(shift(u[eta],FORWARD,nu)) ;

    ans += sum(real(trace(staple_0 * adj(staple_1) )));


    LatticeColorMatrix tmp0 = shift(u[mu],BACKWARD,mu) ;
    LatticeColorMatrix staple_2 = adj(tmp0)*shift(u[nu],BACKWARD,mu) * 
      (shift(tmp0,FORWARD,nu)) ;

    ans += sum(real(trace(staple_2 * adj(staple_1) )));


    tmp0 = shift(u[eta],BACKWARD,eta) ;
    staple_1 =  adj(tmp0)*shift(u[nu],BACKWARD,eta) *
      (shift(tmp0,FORWARD,nu)) ;

    ans += sum(real(trace(staple_2 * adj(staple_1) )));
    ans += sum(real(trace(staple_0 * adj(staple_1) )));

    // Normalize the planes 
    ans /= Double(Layout::vol()*Nc);
    ans /= 4.0 ; 

    END_CODE();
}




  //! Return the value of the really bent Wilson loops of length 6
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param mu1         direction 1 (Read)
   * \param mu2         direction 2 (Read)
   * \param mu3         direction 3 (Read)
   * \param ans         really bent rectangle  (Write)
   */

  void Wloop_really_bent(const multi1d<LatticeColorMatrix>& u, 
			 int mu1, int mu2, int mu3, 
			 Double & ans)
  {
    START_CODE();


    LatticeColorMatrix tmp0 = shift(u[mu3],FORWARD,mu1)  ;
    LatticeColorMatrix third = shift(tmp0,FORWARD,mu2) ;

    LatticeColorMatrix tmp1   = shift(u[mu1],FORWARD,mu3) ;
    LatticeColorMatrix fourth = shift(tmp1,FORWARD,mu2) ;

    LatticeColorMatrix fifth   = shift(u[mu2],FORWARD,mu3) ;


    ans = sum(real(trace(u[mu1] * shift(u[mu2],FORWARD,mu1)*third*adj(fourth)*adj(fifth)*adj(u[mu3] ))));

    ans /= Double(Layout::vol()*Nc);

    END_CODE();
}




  //! Return the value of the really bent Wilson loops of length 6
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param ans         really bent rectangle  (Write)
   */

  void Wloop_really_bent(const multi1d<LatticeColorMatrix>& u, 
			 Double & ans)
  {
    START_CODE();
    Double ans_local ;  
    
    ans = zero ; 

    Wloop_really_bent(u,0,1,2,ans_local ) ; 
    ans += ans_local ;

    Wloop_really_bent(u,0,1,3,ans_local ) ; 
    ans += ans_local ;

    Wloop_really_bent(u,0,2,3,ans_local ) ; 
    ans += ans_local ;

    Wloop_really_bent(u,1,2,3,ans_local ) ; 
    ans += ans_local ;


    ans /= Double(4.0) ; 

    END_CODE();
}




  //! Return the value of the bent 12 Wilson loops
  /*! averaged over directions.
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param ans         bent rectangle  (Write)
   */

  void Wloop_bent(const multi1d<LatticeColorMatrix>& u, 
		  Double & ans)
  {
    START_CODE();

    ans = 0.0 ; 
    Double local_ans ; 

    Wloop_bent(u,0,1,2,local_ans) ; 
    ans += local_ans ; 

    Wloop_bent(u,0,2,3,local_ans) ; 
    ans += local_ans ; 

    Wloop_bent(u,0,1,3,local_ans) ; 
    ans += local_ans ; 

    ans /= 3 ;

    END_CODE();
}




  //! Return the value of the 22, 23,  Wilson loops
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param plane_plaq_22  plane 2x2 Wilson loops  average (Write)
   * \param plane_plaq_23  plane 2x3 Wilson loops  average (Write)
   */

  void Wloop(const multi1d<LatticeColorMatrix>& u, 
	     multi2d<Double>& plane_plaq_22,
	     multi2d<Double>& plane_plaq_23)
  {
    START_CODE();

    plane_plaq_22.resize(Nd,Nd);
    plane_plaq_23.resize(Nd,Nd);

    // Compute the average Wilson loops
    for(int mu=1; mu < Nd; ++mu)
    {
      for(int nu=0; nu < mu; ++nu)
      {
	/** compute the 2x2 Wilson loop **/
	/* staple_0 = u(x,mu) * u(x+mu,nu)*u_dag(x+nu,mu) */

	LatticeColorMatrix tmp1 = shift(u[nu],FORWARD,mu) ; 
	LatticeColorMatrix tmp2 = shift(tmp1,FORWARD,nu) ; 

	LatticeColorMatrix tmp3 = shift(u[mu],FORWARD,nu) ; 
	LatticeColorMatrix tmp4 = shift(tmp3 ,FORWARD,nu) ; 

	LatticeColorMatrix staple_0 = u[mu]*tmp1*tmp2*adj(tmp4) ;
     	LatticeColorMatrix staple_1 = u[mu]*shift(staple_0,FORWARD,mu)*adj(tmp4) ; 

	LatticeColorMatrix tmp5 = shift(u[nu],FORWARD,nu) ; 

	plane_plaq_22[mu][nu] = sum(real(trace(staple_1 * adj(u[nu] * tmp5) )) );


	/** compute the 2x3 Wilson loop **/
	LatticeColorMatrix tmp2_ = shift(tmp2,FORWARD,nu) ; 
	LatticeColorMatrix tmp4_ = shift(tmp4,FORWARD,nu) ; 
	staple_0 = u[mu] * tmp1*tmp2* tmp2_ * adj(tmp4_) ;
	staple_1 = u[mu] * shift(staple_0,FORWARD,mu) * adj(tmp4_) ;

	LatticeColorMatrix tmp5_ = shift(tmp5,FORWARD,nu) ; 

	plane_plaq_23[mu][nu] = sum(real(trace(staple_1 * adj(u[nu] * tmp5 * tmp5_) )) );

      }
    }

    // Normalize the planes
    for(int mu=1; mu < Nd; ++mu)
      for(int nu=0; nu < mu; ++nu)
      {
	plane_plaq_22[mu][nu] /= Double(Layout::vol()*Nc);
	plane_plaq_22[nu][mu] = plane_plaq_22[mu][nu];

	plane_plaq_23[mu][nu] /= Double(Layout::vol()*Nc);
	plane_plaq_23[nu][mu] = plane_plaq_23[mu][nu];

      }

    END_CODE();
}


  //! Return the value of the average Wilson loops normalized to 1
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param w_plaq      plaquette average (Write)
   * \param s_plaq      space-like plaquette average (Write)
   * \param t_plaq      time-like plaquette average (Write)
   * \param plane_plaq  plane plaquette average (Write)
   */

  void Wloop(const multi1d<LatticeColorMatrix>& u, 
	     Double& w_11, Double& s_11, Double& t_11, multi2d<Double>& plane_plaq_11,
	     Double& w_12, Double& s_12, Double& t_12, multi2d<Double>& plane_plaq_12,
	     Double& w_13, Double& s_13, Double& t_13, multi2d<Double>& plane_plaq_13,
	     Double& w_14, Double& s_14, Double& t_14, multi2d<Double>& plane_plaq_14 )
  {
    START_CODE();

    // Compute plane plaquettes 
    Wloop(u,plane_plaq_11,plane_plaq_12,plane_plaq_13,plane_plaq_14);

    // Compute basic Wilson loops
    w_11 = s_11 = t_11 = zero;
    w_12 = s_12 = t_12 = zero;
    w_13 = s_13 = t_13 = zero;
    w_14 = s_14 = t_14 = zero;

    for(int mu=1; mu < Nd; ++mu)
    {
      for(int nu=0; nu < mu; ++nu)
      {
	Double tmp11 = plane_plaq_11[mu][nu];
	Double tmp12 = plane_plaq_12[mu][nu];
	Double tmp13 = plane_plaq_13[mu][nu];
	Double tmp14 = plane_plaq_14[mu][nu];

	w_11 += tmp11;
	w_12 += tmp12;
	w_13 += tmp13;
	w_14 += tmp14;

	if (mu == tDir() || nu == tDir())
	  {
	    t_11 += tmp11;
	    t_12 += tmp12;
	    t_13 += tmp13;
	    t_14 += tmp14;
	  }
	else 
	  {
	    s_11 += tmp11;
	    s_12 += tmp12;
	    s_13 += tmp13;
	    s_14 += tmp14;
	  }

      }
    }
  
    // Normalize
    w_11 *= 2.0 / Double(Nd*(Nd-1));
    w_12 *= 2.0 / Double(Nd*(Nd-1));
    w_13 *= 2.0 / Double(Nd*(Nd-1));
    w_14 *= 2.0 / Double(Nd*(Nd-1));
  
    if (Nd > 2) 
      {
	s_11 *= 2.0 / Double((Nd-1)*(Nd-2));
	s_12 *= 2.0 / Double((Nd-1)*(Nd-2));
	s_13 *= 2.0 / Double((Nd-1)*(Nd-2));
	s_14 *= 2.0 / Double((Nd-1)*(Nd-2));
      }


    t_11 /= Double(Nd-1);
    t_12 /= Double(Nd-1);
    t_13 /= Double(Nd-1);
    t_14 /= Double(Nd-1);
  
    END_CODE();
  }



  //! Return the value of the average Wilson loops normalized to 1
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param w_plaq      plaquette average (Write)
   * \param s_plaq      space-like plaquette average (Write)
   * \param t_plaq      time-like plaquette average (Write)
   * \param plane_plaq  plane plaquette average (Write)
   */

  void Wloop(const multi1d<LatticeColorMatrix>& u, 
	     Double& w_11, Double& s_11, Double& t_11, multi2d<Double>& plane_plaq_11,
	     Double& w_12, Double& s_12, Double& t_12, multi2d<Double>& plane_plaq_12)
  {
    START_CODE();

    // Compute plane plaquettes 
    Wloop(u,plane_plaq_11,plane_plaq_12);

    // Compute basic Wilson loops
    w_11 = s_11 = t_11 = zero;
    w_12 = s_12 = t_12 = zero;

    for(int mu=1; mu < Nd; ++mu)
    {
      for(int nu=0; nu < mu; ++nu)
      {
	Double tmp11 = plane_plaq_11[mu][nu];
	Double tmp12 = plane_plaq_12[mu][nu];

	w_11 += tmp11;
	w_12 += tmp12;

	if (mu == tDir() || nu == tDir())
	  {
	    t_11 += tmp11;
	    t_12 += tmp12;
	  }
	else 
	  {
	    s_11 += tmp11;
	    s_12 += tmp12;
	  }

      }
    }
  
    // Normalize
    w_11 *= 2.0 / Double(Nd*(Nd-1));
    w_12 *= 2.0 / Double(Nd*(Nd-1));

    if (Nd > 2) 
      {
	s_11 *= 2.0 / Double((Nd-1)*(Nd-2));
	s_12 *= 2.0 / Double((Nd-1)*(Nd-2));
      }


    t_11 /= Double(Nd-1);
    t_12 /= Double(Nd-1);

    END_CODE();
  }



  //! Print the value of the average Wilso loops normalized to 1
  /*!
   * \ingroup glue
   *
   * \param xml        plaquette average (Write)
   * \param xml_group  xml file object ( Read )
   * \param u          gauge field (Read)
   */
  void Wloop(XMLWriter& xml, 
	      const string& xml_group,
	      const multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();

    Double w_11, s_11, t_11 ; multi2d<Double> plane_plaq_11;
    Double w_12, s_12, t_12 ; multi2d<Double> plane_plaq_12;
    Double w_13, s_13, t_13 ; multi2d<Double> plane_plaq_13;
    Double w_14, s_14, t_14 ; multi2d<Double> plane_plaq_14;

    Double w_22, s_22, t_22 ; multi2d<Double> plane_plaq_22;
    Double w_23, s_23, t_23 ; multi2d<Double> plane_plaq_23;

    StopWatch swatch_a;

    swatch_a.start();
    Wloop(u, 
	  w_11, s_11, t_11, plane_plaq_11,
	  w_12, s_12, t_12, plane_plaq_12,
	  w_13, s_13, t_13, plane_plaq_13,
	  w_14, s_14, t_14, plane_plaq_14);

    swatch_a.stop() ;

    double time_in_sec_a  = swatch_a.getTimeInSeconds();
    QDPIO::cout << "Time for 11,12,13, and 14 Wilson loops " << time_in_sec_a << " sec\n" ; 



    StopWatch swatch_b;
    swatch_b.start();
    Wloop(u, 
	  w_22, s_22, t_22, plane_plaq_22,
	  w_23, s_23, t_23, plane_plaq_23) ; 

    swatch_b.stop() ;
    double time_in_sec_b  = swatch_b.getTimeInSeconds();
    QDPIO::cout << "Time for 22 and 23 Wilson loops " <<  time_in_sec_b <<  " sec\n" ; 

    StopWatch swatch_c;
    swatch_c.start();
    Double w_bent1x2 ;
    Wloop_bent(u, w_bent1x2) ;
    swatch_c.stop() ;
    double time_in_sec_c  = swatch_c.getTimeInSeconds();
    QDPIO::cout << "Time for bent 1x2 Wilson loops " <<  time_in_sec_c <<  " sec\n" ; 


    StopWatch swatch_d;
    swatch_d.start();
    Double w_really_bent ;
    Wloop_really_bent(u,w_really_bent) ;
    swatch_d.stop() ;
    double time_in_sec_d  = swatch_d.getTimeInSeconds();
    QDPIO::cout << "Time for really bent Wilson loops " <<  time_in_sec_d <<  " sec\n" ; 



    push(xml, xml_group);


    // --------------------------------------------------

    push(xml, "Wilson_loop_11");
    write(xml, "all_av", w_11);
    write(xml, "space_av", s_11);
    write(xml, "time_av", t_11);

    write(xml, "plane_01_plaq", plane_plaq_11[0][1]);
    write(xml, "plane_02_plaq", plane_plaq_11[0][2]);
    write(xml, "plane_12_plaq", plane_plaq_11[1][2]);
    write(xml, "plane_03_plaq", plane_plaq_11[0][3]);
    write(xml, "plane_13_plaq", plane_plaq_11[1][3]);
    write(xml, "plane_23_plaq", plane_plaq_11[2][3]);

    pop(xml);  // Wilson loop 11

    // --------------------------------------------------

    push(xml, "Wilson_loop_12");
    write(xml, "all_av", w_12);
    write(xml, "space_av", s_12);
    write(xml, "time_av", t_12);

    write(xml, "plane_01_plaq", plane_plaq_12[0][1]);
    write(xml, "plane_02_plaq", plane_plaq_12[0][2]);
    write(xml, "plane_12_plaq", plane_plaq_12[1][2]);
    write(xml, "plane_03_plaq", plane_plaq_12[0][3]);
    write(xml, "plane_13_plaq", plane_plaq_12[1][3]);
    write(xml, "plane_23_plaq", plane_plaq_12[2][3]);

    pop(xml);  // Wilson loop 12

    // --------------------------------------------------

    push(xml, "Wilson_loop_13");
    write(xml, "all_av", w_13);
    write(xml, "space_av", s_13);
    write(xml, "time_av", t_13);

    write(xml, "plane_01_plaq", plane_plaq_13[0][1]);
    write(xml, "plane_02_plaq", plane_plaq_13[0][2]);
    write(xml, "plane_12_plaq", plane_plaq_13[1][2]);
    write(xml, "plane_03_plaq", plane_plaq_13[0][3]);
    write(xml, "plane_13_plaq", plane_plaq_13[1][3]);
    write(xml, "plane_23_plaq", plane_plaq_13[2][3]);

    pop(xml);  // Wilson loop 13

    // --------------------------------------------------

    push(xml, "Wilson_loop_14");
    write(xml, "all_av", w_14);
    write(xml, "space_av", s_14);
    write(xml, "time_av", t_14);

    write(xml, "plane_01_plaq", plane_plaq_14[0][1]);
    write(xml, "plane_02_plaq", plane_plaq_14[0][2]);
    write(xml, "plane_12_plaq", plane_plaq_14[1][2]);
    write(xml, "plane_03_plaq", plane_plaq_14[0][3]);
    write(xml, "plane_13_plaq", plane_plaq_14[1][3]);
    write(xml, "plane_23_plaq", plane_plaq_14[2][3]);

    pop(xml);  // Wilson loop 14

    // --------------------------------------------------

    push(xml, "Wilson_loop_22");
    write(xml, "all_av", w_22);
    write(xml, "space_av", s_22);
    write(xml, "time_av", t_22);

    write(xml, "plane_01_plaq", plane_plaq_22[0][1]);
    write(xml, "plane_02_plaq", plane_plaq_22[0][2]);
    write(xml, "plane_12_plaq", plane_plaq_22[1][2]);
    write(xml, "plane_03_plaq", plane_plaq_22[0][3]);
    write(xml, "plane_13_plaq", plane_plaq_22[1][3]);
    write(xml, "plane_23_plaq", plane_plaq_22[2][3]);

    pop(xml);  // Wilson loop 22

    // --------------------------------------------------



    push(xml, "Wilson_loop_23");
    write(xml, "all_av", w_23);
    write(xml, "space_av", s_23);
    write(xml, "time_av", t_23);

    write(xml, "plane_01_plaq", plane_plaq_23[0][1]);
    write(xml, "plane_02_plaq", plane_plaq_23[0][2]);
    write(xml, "plane_12_plaq", plane_plaq_23[1][2]);
    write(xml, "plane_03_plaq", plane_plaq_23[0][3]);
    write(xml, "plane_13_plaq", plane_plaq_23[1][3]);
    write(xml, "plane_23_plaq", plane_plaq_23[2][3]);

    pop(xml);  // Wilson loop 23

    // --------------------------------------------------

    push(xml, "Wilson_loop_bent_rectangle");
    write(xml, "all_av", w_bent1x2);
    pop(xml);  // Wilson loop bent rectangle


    // --------------------------------------------------

    push(xml, "Wilson_loop_really_bent");
    write(xml, "all_av", w_really_bent);
    pop(xml);  // Wilson loop bent rectangle

    // --------------------------------------------------

    pop(xml);  // xml_group

    END_CODE();
  }

}  // end namespace Chroma
