// $Id: mesplq.cc,v 3.1 2006-08-24 21:04:31 edwards Exp $
/*! \file
 *  \brief Plaquette measurement
 */

#include "chromabase.h"
#include "meas/glue/mesplq.h"
#include "meas/glue/polylp.h"


namespace Chroma 
{

  // Primitive way for now to indicate the time direction
  static int tDir() {return Nd-1;}

  //! Return the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param plane_plaq  plane plaquette average (Write)
   * \param link        space-time average link (Write)
   */
   

  template<typename Q>
  void MesPlq_t(const multi1d<Q>& u, 
	      multi2d<Double>& plane_plaq, Double& link)
  {
    START_CODE();

    plane_plaq.resize(Nd,Nd);
    link = zero;

    // Compute the average plaquettes
    for(int mu=1; mu < Nd; ++mu)
    {
      for(int nu=0; nu < mu; ++nu)
      {
#if 0
	// This is the the longer way to write the 1-liner in the else clause

	/* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	LatticeColorMatrix tmp_0 = shift(u[nu],FORWARD,mu) * adj(shift(u[mu],FORWARD,nu));

	/* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	LatticeColorMatrix tmp_1 = tmp_0 * adj(u[nu]);

	/* tmp = sum(tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu))) */
	Double tmp = sum(real(trace(u[mu]*tmp_1)));

#else
	// This is the the short way to write the clause in the above block

	/* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	/* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	/* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
	Double tmp = 
	  sum(real(trace(u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu]))));
#endif

	plane_plaq[mu][nu] = tmp;
      }
    }

    // Normalize the planes
    for(int mu=1; mu < Nd; ++mu)
      for(int nu=0; nu < mu; ++nu)
      {
	plane_plaq[mu][nu] /= Double(Layout::vol()*Nc);
	plane_plaq[nu][mu] = plane_plaq[mu][nu];
      }

    // Compute the average link
    for(int mu=0; mu < Nd; ++mu)
      link += sum(real(trace(u[mu])));

    link /= Double(Layout::vol()*Nd*Nc);

    END_CODE();
  }

  void MesPlq(const multi1d<LatticeColorMatrixF3>& u, 
	      multi2d<Double>& plane_plaq, Double& link) 
  {
      MesPlq_t(u,plane_plaq, link);
  }

  void MesPlq(const multi1d<LatticeColorMatrixD3>& u, 
	      multi2d<Double>& plane_plaq, Double& link)
  {
      MesPlq_t(u,plane_plaq, link);
  }

  //! Return the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param w_plaq      plaquette average (Write)
   * \param s_plaq      space-like plaquette average (Write)
   * \param t_plaq      time-like plaquette average (Write)
   * \param plane_plaq  plane plaquette average (Write)
   * \param link        space-time average link (Write)
   */
  template<typename Q>
  void MesPlq_t(const multi1d<Q>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, 
	      multi2d<Double>& plane_plaq,
	      Double& link)
  {
    START_CODE();

    // Compute plane plaquettes and link
    MesPlq(u, plane_plaq, link);

    // Compute basic plaquettes
    w_plaq = s_plaq = t_plaq = zero;

    for(int mu=1; mu < Nd; ++mu)
    {
      for(int nu=0; nu < mu; ++nu)
      {
	Double tmp = plane_plaq[mu][nu];

	w_plaq += tmp;

	if (mu == tDir() || nu == tDir())
	  t_plaq += tmp;
	else 
	  s_plaq += tmp;
      }
    }
  
    // Normalize
    w_plaq *= 2.0 / Double(Nd*(Nd-1));
  
    if (Nd > 2) 
      s_plaq *= 2.0 / Double((Nd-1)*(Nd-2));
  
    t_plaq /= Double(Nd-1);
  
    END_CODE();
  }

  void MesPlq(const multi1d<LatticeColorMatrixF3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, 
	      multi2d<Double>& plane_plaq,
	      Double& link) 
  {
     MesPlq_t(u,w_plaq,s_plaq,t_plaq, plane_plaq, link);
  }

  void MesPlq(const multi1d<LatticeColorMatrixD3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, 
	      multi2d<Double>& plane_plaq,
	      Double& link)
  {
     MesPlq_t(u,w_plaq,s_plaq,t_plaq, plane_plaq, link);
  }

  //! Return the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param u           gauge field (Read)
   * \param w_plaq      plaquette average (Write)
   * \param s_plaq      space-like plaquette average (Write)
   * \param t_plaq      time-like plaquette average (Write)
   * \param link        space-time average link (Write)
   */

  void MesPlq(const multi1d<LatticeColorMatrixF3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, Double& link)
  {
    START_CODE();

    multi2d<Double> plane_plaq;

    MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);

    END_CODE();
  }
 
  void MesPlq(const multi1d<LatticeColorMatrixD3>& u, 
	      Double& w_plaq, Double& s_plaq, Double& t_plaq, Double& link)
  {
    START_CODE();

    multi2d<Double> plane_plaq;

    MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);

    END_CODE();
  }

  //! Print the value of the average plaquette normalized to 1
  /*!
   * \ingroup glue
   *
   * \param xml        plaquette average (Write)
   * \param xml_group  xml file object ( Read )
   * \param u          gauge field (Read)
   */
  template<typename Q>
  void MesPlq_t(XMLWriter& xml, 
	        const string& xml_group,
	        const multi1d<Q>& u)
  {
    START_CODE();

    Double w_plaq, s_plaq, t_plaq, link;
    multi2d<Double> plane_plaq;
    multi1d<DComplex> pollp;

    MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
    polylp(u, pollp);

    push(xml, xml_group);
    write(xml, "w_plaq", w_plaq);
    write(xml, "s_plaq", s_plaq);
    write(xml, "t_plaq", t_plaq);

    if (Nd >= 2)
    {
      write(xml, "plane_01_plaq", plane_plaq[0][1]);
    }

    if (Nd >= 3)
    {
      write(xml, "plane_02_plaq", plane_plaq[0][2]);
      write(xml, "plane_12_plaq", plane_plaq[1][2]);
    }

    if (Nd >= 4)
    {
      write(xml, "plane_03_plaq", plane_plaq[0][3]);
      write(xml, "plane_13_plaq", plane_plaq[1][3]);
      write(xml, "plane_23_plaq", plane_plaq[2][3]);
    }

// This is commented out because it is redundant and takes up space
// in what can be huge XML output files. However, if the info is
// really useful it is trivial to turn the output back on.
//    push(xml, "PlanePlaq");
//    for(int mu=0; mu < Nd-1; ++mu)
//    {
//      for(int nu=mu+1; nu < Nd; ++nu)
//      {
//	push(xml, "elem");
//	write(xml, "mu", mu);
//	write(xml, "nu", nu);
//	write(xml, "plane_plaq", plane_plaq[mu][nu]);
//	pop(xml);  // elem
//      }
//    }
//    pop(xml);  // PlanePlaq

    write(xml, "link", link);
    write(xml, "pollp", pollp);

    pop(xml);  // xml_group

    END_CODE();
  }

 void MesPlq(XMLWriter& xml, 
	     const string& xml_group,
	     const multi1d<LatticeColorMatrixF3>& u)
 {
   MesPlq_t(xml, xml_group, u);
 }

 void MesPlq(XMLWriter& xml, 
	        const string& xml_group,
	        const multi1d<LatticeColorMatrixD3>& u)
 {
   MesPlq_t(xml, xml_group, u);
 }

}  // end namespace Chroma
