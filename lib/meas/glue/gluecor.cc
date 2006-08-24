// $Id: gluecor.cc,v 3.1 2006-08-24 02:33:52 edwards Exp $
/*! \file
 *  \brief Construct 0++, 2++ and 1+- glueball correlation functions from fuzzy links
 */

#include "chromabase.h"
#include "util/gauge/shift2.h"
#include "meas/glue/gluecor.h"

namespace Chroma 
{ 

  //! Construct 0++, 2++ and 1+- glueball correlation functions from fuzzy links
  /*! 
   * \ingroup glue
   *
   * Construct 0++, 2++ and 1+- glueball correlation functions from
   * fuzzy links at blocking level bl_level and Write them in 
   * XML format.
   *
   * Warning: this works only for Nd = 4 !
   *
   * \param xml_out       xml file object ( Write )
   * \param xml_group     string used for writing xml data ( Read )
   * \param u             (blocked) gauge field ( Read )
   * \param bl_level      blocking level ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   */

  void gluecor(XMLWriter& xml_out, const string& xml_group,
	       const multi1d<LatticeColorMatrix>& u, 
	       const SftMom& phases,
	       int bl_level)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length = phases.numSubsets();
    int j_decay = phases.getDir();

    int half_length = length/2 + 1;

    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;
    LatticeComplex cplaq_tmp;

    multi1d<DComplex> cplaq(length);
    multi1d<Double> glue0(half_length);
    multi1d<Double> glue1(half_length);
    multi1d<Double> glue2(half_length);
    multi1d<Double> op0(length);
    multi3d<Double> op1(length, 3, 2);
    multi2d<Double> op2(length, 3);
    Double vac0;
    Double dummy;
    multi1d< multi1d<DComplex> > plaq(3);
    int mu;
    int nu;
    int plane;
    int t;
    int t0;
    int t1;

    if( Nd != 4 )
      QDP_error_exit("Nd for glueball construction has to be 4 but: ", Nd);

    for(int i=0; i < 3; ++i)
    {
      plaq[i].resize(length);
      plaq[i] = 0;
    }

    /* Initialize glueball operators to zero */
    op0 = 0;
    op1 = 0;
    op2 = 0;

    plane = 0;

    /* Construct glueball operators */
    /* Loop over first direction (!= j_decay) */
    for(mu = 0; mu <= Nd-2; ++mu)
    {
      if( mu != j_decay )
      {
	/* Loop over second direction (!= j_decay) */
	for(nu = mu+1; nu < Nd; ++nu)
	{
	  if( nu != j_decay )
	  {
	    plane = plane + 1;

	    /* Construct (block) plaquettes */
	    /* tmp_1(x) = u(x+mu*2^bl_level,nu) */
	    tmp_1 = shift2(u[nu], FORWARD, mu, bl_level);

	    /* tmp_2(x) = u(x+nu*2^bl_level,mu) */
	    tmp_2 = shift2(u[mu], FORWARD, nu, bl_level);

	    /* tmp_3 = tmp_1 * tmp_2_dagger */
	    tmp_3 = tmp_1 * adj(tmp_2);

	    /* tmp_1 = tmp_3 * u_dagger(x,nu) */
	    tmp_1 = tmp_3 * adj(u[nu]);

	    /* cplaq_tmp = Tr(u(x,mu) * tmp_1) */
	    cplaq_tmp = trace(u[mu] * tmp_1);

	    /* cplaq = slice-wise sum of cplaq_tmp */
	    cplaq = sumMulti(cplaq_tmp, phases.getSet());

	    /* Make the glueball operators from the (block) plaquettes */
	    switch(plane)
	    {
	    case 1:
	      for(t = 0;t  < ( length); ++t )
	      {
		op0[t] += real(cplaq[t]);
		op1[0][0][t] += imag(cplaq[t]);
		op1[1][1][t] += imag(cplaq[t]);
		op2[0][t] += real(cplaq[t]);
		op2[1][t] -= real(cplaq[t]);
	      }

	      break;
	    case 2:
	      for(t = 0;t  < ( length); ++t )
	      {
		op0[t] += real(cplaq[t]);
		op1[0][2][t] += imag(cplaq[t]);
		op1[1][0][t] += imag(cplaq[t]);
		op2[0][t] -= real(cplaq[t]);
		op2[2][t] -= real(cplaq[t]);
	      }

	      break;
	    case 3:
	      for(t = 0;t  < ( length); ++t )
	      {
		op0[t] += real(cplaq[t]);
		op1[0][1][t] += imag(cplaq[t]);
		op1[1][2][t] += imag(cplaq[t]);
		op2[1][t] += real(cplaq[t]);
		op2[2][t] += real(cplaq[t]);
	      }

	      break;
	    default:
	      QDP_error_exit("Too many orthogonal planes: ", plane);
	    }

	    for(t = 0;t  < ( length); ++t )
	    {
	      plaq[plane-1][t] += cplaq[t];
	    }

	  }     /* end if( nu != j_decay ) */
	}       /* end loop over nu */
      }         /* end if( mu != j_decay ) */
    }           /* end loop over mu */



    /* Initialize glueball correlations to zero */
    vac0 = 0;
    glue0 = 0;
    glue1 = 0;
    glue2 = 0;

    /* 0++ vacuum expectation value */
    for(t = 0; t < length; ++t)
      vac0 += op0[t];

    /* And now the glueball correlation functions */
    for(t0 = 0; t0 < length; ++t0)
      for(t = 0; t < half_length; ++t)
      {
	t1 = (t + t0) % length;
	glue0[t] += op0[t0] * op0[t1];

	for(nu = 0;nu  <= ( 2); ++nu )
	{
	  glue1[t] += op1[0][nu][t0] * op1[0][nu][t1];
	  glue1[t] += op1[1][nu][t0] * op1[1][nu][t1];
	  glue2[t] += op2[nu][t0] * op2[nu][t1];
	}
      }

    /* Normalize */
    vac0 /= Double(QDP::Layout::vol());
    dummy = Double(1) / Double(QDP::Layout::vol());
    for(t = 0;t  < ( half_length); ++t )
    {
      glue0[t] = glue0[t] * dummy;
      glue1[t] = glue1[t] * dummy;
      glue2[t] = glue2[t] * dummy;
    }

    /* Finally Write glueball correlations in NAMELIST format */
    push(xml_out, xml_group);

    push(xml_out,"Glueball_0pp");
    write(xml_out, "j_decay", j_decay);
    write(xml_out, "bl_level", bl_level);
    write(xml_out, "vac0", vac0);
    write(xml_out, "glue0", glue0);
    pop(xml_out);

    push(xml_out,"Glueball_1pm");
    write(xml_out, "j_decay", j_decay);
    write(xml_out, "bl_level", bl_level);
    write(xml_out, "glue1", glue1);
    pop(xml_out);

    push(xml_out,"Glueball_2pp");
    write(xml_out, "j_decay", j_decay);
    write(xml_out, "bl_level", bl_level);
    write(xml_out, "glue2", glue2);
    pop(xml_out);

    push(xml_out,"Glueball_plaq");
    write(xml_out, "j_decay", j_decay);
    write(xml_out, "bl_level", bl_level);
    push(xml_out, "Planes");
    for(int i=0; i < plaq.size(); ++i)
    {
      push(xml_out, "elem");
      write(xml_out, "plane", i);
      write(xml_out, "plaq", plaq[i]);
      pop(xml_out);
    }
    pop(xml_out);
    pop(xml_out);

    pop(xml_out);

    END_CODE();
  }

}
