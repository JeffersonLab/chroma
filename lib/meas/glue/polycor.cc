// $Id: polycor.cc,v 3.1 2006-08-24 02:33:52 edwards Exp $
/*! \file
 *  \brief Construct Polyakov loop correlation functions from fuzzy links
 */

#include "chromabase.h"
#include "util/gauge/shift2.h"
#include "meas/glue/polycor.h"

namespace Chroma 
{ 

  //! Construct Polyakov loop correlation functions from fuzzy links
  /*!
   * \ingroup glue
   *
   * Construct Polyakov loop correlation functions from fuzzy links at 
   * blocking level bl_level in the directions orthogonal to j_decay
   * and Write them in (pseudo) XML format.
   *
   * \param xml_out      xml file object ( Write )
   * \param xml_group    string used for writing xml data ( Read )
   * \param u            (blocked) gauge field ( Read )
   * \param block_latt    block lattice size ( Read )
   * \param bl_level      blocking level ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   */

  void polycor(XMLWriter& xml_out, const string& xml_group,
	       const multi1d<LatticeColorMatrix>& u,
	       const SftMom& phases,
	       int bl_level)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length = phases.numSubsets();
    int j_decay = phases.getDir();

    int half_length = length/2 + 1;

    LatticeColorMatrix u_dble;
    LatticeColorMatrix pol;
    LatticeColorMatrix tmp_1;

    multi1d< multi1d<DComplex> > poly_loop_slice(Nd-1);
    multi1d<DComplex> poly_slice(length);
    multi1d<Double> poly_corr(half_length);
    DComplex poly_loop;
    DComplex cdummy;
    Double dummy;
    int block_latt;
    int j_poly;
    int mu;
    int n;
    int t;
    int t0;
    int t1;
    int mu_inc;

    push(xml_out, xml_group);

    for(int i=0; i < poly_loop_slice.size(); ++i)
    {
      poly_loop_slice[i].resize(length);
    }

    /* Loop over directions orthogonal to j_decay */
    mu_inc = 0;
    for(mu = 0; mu < Nd; ++mu)
    {
      if( mu != j_decay )
      {
	j_poly = mu;
	mu_inc = mu_inc + 1;

	/* Write fine size of blocked lattice in direction mu */
	block_latt = QDP::Layout::lattSize()[j_poly] / (1 << bl_level);

	/* Construct double links */
	tmp_1 = shift2(u[mu], FORWARD, mu, bl_level);
	u_dble = u[mu] * tmp_1;

	/* Multiply together to get the Polyakov loop */
	pol = u_dble;

	for(n = 2; n <= block_latt/2; ++n)          /* runs over half the linear block size */
	{
	  tmp_1 = shift2(pol, FORWARD, mu, bl_level+1);
	  pol = u_dble * tmp_1;
	}

	/* Multiply last link for odd blocked lattice size */
	if ( (block_latt & 1) != 0 )
	{
	  tmp_1 = shift2(pol, FORWARD, mu, bl_level);
	  pol = u[mu] * tmp_1;
	}

	/* Take the trace and sum up */
	poly_slice = sumMulti(trace(pol), phases.getSet());

	/* Keep a copy of poly_slice */
	for(t = 0; t < length; ++t)
	  poly_loop_slice[mu_inc-1][t] = poly_slice[t];

	/* Initialize Polyakov loop correlations to zero */
	poly_loop = 0;
	poly_corr = 0;

	/* vacuum expectation value */
	for(t = 0;t  < ( length); ++t )
	  poly_loop += poly_slice[t];

	/* And now the Polyakov loop correlation function */
	for(t0 = 0;t0  < ( length); ++t0 )
	  for(t = 0;t  < ( half_length); ++t )
	  {
	    t1 = (t + t0) % length;
	    cdummy = poly_slice[t0] * adj(poly_slice[t1]);

	    dummy = real(cdummy);
	    poly_corr[t] += dummy;
	  }

	/* Normalize */
	dummy = Double(1) / Double(QDP::Layout::vol());
	poly_loop = poly_loop * dummy;

	for(t = 0;t  < ( half_length); ++t )
	  poly_corr[t] = poly_corr[t] * dummy;

	/* Finally Write Polyakov loop correlations in XML format */
	push(xml_out,"Poly_loop_correlation");
	write(xml_out, "j_decay", j_decay);
	write(xml_out, "bl_level", bl_level);
	write(xml_out, "j_poly", j_poly);
	write(xml_out, "poly_loop", poly_loop);
	write(xml_out, "poly_corr", poly_corr);
	pop(xml_out);
      }       /* end if( mu != j_decay ) */
    }         /* end loop over mu */

    push(xml_out,"Poly_loop_slice");
    write(xml_out, "j_decay", j_decay);
    write(xml_out, "bl_level", bl_level);
    push(xml_out, "Slices");
    for(int i=0; i < poly_loop_slice.size(); ++i)
    {
      push(xml_out, "elem");
      write(xml_out, "slice", i);
      write(xml_out, "poly_loop_slice", poly_loop_slice[i]);
      pop(xml_out);
    }
    pop(xml_out);
    pop(xml_out);

    pop(xml_out);

    END_CODE();
  }

}
