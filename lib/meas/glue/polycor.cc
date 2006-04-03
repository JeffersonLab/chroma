// $Id: polycor.cc,v 3.0 2006-04-03 04:58:58 edwards Exp $

#error "NOT FULLY CONVERTED"


/* Construct Polyakov loop correlation functions from fuzzy links at  */
/* blocking level bl_level in the directions orthogonal to j_decay  */
/* and Write them in (pseudo) NAMELIST format. */

/* u -- (blocked) gauge field ( Read ) */
/* block_latt -- block lattice size ( Read ) */
/* bl_level -- blocking level ( Read ) */
/* j_decay -- direction of exponential decay ( Read ) */

include(types.mh)
SUBROUTINE(polycor, u, bl_level, j_decay)


multi1d<LatticeColorMatrix> `u');
int bl_level;
int j_decay;

{ /* Local Variables */
LatticeColorMatrix u_dble;
LatticeColorMatrix pol;
LatticeColorMatrix tmp_1;
LatticeComplex pol_trace;

multi2d<DComplex> poly_loop_slice(length, Ndm1);
multi1d<DComplex> poly_slice(length);
multi1d<Double> poly_corr(half_length);
DComplex poly_loop;
DComplex cdummy;
Double dummy;
int block_latt;
int length;
int half_length;
int j_poly;
int cb;
int cbs;
int mu;
int n;
int t;
int t0;
int t1;
int mu_inc;
int Ndm1;

include(COMMON_DECLARATIONS)
START_CODE();

length = nrow[j_decay];
half_length = length/2 + 1;
Ndm1 = Nd - 1;




/* Loop over directions orthogonal to j_decay */
mu_inc = 0;
for(mu = 0;mu  < ( Nd); ++mu )
{
  if( mu != j_decay )
  {
    j_poly = mu;
    mu_inc = mu_inc + 1;

/* Initialize poly_slice */
    poly_slice = 0;

/* Write fine size of blocked lattice in direction mu */
    block_latt = nrow[j_poly] / INTEGER_LSHIFT_FUNCTION(1,bl_level);

/* Loop over checkerboards */
    for(cb = 0;cb  <= ( 1); ++cb )
    {

/* Determine checkboard of nearest neighbors at blocking level bl_level */
      if( bl_level == 0 )
        cbs = 1 - cb;
      else
        cbs = cb;

      /* Construct double links */
      NEIGHBOUR(u(cbs,mu), tmp_1, cb, FORWARD, mu, bl_level);
      u_dble = u[mu][cb] * tmp_1;

/* Multiply together to get the Polyakov loop */
      pol = u_dble;

      for(n = 2;n  <= ( block_latt/2 ); ++n )          /* runs over half the linear block size */
      {
        NEIGHBOUR(pol, tmp_1, cb, FORWARD, mu, bl_level+1);
        pol = u_dble * tmp_1;
      }

      /* Multiply last link for odd blocked lattice size */
      if ( INTEGER_MOD_FUNCTION(block_latt,2) != 0 )
      {
        NEIGHBOUR(pol, tmp_1, cbs, FORWARD, mu, bl_level);
        pol = u[mu][cbs] * tmp_1;
      }

      /* Take the trace and sum up */
      pol_trace = trace(pol);
      poly_slice += sumMulti(pol_trace, timeslice);

    }     /* end loop over cb */

/* Keep a copy of poly_slice */
    for(t = 0;t  < ( length); ++t )
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

      t1 = INTEGER_MOD_FUNCTION(t + t0, length);
      cdummy = poly_slice[t0] * adj(poly_slice[t1]);

      dummy = real(cdummy);
      poly_corr[t] += dummy;
    }

/* Normalize */
    dummy = TO_DOUBLE(1) / TO_DOUBLE(vol);
    poly_loop = poly_loop * dummy;

    for(t = 0;t  < ( half_length); ++t )
      poly_corr[t] = poly_corr[t] * dummy;

/* Finally Write Polyakov loop correlations in NAMELIST format */
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
write(xml_out, "poly_loop_slice", poly_loop_slice);
pop(xml_out);




END_CODE();
}(Nd)
