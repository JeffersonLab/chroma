// $Id: gluecor.cc,v 1.1 2005-02-17 02:52:11 edwards Exp $

#error "NOT FULLY CONVERTED"

/* Construct 0++, 2++ and 1+- glueball correlation functions from */
/* fuzzy links at blocking level bl_level and Write them in  */
/* NAMELIST format. */

/* Warning: this works only for Nd = 4 ! */

/* u -- (blocked) gauge field ( Read ) */
/* bl_level -- blocking level ( Read ) */
/* j_decay -- direction of exponential decay ( Read ) */

include(types.mh)

SUBROUTINE(gluecor, u, bl_level, j_decay)

multi1d<LatticeColorMatrix> u(Nd);
int bl_level;
int j_decay;

{ /* Local Variables */
include(COMMON_DECLARATIONS)

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
multi1d<Double> plaq_r(length);
multi1d<Double> plaq_i(length);
Double vac0;
Double dummy;
multi2d<DComplex> plaq(length, 3);
int length;
int half_length;
int cb;
int cbs;
int mu;
int nu;
int plane;
int t;
int t0;
int t1;

START_CODE();

if( Nd != 4 )
  QDP_error_exit("Nd for glueball construction has to be 4 but: ", Nd);

length = nrow[j_decay];
half_length = length/2 + 1;


/* Initialize glueball operators to zero */
op0 = 0;
op1 = 0;
op2 = 0;
plaq = 0;



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

/* Loop over checkerboards */
        for(cb = 0; cb < 2; ++cb)
        {

/* Determine checkboard of nearest neighbors at blocking level bl_level */
          if( bl_level == 0 )
            cbs = 1 - cb;
          else
            cbs = cb;

/* Construct (block) plaquettes */
          /* tmp_1(x) = u(x+mu*2^bl_level,nu) */
          NEIGHBOUR(u(cbs,nu), tmp_1, cb, FORWARD, mu, bl_level);

          /* tmp_2(x) = u(x+nu*2^bl_level,mu) */
          NEIGHBOUR(u(cbs,mu), tmp_2, cb, FORWARD, nu, bl_level);

          /* tmp_3 = tmp_1 * tmp_2_dagger */
          tmp_3 = tmp_1 * adj(tmp_2);

          /* tmp_1 = tmp_3 * u_dagger(x,nu) */
          tmp_1 = tmp_3 * adj(u[nu][cb]);

          /* cplaq_tmp = Tr(u(x,mu) * tmp_1) */
          cplaq_tmp = trace(u[mu][cb] * tmp_1);

          /* cplaq = slice-wise sum of cplaq_tmp */
          cplaq = sumMulti(cplaq_tmp, timeslice);

          /* plaq_r = Re(Tr(cplaq_tmp)) */
          plaq_r = real(cplaq);

          /* plaq_i = Im(Tr(cplaq_tmp)) */
          plaq_i = imag(cplaq);

/* Make the glueball operators from the (block) plaquettes */
          switch(plane)
          {
            case 1:
              for(t = 0;t  < ( length); ++t )
              {
                op0[t] += plaq_r[t];
                op1[0][0][t] += plaq_i[t];
                op1[1][1][t] += plaq_i[t];
                op2[0][t] += plaq_r[t];
                op2[1][t] -= plaq_r[t];
              }

            break;
            case 2:
              for(t = 0;t  < ( length); ++t )
              {
                op0[t] += plaq_r[t];
                op1[0][2][t] += plaq_i[t];
                op1[1][0][t] += plaq_i[t];
                op2[0][t] -= plaq_r[t];
                op2[2][t] -= plaq_r[t];
              }

            break;
            case 3:
              for(t = 0;t  < ( length); ++t )
              {
                op0[t] += plaq_r[t];
                op1[0][1][t] += plaq_i[t];
                op1[1][2][t] += plaq_i[t];
                op2[1][t] += plaq_r[t];
                op2[2][t] += plaq_r[t];
              }

            break;
            default:
              QDP_error_exit("Too many orthogonal planes: ", plane);
          }

          for(t = 0;t  < ( length); ++t )
          {
            plaq[plane-1][t] += cmplx(plaq_r[t],plaq_i[t]);
          }

        }   /* end loop over cb */
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
for(t = 0;t  < ( length); ++t )
  vac0 += op0[t];

/* And now the glueball correlation functions */
for(t0 = 0;t0  < ( length); ++t0 )
for(t = 0;t  < ( half_length); ++t )
{
  t1 = INTEGER_MOD_FUNCTION(t + t0, length);
  glue0[t] += op0[t0] * op0[t1];

  for(nu = 0;nu  <= ( 2); ++nu )
  {
    glue1[t] += op1[0][nu][t0] * op1[0][nu][t1];
    glue1[t] += op1[1][nu][t0] * op1[1][nu][t1];
    glue2[t] += op2[nu][t0] * op2[nu][t1];
  }
}

/* Normalize */
vac0 = vac0 / (TO_DOUBLE(vol));
dummy = TO_DOUBLE(1) / TO_DOUBLE(vol);
for(t = 0;t  < ( half_length); ++t )
{
  glue0[t] = glue0[t] * dummy;
  glue1[t] = glue1[t] * dummy;
  glue2[t] = glue2[t] * dummy;
}

/* Finally Write glueball correlations in NAMELIST format */
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
write(xml_out, "plaq", plaq);
pop(xml_out);



END_CODE();
}
