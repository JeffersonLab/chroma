/*! \file
 *  \brief Construct an instanton or anti-instanton configuration in singular gauge. 
 *
 * Arguments:
 *
 * \param u_inst          Gauge field                   (Write)
 * \param center          location of instanton center ( Read )
 * \param rho             size of instanton ( Read )
 * \param su2_index       SU(2) subgroup index ( Read )
 * \param sign            instanton (1) or anti-instanton ( Read ) 
*/

#include "chromabase.h"
#include "util/gauge/sunfill.h"

namespace Chroma
{
  void instanton(multi1d<LatticeColorMatrix>& u_inst, const multi1d<Real>& center, Real rho, int su2_index, int sign)
  {
    multi1d<LatticeReal> rel_coord(Nd);
    LatticeReal r2;
    LatticeReal phimu;
    LatticeReal phimus;
    LatticeReal a_0;
    LatticeReal a_1;
    LatticeReal a_2;
    LatticeReal a_3;
    LatticeReal tmp1;
    LatticeReal tmp2;
    LatticeReal rel_tmp;
    LatticeInteger my_coord;
    LatticeInteger coord_sum;
    LatticeInteger litmp;
    LatticeBoolean lbtmp;
    int mu;
    int itmp;

    Real ftmp;
    Real ftmp2;

    /* Include any hooks in the start macro */
    START_CODE();			/*# May include declarations and/or code. */

    /* Start executable code */

    if ( Nd != 4 )
      QDP_error_exit("Instanton construction requires Nd=4", Nd);
    if (Nc == 1)
      QDP_error_exit("Instanton construction requires Nc>1", Nc);

    Real rho2 = rho * rho;

                              
    for(int cb = 0; cb < 2; cb++)
    {
      r2 = zero;
      coord_sum = zero;

      /* Compute coordintes relative to the center */
      for(mu = 1; mu < Nd; mu++)
      {
	my_coord = Layout::latticeCoordinate(mu);
	coord_sum += my_coord;
	rel_coord[mu] = my_coord;
	ftmp = center[mu];
	tmp1 = ftmp;
	rel_coord[mu] -= tmp1;
	ftmp = QDP::Layout::lattSize()[mu];
	ftmp2 = 0.5 * ftmp;
	tmp2 = ftmp;

//	lbtmp = rel_coord[mu] > ftmp2;
//	copymask(rel_coord[mu], lbtmp, tmp2, SUBTRACT);
//	tmp1 = -rel_coord[mu];
//	lbtmp = tmp1 > ftmp2;
//	copymask(rel_coord[mu], lbtmp, tmp2, ADD);

	// New version
	lbtmp = rel_coord[mu] > ftmp2;
	rel_tmp = where(lbtmp, tmp2, LatticeReal(zero));
	rel_coord[mu] -= rel_tmp;
	tmp1 = -rel_coord[mu];
	lbtmp = tmp1 > ftmp2;
	rel_tmp = where(lbtmp, tmp2, LatticeReal(zero));
	rel_coord[mu] += rel_tmp;

	r2 += rel_coord[mu] * rel_coord[mu];
      }

      /* For the zeroth direction, we have to determine the real coordinate */
      litmp = cb;
      coord_sum += litmp;
      itmp = 2;
      my_coord = coord_sum % itmp;
      mu = 0;
      litmp = Layout::latticeCoordinate(mu);
      my_coord += litmp * itmp;
      rel_coord[mu] = my_coord;
      ftmp = center[mu];
      tmp1 = ftmp;
      rel_coord[mu] -= tmp1;
      ftmp = QDP::Layout::lattSize()[mu];
      ftmp2 = 0.5 * ftmp;
      tmp2 = ftmp;
      lbtmp = rel_coord[mu] > ftmp2;
      rel_tmp = where(lbtmp, tmp2, LatticeReal(zero));
      rel_coord[mu] -= rel_tmp;
      tmp1 = -rel_coord[mu];
      lbtmp = tmp1 > ftmp2;
      rel_tmp = where(lbtmp, tmp2, LatticeReal(zero));
      rel_coord[mu] += rel_tmp;

      r2 += rel_coord[mu] * rel_coord[mu];

      for(mu = 0; mu < Nd; mu++)
      {
	tmp1 = r2;
	tmp2 = tmp1;
	tmp2 += rel_coord[mu];
	tmp1 -= rel_coord[mu] * rel_coord[mu];
	tmp1 = sqrt(tmp1);
	tmp2 = tmp1 / tmp2;
	tmp2 = atan(tmp2);
	phimu = tmp2 / tmp1;

	tmp1 = rho2;
	tmp1 += r2;
	tmp2 = tmp1;
	tmp2 += rel_coord[mu];
	tmp1 -= rel_coord[mu] * rel_coord[mu];
	tmp1 = sqrt(tmp1);
	tmp2 = tmp1 / tmp2;
	tmp2 = atan(tmp2);
	tmp2 = tmp2 / tmp1;
	phimu -= tmp2;

	if (sign == PLUS)
	{
	  phimus = phimu;
	}
	else
	{
	  phimus = -phimu;
	}

	switch (mu)
	{
	case 0:
	  a_1 = -(rel_coord[3] * phimus);
	  a_2 = rel_coord[2] * phimu;
	  a_3 = -(rel_coord[1] * phimu);
	  break;

	case 1:
	  a_1 = -(rel_coord[2] * phimu);
	  a_2 = -(rel_coord[3] * phimus);
	  a_3 = rel_coord[0] * phimu;
	  break;

	case 2:
	  a_1 = rel_coord[1] * phimu;
	  a_2 = -(rel_coord[0] * phimu);
	  a_3 = -(rel_coord[3] * phimus);
	  break;

	case 3:
	  a_1 = rel_coord[0] * phimus;
	  a_2 = rel_coord[1] * phimus;
	  a_3 = rel_coord[2] * phimus;
	  break;

	default:
	  QDP_error_exit("mu out of range", mu);
	  break;
	}

	tmp1  = a_1 * a_1;
	tmp1 += a_2 * a_2;
	tmp1 += a_3 * a_3;
	tmp1 = sqrt(tmp1);
	tmp2 = 1;
	ftmp = 1.0e-5;
	lbtmp = tmp1 > ftmp;
	copymask(tmp2, lbtmp, tmp1);
	a_0 = cos(tmp2);
	tmp1 = sin(tmp2);
	tmp1 = tmp1 / tmp2;
	a_1 = a_1 * tmp1;
	a_2 = a_2 * tmp1;
	a_3 = a_3 * tmp1;

	/* Where tmp1 < FUZZ above, fill a_i with (1,0,0,0) */
	tmp1 = 1;
	tmp2 = zero;
	lbtmp = !lbtmp;
	copymask(a_0, lbtmp, tmp1);
	copymask(a_1, lbtmp, tmp2);
	copymask(a_2, lbtmp, tmp2);
	copymask(a_3, lbtmp, tmp2);

/*      SUN_FILL(a_0, a_1, a_2, a_3, su2_index, u_inst(cb,mu)); */
	multi1d<LatticeReal> a(4);
	a[0] = a_0;
	a[1] = a_1;
	a[2] = a_2;
	a[3] = a_3;

	sunFill(u_inst[mu], a, su2_index, QDP::rb[cb]);
      }
    }
                              
    /* Close out any other code */
    END_CODE();
  }
}
