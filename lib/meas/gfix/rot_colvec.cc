// $Id: rot_colvec.cc,v 3.1 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Rotate a color vector
 */

#include "chromabase.h"
#include "meas/gfix/polar_dec.h"
#include "util/gauge/sunfill.h"

namespace Chroma {

//! Rotate a color vector
/*!
 * \ingroup gfix
 *
 * Rotate a color vector into the form where the component with index
 * s_index is real and the component with larger index are zero.
 * We do this by a series of SU(2) gauge rotations.
 *
 * \param g          Gauge transformation                              (Write)
 * \param psi        Input color vector field                          (Read)
 * \param chi        Output color vector field                         (Write)
 * \param s_index    color index                                       (Read)
 */

void rot_colvec(LatticeColorMatrix& g, 
		const LatticeColorVector& psi,
		LatticeColorVector& chi,
		int s_index)
{
  multi2d<LatticeComplex> g_e(Nc, Nc);
  multi1d<LatticeComplex> psi_e(Nc);
  multi1d<LatticeComplex> chi_e(Nc);

  START_CODE();

  if ( s_index < 0 || (s_index > (Nc-2) && Nc > 1) ||
       (s_index > 0 && Nc == 1) )
    QDP_error_exit("Illigal value for color index", s_index);

            
  if (Nc > 1)
  {
    /* Make each site color vector unit length, to avoid round-off problems
       and do the opposite rescaling at the end. */
    LatticeReal length = sqrt(localNorm2(psi));
    LatticeReal lr1 = 1 / length;

    for(int i=0; i < Nc; ++i)
      psi_e[i] = peekColor(psi, i);

    for(int i=s_index+1; i<Nc; i++)
    {
      chi_e[i] = 0;
      psi_e[i] *= lr1;
    }
    psi_e[s_index] *= lr1;

    for(int i=0; i<s_index; i++)
      chi_e[i] = psi_e[i];

    multi1d<LatticeReal> a(4);
    LatticeComplex lc1;
    LatticeReal lr2;

    for(int i=Nc-2; i>=s_index; i--)
    {
      lr1 = sqrt(localNorm2(psi_e[i]) + localNorm2(psi_e[i+1]));
      lr2 = 1 / lr1;

      lc1 = LatticeComplex(adj(psi_e[i])) * lr2;
      a[0] = real(lc1);
      a[3] = imag(lc1);

      lc1 = LatticeComplex(adj(psi_e[i+1])) * lr2;
      a[2] = real(lc1);
      a[1] = imag(lc1);

      psi_e[i] = cmplx(lr1,0);

      if (i == Nc-2)
      {
	sunFill(g, a, i, all);
      }
      else
      {
	LatticeColorMatrix v1, v2;
	sunFill(v1, a, i, all);
	v2 = v1 * g;
	g = v2;
      }
    }

    chi_e[s_index] = psi_e[s_index] * length;

    for(int i=0; i < Nc; ++i)
      pokeColor(chi, chi_e[i], i);

  }
  else     /* ! Nc > 1 */
  {
    psi_e[0] = peekColor(psi, 0);
    LatticeReal lr1 = sqrt(localNorm2(psi_e[0]));
    chi_e[0] = cmplx(lr1,0);
    pokeColor(chi, chi_e[0], 0);

    g_e[0][0] = adj(psi_e[0]) / lr1;

    pokeColor(g, g_e[0][0], 0, 0);
  }
           
  END_CODE();
}

}  // end namespace Chroma
